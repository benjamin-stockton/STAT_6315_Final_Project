library(mvtnorm)
library(circular)
library(MASS)
library(latex2exp)
library(truncnorm)
library(MCMCpack)
library(RColorBrewer)


# Abe-Ley Distribution random variable draw

dabeley <- function(x,theta, mu, kappa, lambda, alpha, beta) {
    const <- alpha * beta^alpha / (2*pi * cosh(kappa))
    wc <- (1 + lambda * sin(theta - mu))
    weib <- x^(alpha - 1) * exp(-(beta * x)^alpha * (1 - tanh(kappa) * cos(theta - mu)))
    return(const * wc * weib)
}

rabeley <- function(n, mu, kappa, lambda, alpha, beta) {
    dat <- matrix(numeric(2*n), ncol = 2)
    colnames(dat) <- c("theta", "x")
    for (i in 1:n) {
        t1 <- rwrappedcauchy(1, mu = circular(mu), tanh(kappa/2))
        t <- ifelse(t1 < (1 + lambda * sin(t1 - mu)) / 2, t1, -t1) %% (2*pi)
        scale_par <- beta * (1 - tanh(kappa) * cos(t - mu))^(-1/alpha)
        x <- rweibull(1, shape = alpha, scale = scale_par)
        dat[i,] <- c(t,x)
    }
    return(dat)
}

rssvm <- function(n, mu, kappa, lambda) {
    t1 <- rvonmises(n, circular(mu), kappa)
    # u <- runif(n)
    t <- sapply(t1, FUN = function(t) {return(ifelse(runif(1) < (1 + lambda * sin(t-mu))/2, t, -t))}) %% (2*pi)
    return(t)
}

dbeta_rescale <- function(x, shape1 = 1, shape2 = 1) {
    return(dbeta((x/2 + 1/2), shape1 = shape1, shape2 = shape2)/2)
}

dssvm <- function(theta, mu, kappa, lambda) {
    const <- (2*pi * I.0(kappa))^-1
    skew <- (1 + lambda * sin(theta - mu))
    kern <- exp(kappa*cos(theta - mu))
    return(const * skew * kern)
}

dsswc <- function(theta, mu, kappa, lambda) {
    const <- (1 - tanh(kappa/2)^2) / (2*pi)
    num <- 1 + lambda * sin(theta - mu)
    denom <- 1 + tanh(kappa/2)^2 - 2*tanh(kappa/2) * cos(theta - mu)
    return(const * num / denom)
}

dmodweib <- function(x, alpha, beta, kappa) {
    const <- alpha * beta^alpha * I.0(x^alpha * beta^alpha * tanh(kappa)) / cosh(kappa) 
    kern <- x^(alpha-1) * exp(-(beta*x)^alpha)
    return(const * kern)
}


# Mixture Abe-Ley Simulated Data
dabeley_mixture <- function(x,theta, mix_prop, mu, kappa, lambda, alpha, beta) {
    K <- length(mix_prop)
    
    d <- 0
    
    for (k in 1:K) {
        d <- d + mix_prop[k] * dabeley(x = x, theta = theta, mu = mu[k], kappa = kappa[k],
                                       lambda = lambda[k], alpha = alpha[k], beta = beta[k])
    }
    return(d)
}

rabeley_mixture <- function(N, mix_prop, mu, kappa, lambda, alpha, beta) {
    K <- length(mix_prop)
    dat <- matrix(numeric(2*N), ncol = 2)
    colnames(dat) <- c("theta", "x")
    for (i in 1:N) {
        k <- sample(1:K, size = 1, prob = mix_prop, replace = T)
        tmp <- rabeley(1, mu[k], kappa[k], lambda[k], alpha[k], beta[k])
        dat[i,] <- tmp
    }
    return(dat)
}

rssvm_mixture <- function(N, mix_prop, mu, kappa, lambda) {
    K <- length(mix_prop)
    t <- numeric(N)
    for (i in 1:N) {
        k <- sample(1:K, size = 1, prob = mix_prop, replace = T)
        t[i] <- rssvm(1, mu[k], kappa[k], lambda[k])
    }
    return(t)
}

dsswc_mixture <- function(theta, mu, kappa, lambda, tau) {
    K <- length(tau)
    d <- 0
    
    for (k in 1:K) {
        d <- d + tau[k] * dsswc(theta, mu = mu[k], kappa = kappa[k], lambda[k])
    }
    return(d)
}

dmodweibull_mixture <- function(x, alpha, beta, kappa, tau) {
    K <- length(tau)
    d <- 0
    
    for (k in 1:K) {
        d <- d + tau[k] * dmodweib(x, alpha = alpha[k], beta = beta[k], kappa = kappa[k])
    }
    return(d)
}

# Plotting Functions

cyl_density_plot <- function(dat, mean_pars, mix_prop, main = "Joint Distibution on the Cylinder", xlab = "$\\theta$", ylab = "X") {
    x_seq <- seq(from = min(dat[,2]), to = max(dat[,2]), length.out = 50)
    t_seq <-  seq(from = min(dat[,1]), to = max(dat[,1]), length.out = 50)
    
    z <- t(outer(X = x_seq, Y = t_seq, FUN = dabeley_mixture,
               mix_prop = mix_prop,
               alpha = mean_pars[,1], beta = mean_pars[,2],
               kappa = mean_pars[,3], mu = mean_pars[,4], lambda = mean_pars[,5]))
    
    plot(dat[,1], dat[,2], pch = 16, col = "black",
         main = TeX(main), xlab = TeX(xlab), ylab = TeX(ylab))
    contour(x = t_seq, y = x_seq, z = z, lwd = 2, add = TRUE, nlevels = 10,
            col = hcl.colors(10, "Spectral"))
}

hist_density <- function(dat, main = "Distribution of X", xlab = "", ylab = "Density", lhist = 10,...) {
    dx <- density(dat)
    hx <- hist(dat, plot = F, breaks = seq(from = min(dat), to = max(dat),
                                           length.out = lhist))
    plot(hx, col="grey", main = TeX(main), xlab = TeX(xlab), ylab = TeX(ylab),...)
    lines(x = dx$x, y = dx$y * length(dat) * diff(hx$breaks)[1], lwd = 2)
}

joint_dist_plot <- function(dat, mean_pars, mix_prop, main = "", xlab = "$\\theta$", ylab = "X", lhist = 20) {
    layout(matrix(c(2,0,1,3), nrow = 2, byrow = T),
           widths = c(3,1), heights = c(1,3), respect = T)
    
    par. <- par(mar = c(4, 4, 1,1), oma = rep(.5, 4))
    
    cyl_density_plot(dat = dat, mean_pars = mean_pars, mix_prop = mix_prop, main = main, xlab = xlab, ylab = ylab)
    
    t_hist <- hist(dat[,1], plot = F, breaks = seq(from = min(dat[,1]), to = max(dat[,1]),
                                                   length.out = lhist))
    x_hist <- hist(dat[,2], plot = F, breaks = seq(from = min(dat[,2]), to = max(dat[,2]),
                                                   length.out = lhist))
    
    td <- density(dat[,1])
    xd <- density(dat[,2])
    
    fit_td <- dsswc_mixture(td$x, mu = mean_pars[,4], kappa = mean_pars[,3], lambda = mean_pars[,5], tau = mix_prop)
    
    fit_xd <- dmodweibull_mixture(xd$x, alpha = mean_pars[,1], beta = mean_pars[,2], kappa = mean_pars[,3], tau = mix_prop)
    
    par(mar = c(0, 4, 0,0))
    plot(t_hist, col="grey", main = "", xlab = "", ylab = "", axes = F, add = F)
    lines(x = td$x, y = td$y * length(dat[,1]) * diff(t_hist$breaks)[1], lwd = 2)
    lines(x = td$x, y = fit_td * length(dat[,1]) * diff(t_hist$breaks)[1], lwd = 1, col = "blue", lty = "dashed")
    # * length(dat[,1]) * diff(t_hist$breaks)[1]
    
    par(mar = c(4,0,0,0))
    barplot(x_hist$density, axes = F,
            xlim = c(0, max(x_hist$density)),
            space = 0, horiz = T)
    lines(x = xd$y * length(dat[,2]) * diff(x_hist$breaks)[1], y = xd$x, lwd = 2)
    lines(x = fit_xd * length(dat[,2]) * diff(x_hist$breaks)[1], y = xd$x, lwd = 1, col = "blue", lty = "dashed")
    # * length(dat[,2]) * diff(x_hist$breaks)[1]
    # legend("topright", legend = c("Generic Density", "Fitted Abe-Ley Density"), col = c("black", "blue"), lty = c("solid", "dashed"), lwd = c(2,1))
    par(par.)
    
    par(mfrow = c(1,1))
}

plot_tracestack <- function(mcmc_pars, K, ctrl) {
    if (K == 1) {
        m_props <- matrix(mcmc_pars[,6], ncol = K)
        mix_props <- apply(m_props, 2, mean)
        # print(mix_props)
        iter <- 1:nrow(m_props)
        
        param <- c("alpha", "beta", "kappa", "mu", "lambda", "tau")
        post_means <- matrix(numeric(5*K), ncol = 5)
        for (i in 1:5) {
            post_means[,i] <- apply(matrix(mcmc_pars[,i], ncol = K), 2, mean)
        }
        # print(post_means)
        par(mfrow = c(2,3))
        cols <- brewer.pal(max(3,K), "Dark2")
        for (j in 1:6) {
            if (j == 6) {
                print(paste0("$\\", param[j], "_", 1, "$"))
                print(quantile(m_props[,1]), probs = c(.025, .05, .25, .5, .75, .95, .975))
                plot(iter, m_props[, 1], type = "l",
                     main = TeX(paste0("Chain for $\\", param[j], "$")),
                     xlab = "t", ylab = TeX(paste0("$\\", param[j], "$")))
                abline(h = mix_props[1], col = "red", lty = "dashed")
            } else {
                print(paste0("$\\", param[j], "$"))
                print(quantile(mcmc_pars[,j]), probs = c(.025, .05, .25, .5, .75, .95, .975))
                plot(iter, mcmc_pars[,j], type = "l",
                     main = TeX(paste0("Chain for $\\", param[j], "$")),
                     xlab = "t", ylab = TeX(paste0("$\\", param[j], "$")))
                if (j == 4) {
                    mu_hat <- as.numeric(mean.circular(circular(mcmc_pars[,j])))
                    abline(h = mu_hat, col = "red", lty = "dashed")
                } else {
                    abline(h = post_means[,j], col = "red", lty = "dashed")
                }
            }
        }
    } else {
        m_props <- matrix(mcmc_pars[,,6], ncol = K)
        mix_props <- apply(m_props, 2, mean)
        # print(mix_props)
        iter <- 1:nrow(m_props)
        
        param <- c("alpha", "beta", "kappa", "mu", "lambda", "tau")
        post_means <- matrix(numeric(5*K), ncol = 5)
        for (i in 1:5) {
            post_means[,i] <- apply(matrix(mcmc_pars[,,i], ncol = K), 2, mean)
        }
        # print(post_means)
        par(mfrow = c(2,3))
        cols <- brewer.pal(max(3,K), "Dark2")
        for (j in 1:6) {
            if (j == 6) {
                print(paste0("$\\", param[j], "_", 1, "$"))
                print(quantile(m_props[,1]), probs = c(.025, .05, .25, .5, .75, .95, .975))
                plot(iter, m_props[, 1], type = "l",
                     main = TeX(paste0("Chain for $\\", param[j], "$")),
                     xlab = "t", ylab = TeX(paste0("$\\", param[j], "$")))
                for (k in 2:K) {
                    print(paste0("$\\", param[j], "_", k, "$"))
                    print(quantile(m_props[,k]), probs = c(.025, .05, .25, .5, .75, .95, .975))
                    lines(iter, m_props[, k], col = cols[k])
                    abline(h = mix_props[k], col = "red", lty = "dashed")
                }
            } else {
                print(paste0("$\\", param[j], "$"))
                print(quantile(mcmc_pars[,1,j]), probs = c(.025, .05, .25, .5, .75, .95, .975))
                plot(iter, mcmc_pars[,1,j], type = "l",
                     main = TeX(paste0("Chain for $\\", param[j], "$")),
                     xlab = "t", ylab = TeX(paste0("$\\", param[j], "$")))
                if (j == 4) {
                    mu_hat <- as.numeric(mean.circular(circular(mcmc_pars[,1,j])))
                    abline(h = mu_hat, col = "red", lty = "dashed")
                } else {
                    abline(h = post_means[1,j], col = "red", lty = "dashed")
                }
                for (k in 2:K) {
                    print(paste0("$\\", param[j], "_", k, "$"))
                    print(quantile(mcmc_pars[,k,j]), probs = c(.025, .05, .25, .5, .75, .95, .975))
                    lines(iter, mcmc_pars[,k,j], col = cols[k])
                    if (j == 4) {
                        mu_hat <- as.numeric(mean.circular(circular(mcmc_pars[,k,j])))
                        abline(h = mu_hat, col = "red", lty = "dashed")
                    } else {
                        abline(h = post_means[k,j], col = "red", lty = "dashed")
                    }
                }
            }
        }
    }
    
    
    
    
    par(mfrow = c(1,1))
}

# Metropolis-Hastings Sampler for the Abe-Ley Mixture model
# priors are as specified in the Sadeghianpourihamami paper
# build_colnames <- function(K) {
#     pars <- c("alpha_", "beta_", "kappa_", "mu_", "lambda_")
#     ind <- numeric(5*K)
#     for (k in 1:K) {
#         ind[(5*k - 4):(5*k)] <- rep(k, 5)
#     }
#     return(paste0(pars, ind))
# }

sample_proposal_dist <- function(mu_prev, sd_prev) {
    
    alpha <- rtruncnorm(1, a = 0, mean = mu_prev[1], sd = sd_prev[1])
    beta <- rtruncnorm(1, a = 0, mean = mu_prev[2], sd = sd_prev[2])
    kappa <- rtruncnorm(1, a = 0, mean = mu_prev[3], sd = sd_prev[3])
    mu <- rwrappednormal(1, mu = circular(mu_prev[4]), sd = sd_prev[4])
    lambda <- rtruncnorm(1, a = -1, b = 1, mean = mu_prev[5], sd = sd_prev[5])
    return(c(alpha, beta, kappa, mu, lambda))
}

calc_prop_dist_means <- function(y, sd) {
    mu_star <- numeric(5)
    for (i in 1:5) {
        if (i <= 3) {
            mu_star[i] <- y[i] + sd[i]*(1 - pnorm(-y[i]/sd[i])^(-1) * dnorm(-y[i]/sd[i]))
        } else if (i == 4) {
            mu_star[i] <- y[i]
        } else {
            mu_star[i] <- y[i] - sd[i]*(dnorm((1 - y[i])/sd[i]) - dnorm((-1 - y[i])/sd[i])) / (pnorm((1 - y[i])/sd[i]) - pnorm((-1 - y[i])/sd[i]))
        }
    }
    return(mu_star)
}

accept_ratio <- function(dat, y_star, y_prev, sd_prev) {
    # parameter order: alpha, beta, kappa, mu, lambda
    t <- dat[,1]; x <- dat[,2]
    px_prev <- sum(log(dabeley(x = x, theta = t, alpha = y_prev[1], beta = y_prev[2], 
                           kappa = y_prev[3], mu = y_prev[4], lambda = y_prev[5])))
    A <- numeric(5)
    for (i in 1:5) {
        y_tmp <- y_prev
        y_tmp[i] <- y_star[i]
        px_prop <- sum(log(dabeley(x = x, theta = t, alpha = y_tmp[1], beta = y_tmp[2], 
                                   kappa = y_tmp[3], mu = y_tmp[4], lambda = y_tmp[5])))
        if (i <= 3) {
            p_ystar <- dgamma(y_star[i], shape = 0.001, scale = 1000)
            p_prev <- dgamma(y_prev[i], shape = 0.001, scale = 1000)
            J_prev <- dtruncnorm(y_prev[i], a = 0, mean = y_star[i], sd = sd_prev[i])
            J_star <- dtruncnorm(y_star[i], a = 0, mean = y_prev[i], sd = sd_prev[i])
        }
        else if (i == 4) {
            p_ystar <- dvonmises(y_star[i], mu = circular(0), kappa = 0.001)
            p_prev <- dvonmises(y_prev[i], mu = circular(0), kappa = 0.001)
            J_star <- dwrappednormal(y_star[i], mu = circular(y_prev[i]), sd = sd_prev[i])
            J_prev <- dwrappednormal(y_prev[i], mu = circular(y_star[i]), sd = sd_prev[i])
        }
        else {
            p_ystar <- dbeta_rescale(y_star[i], 1,1)
            p_prev <- dbeta_rescale(y_prev[i], 1,1)
            J_star <- dtruncnorm(y_star[i], a = -1, b = 1, mean = y_prev[i], sd = sd_prev[i])
            J_prev <- dtruncnorm(y_prev[i], a = -1, b = 1, mean = y_star[i], sd = sd_prev[i])
        }
        r <- exp(px_prop + log(p_ystar) + log(J_prev) - px_prev - log(p_prev) - log(J_star))
        A[i] <- min(c(1,r))
    }
    
    return(A)
}

sample_class_labels <- function(dat, nu_t, K) {
    n <- nrow(dat)
    probs <- numeric(K)
    l_sample <- numeric(n)
    for (i in 1:n) {
        for (k in 1:K) {
            tmp <- nu_t[k,]
            # print(tmp)
            # print(dat[i,])
            probs[k] <- dabeley(theta = dat[i,1], x = dat[i,2], alpha = tmp[1], beta = tmp[2], 
                                 kappa = tmp[3], mu = tmp[4], lambda = tmp[5])
        }
        # print(probs)
        l_sample[i] <- sample(1:K, size = 1, prob = probs, replace = T)
    }
    return(l_sample)
}

sample_tau <- function(l_t, K, alpha0 = 1) {
    alpha <- rep(alpha0, K)
    for (k in 1:K) {
        alpha[k] <- alpha[k] + sum(l_t == k)
    }
    tau_star <- rdirichlet(1, alpha)
    return(tau_star)
}

adjust_sd <- function(iter, sd_prev, R, accept_probs) {
    K <- length(sd_prev)
    sd_new <- sd_prev
    if ((iter %% 50) == 0){
        for (k in 1:5) {
            s <- min(c(0.01, sqrt(R / iter)))
            if (accept_probs[k] > 0.44) {
                sd_new[k] <- sd_prev[k] + s
            }
            else {
                sd_new[k] <- ifelse(sd_prev[k] - s > 0, sd_prev[k] - s, sd_prev[k])
            }
        }
    }
    
    return(sd_new)
}

MH_posterior_estimation <- function(dat, K, x.obs = NULL, control = list(Q = 2000, burnin = 1000, sd_init = 1)) {
    # data should be:
    # theta     x
    # t1        x1
    # t2        x2
    # ...       ...
    colnames(dat) <- c("theta", "x")
    start <- Sys.time()
    N <- nrow(dat); sd_init <- control$sd_init; Q <- control$Q; burnin <- control$burnin
    # MCMC Draws will be stored in a QxKx6 array
    # The Q dim is the rows for each iteration of MCMC
    # The K dim is the columns for each mixture component
    # The J dim = 6 is the depth for each parameter
    # J = 1: alpha
    # J = 2: beta
    # J = 3: kappa
    # J = 4: mu
    # J = 5: lambda
    # J = 6: tau (mixing props)
    chains <- array(data = NA, dim = c(Q, K, 6))
    l_chains <- matrix(numeric(Q*N), ncol = N)
    if (!is.null(x.obs)) {
        M <- length(x.obs) 
        theta_pred <- matrix(numeric(Q * M), ncol = M)
        colnames(theta_pred) <- paste0("theta_pred_", 1:M)
    } else {
        theta_pred <- NULL
    }
    
    # Initialization
    sd_prev <- matrix(rep(sd_init, 5*K), ncol = 5)
    nu_0 <- matrix(numeric(5*K), ncol = 5)
    for (k in 1:K) {
        # mu_prev <- calc_prop_dist_means(c(1,1,1, 0, .5), sd_prev[(5*k-4):(5*k)])
        # nu_0[(5*k-4):(5*k)] <- sample_proposal_dist(mu_prev = mu_prev, sd_prev = sd_prev[(5*k-4):(5*k)])
        nu_0[k,] <- sample_proposal_dist(mu_prev = c(1,1,1,0,0), sd_prev = sd_prev[k,])
    }
    
    tau_0 <- rep(1/K, K)
    l_t <- sample(1:K, size = N, prob = tau_0, replace = T)
    
    if (!is.null(x.obs)) {
        M <- length(x.obs)
        theta_new <- numeric(M)
        for (m in 1:M) {
            nu_tmp <- nu_0
            conc_pred <- (nu_tmp[,2] * x.obs[m])^nu_tmp[,1] * tanh(nu_tmp[,3])
            theta_new[m] <- rssvm_mixture(1, mix_prop = tau_0,
                                          mu = nu_tmp[,4], kappa = conc_pred, lambda = nu_tmp[,5])
        }
        theta_pred[1,] <- theta_new
        pred_dat <- cbind(theta_new, x.obs)
        colnames(pred_dat) <- c("theta", "x")
    }
    
    chains[1,,1:5] <- nu_0
    chains[1,,6] <- tau_0
    l_chains[1,] <- l_t
    a_probs <- matrix(rep(1, 5*K), nrow = K)

    for (i in 1:(Q - 1)) {
        iter_start <- Sys.time()
        
        # Sample parameters for each mixture
        nu_prev <- matrix(chains[i,,1:5], nrow = K)
        nu_t <- matrix(numeric(5*K), nrow = K)
        for (k in 1:K) {
            # comp_k_obs <- which(l_t == k)
            # dat_k <- dat[comp_k_obs,]
            if (!is.null(x.obs)) {
                dat_k <- rbind(dat, pred_dat)
            } else {
                dat_k <- dat
            }
            
            nu_prev_k <- nu_prev[k,]
            sd_prev[k,] <- adjust_sd(iter = i, R = 50,
                                    sd_prev = sd_prev[k,],
                                    accept_probs = a_probs[k,]/50)
            
            # mu_prev <- calc_prop_dist_means(y_prev, sd_prev[(5*k-4):(5*k)])
            mu_prev <- nu_prev_k
            nu_star_k <- sample_proposal_dist(mu_prev, sd_prev[k,])
            mu_star <- nu_star_k
            # mu_star <- calc_prop_dist_means(y_star, sd_prev[(5*k-4):(5*k)])
            A <- accept_ratio(dat_k, mu_star, mu_prev, sd_prev[k,])
            U <- runif(5)
            for (j in 1:5) {
                if (U[j] >= A[j]) {
                    nu_star_k[j] <- nu_prev_k[j]
                } else {
                    a_probs[k,j] <- a_probs[k,j] + 1
                    nu_star_k[j] <- nu_star_k[j]
                }
                
            }
            if ((i %% 50) == 0) {
                a_probs <- matrix(rep(0, 5*K), nrow = K)
            }
            nu_t[k,] <- nu_star_k
        }
        chains[i+1,,1:5] <- nu_t
        # Sample the class labels
        l_t <- sample_class_labels(dat, nu_t, K)
        l_chains[i+1,] <- l_t
        # Sample the class probabilities
        chains[i+1,,6] <- sample_tau(l_t, K, alpha0 = 1)
        
        # Sample the missing data 
        # Theta first
        if (!is.null(x.obs)) {
            M <- length(x.obs)
            theta_new <- numeric(M)
            for (m in 1:M) {
                nu_tmp <- matrix(nu_chains[i+1,], byrow = T, ncol = 5)
                conc_pred <- (nu_tmp[,2] * x.obs[m])^nu_tmp[,1] * tanh(nu_tmp[,3])
                theta_new[m] <- rssvm_mixture(1, mix_prop = tau_chains[i+1,],
                                                    mu = nu_tmp[,4], kappa = conc_pred, lambda = nu_tmp[,5])
            }
            theta_pred[i+1,] <- theta_new
            pred_dat <- cbind(theta_new, x.obs)
            colnames(pred_dat) <- c("theta", "x")
        }
        
        iter_end <- Sys.time()
        if (i %% 100 == 0)
            print(paste0("Iteration: ", i, "; Iteration time: ", round(iter_end - iter_start, 3), "; Total time: ", round(iter_end - start, 3)))
    }
    end <- Sys.time()
    print(paste0("Total time elapsed: ", round(end - start, 3)))
    return(list(mcmc_pars = chains, class_labels = l_chains, theta_pred = theta_pred))
}





