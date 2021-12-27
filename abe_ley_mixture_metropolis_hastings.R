library(mvtnorm)
library(circular)
library(MASS)
library(latex2exp)
library(truncnorm)
library(MCMCpack)
library(RColorBrewer)


# Abe-Ley Distribution random variable draw

dabeley <- function(x, theta, mu, kappa, lambda, alpha, beta) {
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
        U <- runif(1)
        t <- ifelse(U < (1 + lambda * sin(t1 - mu)) / 2, t1, -t1)
        scale_par <- beta * (1 - tanh(kappa) * cos(t - mu))^(1/alpha)
        x <- rweibull(1, shape = alpha, scale = 1/scale_par)
        dat[i,] <- c(t %% (2*pi),x)
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
        tmp <- rabeley(1, mu = mu[k], kappa = kappa[k], lambda = lambda[k], alpha = alpha[k], beta = beta[k])
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
        nu_0[k,] <- sample_proposal_dist(mu_prev = c(1,1,1,(k-1)*pi/3,0), sd_prev = sd_prev[k,])
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
                nu_tmp <- nu_t
                tau_tmp <- chains[i+1, , 6]
                conc_pred <- (nu_tmp[,2] * x.obs[m])^nu_tmp[,1] * tanh(nu_tmp[,3])
                theta_new[m] <- rssvm_mixture(1, mix_prop = tau_tmp,
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





