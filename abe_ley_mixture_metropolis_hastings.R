library(mvtnorm)
library(circular)
library(MASS)
library(latex2exp)
library(truncnorm)
library(MCMCpack)
library(RColorBrewer)

##########################################################
# Density Functions
##########################################################

dabeley <- function(x, theta, mu, kappa, lambda, alpha, beta) {
    const <- alpha * beta^alpha / (2*pi * cosh(kappa))
    wc <- (1 + lambda * sin(theta - mu))
    weib <- x^(alpha - 1) * exp(-(beta * x)^alpha * (1 - tanh(kappa) * cos(theta - mu)))
    return(const * wc * weib)
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

##########################################################
## RNG Functions
##########################################################

rabeley <- function(n, mu, kappa, lambda, alpha, beta) {
    dat <- matrix(numeric(2*n), ncol = 2)
    colnames(dat) <- c("theta", "x")
    for (i in 1:n) {
        t1 <- rwrappedcauchy(1, mu = circular(0), tanh(kappa/2))
        U <- runif(1)
        t <- ifelse(U < (1 + lambda * sin(t1)) / 2, t1, t1) + mu
        scale_par <- beta * (1 - tanh(kappa) * cos(mu-t))^(1/alpha)
        x <- rweibull(1, shape = alpha, scale = 1/scale_par)
        dat[i,] <- c(t %% (2*pi),x)
    }
    return(dat)
}

rssvm <- function(n, mu, kappa, lambda) {
    t1 <- rvonmises(n, circular(0), kappa)
    t <- sapply(t1, FUN = function(t) {return(ifelse(runif(1) < (1 + lambda * sin(t))/2, t, -t))}) + mu
    return(t %% (2*pi))
}

rsswc <- function(n, mu, kappa, lambda) {
    t1 <- rwrappedcauchy(n, circular(0), tanh(kappa/2))
    t <- sapply(t1, FUN = function(t) {return(ifelse(runif(1) < (1 + lambda * sin(t))/2, t, -t))}) + mu
    return(t %% (2*pi))
}

rbeta_rescale <- function(n, shape1 = 1, shape2 = 1) {
    U <- rbeta(n, shape1 = shape1, shape2 = shape2)
    V <- 2 * U - 1
    return(V)
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

############################################################
############################################################
# Metropolis-Hastings Sampler for the Abe-Ley Mixture model
# priors are as specified in the Sadeghianpourihamami paper
############################################################
############################################################

sample_proposal_dist <- function(mu_prev, sd_prev) {
    alpha <- rtruncnorm(1, a = 0, mean = mu_prev[1], sd = sd_prev[1])
    beta <- rtruncnorm(1, a = 0, mean = mu_prev[2], sd = sd_prev[2])
    kappa <- rtruncnorm(1, a = 0, mean = mu_prev[3], sd = sd_prev[3])
    mu <- rwrappednormal(1, mu = circular(mu_prev[4]), sd = sd_prev[4])
    lambda <- rtruncnorm(1, a = -1, b = 1, mean = mu_prev[5], sd = sd_prev[5])
    return(c(alpha, beta, kappa, mu, lambda))
}

# calc_prop_dist_means <- function(y, sd) {
#     mu_star <- numeric(5)
#     for (i in 1:5) {
#         if (i <= 3) {
#             mu_star[i] <- y[i] + sd[i]*(1 - pnorm(-y[i]/sd[i])^(-1) * dnorm(-y[i]/sd[i]))
#         } else if (i == 4) {
#             mu_star[i] <- y[i]
#         } else {
#             mu_star[i] <- y[i] - sd[i]*(dnorm((1 - y[i])/sd[i]) - dnorm((-1 - y[i])/sd[i])) / (pnorm((1 - y[i])/sd[i]) - pnorm((-1 - y[i])/sd[i]))
#         }
#     }
#     return(mu_star)
# }

accept_ratio <- function(dat, y_star, y_prev, sd_prev, hyperpar = list(alpha_shape = 1000, alpha_scale = 0.001, beta_shape = 1000, beta_scale = 0.001, kappa_shape = 1000, kappa_scale = 0.001, mu_mean = 0, mu_conc = 0.001, lambda_a = 1, lambda_b = 1, alpha0 = 1)) {
    # print(str(dat))
    # parameter order: alpha, beta, kappa, mu, lambda
    t <- dat[,1]; x <- dat[,2]
    n <- length(t)
    px_prev <- sum(log(dabeley(x = x, theta = t, alpha = y_prev[1], beta = y_prev[2], 
                           kappa = y_prev[3], mu = y_prev[4], lambda = y_prev[5])))
    A <- numeric(5)
    for (i in 1:5) {
        y_tmp <- y_prev
        y_tmp[i] <- y_star[i]
        px_prop <- sum(log(dabeley(x = x, theta = t, alpha = y_tmp[1], beta = y_tmp[2], 
                                   kappa = y_tmp[3], mu = y_tmp[4], lambda = y_tmp[5])))
        if (i == 1) {
            p_ystar <- dgamma(y_star[i], shape = hyperpar$alpha_shape, scale = hyperpar$alpha_scale)
            p_prev <- dgamma(y_prev[i], shape = hyperpar$alpha_shape, scale = hyperpar$alpha_scale)
            J_prev <- dtruncnorm(y_prev[i], a = 0, mean = y_star[i], sd = sd_prev[i])
            J_star <- dtruncnorm(y_star[i], a = 0, mean = y_prev[i], sd = sd_prev[i])
        }
        else if (i == 2) {
            p_ystar <- dgamma(y_star[i], shape = hyperpar$beta_shape, scale = hyperpar$beta_scale)
            p_prev <- dgamma(y_prev[i], shape = hyperpar$beta_shape, scale = hyperpar$beta_scale)
            J_prev <- dtruncnorm(y_prev[i], a = 0, mean = y_star[i], sd = sd_prev[i])
            J_star <- dtruncnorm(y_star[i], a = 0, mean = y_prev[i], sd = sd_prev[i])
        }
        else if (i == 3) {
            p_ystar <- dgamma(y_star[i], shape = hyperpar$kappa_shape, scale = hyperpar$kappa_scale)
            p_prev <- dgamma(y_prev[i], shape = hyperpar$kappa_shape, scale = hyperpar$kappa_scale)
            J_prev <- dtruncnorm(y_prev[i], a = 0, mean = y_star[i], sd = sd_prev[i])
            J_star <- dtruncnorm(y_star[i], a = 0, mean = y_prev[i], sd = sd_prev[i])
        }
        else if (i == 4) {
            p_ystar <- dvonmises(circular(y_star[i]), mu = circular(hyperpar$mu_mean), kappa = hyperpar$mu_conc)
            p_prev <- dvonmises(circular(y_prev[i]), mu = circular(hyperpar$mu_mean), kappa = hyperpar$mu_conc)
            J_star <- dwrappednormal(circular(y_star[i]), mu = circular(y_prev[i]), sd = sd_prev[i])
            J_prev <- dwrappednormal(circular(y_prev[i]), mu = circular(y_star[i]), sd = sd_prev[i])
        }
        else {
            p_ystar <- dbeta_rescale(y_star[i], hyperpar$lambda_a, hyperpar$lambda_b)
            p_prev <- dbeta_rescale(y_prev[i], hyperpar$lambda_a, hyperpar$lambda_b)
            J_star <- dtruncnorm(y_star[i], a = -1, b = 1, mean = y_prev[i], sd = sd_prev[i])
            J_prev <- dtruncnorm(y_prev[i], a = -1, b = 1, mean = y_star[i], sd = sd_prev[i])
        }
        # print(paste0("y_prev = ", y_prev[i]))
        # print(paste0("y_star = ", y_star[i]))
        # print(paste0("p = ", round(c(px_prev, px_prop, p_prev, p_ystar, J_prev, J_star), 5)))
        r <- px_prop + log(p_ystar) + log(J_prev) - px_prev - log(p_prev) - log(J_star)
        A[i] <- min(c(0,r))
        # print(paste0("A = ", A[i]))
    }
    
    return(A)
}

sample_class_labels <- function(dat, nu_t, K, tau) {
    n <- nrow(dat)
    probs <- numeric(K)
    l_sample <- numeric(n)
    for (i in 1:n) {
        for (k in 1:K) {
            tmp <- nu_t[k,]
            # print(tmp)
            # print(dat[i,])
            # print(tau[k])
            probs[k] <- tau[k] * dabeley(theta = dat[i,1], x = dat[i,2], alpha = tmp[1], beta = tmp[2], 
                                 kappa = tmp[3], mu = tmp[4], lambda = tmp[5])
        }
        # print(max(probs/sum(probs)))
        probs <- probs / sum(probs)
        # print("Class probs:")
        # print(probs)
        l_sample[i] <- sample(1:K, size = 1, prob = probs, replace = T)
    }
    return(l_sample)
}

sample_tau <- function(l_t, K, alpha0 = 1) {
    alpha <- rep(alpha0, K)
    for (k in 1:K) {
        alpha[k] <- alpha0 + sum(l_t == k)
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
                sd_new[k] <- sqrt(sd_prev[k]^2 + s)
            }
            else {
                sd_new[k] <- ifelse(sqrt(sd_prev[k]^2 - s) > 0, sqrt(sd_prev[k]^2 - s), sd_prev[k])
            }
        }
    }
    
    return(sd_new)
}

MH_posterior_estimation <- function(dat, K, x.obs = NULL, control = list(Q = 2000, burnin = 1000, sd_init = 1, thin = 1), hyperpar = list(alpha_shape = 1000, alpha_scale = 0.001, beta_shape = 1000, beta_scale = 0.001, kappa_shape = 1000, kappa_scale = 0.001, mu_mean = 0, mu_conc = 0.001, lambda_a = 1, lambda_b = 1, alpha0 = 1)) {
    # data should be:
    # theta     x
    # t1        x1
    # t2        x2
    # ...       ...
    colnames(dat) <- c("theta", "x")
    start <- Sys.time()
    N <- nrow(dat); sd_init <- control$sd_init; Q <- control$Q; burnin <- control$burnin
    # MCMC Draws will be stored in a QxKx6 array
    # The Q dim is the rows; for each iteration of MCMC
    # The K dim is the columns; for each mixture component
    # The J dim = 6 is the depth; for each parameter
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
        nu_0[k,] <- c(rgamma(1, shape = hyperpar$alpha_shape, scale = hyperpar$alpha_scale),
                      rgamma(1, shape = hyperpar$beta_shape, scale = hyperpar$beta_scale),
                      rgamma(1, shape = hyperpar$kappa_shape, scale = hyperpar$kappa_scale),
                      rvonmises(1, circular(hyperpar$mu_mean), kappa = hyperpar$mu_conc),
                      rbeta_rescale(1, shape1 = hyperpar$lambda_a, shape2 = hyperpar$lambda_b))
        # nu_0[k,] <- sample_proposal_dist(c(1,1,1,0,0), rep(sd_init, 5))
    }
    # print(nu_0)
    tau_0 <- rdirichlet(1, rep(hyperpar$alpha0, K))
    l_t <- sample(1:K, size = N, prob = tau_0, replace = T)
    
    if (!is.null(x.obs)) {
        M <- length(x.obs)
        theta_new <- numeric(M)
        for (m in 1:M) {
            conc_pred <- (nu_0[,2] * x.obs[m])^nu_0[,1] * tanh(nu_0[,3])
            theta_new[m] <- rssvm_mixture(1, mix_prop = tau_0,
                                          mu = nu_0[,4], kappa = conc_pred, lambda = nu_0[,5])
        }
        theta_pred[1,] <- theta_new
        pred_dat <- cbind(theta_new, x.obs)
        colnames(pred_dat) <- c("theta", "x")
    }
    
    chains[1,,1:5] <- nu_0
    chains[1,,6] <- tau_0
    l_chains[1,] <- l_t
    a_probs <- matrix(rep(1, 5*K), nrow = K)

    for (q in 1:(Q - 1)) {
        iter_start <- Sys.time()
        
        # Sample parameters for each mixture
        nu_prev <- matrix(chains[q,,1:5], nrow = K)
        nu_t <- matrix(numeric(5*K), nrow = K)
        for (k in 1:K) {
            ind.k <- which(l_chains[q,] == k)
            dat_k <- matrix(dat[ind.k,], ncol = 2)
            
            if (!is.null(x.obs)) {
                dat_k <- rbind(dat_k, pred_dat)
            }
            
            nu_prev_k <- nu_prev[k,]
            sd_prev[k,] <- adjust_sd(iter = q, R = 50,
                                    sd_prev = sd_prev[k,],
                                    accept_probs = a_probs[k,]/50)
            
            nu_star_k <- sample_proposal_dist(nu_prev_k, sd_prev[k,])
            
            A <- accept_ratio(dat_k, nu_star_k, nu_prev_k, sd_prev[k,], hyperpar = hyperpar)
            U <- runif(5)
            for (j in 1:5) {
                if (log(U[j]) >= A[j]) {
                    nu_star_k[j] <- nu_prev_k[j]
                } else {
                    a_probs[k,j] <- a_probs[k,j] + 1
                    nu_star_k[j] <- nu_star_k[j]
                }
            }
            if ((q %% 50) == 0) {
                a_probs <- matrix(rep(0, 5*K), nrow = K)
            }
            nu_t[k,] <- nu_star_k
        }
        chains[q+1,,1:5] <- nu_t
        
        # Sample the class labels
        l_t <- sample_class_labels(dat, nu_t, K, tau = chains[q, , 6])
        l_chains[q+1,] <- l_t
        
        # Sample the class probabilities
        chains[q+1,,6] <- sample_tau(l_t, K, alpha0 = hyperpar$alpha0)
        
        # Sample the missing data 
        # Theta first
        if (!is.null(x.obs)) {
            M <- length(x.obs)
            theta_new <- numeric(M)
            for (m in 1:M) {
                nu_tmp <- nu_t
                tau_tmp <- chains[q+1, , 6]
                conc_pred <- (nu_tmp[,2] * x.obs[m])^nu_tmp[,1] * tanh(nu_tmp[,3])
                theta_new[m] <- rssvm_mixture(1, mix_prop = tau_tmp,
                                                    mu = nu_tmp[,4], kappa = conc_pred, lambda = nu_tmp[,5])
            }
            theta_pred[q+1,] <- theta_new
            pred_dat <- cbind(theta_new, x.obs)
            colnames(pred_dat) <- c("theta", "x")
        }
        
        iter_end <- Sys.time()
        if (q %% 100 == 0)
            print(paste0("Iteration ", q, "; Iteration time: ", round(iter_end - iter_start, 3), "; Total time: ", round(iter_end - start, 3)))
    }
    end <- Sys.time()
    print(paste0("Total time elapsed: ", round(end - start, 3)))
    
    # thin_index <- 1:floor(Q / thin) * ctrl$thin
    # return(list(mcmc_pars = chains[thin_index,,], class_labels = l_chains[thin_index,], theta_pred = theta_pred[thin_index,]))
    return(list(mcmc_pars = chains, class_labels = l_chains, theta_pred = theta_pred))
}





