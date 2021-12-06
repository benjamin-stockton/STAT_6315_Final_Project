library(mvtnorm)
library(circular)
library(MASS)
library(latex2exp)
library(truncnorm)
library(MCMCpack)


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
        t <- ifelse(t1 < (1 + lambda * sin(t1 - mu))/2, t1, -t1) %% (2*pi)
        shape_par <- beta * (1 - tanh(kappa) * cos(t - mu))^(-1/alpha)
        x <- rweibull(1, alpha, shape = shape_par)
        dat[i,] <- c(t,x)
    }
    return(dat)
}

dbeta_rescale <- function(x, shape1 = 1, shape2 = 1) {
    return(dbeta((x/2 + 1/2), shape1 = shape1, shape2 = shape2)/2)
}

cyl_density_plot <- function(dat, main = "Joint Distibution on the Cylinder", xlab = "$\\theta$", ylab = "X") {
    z <- kde2d(dat[,1], dat[,2], n = 50)
    
    plot(dat[,1], dat[,2],
         pch = 16, main = TeX(main), xlab = TeX(xlab), ylab = TeX(ylab))
    contour(z, lwd = 2, add = TRUE, nlevels = 10,
            col = hcl.colors(10, "Spectral"))
}

pred_density_plot <- function(dat, post_means, mix_prop, main = "Fitted Distibution on the Cylinder", xlab = "$\\theta$", ylab = "X") {
    fit_dat <- rabeley_mixture(1500, mix_prop = mix_prop, mu = post_means[,4], kappa = post_means[,3],
                               alpha = post_means[,1], beta = post_means[,2], lambda = post_means[,5])
    
    z <- kde2d(fit_dat[,1], fit_dat[,2], n = 50)
    
    plot(dat[,1], dat[,2],
         pch = 16, main = TeX(main), xlab = TeX(xlab), ylab = TeX(ylab),
         xlim = c(min(z$x), max(z$x)), ylim = c(min(z$y), max(z$y)))
    points(fit_dat[,1], fit_dat[,2], pch = 1, col = "lightgrey")
    contour(z, lwd = 2, add = TRUE, nlevels = 10,
            col = hcl.colors(10, "Spectral"))
    
}

hist_density <- function(dat, main = "Distribution of X", xlab = "", ylab = "Density",...) {
    dx <- density(dat)
    hx <- hist(dat, plot = F)
    plot(hx, col="grey", main = TeX(main), xlab = TeX(xlab), ylab = TeX(ylab),...)
    lines(x = dx$x, y = dx$y * length(dat) * diff(hx$breaks)[1], lwd = 2)
}

joint_dist_plot <- function(dat, main = "", xlab = "$\\theta$", ylab = "X", lhist = 20) {
    layout(matrix(c(2,0,1,3), nrow = 2, byrow = T),
           widths = c(3,1), heights = c(1,3), respect = T)
    
    par. <- par(mar = c(4, 4, 1,1), oma = rep(.5, 4))
    
    cyl_density_plot(dat = dat, main = main, xlab = xlab, ylab = ylab)
    
    t_hist <- hist(dat[,1], plot = F, breaks = seq(from = min(dat[,1]), to = max(dat[,1]),
                                                   length.out = lhist))
    x_hist <- hist(dat[,2], plot = F, breaks = seq(from = min(dat[,2]), to = max(dat[,2]),
                                                   length.out = lhist))
    
    td <- density(dat[,1])
    xd <- density(dat[,2])
    
    par(mar = c(0, 4, 0,0))
    barplot(t_hist$density, axes = F,
            ylim = c(0, max(t_hist$density)), space = 0)
    lines(x = td$x, y = td$y, lwd = 2)
    # * length(dat[,1]) * diff(t_hist$breaks)[1]
    
    par(mar = c(4,0,0,0))
    barplot(x_hist$density, axes = F,
            xlim = c(0, max(x_hist$density)),
            space = 0, horiz = T)
    lines(x = xd$x, y = xd$y, lwd = 2)
    # * length(dat[,2]) * diff(x_hist$breaks)[1]
    
    par(par.)
}


# Mixture Abe-Ley Simulated Data
# dabeley_mixture <- function(x,theta, mix_prop, mu, kappa, lambda, alpha, beta) {
#     
# }
rabeley_mixture <- function(N, mix_prop, mu, kappa, lambda, alpha, beta) {
    dat <- matrix(numeric(2*N), ncol = 2)
    colnames(dat) <- c("theta", "x")
    for (i in 1:N) {
        k <- which(rmultinom(1, 1, mix_prop) == 1)
        tmp <- rabeley(1, mu[k], kappa[k], lambda[k], alpha[k], beta[k])
        dat[i,] <- tmp
    }
    return(dat)
}


# priors are as specified in the Sadeghianpourihamami paper
build_colnames <- function(K) {
    pars <- c("alpha_", "beta_", "kappa_", "mu_", "lambda_")
    ind <- numeric(5*K)
    for (k in 1:K) {
        ind[(5*k - 4):(5*k)] <- rep(k, 5)
    }
    return(paste0(pars, ind))
}

sample_proposal_dist <- function(mu_prev, sd_prev) {
    alpha <- rtruncnorm(1, a = 0, mean = mu_prev[1], sd = sd_prev[1])
    beta <- rtruncnorm(1, a = 0, mean = mu_prev[2], sd = sd_prev[2])
    kappa <- rtruncnorm(1, a = 0, mean = mu_prev[3], sd = sd_prev[3])
    mu <- rwrappednormal(1, mu = circular(mu_prev[4]), sd = sd_prev[4])
    lambda <- rtruncnorm(1, a = -1, b = 1, mean = mu_prev[5], sd = sd_prev[5])
    return(c(alpha, beta, kappa, mu, lambda))
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
        # print(paste0("px_prev: ", round(px_prev, 4), ", px_prop: ", round(px_prop, 4), ", p_prev: ", round(p_prev, 4), ", p_ystar: ", round(p_ystar, 4), ", J_prev: ", round(J_prev, 4), ", J_star: ", round(J_star, 4)))
        r <- exp(px_prop + log(p_ystar) + log(J_prev) - px_prev - log(p_prev) - log(J_star))
        # print(paste0("r = ", round(r, 4)))
        # if (!is.nan(r)) {
        #     A[i] <- min(c(1, r))
        # } else {
        #     A[i] <- -999
        # }
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
            probs[k] <- dabeley(dat[,1], dat[,2], alpha = nu_t[1], beta = nu_t[2], 
                                 kappa = nu_t[3], mu = nu_t[4], lambda = nu_t[5])
        }
        l_sample[i] <- sample(1:K, size = n, prob = probs, replace = T)
    }
    return(l_sample)
}

sample_tau <- function(l_t, K, alpha0 = 1) {

    alpha <- rep(alpha0, K)
    for (k in 1:K) {
        alpha[k] <- alpha[k] + sum(l_t == k)
    }
    tau_star <- rdirichlet(1, alpha)
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

MH_posterior_estimation <- function(dat, K, control = list(Q = 2000, burnin = 1000, sd_init = 1)) {
    start <- Sys.time()
    N <- nrow(dat); sd_init <- control$sd_init; Q <- control$Q; burnin <- control$burnin
    # Each row of the chains is of the form:
    # (rep(alpha_k, beta_k, kappa_k, mu_k, lambda_k, K), tau_1,..., tau_K)
    nu_chains <- matrix(numeric(Q * K*5), nrow = Q)
    colnames(nu_chains) <- build_colnames(K)
    tau_chains <- matrix(numeric(Q * K), nrow = Q)
    colnames(tau_chains) <- paste0("tau_", 1:K)
    l_chains <- matrix(numeric(Q*N), ncol = Q)
    # Initialization
    sd_prev <- rep(sd_init, 5*K)
    nu_0 <- numeric(5*K)
    for (k in 1:K) {
        nu_0[(5*k-4):(5*k)] <- sample_proposal_dist(mu_prev = c(1,1,1, 0, .5), sd_prev = sd_prev[(5*k-4):(5*k)])
    }
    
    tau_0 <- rep(1/K, K)
    l_t <- sample(1:K, size = N, prob = tau_0, replace = T)
    
    nu_chains[1,] <- nu_0
    tau_chains[1,] <- tau_0
    l_chains[,1] <- l_t
    a_probs <- rep(1, 5*K)

    for (i in 1:(Q - 1)) {
        iter_start <- Sys.time()
        
        nu_prev <- nu_chains[i,]
        nu_t <- numeric(5*K)
        for (k in 1:K) {
            # comp_k_obs <- which(l_t == k)
            # dat_k <- dat[comp_k_obs,]
            dat_k <- dat
            # print(paste0("Mixture Component: ", k))
            y_prev <- nu_prev[(5*k-4):(5*k)]
            # print("y_prev: ")
            # print(y_prev, digits = 4)
            sd_prev[(5*k-4):(5*k)] <- adjust_sd(iter = i, R = 50,
                                                sd_prev = sd_prev[(5*k-4):(5*k)],
                                                accept_probs = a_probs[(5*k-4):(5*k)]/50)
            # print("sd_prev: ")
            # print(sd_prev, digits = 4)
            y_star <- sample_proposal_dist(y_prev, sd_prev)
            # print("y_star: ")
            # print(y_star, digits = 4)
            A <- accept_ratio(dat_k, y_star, y_prev, sd_prev)
            # print("Acceptance Ratios: ")
            # print(A, digits = 4)
            U <- runif(5)
            for (j in 1:5) {
                if (U[j] >= A[j]) {
                    y_star[j] <- y_prev[j]
                } else {
                    a_probs[(5*k-4):(5*k)][j] <- a_probs[(5*k-4):(5*k)][j] + 1
                    y_star[j] <- y_star[j]
                }
                
            }
            if ((i %% 50) == 0) {
                a_probs <- rep(0, 5*K)
            }
            # print("y_star: ")
            # print(y_star, digits = 4)
            nu_t[(5*k-4):(5*k)] <- y_star
        }
        # print(i)
        # print(paste0("nu_", i, ": "))
        # print(nu_t, digits = 4)
        nu_chains[i+1,] <- nu_t
        l_t <- sample_class_labels(dat, nu_t, K)
        l_chains[,i+1] <- l_t
        tau_chains[i+1,] <- sample_tau(l_t, K, alpha0 = 1)
        iter_end <- Sys.time()
        if (i %% 100 == 0)
            print(paste0("Iteration: ", i, "; Iteration time: ", iter_end - iter_start, "; Total time: ", iter_end - start))
    }
    end <- Sys.time()
    print(paste0("Total time elapsed: ", end - start))
    return(list(dist_pars = nu_chains, mix_props = tau_chains, class_labels = l_chains))
}





