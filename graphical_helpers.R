# graphical_helpers.R
library(mvtnorm)
library(circular)
library(MASS)
library(latex2exp)
library(truncnorm)
library(MCMCpack)
library(RColorBrewer)
source("abe_ley_mixture_metropolis_hastings.R")

# Plotting Functions

cyl_density_plot <- function(dat, mean_pars, mix_prop, main = "Joint Distibution on the Cylinder", xlab = "$\\theta$", ylab = "X", nlevels = 10) {
    x_seq <- seq(from = min(dat[,2]), to = max(dat[,2]), length.out = 50)
    t_seq <-  seq(from = min(dat[,1]), to = max(dat[,1]), length.out = 50)
    
    z <- t(outer(X = x_seq, Y = t_seq, FUN = dabeley_mixture,
                 mix_prop = mix_prop,
                 alpha = mean_pars[,1], beta = mean_pars[,2],
                 kappa = mean_pars[,3], mu = mean_pars[,4], lambda = mean_pars[,5]))
    
    plot(dat[,1], dat[,2], pch = 16, col = "black",
         main = TeX(main), xlab = TeX(xlab), ylab = TeX(ylab), xlim = c(0, 2*pi))
    contour(x = t_seq, y = x_seq, z = z, lwd = 2, add = TRUE, nlevels = nlevels,
            col = hcl.colors(nlevels, "Spectral"))
}

hist_density <- function(dat, main = "Distribution of X", xlab = "", ylab = "Density", lhist = 10,...) {
    dx <- density(dat)
    hx <- hist(dat, plot = F, breaks = seq(from = min(dat), to = max(dat),
                                           length.out = lhist))
    plot(hx, col="grey", main = TeX(main), xlab = TeX(xlab), ylab = TeX(ylab),...)
    lines(x = dx$x, y = dx$y * length(dat) * diff(hx$breaks)[1], lwd = 2)
}

joint_dist_plot <- function(dat, mean_pars, mix_prop, main = "", xlab = "$\\theta$", ylab = "X", lhist = 20, nlevels = 10) {
    layout(matrix(c(2,0,1,3), nrow = 2, byrow = T),
           widths = c(3,1), heights = c(1,3), respect = T)
    
    par. <- par(mar = c(4, 4, 1,1), oma = rep(.5, 4))
    
    cyl_density_plot(dat = dat, mean_pars = mean_pars, mix_prop = mix_prop, main = main, xlab = xlab, ylab = ylab, nlevels = nlevels)
    
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
        iter <- 1:(dim(mcmc_pars)[1])
        
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
                print(quantile(mcmc_pars[,1,j]), probs = c(.025, .05, .25, .5, .75, .95, .975))
                plot(iter, mcmc_pars[,1,j], type = "l",
                     main = TeX(paste0("Chain for $\\", param[j], "$")),
                     xlab = "t", ylab = TeX(paste0("$\\", param[j], "$")), 
                     ylim = c(min(m_props), max(m_props)))
                for (k in 2:K) {
                    print(paste0("$\\", param[j], "_", k, "$"))
                    print(quantile(m_props[,k]), probs = c(.025, .05, .25, .5, .75, .95, .975))
                    lines(iter, mcmc_pars[,k,j], col = cols[k])
                    abline(h = mix_props[k], col = "red", lty = "dashed")
                }
            } else {
                print(paste0("$\\", param[j], "$"))
                print(quantile(mcmc_pars[,1,j]), probs = c(.025, .05, .25, .5, .75, .95, .975))
                plot(iter, mcmc_pars[,1,j], type = "l",
                     main = TeX(paste0("Chain for $\\", param[j], "$")),
                     xlab = "t", ylab = TeX(paste0("$\\", param[j], "$")), 
                     ylim = c(min(mcmc_pars[,,j]), max(mcmc_pars[,,j])))
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

plot_post_density <- function(mcmc_pars, K, ctrl) {
    pars <- c("alpha", "beta", "kappa", "mu", "lambda", "tau")
    par(mfrow = c(2,3))
    cols <- brewer.pal(K, "Dark2")
    
    for (j in 1:6) {
        d <- density(mcmc_pars[,1,j])
        plot(d, main = TeX(paste0("Posterior Density of $\\", pars[j], "$")),
             xlab = TeX(paste0("$\\", pars[j], "$")), lwd = 2, col = cols[1], 
             xlim = c(min(mcmc_pars[,,j]), max(mcmc_pars[,,j])))
        # polygon(d, lwd = 2, border = cols[j])
        if (K >= 2) {
            for (k in 2:K) {
                d <- density(mcmc_pars[,k,j])
                lines(d$x, d$y, lwd = 2, col = cols[k])
            }
        }
        
    }
    par(mfrow = c(1,1))
}
