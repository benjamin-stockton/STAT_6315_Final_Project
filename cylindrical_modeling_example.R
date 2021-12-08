library(mvtnorm)
library(circular)
library(MASS)
library(latex2exp)
library(truncnorm)
library(MCMCpack)
source("abe_ley_mixture_metropolis_hastings.R")

# Single component Abe-Ley example from original Abe-Ley paper
true_par <- cbind( c(2), c(1), c(1.5), c(0), c(0))
true_mix <- c(1)
# true_par <- cbind(c(1,5,2), c(1,2,3), c(2,1,3), c(0, 2*pi/3, 4*pi/3), c(-.5, 0, .5))
# true_mix <- c(.33, .33, .33)
colnames(true_par) <- c("alpha", "beta", "kappa", "mu", "lambda")
sample_dat2 <- rabeley_mixture(150, mix_prop = true_mix,
                               mu = true_par[,4], kappa = true_par[,3], lambda = true_par[,5],
                               alpha = true_par[,1], beta = true_par[,2])

joint_dist_plot(sample_dat2, mean_pars = true_par, mix_prop = true_mix, main = "True Abe-Ley Mixture Density")


ctrl <- list(Q = 2500, burnin = 1500, sd_init = 1)
K <- 1
fit <- MH_posterior_estimation(sample_dat2, K = K, control = ctrl)

iters <- (ctrl$burnin+1):ctrl$Q
mix_props <- apply(fit$mix_props[iters,], 2, mean)
mix_props
post_means <- matrix(apply(fit$dist_pars[iters,], 2, mean), ncol = 5, byrow = T)
colnames(post_means) <- c("alpha", "beta", "kappa", "mu", "lambda")
post_means

plot_tracestack(fit, ctrl)

joint_dist_plot(sample_dat2, mean_pars = true_par, mix_prop = true_mix, main = "True Abe-Ley Mixture Density")
joint_dist_plot(sample_dat2, mean_pars = post_means, mix_prop = 1, main = "Fitted Abe-Ley Mixture Density")



