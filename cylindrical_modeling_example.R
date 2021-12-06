library(mvtnorm)
library(circular)
library(MASS)
library(latex2exp)
library(truncnorm)
library(MCMCpack)

source("abe_ley_mixture_metropolis_hastings.R")

t <- rwrappedcauchy(100)
plot(t, pch = 16, shrink = 1.5)
lines(density.circular(t, bw = 5))

true_params <- c(10, 5, 1/5, 0, 1/2)
sample_dat <- rabeley(150, mu = 0, kappa = 1, lambda = 0, alpha = 10, beta = 5)

joint_dist_plot(sample_dat, xlab = "$\\theta$", ylab = "X")

sample_dat2 <- rabeley_mixture(150, mix_prop = c(.5, .05, .45),
                               mu = c(0, pi/2, 5*pi/4), kappa = c(2,15,5), lambda = c(-.5, 0, .5),
                               alpha = c(1,1,10), beta = c(1,15,5))

joint_dist_plot(sample_dat2, xlab = "$\\theta$", ylab = "X")

ctrl <- list(Q = 2500, burnin = 1250, sd_init = 1)
K <- 1
fit <- MH_posterior_estimation(sample_dat2, K = K, control = ctrl)

m_props <- fit$mix_props
mix_props <- apply(m_props, 2, mean)
iter <- (ctrl$burnin+1):ctrl$Q

par(mfrow = c(K, 1))
for (k in 1:K) {
    plot(iter, m_props[iter, k], type = "l",
         main = paste0("Chain for Mixing Prop. ", k),
         xlab = "t", ylab = paste0("Prop ", k))
}
mix_props

param_post <- fit$dist_pars
param <- c("alpha", "beta", "kappa", "mu", "lambda")
post_means <- matrix(apply(param_post, 2, mean), byrow = T, ncol = 5)
colnames(post_means) <- param
par(mfrow = c(K, 5))
for (k in 1:K) {
    for (j in 1:5) {
        print(paste0("$\\", param[j], "_", k, "$"))
        print(quantile(param_post[iter, (5*k-4):(5*k)][,j]), probs = c(.025, .05, .25, .5, .75, .95, .975))
        plot(iter, param_post[iter, (5*k-4):(5*k)][,j], type = "l",
             main = TeX(paste0("Chain for $\\", param[j], "_", k, "$")),
             xlab = "t", ylab = TeX(paste0("$\\", param[j], "_", k, "$")))
    }
}
dev.off()


pred_density_plot(sample_dat2, post_means = post_means, mix_props)
