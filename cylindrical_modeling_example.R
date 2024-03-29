library(mvtnorm)
library(circular)
library(MASS)
library(latex2exp)
library(truncnorm)
library(MCMCpack)
library(label.switching)
library(foreach)
library(doParallel)
library(abind)
# numCores <- detectCores() - 1
# numCores
source("abe_ley_mixture_metropolis_hastings.R")
source("graphical_helpers.R")

# true_par <- cbind(c(2,2,10), c(1, 0.2, 0.1), c(1,3,2), c(0, pi, pi/2), c(0.5,-1,.8), c(0,0,1))
true_par <- cbind(c(1,2,10), c(0.07, 0.2, 0.01), c(1,3,2), c(0, pi, pi/2), c(0,-1,.8), c(1/3,1/3,1/3))
# true_par <- cbind(c(1,2,10,3), c(0.07, 0.2, 0.04,1), c(1,3,2,3), c(0, pi, pi/2,2*pi/3), c(0,-1,.8,.5), c(.31, .18, .24, .27))
pars <- c("alpha", "beta", "kappa", "mu", "lambda", "tau")
colnames(true_par) <- pars; N <- 500
sample_dat2 <- rabeley_mixture(N, mix_prop = true_par[,"tau"],
                               mu = true_par[,"mu"], kappa = true_par[,"kappa"], lambda = true_par[,"lambda"],
                               alpha = true_par[,"alpha"], beta = true_par[,"beta"])
joint_dist_plot(sample_dat2, mean_pars = true_par[,1:5], mix_prop = true_par[,6],
                main = "True Abe-Ley Mixture Density",
                nlevels = 15)
true_par
# log(dabeley(x = sample_dat2[,2], t = sample_dat2[,1], alpha = 2, beta = 0.2, kappa = 3, mu = pi, lambda = -1))
Q <- 5000
K <- 3
hp <- list(alpha_shape = 1000, alpha_scale = 0.001,
           beta_shape = 1000, beta_scale = 0.001,
           kappa_shape = 1000, kappa_scale = 0.001,
           mu_mean = 0, mu_conc = 0.001,
           lambda_a = 1, lambda_b = 1, alpha0 = 1/K)
ctrl <- list(Q = Q, burnin = floor(Q/2), sd_init = 1, thin = 5)
source("abe_ley_mixture_metropolis_hastings.R")
fit <- foreach::foreach(i = 1:2, .packages = c("circular", "mvtnorm", "MASS", "MCMCpack", "truncnorm")) %dopar% {
    try(MH_posterior_estimation(sample_dat2, K = K, control = ctrl, hyperpar = hp), silent = F)
}
# fit <- MH_posterior_estimation(sample_dat2, K = K, control = ctrl, hyperpar = hp)

iters <- floor((ctrl$burnin+1)/ctrl$thin):floor(ctrl$Q/ctrl$thin)

for (i in 1:2) {
    if (is.vector(fit[[i]])) {
        guaranteed_fit <- i
    }
}
mcmc_pars <- fit[[guaranteed_fit]]$mcmc_pars[iters,,]
for (i in 1:2) {
    if (is.vector(fit[[i]]) & i != guaranteed_fit) {
        mcmc_pars <- abind(mcmc_pars, fit[[i]]$mcmc_pars[iters,,], along = 1)
    }
        
}

mix_props <- apply(matrix(mcmc_pars[,,6], ncol = K), 2, mean)
mix_props
post_means <- matrix(numeric(5*K), ncol = 5)
for (i in 1:5) {
    post_means[,i] <- apply(matrix(mcmc_pars[,,i], ncol = K), 2, median)
}
colnames(post_means) <- pars[1:5]
print(post_means)

plot_tracestack(mcmc_pars = mcmc_pars, K = K, ctrl)
plot_post_density(mcmc_pars = mcmc_pars, K = K, ctrl)

joint_dist_plot(sample_dat2, mean_pars = true_par, mix_prop = true_par[,6], main = "True Abe-Ley Mixture Density", nlevels = 15)
joint_dist_plot(sample_dat2, mean_pars = post_means, mix_prop = mix_props, main = "Fitted Abe-Ley Mixture Density", nlevels = 15)

# Label-switching Fix

# mcmc.pars <- array(data = NA, dim = c(length(iters), K, 6))
# for (j in 1:6) {
#     # par_names <- colnames(fit$dist_pars[iters,0:(K-1) * 5 +j])
#     mcmc.pars[,,j] <- mcmc_pars[,,j]
#     # dimnames(mcmc.pars[,,j]) <- par_names
# }
# mcmc.pars[1:10,,]
z <- fit[[guaranteed_fit]]$class_labels[iters,]
for (i in 1:2) {
    if (is.vector(fit[[i]]) & i != guaranteed_fit) {
        z <- abind(z, fit[[i]]$class_labels[iters,], along = 1)
    }
    
}
# p <- array(data = NA, dim = c(nrow(z), ncol(z), K))
# 
# for (q in 1:length(iters)) {
#     for (k in 1:K) {
#         par_1 <- fit$mcmc_pars[q,k,]
#         par_mix <- fit$mcmc_pars[q,,]
#         p[q,,k] <- par_1[6] * dabeley(x = sample_dat2[,"x"], theta = sample_dat2[,"theta"],
#                            mu = par_1[4], kappa = par_1[3], lambda = par_1[5],
#                            alpha = par_1[1], beta = par_1[2]) / (dabeley_mixture(x = sample_dat2[,"x"], theta = sample_dat2[,"theta"], mix_prop = par_mix[,6], mu = par_mix[,4], kappa = par_mix[,3], lambda = par_mix[,5], alpha = par_mix[,1], beta = par_mix[,2]))
#     }
# }
# 
# complete.abeley.loglikelihood <- function(dat, z, pars) {
#     # dat should be 2xn with col1 = theta and col2 = x
#     # z should be Qxn allocation vector matrix
#     # pars should be an individual iterations mcmc draw
#     K <- dim(pars)[1]
#     n <- nrow(dat)
#     logl <- rep(0, n)
# 
#     log_theta <- log(pars[,6])
# 
#     alpha <- pars[,1]
#     beta <- pars[,2]
#     kappa <- pars[,3]
#     mu <- pars[,4]
#     lambda <- pars[,5]
#     logl <- log_theta[z] + log(dabeley(x = dat[,2], theta = dat[,1],
#                                    mu = mu[z], kappa = kappa[z], lambda = lambda[z],
#                                    alpha = alpha[z], beta = beta[z]))
#     return(sum(logl))
# }

# run <- label.switching::sjw(mcmc.pars = mcmc.pars, z = z, complete = complete.abeley.loglikelihood, x = sample_dat2, init = 1) # Doesn't Converge
run <- label.switching::dataBased(x = sample_dat2, K = K, z = z) # Doesn't work
# run <- label.switching::ecr.iterative.1(z = z, K = K) # Doesn't work
# run <- label.switching::ecr.iterative.2(z = z, K = K, p = p) # Doesn't converge
# print(run$status)
relabeled.mcmc <- label.switching::permute.mcmc(mcmc_pars, run$permutations)$output
# for (j in 1:5) {
#     par_names <- colnames(fit$dist_pars[iters,0:(K-1) * 5 +j])
#     dimnames(relabeled.mcmc[,,j]) <- par_names
# }

relabeled.mcmc[1:10,,4]

mix_props <- apply(matrix(relabeled.mcmc[,,6], ncol = K), 2, mean)
mix_props
post_means <- matrix(numeric(5*K), ncol = 5)
for (i in 1:5) {
    post_means[,i] <- apply(matrix(relabeled.mcmc[,,i], ncol = K), 2, mean)
}
colnames(post_means) <- pars[1:5]
print(post_means)

plot_tracestack(mcmc_pars = relabeled.mcmc[,,], K = K, ctrl)
plot_post_density(mcmc_pars = relabeled.mcmc[,,], K = K, ctrl)

joint_dist_plot(sample_dat2, mean_pars = true_par, mix_prop = true_par[,6], main = "True Abe-Ley Mixture Density", nlevels = 15)
joint_dist_plot(sample_dat2, mean_pars = post_means, mix_prop = mix_props, main = "Fitted Abe-Ley Mixture Density", nlevels = 15)

# Now trying it with missing data
# 10% missing in theta
# Single component Abe-Ley example from original Abe-Ley paper
true_par <- cbind( c(2), c(1), c(1), c(pi), c(0.5), c(1))
# true_par <- cbind(c(2,2,2), c(1,1,1), c(2,2,2), c(0, 2*pi/3, 4*pi/3)+pi/6, c(-.5, 0, .5), c(.33, .33, .33))
pars <- c("alpha", "beta", "kappa", "mu", "lambda", "tau")
colnames(true_par) <- pars
sample_dat2 <- rabeley_mixture(150, mix_prop = true_par[,6],
                               mu = true_par[,4], kappa = true_par[,3], lambda = true_par[,5],
                               alpha = true_par[,1], beta = true_par[,2])
mis <- sample(1:nrow(sample_dat2), size = 15, replace = F)

comp_sample <- sample_dat2[!mis, ]
obs.x <- sample_dat2[mis,"x"]

source("abe_ley_mixture_metropolis_hastings.R")
ctrl <- list(Q = 2500, burnin = 1500, sd_init = 1)
K <- 1
fit <- MH_posterior_estimation(sample_dat2, K = K, x.obs = obs.x, control = ctrl)

iters <- (ctrl$burnin+1):ctrl$Q
mix_props <- apply(matrix(fit$mcmc_pars[iters,,6], ncol = K), 2, mean)
mix_props
post_means <- matrix(numeric(5*K), ncol = 5)
for (i in 1:5) {
    post_means[,i] <- apply(matrix(fit$mcmc_pars[iters,,i], ncol = K), 2, mean)
}
colnames(post_means) <- pars[1:5]
print(post_means)

plot_tracestack(mcmc_pars = fit$mcmc_pars[iters,,], K = K, ctrl)

joint_dist_plot(sample_dat2, mean_pars = matrix(true_par[,1:5], nrow = K), mix_prop = matrix(true_par[,6], nrow = K), main = "True Abe-Ley Mixture Density")
joint_dist_plot(sample_dat2, mean_pars = post_means, mix_prop = 1, main = "Fitted Abe-Ley Mixture Density")

s <- sample(iters, size = 1)
theta_pred <- fit$theta_pred[s,] %% (2*pi)

par(mfrow = c(1,1))
plot(sample_dat2[mis, "theta"], obs.x, pch = 16,
     main = "Distribution of Predictive Posterior Draws",
     xlab = TeX("$\\theta_{pred}$"), ylab = "X")
legend("topright", legend = c("Pred. Post. Draws", "Obs X"),
       col = c("black", "red"),
       pch = c(1,16))
for (i in 1:100) {
    s <- sample(iters, size = 1)
    points(fit$theta_pred[s,] %% (2*pi), obs.x)
}
points(sample_dat2[mis, "theta"], obs.x, pch = 4, col = "red", lwd = 2)


