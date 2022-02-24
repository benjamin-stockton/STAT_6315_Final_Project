# Testing_Rcpp_sampler.R

library(Rcpp)
library(RcppGSL)
library(RcppArmadillo)
library(circular)
library(MCMCpack)
library(latex2exp)
library(RColorBrewer)
library(label.switching)

source("graphical_helpers.R")
source("abe_ley_mixture_metropolis_hastings.R")
sourceCpp("abe_ley_MM_sampler.cpp")

# true_par <- cbind(c(2,2,10), c(1, 0.2, 0.1), c(1,3,2), c(0, pi, pi/2), c(0.5,-1,.8), c(0,0,1))
# true_par <- cbind(c(1,2,10), c(0.07, 0.2, 0.01), c(1,3,2), c(0, pi, pi/2), c(0,-1,.8), c(1/3,1/3,1/3))
true_par <- cbind(c(1,2,10,3), c(0.07, 0.2, 0.04,1), c(1,3,2,3), c(0, pi, pi/2,2*pi/3), c(0,-1,.8,.5), c(.31, .18, .24, .27))
pars <- c("alpha", "beta", "kappa", "mu", "lambda", "tau")
colnames(true_par) <- pars; N <- 500
dat <- rabeley_mixture(N, mix_prop = true_par[,"tau"],
                               mu = true_par[,"mu"], kappa = true_par[,"kappa"], lambda = true_par[,"lambda"],
                               alpha = true_par[,"alpha"], beta = true_par[,"beta"])
source("graphical_helpers.R")
joint_dist_plot(dat, mean_pars = true_par[,1:5], mix_prop = true_par[,6],
                # main = "True Abe-Ley Mixture Density",
                main = "",
                nlevels = 20)
dev.off()
true_par

#########################################################################
# Individual Function Testing
#########################################################################

# Density Functions
sourceCpp("abe_ley_MM_sampler.cpp")

# dabeley()
dabeley_cpp(pi, 20, 2, 0.2, 3, pi, -2)
dabeley(20, pi, pi, 3, -1, 2, 0.2)
all.equal(dabeley_cpp(pi, 20, 2, 0.2, 3, pi, -2),
          dabeley(20, pi, pi, 3, -1, 2, 0.2))

# dbeta_rescale()
dbeta_rescale(0, 5, .5)
dbeta_rescale_cpp(0, 5, 0.5)
all.equal(dbeta_rescale(0, 5, .5),
          dbeta_rescale_cpp(0, 5, 0.5))

# dwrappednormal()
sigma <- c(.5, 1, 2, 2.5, 4, 8)
cpp_probs <- numeric(length(sigma))
r_probs <- numeric(length(sigma))
for (i in seq_len(length(sigma))) {
    cpp_probs[i] <- dwrappednorm_cpp(pi, pi/2, sigma[i])
    r_probs[i] <- dwrappednormal(circular(pi), circular(pi/2), sd = sigma[i])
}
cpp_probs; r_probs
all.equal(cpp_probs, r_probs)

# dtruncnorm()
dtruncnorm(10, 0, Inf, 5, 10)
dtruncnorm_cpp(10, 0, Inf, 5, 10)
all.equal(dtruncnorm(10, 0, Inf, 5, 10),
          dtruncnorm_cpp(10, 0, Inf, 5, 10))

dtruncnorm(.5, -1,1, 0, 2)
dtruncnorm_cpp(.5, -1,1, 0, 2)
all.equal(dtruncnorm(.5, -1,1, 0, 2),
          dtruncnorm_cpp(.5, -1,1, 0, 2))

# RNG Functions
sourceCpp("abe_ley_MM_sampler.cpp")

# rtruncnorm()
set.seed(1099)
tn_cpp <- rtruncnorm_cpp(1000, 0, Inf, 5, 10)
set.seed(1099)
tn_r <- rtruncnorm(1000, 0, Inf, 5, 10)
tn_cpp[1:10]; tn_r[1:10]
all.equal(tn_cpp, tn_r)

plot(density(tn_cpp), xlim = c(min(c(tn_cpp, tn_r)), max(c(tn_cpp, tn_r))))
lines(density(tn_r), col = "red")
x_seq <- seq(min(c(tn_cpp, tn_r)), max(c(tn_cpp, tn_r)), length.out = 1000)
lines(x_seq, dtruncnorm(x_seq, 0, Inf, 5,10), col = "blue", lty = "dashed")

set.seed(1099)
tn_cpp <- rtruncnorm_cpp(1000, -1, 1, 0, 3)
set.seed(1099)
tn_r <- rtruncnorm(1000, -1, 1, 0, 3)
tn_cpp[1:10]; tn_r[1:10]
all.equal(tn_cpp, tn_r)

plot(density(tn_cpp))
lines(density(tn_r), col = "red")
x_seq <- seq(-1,1, length.out = 1000)
lines(x_seq, dtruncnorm(x_seq, -1,1,0, 3), col = "blue", lty = "dashed")

# rbeta_rescale()
set.seed(899)
b_cpp <- rbeta_rescale_cpp(1000, 5, .5)[,1]
set.seed(899)
b_r <- rbeta_rescale(1000, 5, .5)
b_cpp[1:10]; b_r[1:10]
all.equal(b_cpp, b_r)

# rwrappednormal()
set.seed(7677)
wn_cpp <- rwrappednorm_cpp(1000, pi, 3)[,1] %% (2*pi)
set.seed(7677)
wn_r <- as.numeric(rwrappednormal(1000, circular(pi), sd = 3))
wn_cpp[1:10]; wn_r[1:10]
all.equal(wn_cpp, wn_r)

plot(density(wn_cpp))
lines(density(wn_r), col = "red")
x_seq <- seq(0,2*pi, length.out = 1000)
lines(x_seq, dwrappednormal(circular(x_seq), circular(pi), sd = 3), col = "blue", lty = "dashed")

# rwrappedcauchy()
set.seed(7677)
wc_cpp <- rwrappedcauchy_cpp(1000, pi, 1)[,1]
set.seed(7677)
wc_r <- as.numeric(rwrappedcauchy(1000, circular(pi), rho = exp(-1)))
wc_cpp[1:10]; wc_r[1:10]
all.equal(wn_cpp, wc_r)

plot(density(wc_cpp))
lines(density(wc_r), col = "red")
x_seq <- seq(0,2*pi, length.out = 1000)
lines(x_seq, dwrappedcauchy(circular(x_seq), circular(pi), rho = exp(-1)), col = "blue", lty = "dashed")

# rssvm()
set.seed(7677)
ssvm_cpp <- rssvm_cpp(1000, pi, 3, .5)
set.seed(7677)
ssvm_r <- rssvm(1000, pi, 3, .5)
ssvm_cpp[1:10]; ssvm_r[1:10]
all.equal(ssvm_cpp, ssvm_r)

plot(density(ssvm_cpp))
lines(density(ssvm_r), col = "red")
x_seq <- seq(0,2*pi, length.out = 1000)
lines(x_seq, dssvm(x_seq, pi, 3, .5), col = "blue", lty = "dashed")

# rabeley()
set.seed(494)
ab_cpp <- rabeley_cpp(1000, 2, .2, 3, pi, -1)
set.seed(494)
ab_r <- rabeley(1000, pi, 3, -1, 2, .2)
ab_cpp[1:10,]; ab_r[1:10,]
all.equal(ab_cpp, ab_r)

plot(ab_cpp[,1] %% (2*pi), ab_cpp[,2], col = "blue")
points(ab_r[,1], ab_r[,2], col = "red")

# rssvm_mixture()
set.seed(789)
ssvmm_cpp <- rssvm_mixture_cpp(1000, c(.5, .5), c(pi, pi/2), c(25, 15), c(-.75,-.75))[,1] %% (2*pi)
set.seed(789)
ssvmm_r <- rssvm_mixture(1000, c(.5, .5), c(pi, pi/2), c(25, 15), c(-.75,-.75))
all.equal(ssvmm_cpp, ssvmm_r)

plot(density(ssvmm_cpp))
lines(density(ssvmm_r), col = "red")
x_seq <- seq(0,2*pi, length.out = 1000)
lines(x_seq, dssvm(x_seq, pi, 3, .5), col = "blue", lty = "dashed")

# Sampler Helper Functions
sourceCpp("abe_ley_MM_sampler.cpp")

# count_cluster_sizes()
cnts <- rep(1:4, 25)
count_cluster_size(cnts, 4)

# sample_proposal_dist()
set.seed(5677)
sp_cpp <- sample_proposal_dist_cpp(c(1,1,1,pi,0), c(1,1,1,.5,.5))[,1]
set.seed(5677)
sp_r <- sample_proposal_dist(c(1,1,1,pi,0), c(1,1,1,.5,.5))
sp_cpp; sp_r
all.equal(sp_cpp, sp_r)

# accept_ratio()
accept_ratio_cpp(dat, sp_r, sp_cpp, c(1,1,1,.5,.25), list(alpha_shape = 1000, alpha_scale = 0.001,
                                                          beta_shape = 1000, beta_scale = 0.001, 
                                                          kappa_shape = 1000, kappa_scale = 0.001,
                                                          mu_mean = 0, mu_conc = 0.001, 
                                                          lambda_a = 1, lambda_b = 1, alpha0 = 1))[,1]
accept_ratio(dat, sp_r, sp_cpp, c(1,1,1,.5,.25))

# sample_class_labels()
nu_t <- true_par
set.seed(16897)
s_cpp <- sample_class_labels_cpp(dat, nu_t, 3)[,1]
set.seed(16897)
s_r <- sample_class_labels(dat, nu_t[,1:5], 3, nu_t[,6])
all.equal(s_cpp, s_r)

# adjust_sd()
a_probs <- cbind(c(.8, .8, .8), c(.4,.6,.8), c(.2,.2,.2), c(0,.2,.2), c(1,1,1))
adjust_sd_cpp(5, rep(5,5), 50, a_probs)[,1]
adjust_sd(5, rep(5,5), 50, a_probs)

adjust_sd_cpp(50, rep(.5,5), 50, a_probs)[,1]
adjust_sd(50, rep(.5,5), 50, a_probs)


# Sampler
Q <- 50000
K <- 4
hp <- list(alpha_shape = .5, alpha_scale = 2,
           beta_shape = .5, beta_scale = 2,
           kappa_shape = .5, kappa_scale = 2,
           mu_mean = 0, mu_conc = 0.001,
           lambda_a = 1, lambda_b = 1, alpha0 = 1/K)
ctrl <- list(Q = Q, burnin = floor(Q/2), sd_init = 5, thin = 5)
# sourceCpp("abe_ley_MM_sampler.cpp")
fit <- abeley_MM_sampler_MH_cpp(dat, K = K, control = ctrl, hyperpar = hp)

# plot_tracestack(mcmc_pars = fit$chains, K = K, ctrl)

iters <- floor((ctrl$Q-15000)/ctrl$thin):floor(ctrl$Q/ctrl$thin) * ctrl$thin

mcmc_pars <- fit$chains[iters,,]
mix_props <- apply(matrix(mcmc_pars[,,6], ncol = K), 2, mean)
post_means <- matrix(numeric(5*K), ncol = 5)
for (i in 1:5) {
    post_means[,i] <- apply(matrix(mcmc_pars[,,i], ncol = K), 2, median)
}
colnames(post_means) <- pars[1:5]
print(cbind(post_means, mix_props)); true_par

plot_tracestack(mcmc_pars = mcmc_pars, K = K, ctrl)
plot_post_density(mcmc_pars = mcmc_pars, K = K, ctrl)

source("graphical_helpers.R")
joint_dist_plot(dat, mean_pars = post_means, mix_prop = mix_props, main = "Fitted Abe-Ley Mixture Density", nlevels = 15)
joint_dist_plot(dat, mean_pars = true_par[,1:5], mix_prop = true_par[,6],
                main = "True Abe-Ley Mixture Density", nlevels = 15)

z <- fit$class_labels[iters,]
run <- label.switching::dataBased(x = dat, K = K, z = z) 
relabeled.mcmc <- label.switching::permute.mcmc(mcmc_pars, run$permutations)$output

mix_props <- apply(matrix(relabeled.mcmc[,,6], ncol = K), 2, mean)
post_means <- matrix(numeric(5*K), ncol = 5)
for (i in 1:5) {
    post_means[,i] <- apply(matrix(relabeled.mcmc[,,i], ncol = K), 2, mean)
}
colnames(post_means) <- pars[1:5]
print(cbind(post_means, mix_props)); true_par

plot_tracestack(mcmc_pars = relabeled.mcmc[,,], K = K, ctrl)
plot_post_density(mcmc_pars = relabeled.mcmc[,,], K = K, ctrl)

joint_dist_plot(dat, mean_pars = true_par, mix_prop = true_par[,6], main = "True Abe-Ley Mixture Density", nlevels = 15)
joint_dist_plot(dat, mean_pars = post_means, mix_prop = mix_props, main = "Fitted Abe-Ley Mixture Density", nlevels = 15)
