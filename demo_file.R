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

# Single component Abe-Ley example from original Abe-Ley paper
true_par <- cbind(c(1,1,1), c(2,1,.5), c(1.5,1,0.5), c(6*pi/3,2*pi/3,4*pi/3)+pi/3, c(0,0,0), c(.5,.25,.25))
pars <- c("alpha", "beta", "kappa", "mu", "lambda", "tau")
colnames(true_par) <- pars
sample_dat2 <- rabeley_mixture(300, mix_prop = true_par[,"tau"],
                               mu = true_par[,"mu"], kappa = true_par[,"kappa"], lambda = true_par[,"lambda"],
                               alpha = true_par[,"alpha"], beta = true_par[,"beta"])

joint_dist_plot(sample_dat2, mean_pars = true_par[,1:5], mix_prop = true_par[,6], main = "True Abe-Ley Mixture Density")
true_par