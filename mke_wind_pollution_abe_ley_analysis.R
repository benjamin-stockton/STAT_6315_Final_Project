library(circular)
library(mvtnorm)
library(circglmbayes)
library(kableExtra)
library(magrittr)
library(dplyr)
library(lubridate)
library(astsa)
library(ggplot2)
library(xtable)
library(foreach)
library(abind)
library(Rcpp)
library(RcppGSL)
library(RcppArmadillo)
source("~/Documents/Research/CL-Regression_MI/weighted_circular_regression_cl.R")
source("~/Documents/Research/CL-Regression_MI/imputation_funcs.R")
source("abe_ley_mixture_metropolis_hastings.R")
source("graphical_helpers.R")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
theme_set(theme_classic())

dat <- read.csv("~/Documents/Research/Data/aqi_sample_data_42101_61103_61104_bdate_20140101_edate_20201231_state_55_county_079_site_0056.csv", header = T, stringsAsFactors = F) %>% dplyr::select(!X)

# Data Overview

head(dat)

# data type conversions for dates and circulars
dat$datetime_gmt %<>% ymd_hm()
dat$sample_measurement_61104 %<>% circular(., type = "directions", units = "degrees") %>% conversion.circular("radians")

# sort by datetime
dat %<>% arrange(datetime_gmt)
dat %>% tail(n= 10)

# summaries to see how many missing values there are
summary(dat$sample_measurement_42101)
##  1577 missing values for CO
summary(dat$sample_measurement_61103)
##  8050 missing values for both wind measurements
summary(dat$sample_measurement_61104)

# All of the missing wind observations are missing in the same rows. 
## WS is missing iff WD is missing
sum(is.na(dat[is.na(dat$sample_measurement_61104),"sample_measurement_61103"]))

# Now only look at 2017

int_2017 <- interval(start = "2017-01-01", end = "2017-12-31")
dat.2017 <- dat %>% filter(datetime_gmt %within% int_2017)

head(dat.2017)

# summaries to see how many missing values there are in 2017
summary(dat.2017$sample_measurement_42101)
##  475 missing values for CO
summary(dat.2017$sample_measurement_61103)
##  75 missing values for both wind measurements
summary(dat.2017$sample_measurement_61104)

# Then I'll make a 1-hour lagged variable for CO
dat.2017$co_lag1 <- dat.2017$sample_measurement_42101 %>% lag()

plot(dat.2017$sample_measurement_42101, dat.2017$co_lag1)
# There's a strong linear relationship between CO from one hour to the next, although this tends to be heteroscedastic (higher variance for higher earlier measurements)
cor(dat.2017$sample_measurement_42101, dat.2017$co_lag1, use = "complete.obs")

# Complete Data exploration
mis <- which(is.na(dat.2017$sample_measurement_42101))
table(dat.2017[mis, "date_gmt"])
mis <- which(is.na(dat.2017$sample_measurement_61103))
table(dat[mis, "date_gmt"])

cdat.2017 <- dat.2017 %>% filter(!is.na(sample_measurement_42101) & !is.na(co_lag1) & !is.na(sample_measurement_61103))
cdat.2017 %>% head(10)

#####################################################
# Annual Data
## Taking the measurements from Jan 1 and Jul 1 each year at noon
dat$co_lag1 <- dat$sample_measurement_42101 %>% lag()
adat <- dat %>% filter(day(datetime_gmt) == 1,
                       month(datetime_gmt) %in% c(1,7),
                       hour(datetime_gmt) == 12,
                       # !is.na(sample_measurement_61104),
                       !is.na(sample_measurement_42101)) %>% dplyr::select(sample_measurement_61104, sample_measurement_42101)
colnames(adat) <- c("theta", "x")
adat$theta %<>% as.numeric()
adat

comp.adat <- adat[3:13,]
obs.x <- adat[1:2,"x"]

cdat.apr.2017 <- dat.2017 %>% filter(!is.na(sample_measurement_42101),
                                     !is.na(sample_measurement_61103),
                                     month(datetime_gmt) == 4) %>% dplyr::select(sample_measurement_61104, sample_measurement_42101)
x.obs <- dat.2017 %>% filter(is.na(sample_measurement_61104), 
                             !is.na(sample_measurement_42101), 
                             month(datetime_gmt) == 4) %>% dplyr::select(sample_measurement_42101)
x.obs <- x.obs$sample_measurement_42101
colnames(cdat.apr.2017) <- c("theta", "x")
summary(cdat.apr.2017)
cdat.apr.2017$theta %<>% as.numeric()

par(mfrow = c(1,1))
plot(cdat.apr.2017$theta, cdat.apr.2017$x, pch = 16,
     main = "Wind Directions and Carbon Monoxide Distribution for April 2017",
     xlab = TeX("Wind Direction $\\theta$ (rad)"), ylab = "Carbon Monoxide")
# plot(comp.adat$theta, comp.adat$x, pch = 16, xlab = TeX("Wind Direction $\\theta$ (rad)"), ylab = "Carbon Monoxide")


source("abe_ley_mixture_metropolis_hastings.R")
sourceCpp("abe_ley_MM_sampler.cpp")
# Q <- 5000
# ctrl <- list(Q = Q, burnin = floor(Q/2), sd_init = 1, thin = 5)
# K <- 4
# fit <- foreach::foreach(i = 1:4, .packages = c("circular", "mvtnorm", "MASS", "MCMCpack", "truncnorm")) %dopar% {
#     # try(MH_posterior_estimation(as.matrix(comp.adat), K = K, control = ctrl), silent = F)
#     try(MH_posterior_estimation(as.matrix(cdat.apr.2017), K = K, control = ctrl), silent = F)
# }

Q <- 50000
K <- 2
hp <- list(alpha_shape = .5, alpha_scale = 2,
           beta_shape = .5, beta_scale = 2,
           kappa_shape = .5, kappa_scale = 2,
           mu_mean = 0, mu_conc = 0.001,
           lambda_a = 1, lambda_b = 1, alpha0 = 1/K)
ctrl <- list(Q = Q, burnin = floor(Q/2), sd_init = 5, thin = 5)
# sourceCpp("abe_ley_MM_sampler.cpp")
fit <- abeley_MM_sampler_MH_cpp(as.matrix(cdat.apr.2017), K = K, control = ctrl, hyperpar = hp)

# plot_tracestack(mcmc_pars = fit$chains, K = K, ctrl)

iters <- floor((ctrl$Q-15000)/ctrl$thin):floor(ctrl$Q/ctrl$thin) * ctrl$thin

# iters <- floor((ctrl$burnin+1)/thin):floor(ctrl$Q/thin)
# for (i in 1:2) {
#     if (is.vector(fit[[i]])) {
#         guaranteed_fit <- i
#     }
# }
# mcmc_pars <- fit[[guaranteed_fit]]$mcmc_pars[iters,,]
# for (i in 1:2) {
#     if (is.vector(fit[[i]]) & i != guaranteed_fit) {
#         mcmc_pars <- abind(mcmc_pars, fit[[i]]$mcmc_pars[iters,,], along = 1)
#     }
# }
mcmc_pars <- fit$chains[iters,,]
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

# joint_dist_plot(comp.adat, mean_pars = post_means, mix_prop = mix_props, main = "Fitted Abe-Ley Mixture Density", nlevels = 15)

joint_dist_plot(cdat.apr.2017, mean_pars = post_means, mix_prop = mix_props, main = "Fitted Abe-Ley Mixture Density", nlevels = 15)

# Label-switching Fix

z <- fit$class_labels[iters,]
# z <- fit[[guaranteed_fit]]$class_labels[iters,]
# for (i in 1:2) {
#     if (is.vector(fit[[i]]) & i != guaranteed_fit) {
#         z <- abind(z, fit[[i]]$class_labels[iters,], along = 1)
#     }
#     
# }
run <- label.switching::dataBased(x = as.matrix(cdat.apr.2017), K = K, z = z) 
relabeled.mcmc <- label.switching::permute.mcmc(mcmc_pars, run$permutations)$output

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

# joint_dist_plot(comp.adat, mean_pars = post_means, mix_prop = mix_props, main = "Fitted Abe-Ley Mixture Density", nlevels = 15)
plot(cdat.apr.2017$theta, cdat.apr.2017$x)
fl <- kde2d(cdat.apr.2017$theta, cdat.apr.2017$x, lims = c(0, 2*pi, 0, max(cdat.apr.2017$x)))
contour(fl, xlab = "Angle", ylab = "CO2", nlevels = 7, add = T, lwd = 1.5)
joint_dist_plot(cdat.apr.2017, mean_pars = post_means, mix_prop = mix_props, main = "Fitted Abe-Ley Mixture Density", nlevels = 15)

