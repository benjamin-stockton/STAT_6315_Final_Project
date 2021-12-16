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
source("~/Documents/Research/RCode/Schafer_Graham_Example/weighted_circular_regression_cl.R")
source("~/Documents/Research/RCode/Schafer_Graham_Example/imputation_funcs.R")
source("abe_ley_mixture_metropolis_hastings.R")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
theme_set(theme_classic())

dat <- read.csv("~/Documents/Research/RCode/Data/aqi_sample_data_42101_61103_61104_bdate_20140101_edate_20201231_state_55_county_079_site_0056.csv", header = T, stringsAsFactors = F) %>% dplyr::select(!X)

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

comp.adat <- adat[2:13,]
obs.x <- adat[1:2,"x"]

cdat.jan.2017 <- dat.2017 %>% filter(!is.na(sample_measurement_42101),
                                     !is.na(sample_measurement_61103),
                                     month(datetime_gmt) == 1) %>% dplyr::select(sample_measurement_61104, sample_measurement_42101)
x.obs <- dat.2017 %>% filter(is.na(sample_measurement_61104), 
                             !is.na(sample_measurement_42101), 
                             month(datetime_gmt) == 1) %>% dplyr::select(sample_measurement_42101)
x.obs <- x.obs$sample_measurement_42101
colnames(cdat.jan.2017) <- c("theta", "x")
summary(cdat.jan.2017)
cdat.jan.2017$theta %<>% as.numeric()

par(mfrow = c(1,1))
plot(cdat.jan.2017$theta, cdat.jan.2017$x, pch = 16,
     main = "Wind Directions and Carbon Monoxide Distribution for Jan. 2017",
     xlab = TeX("Wind Direction $\\theta$ (rad)"), ylab = "Carbon Monoxide")
# plot(comp.adat$theta, comp.adat$x, pch = 16, xlab = TeX("Wind Direction $\\theta$ (rad)"), ylab = "Carbon Monoxide")


source("abe_ley_mixture_metropolis_hastings.R")
ctrl <- list(Q = 2500, burnin = 1500, sd_init = .5)
K <- 1
# fit <- MH_posterior_estimation(comp.adat, K = K, control = ctrl)
fit <- MH_posterior_estimation(cdat.jan.2017, K = K, x.obs = x.obs, control = ctrl)

iters <- (ctrl$burnin+1):ctrl$Q
mix_props <- apply(fit$mix_props, 2, mean)
mix_props
post_means <- matrix(apply(fit$dist_pars[iters,], 2, mean), ncol = 5, byrow = T)
colnames(post_means) <- c("$\\alpha$", "$\\beta$", "$\\kappa$", "$\\mu$", "$\\lambda$")
post_means

plot_tracestack(fit, ctrl)

joint_dist_plot(cdat.jan.2017, mean_pars = post_means, mix_prop = mix_props, main = "Fitted Abe-Ley Mixture Density")
# joint_dist_plot(adat, mean_pars = post_means, mix_prop = mix_props, main = "Fitted Abe-Ley Mixture Density")

fit$dist_pars[,4] <- circular(fit$dist_pars[,4])
out_tab <- round(t(apply(fit$dist_pars, 2, quantile, probs = c(0.025, 0.05, .25,.5,.75,.95,.975))), 4)
rownames(out_tab) <- c("$\\alpha$", "$\\beta$", "$\\kappa$", "$\\mu$", "$\\lambda$")
xtable(out_tab, caption = "Posterior Quantiles of the Abe-Ley distribution from MCMC samples")

# alpha <- 2.3
# beta <- 4.94
# kappa <- post_means[1,3]
# mu <- post_means[1,4]
# lambda <- post_means[1,5]
# scale_par <- 1/(beta * (1 - tanh(kappa) * cos(0))^(1/alpha))
# scale_par * log(2)^(1/alpha)

s <- sample(iters, size = 1)
plot(fit$theta_pred[s,] %% (2*pi), x.obs, pch = 16,
     main = "Distribution of Draws from the Predictive Dist.",
     xlab = TeX("$\\theta_{pred}$"), ylab = "X")
for (i in 1:100) {
    s <- sample(iters, size = 1)
    points(fit$theta_pred[s,] %% (2*pi), x.obs)
}