library(circular)
library(mvtnorm)
library(circglmbayes)
library(kableExtra)
library(magrittr)
library(dplyr)
library(lubridate)
library(astsa)
library(ggplot2)
source("~/Documents/Research/RCode/Schafer_Graham_Example/weighted_circular_regression_cl.R")
source("~/Documents/Research/RCode/Schafer_Graham_Example/imputation_funcs.R")
source("abe_ley_mixture_metropolis_hastings.R")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
theme_set(theme_classic())

dat <- read.csv("~/Documents/Research/RCode/Data/aqi_sample_data_42101_61103_61104_bdate_20140101_edate_20201231_state_55_county_079_site_0056.csv", header = T, stringsAsFactors = F) %>% select(!X)

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
adat <- dat %>% filter(day(datetime_gmt) == 1 & month(datetime_gmt) %in% c(1,7) & hour(datetime_gmt) == 12)
adat

cdat.2017.2 <- cdat.2017 %>% select(sample_measurement_61104, sample_measurement_42101) %>% filter(!is.na(sample_measurement_61104), !is.na(sample_measurement_42101))
colnames(cdat.2017.2) <- c("theta", "x")

cdat.2017.2$theta %<>% as.numeric()
joint_dist_plot(cdat.2017.2, xlab = "$\\theta$", ylab = "X")

ctrl <- list(Q = 2500, burnin = 1250, sd_init = 1)
K <- 1
fit <- MH_posterior_estimation(cdat.2017.2, K = K, control = ctrl)

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


pred_density_plot(adat, post_means = post_means, mix_props)

