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

###################################################
# SSVM sampling test
###################################################

mu <- pi/2; kappa <- 4; lambda <- 1
y_seq <- seq(0,2*pi, length.out = 200)
y <- rssvm(100, mu = mu, kappa = kappa, lambda = lambda) 
par(mfrow = c(1,2))
plot(density(circular(y), bw = 10), shrink = 1.25, main = paste0("SSVM; mu = ", round(mu, 3)))
points(circular(y))

t_seq <- seq(0, 2*pi, length.out = 50)
d <- density(y)

hist(y, breaks = 25, freq = F, main = paste0("SSVM; mu = ", round(mu, 3)))
lines(d$x, d$y, lwd = 2)
lines(t_seq, dssvm(t_seq, mu, kappa, lambda), col = "red", lwd = 2)
abline(v = mu * c(-1,1), col = 'blue', lty = "dashed")

###################################################
# SSWC
###################################################

par(mfrow = c(1,2))
M <- 1000
mu <- runif(1, 0, 2*pi); kappa <- 2; lambda <- 1
t <- rsswc(M, mu, kappa, lambda)

plot(density(circular(t, modulo = "asis"), bw = tanh(kappa/2)), shrink = 1.5, main = paste0("SSWC; mu = ", round(mu, 3)))
points(circular(t), stack = T)

t_seq <- seq(0, 2*pi, length.out = 50)
d <- density(t)

hist(t, breaks = 25, freq = F, main = paste0("SSWC; mu = ", round(mu, 3)))
lines(d$x, d$y, lwd = 2)
lines(t_seq, dsswc(t_seq, mu, kappa, lambda), col = "red", lwd = 2)
abline(v = mu * c(-1,1), col = 'blue', lty = "dashed")

#########################################
# Abe-Ley Distribution
#########################################

true_par <- cbind(c(1,2,10,3), c(0.07, 0.2, 0.04,1), c(1,3,2,3), c(0, pi, pi/2,2*pi/3), c(0,-1,.8,.5), c(.31, .18, .24, .27))
pars <- c("alpha", "beta", "kappa", "mu", "lambda", "tau")
colnames(true_par) <- pars; N <- 3000

dat <- rbind(rabeley(N/4, mu = true_par[1, "mu"], kappa = true_par[1,"kappa"], lambda = true_par[1,"lambda"], alpha = true_par[1,"alpha"], beta = true_par[1,"beta"]),
             rabeley(N/4, mu = true_par[2, "mu"], kappa = true_par[2,"kappa"], lambda = true_par[2,"lambda"], alpha = true_par[2,"alpha"], beta = true_par[2,"beta"]),
             rabeley(N/4, mu = true_par[3, "mu"], kappa = true_par[3,"kappa"], lambda = true_par[3,"lambda"], alpha = true_par[3,"alpha"], beta = true_par[3,"beta"]),
             rabeley(N/4, mu = true_par[4, "mu"], kappa = true_par[4,"kappa"], lambda = true_par[4,"lambda"], alpha = true_par[4,"alpha"], beta = true_par[4,"beta"]))


x_seq <- seq(from = min(dat[,2]), to = max(dat[,2])+10, length.out = 50)
t_seq <-  seq(from = min(dat[,1]), to = max(dat[,1]), length.out = 50)

z1 <- t(outer(X = x_seq, Y = t_seq, FUN = dabeley_mixture,
              mix_prop = true_par[1,6],
              mu = true_par[1, "mu"], kappa = true_par[1,"kappa"], lambda = true_par[1,"lambda"],
              alpha = true_par[1,"alpha"], beta = true_par[1,"beta"]))
z2 <- t(outer(X = x_seq, Y = t_seq, FUN = dabeley_mixture,
              mix_prop = true_par[2,6],
              mu = true_par[2, "mu"], kappa = true_par[2,"kappa"], lambda = true_par[2,"lambda"],
              alpha = true_par[2,"alpha"], beta = true_par[2,"beta"]))
z3 <- t(outer(X = x_seq, Y = t_seq, FUN = dabeley_mixture,
              mix_prop = true_par[3,6],
              mu = true_par[3, "mu"], kappa = true_par[3,"kappa"], lambda = true_par[3,"lambda"],
              alpha = true_par[3,"alpha"], beta = true_par[3,"beta"]))
z4 <- t(outer(X = x_seq, Y = t_seq, FUN = dabeley_mixture,
              mix_prop = true_par[4,6],
              mu = true_par[4, "mu"], kappa = true_par[4,"kappa"], lambda = true_par[4,"lambda"],
              alpha = true_par[4,"alpha"], beta = true_par[4,"beta"]))
z.list <- list(z1, z2, z3, z4)

par(mfrow = c(2,2))
cols <- c("red", "blue", "green", "orange")
for (i in 1:nrow(true_par)) {
    nlvls <- 20
    bounds <- c(min(dat[(N*(i-1)/4+1):(N*(i)/4),"x"]), max(dat[(N*(i-1)/4+1):(N*(i)/4),"x"])+10)
    contour(x = t_seq, y = x_seq, z = z.list[[i]], lwd = 2, add = FALSE, nlevels = nlvls,
            col = hcl.colors(nlvls, "Spectral"),
            ylim = bounds)
    points(dat[(N*(i-1)/4+1):(N*(i)/4),"theta"], dat[(N*(i-1)/4+1):(N*(i)/4),"x"], col = cols[i], cex = .75)
    print(paste0(c("alpha = ", "beta = ", "kappa = ", "mu = ", "lambda = ", "tau = "), round(true_par[i,], 2)))
}

