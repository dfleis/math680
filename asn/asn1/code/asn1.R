#==========
#
# Main file for Assignment 1
# 
#==========
#===== Load external functions/libraries =====#
source("code/asn1_functions.R")
source("code/asn1_functions_tests.R")
library(doParallel) # for Q4

#===== Q2 =====#

## (e)
set.seed(124)

# set parameters
nsims <- 1e3
n <- 25
p <- 7
lam <- 4
beta_star <- 1:p
sigma_star <- 1

# generate data
X <- cbind(1, matrix(rnorm(n * (p - 1)), nrow = n))

# compute theoretical mean and variance
par_true <- ridge_coef_params(X, lam, beta_star, sigma_star)
b_true <- as.vector(par_true$b)
vcv_true <- par_true$vcv

# simulate ridge coefficients nsims times
# outputs a matrix with rows corresponding to coefficients # and columns correspond to simulation number
b_hat <- replicate(nsims, {
  y <- X %*% beta_star + rnorm(n, 0, sigma_star)
  as.vector(ridge_coef(X, y, lam)$b) 
})

# estimate variance of b1, ..., b_p estimates
vcv_hat <- var(t(b_hat))

# print estimated fused ridge coefficients vs. expected values
b <- rbind(rowMeans(b_hat), b_true) 
rownames(b) <- c("b_hat", "b_true") 
round(b, 4)

# print absolute error between estimated and true fused ridge variances
round(abs(vcv_true - vcv_hat), 4)

#===== Q4 =====#
set.seed(124)
  
# global parameters
nsims <- 4
n <- 100
lams <- 10^seq(-8, 8, 0.5)
Ks <- c(5, 10, n)
sigma_star <- sqrt(1/2)

## (a)
# set parameters
p <- 50
theta <- 0.5

# generate data
beta_star <- rnorm(p, 0, sigma_star)
SIGMA <- outer(1:(p - 1), 1:(p - 1), FUN = function(a, b) theta^abs(a - b))
X <- cbind(1, rmvn(n, p - 1, 0, SIGMA))

coef_list <- vector(mode = 'list', length = length(Ks) + 1)
names(coef_list) <- c("OLS", "K5", "K10", "Kn")

# simulation
pt <- proc.time()
registerDoParallel(cores = 4)
sim <- foreach(1:nsims, .combine = cbind) %dopar% {
  y <- X %*% beta_star + rnorm(n, 0, sigma_star)
  
  ols_fit <- ridge_coef(X, y, 0)
  coef_list[[1]] <- c(ols_fit$b1, ols_fit$b)

  coef_list[2:(length(Ks) + 1)] <- sapply(Ks, function(k) {
    rcv <- ridge_cv(X, y, lam.vec = lams, K = k)
    list(coefs = c(rcv$b1, rcv$b))
  })
  
  l1 <- sapply(coef_list, function(b) loss1(beta_star, b))
  l2 <- sapply(coef_list, function(b) loss2(X, beta_star, b))
  list(l1, l2)
}
sim_loss <- lapply(1:nrow(sim),
                   function(i) sapply(sim[i,], function(s) s))
names(sim_loss) <- c("Loss 1", "Loss 2")

sim_means <- t(sapply(sim_loss, function(s) rowMeans(s)))
sim_se <- t(sapply(sim_loss,
                   function(s) apply(s, 1, function(x) sd(x)/sqrt(length(x)))))
proc.time() - pt

# report results
round(sim_means, 4)
round(sim_se, 4)


  





































