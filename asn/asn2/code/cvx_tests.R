#===== libraries =====#
library(CVXR)
library(fields)

#===== generate data =====#
set.seed(124)

# Problem data
n <- 20
m <- 1000
DENSITY <- 0.25    # Fraction of non-zero beta
beta_true <- matrix(rnorm(n), ncol = 1)
idxs <- sample.int(n, size = floor((1 - DENSITY) * n), replace = FALSE)
beta_true[idxs] <- 0

sigma <- 45
X <- matrix(rnorm(m * n, sd = 5), nrow = m, ncol = n)
eps <- matrix(rnorm(m, sd = sigma), ncol = 1)
y <- X %*% beta_true + eps

elastic_reg <- function(beta, lambda = 0, alpha = 0) {
  #ridge <- (1 - alpha) * sum(beta^2)
  #lasso <- alpha * p_norm(beta, 1)
  #lambda * (lasso + ridge)
  lambda * p_norm(beta, 1.5)
}

TRIALS <- 20
beta_vals <- matrix(0, nrow = n, ncol = TRIALS)
lambda_vals <- 10^seq(0, log10(50), length.out = TRIALS)
beta <- Variable(n)  
loss <- sum((y - X %*% beta)^2)/(2*m)

## Elastic-net regression
alpha <- 1
pt <- proc.time()
for(i in 1:TRIALS) {
  lambda <- lambda_vals[i]
  obj <- loss + elastic_reg(beta, lambda, alpha)
  prob <- Problem(Minimize(obj))
  result <- solve(prob)
  beta_vals[,i] <- result$getValue(beta)
}
proc.time() - pt


plot(NA, type = "n", main = "Regularization Path for Elastic-net Regression",
     log = 'x',
     xlab = expression(lambda), ylab = expression(beta),
     ylim = range(beta_vals), xlim = range(lambda_vals))
matlines(lambda_vals, t(beta_vals))
