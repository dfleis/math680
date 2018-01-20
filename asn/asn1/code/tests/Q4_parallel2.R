#===== libraries =====#
library(doParallel)

#===== functions =====#
ridge_coef <- function(X, y, lam) {
  Xm1 <- X[,-1] # remove leading column of 1's marking the intercept
  
  ytilde <- y - mean(y) # center response
  xbar <- colMeans(Xm1) # find predictor means
  Xtilde <- sweep(Xm1, 2, xbar) # center each predictor according to its mean
  
  # compute the SVD on the centered design matrix
  Xtilde_svd <- svd(Xtilde)
  U <- Xtilde_svd$u
  d <- Xtilde_svd$d
  V <- Xtilde_svd$v
  
  # compute the inverse (D^T D + lambda I_{p-1})^{-1} D^T
  Dstar <- diag(d/(d^2 + lam))
  
  # compute ridge coefficients
  b <- V %*% (Dstar %*% crossprod(U, ytilde)) # slopes
  b1 <- mean(y) - crossprod(xbar, b) # intercept
  list(b1 = b1, b = b)
}
ridge_fit <- function(X, y, lam) {
  # fully fit a ridge regression model given predictors, response, and penalty
  
  b <- unlist(ridge_coef(X, y, lam)) # extract coefficient estimates
  yhat <- X %*% b # fit a response estimate given fitted coefficients
  res <- sum((y - yhat)^2) # find prediction error
  
  list(X = X, y = y, lam = lam, coef = b, fit = yhat, res = res)
}
ridge_cv_lam <- function(X, y, lam, K) {
  # Helper function for ridge_cv()
  # perform K-fold cross-validation on the ridge regression 
  # estimation problem over a single tuning parameter lam
  n <- nrow(X)
  
  if (K > n) { 
    stop(paste0("K > ", n, "."))
  } else if (K < 2) {
    stop("K < 2.")
  }
  
  # groups to cross-validate over
  folds <- cut(1:nrow(X), breaks = K, labels = F)
  # get indices of training subset
  train_idxs <- lapply(1:K, function(i) !(folds %in% i))
  
  cv_err <- sapply(train_idxs, function(tis) {
    # train our model
    train_fit <- ridge_fit(X[tis,], y[tis], lam)
    
    # find observations needed for testing fits
    test_idx <- !((1:n) %in% tis)
    
    # extract fitted coefficients
    b <- train_fit$coef
    # fit data
    yhat <- X[test_idx,] %*% b
    # compute test error
    sum((y[test_idx] - yhat)^2)
  })
  # weighted average (according to group size, some groups may have
  # +/- 1 member depending on whether sizes divided unevenly) of
  # cross validation error for a fixed lambda
  sum((cv_err * table(folds)))/n
}
ridge_cv <- function(X, y, lam.vec, K) {
  # perform K-fold cross-validation on the ridge regression 
  # estimation problem over tuning parameters given in lam.vec
  n <- nrow(X); p <- ncol(X); L <- length(lam.vec)
  
  cv.error <- sapply(1:L, function(i) ridge_cv_lam(X, y, lam.vec[i], K))
  
  # extract best tuning parameter and corresponding coefficient estimates
  best.lam <- lam.vec[cv.error == min(cv.error)]
  best.fit <- ridge_fit(X, y, best.lam)
  b1 <- best.fit$coef[1]
  b <- best.fit$coef[-1]
  
  list(b1 = b1, b = b, best.lam = best.lam, cv.error = cv.error)
}
rmvn <- function(n, p, mu = 0, S = diag(p)) {
  # generates n (potentially correlated) p-dimensional normal deviates
  # given mean vector mu and variance-covariance matrix S
  # NOTE: S must be a positive-semidefinite matrix
  Z <- matrix(rnorm(n * p), nrow = n, ncol = p) # generate iid normal deviates
  C <- chol(S)
  mu + Z %*% C # compute our correlated deviates
}

loss1 <- function(beta, b) sum((b - beta)^2)
loss2 <- function(X, beta, b) sum((X %*% (beta - b))^2)

#===== global parameters =====#
set.seed(125)

nsims <- 3
n <- 100
Ks <- c(5, 10, n)
lams <- 10^seq(-8, 8, 0.5)
sigma_star <- sqrt(1/2)

#===== (a) =====#
# set parameters
p <- 50
theta <- 0.5

# generate data
beta_star <- rnorm(p, 0, sigma_star)
SIGMA <- outer(1:(p - 1), 1:(p - 1), FUN = function(a, b) theta^abs(a - b))
X <- cbind(1, rmvn(n, p - 1, 0, SIGMA))

pt <- proc.time()
registerDoParallel(cores = 4)

sim <- foreach(1:nsims, .combine = cbind) %dopar% {
  y <- X %*% beta_star + rnorm(n, 0, sigma_star)  
  
  ols_fit <- ridge_fit(X, y, 0)
  k5_fit <- ridge_cv(X, y, lam.vec = lams, K = 5)
  k10_fit <- ridge_cv(X, y, lam.vec = lams, K = 10)
  kn_fit <- ridge_cv(X, y, lam.vec = lams, K = n)
  
  coef_list <- list(OLS = ols_fit$coef, 
                    k5  = c(k5_fit$b1, k5_fit$b), 
                    k10 = c(k10_fit$b1, k10_fit$b),
                    kn  = c(kn_fit$b1, kn_fit$b))
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
sim_means
