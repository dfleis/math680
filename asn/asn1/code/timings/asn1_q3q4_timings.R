#==========
#
# Test speed of different solutions to Q2
#
#==========
#===== libraries =====#
library(microbenchmark)
source("code/asn1_functions.R")

#===== functions =====#
ridge_cv1 <- function(X, y, lam.vec, K) {
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
ridge_cv2 <- function(X, y, lam.vec, K) {
  # perform K-fold cross-validation on the ridge regression 
  # estimation problem over tuning parameters given in lam.vec
  n <- nrow(X); p <- ncol(X); L <- length(lam.vec)
  
  cv.error <- vector(mode = 'numeric', length = L)
  for (i in 1:L) 
    cv.error[i] <- ridge_cv_lam(X, y, lam.vec[i], K)
  
  # extract best tuning parameter and corresponding coefficient estimates
  best.lam <- lam.vec[cv.error == min(cv.error)]
  best.fit <- ridge_fit(X, y, best.lam)
  b1 <- best.fit$coef[1]
  b <- best.fit$coef[-1]
  
  list(b1 = b1, b = b, best.lam = best.lam, cv.error = cv.error)
}
ridge_cv3 <- function(X, y, lam.vec, K) {
  # perform K-fold cross-validation on the ridge regression 
  # estimation problem over tuning parameters given in lam.vec
  n <- nrow(X); p <- ncol(X); L <- length(lam.vec)
  
  cv.error <- sapply(1:L, function(i) ridge_cv_lam(X, y, lam.vec[i], K))
  
  # extract best tuning parameter and corresponding coefficient estimates
  best.lam <- lam.vec[cv.error == min(cv.error)]
  best.fit <- ridge_coef(X, y, best.lam)
  b1 <- best.fit$b1
  b <- best.fit$b
  
  list(b1 = b1, b = b, best.lam = best.lam, cv.error = cv.error)
}


#===== set parameters/generate data =====#
set.seed(124)
n <- 10
p <- 50
lams <- 10^seq(-8, 8, 0.5)
K <- 5

X <- cbind(1, matrix(rnorm(n * (p - 1)), nrow = n))
beta <- rnorm(p)
y <- X %*% beta + rnorm(n)

pt <- proc.time()
mb <- microbenchmark(ridge_cv1(X, y, lams, K), 
                     ridge_cv2(X, y, lams, K),
                     ridge_cv3(X, y, lams, K), times = 100, unit = "ms")
boxplot(mb, outline = F)
proc.time() - pt
mb

