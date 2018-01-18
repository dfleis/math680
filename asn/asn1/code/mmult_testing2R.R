library(microbenchmark)
set.seed(124)

#===== functions =====#
# define functions to benchmark
f1 <- function(X, lam) {
  I <- diag(ncol(X))
  
  solve(I + lam * solve(t(X) %*% X)) %*% 
    solve(t(X) %*% X) %*% 
    solve(I + lam * solve(t(X) %*% X))
}
f2 <- function(X, lam) {
  X_svd <- svd(X)
  d <- X_svd$d
  V <- X_svd$v
  
  Dxx <- diag(1/(d + lam/d)^2)
  V %*% tcrossprod(Dxx, V)
}
f3 <- function(X, lam) {
  X_svd <- svd(X)
  d <- X_svd$d
  V <- X_svd$v
  
  V %*% tcrossprod(diag(d^2/(d^2 + lam)^2), V)
}

#===== data =====#
# parameters
n <- 1e1
p <- 1e2
lam <- 1

# generate data
X <- matrix(rnorm(n * p), nrow = n, ncol = p)

#===== work =====#

all.equal(f1(X, lam), f2(X, lam))
mb <- microbenchmark(f1(X, lam), f2(X, lam), f3(X, lam), times = 1e2, unit = "us")
boxplot(mb, outline = F)
mb


X_svd <- svd(X)
d <- X_svd$d
V <- X_svd$v
U <- X_svd$u

Dxx <- diag(1/(d + lam/d)^2)
V %*% Dxx %*% t(V)

x <- 
y <- f1(X, lam)
x - y



