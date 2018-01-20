library(microbenchmark)
set.seed(124)
n <- 10
p <- 11
lam <- 2.5

X <- matrix(rnorm(n * p), nrow = n)
I <- diag(p)

f1 <- function() {
  solve(t(X) %*% X + lam * I) %*% t(X) %*% X %*% solve(t(X) %*% X + lam * I)
}
f2 <- function() {
  inv <- solve(t(X) %*% X + lam * I)
  inv %*% t(X) %*% X %*% inv
}
f3 <- function() {
  inv <- solve(t(X) %*% X + lam * I)
  (inv %*% t(X)) %*% (X %*% inv)
}
f4 <- function() {
  inv <- solve(crossprod(X) + lam * I)
  inv %*% crossprod(X) %*% inv
}
f5 <- function() {
  X_svd <- svd(X)
  V <- X_svd$v
  d <- X_svd$d
  Dstar <- diag(d^2/(d^2 + lam)^2) # solve(t(D) %*% D + lam * I)
  V %*% Dstar %*% t(V)
}
f6 <- function() {
  X_svd <- svd(X)
  V <- X_svd$v
  d <- X_svd$d
  Dstar <- diag(d^2/(d^2 + lam)^2) # solve(t(D) %*% D + lam * I)
  V %*% tcrossprod(Dstar, V)
}

pt <- proc.time()
mb <- microbenchmark(f1(), f2(), f3(), f4(), 
                     f5(), f6(), times = 1e2, unit = "us")
proc.time() - pt
boxplot(mb, outline = F)
mb






