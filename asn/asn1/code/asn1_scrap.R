set.seed(124)

n <- 100
p <- 2e4

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
beta <- rnorm(p)
eps <- rnorm(n)

y <- X %*% beta + eps

ridge_coef <- function(X, y, lam) {
  ytilde <- y - mean(y)
  xbar <- colMeans(X)
  Xtilde <- sweep(X, 2, xbar)
  
  Xtilde_svd <- svd(Xtilde)
  U <- Xtilde_svd$u
  d <- Xtilde_svd$d
  V <- Xtilde_svd$v
  
  Dstar <- diag(d/(d^2 + lam))
  
  b1 <- mean(y) - crossprod(xbar, b)
  b <- V %*% (Dstar %*% crossprod(U, ytilde))
  return (list(b1 = b1, b = b))
}

f1 <- function() V %*% Dstar %*% t(U) %*% ytilde
f2 <- function() V %*% Dstar %*% (t(U) %*% ytilde)
f3 <- function() V %*% (Dstar %*% (t(U) %*% ytilde))
f4 <- function() V %*% (Dstar %*% crossprod(U, ytilde))
f5 <- function() V %*% crossprod(Dstar, crossprod(U, ytilde))

microbenchmark(f1(), f2(), f3(), f4(), f5(), times = 1e2, unit = "us")




