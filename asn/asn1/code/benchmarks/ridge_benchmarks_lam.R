#==========
#
# Benchmarking different ways of solving (X^T X + lam * I) beta = X^T y
# over a vector of lambdas
#
#==========
library(microbenchmark)

#===== functions =====#
make_data <- function(n, p, sigma = 1) {
  X <- cbind(1, matrix(rnorm(n * (p - 1)), nrow = n))
  beta <- rnorm(p)
  y <- X %*% beta + rnorm(n, 0, sigma)
  list(X = X, beta = beta, y = y)
}
f01 <- function(n, p, lams, sigma = 1) { # naive 1
  data <- make_data(n, p, sigma)
  lapply(lams, function(lam) 
    solve(t(data$X) %*% data$X + lam * diag(p)) %*% t(data$X) %*% data$y
  )
}  
f02 <- function(n, p, lams, sigma = 1) { # naive 2
  data <- make_data(n, p, sigma)
  Xty <- crossprod(data$X, data$y)
  XtX <- crossprod(data$X)
  I <- diag(p)
  lapply(lams, function(lam) 
    solve(XtX + lam * I, Xty)
  )
} 
f03 <- function(n, p, lams, sigma = 1) { # svd
  data <- make_data(n, p, sigma)
  Xsvd <- svd(data$X)
  U <- Xsvd$u; d <- Xsvd$d; V <- Xsvd$v
  Uty <- crossprod(U, data$y)
  d2 <- d^2
  lapply(lams, function(lam)
    V %*% (diag(d/(d2 + lam)) %*% Uty)
  )
}
f04 <- function(n, p, lams, sigma = 1) { # cholesky 1
  data <- make_data(n, p, sigma)
  XtX <- crossprod(data$X)
  Xty <- crossprod(data$X, data$y)
  I <- diag(p)
  lapply(lams, function(lam) {
    C <- chol(XtX + lam * I)
    chol2inv(C) %*% Xty
  })
}

#===== benchmark =====#
set.seed(124)
n <- 50
p <- 5
lams <- seq(-8, 8, 0.5)

pt <- proc.time()
mb <- microbenchmark(f01(n, p, lams), f02(n, p, lams), f03(n, p, lams), f04(n, p, lams),
                     times = 250, unit = "us")
tm <- proc.time() - pt
boxplot(mb, outline = F, xaxt = 'n', 
        main = paste0("(n, p) = ", "(", n, ", ", p, ")"),
        xlab = "Method")
axis(1, at = 1:nrow(summary(mb)), labels = 1:nrow(summary(mb)))
mb
tm
























