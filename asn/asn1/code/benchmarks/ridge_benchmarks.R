#==========
#
# Benchmarking different ways of solving (X^T X + lam * I) beta = X^T y
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
f01 <- function(n, p, lam = 1, sigma = 1) { # naive 1
  data <- make_data(n, p, sigma)
  solve(t(data$X) %*% data$X + lam * diag(p)) %*% t(data$X) %*% data$y
}  
f02 <- function(n, p, lam = 1, sigma = 1) { # naive 2
  data <- make_data(n, p, sigma)
  solve(crossprod(data$X) + lam * diag(p)) %*% crossprod(data$X, data$y)
}  
f03 <- function(n, p, lam = 1, sigma = 1) { # naive 3
  data <- make_data(n, p, sigma)
  solve(crossprod(data$X) + lam * diag(p), crossprod(data$X, data$y))
} 
f04 <- function(n, p, lam = 1, sigma = 1) { # naive 4
  data <- make_data(n, p, sigma)
  solve(.Internal(crossprod(data$X, data$X)) + lam * diag(p),
        .Internal(crossprod(data$X, data$y)))
} 
f05 <- function(n, p, lam = 1, sigma = 1) { # svd 1
  data <- make_data(n, p, sigma)
  Xsvd <- svd(data$X)
  U <- Xsvd$u; d <- Xsvd$d; V <- Xsvd$v
  V %*% diag(d/(d^2 + lam)) %*% t(U) %*% data$y
}
f06 <- function(n, p, lam = 1, sigma = 1) { # svd 2
  data <- make_data(n, p, sigma)
  Xsvd <- svd(data$X)
  U <- Xsvd$u; d <- Xsvd$d; V <- Xsvd$v
  V %*% (diag(d/(d^2 + lam)) %*% .Internal(crossprod(U, data$y)))
}
f07 <- function(n, p, sigma = 1) { # cholesky 1
  data <- make_data(n, p, sigma)
  C <- chol(crossprod(data$X) + lam * diag(p))
  chol2inv(C) %*% crossprod(data$X, data$y)
}
f08 <- function(n, p, sigma = 1) { # cholesky 2
  data <- make_data(n, p, sigma)
  C <- chol(.Internal(crossprod(data$X, data$X)) + lam * diag(p))
  chol2inv(C) %*% .Internal(crossprod(data$X, data$y))
}

#===== benchmark =====#
set.seed(124)
n <- 50
p <- 10

pt <- proc.time()
mb <- microbenchmark(f01(n, p), f02(n, p), f03(n, p), f04(n, p), f05(n, p),
                     f06(n, p), f07(n, p), f08(n, p), 
                     times = 250, unit = "us")
tm <- proc.time() - pt
boxplot(mb, outline = F, xaxt = 'n', 
        main = paste0("(n, p) = ", "(", n, ", ", p, ")"),
        xlab = "Method")
axis(1, at = 1:nrow(summary(mb)), labels = 1:nrow(summary(mb)))
mb
tm
























