#==========
#
# Benchmarking different ways of solving X^T X beta = X^T y
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
f01 <- function(n, p, sigma = 1) { # base R 1
  data <- make_data(n, p, sigma)
  lm(data$y ~ data$X[,-1])$coef 
}
f02 <- function(n, p, sigma = 1) { # base R 2
  data <- make_data(n, p, sigma)
  coef(lm(data$y ~ data$X[,-1]))
}  
f03 <- function(n, p, sigma = 1) { # naive 1
  data <- make_data(n, p, sigma)
  solve(t(data$X) %*% data$X) %*% t(data$X) %*% data$y
}  
f04 <- function(n, p, sigma = 1) { # naive 2
  data <- make_data(n, p, sigma)
  solve(crossprod(data$X)) %*% crossprod(data$X, data$y)
}  
f05 <- function(n, p, sigma = 1) { # naive 3
  data <- make_data(n, p, sigma)
  solve(crossprod(data$X), crossprod(data$X, data$y))
} 
f06 <- function(n, p, sigma = 1) { # naive 4
  data <- make_data(n, p, sigma)
  solve(.Internal(crossprod(data$X, data$X)),
        .Internal(crossprod(data$X, data$y)))
} 
f07 <- function(n, p, sigma = 1) { # svd 1
  data <- make_data(n, p, sigma)
  Xsvd <- svd(data$X)
  U <- Xsvd$u; d <- Xsvd$d; V <- Xsvd$v
  V %*% diag(1/d) %*% t(U) %*% data$y
}
f08 <- function(n, p, sigma = 1) { # svd 2
  data <- make_data(n, p, sigma)
  Xsvd <- svd(data$X)
  U <- Xsvd$u; d <- Xsvd$d; V <- Xsvd$v
  V %*% (diag(1/d) %*% .Internal(crossprod(U, data$y)))
}


#===== benchmark =====#
set.seed(124)
n <- 50
p <- 50

pt <- proc.time()
mb <- microbenchmark(f01(n, p), f02(n, p), f03(n, p), f04(n, p), f05(n, p),
                     f06(n, p), f07(n, p), f08(n, p), times = 250, unit = "us")
tm <- proc.time() - pt
boxplot(mb, outline = F, xaxt = 'n', 
        main = paste0("(n, p) = ", "(", n, ", ", p, ")"),
        xlab = "Method")
axis(1, at = 1:nrow(summary(mb)), labels = 1:nrow(summary(mb)))
mb
tm
























