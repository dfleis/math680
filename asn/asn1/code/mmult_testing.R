#==========
#
# Testing the speed of different matrix multiplication routines
# for the matrix product M X^T X M for X (n * p) dimensional and 
# M (p * p) dimensional
# 
#==========
library(microbenchmark)
set.seed(124)

#===== Functions =====#
# define our matrix multiplication functions
f1 <- function(x, m) m %*% t(x) %*% x %*% t(m)
f2 <- function(x, m) m %*% t(x) %*% (x %*% t(m))
f3 <- function(x, m) m %*% (t(x) %*% (x %*% t(m)))
f4 <- function(x, m) (m %*% t(x)) %*% (x %*% t(m))
f5 <- function(x, m) m %*% (t(x) %*% x) %*% t(m)
f6 <- function(x, m) m %*% ((t(x) %*% x) %*% t(m))
f7 <- function(x, m) m %*% t(x) %*% tcrossprod(x, m)
f8 <- function(x, m) m %*% crossprod(x, tcrossprod(x, m))
f9 <- function(x, m) tcrossprod(m, x) %*% tcrossprod(x, m)
f10 <- function(x, m) m %*% crossprod(x) %*% t(m)
f11 <- function(x, m) m %*% tcrossprod(crossprod(x), m)

#===== Large n case, n >> p =====#
# define parameters
n <- 100
p <- 10

# generate data
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
M <- matrix(rnorm(p^2), ncol = p, nrow = p)

# compute tests
mb1 <- microbenchmark(f1(X, M), f2(X, M), f3(X, M), 
                      f4(X, M), f5(X, M), f6(X, M), 
                      f7(X, M), f8(X, M), f9(X, M),
                      f10(X, M), f11(X, M), times = 1e2, unit = "us")

#===== Large p case, p >> n =====#
# define parameters
n <- 10
p <- 100

# generate data
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
M <- matrix(rnorm(p^2), ncol = p, nrow = p)

# compute tests
mb2 <- microbenchmark(f1(X, M), f2(X, M), f3(X, M), 
                      f4(X, M), f5(X, M), f6(X, M), 
                      f7(X, M), f8(X, M), f9(X, M),
                      f10(X, M), f11(X, M), times = 1e2, unit = "us")

#===== n ~ p case =====#
# define parameters
n <- 32
p <- 32

# generate data
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
M <- matrix(rnorm(p^2), ncol = p, nrow = p)

# compute tests
mb3 <- microbenchmark(f1(X, M), f2(X, M), f3(X, M), 
                      f4(X, M), f5(X, M), f6(X, M), 
                      f7(X, M), f8(X, M), f9(X, M),
                      f10(X, M), f11(X, M), times = 1e2, unit = "us")

#===== Plots =====#

boxplot(mb1, outline = F, main = "n >> p")
boxplot(mb2, outline = F, main = "p >> n")
boxplot(mb3, outline = F, main = "n ~ p")










