set.seed(124)

# parameters
n <- 2
p <- 5
lam1 <- 1
lam2 <- 1

# build matrices
X <- matrix(rnorm(n * (p - 1)), nrow = n, ncol = p - 1)
O <- cbind(rep(0, p - 2), diag(p - 2)) # off-diag matrix along the upper diagonal
J <- -1 * cbind(diag(p - 2), rep(0, p - 2)) # diag (p - 2)*(p - 1) matrix

A <- J + O
A2 <- t(A) %*% A

q <- qr.Q(qr(A))
r <- qr.R(qr(A))

r

chol(A2)
lu.decomposition(A2)

eig <- eigen(A2)
P <- eig$vectors
L <- diag(eig$values)
all.equal(P %*% L %*% t(P), A2)




set.seed(124)

# parameters
n <- 10
p <- 10
lam1 <- 1
lam2 <- 1

# build data
X <- cbind(rep(1, n), matrix(rnorm(n * (p - 1)), nrow = n, ncol = p - 1))

# compute
Xm1 <- X[,-1]
xbar <- colMeans(Xm1)
Xtilde <- sweep(Xm1, 2, xbar)

Xtilde_svd <- svd(Xtilde)
V <- Xtilde_svd$v
d <- Xtilde_svd$d
D <- diag(d)

I <- diag(p - 1)
O <- cbind(rep(0, p - 2), diag(p - 2)) # off-diag matrix along the upper diagonal
J <- -1 * cbind(diag(p - 2), rep(0, p - 2)) # diag (p - 2)*(p - 1) matrix
A <- J + O

D2lam1 <- diag(d^2 + lam1)

f1 <- function() solve(crossprod(Xtilde) + lam1 * I + lam2 * crossprod(A))
f2 <- function() solve(V %*% tcrossprod(D2lam1, V) + lam2 * crossprod(A))

mb <- microbenchmark(f1(), f2(), times = 1e2, unit = "us")
mb

Z <- t(V) %*% crossprod(A) %*% V
all.equal(V %*% Z %*% t(V), crossprod(A))
Z
image(Z)
crossprod(A)

Z <- lam1 * I + lam2 * crossprod(A)
qr.Q(qr(Z))

























