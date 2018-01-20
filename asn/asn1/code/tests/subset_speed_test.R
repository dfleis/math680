set.seed(124)

# parameters
n <- 1e2
p <- 1e3

# build data
X <- matrix(rnorm(n * p), nrow = n, ncol = p)

# test
mb <- microbenchmark(X[,-1], X[,1:p], X[,1:ncol(X)], times = 1e2, unit = "us")
boxplot(mb)
