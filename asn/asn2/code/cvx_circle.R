#===== libraries =====#
library(CVXR)
library(fields)

#===== load data =====#
circle <- as.matrix(read.csv("data/circle.csv", header = F))


nr <- nrow(circle)
nc <- ncol(circle)
n <- length(circle)

#==== define fused penalty ====#
fused_lasso_2d <- function(theta, lambda = 0) {
  nr <- nrow(theta); nc <- ncol(theta)
  S <- theta[1:(nr - 1),] - theta[2:nr,] # SOUTH
  N <- theta[2:nr,] - theta[1:(nr - 1),] # NORTH
  E <- theta[,1:(nc - 1)] - theta[,2:nc] # EAST
  W <- theta[,2:nc] - theta[,1:(nc - 1)] # WEST
  lambda * (sum(abs(S)) + sum(abs(N)) + sum(abs(E)) + sum(abs(W)))
}

#==== tests =====#
lambda_vals <- c(1, 2)
theta_vals <- vector(mode = 'list', length = length(lambda_vals))
theta <- Variable(nr, nc)
loss <- sum(0.5 * (circle - theta)^2)

for (i in 1:length(lambda_vals)) {
  lambda <- lambda_vals[i]
  obj <- loss + fused_lasso_2d(theta, lambda)
  prob <- Problem(Minimize(obj))
  res <- solve(prob)
  theta_vals[[i]] <- res$getValue(theta)
}

#===== plots =====#
cols <- colorRampPalette(c("black", "white"))
image.plot(circle, col = cols(10))
image.plot(theta_vals[[2]])





