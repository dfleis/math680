#===== libraries =====#
library(CVXR)

#===== load data =====#
lenna <- as.matrix(read.csv("data/lenna_64.csv", header = F))
lenna <- lenna[seq(4, 64, 4), seq(4, 64, 4)]

nr <- nrow(lenna)
nc <- ncol(lenna)
n <- length(lenna)

#==== define fused penalty ====#
fused_lasso_2d <- function(theta, lambda = 0) {
  nr <- nrow(theta); nc <- ncol(theta)
  S <- theta[1:(nr - 1),] - theta[2:nr,] # SOUTH
  N <- theta[2:nr,] - theta[1:(nr - 1),] # NORTH
  E <- theta[,1:(nc - 1)] - theta[,2:nc] # EAST
  W <- theta[,2:nc] - theta[,1:(nc - 1)] # WEST
  #lambda * p_norm(c(S, N, E, W), 1)
  lambda * (p_norm(S, 1) + p_norm(N, 1) + p_norm(E, 1) + p_norm(W, 1))
  #lambda * (sum(abs(S)) + sum(abs(N)) + sum(abs(E)) + sum(abs(W)))
}

#==== tests =====#
lambda_vals <- 10^(-(0:8)/4)
theta_vals <- vector(mode = 'list', length = length(lambda_vals))
obj_val <- vector(mode = 'numeric', length = length(lambda_vals))
theta <- Variable(nr, nc)
loss <- sum(0.5 * (lenna - theta)^2)

for (i in 1:length(lambda_vals)) {
  pt <- proc.time()
  lambda <- lambda_vals[i]
  obj <- loss + fused_lasso_2d(theta, lambda)
  prob <- Problem(Minimize(obj))
  res <- solve(prob)
  theta_vals[[i]] <- res$getValue(theta)
  obj_val[i] <- res$value
  print(proc.time() - pt)
}

#===== plots =====#
cols <- colorRampPalette(c("black", "white"))
lenna2 <- t(apply(lenna, 2, rev))

image(lenna2, col = cols(10), xaxt = 'n', yaxt = 'n')

k <- 9
theta_hat <- t(apply(theta_vals[[k]], 2, rev))
image(theta_hat, col = cols(10), xaxt = 'n', yaxt = 'n',
      main = substitute(paste(lambda, " = ", lam), list(lam = lambda_vals[k])))
hist(theta_hat, breaks = 1e2, xlim = c(0, 1),
     main = substitute(paste(lambda, " = ", lam), list(lam = lambda_vals[k])))


