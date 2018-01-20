ridge_coef <- function(X, y, lam) {
  Xm1 <- X[,-1] # remove leading column of 1's marking the intercept
  
  ytilde <- y - mean(y) # center response
  xbar <- colMeans(Xm1) # find predictor means
  Xtilde <- sweep(Xm1, 2, xbar) # center each predictor according to its mean
  
  # compute the SVD on the centered design matrix
  Xtilde_svd <- svd(Xtilde)
  U <- Xtilde_svd$u
  d <- Xtilde_svd$d
  V <- Xtilde_svd$v
  
  # compute the inverse (D^T D + lambda I_{p-1})^{-1} D^T
  Dstar <- diag(d/(d^2 + lam))
  
  # compute ridge coefficients
  b <- V %*% (Dstar %*% crossprod(U, ytilde)) # slopes
  b1 <- mean(y) - crossprod(xbar, b) # intercept
  return (list(b1 = b1, b = b))
}
ridge_coef_params <- function(X, lam, beta, sigma) {
  n <- nrow(X); p <- ncol(X)
  betam1 <- beta[-1] # remove intercept term
  Xm1 <- matrix(X[,-1]) # remove leading column of 1's in our design matrix
  
  xbar <- colMeans(Xm1) # find prector means
  Xtilde <- sweep(Xm1, 2, xbar) # center each predictor according to its mean
  
  if (n >= p) {
    I <- diag(p - 1)
    inv <- solve(crossprod(Xtilde) + lam * I)
    
    b <- solve(crossprod(Xtilde) + lam * I) %*% (crossprod(Xtilde) %*% betam1)
    vcv <- sigma^2 * inv %*% crossprod(Xtilde) %*% inv
    return (list(b = b, vcv = vcv)) 
    
  } else {
    # compute SVD on the centered design matrix
    Xtilde_svd <- svd(Xtilde)
    d <- Xtilde_svd$d
    V <- Xtilde_svd$v
    
    Dstar <- diag(d^2/(d^2 + lam))
    Dstar2 <- diag(d^2/(d^2 + lam)^2)
    
    b <- V %*% (Dstar %*% crossprod(V, betam1))
    vcv <- V %*% tcrossprod(Dstar2, V)
    return (list(b = b, vcv = vcv))
  }
}
ridge_fit <- function(X, y, lam) {
  b <- unlist(ridge_coef(X, y, lam)) # extract coefficient estimates
  yhat <- X %*% b # fit a response estimate given fitted coefficients
  res <- sum((y - yhat)^2) # find prediction error
  
  return (list(X = X, y = y, lam = lam, coef = b, fit = yhat, res = res))
}
ridge_cv <- function(X, y, lam.vec, K) {
  # perform K-fold cross-validation on the ridge regression 
  # estimation problem over tuning parameters given in lam.vec
  n <- nrow(X); p <- ncol(X); L <- length(lam.vec)
  
  cv.error <- vector(mode = "numeric", length = L)
  for (i in 1:L) {
    cv.error[i] <- ridge_cv_lam(X, y, lam.vec[i], K)  
  }

  best.lam <- lam.vec[which(cv.error == min(cv.error))]
  best.fit <- ridge_fit(X, y, best.lam)
  b1 <- best.fit$coef[1]
  b <- best.fit$coef[-1]
  
  return (list(b1 = b1, b = b, best.lam = best.lam, cv.error = cv.error))
}
ridge_cv_lam <- function(X, y, lam, K) {
  # Helper function for ridge_cv()
  # perform K-fold cross-validation on the ridge regression 
  # estimation problem over a single tuning parameter lam
  
  if (K > n) { 
    stop(paste0("K > ", n, "."))
  } else if (K < 2) {
    stop("K < 2.")
  }
  
  # groups to cross-validate over
  folds <- cut(1:nrow(X), breaks = K, labels = F)
  train_grps <- lapply(1:K, function(i) which(!(1:K %in% i)))
  # get indices of our training subsets
  train_idxs <- lapply(train_grps, function(tgs) which(folds %in% tgs))
  
  cv_err <- sapply(train_idxs, function(tis) {
    # train our model
    train_fit <- ridge_fit(X[tis,], y[tis], lam)
    
    # find observations needed for testing fits
    test_idx <- which(!((1:n) %in% tis))
    
    # extract fitted coefficients
    b <- train_fit$coef
    # fit data
    yhat <- X[test_idx,] %*% b
    # compute test error
    sum((y[test_idx] - yhat)^2)
  })
  # weighted average (according to group size, some groups may have
  # +/- 1 member depending on whether sizes divided unevenly) of
  # cross validation error for a fixed lambda
  sum((cv_err * table(folds)))/n
}


set.seed(124)
n <- 1e3
p <- 3
beta <- 1:p
sigma <- 1
K <- 5
lams <- seq(0.001, 1, length.out = 100)

X <- cbind(1, matrix(rnorm(n * (p - 1)), nrow = n))
y <- X %*% beta + rnorm(n, 0, sigma)

pt <- proc.time()
cv <- ridge_cv(X, y, lams, K)
proc.time() - pt

plot(cv$cv.error ~ lams, type = 'l')
cv$best.lam
cv$b1
cv$b

