#==========
#
# R file containing all functions for use in Assignment 1
#
#==========

#===== Q2 =====#

ridge_cv_full <- function(X, y, lam.vec, K) {
  # perform K-fold cross-validation on the ridge regression 
  # estimation problem over tuning parameters given in lam.vec
  n <- nrow(X); p <- ncol(X)
  
  Xm1 <- X[,-1] # remove leading column of 1's marking the intercept
  ytilde <- y - mean(y) # center response
  xbar <- colMeans(Xm1) # find predictor means
  Xtilde <- Xm1 - tcrossprod(rep(1, nrow(Xm1)), xbar) # center each predictor according to its mean
  
  cv.error <- sapply(lam.vec, function(l) {
    # groups to cross-validate over
    folds <- cut(1:n, breaks = K, labels = F)
    # get indices of training subset
    train_idxs <- lapply(1:K, function(i) !(folds %in% i))
    
    cv_err <- sapply(train_idxs, function(tis) {
      # train our model, extract fitted coefficients
      Xsub <- Xtilde[tis,]
      ysub <- ytilde[tis]
      
      # compute the SVD on the centered design matrix
      Xsub_svd <- svd(Xsub)
      U <- Xsub_svd$u
      d <- Xsub_svd$d
      V <- Xsub_svd$v
      
      # compute the inverse (D^T D + lambda I_{p-1})^{-1} D^T
      Dstar <- diag(d/(d^2 + l))
      
      # compute ridge coefficients
      b <- V %*% (Dstar %*% crossprod(U, ysub)) # slopes
      b1 <- mean(y) - crossprod(xbar, b) # intercept
      b_train  <- c(b1, b)
        
      # find observations needed for testing fits
      test_idx <- !((1:n) %in% tis)
      # fit data
      yhat <- X[test_idx,] %*% b_train
      # compute test error
      sum((y[test_idx] - yhat)^2)
    })
    # weighted average (according to group size, some groups may have
    # +/- 1 member depending on whether sizes divided unevenly) of
    # cross validation error for a fixed lambda
    sum((cv_err * table(folds)))/n
  })
  
  # extract best tuning parameter and corresponding coefficient estimates
  best.lam <- lam.vec[cv.error == min(cv.error)]
  best.fit <- ridge_coef(X, y, best.lam)
  b1 <- best.fit$b1
  b <- best.fit$b
  
  list(b1 = b1, b = b, best.lam = best.lam, cv.error = cv.error)
}
