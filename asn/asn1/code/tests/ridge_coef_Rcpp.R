library(Rcpp)
sourceCpp("code/ridge_svd.cpp")

ridge_coef <- function(X, y, lam) {
  # Commented-out scaling parameters represent the transformations
  # used to make the output identical to that of 
  # coef(MASS::lm.ridge(y ~ X[,-1], lambda = lam))
  Xm1 <- X[,-1]; xbar <- colMeans(Xm1); ytilde <- y - mean(y)
  
  # center each predictor according to its mean
  # Xtilde <- (Xm1 - tcrossprod(rep(1, nrow(Xm1)), xbar)) 
  Xtilde <- scale(Xm1) * sqrt(n/(n - 1))

  # compute the SVD on the centered design matrix
  Xtilde_svd <- svd(Xtilde)
  U <- Xtilde_svd$u; d <- Xtilde_svd$d; V <- Xtilde_svd$v
  
  # compute the inverse (D^T D + lambda I_{p-1})^{-1} D^T
  Dstar <- diag(d/(d^2 + lam))
  
  # compute ridge coefficients
  b <- V %*% (Dstar %*% crossprod(U, ytilde)) * 
    1/attr(Xtilde, "scaled:scale") * sqrt(n/(n - 1))
  b1 <- mean(y) - crossprod(xbar, b) 
  list(b1 = b1, b = b)
}
ridge_coef_cpp <- function(X, y, lam) {
  # Commented-out scaling parameters represent the transformations
  # used to make the output identical to that of 
  # coef(MASS::lm.ridge(y ~ X[,-1], lambda = lam))
  Xm1 <- X[,-1]; xbar <- colMeans(Xm1); ytilde <- y - mean(y)
  
  # center each predictor according to its mean
  # Xtilde <- (Xm1 - tcrossprod(rep(1, nrow(Xm1)), xbar)) 
  Xtilde <- scale(Xm1) * sqrt(n/(n - 1))
  
  # compute the SVD on the centered design matrix
  Xtilde_svd <- svd_cpp(Xtilde)
  U <- Xtilde_svd$u; d <- Xtilde_svd$d; V <- Xtilde_svd$v
  
  # compute the inverse (D^T D + lambda I_{p-1})^{-1} D^T
  Dstar <- diag(d/(d^2 + lam))
  
  # compute ridge coefficients
  b <- V %*% (Dstar %*% crossprod(U, ytilde)) * 
    1/attr(Xtilde, "scaled:scale") * sqrt(n/(n - 1))
  b1 <- mean(y) - crossprod(xbar, b) 
  list(b1 = b1, b = b)
}

ridge_cv <- function(X, y, lam.vec, K) {
  # perform K-fold cross-validation on the ridge regression 
  # estimation problem over tuning parameters given in lam.vec
  n <- nrow(X); p <- ncol(X); L <- length(lam.vec)
  
  # groups to cross-validate over
  folds <- cut(1:n, breaks = K, labels = F)
  # get indices of training subset
  train_idxs <- lapply(1:K, function(i) !(folds %in% i))
  
  # preallocate empty data structure to store our CV errors 
  # for each lambda & fold
  cv_errs_lams_folds <- matrix(0, nrow = L, ncol = K)
  cv_errs_lams_folds <- sapply(train_idxs, function(trn_idx) {
    tst_idx <- !trn_idx # find test subset indices
    
    # subset data and remove intercept column in the design matrix
    Xm1_trn <- X[trn_idx, -1]; y_trn <- y[trn_idx]; n_trn <- nrow(Xm1_trn) 
    xbar_trn <- colMeans(Xm1_trn); ytilde_trn <- y_trn - mean(y_trn) 
    
    # center each predictor according to its mean
    # Xtilde_trn <- (Xm1_trn - tcrossprod(rep(1, nrow(Xm1_trn)), xbar_trn)) 
    Xtilde_trn <- scale(Xm1_trn) * sqrt(n_trn/(n_trn - 1))
    
    # compute the SVD on the centered design matrix
    # Note: This is done before computing our coefficients for each
    # lambda since the SVD will be identical across the vector of
    # tuning parameters
    Xtilde_trn_svd <- svd(Xtilde_trn)
    U <- Xtilde_trn_svd$u; d <- Xtilde_trn_svd$d; V <- Xtilde_trn_svd$v
    d2 <- d^2
    
    # preallocate empty data structure to store CV errors 
    # for each value of lambda
    cv_errs_lams <- vector(mode = 'numeric', length = L)
    cv_errs_lams <- sapply(lam.vec, function(lam) {
      # compute the inverse (D^T D + lambda I_{p-1})^{-1} D^T
      Dstar <- diag(d/(d2 + lam))
      # compute ridge coefficients
      b_trn <- V %*% (Dstar %*% crossprod(U, ytilde_trn)) * 
        1/attr(Xtilde_trn, "scaled:scale") * sqrt(n_trn/(n_trn - 1))
      b1_trn <- mean(y_trn) - crossprod(xbar_trn, b_trn)
      
      # fit test data
      yhat_tst <- X[tst_idx,] %*% c(b1_trn, b_trn)
      # compute test error
      sum((y[tst_idx] - yhat_tst)^2)
    })
    cv_errs_lams
  })
  
  #  weighted average according to group size (some groups may have 
  # +/- 1 member depending on whether K can't divide sizes evenly) of
  # cross validation error for each value of lambda
  cv.error <- apply(cv_errs_lams_folds, 1, 
                    function(cv_errs_folds) {
                      sum(cv_errs_folds * tabulate(folds))
                    })/n
  
  # extract the optimal value of our tuning parameter lambda
  # and (re)compute the corresponding coefficient estimates
  best.lam <- lam.vec[cv.error == min(cv.error)]
  best.fit <- ridge_coef(X, y, best.lam)
  b1 <- best.fit$b1
  b <- best.fit$b
  list(b1 = b1, b = b, best.lam = best.lam, cv.error = cv.error)
}

ridge_cv_cpp <- function(X, y, lam.vec, K) {
  # perform K-fold cross-validation on the ridge regression 
  # estimation problem over tuning parameters given in lam.vec
  n <- nrow(X); p <- ncol(X); L <- length(lam.vec)
  
  # groups to cross-validate over
  folds <- cut(1:n, breaks = K, labels = F)
  # get indices of training subset
  train_idxs <- lapply(1:K, function(i) !(folds %in% i))
  
  # preallocate empty data structure to store our CV errors 
  # for each lambda & fold
  cv_errs_lams_folds <- matrix(0, nrow = L, ncol = K)
  cv_errs_lams_folds <- sapply(train_idxs, function(trn_idx) {
    tst_idx <- !trn_idx # find test subset indices
    
    # subset data and remove intercept column in the design matrix
    Xm1_trn <- X[trn_idx, -1]; y_trn <- y[trn_idx]; n_trn <- nrow(Xm1_trn) 
    xbar_trn <- colMeans(Xm1_trn); ytilde_trn <- y_trn - mean(y_trn) 
    
    # center each predictor according to its mean
    # Xtilde_trn <- (Xm1_trn - tcrossprod(rep(1, nrow(Xm1_trn)), xbar_trn)) 
    Xtilde_trn <- scale(Xm1_trn) * sqrt(n_trn/(n_trn - 1))
    
    # compute the SVD on the centered design matrix
    # Note: This is done before computing our coefficients for each
    # lambda since the SVD will be identical across the vector of
    # tuning parameters
    Xtilde_trn_svd <- svd_cpp(Xtilde_trn)
    U <- Xtilde_trn_svd$u; d <- Xtilde_trn_svd$d; V <- Xtilde_trn_svd$v
    d2 <- d^2
    
    # preallocate empty data structure to store CV errors 
    # for each value of lambda
    cv_errs_lams <- vector(mode = 'numeric', length = L)
    cv_errs_lams <- sapply(lam.vec, function(lam) {
      # compute the inverse (D^T D + lambda I_{p-1})^{-1} D^T
      Dstar <- diag(d/(d2 + lam))
      # compute ridge coefficients
      b_trn <- V %*% (Dstar %*% crossprod(U, ytilde_trn)) * 
        1/attr(Xtilde_trn, "scaled:scale") * sqrt(n_trn/(n_trn - 1))
      b1_trn <- mean(y_trn) - crossprod(xbar_trn, b_trn)
      
      # fit test data
      yhat_tst <- X[tst_idx,] %*% c(b1_trn, b_trn)
      # compute test error
      sum((y[tst_idx] - yhat_tst)^2)
    })
    cv_errs_lams
  })
  
  #  weighted average according to group size (some groups may have 
  # +/- 1 member depending on whether K can't divide sizes evenly) of
  # cross validation error for each value of lambda
  cv.error <- apply(cv_errs_lams_folds, 1, 
                    function(cv_errs_folds) {
                      sum(cv_errs_folds * tabulate(folds))
                    })/n
  
  # extract the optimal value of our tuning parameter lambda
  # and (re)compute the corresponding coefficient estimates
  best.lam <- lam.vec[cv.error == min(cv.error)]
  best.fit <- ridge_coef(X, y, best.lam)
  b1 <- best.fit$b1
  b <- best.fit$b
  list(b1 = b1, b = b, best.lam = best.lam, cv.error = cv.error)
}

#===== start tests =====#
set.seed(124)

# set parameters
n <- 250
p <- 20
sigma <- 2
lams <- seq(0, 1, length.out = 500); lams <- lams[lams != 0]
K <- 5

# generate data (polynomial curve)
x <- runif(n, -1, 1); X <- sapply(0:p, function(k) x^k)
beta <- sample((-1)^(1:(p + 1)) * (1:(p + 1)))
y <- X %*% beta + rnorm(n, 0, sigma)

# run cross-validated ridge regression
pt <- proc.time()
rcv <- ridge_cv(X = X, y = y, lam.vec = lams, K = K)
proc.time() - pt

# fit estimates
bhat <- c(rcv$b1, rcv$b)
yhat <- X %*% bhat

plot(y ~ x, cex = 0.75, pch = 21, bg = 'white', xlab = "X", ylab = "Y",
     sub = substitute(paste(lambda[min], " = ", lam_min), 
                      list(lam_min = round(rcv$best.lam, 4))))
lines(yhat[order(x)] ~ sort(x), lwd = 2)

plot(rcv$cv.error ~ lams, type = 'l', lwd = 2, main = "Test Error", 
     xlab = expression(lambda), ylab = "Average K-fold MSE")

mb <- microbenchmark(ridge_cv(X, y, c(1, 2), K = 5), ridge_cv_cpp(X, y, c(1, 2), K = 5), times = 1e2, unit = 'us')
boxplot(mb, outline = F)
mb








