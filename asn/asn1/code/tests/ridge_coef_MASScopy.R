library(Rcpp)
sourceCpp("code/ridge_coef.cpp")
ridge_coef <- function(X, y, lam) {
  # Commented-out scaling parameters represent the transformations
  # used to make the output identical to that of 
  # coef(MASS::lm.ridge(y ~ X[,-1], lambda = lam))
  Xm1 <- X[,-1]; ytilde <- y - mean(y)
  # center each predictor according to its mean
  Xtilde <- scale(Xm1) * sqrt(n/(n - 1))

  # compute the SVD on the centered design matrix
  Xtilde_svd <- svd(Xtilde)
  U <- Xtilde_svd$u; d <- Xtilde_svd$d; V <- Xtilde_svd$v
  
  # compute the inverse (D^T D + lambda I_{p-1})^{-1} D^T
  Dstar <- diag(d/(d^2 + lam))
  
  # compute ridge coefficients
  b <- V %*% (Dstar %*% crossprod(U, ytilde)) * 1/attr(Xtilde, "scaled:scale") * sqrt(n/(n - 1))
  b0 <- mean(y) - crossprod(attr(Xtilde, "scaled:center"), b) 
  c(b0 = b0, b = b)
}

ridge_cv <- function(X, y, lam_vec, K) {
  # perform K-fold cross-validation on the ridge regression 
  # estimation problem over tuning parameters given in lam.vec
  n <- nrow(X); p <- ncol(X); L <- length(lam_vec)
  
  # groups to cross-validate over
  folds <- cut(1:n, breaks = K, labels = F)
  # get indices of training subset
  train_idxs <- lapply(1:K, function(i) !(folds %in% i))
  
  # data structure to store our CV errors 
  cv_errs_lams_folds <- lapply(train_idxs, function(trn_idx) {
    tst_idx <- !trn_idx # find test subset indices
    
    # subset data and remove intercept column in the design matrix
    Xm1_trn <- X[trn_idx, -1]; y_trn <- y[trn_idx]; ytilde_trn <- y_trn - mean(y_trn) 
    n_trn <- nrow(Xm1_trn) 
    
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
    cv_errs_lams <- sapply(lam_vec, function(lam) {
      # compute the inverse (D^T D + lambda I_{p-1})^{-1} D^T
      Dstar <- diag(d/(d2 + lam))
      
      # compute ridge coefficients
      b_trn <- V %*% (Dstar %*% crossprod(U, ytilde_trn)) * 
        1/attr(Xtilde_trn, "scaled:scale") * sqrt(n_trn/(n_trn - 1))
      b0_trn <- mean(y_trn) - crossprod(attr(Xtilde_trn, "scaled:center"), b_trn)
      bhat_trn <- c(b0_trn, b_trn)
      
      # fit test data
      yhat_trn <- X[trn_idx,] %*% bhat_trn
      yhat_tst <- X[tst_idx,] %*% bhat_trn
      
      # compute test error
      c(trn_err = sum((y[trn_idx] - yhat_trn)^2)/sum(trn_idx),
        tst_err = sum((y[tst_idx] - yhat_tst)^2)/sum(tst_idx))
    })
    cv_errs_lams
  })
  cv_errs_lams_folds <- simplify2array(cv_errs_lams_folds)
  
  #  weighted average according to group size (some groups may have 
  # +/- 1 member depending on whether K can't divide sizes evenly) of
  # cross validation error for each value of lambda
  cv_errs_trn <- apply(cv_errs_lams_folds["trn_err",,], 1, 
                       function(cv_errs_folds) sum(cv_errs_folds * tabulate(folds)))/n
  cv_errs_tst <- apply(cv_errs_lams_folds["tst_err",,], 1, 
                       function(cv_errs_folds) sum(cv_errs_folds * tabulate(folds)))/n
  
  # extract the optimal value of our tuning parameter lambda
  # and (re)compute the corresponding coefficient estimates
  best_test_err <- min(cv_errs_tst)
  best_lam <- lam_vec[cv_errs_tst == best_test_err]
  best_b <- ridge_coef(X, y, best_lam)
  
  # overall loss, not just train/testing subset
  yhat <- X %*% best_b
  loss <- (yhat - y)^2 
  loss_avg <- mean(loss)
  loss_se <- sd(loss)/length(loss)
  
  list(best_b = best_b, best_lam = best_lam, cv_errs_trn = cv_errs_trn, 
       cv_errs_tst = cv_errs_tst, best_test_err = best_test_err, loss_avg = loss_avg, loss_se = loss_se)
}

#===== start tests =====#
set.seed(124)

# set parameters
n <- 250
p <- 15
sigma <- 1/2
lams <- seq(0, 2, length.out = 250); lams <- lams[lams != 0]
K <- 5; Ks <- seq(2, 20, 2)

# generate data (polynomial curve)
x <- runif(n, -1, 1); X <- sapply(0:p, function(k) x^k)
beta <- sample((-1)^(1:(p + 1)) * (1:(p + 1)))
y <- X %*% beta + rnorm(n, 0, sigma)

# run cross-validated ridge regression
pt <- proc.time()
rcv <- ridge_cv(X = X, y = y, lam_vec = lams, K = K)
proc.time() - pt

# fit estimates
yhat <- X %*% rcv$best_b

plot(y ~ x, cex = 0.75, pch = 21, bg = 'white', xlab = "X", ylab = "Y", 
     sub = substitute(paste(lambda[min], " = ", lam_min), 
                      list(lam_min = round(rcv$best_lam, 4))))
lines(yhat[order(x)] ~ sort(x), lwd = 2)

plot(rcv$cv_errs_trn ~ lams, type = 'l', lwd = 2, lty = 'dotdash', 
     ylim = range(c(rcv$cv_errs_trn, rcv$cv_errs_tst)),
     main = "CV Errors", xlab = expression(lambda), ylab = "Average K-fold Error")
lines(rcv$cv_errs_tst ~ lams, lwd = 2)


pt <- proc.time()
rcv_ks <- lapply(Ks, function(k) ridge_cv(X = X, y = y, lam_vec = lams, K = k))
proc.time() - pt

rcv_ks_best <- sapply(rcv_ks, function(r) 
  c(best_lam = r$best_lam, loss_avg = r$loss_avg, loss_se = r$loss_se))
plot(rcv_ks_best["best_lam",] ~ Ks, type = 'o', pch = 21, bg = 'white',
     xlab = "K", ylab = expression(lambda[min]))
plot(rcv_ks_best["loss_avg",] ~ Ks, type = 'o', pch = 21, bg = 'white',
     xlab = "K", ylab = "Mean K-fold Loss")
plot(rcv_ks_best["loss_se",] ~ Ks, type = 'o', pch = 21, bg = 'white',
     xlab = "K", ylab = "Std Err. K-fold Loss")





