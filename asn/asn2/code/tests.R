a <- c(0, 0.25, 0.5, 0.75)
f <- function(x, y, a) {
  abs(x * y) + a * (x^2 + y^2)
}

x <- y <- seq(-1, 1, length.out = 50)
Za <- lapply(a, function(ai) outer(x, y, FUN = function(x1, x2) f(x1, x2, ai)))

i <- 3

persp(x, y, Za[[i]], zlab = expression(f(x, y)), theta = 15, phi = 30, 
      main = paste0("f(x, y; a), a = ", a[i]))











n <- 100
x <- runif(n, 0, 10)
xnorm <- x/sum(x)
xbar <- mean(x)

-sum(xnorm * log(xnorm))
-1/(n * xbar) * sum(x * log(x)) + log(n * xbar)


















