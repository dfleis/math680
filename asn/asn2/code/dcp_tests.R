
n <- 200
x <- seq(-4, 4, length.out = n)
y <- seq(-4, 4, length.out = n)

M1 <- outer(x, y, function(x, y) (abs(x) + abs(y))^3 <= 5 * x + 7)
M3 <- outer(x, y, function(x, y) sqrt(x^2 + 4) + 2 * y <= -5 * x)
M7 <- outer(x, y, function(x, y) x * y >= 1)
#M8 <- outer(x, y, function(x, y) (y * log(y/(2 * x)) <= y + x - 30) & (x > 0) & (y > 0))
image.plot(M7)
sum(M7)












