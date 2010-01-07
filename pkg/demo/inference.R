
library("mboost")

data("bodyfat", package = "mboost")

object <- glmboost(DEXfat ~ ., data = bodyfat, center = TRUE)
a <- mboost:::basesel(object, folds = cv(model.weights(object), type = "subsampling", B = 100))

object <- gamboost(DEXfat ~ ., data = bodyfat)
a <- mboost:::basesel(object)
matplot(t(a), type = "l")

apply(a, 1, max)

nd <- data.frame(hipcirc = 90:120)
b <- mboost:::fitsel(object, newdata = nd, which = "hipcirc")

apply(b, 1, max)

n <- 500
x <- seq(from = 0, to = 6 * pi, length = n)
y <- 2 * sin(x) + rnorm(n)
object <- gamboost(y ~ x)
d <- mboost:::fitsel(object, which = 1)
a <- apply(d, 1, max)
plot(x, fitted(object), col = (a > 0.99) + 1)
lines(x, 2 * sin(x))
abline(h = 0)
