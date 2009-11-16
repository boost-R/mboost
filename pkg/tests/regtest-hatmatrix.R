
library("mboost")

set.seed(290875)

df <- function(n = 100) {
    x <- matrix(runif(n * 10), ncol = 10)
    y <- x %*% c(1:3, rep(0, 7)) + rnorm(n)
    data.frame(y = y, x)
}

mydf <- df()

w <- c(rep(0, 50), rep(1, 50))
mod1 <- glmboost(y ~ ., data = mydf, weights = w)
cf1 <- coef(mod1)

### hat matrix: fast for linear models
H <- attr(hatvalues(mod1), "hatmatrix")
stopifnot(max(abs(H %*% (mydf$y - weighted.mean(mydf$y, w)) - fitted(mod1) + weighted.mean(mydf$y, w))) <
          sqrt(.Machine$double.eps))

mb <- function(x) bols(x, intercept = FALSE)
mydf$int <- rep(1, nrow(mydf))
## used to be: mod2 <- gamboost(y ~ ., data = mydf, weights = w, baselearner = mb) ## but doesn't work as .mboost functions deprecated.
mod2 <- gamboost(y ~ ., data = mydf, weights = w, baselearner = mb)
cf2 <- coef(mod2)
stopifnot(all(round(unlist(cf2), 5) %in% round(cf1, 5)))

### hat matrix: less faster but in C
H <- attr(hatvalues(mod2), "hatmatrix")
stopifnot(max(abs(H %*% (mydf$y - weighted.mean(mydf$y, w)) - fitted(mod2) + weighted.mean(mydf$y, w))) <
          sqrt(.Machine$double.eps))

### hat matrix: directly in R
mod2$family <- Laplace()
H <- attr(hatvalues(mod2), "hatmatrix")
stopifnot(max(abs(H %*% (mydf$y - weighted.mean(mydf$y, w)) - fitted(mod2) + weighted.mean(mydf$y, w))) <
          sqrt(.Machine$double.eps))
