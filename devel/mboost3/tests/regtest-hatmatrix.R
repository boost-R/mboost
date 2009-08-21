
library("mboost3")

set.seed(290875)

df <- function(n = 100) {
    x <- matrix(runif(n * 10), ncol = 10)
    y <- x %*% c(1:3, rep(0, 7)) + rnorm(n)
    data.frame(y = y, x)
}

mydf <- df()

w <- c(rep(0, 50), rep(1, 50))
mod <- glmboost(y ~ ., data = mydf, weights = w)

### hat matrix: fast for linear models
H <- attr(hatvalues(mod), "hatmatrix")
stopifnot(max(abs(H %*% (mydf$y - weighted.mean(mydf$y, w)) - fitted(mod) + weighted.mean(mydf$y, w))) < 
          sqrt(.Machine$double.eps))

### hat matrix: less faster but in C
H <- attr(mboost3:::hatvalues.mboost(mod), "hatmatrix")
stopifnot(max(abs(H %*% (mydf$y - weighted.mean(mydf$y, w)) - fitted(mod) + weighted.mean(mydf$y, w))) < 
          sqrt(.Machine$double.eps))

### hat matrix: directly in R
mod$family <- Laplace()
H <- attr(mboost3:::hatvalues.mboost(mod), "hatmatrix")
stopifnot(max(abs(H %*% (mydf$y - weighted.mean(mydf$y, w)) - fitted(mod) + weighted.mean(mydf$y, w))) < 
          sqrt(.Machine$double.eps))
