
library("splines")

get_vars <- function(data, index, classes = c("numeric", "factor")) {

    x <- data[, index, drop = FALSE]
    stopifnot(all(sapply(x, class) %in% classes))
    x
}

bl_pLS <- function(data, index, weights = rep.int(1, nrow(data)), 
                   nknots = 20, degree = 3, df = 4, differences = 2, ...) {

    x <- get_vars(data, index, classes = "numeric")
    stopifnot(length(x) == 1)
    x <- x[[1]]

    offset <- diff(range(x, na.rm = TRUE)) / 10000
    diffs <- (diff(range(x, na.rm = TRUE)) + 2 * offset) / (nknots - 1)
    knots <- seq(from = min(x, na.rm = TRUE) - offset - degree * diffs, 
                 to = max(x, na.rm = TRUE) + offset + degree * diffs,
                 by = diffs)

    xu <- unique(x)
    inx <- match(x, xu)
    wu <- as.vector(tapply(weights, inx, sum))

    design <- function(x)
        splineDesign(knots, x, ord = degree + 1, outer.ok = TRUE) 

    Xu <- design(xu)
    K <- diff(diag(ncol(Xu)), differences = differences)
    K <- crossprod(K, K)
    K <- K * mboost:::df2lambda(Xu, df = df, dmat = K, weights = wu)

    Xw <- Xu * wu
    XtX <- crossprod(Xw, Xu)
    Xsolve <- solve(XtX + K, t(Xu))
    
    fit <- function(y) {
        yu <- tapply(y * weights, inx, sum)
        as.vector(Xsolve %*% yu)
    }
    list(fit = fit, design = design)
}

x <- round(runif(10000, max = pi), 1)
y <- sin(x) + rnorm(length(x), sd = 0.1)
w <- rpois(length(x), lambda = 1)

tmp <- bl_pLS(data.frame(x), 1, nknots = 5, weights = w)
lm.wfit(x = tmp$design(x), y = y, w = w)$coef

tmp$fit(y)

