
library("splines")

tensor <- function(x, y) {

    xi <- attr(x, "index")
    yi <- attr(y, "index")
    if (is.null(xi)) xi <- 1:NROW(x)
    if (is.null(yi)) yi <- 1:NROW(y)
    stopifnot(length(xi) == length(yi))

    xiyi <- paste(xi, yi, sep = "_")
    d <- !duplicated(xiyi)
    zi <- match(xiyi, xiyi[d])

    if (!is.matrix(x)) x <- matrix(x, ncol = 1)
    if (!is.matrix(y)) y <- matrix(y, ncol = 1)

    X <- x[xi[d],,drop = FALSE] 
    Y <- y[yi[d],,drop = FALSE]
    ret <- NULL
    for (j in 1:ncol(y))
        ret <- cbind(ret, X * Y[, j])
    attr(ret, "index") <- zi   
    return(ret)
}

Tensor <- function(...) {

    X <- list(...)[[1]]  
    if (length(X) == 1) 
        return(X[[1]])
    ret <- X[[1]]
    for (i in 2:length(X))
        ret <- tensor(ret, X[[i]])
    ret
}

get_var <- function(data, varid, classes = c("numeric", "factor")) {

    stopifnot(length(varid) == 1 && varid %in% 1:ncol(data))
    x <- data[, varid, drop = TRUE]
    stopifnot(class(x) %in% classes)
    xu <- unique(x)
    index <- match(x, xu)     
    attr(xu, "index") <- index
    xu
}

spDes <- function(data, varid, nknots = 20, degree = 3, differences = 2) {

    x <- get_var(data, varid)
    offset <- diff(range(x, na.rm = TRUE)) / 10000
    diffs <- (diff(range(x, na.rm = TRUE)) + 2 * offset) / (nknots - 1)
    knots <- seq(from = min(x, na.rm = TRUE) - offset - degree * diffs, 
                 to = max(x, na.rm = TRUE) + offset + degree * diffs,
                 by = diffs)

    function(data) {
        x <- get_var(data, varid)
        Xnp <- splineDesign(knots, x, ord = degree + 1, outer.ok = TRUE) 
        D <- diff(diag(ncol(Xnp)), differences = differences)
        ret <- tcrossprod(Xnp, D) %*% solve(tcrossprod(D))
        attr(ret, "index") <- attr(x, "index")
        class(ret) <- "spDes"
        ret
    }
}

polyDes <- function(data, varid, degree = 1) {

    x <- get_var(data, varid)
    if (degree > 0)
        X <- poly(x, degree = degree)
    function(data) {
        if (degree == 0) {
            ret <- matrix(1, nrow = 1, ncol = 1)
            attr(ret, "index") <- rep(1, nrow(data))
            return(ret)
        }
        x <- get_var(data, varid)
        ret <- predict(X, newdata = x)[, degree, drop = FALSE]
        attr(ret, "index") <- attr(x, "index")
        class(ret) <- "polyDes"
        ret
    }
}

factorDes <- function(data, varid) {

    x <- get_var(data, varid, class = "factor")
    function(data) {
        ret <- model.matrix(~ x - 1)
        attr(ret, "index") <- attr(x, "index")
        class(ret) <- "factorDes"
        return(ret)
    }
}

PSpline <- function(data, varid, differences = 2, ...) {

    stopifnot(length(varid) == 1)
    ret <- vector(mode = "list", length = differences + 1)
    foo <- function(d) {
            degree <- d - 1
            function(data) polyDes(data, varid, degree = degree)
    }
    for (d in 1:differences)
        ret[[d]] <- foo(d)
    ret[[differences + 1]] <- function(data) spDes(data, varid, 
       differences = differences, ...)
    ret
}
class(PSpline) <- "design_generator"

Factor <- function(data, varid)
    return(function(data) factorDes(data, varid))
class(Factor) <- "design_generator"

### koennen wir auch baselearner haben, die aus
### interactions und haupteffekten bestehen (so dass immer
### die Haupteffekte mit ausgewaehlt werden) -> GT?
Terms <- function(data, varids, FUN = PSpline,  ...) {

    stopifnot(inherits(FUN, "design_generator"))
    base <- lapply(varids, FUN, data = data, ...)
    if (length(varids) == 1) {
        ret <- lapply(base[[1]], function(f) f(data))
        class(ret) <- "Terms"
        return(ret)
    }
    
    ind <- vector(mode = "list", length = length(base))
    for (i in 1:length(base))
        ind[[i]] <- 1:length(base[[i]])
    ind <- as.matrix(do.call("expand.grid", ind))

    baseX <- vector(mode = "list", length = nrow(ind))
    for (i in 1:nrow(ind)) {
    
        foo <- function(i) {
            myind <- ind[i,]
            function(data) {
                tmp <- vector(mode = "list", length = ncol(ind))
                for (j in 1:ncol(ind))
                    tmp[[j]] <- base[[j]][[myind[j]]](data)(data)
                ret <- Tensor(tmp)
                class(ret) <- "polyDes"
                if (any(sapply(tmp, class) == "spDes"))
                    class(ret) <- "spDes"
                ret
            }
        }
        baseX[[i]] <- foo(i)
    }
    class(baseX) <- "Terms"
    baseX
}

model.matrix.Terms <- function(object, data, ...)
    lapply(object, function(f) f(data))

mm2fit <- function(mm, weights, df = 4) {

    ret <- lapply(mm, function(Xu) {

        index <- attr(Xu, "index")
        if (is.null(index)) index <- 1:NROW(Xu)
        wu <- as.vector(tapply(weights, index, sum))

        Xw <- Xu * wu
        XtX <- crossprod(Xw, Xu)
        if (inherits(Xu, "spDes")) {
            lambda <- mboost:::df2lambda(Xu, df = df, dmat = diag(ncol(Xu)),
                          weights = wu)
            Xsolve <- solve(XtX + lambda * diag(ncol(Xu)), t(Xu))
        } else {
            Xsolve <- solve(XtX, t(Xu))
        }

        fit <- function(y) {
            yu <- tapply(y * weights, index, sum)
            as.vector(Xsolve %*% yu)
        }
        return(fit)
    })
    names(ret) <- names(mm)
    ret
}

mult <- function(x, y) {

    xi <- attr(x, "index")
    if (is.null(xi)) xi <- 1:NROW(x)
    
    (x %*% y)[xi]
}

boost <- function(y, mm, weights = NULL, mstop = 100, ...) {

    if (is.null(weights)) weights <- rep(1, length(y))
    bf <- mm2fit(mm, weights, ...)

    offset <- mean(y)
    f <- offset
    u <- y - f  

    sel <- numeric(mstop)
    for (m in 1:mstop) {

       mse <- numeric(length(bf)) 
       for (i in 1:length(bf)) {
            beta <- bf[[i]](u)
            uhat <- mult(mm[[i]], beta)
            mse[i] <- sum((uhat - u)^2)
       }
       istar <- sel[m] <- which.min(mse)

       f <- f + 0.1 * mult(mm[[istar]], bf[[istar]](u))

       u <- y - f
    }
    list(fhat = f, sel = names(mm)[sel])
}

data("bodyfat", package = "mboost")
y <- bodyfat$DEXfat
inputs <- which(colnames(bodyfat) != "DEXfat")
tm <- lapply(inputs, Terms, data = bodyfat)
names(tm) <- colnames(bodyfat)[inputs]
tm <- unlist(tm)
class(tm) <- "Terms"
b <- boost(y, model.matrix(tm, bodyfat), df = 4)
