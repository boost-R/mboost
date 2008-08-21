
library("splines")

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
        if (degree == 0) return(matrix(1, nrow = nrow(data)))
        x <- get_var(data, varid)
        ret <- predict(X, newdata = x)[, degree, drop = FALSE]
        attr(ret, "index") <- attr(x, "index")
        class(ret) <- "polyDes"
        ret
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

tm <- Terms(iris, 1:2)
mm <- model.matrix(tm, iris)
w <- rep(1, nrow(iris))
fit <- mm2fit(mm, w, df = 3)
