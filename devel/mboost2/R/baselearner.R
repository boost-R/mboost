
library("splines")

get_vars <- function(data, index, classes = c("numeric", "factor")) {

    x <- data[, index, drop = FALSE]
    stopifnot(all(sapply(x, class) %in% classes))
    x
}

spDes <- function(x, nknots = 20, degree = 3) {

    offset <- diff(range(x, na.rm = TRUE)) / 10000
    diffs <- (diff(range(x, na.rm = TRUE)) + 2 * offset) / (nknots - 1)
    knots <- seq(from = min(x, na.rm = TRUE) - offset - degree * diffs, 
                 to = max(x, na.rm = TRUE) + offset + degree * diffs,
                 by = diffs)

    design <- function(x)
        splineDesign(knots, x, ord = degree + 1, outer.ok = TRUE) 

    design
}

npp <- function(data, index, setup = spDes, differences = 2, ...) {

    x <- get_vars(data, index, classes = "numeric")
    stopifnot(length(x) == 1)
    nam <- names(x)
    x <- x[[1]]

    xu <- unique(x)
    index <- match(x, xu)

    design <- setup(xu, ...)
    Xnp <- design(xu)
    ### FIXME: spaltenweise Mittelwert abziehen?
    D <- diff(diag(ncol(Xnp)), differences = differences)
    Xnp <- tcrossprod(Xnp, D) %*% solve(tcrossprod(D))
    attr(Xnp, "index") <- index

    Xp <- poly(xu, degree = differences - 1)
    Xp <- lapply(1:ncol(Xp), function(i) {
        ret <- Xp[,i , drop = FALSE]
        attr(ret, "index") <- index
        ret
    })
 
    Xone <- matrix(1, ncol = 1, nrow = 1)

    ret <- c(Xone, Xp, list(Xnp))
    ret[[1]] <- as.matrix(ret[[1]])
    attr(ret[[1]], "index") <- rep.int(1, length(x))
    names(ret) <- c("(Intercept)", 
        paste("poly(", nam, ", ", 1:(length(ret) - 2), ")", sep = ""), 
        paste("s(", nam, ")", sep = ""))
    return(ret)
}

prebase <- function(data, index, by = NULL, ...) {

     if (length(index) == 1) {
         baseX <- npp(data, index, ...)
         if (!is.null(by))
             baseX <- do_by(baseX, data, by)
         return(baseX)
     }

     x <- get_vars(data, index, classes = "numeric")

     dmat <- vector(mode = "list", length = length(x))
     for (i in 1:length(x)) 
         dmat[[i]] <- npp(data, index[i], ...)

     ind <- vector(mode = "list", length = length(x))
     for (i in 1:length(x))
         ind[[i]] <- 1: length(dmat[[i]])
     ind <- as.matrix(do.call("expand.grid", ind))[-1,, drop = FALSE]

     tmp <- vector(mode = "list", length = ncol(ind))
     baseX <- vector(mode = "list", length = nrow(ind))

     for (i in 1:nrow(ind)) {
         for (j in 1:ncol(ind)) {
             tmp[[j]] <- dmat[[j]][[ind[i,j]]]
             names(tmp)[j] <- names(dmat[[j]])[ind[i,j]]
         }
         baseX[[i]] <- Tensor(tmp)
         names(baseX)[i] <- paste(names(tmp), collapse = ":")
     }
     if (!is.null(by))
         baseX <- do_by(baseX, data, by)
     return(baseX)
}

do_by <- function(base, data, by) {

    x <- get_vars(data, by, classes = "numeric")
    stopifnot(length(x) == 1)
    nam <- names(x)
    x <- x[[1]]

    xu <- unique(x)
    index <- match(x, xu)

    attr(xu, "index") <- index

    ret <- lapply(base, tensor, y = xu)
    names(ret) <- paste(names(base), "*", nam, sep = "")
    ret
}

base2fit <- function(base, weights, df = 4) {

    ret <- lapply(base, function(Xu) {

        index <- attr(Xu, "index")
        if (is.null(index)) index <- 1:NROW(Xu)
        wu <- as.vector(tapply(weights, index, sum))

        Xw <- Xu * wu
        XtX <- crossprod(Xw, Xu)
        if (NCOL(Xu) > 1) {
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
    names(ret) <- names(base)
    ret
}

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

boost <- function(y, data, index, weights = NULL, mstop = 100, ...) {

    if (is.null(weights)) weights <- rep(1, nrow(data))

    pb <- list()
    if (is.list(index)) {
        for (i in 1:length(index))
            pb <- c(pb, prebase(data, index[[i]], ...))
    } else {
        pb <- prebase(data, index, ...)
    }

    bf <- base2fit(pb, weights)

    offset <- mean(y)
    f <- offset
    u <- y - f

    sel <- numeric(mstop)
    for (m in 1:mstop) {

       mse <- numeric(length(bf)) 
       for (i in 1:length(bf)) {
            beta <- bf[[i]](u)
            uhat <- pb[[i]] %*% beta
            mse[i] <- sum((uhat[attr(pb[[i]], "index")] - u)^2)
       }
       istar <- sel[m] <- which.min(mse)

       f <- f + 0.1 * (pb[[istar]] %*% bf[[istar]](u))[attr(pb[[istar]], "index")]

       u <- y - f
    }
    list(fhat = f, sel = names(pb)[sel])
}

data("bodyfat", package = "mboost")
y <- bodyfat$DEXfat
b <- boost(y, bodyfat, list(1, 3:4, 5, 6))

