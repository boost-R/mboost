### what happens to weights
### when calculating knots etc?
bbs <- function(x, z = NULL, df = 4, knots = 20, degree = 3, differences = 2,
                center = FALSE, xname = NULL, zname = NULL) {

    if (is.null(xname)) xname <- deparse(substitute(x))
    if (is.null(zname)) zname <- deparse(substitute(z))

    if (!is.numeric(z) && (is.factor(z) && length(unique(z)) != 2))
        stop(sQuote("z"), " must be binary or numeric")

    if(is.factor(z) && length(unique(z)) == 2)
        ## FIXME is there a more elegant way to produce a binary with 0/1?
        z <- as.numeric(z[, drop = TRUE]) - 1

    if (all(x %in% c(0, 1)))
        return(bols(x = x, z = z, xname = xname, zname = zname,
                    center = center || all(x == 1)))

    if (is.factor(x) || (df <= 2 && !center))
        return(bols(x = x, z = z, xname = xname, zname = zname))


    if (!differences %in% 1:3)
        stop(sQuote("differences"), " are not in 1:3")
    if ((!center) && (df < differences))
        stop(sQuote("df"), " is less than ", sQuote("differences"))
    if(center && (degree < (differences-1)))
        stop(sQuote("degree"), " is less than ", sQuote("differences"), "-1")
    if (length(unique(x)) < 6)
        stop(sQuote(xname), " has less than 6 unique values")

    if (length(unique(round(diff(knots), 10))) > 1)
            warning("non-equidistant ", sQuote("knots"),
                    " might be inappropriate")

    X <- matrix(x, ncol = 1)

    if (is.null(z)) {
        ux <- sort(unique(x), na.last = TRUE)
        index <- match(x, ux)
        x <- ux
    }

    cc <- complete_cases(x = x, z = z)
    if (any(!cc)) {
        if (is.null(index)) index <- 1:length(x)
        x <- x[cc]
        if (!is.null(z)) z <- z[cc]
        index[index %in% which(!cc)] <- NA
        
    }

    dpp <- function(weights) {

        oweights <- weights
        if (!is.null(index)) weights <- as.vector(tapply(weights, index, sum))

        ### knots may depend on weights
        boundary.knots <- range(x, na.rm = TRUE)
        bnw <- range(x[weights > 0], na.rm = TRUE)
        if (!isTRUE(all.equal(boundary.knots, bnw)))
            warning("knots (and therefore model) depend on observations with zero weight")

        if (length(knots) == 1) {
            knots <- seq(from = boundary.knots[1], to = boundary.knots[2], length = knots+2)
            knots <- knots[2:(length(knots) - 1)]
        }

        newX <- function(x, z = NULL) {
            X <- bs(x, knots = knots, degree = degree, intercept = TRUE,
                    Boundary.knots = boundary.knots)
            class(X) <- "matrix"
            X <- as(X, "CsparseMatrix")
            if (!is.null(z))
                X <- X * z
            if (center) {
                K <- diff(Diagonal(ncol(X)), differences = differences)
                X <- tcrossprod(X, K) %*% solve(tcrossprod(K))
            }
            return(X)
        }
        X <- newX(x, z)

        if (center) {
            K <- Diagonal(ncol(X))
        } else {
            K <- diff(Diagonal(ncol(X)), differences = differences)
            K <- crossprod(K, K)
        }

        lambda <- df2lambda(X, df = df, dmat = K, weights = weights)

        XtX <- crossprod(X * weights, X) + lambda * K
        if (isSymmetric(XtX)) XtX <- forceSymmetric(XtX)
	
        modelmatrixfun <- function(newdata = NULL) {
            if (is.null(newdata)) {
                if (!is.null(index)) return(X[index,])
                return(X)
            }
            return(newX(x = newdata[[xname]], z = newdata[[zname]]))
        }

        fitfun <- function(y) {

            if (!is.null(index)) y <- as.vector(tapply(oweights * y, index, sum))
            coef <- solve(XtX, crossprod(X, y))

            predictfun <- function(newdata = NULL) {
                XX <- modelmatrixfun(newdata = newdata)
                return(as.vector(XX %*% coef))
            }
            ret <- list(model = coef, predict = predictfun, 
                        fitted = function() {
                            if (!is.null(index)) return(as.vector(X %*% coef)[index])
                            return(as.vector(X %*% coef))
                        })
            class(ret) <- c("basefit", "baselm")
            ret
        }
        ret <- list(fit = fitfun, modelmatrix = modelmatrixfun, 
                    hatmatrix = function() {
            if (!is.null(index)) {
                X <- as.matrix(X)
                return(as.matrix(tcrossprod(X[index,] %*% solve(XtX), X[index,] * oweights)))
            }
            return(as.matrix(tcrossprod(X %*% solve(XtX), X * oweights)))
            })
        class(ret) <- "basisdpp"
        ret
    }
    attr(X, "dpp") <- dpp
    return(X)
}

bspatial <- function(x, y, z = NULL, df = 5, xknots = 20, yknots = 20,
                     degree = 3, differences = 2, center = FALSE, xname = NULL,
                     yname = NULL, zname = NULL) {

    if (!is.numeric(x) || !is.numeric(y))
        stop(sQuote("x"), " and ", sQuote("y"), " must be numeric")
    if (!is.numeric(z) && (is.factor(z) && length(unique(z)) != 2))
        stop(sQuote("z"), " must be binary or numeric")

    if(is.factor(z) && length(unique(z)) == 2)
        ## FIXME is there a more elegant way to produce a binary with 0/1?
        z <- as.numeric(z[, drop=T]) - 1

    if (is.null(xname)) xname = deparse(substitute(x))
    if (is.null(yname)) yname = deparse(substitute(y))
    if (is.null(zname)) zname = deparse(substitute(z))

    cc <- complete_cases(x = x, y = y, z = z)

#    if (df <= 2) stop(sQuote("df"), " must be greater two")
    if (!differences %in% 1:3)
        stop(sQuote("differences"), " are not in 1:3")
    if ((!center) && (df < differences^2))
        stop(sQuote("df"), " is less than ", sQuote("differences^2"))
    if (length(unique(x)) < 6)
        stop(sQuote(xname), " has less than 6 unique values")
    if (length(unique(y)) < 6)
        stop(sQuote(yname), " has less than 6 unique values")

    if (length(unique(diff(xknots))) > 1)
            warning("non-equidistant ", sQuote("xknots"),
                    " might be inappropriate")
    if (length(unique(diff(yknots))) > 1)
            warning("non-equidistant ", sQuote("yknots"),
                    " might be inappropriate")

    if (length(xknots) == 1) {
        xknots <- seq(from = min(x, na.rm = TRUE),
                     to = max(x, na.rm = TRUE), length = xknots + 2)
        xknots <- xknots[2:(length(xknots) - 1)]
    }
    if (length(yknots) == 1) {
        yknots <- seq(from = min(y, na.rm = TRUE),
                     to = max(y, na.rm = TRUE), length = yknots + 2)
        yknots <- yknots[2:(length(yknots) - 1)]
    }

    X <- matrix(x, ncol = 1)

    dpp <- function(weights) {

        if (any(!cc)) weights <- weights[cc]

        ### knots may depend on weights
        xboundary.knots <- range(x[cc], na.rm = TRUE)
        xbnw <- range(x[cc][weights > 0], na.rm = TRUE)
        if (!isTRUE(all.equal(xboundary.knots, xbnw)))
            warning("knots (and therefore model) depend on observations with zero weight")

        if (length(xknots) == 1) {
            xknots <- seq(from = xboundary.knots[1], to = xboundary.knots[2], length = xknots+2)
            xknots <- xknots[2:(length(xknots) - 1)]
        }

        yboundary.knots <- range(y[cc], na.rm = TRUE)
        ybnw <- range(y[cc][weights > 0], na.rm = TRUE)
        if (!isTRUE(all.equal(yboundary.knots, ybnw)))
            warning("knots (and therefore model) depend on observations with zero weight")

        if (length(yknots) == 1) {
            yknots <- seq(from = yboundary.knots[1], to = yboundary.knots[2], length = yknots+2)
            yknots <- yknots[2:(length(yknots) - 1)]
        }

        newX <- function(x, y, z = NULL, weights = NULL, na.rm = TRUE) {

          if (na.rm) {
                x <- x[cc]
                y <- y[cc]
                if (!is.null(z))
                    z <- z[cc]
                if (!is.null(weights))
                    weights <- weights[cc]
            }

            Xx <- bs(x, knots = xknots, degree = degree, intercept = TRUE,
                     Boundary.knots = xboundary.knots)
            class(Xx) <- "matrix"
            Xx <- Matrix(Xx)
            Xy <- bs(y, knots = yknots, degree = degree, intercept = TRUE,
                     Boundary.knots = yboundary.knots)
            class(Xy) <- "matrix"
            Xy <- Matrix(Xy)

            X <- kronecker(Xx, Matrix(1, nc = ncol(Xy))) * kronecker(Matrix(1, nc = ncol(Xx)), Xy)
            if (!is.null(z))
                X <- X * z
            return(X)
        }
        X <- newX(x, y, z)
        Xna <- X
        if (any(!cc))
            Xna <- newX(x, y, z, weights = weights, na.rm = FALSE)

        xd <- length(xknots) + degree + 1
        yd <- length(yknots) + degree + 1

        Kx <- diff(Diagonal(xd), differences = differences)
        Kx <- crossprod(Kx, Kx)
        Ky <- diff(Diagonal(yd), differences = differences)
        Ky <- crossprod(Ky, Ky)
        K <- kronecker(Kx, Diagonal(yd)) + kronecker(Diagonal(xd), Ky)

        L <- 0
        if(center) {
            L <- eigen(K, symmetric=TRUE, EISPACK=TRUE)
            L$vectors <- L$vectors[,1:(ncol(X)-differences^2)]
            L$values <- sqrt(L$values[1:(ncol(X)-differences^2)])
            L <- L$vectors%*%Diagonal(1/L$values)
            X <- X%*%L
            K <- Diagonal(ncol(X))
        }

        lambda <- df2lambda(X, df = df, dmat = K, weights = weights)
        Xw <- X * weights
        XtX <- crossprod(Xw, X) + lambda * K
        if (isSymmetric(XtX)) XtX <- forceSymmetric(XtX)
        #### if (any(!cc)) rm(X)

            modelmatrixfun <- function(newdata = NULL) {
                if (is.null(newdata)) return(Xna)
                nX <- newX(x = newdata[[xname]], y = newdata[[yname]],
                           z = newdata[[zname]], na.rm = FALSE)
                if(center) {
                    nX <- nX%*%L
                }
            }


        fitfun <- function(y) {
            coef <- solve(XtX, crossprod(Xw, y))

            predictfun <- function(newdata = NULL) {
                XX <- modelmatrixfun(newdata = newdata)
                as.vector(XX %*% coef)
            }

            ret <- list(model = coef, predict = predictfun,
                        fitted = function() as.vector(Xna %*% coef))
            class(ret) <- c("basefit", "baselm")
            ret
        }
        ret <- list(fit = fitfun, modelmatrix = modelmatrixfun, 
                    hatmatrix = function() 
                    as.matrix(tcrossprod(X %*% solve(XtX), Xw)))
        class(ret) <- "basisdpp"
        ret
    }
    attr(X, "dpp") <- dpp
    return(X)
}

model.matrix.basisdpp <- function(object, newdata = NULL)
    object$modelmatrixfun(newdata)
