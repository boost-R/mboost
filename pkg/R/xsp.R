
hatMatTH <- function(x, w = NULL, df = 4) {
    n <- NROW(x)
    indx <- diag(n)
    x <- signif(x, 10)
    apply(indx, 2, function(y)
        predict(smoothbase(x = x, ux = unique(sort(x)), y = y, w = w, df = df),
                x = x)$y)
}

complete_cases <- function(x, y = NULL, z = NULL) {

    tmp <- list(x = x, y = y, z = z)
    tmp <- tmp[!sapply(tmp, is.null)]
    rowSums(sapply(tmp, is.na)) == 0
}

predict.baselist <- function(object, ...) {

    pr <- predict(object[[1]], ...)
    if (length(object) == 1) return(pr)
    for (i in 2:length(object)) {
        if (!any(is.na(pr)))
            break
        tmp <- predict(object[[i]], ...)
        pr[is.na(pr)] <- tmp[is.na(pr)]
    }
    return(pr)
}

fitted.baselist <- function(object) {

    pr <- fitted(object[[1]])
    if (length(object) == 1) return(pr)
    for (i in 2:length(object)) {
        if (!any(is.na(pr)))
            break
        tmp <- fitted(object[[i]])
        pr[is.na(pr)] <- tmp[is.na(pr)]
    }
    return(pr)
}

### what happens to weights
### when calculating knots etc?
bbs <- function(x, z = NULL, df = 4, knots = 20, degree = 3, differences = 2,
                center = FALSE, xname = NULL, zname = NULL) {

    cc <- complete_cases(x = x, z = z)

    if (is.null(xname)) xname <- deparse(substitute(x))
    if (is.null(zname)) zname <- deparse(substitute(z))

    if (all(x %in% c(0, 1)))
        return(bols(x = x, z = z, xname = xname, zname = zname, 
                    center = center || all(x == 1)))

    if (is.factor(x) || (df <= 2 && !center))
        return(bols(x = x, z = z, xname = xname, zname = zname))

    if (!is.numeric(z) && (is.factor(z) && length(unique(z)) != 2))
        stop(sQuote("z"), " must be binary or numeric")

    if(is.factor(z) && length(unique(z)) == 2)
        ## FIXME is there a more elegant way to produce a binary with 0/1?
        z <- as.numeric(z[, drop = TRUE]) - 1

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


    dpp <- function(weights) {

        if (any(!cc)) weights <- weights[cc]

        ### knots may depend on weights
        boundary.knots <- range(x[cc], na.rm = TRUE)
        bnw <- range(x[cc][weights > 0], na.rm = TRUE)
        if (!isTRUE(all.equal(boundary.knots, bnw)))
            warning("knots depend on weights")

        if (length(knots) == 1) {
            knots <- seq(from = bn[1], to = bn[2], length = knots+2)
            knots <- knots[2:(length(knots) - 1)]
        }

        newX <- function(x, z = NULL, weights = NULL, na.rm = TRUE) {
            if (na.rm) {
                x <- x[cc]
                if (!is.null(z))
                    z <- z[cc]
                if (!is.null(weights))
                    weights <- weights[cc]
            }
            X <- bs(x, knots = knots, degree = degree, intercept = TRUE,
                    Boundary.knots = boundary.knots)
            if (!is.null(z))
                X <- X * z
            if (center) {
                K <- diff(diag(ncol(X)), differences = differences)
                X <- tcrossprod(X, K) %*% solve(tcrossprod(K))
            }
            return(X)
        }
        X <- newX(x, z, weights = weights)
        Xna <- X
        if (any(!cc))
            Xna <- newX(x, z, weights = weights, na.rm = FALSE)

        if (center) {
            K <- diag(ncol(X))
        } else {
            K <- diff(diag(ncol(X)), differences = differences)
            K <- crossprod(K, K)
        }

        lambda <- df2lambda(X, df = df, dmat = K, weights = weights)

        Xw <- X * weights
        XtX <- crossprod(Xw, X)
        Xsolve <- tcrossprod(solve(XtX + lambda * K), Xw)

        fitfun <- function(y) {

            if (any(!cc)) y <- y[cc]
            coef <- Xsolve %*% y

            predictfun <- function(newdata = NULL) {
                if (is.null(newdata)) return(Xna %*% coef)
                nX <- newX(x = newdata[[xname]], z = newdata[[zname]], na.rm = FALSE)
                nX %*% coef
            }
            ret <- list(model = coef, predict = predictfun, fitted = function() Xna %*% coef)
            class(ret) <- c("basefit", "baselm")
            ret
        }
        ret <- list(fit = fitfun, hatmatrix = function() X %*% Xsolve)
        class(ret) <- "basisdpp"
        ret
    }
    attr(X, "dpp") <- dpp
    return(X)
}

predict.basefit <- function(object, newdata = NULL)
    object$predict(newdata)

fitted.basefit <- function(object)
    object$fitted()

coef.baselm <- function(object)
    object$model

coef.bssfit <- function(object)
    object$basemodel$coef

df2lambda <- function(X, df = 4, dmat = NULL, weights) {

#   if (df <= 2) stop(sQuote("df"), " must be greater than two")

    if (is.null(dmat)) {
        dmat <- diff(diag(ncol(X)), differences = 2)
        dmat <- crossprod(dmat, dmat)
    }

    # Cholesky decomposition

    A <- crossprod(X * weights, X) + dmat*10e-10
    Rm <- solve(chol(A))

    decomp <- svd(crossprod(Rm,dmat)%*%Rm)
    d <- decomp$d

    # df2lambda
    df2l <- function(lambda)
        (sum( 1/(1+lambda*d) ) - df)^2

    lower.l <- 0
    upper.l <- 5000
    lambda <- upper.l

    while (lambda >= upper.l - 200 ) {
        upper.l <- upper.l * 1.5

        tl <- try(lambda <- optimize(df2l, interval=c(lower.l,upper.l))$minimum,
        silent=T)
        if (class(tl)=="try-error") stop("problem of
        converting df into lambda cannot be solved - please increase value of
        df")
        lower.l <- upper.l-200
        if (lower.l > 1e+06){
            lambda <- 1e+06
            warning("lambda needs to be larger than 1e+06 for given value of df,
            setting lambda = 1e+06 \n trace of hat matrix differs from df by ",
            round(sum( 1/(1+lambda*d) )-df,6))
            break
            }
    }

    ### tmp <- sum(diag(X %*% solve(crossprod(X * weights, X) +
    ###                   lambda*dmat) %*% t(X * weights))) - df
    ### if (abs(tmp) > sqrt(.Machine$double.eps))
    ###   warning("trace of hat matrix is not equal df with difference", tmp)

    lambda
}

bns <- function(x, z = NULL, df = 4, knots = 20, differences = 2,
                xname = NULL, zname = NULL) {

    if (is.null(xname)) xname <- deparse(substitute(x))
    if (is.null(zname)) zname <- deparse(substitute(z))

    if (is.factor(x) || df <= 2)
        return(bols(x = x, z = z, xname = xname, zname = zname))

    if (!is.numeric(z) && (is.factor(z) && length(unique(z)) != 2))
        stop(sQuote("z"), " must be binary or numeric")

    if(is.factor(z) && length(unique(z)) == 2)
        ## FIXME is there a more elegant way to produce a binary with 0/1?
        z <- as.numeric(z[, drop=T]) - 1

    if (!differences %in% 1:3)
        stop(sQuote("differences"), " are not in 1:3")
    if (df < differences)
        stop(sQuote("df"), " is less than ", sQuote("differences"))
    if (length(unique(x)) < 6)
        stop(sQuote(xname), " has less than 6 unique values")

    if (length(unique(diff(knots))) > 1)
            warning("non-equidistant ", sQuote("knots"),
                    " might be inappropriate")

    if (length(knots) == 1) {
        knots <- seq(from = min(x, na.rm = TRUE),
                     to = max(x, na.rm = TRUE), length = knots + 2)
        knots <- knots[2:(length(knots) - 1)]
        #knots <- c(min(x)-sd(x),knots,max(x)+sd(x))
    }
    newX <- function(x, z = NULL) {
        epsilon <- diff(range(x)) / 10
        X <- ns(x, knots = knots, intercept = TRUE, Boundary.knots = c(min(x)-epsilon,
        max(x)+epsilon) )
        if (!is.null(z))
            X <- X * z
        return(X)
    }
    X <- newX(x, z)

    K <- diff(diag(ncol(X)), differences = differences)
    K <- crossprod(K, K)

    dpp <- function(weights) {

        lambda <- df2lambda(X, df = df, dmat = K, weights = weights)

        Xw <- X * weights
        XtX <- crossprod(Xw, X)
        Xsolve <- tcrossprod(solve(XtX + lambda * K), Xw)

        fitfun <- function(y) {
            coef <- Xsolve %*% y

            predictfun <- function(newdata = NULL) {
                if (is.null(newdata)) return(X %*% coef)
                nX <- newX(x = newdata[[xname]], z = newdata[[zname]])
                nX %*% coef
            }
            ret <- list(model = coef, predict = predictfun, fitted = function() X %*% coef)
            class(ret) <- c("basefit", "baselm")
            ret
        }
        ret <- list(fit = fitfun, hatmatrix = function() X %*% Xsolve)
        class(ret) <- "basisdpp"
        ret
    }
    attr(X, "dpp") <- dpp
    return(X)
}

bss <- function(x, df = 4, xname = NULL) {

    if (is.null(xname)) xname = deparse(substitute(x))
    if (is.factor(x) || df <= 2)
        return(bols(x = x, xname = xname))

    xs <- signif(x, 10)
    ux <- unique(sort(xs))

    dpp <- function(weights) {

        fitfun <- function(y) {
            object <- smoothbase(x = xs, ux = ux, y = y, w = weights, df = df)

            predictfun <- function(newdata = NULL) {
                if (is.null(newdata)) return(stats:::predict.smooth.spline.fit(object, x = xs)$y)
                stats:::predict.smooth.spline.fit(object, x = newdata[[xname]])$y
            }

            ret <- list(basemodel = object, predict = predictfun, fitted = predictfun)
            class(ret) <- c("basefit", "bssfit")
            ret
        }

        ret <- list(fit = fitfun, hatmatrix = function()
                                      hatMatTH(x = x, w = weights, df = df))
        class(ret) <- "basisdpp"
        ret
    }
    x <- matrix(x, nc = 1)
    attr(x, "dpp") <- dpp
    return(x)
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
    newX <- function(x, y, z) {
        Xx <- bs(x, knots = xknots, degree = degree, intercept = TRUE)
        Xy <- bs(y, knots = yknots, degree = degree, intercept = TRUE)
        X <- kronecker(Xx, matrix(1, nc = ncol(Xy))) * kronecker(matrix(1, nc = ncol(Xx)), Xy)
        if (!is.null(z))
            X <- X * z
        return(X)
    }
    X <- newX(x, y, z)

    xd <- length(xknots) + degree + 1
    yd <- length(yknots) + degree + 1

    Kx <- diff(diag(xd), differences = differences)
    Kx <- crossprod(Kx, Kx)
    Ky <- diff(diag(yd), differences = differences)
    Ky <- crossprod(Ky, Ky)
    K <- kronecker(Kx, diag(yd)) + kronecker(diag(xd), Ky)

    L <- 0
    if(center) {
        L <- eigen(K, symmetric=TRUE, EISPACK=TRUE)
        L$vectors <- L$vectors[,1:(ncol(X)-differences^2)]
        L$values <- sqrt(L$values[1:(ncol(X)-differences^2)])
        L <- L$vectors%*%diag(1/L$values)
        X <- X%*%L
        K <- diag(ncol(X))
    }

    dpp <- function(weights) {

        lambda <- df2lambda(X, df = df, dmat = K, weights = weights)
        Xw <- X * weights
        XtX <- crossprod(Xw, X)
        Xsolve <- tcrossprod(solve(XtX + lambda * K), Xw)

        fitfun <- function(y) {
            coef <- Xsolve %*% y

            predictfun <- function(newdata = NULL) {
                if (is.null(newdata)) return(X %*% coef)
                nX <- newX(x = newdata[[xname]], y = newdata[[yname]],
                           z = newdata[[zname]])
                if(center) {
                    nX <- nX%*%L
                }
                nX %*% coef
            }
            ret <- list(model = coef, predict = predictfun,
                        fitted = function() X %*% coef)
            class(ret) <- c("basefit", "baselm")
            ret
        }
        ret <- list(fit = fitfun, hatmatrix = function() X %*% Xsolve)
        class(ret) <- "basisdpp"
        ret
    }
    attr(X, "dpp") <- dpp
    return(X)
}

bols <- function(x, z = NULL, xname = NULL, zname = NULL, center = FALSE,
                 df = NULL, contrasts.arg = "contr.treatment") {

     if (is.null(xname)) xname = deparse(substitute(x))
     if (is.null(zname)) zname = deparse(substitute(z))

     cc <- complete_cases(x = x, z = z)

     newX <- function(x, z = NULL, na.rm = TRUE) {
         if (na.rm) {
             x <- x[cc]
             if (!is.null(z))
                 z <- z[cc]
         }
         
         if (is.factor(x)) {
             X <- model.matrix(~ x, contrasts.arg = list(x = contrasts.arg))
         } else {
             X <- model.matrix(~ x)
         }

         if (center)
            X <- X[, -1, drop = FALSE]

         if (any(!cc) & !na.rm) {
             Xtmp <- matrix(NA, ncol = ncol(X), nrow = length(cc))
             Xtmp[cc,] <- X
             X <- Xtmp
         }
         if (!is.null(z)) X <- X * z
         X
     }
     X <- newX(x, z)
     Xna <- X
     if (any(!cc))
         Xna <- newX(x, z, na.rm = FALSE)

     K <- diag(ncol(X))

     dpp <- function(weights) {

         if (any(!cc)) weights <- weights[cc]
         Xw <- X * weights
         XtX <- crossprod(Xw, X)

         if (is.null(df) || df >= ncol(K)) {
             Xsolve <- tcrossprod(solve(crossprod(Xw, X)), Xw)
         } else {
             lambda <- df2lambda(X, df = df, dmat = K, weights = weights)
             Xsolve <- tcrossprod(solve(XtX + lambda * K), Xw)
         }

         fitfun <- function(y) {

             if (any(!cc)) y <- y[cc]
             coef <- Xsolve %*% y

             predictfun <- function(newdata = NULL) {
                 if (is.null(newdata)) return(Xna %*% coef)
                 nX <- newX(x = newdata[[xname]], z = newdata[[zname]],
                            na.rm = FALSE)
                 nX %*% coef
             }
             ret <- list(model = coef, predict = predictfun,
                         fitted = function() Xna %*% coef)
             class(ret) <- c("basefit", "baselm")
             ret
         }
         ret <- list(fit = fitfun, hatmatrix = function() X %*% Xsolve)
         class(ret) <- "basisdpp"
         ret
     }
     attr(X, "dpp") <- dpp
     return(X)
}

brandom <- function(x, z = NULL, df = 4, xname = NULL,
                    zname = NULL) {

    if (is.null(xname)) xname = deparse(substitute(x))
    if (is.null(zname)) zname = deparse(substitute(z))

    if (!is.numeric(z) && (is.factor(z) && length(unique(z)) != 2))
        stop(sQuote("z"), " must be binary or numeric")

    if(is.factor(z) && length(unique(z)) == 2)
        ## FIXME is there a more elegant way to produce a binary with 0/1?
        z <- as.numeric(z[, drop=T]) - 1

    newX <- function(x, z = NULL) {
        if (!is.factor(x)) stop(sQuote("x"), " is not a factor")
        X <- model.matrix(~ x - 1)
        if (!is.null(z))
            X <- X * z
        return(X)
    }
    X <- newX(x, z)

    K <- diag(ncol(X))

    dpp <- function(weights) {

        lambda <- df2lambda(X, df = df, dmat = K, weights = weights)

        Xw <- X * weights
        XtX <- crossprod(Xw, X)
        ### XtX and K are diagonal matrices
        Xsolve <- tcrossprod(solve(XtX + lambda * K), Xw)

        fitfun <- function(y) {
            coef <- Xsolve %*% y

            predictfun <- function(newdata = NULL) {
                if (is.null(newdata)) return(X %*% coef)
                nX <- newX(x = newdata[[xname]], z = newdata[[zname]])
                nX %*% coef
            }
            ret <- list(model = coef, predict = predictfun, fitted = function() X %*% coef)
            class(ret) <- c("basefit", "baselm")
            ret
        }
        ret <- list(fit = fitfun, hatmatrix = function() X %*% Xsolve)
        class(ret) <- "basisdpp"
        ret
    }
    attr(X, "dpp") <- dpp
    return(X)
}

btree <- function(..., tree_controls = ctree_control(stump = TRUE,
    mincriterion = 0), xname = NULL) {

    x <- as.data.frame(list(...))

    if (is.null(xname)) {
        cl <- as.list(match.call(expand.dots = FALSE))[2][[1]]
        xname <- sapply(cl, function(x) as.character(x))
        colnames(x) <- xname
    } else {
        colnames(x) <- xname
    }

    X <- matrix(numeric(nrow(x)))

    dpp <- function(weights) {

        ### construct design matrix etc.
        y <- vector(length = nrow(x), mode = "numeric")
        ### name for working response (different from any x)
        rname <- paste("R", paste(colnames(x), collapse = "_"), sep = "_")
        fm <- as.formula(paste(rname, " ~ ", paste(xname, collapse = "+")))
        df <- x
        df[[rname]] <- y
        object <- party:::ctreedpp(fm, data = df)
        fitmem <- ctree_memory(object, TRUE)
        where <- rep.int(0, nrow(x))
        storage.mode(where) <- "integer"
        storage.mode(weights) <- "double"

        fitfun <- function(y) {

            .Call("R_modify_response", as.double(y), object@responses,
                 PACKAGE = "party")
            tree <- .Call("R_TreeGrow", object, weights, fitmem, tree_controls,
                          where, PACKAGE = "party")
            .Call("R_remove_weights", tree, package = "party")

            predictfun <- function(newdata = NULL) {
                if (is.null(newdata)) {
                    wh <- .Call("R_get_nodeID", tree, object@inputs, 0.0, PACKAGE = "party")
                    return(unlist(.Call("R_getpredictions", tree, wh, PACKAGE = "party")))
                }
                newinp <- party:::newinputs(object, newdata)
                wh <- .Call("R_get_nodeID", tree, newinp, 0.0,
                        PACKAGE = "party")
                unlist(.Call("R_getpredictions", tree, wh, PACKAGE = "party"))
            }
            ret <- list(model = tree, predict = predictfun, fitted = predictfun)
            class(ret) <- "basefit"
            ret
        }
        ret <- list(fit = fitfun, hatmatrix = function() NA)
        class(ret) <- "basisdpp"
        ret
    }
    attr(X, "dpp") <- dpp
    return(X)
}
