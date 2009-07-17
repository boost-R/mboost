
df2lambda <- function(X, df = 4, lambda = NULL, dmat = diag(ncol(X)), weights) {

    if (!is.null(df))
        if (df > ncol(X)) return(c(df = df, lambda = 0))

    # Demmler-Reinsch Orthogonalization (cf. Ruppert et al., 2003, 
    # Semiparametric Regression, Appendix B.1.1).

    ### option
    A <- crossprod(X * weights, X) + dmat * 10e-10
    Rm <- solve(chol(A))

    decomp <- svd(crossprod(Rm, dmat) %*% Rm)
    d <- decomp$d[decomp$d > .Machine$double.eps]

    if (!is.null(lambda)) 
        return(c(df = sum(1 / (1 + lambda * d)), lambda = lambda))
    if (df >= length(d)) return(c(df = df, lambda = 0))

    # search for appropriate lambda using uniroot
    df2l <- function(lambda)
        sum(1/(1 + lambda * d)) - df

    ### option
    if (df2l(1e+10) > 0) return(c(df = df, lambda = 1e+10))
    return(c(df = df, 
             lambda = uniroot(df2l, c(0, 1e+10), 
                              tol = sqrt(.Machine$double.eps))$root))
}

hyper_ols <- function(df = NULL, lambda = NULL, intercept = TRUE, 
                      contrasts.arg = "contr.treatment") 
    list(pen = !is.null(df) || !is.null(lambda),
         df = df, lambda = lambda,
         intercept = intercept, 
         contrasts.arg = contrasts.arg)

X_ols <- function(mf, vary, args) {

    if (isMATRIX(mf)) {
        X <- mf
        contr <- NULL
    } else {
        ### set up model matrix
        fm <- paste("~ ", paste(colnames(mf)[colnames(mf) != vary], 
                    collapse = "+"), sep = "")
        if (!args$intercept)
            fm <- paste(fm, "-1", collapse = "")
        X <- model.matrix(as.formula(fm), data = mf, contrasts.arg = args$contrasts.arg)
        contr <- attr(X, "contrasts")
        if (vary != "") {
            ### <FIXME> is this really what we want?
            z <- model.matrix(as.formula(paste("~", vary, collapse = "")), data = mf)[,2]
            X <- X * z
            ### </FIXME>
        }
    }
    K <- NULL
    if (args$pen) {
        ### set up penalty matrix
        ANCOVA <- !is.null(contr)
        if (ANCOVA) { 
            diag <- Diagonal
            X <- Matrix(X)
        }
        K <- diag(ncol(X))
        ### for ordered factors use difference penalty
        if (ANCOVA && any(sapply(mf[, names(contr), drop = FALSE], is.ordered))) {
            K <- diff(diag(ncol(X)), differences = 2)
            K <- crossprod(K)
        }
    }
    list(X = X, K = K)
}

hyper_bbs <- function(mf, vary, knots = 20, degree = 3, differences = 2, df = 4, lambda = NULL,
                      center = FALSE) {

    knotf <- function(x, knots) {	
        boundary.knots <- range(x, na.rm = TRUE)
        if (length(knots) == 1) {
            knots <- seq(from = boundary.knots[1], 
                         to = boundary.knots[2], length = knots + 2)
            knots <- knots[2:(length(knots) - 1)]
        }
        list(knots = knots, boundary.knots = boundary.knots)
    }
    nm <- colnames(mf)[colnames(mf) != vary]
    if (is.list(knots)) stopifnot(all(names(knots) %in% nm))
    ret <- vector(mode = "list", length = length(nm))
    names(ret) <- nm
    for (n in nm)
        ret[[n]] <- knotf(mf[[n]], if (is.list(knots)) knots[[n]] else knots)
    list(knots = ret, degree = degree, differences = differences, pen = TRUE,
         df = df, lambda = lambda, center = center)
}

X_bbs <- function(mf, vary, args) {

    mm <- lapply(which(colnames(mf) != vary), function(i) {
        X <- bs(mf[[i]], knots = args$knots[[i]]$knots, degree = args$degree,
           Boundary.knots = args$knots[[i]]$boundary.knots, intercept = TRUE)
        class(X) <- "matrix"
        return(X)
    })
    ### options
    MATRIX <- any(sapply(mm, dim) > c(500, 50)) || (length(mm) > 1)
    if (MATRIX) {
        diag <- Diagonal
        for (i in 1:length(mm)) mm[[i]] <- Matrix(mm[[i]])
    }
    if (length(mm) == 1) {
        X <- mm[[1]]
        K <- diff(diag(ncol(X)), differences = args$differences)
        K <- crossprod(K)
    }
    if (length(mm) == 2) {
        X <- kronecker(mm[[1]], matrix(1, nc = ncol(mm[[2]]))) * 
             kronecker(matrix(1, nc = ncol(mm[[1]])), mm[[2]])
        Kx <- diff(diag(ncol(mm[[1]])), differences = args$differences)
        Kx <- crossprod(Kx)
        Ky <- diff(diag(ncol(mm[[2]])), differences = args$differences)
        Ky <- crossprod(Ky)
        K <- kronecker(Kx, diag(ncol(mm[[2]]))) + 
             kronecker(diag(ncol(mm[[1]])), Ky)
    }
    ### <FIXME>
    if (vary != "") {
        z <- model.matrix(as.formula(paste("~", vary, collapse = "")), data = mf)[,2]
        X <- X * z
    }
    ### </FIXME>
    if (args$center) {
        L <- eigen(K, symmetric = TRUE, EISPACK = TRUE)
        L$vectors <- L$vectors[,1:(ncol(X) - args$differences^2)]
        L$values <- sqrt(L$values[1:(ncol(X) - args$differences^2)])
        L <- L$vectors %*% (diag(length(L$values)) * (1/L$values))
        X <- as(X %*% L, "matrix")
        K <- as(diag(ncol(X)), "matrix")
    }
    return(list(X = X, K = K))
}

bols3 <- function(..., z = NULL, index = NULL, intercept = TRUE, df = NULL, lambda = NULL,
                  contrasts.arg = "contr.treatment") {

    mf <- list(...)
    if (length(mf) == 1 && (isMATRIX(mf[[1]]) || is.data.frame(mf[[1]]))) {
        mf <- mf[[1]]
        ### spline bases should be matrices
        if (isMATRIX(mf) && !is(mf, "Matrix"))
            class(mf) <- "matrix"
    } else {
        mf <- as.data.frame(mf)
        cl <- as.list(match.call(expand.dots = FALSE))[2][[1]]
        colnames(mf) <- sapply(cl, function(x) as.character(x))
    }
    vary <- ""
    if (!is.null(z)) {
        stopifnot(is.data.frame(mf))
        stopifnot(is.numeric(z) || (is.factor(z) && nlevels(z) == 2))
        mf <- cbind(mf, z)
        colnames(mf)[ncol(mf)] <- vary <- deparse(substitute(z))
    }

    CC <- all(Complete.cases(mf))
    ### option
    DOINDEX <- is.data.frame(mf) && (nrow(mf) > 10000 || is.factor(mf[[1]]))
    if (is.null(index)) {
        if (!CC || DOINDEX) {
            index <- get_index(mf)
            mf <- mf[index[[1]],,drop = FALSE]
            index <- index[[2]]
        }
    }

    ret <- list(model.frame = function() 
                    if (is.null(index)) return(mf) else return(mf[index,,drop = FALSE]),
                get_names = function() colnames(mf),
                get_vary = function() vary,
                set_names = function(value) attr(mf, "names") <<- value)
    class(ret) <- "blg"

    ret$dpp <- bl_lin(mf, vary, index = index, Xfun = X_ols, args = hyper_ols(
                      df = df, lambda = lambda, 
                      intercept = intercept, contrasts.arg = contrasts.arg))
    return(ret)
}

bbs3 <- function(..., z = NULL, index = NULL, knots = 20, degree = 3, 
                 differences = 2, df = 4, lambda = NULL, center = FALSE) {

    mf <- list(...)
    if (length(mf) == 1 && (is.matrix(mf[[1]]) || is.data.frame(mf[[1]]))) {
        mf <- as.data.frame(mf[[1]])
    } else {
        mf <- as.data.frame(mf)
        cl <- as.list(match.call(expand.dots = FALSE))[2][[1]]
        colnames(mf) <- sapply(cl, function(x) as.character(x))
    }
    stopifnot(is.data.frame(mf))
    stopifnot(all(sapply(mf, is.numeric)))
    vary <- ""
    if (!is.null(z)) {
        stopifnot(is.numeric(z) || (is.factor(z) && nlevels(z) == 2))
        mf <- cbind(mf, z)
        colnames(mf)[ncol(mf)] <- vary <- deparse(substitute(z))
    }

    CC <- all(Complete.cases(mf))
    ### option
    DOINDEX <- (nrow(mf) > 10000)
    if (is.null(index)) {
        if (!CC || DOINDEX) {
            index <- get_index(mf)
            mf <- mf[index[[1]],,drop = FALSE]
            index <- index[[2]]
        }
    }

    ret <- list(model.frame = function() 
                    if (is.null(index)) return(mf) else return(mf[index,,drop = FALSE]),
                get_vary = function() vary,
                get_names = function() colnames(mf),
                set_names = function(value) attr(mf, "names") <<- value)
    class(ret) <- "blg"

    ret$dpp <- bl_lin(mf, vary, index = index, Xfun = X_bbs, 
                      args = hyper_bbs(mf, vary, knots = knots,
                      degree = degree, differences = differences, 
                      df = df, lambda = lambda, center = center))
    return(ret)
}

bl_lin <- function(mf, vary, index = NULL, Xfun, args) {

    newX <- function(newdata = NULL) {
        if (!is.null(newdata)) {
            stopifnot(all(names(newdata) == names(mf)))
            mf <- newdata[,colnames(mf),drop = FALSE]
        }
        return(Xfun(mf, vary, args))
    }
    X <- newX()
    K <- X$K
    X <- X$X

    dpp <- function(weights) {

        weights[!Complete.cases(mf)] <- 0
        w <- weights
        if (!is.null(index)) w <- as.vector(tapply(weights, index, sum))
        XtX <- crossprod(X * w, X)
        if (args$pen) {
            if (is.null(args$lambda)) {
                lambda <- df2lambda(X, df = args$df, dmat = K, weights = w)["lambda"]
            } else {
                lambda <- args$lambda
            }
            XtX <- XtX + lambda * K
        }

        if (is(X, "Matrix")) {
            ### chol benutzen
            XtXC <- Cholesky(forceSymmetric(XtX))
            mysolve <- function(y) solve(XtXC, crossprod(X, y))
        } else {
            mysolve <- function(y) 
                .Call("La_dgesv", XtX, crossprod(X, y), .Machine$double.eps, 
                      PACKAGE = "base")
        }	

        fit <- function(y) {
            if (!is.null(index)) {
                y <- as.vector(tapply(weights * y, index, sum))
            } else {
                y <- y * weights
            }
            coef <- mysolve(y)
            ret <- list(model = coef, 
                        fitted = function() {
                            ret <- as.vector(X %*% coef)
                            if (is.null(index)) return(ret)
                            return(ret[index])
                        })
            class(ret) <- c("bm_lin", "bm")
            ret
        }

        ### check for n
        hatvalues <- function() {
            ret <- as.matrix(tcrossprod(X %*% solve(XtX), X * w))
            if (is.null(index)) return(ret)
            return(ret[index,index])
        }

        df <- function() {
            if (args$pen) 
                return(df2lambda(X, df = NULL, lambda = lambda, 
                                 dmat = K, weights = w))
            return(ncol(X))
        }

        predict <- function(bm, newdata = NULL, aggregate = c("sum", "cumsum", "none")) {
            cf <- sapply(bm, coef)
            if(!is.null(newdata)) {
                index <- NULL
                nm <- colnames(mf)
                newdata <- newdata[,nm, drop = FALSE]
                ### option
                if (nrow(newdata) > 10000) {
                    index <- get_index(newdata)
                    newdata <- newdata[index[[1]],,drop = FALSE]
                    index <- index[[2]]
                }
                X <- newX(newdata)$X
            }
            aggregate <- match.arg(aggregate)
            pr <- switch(aggregate, "sum" = 
                as(X %*% rowSums(cf), "matrix"),
            "cumsum" = {
                M <- triu(crossprod(Matrix(1, nc = ncol(cf))))
                as(X %*% (cf %*% M), "matrix")
            },
            "none" = as(X %*% cf, "matrix"))
            if (is.null(index)) return(pr[,,drop = FALSE])
            return(pr[index,,drop = FALSE])
        }

        ret <- list(fit = fit, hatvalues = hatvalues, 
                    predict = predict, df = df, 
                    Xnames = colnames(X))
        class(ret) <- c("bl_lin", "bl")
        return(ret)

    }
    return(dpp)
}

bspatial3 <- function(...) {
    cl <- match.call()
    cl[[1L]] <- as.name("bbs3")
    eval(cl, parent.frame())
}

brandom3 <- function(..., df = 4) {
    cl <- match.call()
    if (is.null(cl$df)) cl$df <- df
    cl$intercept <- FALSE
    cl[[1L]] <- as.name("bols3")
    eval(cl, parent.frame())
}


names.blg <- function(x)
    x$get_names()

model.frame.blg <- function(formula)
    formula$model.frame()

coef.bm_lin <- function(object) {
    ret <- as.vector(object$model)
    names(ret) <- object$Xnames
    ret
}

fitted.bm <- function(object)
    object$fitted()

hatvalues.bl_lin <- function(model)
    model$hatvalues()

dpp <- function(object, weights)
    UseMethod("dpp", object)

fit <- function(object, y)
    UseMethod("fit", object)

dpp.blg <- function(object, weights)
    object$dpp(weights)

fit.bl <- function(object, y)
    object$fit(y)

