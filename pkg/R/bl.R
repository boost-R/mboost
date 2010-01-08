
### compute Ridge shrinkage parameter lambda from df
### or the other way round
df2lambda <- function(X, df = 4, lambda = NULL, dmat = diag(ncol(X)), weights) {


    stopifnot(xor(is.null(df), is.null(lambda)))
    if (!is.null(df))
        if (df >= ncol(X)) return(c(df = df, lambda = 0))
    if (!is.null(lambda))
        if (lambda == 0) return(c(df = ncol(X), lambda = 0))

    # Demmler-Reinsch Orthogonalization (cf. Ruppert et al., 2003,
    # Semiparametric Regression, Appendix B.1.1).

    ### option
    A <- crossprod(X * weights, X) + dmat * 10e-10
    Rm <- solve(chol(A))

    decomp <- svd(crossprod(Rm, dmat) %*% Rm)
    d <- decomp$d[decomp$d > sqrt(.Machine$double.eps)]

    if (!is.null(lambda))
        return(c(df = sum(1 / (1 + lambda * d)), lambda = lambda))
    if (df >= length(d)) return(c(df = df, lambda = 0))

    ### <FIXME> Buja definition of df
    stopifnot(options("mboost_dftraceS")[[1]])
    ### </FIXME>

    # search for appropriate lambda using uniroot
    df2l <- function(lambda)
        sum(1/(1 + lambda * d)) - df

    ### option
    if (df2l(1e+10) > 0) return(c(df = df, lambda = 1e+10))
    return(c(df = df,
             lambda = uniroot(df2l, c(0, 1e+10),
                              tol = sqrt(.Machine$double.eps))$root))
}

### hyper parameters for ols baselearner
hyper_ols <- function(df = NULL, lambda = 0, intercept = TRUE,
                      contrasts.arg = "contr.treatment")
    list(df = df, lambda = lambda,
         intercept = intercept,
         contrasts.arg = contrasts.arg)

### model.matrix for ols baselearner
X_ols <- function(mf, vary, args) {

    if (isMATRIX(mf)) {
        X <- mf
        contr <- NULL
    } else {
        ### set up model matrix
        fm <- paste("~ ", paste(colnames(mf)[colnames(mf) != vary],
                    collapse = "+"), sep = "")
        X <- model.matrix(as.formula(fm), data = mf, contrasts.arg = args$contrasts.arg)
        if (!args$intercept)
            X <- X[ , -1, drop=FALSE]
        contr <- attr(X, "contrasts")
        ### <FIXME>
        if (vary != "") {
            by <- model.matrix(as.formula(paste("~", vary, collapse = "")), data = mf)[,2]
            X <- X * by
        }
        ### </FIXME>
    }
    ### <FIXME> penalize intercepts???
    ### set up penalty matrix
    ANOVA <- (!is.null(contr) && (length(contr) == 1)) && (ncol(mf) == 1)
    K <- diag(ncol(X))
    ### for ordered factors use difference penalty
    if (ANOVA && any(sapply(mf[, names(contr), drop = FALSE], is.ordered))) {
        K <- diff(diag(ncol(X) + 1), differences = 1)[, -1, drop = FALSE]
        K <- crossprod(K)
    }
    ### </FIXME>
    list(X = X, K = K)
}

### hyper parameters for P-splines baselearner (including tensor product P-splines)
hyper_bbs <- function(mf, vary, knots = 20, degree = 3, differences = 2, df = 4,
                      lambda = NULL, center = FALSE) {

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
    if (is.list(knots)) if(!all(names(knots) %in% nm)) 
        stop("variable names and knot names must be the same")
    ret <- vector(mode = "list", length = length(nm))
    names(ret) <- nm
    for (n in nm)
        ret[[n]] <- knotf(mf[[n]], if (is.list(knots)) knots[[n]] else knots)
    list(knots = ret, degree = degree, differences = differences,
         df = df, lambda = lambda, center = center)
}

### model.matrix for P-splines baselearner (including tensor product P-splines)
X_bbs <- function(mf, vary, args) {

    stopifnot(is.data.frame(mf))
    mm <- lapply(which(colnames(mf) != vary), function(i) {
        X <- bs(mf[[i]], knots = args$knots[[i]]$knots, degree = args$degree,
           Boundary.knots = args$knots[[i]]$boundary.knots, intercept = TRUE)
        class(X) <- "matrix"
        return(X)
    })
    ### options
    MATRIX <- any(sapply(mm, dim) > c(500, 50)) || (length(mm) > 1)
    MATRIX <- MATRIX && options("mboost_useMatrix")$mboost_useMatrix
    if (MATRIX) {
        diag <- Diagonal
        for (i in 1:length(mm)) mm[[i]] <- Matrix(mm[[i]])
    }
    if (length(mm) == 1) {
        X <- mm[[1]]
        ### <FIXME>
        if (vary != "") {
            by <- model.matrix(as.formula(paste("~", vary, collapse = "")), data = mf)[,2]
            X <- X * by
        }
        ### </FIXME>
        K <- diff(diag(ncol(X)), differences = args$differences)
        if (args$center) {
            X <- tcrossprod(X, K) %*% solve(tcrossprod(K))
            K <- diag(ncol(X))
        } else {
            K <- crossprod(K)
        }
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
        ### <FIXME>
        if (vary != "") {
            by <- model.matrix(as.formula(paste("~", vary, collapse = "")), data = mf)[,2]
            X <- X * by
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
    }
    if (length(mm) > 2)
        stop("not possible to specify more than two variables in ", 
             sQuote("..."), " argument of smooth base-learners")
    return(list(X = X, K = K))
}

### Linear baselearner, potentially Ridge-penalized (but not by default)
bols <- function(..., by = NULL, index = NULL, intercept = TRUE, df = NULL, 
                 lambda = 0, contrasts.arg = "contr.treatment") {

    if (!is.null(df)) lambda <- NULL

    cll <- match.call()
    cll[[1]] <- as.name("bols")
    cll <- deparse(cll)
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
    if (!is.null(by)){
        stopifnot(is.data.frame(mf))
        stopifnot(is.numeric(by) || (is.factor(by) && nlevels(by) == 2))
        mf <- cbind(mf, by)
        colnames(mf)[ncol(mf)] <- vary <- deparse(substitute(by))
    }

    CC <- all(Complete.cases(mf))
    ### option
    DOINDEX <- is.data.frame(mf) && 
        (nrow(mf) > options("mboost_indexmin")[[1]] || is.factor(mf[[1]]))
    if (is.null(index)) {
        ### try to remove duplicated observations or
        ### observations with missings
        if (!CC || DOINDEX) {
            index <- get_index(mf)
            mf <- mf[index[[1]],,drop = FALSE]
            index <- index[[2]]
        }
    }

    ret <- list(model.frame = function()
                    if (is.null(index)) return(mf) else return(mf[index,,drop = FALSE]),
                get_call = function() cll,
                get_data = function() mf,
                get_index = function() index,
                get_names = function() colnames(mf),
                get_vary = function() vary,
                set_names = function(value) {
                    attr(mf, "names") <<- value
                    cll <<- paste("bols", "(", paste(colnames(mf), 
                        collapse = ", "), ")", sep = "")
                })
    class(ret) <- "blg"

    ret$dpp <- bl_lin(ret, Xfun = X_ols, args = hyper_ols(
                      df = df, lambda = lambda,
                      intercept = intercept, contrasts.arg = contrasts.arg))
    return(ret)
}

### P-spline (and tensor-product spline) baselearner
bbs <- function(..., by = NULL, index = NULL, knots = 20, degree = 3,
                 differences = 2, df = 4, lambda = NULL, center = FALSE) {

    if (!is.null(lambda)) df <- NULL

    cll <- match.call()
    cll[[1]] <- as.name("bbs")
    cll <- deparse(cll)

    mf <- list(...)
    if (length(mf) == 1 && (is.matrix(mf[[1]]) || is.data.frame(mf[[1]]))) {
        mf <- as.data.frame(mf[[1]])
    } else {
        mf <- as.data.frame(mf)
        cl <- as.list(match.call(expand.dots = FALSE))[2][[1]]
        colnames(mf) <- sapply(cl, function(x) deparse(x))
    }
    stopifnot(is.data.frame(mf))
    if(!(all(sapply(mf, is.numeric)))) {
        if (ncol(mf) == 1) return(bols(..., by = by, index = index))
        stop("cannot compute bbs for non-numeric variables")
    }
    ### use bols when appropriate
    if (!is.null(df) & !center) {
        if (df <= ncol(mf))
            return(bols(..., by = by, index = index))
    }
    vary <- ""
    if (!is.null(by)){
        stopifnot(is.numeric(by) || (is.factor(by) && nlevels(by) == 2))
        mf <- cbind(mf, by)
        colnames(mf)[ncol(mf)] <- vary <- deparse(substitute(by))
    }

    CC <- all(Complete.cases(mf))
    ### option
    DOINDEX <- (nrow(mf) > options("mboost_indexmin")[[1]])
    if (is.null(index)) {
        if (!CC || DOINDEX) {
            index <- get_index(mf)
            mf <- mf[index[[1]],,drop = FALSE]
            index <- index[[2]]
        }
    }

    ret <- list(model.frame = function()
                    if (is.null(index)) return(mf) else return(mf[index,,drop = FALSE]),
                get_call = function() cll,
                get_data = function() mf,
                get_index = function() index,
                get_vary = function() vary,
                get_names = function() colnames(mf),
                set_names = function(value) {
                    attr(mf, "names") <<- value
                    cll <<- paste("bbs", "(", paste(colnames(mf), 
                        collapse = ", "), ")", sep = "")
                })
    class(ret) <- "blg"

    ret$dpp <- bl_lin(ret, Xfun = X_bbs,
                      args = hyper_bbs(mf, vary, knots = knots,
                      degree = degree, differences = differences,
                      df = df, lambda = lambda, center = center))
    return(ret)
}

### workhorse for fitting (Ridge-penalized) baselearners
bl_lin <- function(blg, Xfun, args) {

    mf <- blg$get_data()
    index <- blg$get_index()
    vary <- blg$get_vary()

    newX <- function(newdata = NULL) {
        if (!is.null(newdata)) {
            stopifnot(all(names(newdata) == names(blg)))
            stopifnot(all(class(newdata) == class(mf)))
            mf <- newdata[,names(blg),drop = FALSE]
        }
        return(Xfun(mf, vary, args))
    }
    X <- newX()
    K <- X$K
    X <- X$X

    dpp <- function(weights) {

        weights[!Complete.cases(mf)] <- 0
        w <- weights
        if (!is.null(index)) 
            w <- .Call("R_ysum", as.double(weights), as.integer(index), PACKAGE = "mboost")
        XtX <- crossprod(X * w, X)
        lambdadf <- df2lambda(X, df = args$df, lambda = args$lambda, 
                              dmat = K, weights = w)
        lambda <- lambdadf["lambda"]
        XtX <- XtX + lambda * K

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
                y <- .Call("R_ysum", as.double(weights * y), as.integer(index), PACKAGE = "mboost")
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

        ### <FIXME> check for large n, option?
        hatvalues <- function() {
            ret <- as.matrix(tcrossprod(X %*% solve(XtX), X * w))
            if (is.null(index)) return(ret)
            return(ret[index,index])
        }
        ### </FIXME>

        ### actually used degrees of freedom (trace of hat matrix)
        df <- function() lambdadf

        ### prepare for computing predictions
        predict <- function(bm, newdata = NULL, aggregate = c("sum", "cumsum", "none")) {
            cf <- sapply(bm, coef)
            if (!is.matrix(cf)) cf <- matrix(cf, nrow = 1)
            if(!is.null(newdata)) {
                index <- NULL
                nm <- names(blg)
                newdata <- newdata[,nm, drop = FALSE]
                ### option
                if (nrow(newdata) > options("mboost_indexmin")[[1]]) {
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
                as(X %*% .Call("R_mcumsum", as(cf, "matrix")), "matrix")
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

### tensor-product spline baselearner
bspatial <- function(...) {
    cl <- match.call()
    cl[[1L]] <- as.name("bbs")
    eval(cl, parent.frame())
}

### random-effects (Ridge-penalized ANOVA) baselearner
brandom <- function(..., df = 4) {
    cl <- match.call()
    if (is.null(cl$df)) cl$df <- df
    cl$intercept <- FALSE
    cl[[1L]] <- as.name("bols")
    eval(cl, parent.frame())
}

### extract variables names from baselearner
names.blg <- function(x)
    x$get_names()

### extract data from baselearner
model.frame.blg <- function(formula, ...)
    formula$model.frame(...)

### extract coefficients
coef.bm_lin <- function(object, ...) {
    ret <- as.vector(object$model)
    names(ret) <- object$Xnames
    ret
}

### extract fitted values
fitted.bm <- function(object)
    object$fitted()

### extract hatmatrix
hatvalues.bl_lin <- function(model)
    model$hatvalues()

### data preprocessing (plug in weights)
dpp <- function(object, weights)
    UseMethod("dpp", object)

dpp.blg <- function(object, weights)
    object$dpp(weights)

### actually fit a baselearner to response y
fit <- function(object, y)
    UseMethod("fit", object)

fit.bl <- function(object, y)
    object$fit(y)

"%+%" <- function(bl1, bl2) {

    if (is.list(bl1) && !inherits(bl1, "blg"))
        return(lapply(bl1, "%+%", bl2 = bl2))

    if (is.list(bl2) && !inherits(bl2, "blg"))
        return(lapply(bl2, "%+%", bl1 = bl1))

    cll <- paste(deparse(bl1$get_call()), "%+%", 
                 deparse(bl2$get_call()), collapse = "")
    stopifnot(inherits(bl1, "blg"))
    stopifnot(inherits(bl2, "blg"))

    mf <- cbind(model.frame(bl1), model.frame(bl2))
    index1 <- bl1$get_index()
    index2 <- bl2$get_index()
    if (is.null(index1)) index1 <- 1:nrow(mf)
    if (is.null(index2)) index2 <- 1:nrow(mf)

    mfindex <- cbind(index1, index2)
    index <- NULL
    
    CC <- all(Complete.cases(mf))
    ### option
    DOINDEX <- (nrow(mf) > options("mboost_indexmin")[[1]])
    if (is.null(index)) {
        if (!CC || DOINDEX) {
            index <- get_index(mfindex)
            mf <- mf[index[[1]],,drop = FALSE]
            index <- index[[2]]
        }
    }

    vary <- ""

    ret <- list(model.frame = function()
                    if (is.null(index)) return(mf) else return(mf[index,,drop = FALSE]),
                get_call = function() cll,
                get_data = function() mf,
                get_index = function() index,
                get_vary = function() vary,
                get_names = function() colnames(mf),
                set_names = function(value) attr(mf, "names") <<- value)
    class(ret) <- "blg"

    args1 <- environment(bl1$dpp)$args
    args2 <- environment(bl2$dpp)$args
    l1 <- args1$lambda        
    l2 <- args2$lambda
    if (!is.null(l1) && !is.null(l2)) {
        args <- list(lambda = 1, df = NULL)
    } else {
        args <- list(lambda = NULL,
            df = ifelse(is.null(args1$df), 0, args1$df) +
                 ifelse(is.null(args2$df), 0, args2$df))
    }

    Xfun <- function(mf, vary, args) {

        newX1 <- environment(bl1$dpp)$newX
        newX2 <- environment(bl2$dpp)$newX

        X1 <- newX1(mf[, bl1$get_names(), drop = FALSE])
        K1 <- X1$K
        if (!is.null(l1)) K1 <- l1 * K1
        X1 <- X1$X

        X2 <- newX2(mf[, bl2$get_names(), drop = FALSE])
        K2 <- X2$K
        if (!is.null(l2)) K2 <- l2 * K2
        X2 <- X2$X

        K <- matrix(0, ncol = ncol(K1) + ncol(K2),
                    nrow = nrow(K1) + nrow(K2))
        K[1:nrow(K1), 1:ncol(K1)] <- K1
        K[-(1:nrow(K1)), -(1:ncol(K1))] <- K2
        list(X = cbind(X1, X2), K = K)
    }

    ret$dpp <- bl_lin(ret, Xfun = Xfun, args = args)

    return(ret)
}

"%X%" <- function(bl1, bl2) {

    if (is.list(bl1) && !inherits(bl1, "blg"))
        return(lapply(bl1, "%X%", bl2 = bl2))

    if (is.list(bl2) && !inherits(bl2, "blg"))
        return(lapply(bl2, "%X%", bl1 = bl1))

    cll <- paste(deparse(bl1$get_call()), "%X%", 
                 deparse(bl2$get_call()), collapse = "")
    stopifnot(inherits(bl1, "blg"))
    stopifnot(inherits(bl2, "blg"))

    mf <- cbind(model.frame(bl1), model.frame(bl2))
    index1 <- bl1$get_index()
    index2 <- bl2$get_index()
    if (is.null(index1)) index1 <- 1:nrow(mf)
    if (is.null(index2)) index2 <- 1:nrow(mf)

    mfindex <- cbind(index1, index2)
    index <- NULL
    
    CC <- all(Complete.cases(mf))
    ### option
    DOINDEX <- (nrow(mf) > options("mboost_indexmin")[[1]])
    if (is.null(index)) {
        if (!CC || DOINDEX) {
            index <- get_index(mfindex)
            mf <- mf[index[[1]],,drop = FALSE]
            index <- index[[2]]
        }
    }

    vary <- ""

    ret <- list(model.frame = function()
                    if (is.null(index)) return(mf) else return(mf[index,,drop = FALSE]),
                get_call = function() cll,
                get_data = function() mf,
                get_index = function() index,
                get_vary = function() vary,
                get_names = function() colnames(mf),
                set_names = function(value) attr(mf, "names") <<- value)
    class(ret) <- "blg"

    args1 <- environment(bl1$dpp)$args
    args2 <- environment(bl2$dpp)$args
    l1 <- args1$lambda        
    l2 <- args2$lambda
    if (!is.null(l1) && !is.null(l2)) {
        args <- list(lambda = 1, df = NULL)
    } else {
        args <- list(lambda = NULL,
            df = ifelse(is.null(args1$df), 1, args1$df) * 
                 ifelse(is.null(args2$df), 1, args2$df))
    }

    Xfun <- function(mf, vary, args) {

        newX1 <- environment(bl1$dpp)$newX
        newX2 <- environment(bl2$dpp)$newX

        X1 <- newX1(mf[, bl1$get_names(), drop = FALSE])
        K1 <- X1$K
        X1 <- X1$X
        if (!is.null(l1)) K1 <- l1 * K1
        if (options("mboost_useMatrix")$mboost_useMatrix) {
            X1 <- Matrix(X1)
            K1 <- Matrix(K1)
        }

        X2 <- newX2(mf[, bl2$get_names(), drop = FALSE])
        K2 <- X2$K
        X2 <- X2$X
        if (!is.null(l2)) K2 <- l2 * K2
        if (options("mboost_useMatrix")$mboost_useMatrix) {
            X2 <- Matrix(X2)
            K2 <- Matrix(K2)
        }

        X <- kronecker(X1, matrix(1, nc = ncol(X2))) *
             kronecker(matrix(1, nc = ncol(X1)), X2)
        K <- kronecker(K1, diag(ncol(X2))) +
             kronecker(diag(ncol(K1)), K2)
        list(X = X, K = K)
    }

    ret$dpp <- bl_lin(ret, Xfun = Xfun, args = args)

    return(ret)
}

