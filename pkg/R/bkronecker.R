
### workhorse for fitting matrix-response (Ridge-penalized) baselearners
### Y = kronecker(X2, X1)
### see Currie, Durban, Eilers (2006, JRSS B)
bl_lin_matrix <- function(blg, Xfun, args) {

    mf <- blg$get_data()
    index <- blg$get_index()
    vary <- blg$get_vary()

    newX <- function(newdata = NULL) {
        if (!is.null(newdata)) {
            stopifnot(all(names(newdata) == names(blg)))
            # stopifnot(all(class(newdata) == class(mf)))
            mf <- newdata[names(blg)]
        }
        return(Xfun(mf, vary, args))
    }
    X <- newX()
    K <- X$K
    X <- X$X
    c1 <- ncol(X$X1)
    c2 <- ncol(X$X2)
    n1 <- nrow(X$X1)
    n2 <- nrow(X$X2)

    G <- function(x) {
        one <- matrix(rep(1, ncol(x)), nrow = 1)
        kronecker(x, one) * kronecker(one, x)
    }

    dpp <- function(weights) {

        W <- matrix(weights, nrow = n1, ncol = n2)

        ### X = kronecker(X2, X1)
        XtX <- as(crossprod(G(X$X1), W) %*% G(X$X2), "matrix")
        XtX <- array(XtX, c(c1, c1, c2, c2))
        XtX <- matrix(aperm(XtX, c(1, 3, 2, 4)), nrow = c1 * c2)

        ### <FIXME> This does not happen in bl_lin / df2lambda.
        ### For one base learner only, it makes sense to allow
        ### for a direct choice of lambda (regardless of df)
        ### </FIXME>
        if (is.null(args$lambda)) {

            ### <FIXME>: is there a better way to feed XtX into lambdadf?
            lambdadf <- df2lambda(matrix(0, ncol = ncol(X$X1) + ncol(X$X2)), 
                                  df = args$df, lambda = args$lambda,
                                  dmat = K, weights = weights, XtX = XtX)
            ### </FIXME>
            lambda <- lambdadf["lambda"]
        } else {
            lambda <- args$lambda
        }
        XtX <- XtX + lambda * K

        ## matrizes of class dgeMatrix are dense generic matrices; they should
        ## be coerced to class matrix and handled in the standard way
        if (is(XtX, "Matrix") && !extends(class(XtX), "dgeMatrix")) {
            XtXC <- Cholesky(forceSymmetric(XtX))
            mysolve <- function(y) {
                Y <- matrix(y, nrow = n1) * W
                XWY <- as.vector(crossprod(X$X1, Y) %*% X$X2)
                solve(XtXC, XWY)  ## special solve routine from
                                  ## package Matrix
                }
        } else {
            if (is(XtX, "Matrix")) {
                ## coerce Matrix to matrix
                XtX <- as(XtX, "matrix")
            }
            mysolve <- function(y) {
                Y <- matrix(y, nrow = n1) * W
                XWY <- crossprod(X$X1, Y) %*% X$X2
                .Call("La_dgesv", XtX, matrix(as(XWY, "matrix"), ncol = 1), 
                      .Machine$double.eps,
                      PACKAGE = "base")
            }
        }

        cfprod <- function(b) tcrossprod(X$X1 %*% b, X$X2)

        fit <- function(y) {
            coef <- as(mysolve(y), "matrix")
            if (nrow(coef) != c1) coef <- matrix(as.vector(coef), nrow = c1)
            f <- cfprod(coef)
            f <- as(f, "matrix")
            if (options("mboost_Xmonotone")$mboost_Xmonotone) {
                md <- apply(f, 1, function(x) min(diff(x)))
                if (any(md < 0)) {
                    coef <- matrix(0, nrow = nrow(coef), ncol = ncol(coef))
                    f <- matrix(0, nrow = nrow(f), ncol = ncol(f))
                }
            }
            ret <- list(model = coef,
                        fitted = function() as.vector(f))
            class(ret) <- c("bm_lin", "bm")
            ret
        }

        ### <FIXME> check for large n, option?
        hatvalues <- function() {
            return(NULL)
        }
        ### </FIXME>

        ### actually used degrees of freedom (trace of hat matrix)
        df <- function() lambdadf

        ### prepare for computing predictions
        predict <- function(bm, newdata = NULL, aggregate = c("sum", "cumsum", "none")) {
            cf <- lapply(bm, function(x) x$model)
            if(!is.null(newdata)) {
                index <- NULL
                nm <- names(blg)
                newdata <- newdata[nm]
                X <- newX(newdata)$X
            }
            ncfprod <- function(b) 
                as.vector(as(tcrossprod(X$X1 %*% b, X$X2), "matrix"))
            aggregate <- match.arg(aggregate)
            pr <- switch(aggregate, "sum" = {
                cf2 <- 0
                for (b in cf) cf2 <- cf2 + b
                ncfprod(cf2)
            },
            "cumsum" = {
                cf2 <- 0
                ret <- c()
                for (b in cf) {
                    cf2 <- cf2 + b
                    ret <- cbind(ret, ncfprod(cf2))
                }
                ret
            },
            "none" = {
                ret <- c()
                for (b in cf) {   
                    ret <- cbind(ret, ncfprod(b))
                }
                ret
            })
            return(pr)
        }

        ret <- list(fit = fit, hatvalues = hatvalues,
                    predict = predict, df = df,
                    Xnames = c(interaction(colnames(X$X1), colnames(X$X2))))
        class(ret) <- c("bl_lin", "bl")
        return(ret)

    }
    return(dpp)
}

### kronecker product of two baselearners
### bbs(x1) %O% bbs(x2) means that
### X = kronecker(X2, X1) and x1 varies fastest
"%O%" <- function(bl1, bl2) {

    if (is.list(bl1) && !inherits(bl1, "blg"))
        return(lapply(bl1, "%X%", bl2 = bl2))

    if (is.list(bl2) && !inherits(bl2, "blg"))
        return(lapply(bl2, "%X%", bl1 = bl1))

    cll <- paste(bl1$get_call(), "%O%",
                 bl2$get_call(), collapse = "")
    stopifnot(inherits(bl1, "blg"))
    stopifnot(inherits(bl2, "blg"))

    mf1 <- model.frame(bl1)
    mf2 <- model.frame(bl2)
    stopifnot(!any(colnames(mf1) %in%
                   colnames(mf2)))
    mf <- c(mf1, mf2)
    stopifnot(all(complete.cases(mf[[1]])))
    stopifnot(all(complete.cases(mf[[2]])))

    index <- NULL

    vary <- ""

    ret <- list(model.frame = function()
                    return(mf),
                get_call = function(){
                    cll <- deparse(cll, width.cutoff=500L)
                    if (length(cll) > 1)
                        cll <- paste(cll, collapse="")
                    cll
                },
                get_data = function() mf,
                get_index = function() index,
                get_vary = function() vary,
                get_names = function() names(mf),
                ## <FIXME> Is this all we want to change if we set names here?
                set_names = function(value) attr(mf, "names") <<- value)
                ## </FIXME>
    class(ret) <- "blg"

    args1 <- environment(bl1$dpp)$args
    args2 <- environment(bl2$dpp)$args
    l1 <- args1$lambda
    l2 <- args2$lambda
    if (!is.null(l1) && !is.null(l2)) {
        args <- list(lambda = l1 + l2, df = NULL)
    } else {
        args <- list(lambda = NULL,
            df = ifelse(is.null(args1$df), 1, args1$df) *
                 ifelse(is.null(args2$df), 1, args2$df))
    }

    Xfun <- function(mf, vary, args) {

        newX1 <- environment(bl1$dpp)$newX
        newX2 <- environment(bl2$dpp)$newX

        X1 <- newX1(as.data.frame(mf[bl1$get_names()]))
        K1 <- X1$K
        X1 <- X1$X
        if (!is.null(l1)) K1 <- l1 * K1
        MATRIX <- options("mboost_useMatrix")$mboost_useMatrix
        if (MATRIX & !is(X1, "Matrix"))
            X1 <- Matrix(X1)
        if (MATRIX & !is(K1, "Matrix"))
            K1 <- Matrix(K1)

        X2 <- newX2(as.data.frame(mf[bl2$get_names()]))
        K2 <- X2$K
        X2 <- X2$X
        if (!is.null(l2)) K2 <- l2 * K2
        if (MATRIX & !is(X2, "Matrix"))
            X2 <- Matrix(X2)
        if (MATRIX & !is(K2, "Matrix"))
            K2 <- Matrix(K2)

        K <- kronecker(K2, diag(ncol(X1))) +
             kronecker(diag(ncol(X2)), K1)
        list(X = list(X1 = X1, X2 = X2), K = K)
    }

    ret$dpp <- bl_lin_matrix(ret, Xfun = Xfun, args = args)

    return(ret)
}
