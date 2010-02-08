bmono <- function(..., constraint = c("increasing", "decreasing", "convex", "concave"), by = NULL, index = NULL, knots = 20, degree = 3,
                 differences = 2, df = 4, lambda = NULL, lambda2 = 1e6, niter=10) {

    if (!is.null(lambda)) df <- NULL

    cll <- match.call()
    cll[[1]] <- as.name("bmono")
    cll <- deparse(cll)

    mf <- list(...)
    constraint <- match.arg(constraint)
    if (length(mf) == 1 && (is.matrix(mf[[1]]) || is.data.frame(mf[[1]]))) {
        mf <- as.data.frame(mf[[1]])
    } else {
        mf <- as.data.frame(mf)
        cl <- as.list(match.call(expand.dots = FALSE))[2][[1]]
        colnames(mf) <- sapply(cl, function(x) deparse(x))
    }
    stopifnot(is.data.frame(mf))
    if(!(all(sapply(mf, is.numeric)))) {
        if (ncol(mf) == 1) return(bols(as.data.frame(...), by = by, index = index))
        stop("cannot compute ", sQuote("bmono"), " for non-numeric variables")
    }
    ### use bols when appropriate
    if (!is.null(df)) {
        if (df <= (ncol(mf) + 1))
            return(bols(as.data.frame(...), by = by, index = index))
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
                    cll <<- paste("bmono", "(", paste(colnames(mf),
                        collapse = ", "), ")", sep = "")
                })
    class(ret) <- "blg"
    ret$dpp <- bl_mono(ret, Xfun = X_bbs,
                      args = c(hyper_bbs(mf, vary, knots = knots,
                               degree = degree, differences = differences,
                               df = df, lambda = lambda, center = FALSE),
                               constraint = constraint, lambda2 = lambda2, niter=niter))
    return(ret)
}

bl_mono <- function(blg, Xfun, args) {
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

    if (args$constraint %in% c("increasing", "decreasing")){
        diff_order <- 1
    } else { # i.e. args$constraint %in% c("convex", "concave")
        diff_order <- 2
    }

    D <- diff(diag(ncol(X)), differences = diff_order)
    V <- matrix(0, ncol = ncol(X) - diff_order, nrow =  ncol(X) - diff_order)
    lambda2 <- args$lambda2

    dpp <- function(weights) {
        weights[!Complete.cases(mf)] <- 0
        w <- weights
        if (!is.null(index))
            w <- .Call("R_ysum", as.double(weights), as.integer(index), PACKAGE = "mboost")
        lambdadf <- df2lambda(X, df = args$df, lambda = args$lambda,
                              dmat = K, weights = w)
        lambda <- lambdadf["lambda"]
        XtX <- crossprod(X * w, X)
        XtX <- XtX + lambda * K

        if (is(X, "Matrix")) {
            ## use chol
            mysolve <- function(y, V) {
                XtXC <- Cholesky(forceSymmetric(XtX + lambda2 * crossprod(D, V %*% D)))
                solve(XtXC, crossprod(X, y))
            }
        } else {
            mysolve <- function(y, V)
                .Call("La_dgesv", XtX + lambda2 * crossprod(D, V %*% D), crossprod(X, y), .Machine$double.eps,
                      PACKAGE = "base")
        }
        fit <- function(y) {
            if (!is.null(index)) {
                y <- .Call("R_ysum", as.double(weights * y), as.integer(index), PACKAGE = "mboost")
            } else {
                y <- y * weights
            }

            for (i in 1:args$niter){
                coef <- mysolve(y, V)
                if (all( V == (V <- do.call(args$constraint, args=list(as.vector(coef)))) )) # compare old and new V
                    break                                                                    # if both are equal: done!
                if (i == args$niter) warning("no convergence of coef in bmono")
            }

            ret <- list(model = coef,
                        fitted = function() {
                            ret <- as.vector(X %*% coef)
                            if (is.null(index)) return(ret)
                            return(ret[index])
                        })
            class(ret) <- c("bm_lin", "bm")
            ret
        }


        hatvalues <- function() {
            stop("not possible for monotonic base-learners")
        }

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

increasing <- function(coef){
    diag(as.numeric(c(diff(coef, differences=1))) <= 0)
}
decreasing <- function(coef){
    diag(as.numeric(c(diff(coef, differences=1))) >= 0)
}

convex <- function(coef){
    diag(as.numeric(c(diff(coef, differences=2))) <= 0)
}

concave <- function(coef){
    diag(as.numeric(c(diff(coef, differences=2))) >= 0)
}
