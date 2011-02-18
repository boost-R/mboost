### P-spline base-learner with (monotonicity) constraints
bmono <- function(..., constraint = c("increasing", "decreasing",
                                      "convex", "concave"),
                  by = NULL, index = NULL, knots = 20, boundary.knots = NULL,
                  degree = 3, differences = 2, df = 4,
                  lambda = NULL, lambda2 = 1e6, niter = 10,
                  intercept = TRUE, contrasts.arg = "contr.treatment") {

    if (!is.null(lambda)) df <- NULL

    cll <- match.call()
    cll[[1]] <- as.name("bmono")

    mf <- list(...)

    ## change this to use lists as for knots
    if (!is.list(constraint)) {
        constraint <- match.arg(constraint)
    } else {
        c.args <- eval(formals(sys.function(sys.parent()))[["constraint"]])
        constraint <- lapply(constraint, match.arg, choices = c.args)
    }
    if (length(mf) > 1){
        if(length(mf) != 2)
            stop("Only one or two covariates can be specified in",
                 sQuote("bmono"))
        if (length(constraint) == 1)
            constraint <- list(constraint, constraint)
    }
    ##

    if (length(mf) == 1 && (is.matrix(mf[[1]]) || is.data.frame(mf[[1]]))) {
        mf <- as.data.frame(mf[[1]])
    } else {
        mf <- as.data.frame(mf)
        cl <- as.list(match.call(expand.dots = FALSE))[2][[1]]
        colnames(mf) <- sapply(cl, function(x) deparse(x))
    }
    stopifnot(is.data.frame(mf))
    if (!(all(sapply(mf, is.numeric))) && !all(sapply(mf, is.ordered))) {
        stop("cannot compute ", sQuote("bmono"),
             ": variables must be numeric or ordered")
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
    DOINDEX <- (nrow(mf) > options("mboost_indexmin")[[1]] ||
                is.factor(mf[[1]]))
    if (is.null(index)) {
        if (!CC || DOINDEX) {
            ### try to remove duplicated observations or
            ### observations with missings
            index <- get_index(mf)
            mf <- mf[index[[1]],,drop = FALSE]
            index <- index[[2]]
        }
    }

    ret <- list(model.frame = function()
                    if (is.null(index)) return(mf) else
                        return(mf[index,,drop = FALSE]),
                get_call = function(){
                    cll <- deparse(cll, width.cutoff=500L)
                    if (length(cll) > 1)
                        cll <- paste(cll, collapse="")
                    cll
                },
                get_data = function() mf,
                get_index = function() index,
                get_vary = function() vary,
                get_names = function() colnames(mf),
                set_names = function(value) {
                    if(length(value) != length(colnames(mf)))
                        stop(sQuote("value"), " must have same length as ",
                             sQuote("colnames(mf)"))
                    for (i in 1:length(value)){
                        cll[[i+1]] <<- as.name(value[i])
                    }
                    attr(mf, "names") <<- value
                })
    class(ret) <- "blg"
    if (!is.factor(mf[[1]])){
        args <- hyper_bbs(mf, vary, knots = knots,
                          boundary.knots =  boundary.knots,
                          degree = degree, differences = differences,
                          df = df, lambda = lambda, center = FALSE)
        args$constraint <- constraint
        args$lambda2 <- lambda2
        args$niter <- niter
        ret$dpp <- bl_mono(ret, Xfun = X_bbs,
                           args = args)
    } else {
        args <- hyper_ols(df = df, lambda = lambda,
                          intercept = intercept,
                          contrasts.arg = contrasts.arg)
        args$constraint <- constraint
        args$lambda2 <- lambda2
        args$niter <- niter
        ret$dpp <- bl_mono(ret, Xfun = X_ols,
                           args = args)
    }
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

    if (length(args$constraint) == 1) {
        if (args$constraint %in% c("increasing", "decreasing")){
            diff_order <- 1
        } else { # i.e. args$constraint %in% c("convex", "concave")
            diff_order <- 2
        }
        D <- V <- lambda2 <- vector(mode = "list", length =2)
        D[[1]] <- diff(diag(ncol(X)), differences = diff_order)
        V[[1]] <- matrix(0, ncol = ncol(X) - diff_order, nrow =  ncol(X) -
        diff_order)

        lambda2[[1]] <- args$lambda2
        lambda2[[2]] <- 0
    }
    if (length(args$constraint) == 2) {
        diff_order <- lapply(args$constraint, function(x){
            ifelse( x %in% c("increasing", "decreasing"), 1, 2) } )

        ## ncol1 = length(knots[[1]]) + degree + 1
        ## ncol2 = length(knots[[1]]) + degree + 1
        ## ncol(X) = ncol1 * ncol2
        ncoli <- lapply(args$knots, function(x)
                        length(x$knots) + args$degree + 1)
        stopifnot(ncoli[[1]] * ncoli[[2]] == ncol(X))

        D <- V <- lambda2 <- vector(mode = "list", length =2)
        D[[1]] <- kronecker(diff(diag(ncoli[[1]]),
                                 differences = diff_order[[1]]),
                            diag(ncoli[[2]]))
        D[[2]] <- kronecker(diag(ncoli[[1]]),
                            diff(diag(ncoli[[2]]),
                                 differences = diff_order[[2]]))
        V[[1]] <- matrix(0, ncol = nrow(D[[1]]), nrow =  nrow(D[[1]]))
        V[[2]] <- matrix(0, ncol = nrow(D[[2]]), nrow =  nrow(D[[2]]))
        if (length(args$lambda2) == 1) {
            lambda2[[1]] <- lambda2[[2]] <- args$lambda2
        } else {
            lambda2 <- args$lambda2
        }
    }

    dpp <- function(weights) {
        weights[!Complete.cases(mf)] <- 0
        w <- weights
        if (!is.null(index))
            w <- .Call("R_ysum", as.double(weights), as.integer(index),
                       PACKAGE = "mboost")
        lambdadf <- df2lambda(X, df = args$df, lambda = args$lambda,
                              dmat = K, weights = w)
        lambda <- lambdadf["lambda"]
        XtX <- crossprod(X * w, X)
        XtX <- XtX + lambda * K

        if (is(X, "Matrix")) {
            if (lambda2[[2]] == 0){
                mysolve <- function(y, V) {
                    XtXC <- Cholesky(forceSymmetric(XtX +
                        lambda2[[1]] * crossprod(D[[1]], V[[1]] %*% D[[1]])))
                    solve(XtXC, crossprod(X, y))
                }
            } else {
                ## use chol
                mysolve <- function(y, V) {
                    XtXC <- Cholesky(forceSymmetric(XtX +
                        lambda2[[1]] * crossprod(D[[1]], V[[1]] %*% D[[1]])) +
                        lambda2[[2]] * crossprod(D[[2]], V[[2]] %*% D[[2]]))
                    solve(XtXC, crossprod(X, y))
                }
            }
        } else { ## not Matrix
            if (lambda2[[2]] == 0){#
                mysolve <- function(y, V)
                    .Call("La_dgesv", XtX +
                          lambda2[[1]] * crossprod(D[[1]], V[[1]] %*% D[[1]]),
                          crossprod(X, y), .Machine$double.eps,
                          PACKAGE = "base")
            } else {
                mysolve <- function(y, V)
                    .Call("La_dgesv", XtX +
                          lambda2[[1]] * crossprod(D[[1]], V[[1]] %*% D[[1]]) +
                          lambda2[[2]] * crossprod(D[[2]], V[[2]] %*% D[[2]]),
                          crossprod(X, y), .Machine$double.eps,
                          PACKAGE = "base")
            }
        }
        fit <- function(y) {
            if (!is.null(index)) {
                y <- .Call("R_ysum", as.double(weights * y), as.integer(index),
                           PACKAGE = "mboost")
            } else {
                y <- y * weights
            }

            for (i in 1:args$niter){
                coef <- mysolve(y, V)
                ## compare old and new V
                tmp1 <- do.call(args$constraint[[1]],
                                args=list(D[[1]] %*% coef))
                if (lambda2[[2]] != 0)
                    tmp2 <- do.call(args$constraint[[2]],
                                    args=list(D[[2]] %*% coef))

                if ( all( V[[1]] == tmp1 ) &&
                    ( lambda2[[2]] == 0 || all( V[[2]] == tmp2 ) ) )
                    break    # if both are equal: done!
                V[[1]] <- tmp1
                if (lambda2[[2]] != 0)
                    V[[2]] <- tmp2
                if (i == args$niter)
                    warning("no convergence of coef in bmono\n",
                            "You could try increasing ", sQuote("niter"),
                            " or ", sQuote("lambda2"))
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
        predict <- function(bm, newdata = NULL,
                            aggregate = c("sum", "cumsum", "none")) {
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

increasing <- function(diffs)
    diag(c(as.numeric(diffs)) <= 0)

decreasing <- function(coef)
    diag(c(as.numeric(diffs)) >= 0)

convex <- function(coef)
    diag(c(as.numeric(diffs)) <= 0)

concave <- function(coef)
    diag(c(as.numeric(diffs)) >= 0)
