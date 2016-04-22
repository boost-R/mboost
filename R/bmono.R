### P-spline base-learner with (monotonicity) constraints
bmono <- function(..., constraint = c("increasing", "decreasing",
                                      "convex", "concave", "none",
                                      "positive", "negative"),
                  type = c("quad.prog", "iterative"),
                  by = NULL, index = NULL, knots = 20, boundary.knots = NULL,
                  degree = 3, differences = 2, df = 4,
                  lambda = NULL, lambda2 = 1e6, niter = 10,
                  intercept = TRUE, contrasts.arg = "contr.treatment",
                  boundary.constraints = FALSE,
                  cons.arg = list(lambda = 1e6, n = NULL, diff_order = NULL)){
                                  #type = c("none", "constant", "linear"))) {

    if (!is.null(lambda)) df <- NULL

    cll <- match.call()
    cll[[1]] <- as.name("bmono")

    mf <- list(...)

    ## change this to use lists as for knots
    if (!is.list(constraint)) {
        constraint <- match.arg(constraint)
    } else {
        c.args <- eval(formals(sys.function())[["constraint"]])
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

    type = match.arg(type)

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
    if (!CC)
        warning("base-learner contains missing values;\n",
                "missing values are excluded per base-learner, ",
                "i.e., base-learners may depend on different",
                " numbers of observations.")
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
        args$type <- type
        args$lambda2 <- lambda2
        args$niter <- niter
        args$boundary.constraints <- boundary.constraints
        if(boundary.constraints) {
            if (is.null(cons.arg$n)){
                ## use 10% of the knots on each side per default
                cons.arg$n <- sapply(args$knots,
                                     function(x) rep(length(x$knots), 2)) * 0.1
            } else {
                if(length(cons.arg$n) != 2 && !is.list(cons.arg$n)){
                    cons.arg$n <- rep(cons.arg$n, 2)
                }
                ## <fixme> was passiert bei bivariatem bmono? </fixme>
            }
            ## diff_order for boundary constraints
            if(is.null(cons.arg$diff_order)){
                ## use same difference order as defined by "constraint":
                ## <FIXME> args$constraint may be a list of length 2 for spatial effects
                cons.arg$diff_order <- differences(args$constraint)
            }
            if(is.null(cons.arg$lambda)){
                cons.arg$lambda <- 1e6
            }
        }
        args$cons.arg <- cons.arg
        ret$dpp <- bl_mono(ret, Xfun = X_bbs, args = args)
    } else {
        args <- hyper_ols(df = df, lambda = lambda,
                          intercept = intercept,
                          contrasts.arg = contrasts.arg)
        args$constraint <- constraint
        args$type <- type
        args$lambda2 <- lambda2
        args$niter <- niter
        ## <FIXME> Was machen wir bei kateg. Effekten? Da muesste das doch auch gehen!
        args$boundary.constraints <- boundary.constraints
        args$cons.arg$n <- cons.arg$n
        ret$dpp <- bl_mono(ret, Xfun = X_ols, args = args)
    }
    return(ret)
}

bl_mono <- function(blg, Xfun, args) {
    mf <- blg$get_data()
    index <- blg$get_index()
    vary <- blg$get_vary()

    newX <- function(newdata = NULL, prediction = FALSE) {
        if (!is.null(newdata)) {
            mf <- check_newdata(newdata, blg, mf)
        }
        ## this argument is currently only used in X_bbs --> bsplines
        args$prediction <- prediction
        return(Xfun(mf, vary, args))
    }
    X <- newX()
    K <- X$K
    X <- X$X

    if (length(args$constraint) == 1) {
        D <- V <- lambda2 <- vector(mode = "list", length = 2)
        ## set up difference matrix
        D[[1]] <- differences(args$constraint, ncol(X))

        if (is.factor(mf[[1]]) && args$intercept) {
            D[[1]][1,1] <- 0
        }

        if (!is.factor(mf[[1]]) && args$boundary.constraints) {
            ## set up boundary constraints
            cons.arg <- args$cons.arg
            idx <- rep(0, ncol(X) - cons.arg$diff_order)
            if (cons.arg$n[1] == 0) {
                lower <- 0
            } else {
                lower <- 1:cons.arg$n[1]
            }
            if (cons.arg$n[2] == length(idx)) {
                upper <- 0
            } else {
                upper <- length(idx) - 1:cons.arg$n[2] + 1
            }
            idx[c(lower, upper)] <- 1
            V3 <- diag(idx)
            D3 <- V3 %*% diff(diag(ncol(X)), differences = cons.arg$diff_order)
        }

        V[[1]] <- matrix(0, ncol = nrow(D[[1]]), nrow =  nrow(D[[1]]))
        lambda2[[1]] <- ifelse(args$constraint == "none", 0, args$lambda2)
        lambda2[[2]] <- 0
        if (args$boundary.constraints) {
            lambda3 <- cons.arg$lambda
        } else {
            lambda3 <- 0
        }
    }
    if (length(args$constraint) == 2) {
        if (is.factor(mf[[1]]))
            stop(paste("Bivariate monotonic effects currently not",
                       "implemented for ordered factors"))
        ## ncol1 = length(knots[[1]]) + degree + 1
        ## ncol2 = length(knots[[2]]) + degree + 1
        ## ncol(X) = ncol1 * ncol2
        ncoli <- lapply(args$knots, function(x)
                        length(x$knots) + args$degree + 1)
        stopifnot(ncoli[[1]] * ncoli[[2]] == ncol(X))

        D <- V <- lambda2 <- vector(mode = "list", length =2)
        ## set up difference matrices
        D <- differences(args$constraint, ncoli)
        idx <- !sapply(D, is.null)
        V[idx] <- lapply(D[idx], function(m) matrix(0, nrow(m), nrow(m)))

        if (length(args$lambda2) == 1) {
            lambda2[[1]] <- lambda2[[2]] <- args$lambda2
        } else {
            lambda2 <- args$lambda2
        }
        ## set lambda2 = 0 if no constraint is used
        if (any(none <- args$constraint == "none"))
            lambda2[none] <- 0
        if (any(none <- lambda2 == 0))
            args$constraint[none] <- "none"
        ## <FIXME> Boundary constraints for bivariate smooths are currently not
        ## implemented
        if (args$boundary.constraints)
            warning("Boundary constraints for bivariate smooths",
                    "are currently not implemented")
        lambda3 <- 0
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

        if (args$type == "iterative") {
            fct <- define_solver(lambda2, lambda3, X)
            ## deparsing and parsing again needed to tidy-up code.
            mysolve <- eval(parse(text = deparse(eval(parse(text = fct)))))
        }

        fit <- function(y) {
            if (!is.null(index)) {
                y <- .Call("R_ysum", as.double(weights * y), as.integer(index),
                           PACKAGE = "mboost")
            } else {
                y <- y * weights
            }

            if (args$type == "iterative") {
                for (i in 1:args$niter){
                    coef <- mysolve(y, V)
                    if (args$constraint[[1]] == "none")
                        break ## as there is no need to iterate

                    ## compare old and new V
                    tmp1 <- violations(D[[1]] %*% coef)
                    if (lambda2[[2]] != 0)
                        tmp2 <- violations(D[[2]] %*% coef)

                    if ( all( V[[1]] == tmp1 ) &&
                        ( lambda2[[2]] == 0 || all( V[[2]] == tmp2 ) ) )
                        break    # if both are equal: done!
                        #if (args$boundary.constraints &&
                        #    all( V[[1]][-idxB, -idxB] == tmp1[-idxB, -idxB]) )
                        #    break   # if both are equal (without V for boundary
                        #            # constraints): done!
                    V[[1]] <- tmp1
                        #if (args$boundary.constraints) {
                        #    V[[1]][idxFlat, idxFlat] <- diag(rep(1, length(idxFlat)))
                        #}
                    if (lambda2[[2]] != 0)
                        V[[2]] <- tmp2
                    if (i == args$niter)
                        warning("no convergence of coef in bmono\n",
                                "You could try increasing ", sQuote("niter"),
                                " or ", sQuote("lambda2"))
                }
            } else {  ## i.e. type == "quad.prog"
                if (lambda2[[2]] == 0) {
                    coef <- solveLSEI(XtX, crossprod(X, y), D = D[[1]])
                } else {
                    coef <- solveLSEI(XtX, crossprod(X, y), D = D)
                }
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
            warning("hatvalues might be a very poor approximation",
                    "for monotonic base-learners.")
            if(lambda2[[2]] == 0)
                pen2 <- lambda2[[1]] * crossprod(D[[1]], V[[1]] %*% D[[1]])
            else
                pen2 <- lambda2[[1]] * crossprod(D[[1]], V[[1]] %*% D[[1]]) +
                        lambda2[[2]] * crossprod(D[[2]], V[[2]] %*% D[[2]])
            ret <- as.matrix(tcrossprod(X %*% solve(XtX + pen2), X * w))
            if (is.null(index)) return(ret)
            return(ret[index,index])
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
                ## Use sparse data represenation if data set is huge
                ## and a data.frame
                if (is.data.frame(newdata) && nrow(newdata) > options("mboost_indexmin")[[1]]) {
                    index <- get_index(newdata)
                    newdata <- newdata[index[[1]], , drop = FALSE]
                    index <- index[[2]]
                }
                X <- newX(newdata, prediction = TRUE)$X
            }
            aggregate <- match.arg(aggregate)
            pr <- switch(aggregate, "sum" =
                as(X %*% rowSums(cf), "matrix"),
            "cumsum" = {
                as(X %*% .Call("R_mcumsum", as(cf, "matrix"),
                               PACKAGE = "mboost"), "matrix")
            },
            "none" = as(X %*% cf, "matrix"))
            if (is.null(index))
                return(pr[, , drop = FALSE])
            return(pr[index, , drop = FALSE])
        }

        ret <- list(fit = fit, hatvalues = hatvalues,
                    predict = predict, df = df,
                    Xnames = colnames(X))
        class(ret) <- c("bl_lin", "bl")
        return(ret)

    }
    return(dpp)
}

violations <- function(diffs)
    diag(c(as.numeric(diffs)) <= 0)

define_solver <- function(lambda2, lambda3, X) {
    ## define function as text and eval(parse()) later.
    l2txt <- "+ lambda2[[2]] * crossprod(D[[2]], V[[2]] %*% D[[2]])"
    l3txt <- "+ lambda3 * crossprod(D3, V3 %*% D3)"
    fct <- c("function(y, V) {",
             "    XtXC <- Cholesky(forceSymmetric(XtX +",
             "       lambda2[[1]] * crossprod(D[[1]], V[[1]] %*% D[[1]])",
             ## add if lambda2[[2]] != 0
             ifelse(lambda2[[2]] != 0, l2txt,""),
             ## add if lambda3 != 0
             ifelse(lambda3 != 0, l3txt,""),
             "                                   ))",
             "    solve(XtXC, crossprod(X, y), LINPACK = FALSE)",
             "}"
             )

    if (!is(X, "Matrix")) {
        ## some lines must be replaced in order to solve directly
        fct[2] <- '    solve(XtX +'
        fct[6] <- "          , crossprod(X, y),"
        fct[7] <- '          LINPACK = FALSE)'
    }
    fct
}
