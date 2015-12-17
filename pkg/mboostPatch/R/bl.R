
### compute Ridge shrinkage parameter lambda from df
### or the other way round
df2lambda <- function(X, df = 4, lambda = NULL, dmat = NULL, weights,
                      XtX = NULL) {

    stopifnot(xor(is.null(df), is.null(lambda)))
    if (!is.null(df)) {
        rank_X <- rankMatrix(X, method = 'qr', warn.t = FALSE)
        if (df >= rank_X) {
            if (df > rank_X)
                warning(sQuote("df"),
                        " too large:\n  Degrees of freedom cannot be larger",
                        " than the rank of the design matrix.\n",
                        "  Unpenalized base-learner with df = ",
                        rank_X, " used. Re-consider model specification.")
            return(c(df = df, lambda = 0))
        }
    }
    if (!is.null(lambda))
        if (lambda == 0)
            return(c(df = rankMatrix(X), lambda = 0))

    ## check for possible instability
    if (options("mboost_check_df2lambda")[[1]] && max(abs(X)) > 10)
        warning("Some absolute values in design matrix are greater 10. Hence, ",
                sQuote("df2lambda"), " might be numerically instable.\n  ",
                "See documentation of argument ", sQuote("by"),
                " in ?bbs for further information.",
                immediate. = TRUE)
    ## instable df2lambda might for example occure if one uses bbs(x, by = z)
    ## with large values of z

    # Demmler-Reinsch Orthogonalization (cf. Ruppert et al., 2003,
    # Semiparametric Regression, Appendix B.1.1).

    ### there may be more efficient ways to compute XtX, but we do this
    ### elsewhere (e.g. in %O%)
    if (is.null(XtX)) XtX <- crossprod(X * weights, X)
    if (is.null(dmat)) {
        if(is(XtX, "Matrix")) diag <- Diagonal
        dmat <- diag(ncol(XtX))
    }
    A <- XtX + dmat * options("mboost_eps")[[1]]
    ## make sure that A is also numerically symmetric
    if (is(A, "Matrix"))
        A <- forceSymmetric(A)
    Rm <- solve(chol(A))
    ## singular value decomposition without singular vectors
    d <- try(svd(crossprod(Rm, dmat) %*% Rm, nu=0, nv=0)$d)
    ## if unsucessfull try the same computation but compute singular vectors as well
    if (inherits(d, "try-error"))
        d <- svd(crossprod(Rm, dmat) %*% Rm)$d

    ### option
    if (options("mboost_dftraceS")[[1]]){
        ## df := trace(S)
        dfFun <- function(lambda) sum(1 / (1 + lambda * d))
    } else {
        ## df := trace(2S - S'S)
        dfFun <- function(lambda) 2 * sum( 1/(1+lambda*d) ) - sum( 1/(1+lambda*d)^2 )
    }

    if (!is.null(lambda))
        return(c(df = dfFun(lambda), lambda = lambda))
    if (df >= length(d)) return(c(df = df, lambda = 0))

    # search for appropriate lambda using uniroot
    df2l <- function(lambda)
        dfFun(lambda) - df

    ### option
    lambdaMax <- options("mboost_lambdaMax")[[1]]

    if (df2l(lambdaMax) > 0){
        if (df2l(lambdaMax) > sqrt(.Machine$double.eps))
            warning("lambda needs to be larger than ", lambdaMax, " for given ",
                    sQuote("df"), ";\n  setting lambda = ", lambdaMax,
                    " leeds to an deviation from ", sQuote("df"), " of ",
                    df2l(lambdaMax), ";\n  You can increase lambda_max via ",
                    sQuote("options(mboost_lambdaMax = value)"))
        return(c(df = df, lambda = lambdaMax))
    }
    lambda <- uniroot(df2l, c(0, lambdaMax), tol = sqrt(.Machine$double.eps))$root
    if (abs(df2l(lambda)) > sqrt(.Machine$double.eps))
        warning("estimated degrees of freedom differ from ", sQuote("df"),
                " by ", df2l(lambda))
    return(c(df = df, lambda = lambda))
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
        fac <- sapply(mf[colnames(mf) != vary], is.factor)
        if (any(fac)){
            if (!is.list(args$contrasts.arg)){
                ## first part needed to prevent warnings from calls such as
                ## contrasts.arg = contr.treatment(4, base = 1):
                if (is.character(args$contrasts.arg) &&
                    args$contrasts.arg == "contr.dummy"){
                    if (!args$intercept)
                        stop('"contr.dummy" can only be used with ',
                             sQuote("intercept = TRUE"))
                    fm <- paste(fm, "-1")
                } else {
                    txt <- paste("list(", paste(colnames(mf)[colnames(mf) != vary][fac],
                                                "= args$contrasts.arg", collapse = ", "),")")
                    args$contrasts.arg <- eval(parse(text=txt))
                }
            } else {
                ## if contrasts are given as list check if "contr.dummy" is specified
                if (any(args$contrasts.arg == "contr.dummy"))
                    stop('"contr.dummy"',
                         " can only be used for all factors at the same time.\n",
                         "Use ", sQuote('contrasts.arg = "contr.dummy"'),
                         " to achieve this.")
            }
        } else {
            args$contrasts.arg <- NULL
        }
        X <- model.matrix(as.formula(fm), data = mf, contrasts.arg = args$contrasts.arg)
        if (!is.null(args$contrasts.arg) && args$contrasts.arg == "contr.dummy")
            attr(X, "contrasts") <- lapply(attr(X, "contrasts"),
                                           function(x) x <- "contr.dummy")
        contr <- attr(X, "contrasts")
        if (!args$intercept)
            X <- X[ , -1, drop = FALSE]
        MATRIX <- any(dim(X) > c(500, 50)) && any(fac)
        MATRIX <- MATRIX && options("mboost_useMatrix")$mboost_useMatrix
        if (MATRIX) {
            diag <- Diagonal
            cbind <- cBind
            if (!is(X, "Matrix"))
                X <- Matrix(X)
        }
        if (vary != "") {
            by <- model.matrix(as.formula(paste("~", vary, collapse = "")),
                               data = mf)[ , -1, drop = FALSE] # drop intercept
            DM <- lapply(1:ncol(by), function(i) {
                ret <- X * by[, i]
                colnames(ret) <- paste(colnames(ret), colnames(by)[i], sep = ":")
                ret
            })
            if (is(X, "Matrix")) {
                X <- do.call("cBind", DM)
            } else {
                X <- do.call("cbind", DM)
            }
        }
    }
    ### <FIXME> penalize intercepts???
    ### set up penalty matrix
    ANOVA <- (!is.null(contr) && (length(contr) == 1)) && (ncol(mf) == 1)
    K <- diag(ncol(X))
    ### for ordered factors use difference penalty
    if (ANOVA && any(sapply(mf[, names(contr), drop = FALSE], is.ordered))) {
        K <- diff(diag(ncol(X) + 1), differences = 1)[, -1, drop = FALSE]
        if (vary != "" && ncol(by) > 1){       # build block diagonal penalty
            suppressMessages(K <- kronecker(diag(ncol(by)), K))
        }
        K <- crossprod(K)
    }
    ### </FIXME>
    if (is(X, "Matrix") && !is(K, "Matrix"))
        K <- Matrix(K)
    list(X = X, K = K)
}

### hyper parameters for P-splines baselearner (including tensor product P-splines)
hyper_bbs <- function(mf, vary, knots = 20, boundary.knots = NULL, degree = 3,
                      differences = 2, df = 4, lambda = NULL, center = FALSE,
                      cyclic = FALSE, constraint = "none", deriv = 0L) {

    knotf <- function(x, knots, boundary.knots) {
        if (is.null(boundary.knots))
            boundary.knots <- range(x, na.rm = TRUE)
        ## <fixme> At the moment only NULL or 2 boundary knots can be specified.
        ## Knot expansion is done automatically on an equidistand grid.</fixme>
        if ((length(boundary.knots) != 2) || !boundary.knots[1] < boundary.knots[2])
            stop("boundary.knots must be a vector (or a list of vectors) ",
                 "of length 2 in increasing order")
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
    if (is.list(boundary.knots)) if(!all(names(boundary.knots) %in% nm))
        stop("variable names and boundary.knot names must be the same")
    if (!identical(center, FALSE) && cyclic)
        stop("centering of cyclic covariates not yet implemented")
    ret <- vector(mode = "list", length = length(nm))
    names(ret) <- nm
    for (n in nm)
        ret[[n]] <- knotf(mf[[n]], if (is.list(knots)) knots[[n]] else knots,
                          if (is.list(boundary.knots)) boundary.knots[[n]]
                          else boundary.knots)
    if (cyclic & constraint != "none")
        stop("constraints not implemented for cyclic B-splines")
    stopifnot(is.numeric(deriv) & length(deriv) == 1)
    ## prediction is usually set in/by newX()
    list(knots = ret, degree = degree, differences = differences,
         df = df, lambda = lambda, center = center, cyclic = cyclic,
         Ts_constraint = constraint, deriv = deriv, prediction = FALSE)
}

### model.matrix for P-splines baselearner (including tensor product P-splines)
X_bbs <- function(mf, vary, args) {

    stopifnot(is.data.frame(mf))
    mm <- lapply(which(colnames(mf) != vary), function(i) {
        if (!args$cyclic) {
            X <- bsplines(mf[[i]],
                          knots = args$knots[[i]]$knots,
                          boundary.knots = args$knots[[i]]$boundary.knots,
                          degree = args$degree,
                          Ts_constraint = args$Ts_constraint,
                          deriv = args$deriv, extrapolation = args$prediction)
        } else { ## if cyclic spline
            X <- cbs(mf[[i]],
                     knots = args$knots[[i]]$knots,
                     boundary.knots = args$knots[[i]]$boundary.knots,
                     degree = args$degree,
                     deriv = args$deriv)
        }
        class(X) <- "matrix"
        return(X)
    })
    ### options
    MATRIX <- any(sapply(mm, dim) > c(500, 50)) || (length(mm) > 1)
    MATRIX <- MATRIX && options("mboost_useMatrix")$mboost_useMatrix
    if (MATRIX) {
        diag <- Diagonal
        cbind <- cBind
        for (i in 1:length(mm)){
            tmp <- attributes(mm[[i]])[c("degree", "knots", "Boundary.knots")]
            mm[[i]] <- Matrix(mm[[i]])
            attributes(mm[[i]])[c("degree", "knots", "Boundary.knots")] <- tmp
        }
    }
    if (length(mm) == 1) {
        X <- mm[[1]]
        if (vary != "") {
            by <- model.matrix(as.formula(paste("~", vary, collapse = "")),
                               data = mf)[ , -1, drop = FALSE] # drop intercept
            DM <- lapply(1:ncol(by), function(i) {
                ret <- X * by[, i]
                colnames(ret) <- paste(colnames(ret), colnames(by)[i], sep = ":")
                ret
            })
            if (is(X, "Matrix")) {
                X <- do.call("cBind", DM)
            } else {
                X <- do.call("cbind", DM)
            }
        }
        if (args$differences > 0){
            if (!args$cyclic) {
                K <- diff(diag(ncol(mm[[1]])), differences = args$differences)
            } else {
                ## cyclic P-splines
                differences <- args$differences
                K <- diff(diag(ncol(mm[[1]]) + differences),
                          differences = differences)
                tmp <- K[,(1:differences)]   # save first "differences" columns
                K <- K[,-(1:differences)]    # drop first "differences" columns
                indx <- (ncol(mm[[1]]) - differences + 1):(ncol(mm[[1]]))
                K[,indx] <- K[,indx] + tmp   # add first "differences" columns
            }
        } else {
            if (args$differences != 0)
                stop(sQuote("differences"), " must be an non-neative integer")
            K <- diag(ncol(mm[[1]]))
        }

        if (vary != "" && ncol(by) > 1){       # build block diagonal penalty
                suppressMessages(K <- kronecker(diag(ncol(by)), K))
        }
        if (!identical(args$center, FALSE)) {
            tmp <- attributes(X)[c("degree", "knots", "Boundary.knots")]
            center <- match.arg(as.character(args$center),
                                choices = c("TRUE", "differenceMatrix", "spectralDecomp"))
            if (center == "TRUE") center <- "differenceMatrix"
            X <- switch(center,
                ### L = t(D) in Section 2.3. of Fahrmeir et al. (2004, Stat Sinica)
                "differenceMatrix" = tcrossprod(X, K) %*% solve(tcrossprod(K)),
                ### L = \Gamma \Omega^1/2 in Section 2.3. of
                ### Fahrmeir et al. (2004, Stat Sinica)
                "spectralDecomp" = {
                    SVD <- eigen(crossprod(K), symmetric = TRUE)
                    ev <- SVD$vector[, 1:(ncol(X) - args$differences), drop = FALSE]
                    ew <- SVD$values[1:(ncol(X) - args$differences), drop = FALSE]
                    X %*% ev %*% diag(1/sqrt(ew))
                }
            )
            attributes(X)[c("degree", "knots", "Boundary.knots")] <- tmp
            K <- diag(ncol(X))
        } else {
            K <- crossprod(K)
        }
        if (!is.null(attr(X, "Ts_constraint"))) {
            D <- attr(X, "D")
            K <- crossprod(D, K) %*% D
        }
    }
    if (length(mm) == 2) {
        suppressMessages(
            X <- kronecker(mm[[1]], matrix(1, ncol = ncol(mm[[2]]))) *
                 kronecker(matrix(1, ncol = ncol(mm[[1]])), mm[[2]])
            )
        if (vary != "") {
            by <- model.matrix(as.formula(paste("~", vary, collapse = "")),
                               data = mf)[ , -1, drop = FALSE] # drop intercept
            DM <- X * by[,1]
            if (ncol(by) > 1){
                for (i in 2:ncol(by))
                    DM <- cbind(DM, (X * by[,i]))
            }
            X <- DM
            ### <FIXME> Names of X if by is given
        }
        if (args$differences > 0){
            if (!args$cyclic) {
                Kx <- diff(diag(ncol(mm[[1]])), differences = args$differences)
                Ky <- diff(diag(ncol(mm[[2]])), differences = args$differences)
            } else {
                ## cyclic P-splines
                differences <- args$differences
                Kx <- diff(diag(ncol(mm[[1]]) + differences),
                           differences = differences)
                Ky <- diff(diag(ncol(mm[[2]]) + differences),
                           differences = differences)

                tmp <- Kx[,(1:differences)]   # save first "differences" columns
                Kx <- Kx[,-(1:differences)]    # drop first "differences" columns
                indx <- (ncol(mm[[1]]) - differences + 1):(ncol(mm[[1]]))
                Kx[,indx] <- Kx[,indx] + tmp   # add first "differences" columns

                tmp <- Ky[,(1:differences)]   # save first "differences" columns
                Ky <- Ky[,-(1:differences)]    # drop first "differences" columns
                indx <- (ncol(mm[[2]]) - differences + 1):(ncol(mm[[2]]))
                Ky[,indx] <- Ky[,indx] + tmp   # add first "differences" columns
            }
        } else {
            if (args$differences != 0)
                stop(sQuote("differences"), " must be an non-neative integer")
            Kx <- diag(ncol(mm[[1]]))
            Ky <- diag(ncol(mm[[2]]))
        }

        Kx <- crossprod(Kx)
        Ky <- crossprod(Ky)
        suppressMessages(
            K <- kronecker(Kx, diag(ncol(mm[[2]]))) +
                kronecker(diag(ncol(mm[[1]])), Ky)
            )
        if (vary != "" && ncol(by) > 1){       # build block diagonal penalty
            suppressMessages(K <- kronecker(diag(ncol(by)), K))
        }
        if (!identical(args$center, FALSE)) {
            ### L = \Gamma \Omega^1/2 in Section 2.3. of Fahrmeir et al.
            ### (2004, Stat Sinica), always
            L <- eigen(K, symmetric = TRUE)
            L$vectors <- L$vectors[,1:(ncol(X) - args$differences^2), drop = FALSE]
            L$values <- sqrt(L$values[1:(ncol(X) - args$differences^2), drop = FALSE])
            L <- L$vectors %*% (diag(length(L$values)) * (1/L$values))
            X <- as(X %*% L, "matrix")
            K <- as(diag(ncol(X)), "matrix")
        }
    }
    if (length(mm) > 2)
        stop("not possible to specify more than two variables in ",
             sQuote("..."), " argument of smooth base-learners")

    ## compare specified degrees of freedom to dimension of null space
    if (!is.null(args$df)){
        rns <- ncol(K) - qr(as.matrix(K))$rank # compute rank of null space
        if (rns == args$df)
            warning( sQuote("df"), " equal to rank of null space ",
                    "(unpenalized part of P-spline);\n  ",
                    "Consider larger value for ", sQuote("df"),
                    " or set ", sQuote("center != FALSE"), ".", immediate.=TRUE)
        if (rns > args$df)
            stop("not possible to specify ", sQuote("df"),
                 " smaller than the rank of the null space\n  ",
                 "(unpenalized part of P-spline). Use larger value for ",
                 sQuote("df"), " or set ", sQuote("center != FALSE"), ".")
    }
    return(list(X = X, K = K))
}

### Linear baselearner, potentially Ridge-penalized (but not by default)
bols <- function(..., by = NULL, index = NULL, intercept = TRUE, df = NULL,
                 lambda = 0, contrasts.arg = "contr.treatment") {

    if (!is.null(df)) lambda <- NULL

    cll <- match.call()
    cll[[1]] <- as.name("bols")

    mf <- list(...)
    if (length(mf) == 1 && ((isMATRIX(mf[[1]]) || is.data.frame(mf[[1]])) &&
                            ncol(mf[[1]]) > 1 )) {
        mf <- mf[[1]]
        ### spline bases should be matrices
        if (isMATRIX(mf) && !is(mf, "Matrix"))
            class(mf) <- "matrix"
    } else {
        mf <- as.data.frame(mf)
        cl <- as.list(match.call(expand.dots = FALSE))[2][[1]]
        colnames(mf) <- sapply(cl, function(x) as.character(x))
    }
    if(!intercept && !any(sapply(mf, is.factor)) &&
       !any(sapply(mf, function(x){uni <- unique(x);
                                   length(uni[!is.na(uni)])}) == 1)){
        ## if no intercept is used and no covariate is a factor
        ## and if no intercept is specified (i.e. mf[[i]] is constant)
        if (any(sapply(mf, function(x) abs(mean(x, na.rm=TRUE) / sd(x,na.rm=TRUE))) > 0.1))
            ## if covariate mean is not near zero
            warning("covariates should be (mean-) centered if ",
                    sQuote("intercept = FALSE"))
    }
    vary <- ""
    if (!is.null(by)){
        stopifnot(is.data.frame(mf))
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
                get_call = function(){
                    cll <- deparse(cll, width.cutoff=500L)
                    if (length(cll) > 1)
                        cll <- paste(cll, collapse="")
                    cll
                },
                get_data = function() mf,
                get_index = function() index,
                get_names = function() colnames(mf),
                get_vary = function() vary,
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

    ret$dpp <- bl_lin(ret, Xfun = X_ols, args = hyper_ols(
                      df = df, lambda = lambda,
                      intercept = intercept, contrasts.arg = contrasts.arg))
    return(ret)
}

### P-spline (and tensor-product spline) baselearner
bbs <- function(..., by = NULL, index = NULL, knots = 20, boundary.knots = NULL,
                degree = 3, differences = 2, df = 4, lambda = NULL, center = FALSE,
                cyclic = FALSE, constraint = c("none", "increasing", "decreasing"),
                deriv = 0) {

    if (!is.null(lambda)) df <- NULL

    cll <- match.call()
    cll[[1]] <- as.name("bbs")

    mf <- list(...)
    if (length(mf) == 1 && ((is.matrix(mf[[1]]) || is.data.frame(mf[[1]])) &&
                            ncol(mf[[1]]) > 1 )) {
        mf <- as.data.frame(mf[[1]])
    } else {
        mf <- as.data.frame(mf)
        cl <- as.list(match.call(expand.dots = FALSE))[2][[1]]
        colnames(mf) <- sapply(cl, function(x) deparse(x))
    }
    stopifnot(is.data.frame(mf))
    if(!(all(sapply(mf, is.numeric)))) {
        if (ncol(mf) == 1){
            warning("cannot compute ", sQuote("bbs"),
                    " for non-numeric variables; used ",
                    sQuote("bols"), " instead.")
            return(bols(mf, by = by, index = index))
        }
        stop("cannot compute bbs for non-numeric variables")
    }
    vary <- ""
    if (!is.null(by)){
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

    ret$dpp <- bl_lin(ret, Xfun = X_bbs,
                      args = hyper_bbs(mf, vary, knots = knots, boundary.knots =
                      boundary.knots, degree = degree, differences = differences,
                      df = df, lambda = lambda, center = center, cyclic = cyclic,
                      constraint = match.arg(constraint), deriv = deriv))
    return(ret)
}

### cyclic B-splines
### adapted version of mgcv::cSplineDes from S.N. Wood
cbs <- function (x, knots, boundary.knots, degree = 3, deriv = 0L) {

    if (any(x < boundary.knots[1], na.rm = TRUE) |
        any(x > boundary.knots[2], na.rm = TRUE))
        stop("some ", sQuote("x"), " values are beyond ",
             sQuote("boundary.knots"))

    nx <- names(x)
    x <- as.vector(x)
    ## handling of NAs
    nax <- is.na(x)
    if (nas <- any(nax))
        x <- x[!nax]

    knots <- c(boundary.knots[1], knots, boundary.knots[2])
    nKnots <- length(knots)
    ord <- degree + 1
    xc <- knots[nKnots - ord + 1]
    knots <- c(boundary.knots[1] -
               (boundary.knots[2] - knots[(nKnots - ord + 1):(nKnots - 1)]),
               knots)
    ind <- x > xc
    X <- splineDesign(knots, x, ord, derivs = rep(deriv, length(x)), outer.ok = TRUE)
    x[ind] <- x[ind] - boundary.knots[2] + boundary.knots[1]
    if (sum(ind)) {
        Xtmp <- splineDesign(knots, x[ind], ord, derivs = rep(deriv, length(x[ind])),
                             outer.ok = TRUE)
        X[ind, ] <- X[ind, ] + Xtmp
    }
    ## handling of NAs
    if (nas) {
        tmp <- matrix(NA, length(nax), ncol(X))
        tmp[!nax, ] <- X
        X <- tmp
    }
    ## add attributes
    attr(X, "degree") <- degree
    attr(X,"knots") <- knots
    attr(X,"boundary.knots") <- boundary.knots
    if (length(deriv) > 1 || deriv != 0)
        attr(X, "deriv") <- deriv
    dimnames(X) <- list(nx, 1L:ncol(X))
    return(X)
}

bsplines <- function(x, knots, boundary.knots, degree,
                     Ts_constraint = "none", deriv = 0L,
                     extrapolation = FALSE) {

    ## do not allow data beyond boundary knots while fitting
    if (!extrapolation && (any(x < boundary.knots[1], na.rm = TRUE) |
                               any(x > boundary.knots[2], na.rm = TRUE)))
        stop("some ", sQuote("x"), " values are beyond ",
             sQuote("boundary.knots"))

    ## allow extrapolation when predicting
    if (extrapolation <- extrapolation &&
        (any(x < boundary.knots[1], na.rm = TRUE) |
             any(x > boundary.knots[2], na.rm = TRUE))) {
        warning("Some ", sQuote("x"), " values are beyond ",
                sQuote("boundary.knots"), "; Linear extrapolation used.")
    }

    nx <- names(x)
    x <- as.vector(x)
    ## handling of NAs
    nax <- is.na(x)
    if (nas <- any(nax))
        x <- x[!nax]
    ## use equidistant boundary knots
    dx <- diff(boundary.knots)/(length(knots) + 1)
    bk_lower <- seq(boundary.knots[1] - degree * dx, boundary.knots[1],
                    length = degree + 1)
    bk_upper <- seq(boundary.knots[2], boundary.knots[2] + degree * dx,
                    length = degree + 1)
    ## complete knot mesh
    k <- c(bk_lower, knots, bk_upper)
    ## construct design matrix
    X <- splineDesign(k, x, degree + 1, derivs = rep(deriv, length(x)),
                      outer.ok = TRUE)

    ## code along the lines of mgcv::Predict.matrix.pspline.smooth
    if (extrapolation) {
        ## Build matrix to map coeficients to value (deriv = 0) and
        ## slope (deriv = 1) at end points.
        if (deriv != 0L) {
            warning("deriv != 0L; Linear extrapolation overwritten")
        } else {
              deriv <- c(0, 1, 0, 1)
          }
        D <- splineDesign(knots = k, x = rep(boundary.knots, each = 2),
                          ord = degree + 1, deriv)
        ## Add rows for linear extrapolation
        idx <- x < boundary.knots[1]
        if (any(idx, na.rm = TRUE))
            X[idx,] <- cbind(1, x[idx] - boundary.knots[1]) %*% D[1:2, ]
        idx <- x > boundary.knots[2]
        if (any(idx, na.rm = TRUE))
            X[idx,] <- cbind(1, x[idx] - boundary.knots[2]) %*% D[3:4, ]
    }

    ## handling of NAs
    if (nas) {
        tmp <- matrix(NA, length(nax), ncol(X))
        tmp[!nax, ] <- X
        X <- tmp
    }
    ### constraints; experimental
    D <- diag(ncol(X))
    D[lower.tri(D)] <- 1
    X <- switch(Ts_constraint, "none" = X,
                            "increasing" = X %*% D,
                            "decreasing" = -X %*% D)
    ## add attributes
    attr(X, "degree") <- degree
    attr(X, "knots") <- knots
    attr(X, "boundary.knots") <- list(lower = bk_lower, upper = bk_upper)
    if (Ts_constraint != "none")
        attr(X, "Ts_constraint") <- Ts_constraint
    if (Ts_constraint != "none")
        attr(X, "D") <- D
    if (length(deriv) > 1 || deriv != 0)
        attr(X, "deriv") <- deriv
    dimnames(X) <- list(nx, 1L:ncol(X))
    return(X)
}

### workhorse for fitting (Ridge-penalized) baselearners
bl_lin <- function(blg, Xfun, args) {

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

    dpp <- function(weights) {

        if (!is.null(attr(X, "deriv")))
            stop("fitting of derivatives of B-splines not implemented")

        weights[!Complete.cases(mf)] <- 0
        w <- weights
        if (!is.null(index))
            w <- .Call("R_ysum", as.double(weights), as.integer(index), PACKAGE = "mboost")
        XtX <- crossprod(X * w, X)
        lambdadf <- df2lambda(X, df = args$df, lambda = args$lambda,
                              dmat = K, weights = w, XtX = XtX)
        lambda <- lambdadf["lambda"]
        XtX <- XtX + lambda * K

        ## matrizes of class dgeMatrix are dense generic matrices; they should
        ## be coerced to class matrix and handled in the standard way
        if (is(X, "Matrix") && !extends(class(XtX), "dgeMatrix")) {
            XtXC <- Cholesky(forceSymmetric(XtX))
            mysolve <- function(y) {
                if (is.null(attr(X, "Ts_constraint")))
                    return(solve(XtXC, crossprod(X, y)))  ## special solve routine from
                                                          ## package Matrix
                ### non-negative LS only at the moment
                return(nnls1D(as(XtX, "matrix"), as(X, "matrix"), y))
            }
        } else {
            if (is(X, "Matrix")) {
                ## coerce Matrix to matrix
                X <- as(X, "matrix")
                XtX <- as(XtX, "matrix")
            }
            mysolve <- function(y) {
                if (is.null(attr(X, "Ts_constraint")))
                    return(solve(XtX, crossprod(X, y), LINPACK = FALSE))
                ### non-negative LS only at the moment
                return(nnls1D(XtX, X, y))
            }
        }

        fit <- function(y) {
            if (!is.null(index)) {
                y <- .Call("R_ysum", as.double(weights * y), as.integer(index),
                           PACKAGE = "mboost")
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
            return(ret[index, index])
        }
        ### </FIXME>

        ### actually used degrees of freedom (trace of hat matrix)
        df <- function() lambdadf

        ### prepare for computing predictions
        predict <- function(bm, newdata = NULL, aggregate = c("sum", "cumsum", "none")) {
            cf <- sapply(bm, coef)
            if (!is.matrix(cf))
                cf <- matrix(cf, nrow = 1)
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
                return(pr[ , , drop = FALSE])
            return(pr[index, ,drop = FALSE])
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
bspatial <- function(..., df = 6) {
    cl <- cltmp <- match.call()
    if (is.null(cl$df)) cl$df <- df
    cl[[1L]] <- as.name("bbs")
    ret <- eval(cl, parent.frame())
    cltmp[[1]] <- as.name("bspatial")
    assign("cll", cltmp, envir = environment(ret$get_call))
    ret
}

### random-effects (Ridge-penalized ANOVA) baselearner
brandom <- function(..., by = NULL, index = NULL, df = 4, lambda = NULL,
                    contrasts.arg = "contr.dummy") {
    cl <- cltmp <- match.call()
    x <- list(...)
    ## drop further arguments to be passed to bols
    if (!is.null(names(x)))
        x <- x[names(x) == ""]

    if (!all(sapply(x, is.factor) |
             sapply(x, is.matrix) |
             sapply(x, is.data.frame)))
        stop(sQuote("..."), " must be a factor or design matrix in ",
             sQuote("brandom"))

    if (is.null(cl$df) && is.null(cl$lambda))
        cl$df <- df
    if (is.null(cl$contrasts.arg))
        cl$contrasts.arg <- contrasts.arg
    cl[[1L]] <- as.name("bols")
    ret <- eval(cl, parent.frame())
    cltmp[[1]] <- as.name("brandom")
    assign("cll", cltmp, envir = environment(ret$get_call))
    ret
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

    cll <- paste(bl1$get_call(), "%+%",
                 bl2$get_call(), collapse = "")
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
    if (!CC)
        warning("base-learner contains missing values;\n",
                "missing values are excluded per base-learner, ",
                "i.e., base-learners may depend on different",
                " numbers of observations.")
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
                ## <FIXME> Is this all we want to change if we set names here?
                set_names = function(value) attr(mf, "names") <<- value)
                ## </FIXME>
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

        X1 <- newX1(mf[, bl1$get_names(), drop = FALSE],
                    prediction = args$prediction)
        K1 <- X1$K
        if (!is.null(l1)) K1 <- l1 * K1
        X1 <- X1$X

        X2 <- newX2(mf[, bl2$get_names(), drop = FALSE],
                    prediction = args$prediction)
        K2 <- X2$K
        if (!is.null(l2)) K2 <- l2 * K2
        X2 <- X2$X

        K <- matrix(0, ncol = ncol(K1) + ncol(K2),
                    nrow = nrow(K1) + nrow(K2))
        K[1:nrow(K1), 1:ncol(K1)] <- as.matrix(K1)
        K[-(1:nrow(K1)), -(1:ncol(K1))] <- as.matrix(K2)
        X <- cbind(as.matrix(X1), as.matrix(X2))
        MATRIX <- any(dim(X) > c(500, 50)) &&
                  options("mboost_useMatrix")$mboost_useMatrix
        if (MATRIX & !is(X, "Matrix"))
            X <- Matrix(X)
        if (MATRIX & !is(K, "Matrix"))
            K <- Matrix(K)
        list(X = X, K = K)
    }

    ret$dpp <- bl_lin(ret, Xfun = Xfun, args = args)

    return(ret)
}

"%X%" <- function(bl1, bl2) {

    if (is.list(bl1) && !inherits(bl1, "blg"))
        return(lapply(bl1, "%X%", bl2 = bl2))

    if (is.list(bl2) && !inherits(bl2, "blg"))
        return(lapply(bl2, "%X%", bl1 = bl1))

    cll <- paste(bl1$get_call(), "%X%",
                 bl2$get_call(), collapse = "")
    stopifnot(inherits(bl1, "blg"))
    stopifnot(inherits(bl2, "blg"))

    stopifnot(!any(colnames(model.frame(bl1)) %in%
                   colnames(model.frame(bl2))))
    mf <- cbind(model.frame(bl1), model.frame(bl2))
    index1 <- bl1$get_index()
    index2 <- bl2$get_index()
    if (is.null(index1)) index1 <- 1:nrow(mf)
    if (is.null(index2)) index2 <- 1:nrow(mf)

    mfindex <- cbind(index1, index2)
    index <- NULL

    CC <- all(Complete.cases(mf))
    if (!CC)
        warning("base-learner contains missing values;\n",
                "missing values are excluded per base-learner, ",
                "i.e., base-learners may depend on different",
                " numbers of observations.")
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
                ## <FIXME> Is this all we want to change if we set names here?
                set_names = function(value) attr(mf, "names") <<- value)
                ## </FIXME>
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

        X1 <- newX1(mf[, bl1$get_names(), drop = FALSE],
                    prediction = args$prediction)
        K1 <- X1$K
        X1 <- X1$X
        if (!is.null(l1)) K1 <- l1 * K1
        MATRIX <- options("mboost_useMatrix")$mboost_useMatrix
        if (MATRIX & !is(X1, "Matrix"))
            X1 <- Matrix(X1)
        if (MATRIX & !is(K1, "Matrix"))
            K1 <- Matrix(K1)

        X2 <- newX2(mf[, bl2$get_names(), drop = FALSE],
                    prediction = args$prediction)
        K2 <- X2$K
        X2 <- X2$X
        if (!is.null(l2)) K2 <- l2 * K2
        if (MATRIX & !is(X2, "Matrix"))
            X2 <- Matrix(X2)
        if (MATRIX & !is(K2, "Matrix"))
            K2 <- Matrix(K2)
        suppressMessages(
            X <- kronecker(X1, Matrix(1, ncol = ncol(X2),
                                  dimnames = list("", colnames(X2))),
                       make.dimnames = TRUE) *
                kronecker(Matrix(1, ncol = ncol(X1),
                              dimnames = list("", colnames(X1))),
                       X2, make.dimnames = TRUE)
            )
        suppressMessages(
            K <- kronecker(K1, diag(ncol(X2))) +
                 kronecker(diag(ncol(X1)), K2)
            )
        list(X = X, K = K)
    }

    ret$dpp <- bl_lin(ret, Xfun = Xfun, args = args)

    return(ret)
}

bns <- function(...)
    stop("Base-learner ",  sQuote("bns"), " has ben removed. Consider ",
         sQuote("bbs"), " instead.")

bss <- function(...)
    stop("Base-learner ",  sQuote("bss"), " has ben removed. Consider ",
         sQuote("bbs"), " instead.")
