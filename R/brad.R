brad <- function(..., by = NULL, index = NULL, knots = 100, df = 4, lambda = NULL,
                 covFun = fields::stationary.cov,
                 args = list(Covariance = "Matern", smoothness = 1.5, theta = NULL)) {

    cll <- match.call()
    cll[[1]] <- as.name("brad")

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
        stop("cannot compute brad for non-numeric variables")
    }
    ### use bols when appropriate
    ## <FIXME>
    #if (!is.null(df)) {
    #    if (df <= ncol(mf))
    #        return(bols(..., by = by, index = index))
    #}
    ## </FIXME>
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

    ret$dpp <- bl_lin(ret, Xfun = X_brad,
                      args = hyper_brad(mf, vary, knots = knots,
                      df = df, lambda = lambda, covFun = covFun, args = args))
    return(ret)
}


### model.matrix for kriging base-learners
X_brad <- function(mf, vary, args) {
    stopifnot(is.data.frame(mf))

    INDX <- which(colnames(mf) != vary)
    #if (length(INDX) != 2)
    #    stop("only bivariate kriging implemented atm.")
    args$args$x1 <- mf[INDX]
    args$args$x2 <- args$knots
    X <- do.call(args$covFun, args$args)
    args$args$x1 <- args$knots
    PEN <- do.call(args$covFun, args$args)
    e <- eigen(PEN)
    PEN_sqrt_INV <- solve(e$vectors %*% diag(sqrt(e$values)) %*% t(e$vectors))
    X <- X %*% PEN_sqrt_INV

    K <- diag(ncol(X))
    ### <FIXME>
    if (vary != "") {
        by <- model.matrix(as.formula(paste("~", vary, collapse = "")), data = mf)[,2]
        X <- X * by
    }
    ### </FIXME>
    attr(X, "knots") <- args$knots
    return(list(X = X, K = K))
}

### hyper parameters for kriging base-learners
hyper_brad <- function(mf, vary, knots = 100, df = 4, lambda = NULL,
                       covFun = covFun, args = args) {

    ## first we need to build a correct matrix of mf
    x <- as.matrix(mf[which(colnames(mf) != vary)])
    if (length(knots) == 1) {
        if (!requireNamespace("fields"))
            stop("Cannot load package", sQuote("fields"),
                 ", which is needed for the automatic knot placement")
        knots <- fields::cover.design(R = unique(x), nd = knots)$design
    }
    if ("theta" %in% names(args) && is.null(args$theta)){
        ## (try to) compute effective range
        args$theta <- effective_range(x, eps = 0.001, interval = c(0.1, 100),
                                      covFun = covFun, args = args)
    }
    if (is.list(knots))
        stop(sQuote("knots"), " must be an integer defining the number of knots",
             " or a matrix specifiying the location of the knots")
    list(knots = knots, pen = TRUE, df = df, lambda = lambda, covFun = covFun,
         args = args)
}

######
# theta: effective range (per default theta such that
#       rho( max(x_(i) - x_(j)), smoothness, theta = max(x_(i) - x_(j))/c ) = 0.001
#  <==> rho(c, smoothness, theta = 1) = 0.001
effective_range <- function(x, eps = 0.001, interval = c(0.1, 100),
                            covFun = fields::stationary.cov, args = list()){

    if ( !( length(deparse(covFun)) == length(deparse(fields::stationary.cov))
           && all(deparse(covFun) == deparse(fields::stationary.cov)) ) &&
         !( length(deparse(covFun)) == length(deparse(fields::Exp.cov))
           && all(deparse(covFun) == deparse(fields::Exp.cov)) ) ){
        ## if cov.funcion is not one of stationary.cov and Exp.cov
        warning(sQuote("effective_range()"), " is only implemented for ",
                sQuote("stationary.cov"), " and ", sQuote("Exp.cov"),
                " from package ", sQuote("fields"))
        return(NULL)
    }

    args$theta <- 1
    args$x1 <- 0
    maxX <- max(dist(x))
    rho <- function(cval, eps){
        args$x2 <- cval
        # make a call to covFun with the corresponding arguments
        return(c(do.call(covFun, args) - eps))
    }
    cval <- uniroot(rho, interval = interval, eps = eps)$root
    RET <- maxX/cval
    attr(RET, "c_value") <- cval
    return(RET)
}
