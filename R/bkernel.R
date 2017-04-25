bkernel <- function(..., df = 4, lambda = NULL,
                    kernel = c("lin", "sia", "net"),
                    pathway = NULL, knots = NULL, args = list()) {

    if (!requireNamespace("kangar00", quietly = TRUE))
        stop("Please install package 'kangar00' for using bkernel base-learners.")

    cll <- match.call()
    cll[[1]] <- as.name("bkernel")

    mf <- list(...)
    if (length(mf) > 1 || !inherits(mf[[1]], "GWASdata"))
        stop(sQuote("..."), " must be a single object of class ",
             sQuote("GWASdata"), " as defined in package 'kangar00'")
    mf <- mf[[1]]

    kernel <- match.arg(kernel)

    ret <- list(model.frame = function() mf,
                get_call = function(){
                    cll <- deparse(cll, width.cutoff=500L)
                    if (length(cll) > 1)
                        cll <- paste(cll, collapse="")
                    cll
                },
                get_data = function() mf,
                get_index = function() NULL,
                get_vary = function() "",
### <FIXME> This now returns the name of the GWASdata object and not the column
###         names of mf@geno.
                get_names = function() as.character(cll[2]),
                set_names = function(value) {
                    if(length(value) != 1)
                        stop(sQuote("value"), " must have same length 1.")
                    cll[[2]] <<- as.name(value)
                })
## </FIXME>
    class(ret) <- "blg"

    ret$dpp <- bl_lin(ret, Xfun = X_kernel,
                      args = hyper_kernel(df = df, lambda = lambda,
                                          kernel = kernel, pathway = pathway,
                                          knots = if(!is.null(knots)) knots else mf,
                                          args = args))
    return(ret)
}

## - Think about storing all but X (e.g. (ANA')X') in args?
##   As we have different kernels, we would need to return this from
##   calc_kernel(), then we can use this.

## <FIXME> Slot @anno should be an hyper argument

### model.matrix for kernel base-learners
X_kernel <- function(mf, vary, args) {

    if (length(args$args) > 0)
        stop("further arguments are curently not implemented")

    ## compute kernel with arguments args$args
    X <- kangar00::calc_kernel(object = mf,
                               type = args$kernel,
                               pathway = args$pathway,
                               knots = args$knots)@kernel
    
    ## compute penalty matrix for prediction or if X matrix is not quadratic
    if (args$prediction || nrow(X) != ncol(X)) {
        K <- kangar00::calc_kernel(object = args$knots,
                                   type = args$kernel,
                                   pathway = args$pathway,
                                   knots = args$knots)@kernel
    } else {
        K <- X
    }
    ## make sure K is positive semidefinit
    K <- make_psd(K)

    ## compute inverse of the square root of the penalty matrix (always use knots)
    e <- eigen(K)
    if (any(e$values < 0)) {
        warning("Kernel should be positiv semidefinit, yet some eigenvalues are < 0.\n",
                "These eigenvalues were set to 0.")
        e$values[e$values < 0] <- 0
    }

    ## t(e$vectors) can be used instead of solve(e$vectors) as K is symmetric
    sqrtK <- e$vectors %*% diag(sqrt(e$values)) %*% t(e$vectors)
    sqrtK <- forceSymmetric(sqrtK)
    ## compute inverse of sqrt(K):
    ## chol2inv works for symmetric, posdef matrices and is much quicker than
    ## solve. Note that its argument is chol(...).
    chol_sqrtK <- try(chol(sqrtK), silent = TRUE)
    if (inherits(chol_sqrtK, "try-error")) {
        PEN_sqrt_INV <- solve(sqrtK)
    } else {
        PEN_sqrt_INV <- chol2inv(chol_sqrtK)
    }

    X <- X %*% PEN_sqrt_INV
    K <- diag(ncol(X))
    return(list(X = X, K = K))
}

### hyper parameters for kernel base-learners
# knots are used for prediction
hyper_kernel <- function(df = 4, lambda = NULL, kernel = "lin",
                         knots, pathway = NULL, args = list()) {
    ## knots takes the observations which where used at the time of fitting.
    ## needed for predictions.
    list(df = df, lambda = lambda, kernel = kernel, knots = knots,
         pathway = pathway, args = args)
}

### <FIXME> Add possibility to directly specify ANA'
