bkernel <- function(..., df = 4, lambda = NULL, kernel = kernel.lin,
                    pathway = NULL, args = list()) {

    cll <- match.call()
    cll[[1]] <- as.name("bkernel")

    mf <- list(...)
    if (length(mf) > 1 || !inherits(mf[[1]], "databel"))
        stop(sQuote("..."), " must be a single object of class ",
             sQuote("databel"))

    mf <- mf[[1]]
    if (is.null(attr(mf, "anno")))
        stop("SNP data needs annotation as ",
             sQuote('attr(, "anno")'))

    ### all data checks should be done in kangar00!

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
### <FIXME> Do we really want these colnames or simply the name of ...?
                get_names = function() colnames(mf@geno),
                set_names = function(value) {
                    if(length(value) != length(colnames(mf@geno)))
                        stop(sQuote("value"), " must have same length as ",
                             sQuote("colnames(mf)"))
                    for (i in 1:length(value)){
                        cll[[i+1]] <<- as.name(value[i])
                    }
                    attr(mf@geno, "names") <<- value
                })
## </FIXME>@geno
    class(ret) <- "blg"

    ret$dpp <- bl_lin(ret, Xfun = X_kernel,
                      args = hyper_kernel(df = df, lambda = lambda,
                                          kernel = kernel, pathway = pathway,
                                          args = args))
    return(ret)
}


### model.matrix for kriging base-learners
X_kernel <- function(mf, vary, args) {

    ## add GWASdata to args (in order to use do.call)
    args$args$data <- mf
    ## compute kernel with arguments args$args
    X <- K <- do.call(args$kernel, args$args)@kernel

    e <- eigen(X)

    PEN_sqrt_INV <- solve(e$vectors %*% diag(sqrt(e$values)) %*% t(e$vectors))
    X <- X %*% PEN_sqrt_INV
    K <- diag(ncol(X))
    return(list(X = X, K = K))
}

### hyper parameters for kernel base-learners
hyper_kernel <- function(df = 4, lambda = NULL, kernel = kernel.lin,
                       pathway = NULL, args = list()) {
    list(df = df, lambda = lambda, kernel = kernel, args = c(pathway = pathway, args))
}

### <FIXME> Add possibility to specify ANA'
### <FIXME> What about newdata in bl_lin?
