# Methods to orthogonalize predictors
#
# Authors: Christian Lindenlaub, Benjamin Hofner
###############################################################################

# ================================================= #
#  %ortho%-Operator		          	    #
# ================================================= #

### args:
# bl1: Base-learner 1
# bl2: Base-learner 2


"%ortho%" <- function(bl1, bl2){

    if(is.list(bl1) && !inherits(bl1, "blg"))
        return(lapply(bl1, "%ll%", bl2 = bl2))

    if(is.list(bl2) && !inherits(bl2, "blg"))
        return(lapply(bl2, "%ll%", bl1 = bl1))

    ## set baselearner name
    cll <- paste(bl1$get_call(), "%ll%",
                 bl2$get_call(), collapse = "")
    cll <- paste(bl1$get_call())

    ## test if baselearners
    stopifnot(inherits(bl1, "blg"))
    stopifnot(inherits(bl2, "blg"))

    ## build model.frame
    mf <- cbind(model.frame(bl1), model.frame(bl2))

    ## index
    index <- NULL

    ## vary
    vary <- ""

    ## return
    ret <- list(

        ## model.frame
        model.frame = function() mf,
        ## function
        get_call = function(){
            ##cll <- deparse(cll, width.cutoff = 500L)
            if (length(cll) > 1)
                cll <- paste(cll, collapse = "")
            cll
        },

        ## model.frame data
        get_data = function() mf,

        ## index
        get_index = function() index,
        get_vary = function() vary,

        ## get the names of the model.frame
        get_names = function() colnames(mf),

        ## change the names of the model.frame
        set_names = function(value) attr(mf, "names") <<- value
    )
    ## class return
    class(ret) <- "blg"

    ## read arguments
    args1 <- environment(bl1$dpp)$args
    args2 <- environment(bl2$dpp)$args

    ## lambda
    l1 <- args1$lambda
    l2 <- args2$lambda
    if (!is.null(l1) && !is.null(l2)){
        args <- list(lambda = 1, df = NULL)
    }
    else{
        args <- list(lambda = NULL,
                     df = ifelse(is.null(args1$df), 0, args1$df) +
                         ifelse(is.null(args2$df), 0, args2$df))
    }

    ## Xfun
    Xfun <- function(mf, vary, args){
        ## create x and k matrices
        newX1 <- environment(bl1$dpp)$newX
        newX2 <- environment(bl2$dpp)$newX
        ## extract x and k matrices
        X1 <- newX1(mf[, bl1$get_names(), drop = FALSE])
        K1 <- X1$K
        if (!is.null(l1)) K1 <- l1 * K1
        X1 <- X1$X

        X2 <- newX2(mf[, bl2$get_names(), drop = FALSE])
        K2 <- X2$K
        if (!is.null(l2)) K2 <- l2 * K2
        X2 <- X2$X

        ## make x1 orthogonal to x2
        # qr.resid(qr, y)
        # X1orth <- qr.resid(qr(X2), X1)
        # X1orth2 <-  (I - (X2 (X2'X2)^-1 X2') X1)
        # X1 <- qr.resid(qr(X2), X1) = I - (X2 (X2'X2)^-1 X2') X1

        ## only one is of type "Matrix"
        if (xor(is(X1, "Matrix"), is(X2, "Matrix"))) {
            X1 <- Matrix(X1)
            X2 <- Matrix(X2)
        }
        ## new design matrix X
        X <- qr.resid(qr(X2), X1)
        ## new penalty matrix K
        K <- K1

        ## return
        list(X = X, K = K)
    }
    ret$dpp <- bl_lin(ret, Xfun = Xfun, args = args)

    return(ret)
}
