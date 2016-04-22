###
# User-specified base-learner with quadratic penalty
buser <- function(X, K = NULL,  by = NULL, index = NULL, df = 4, lambda = NULL){

    ## TODO:
    ## index should be available
    ## check if CC and index work

    if (!is.null(lambda)) df <- NULL
    cll <- match.call()
    cll[[1]] <- as.name("buser")

    mf <- as.data.frame(X) ##<FIXME> is this correct this way?

    vary <- ""
    if (!is.null(by)){
        #stopifnot(is.data.frame(mf))
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

    ret$dpp <- bl_lin(ret, Xfun = X_user,
                      args = hyper_user(mf, vary, K, df, lambda))
    return(ret)
}


### model.matrix for arbitrary user-specified base-learners
X_user <- function(mf, vary, args) {
    X <- mf
    K <- args$K
    if (vary != "") {
        by <- model.matrix(as.formula(paste("~", vary, collapse = "")),
                           data = as.data.frame(mf))[ , -1, drop = FALSE]
        X <- X[,colnames(mf) != vary]
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
        if (ncol(by) > 1)
            suppressMessages(K <- kronecker(diag(ncol(by)), K))
    }
    X <- as.matrix(X)
    # contr <- NULL ##<FIXME> Do we need this?
    return(list(X = X, K = K))
}

### hyper parameters for arbitrary user-specified base-learners
hyper_user <- function(mf, vary, K, df, lambda) {
    if (is.null(K)){
        K <- diag(ncol(mf[,colnames(mf) != vary]))
        lambda <- 0
        df <- NULL
    }
    list(df = df, lambda = lambda, K = K)
}
