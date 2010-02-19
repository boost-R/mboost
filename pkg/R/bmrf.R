bmrf <-
function (..., by = NULL, index = NULL, bnd = NULL, df = 4, lambda = NULL, 
    center = FALSE) 
{
    if (!require("BayesX"))
        stop("cannot load ", sQuote("BayesX"))

    cll <- match.call()
    cll[[1]] <- as.name("bmrf")
    cll <- deparse(cll)
    mf <- list(...)
    if (length(mf) == 1 && (is.matrix(mf[[1]]) || is.data.frame(mf[[1]]))) {
        mf <- as.data.frame(mf[[1]])
    }
    else {
        mf <- as.data.frame(mf)
        cl <- as.list(match.call(expand.dots = FALSE))[2][[1]]
        colnames(mf) <- sapply(cl, function(x) deparse(x))
    }
    stopifnot(is.data.frame(mf))
    if (!(all(sapply(mf, is.factor)))) {
        stop("cannot compute bmrf for numeric variables")
    }
    vary <- ""
    if (!is.null(by)) {
        stopifnot(is.numeric(by) || (is.factor(by) && nlevels(by) == 
            2))
        mf <- cbind(mf, by)
        colnames(mf)[ncol(mf)] <- vary <- deparse(substitute(by))
    }
    CC <- all(complete.cases(mf))
    DOINDEX <- (nrow(mf) > 10000)
    if (is.null(index)) {
        if (!CC || DOINDEX) {
            index <- get_index(mf)
            mf <- mf[index[[1]], , drop = FALSE]
            index <- index[[2]]
        }
    }
    ret <- list(model.frame = function() if (is.null(index)) return(mf) else return(mf[index, 
        , drop = FALSE]), get_call = function() cll, get_data = function() mf, 
        get_index = function() index, get_vary = function() vary, 
        get_names = function() colnames(mf), set_names = function(value) attr(mf, 
            "names") <<- value)
    class(ret) <- "blg"
    ret$dpp <- mboost:::bl_lin(ret, Xfun = X_bmrf, args = hyper_bmrf(mf, 
        vary, bnd = bnd, df = df, lambda = lambda, center = center))
    return(ret)
}

hyper_bmrf <-
function (mf, vary, bnd = NULL, df = 4, lambda = NULL, center = FALSE) 
{
    if (is.null(bnd)) 
        stop("Neighbourhood must be defined")
    else if (class(bnd) == "bnd") 
        K <- bnd2gra(bnd)
    else if (is.matrix(bnd) && dim(bnd)[1] == dim(bnd)[2] && 
        sum(bnd) == 0 && !is.null(dimnames(bnd)[[1]])) 
        K <- bnd
    else stop("No proper neighbourhood definition given")
    class(K) <- "matrix"
    nm <- colnames(mf)[colnames(mf) != vary]
    list(K = K, bnd = bnd, pen = TRUE, df = df, lambda = lambda, 
        center = center)
}

X_bmrf <-
function (mf, vary, args) 
{
    stopifnot(is.data.frame(mf))
    K = args$K
    districts = dimnames(K)[[1]]
    X = matrix(0, nrow = dim(mf)[1], ncol = dim(K)[2])
    for (i in 1:dim(mf)[1]) X[i, which(districts == mf[i, 1])] = 1
    MATRIX <- any(dim(X) > c(500, 50))
    MATRIX <- MATRIX && options("mboost_useMatrix")$mboost_useMatrix
    if (MATRIX && require(Matrix)) {
        diag <- Diagonal
        X <- Matrix(X)
    }
    if (vary != "") {
        by <- model.matrix(as.formula(paste("~", vary, collapse = "")), 
            data = mf)[, 2]
        X <- X * by
    }
    if (args$center) {
        e <- eigen(K)
        L <- e$vectors[, -dim(e$vectors)[2]] * sqrt(e$values[-length(e$values)])
        Zspathelp <- L %*% solve(t(L) %*% L)
        X <- X %*% Zspathelp
        K <- diag(ncol(X))
    }
    return(list(X = X, K = K))
}
