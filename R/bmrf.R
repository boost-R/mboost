bmrf <- function (..., by = NULL, index = NULL, bnd = NULL, df = 4,
                  lambda = NULL, center = FALSE) {
    if (!requireNamespace("BayesX"))
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
        stop("cannot compute bmrf for numeric variables, region variable must be given as factor.")
    }
    vary <- ""
    if (!is.null(by)) {
        stopifnot(is.numeric(by) || (is.factor(by) && nlevels(by) ==
            2))
        mf <- cbind(mf, by)
        colnames(mf)[ncol(mf)] <- vary <- deparse(substitute(by))
    }
    CC <- all(complete.cases(mf))
    if (!CC)
        warning("base-learner contains missing values;\n",
                "missing values are excluded per base-learner, ",
                "i.e., base-learners may depend on different",
                " numbers of observations.")
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
    ret$dpp <- bl_lin(ret, Xfun = X_bmrf, args = hyper_bmrf(mf,
        vary, bnd = bnd, df = df, lambda = lambda, center = center))
    return(ret)
}

hyper_bmrf <- function (mf, vary, bnd = NULL, df = 4, lambda = NULL,
                        center = FALSE) {
    if (is.null(bnd))
        stop("Neighbourhood relations must be given in matrix or boundary format.")
    else if (inherits(bnd, "bnd"))
        K <- BayesX::bnd2gra(bnd)
    else if (isMATRIX(bnd) &&
             nrow(bnd) == ncol(bnd) &&
             nlevels(mf[[1]]) <= nrow(bnd) &&
             sum(bnd) == 0 &&
             all(levels(mf[[1]]) %in% rownames(bnd)))
        K <- bnd
    else stop("Neighbourhood matrix not defined as stated in manual page.")
    ###    K <- as(K, "matrix")
    if (!is(K, "Matrix"))
        K <- Matrix(unclass(K))
    nm <- colnames(mf)[colnames(mf) != vary]
    list(K = K, bnd = bnd, pen = TRUE, df = df, lambda = lambda,
         center = center)
}

X_bmrf <- function (mf, vary, args) {
    K <- args$K
    districts <- rownames(K)
    X <- Diagonal(nrow(K))
    rownames(X) <- districts
    colnames(X) <- districts
    X <- X[as.character(mf[,1]),]
    # X <- matrix(0, nrow = nrow(mf), ncol = ncol(K))
    # for (i in 1:nrow(mf)) X[i, which(districts == mf[i, 1])] <- 1
    # MATRIX <- any(dim(X) > c(500, 50))
    # MATRIX <- MATRIX && options("mboost_useMatrix")$mboost_useMatrix
    # if (MATRIX) {
    #    diag <- Diagonal
    #     X <- Matrix(X)
    # }
    if (vary != "") {
        by <- model.matrix(as.formula(paste("~", vary, collapse = "")),
            data = mf)[, 2]
        X <- X * by
    }
    if (isTRUE(args$center)) {
        ### L = \Gamma \Omega^1/2 in Section 2.3. of Fahrmeir et al.
        ### (2004, Stat Sinica)
        SVD <- eigen(K)
        ev <- SVD$vectors[, -dim(SVD$vectors)[2]]
        ew <- SVD$values[-length(SVD$values)]
        X <- as(X %*% ev %*% diag(1/sqrt(ew)), "matrix")
        K <- as(diag(ncol(X)), "matrix")
    }
    return(list(X = X, K = K))
}
