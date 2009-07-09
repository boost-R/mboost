
df2lambda <- function(X, df = 4, dmat = NULL, weights) {

#   if (df <= 2) stop(sQuote("df"), " must be greater than two")

    if (is.null(dmat)) {
        dmat <- diff(diag(ncol(X)), differences = 2)
        dmat <- crossprod(dmat, dmat)
    }

    # Cholesky decomposition

    A <- crossprod(X * weights, X) + dmat*10e-10
    Rm <- solve(chol(A))

    decomp <- svd(crossprod(Rm,dmat)%*%Rm)
    d <- decomp$d

    # df2lambda
    df2l <- function(lambda)
        (sum( 1/(1+lambda*d) ) - df)^2

    lower.l <- 0
    upper.l <- 5000
    lambda <- upper.l

    while (lambda >= upper.l - 200 ) {
        upper.l <- upper.l * 1.5

        tl <- try(lambda <- optimize(df2l, interval=c(lower.l,upper.l))$minimum,
        silent=T)
        if (class(tl)=="try-error") stop("problem of
        converting df into lambda cannot be solved - please increase value of
        df")
        lower.l <- upper.l-200
        if (lower.l > 1e+06){
            lambda <- 1e+06
            warning("lambda needs to be larger than 1e+06 for given value of df,
            setting lambda = 1e+06 \n trace of hat matrix differs from df by ",
            round(sum( 1/(1+lambda*d) )-df,6))
            break
            }
    }

    ### tmp <- sum(diag(X %*% solve(crossprod(X * weights, X) +
    ###                   lambda*dmat) %*% t(X * weights))) - df
    ### if (abs(tmp) > sqrt(.Machine$double.eps))
    ###   warning("trace of hat matrix is not equal df with difference", tmp)

    lambda
}

hyper_ols <- function(df = NULL, center = FALSE, contrasts.arg = "contr.treatment", ...) 
    list(df = df, center = FALSE, contrasts.arg = contrasts.arg, ...)

mm_ols <- function(mf, vary, args) {

    ANOVA <- (length(mf) == 1 && is.factor(mf[[1]]))
    fm <- paste("~ ", paste(colnames(mf)[colnames(mf) != vary], 
                collapse = "+"), sep = "")
    if (ANOVA || args$center)
        fm <- paste(fm, "-1", collapse = "")
    X <- model.matrix(as.formula(fm), data = mf, contrasts.arg = args$contrasts.arg)
    if (vary != "") {
        z <- model.matrix(as.formula(paste("~", vary, "-1", collapse = "")), data = mf)
        xtmp <- NULL
        for (i in 1:ncol(z))
            xtmp <- cbind(xtmp, X * z[,i])
        X <- xtmp
    }
    K <- NULL
    if (!is.null(args$df)) {
        if (ANOVA) { 
            diag <- Diagonal
            X <- Matrix(X)
        }
        if (ANOVA && is.ordered(mf[[1]])) K <- diff(diag(ncol(X)), differences = 2)
        K <- diag(ncol(X))
    }
    list(mm = X, K = K)
}

hyper_bbs <- function(mf, vary, knots = 20, degree = 3, differences = 2, df = 4) {

    knotf <- function(x, knots) {	
        boundary.knots <- range(x, na.rm = TRUE)
        if (length(knots) == 1) {
            knots <- seq(from = boundary.knots[1], 
                         to = boundary.knots[2], length = knots + 2)
            knots <- knots[2:(length(knots) - 1)]
        }
        list(knots = knots, boundary.knots = boundary.knots)
    }
    ret <- lapply(which(colnames(mf) != vary), function(i) {
        knotf(mf[[i]], if (is.list(knots)) knots[[i]] else knots)
    })
    list(knots = ret, degree = degree, differences = differences, df = df)
}

mm_bbs <- function(mf, vary, args) {

    mm <- lapply(which(colnames(mf) != vary), function(i) {
        X <- bs(mf[[i]], knots = args$knots[[i]]$knots, degree = args$degree,
           Boundary.knots = args$knots[[i]]$boundary.knots, intercept = TRUE)
        class(X) <- "matrix"
        if (nrow(X) > 500 || (ncol(X) > 50 || sum(colnames(mf) != vary) > 1))
            return(Matrix(X))
        return(X)
    })
    if (nrow(mm[[1]]) > 500 || (ncol(mm[[1]]) > 50 || length(mm) > 1))
        diag <- Diagonal
    if (length(mm) == 1) {
        X <- mm[[1]]
        K <- diff(diag(ncol(X)), differences = args$differences)
        K <- crossprod(K, K)
    }
    if (length(mm) == 2) {
        X <- kronecker(mm[[1]], matrix(1, nc = ncol(mm[[2]]))) * 
             kronecker(matrix(1, nc = ncol(mm[[1]])), mm[[2]])
        Kx <- diff(diag(ncol(mm[[1]])), differences = args$differences)
        Kx <- crossprod(Kx, Kx)
        Ky <- diff(diag(ncol(mm[[2]])), differences = args$differences)
        Ky <- crossprod(Ky, Ky)
        K <- kronecker(Kx, diag(ncol(mm[[2]]))) + 
             kronecker(diag(ncol(mm[[1]])), Ky)
    }
    if (vary != "") {
        z <- model.matrix(as.formula(paste("~", vary, "-1", collapse = "")), data = mf)[,-1]
        X <- X * z
#        xtmp <- NULL
#        for (i in 1:ncol(z))
#            xtmp <- cbind(xtmp, X * z[,i])
#        X <- xtmp
    }
    return(list(mm = X, K = K))
}

bols3 <- function(..., z = NULL, index = NULL, center = FALSE, df = NULL, 
                  contrasts.arg = "contr.treatment") {

    mf <- list(...)
    if (length(mf) == 1 && (is.matrix(mf[[1]]) || is.data.frame(mf[[1]]))) {
        mf <- mf[[1]]
    } else {
        mf <- as.data.frame(mf)
        cl <- as.list(match.call(expand.dots = FALSE))[2][[1]]
        colnames(mf) <- sapply(cl, function(x) as.character(x))
    }
    vary <- ""
    if (!is.null(z)) {
        mf <- cbind(mf, z)
        colnames(mf) <- c(colnames(mf), deparse(substitute(z)))
        vary <- colnames(mf)[ncol(mf)]
    }

    if (is.null(index) & (!is.matrix(mf) & nrow(mf) > 10000)) {
        index <- get_index(mf)
        mf <- mf[index[[1]],,drop = FALSE]
        index <- index[[2]]
    }

    ret <- list(model.frame = function() 
                    if (is.null(index)) return(mf) else return(mf[index,,drop = FALSE]),
                get_names = function() colnames(mf),
                set_names = function(value) attr(mf, "names") <<- value)
    class(ret) <- "blg"

    ret$dpp <- bl_lin(mf, vary, index = index, Xfun = mm_ols, args = hyper_ols(
                      df = df, center = center, contrasts.arg = contrasts.arg))
    return(ret)
}

bbs3 <- function(..., z = NULL, index = NULL, knots = 20, degree = 3, 
                 differences = 2, df = 4) {

    mf <- list(...)
    if (length(mf) == 1 && (is.matrix(mf[[1]]) || is.data.frame(mf[[1]]))) {
        mf <- mf[[1]]
    } else {
        mf <- as.data.frame(mf)
        cl <- as.list(match.call(expand.dots = FALSE))[2][[1]]
        colnames(mf) <- sapply(cl, function(x) as.character(x))
    }
    vary <- ""
    if (!is.null(z)) {
        mf <- cbind(mf, z)
        colnames(mf)[ncol(mf)] <- vary <- deparse(substitute(z))
    }

    if (is.null(index) & (!is.matrix(mf) & nrow(mf) > 100)) {
        index <- get_index(mf)
        mf <- mf[index[[1]],,drop = FALSE]
        index <- index[[2]]
    }

    ret <- list(model.frame = function() 
                    if (is.null(index)) return(mf) else return(mf[index,,drop = FALSE]),
                get_names = function() colnames(mf),
                set_names = function(value) attr(mf, "names") <<- value)
    class(ret) <- "blg"

    ret$dpp <- bl_lin(mf, vary, index = index, Xfun = mm_bbs, 
                      args = hyper_bbs(mf, vary, knots = knots,
                      degree = degree, differences = differences, df = df))
    return(ret)
}

bl_lin <- function(mf, vary, index = NULL, Xfun, args) {

    newX <- function(newdata = NULL) {
        if (!is.null(newdata) && all(names(newdata) == names(mf)))
            mf <- newdata[,colnames(mf),drop = FALSE]        
        if (is.matrix(mf)) {
            X <- mf[,colnames(mf) != vary,drop = FALSE]
        } else {
            X <- Xfun(mf, vary, args)
            K <- X$K
            X <- X$mm
        }
        return(list(X = X, K = K))
    }
    X <- newX()
    K <- X$K
    X <- X$X

    dpp <- function(weights) {

        weights[!complete.cases(mf)] <- 0
        w <- weights
        if (!is.null(index)) w <- as.vector(tapply(weights, index, sum))
        XtX <- crossprod(X * w, X)
        if (!is.null(args$df)) {
            lambda <- df2lambda(X, df = args$df, dmat = K, weights = w)
            XtX <- XtX + lambda * K
        }

        if (inherits(X, "Matrix")) {
            mysolve <- solve
        } else {
            mysolve <- function(a, b) 
                .Call("La_dgesv", a, b, .Machine$double.eps, PACKAGE = "base")
        }
        
        fit <- function(y) {
            if (!is.null(index)) {
                y <- as.vector(tapply(weights * y, index, sum))
            } else {
                y <- y * weights
            }
            coef <- mysolve(XtX, crossprod(X, y))
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
            ret <- as.matrix(tcrossprod(X %*% solve(XtX), X * w))
            if (is.null(index)) return(ret)
            return(ret[index,index])
        }

        predict <- function(bm, newdata = NULL, aggregate = c("sum", "cumsum", "none")) {
            cf <- sapply(bm, coef)
            if(!is.null(newdata)) {
                index <- NULL
                nm <- colnames(mf)
                if (vary != "") nm <- c(nm, vary)
                mf <- newdata[,nm, drop = FALSE]
                if (nrow(newdata) > 1000) {
                    index <- get_index(mf)
                    mf <- mf[index[[1]],,drop = FALSE]
                    index <- index[[2]]
                }
                X <- newX(mf)$X
            }
            aggregate <- match.arg(aggregate)
            pr <- switch(aggregate, "sum" = 
                as(X %*% rowSums(cf), "matrix"),
            "cumsum" = {
                M <- triu(crossprod(Matrix(1, nc = ncol(cf))))
                as(X %*% (cf %*% M), "matrix")
            },
            "none" = as(X %*% cf, "matrix"))
            if (is.null(index)) return(pr[,,drop = FALSE])
            return(pr[index,,drop = FALSE])
        }

        ret <- list(fit = fit, hatvalues = hatvalues, predict = predict)
        class(ret) <- c("bl_lin", "bl")
        return(ret)

    }
    return(dpp)
}

bspatial3 <- function(...) {
    cl <- match.call()
    cl[[1L]] <- as.name("bbs3")
    eval(cl, parent.frame())
}

brandom3 <- function(..., df = 4) {
    cl <- match.call()
    if (is.null(cl$df)) cl$df <- df
    cl$center <- FALSE
    cl[[1L]] <- as.name("bols3")
    eval(cl, parent.frame())
}


names.blg <- function(x)
    x$get_names()

model.frame.blg <- function(formula)
    formula$model.frame()

coef.bm_lin <- function(object)
    as.vector(object$model)

fitted.bm_lin <- function(object)
    object$fitted()

hatvalues.bl_lin <- function(model)
    model$hatvalues()

dpp <- function(object, weights)
    UseMethod("dpp", object)

fit <- function(object, y)
    UseMethod("fit", object)

dpp.blg <- function(object, weights)
    object$dpp(weights)

fit.bl <- function(object, y)
    object$fit(y)

