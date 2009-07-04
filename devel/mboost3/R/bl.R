
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


hyper_ols <- function(mf, vary) list(center = FALSE)

mm_ols <- function(mf, vary, args) {

    fm <- as.formula(paste("~ ", paste(colnames(mf)[colnames(mf) != vary], 
                     collapse = "+"), sep = ""))
    X <- model.matrix(fm, data = mf)
    if (vary != "")
        X <- X * mf[, vary]
    list(mm = X, K = NULL)
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
        Matrix(X)
    })
    if (length(mm) == 1) {
        X <- mm[[1]]
        K <- diff(Diagonal(ncol(X)), differences = args$differences)
        K <- crossprod(K, K)
    }
    if (length(mm) == 2) {
        X <- kronecker(mm[[1]], matrix(1, nc = ncol(mm[[2]]))) * 
             kronecker(matrix(1, nc = ncol(mm[[1]])), mm[[2]])
        Kx <- diff(Diagonal(ncol(mm[[1]])), differences = args$differences)
        Kx <- crossprod(Kx, Kx)
        Ky <- diff(Diagonal(ncol(mm[[2]])), differences = args$differences)
        Ky <- crossprod(Ky, Ky)
        K <- kronecker(Kx, Diagonal(ncol(mm[[2]]))) + 
             kronecker(Diagonal(ncol(mm[[1]])), Ky)
    }
    if (vary != "")
        X <- X * mf[, vary]
    return(list(mm = X, K = K))
}

bl <- function(..., z = NULL, index = NULL, hfun = hyper_ols, Xfun = mm_ols) {

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

    if (is.null(index) & !is.matrix(mf)) {
        index <- get_index(mf)
        mf <- mf[index[[1]],,drop = FALSE]
        index <- index[[2]]
    }

    ret <- list(model.frame = function() 
                    if (is.null(index)) return(mf) else return(mf[index,,drop = FALSE]),
                get_names = function() colnames(mf),
                set_names = function(value) attr(mf, "names") <<- value)
    class(ret) <- "blg"

    args <- hfun(mf, vary)

    newX <- function(newdata = NULL) {
        if (!is.null(newdata) && !all(names(newdata) == names(mf)))
            mf <- newdata[,colnames(mf),drop = FALSE]        
        if (is.matrix(mf)) {
            X <- mf[,colnames(mf) != vary,drop = FALSE]
        } else {
            X <- Xfun(mf, vary, args)
            K <- X$K
            X <- X$mm
        }
        return(list(X = Matrix(X), K = K))
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
        
        fit <- function(y) {
            if (!is.null(index)) {
                y <- as.vector(tapply(weights * y, index, sum))
            } else {
                y <- y * weights
            }
            coef <- solve(XtX, crossprod(X, y))
            ret <- list(model = as.vector(coef), 
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

        predict <- function(bm, newdata = NULL, Sum = TRUE) {
            cf <- sapply(bm, coef)
            if(!is.null(newdata)) {
                index <- get_index(mf)
                mf <- mf[index[[1]],,drop = FALSE]
                index <- index[[2]]
                X <- newX(mf)$X
            }
            if (Sum) {
                pr <- X %*% rowSums(cf)
            } else {
                M <- triu(crossprod(Matrix(1, nc = ncol(cf))))
                pr <- X %*% (cf %*% M)
            }
            if (is.null(index)) return(pr[,,drop = TRUE])
            return(pr[index,,drop = TRUE])
        }

        ret <- list(fit = fit, hatvalues = hatvalues, predict = predict)
        class(ret) <- c("bl_lin", "bl")
        return(ret)

    }
    ret$dpp <- dpp
    return(ret)
}

names.blg <- function(x)
    x$get_names()

model.frame.blg <- function(formula)
    formula$model.frame()

coef.bm_lin <- function(object)
    object$model

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

