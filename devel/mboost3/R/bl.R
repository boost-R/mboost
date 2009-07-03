
library("splines")
library("Matrix")
library("mboost")

get_index <- function(x) {
    
    tmp <- do.call("paste", x)           
    nd <- which(!duplicated(tmp))
    nd <- nd[complete.cases(x[nd,])]
    index <- match(tmp, tmp[nd])
    return(list(nd, index))
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
                    if (is.null(index)) return(mf) else return(mf[index,]),
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
            lambda <- mboost:::df2lambda(X, df = args$df, dmat = K, weights = w)
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
    object$hatvalues()

dpp <- function(object, weights)
    UseMethod("dpp", object)

fit <- function(object, y)
    UseMethod("fit", object)

dpp.blg <- function(object, weights)
    object$dpp(weights)

fit.bl <- function(object, y)
    object$fit(y)


mboost <- function(blg, y, weights = NULL, family = GaussReg(), control = boost_control()) {

    mboost:::check_y_family(y, family)

    ### hyper parameters
    mstop <- control$mstop
    risk <- control$risk
    constraint <- control$constraint
    nu <- control$nu
    trace <- control$trace
    tracestep <- options("width")$width / 2

    ### extract negative gradient and risk functions
    ngradient <- family@ngradient
    riskfct <- family@risk

    ### unweighted problem
    if (is.null(weights)) weights <- rep.int(1, length(y))
    WONE <- (max(abs(weights - 1)) < .Machine$double.eps)
    if (!family@weights && !WONE)
        stop(sQuote("family"), " is not able to deal with weights")

    ### rescale weights (because of the AIC criterion)
    ### <FIXME> is this correct with zero weights??? </FIXME>
    weights <- mboost:::rescale_weights(weights)
    oobweights <- as.numeric(weights == 0)

    bl <- lapply(blg, dpp, weights = weights)

    xselect <- integer(mstop)
    ens <- vector(mode = "list", length = mstop)

    ### vector of empirical risks for all boosting iterations
    ### (either in-bag or out-of-bag)
    mrisk <- numeric(mstop)
    mrisk[1:mstop] <- NA
    tsums <- numeric(length(bl))
    ss <- vector(mode = "list", length = length(bl))

    fit <- offset <- family@offset(y, weights)
    u <- ustart <- ngradient(y, fit, weights)

    ### start boosting iteration
    for (m in 1:mstop) {

        ### fit least squares to residuals _componentwise_
        for (i in 1:length(bl)) {
            tsums[i] <- -1
            ss[[i]] <- try(fit(bl[[i]], y = u))
            if (inherits(ss[[i]], "try-error")) next
            tsums[i] <- mean(weights * (fitted(ss[[i]]) - u)^2, na.rm = TRUE)
        }

        if (all(tsums < 0))
            stop("could not fit base learner in boosting iteration ", m)
        xselect[m] <- order(tsums)[1]
        basess <- ss[[xselect[m]]]

        ### update step
        fit <- fit + nu * fitted(basess)

        ### negative gradient vector, the new `residuals'
        u <- ngradient(y, fit, weights)

        ### evaluate risk, either for the learning sample (inbag)
        ### or the test sample (oobag)
        if (risk == "inbag") mrisk[m] <- riskfct(y, fit, weights)
        if (risk == "oobag") mrisk[m] <- riskfct(y, fit, oobweights)

        ### save the model, i.e., the selected coefficient and variance
        if (control$saveensss)
            ens[[m]] <- basess

        ## free memory
        rm("basess")

        ### print status information
        if (trace)
            do_trace(m, risk = mrisk, step = tracestep, width = mstop)
    }

    RET <- list(xselect = xselect,
                ensemble = ens,         ### list of baselearners
                fit = fit,              ### vector of fitted values
                offset = offset,        ### offset
                ustart = ustart,        ### first negative gradients
                risk = mrisk,           ### empirical risks for m = 1, ..., mstop
                control = control,      ### control parameters
                family = family,        ### family object
                response = y,           ### the response variable
                weights = weights       ### weights used for fitting
    )

    ### function for computing hat matrices of individual predictors
    RET$hat <- function(j) hatvalues(bl[[j]])

    RET$predict <- function(newdata = NULL, which = NULL, components = FALSE, Sum = TRUE) {

        if (is.null(which))
            which <- sort(unique(xselect))

        if (length(which) == 1)
            return(bl[[which]]$predict(ens[xselect == which], newdata = newdata, Sum = Sum))

        pr <- sapply(which, function(w) 
            bl[[w]]$predict(ens[xselect == w], newdata = newdata))
        colnames(pr) <- names(bl)[which]
        if (components) return(nu * pr)
        offset + nu * rowSums(pr)
    }

    class(RET) <- "mboost"
    return(RET)
}


x <- gl(50, 1009) ###rpois(10, lambda = 10)
x[sample(1:length(x), 100)] <- NA
y <- rnorm(length(x))
x[sample(1:length(x), 100)] <- NA
w <- rpois(length(x), lambda = 1)

system.time(c1 <- bl(x)$dpp(w)$fit(y)$model)
system.time(c2 <- coef(lm(y ~ x, weights = w)))
max(abs(c1 - c2))

set.seed(29)
x <- rpois(100, lambda = 7)
y <- rnorm(length(x))
w <- rpois(length(x), lambda = 1)
system.time(a1 <- bl(x, hfun = hyper_bbs, Xfun = mm_bbs)$dpp(w)$fit(y)$model)
system.time(a2 <- as.vector(attr(mboost:::bbs(x), "dpp")(w)$fit(y)$model))
max(abs(a1 - a2))

b1 <- bl(x, hfun = hyper_bbs, Xfun = mm_bbs)$dpp(w)
b1$predict(list(b1$fit(y), b1$fit(y + 2)), newdata = NULL, Sum = FALSE)
fitted(b1$fit(y))

a1 <- bl(x, y, hfun = hyper_bbs, Xfun = mm_bbs)$dpp(w)$fit(y)$model
a2 <- as.vector(attr(mboost:::bspatial(x, y, df = 4), "dpp")(w)$fit(y)$model)
max(abs(a1 - a2))

data("bodyfat", package = "mboost")

attach(bodyfat)
b <- list(blage = bl(age, hfun = hyper_bbs, Xfun = mm_bbs), 
          blhih = bl(hipcirc, hfun = hyper_bbs, Xfun = mm_bbs))
a <- mboost(b, DEXfat)

cc <- gamboost(DEXfat ~ age + hipcirc)

max(abs(a$predict() - cc$fit))

plot(model.frame(b[[1]]), a$predict(which = 1))

matplot(a$predict(which = 1, Sum = FALSE))
