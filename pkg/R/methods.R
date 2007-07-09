
### 
### Methods for gradient boosting objects
### 


### methods: subset
"[.gb" <- function(x, i, ...) {
    mstop <- mstop(x)
    if (i == mstop) return(x)
    if (length(i) != 1)
        stop("not a positive integer")
    if (i < 1 || i > mstop)
        warning("invalid number of boosting iterations")
    indx <- 1:min(max(i, 1), mstop)
    x$ensemble <- x$ensemble[indx, , drop = FALSE]
    x$risk <- x$risk[indx]
    x$fit <- x$predict(mstop = max(indx))
    x
}

### methods: prediction
predict.gb <- function(object, newdata = NULL, type = c("lp", "response"), ...) {
    type <- match.arg(type)
    y <- object$data$y
    lp <- object$predict(newdata = newdata, mstop = mstop(object), ...)
    if (type == "response" && is.factor(y))
       return(factor(levels(y)[(lp > 0) + 1], levels = levels(y)))
    return(lp)
}

### methods: fitted
fitted.gb <- function(object, type = c("lp", "response"), ...) {
    type <- match.arg(type)
    lp <- object$fit
    y <- object$data$y
    if (type == "response" && is.factor(y))
       return(factor(levels(y)[(lp > 0) + 1], levels = levels(y)))
    return(lp)
}

fitted.blackboost <- function(object, type = c("lp", "response"), ...) {
    type <- match.arg(type)
    lp <- object$fit
    y <- object$data@responses@variables[[1]]
    if (type == "response" && is.factor(y))
       return(factor(levels(y)[(lp > 0) + 1], levels = levels(y)))
    return(lp)
}

### methods: resid
resid.gb <- function(object, ...)
    object$family@ngradient(object$yfit, fitted(object))

resid.blackboost <- resid.gb 

### methods: hatvalues, either exact (L2) or approximately
hatvalues.gb <- function(model, ...) {

    n <- nrow(model$data$x)   
    p <- ncol(model$data$x)
    ens <- model$ensemble
    nu <- model$control$nu

    ### list of hat matices
    H <- vector(mode = "list", length = p)

    for (xs in unique(ens[,"xselect"]))
        ### compute hat matrix for this covariate
        if (is.null(H[[xs]]))
            H[[xs]] <- nu * model$hat(xs)

    if (checkL2(model)) {
        op <- .Call("R_trace_gamboost", as.integer(n), H, 
                    as.integer(ens[,"xselect"]), PACKAGE = "mboost")
    } else {
        g <- gm(model)
        fitm <- t(apply(g, 1, function(a) cumsum(a)))
        op <- bhatmat(n, H, ens[,"xselect"], fitm, model$family@fW)
    }
    RET <- diag(op[[1]])
    attr(RET, "hatmatrix") <- op[[1]]
    attr(RET, "trace") <- op[[2]]
    RET
}


### methods: AIC
AIC.gb <- function(object, method = c("corrected", "classical"), ...) {

    if (object$control$risk != "inbag")
        return(NA)
    method <- match.arg(method)

    if (checkL2(object) && method == "classical")
        stop("classical AIC method not implemented for Gaussian family")
    if (!checkL2(object) && method == "corrected")
        stop("corrected AIC method not implemented for non-Gaussian family")

    hatval <- hatvalues(object)
    hatmat <- attr(hatval, "hatmatrix")
    trace <- attr(hatval, "trace")

    sumw <- sum(object$weights)
    if (method == "corrected")
        AIC <- log(object$risk / sumw) + 
               (1 + trace/sumw) / (1 - (trace + 2)/sumw)

    ### loss-function is to be MINIMIZED, take -2 * logLik == 2 * risk
    if (method == "classical")
        AIC <- 2 * object$risk + 2 * trace

    mstop <- which.min(AIC)
    RET <- AIC[mstop]

    attr(RET, "hatmatrix") <- hatmat
    attr(RET, "mstop") <- which.min(AIC)
    attr(RET, "df") <- trace
    attr(RET, "AIC") <- AIC
    attr(RET, "corrected") <- method == "corrected"

    class(RET) <- "gbAIC"
    return(RET)
}

logLik.gb <- function(object, ...)
    -object$family@risk(object$data$yfit, fitted(object), object$weights)

print.gbAIC <- function(x, ...) {
    mstop <- mstop(x)
    df <- attr(x, "df")[mstop]
    attributes(x) <- NULL
    print(x)
    cat("Optimal number of boosting iterations:", mstop, "\n")
    cat("Degrees of freedom", paste("(for mstop = ", mstop, "):", sep = "", 
                                    collapse = ""), df, "\n")
    invisible(x)
}

plot.gbAIC <- function(x, y = NULL, ...) {
    mstop <- mstop(x)
    class(x) <- NULL
    plot(attr(x, "AIC"), xlab = "Number of boosting iterations", 
         ylab = ifelse(attr(x, "corrected"), "Corrected AIC", "AIC"),
         type = "l", ...)
    points(mstop, x)
    ylim <- list(...)$ylim
    if (!is.null(ylim)) {
        ymin <- ylim[1] * ifelse(ylim[1] < 0, 2, 0.5)
    } else {
        ymin <- x - x/2
    }
    lines(c(mstop, mstop), 
          c(ymin, x), lty = 2)
}

mstop <- function(object, ...) UseMethod("mstop")

mstop.gbAIC <- function(object, ...) attr(object, "mstop")

mstop.gb <- function(object, ...) nrow(object$ensemble)

mstop.blackboost <- function(object, ...) length(object$ensemble)

### methods: variance of predictions (!)
vcov.gb <- function(object, ...) {
    aic <- AIC(object)
    Bhat <- attr(aic, "hatmatrix")
    df <- attr(aic, "df")[length(attr(aic, "df"))]
    sigma2 <- sum(resid(object)^2) / (length(object$yfit) - df)
    return(sigma2 * colSums(Bhat^2))
}
