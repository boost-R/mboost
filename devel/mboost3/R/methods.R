
### compute predictions
### <FIXME>: add arguments (for documentation)
###          add link argument (family needs to be touched)
### </FIXME>
predict.mboost <- function(object, type = c("lp", "response"), ...) {
    pr <- object$predict(...)
    type <- match.arg(type)
    if (is.factor(y <- object$response) && type == "response")
        return(factor(levels(y)[(pr > 0) + 1], levels = levels(y)))
    return(pr)
}

### extract coefficients
coef.mboost <- function(object, ...)
    object$coef(...)

### compute boosting hat matrix and its trace
hatvalues.mboost <- function(model, ...) {
    H <- model$hatvalues(...)
    n <- length(model$response)
    if (checkL2(model)) {
        op <- .Call("R_trace_gamboost", as.integer(n), H,
                    as.integer(model$xselect()), PACKAGE = "mboost")
    } else {
        fitm <- predict(model, aggregate = "cumsum")
        op <- bhatmat(n, H, model$xselect(), fitm, model$family@fW)
    }
    RET <- diag(op[[1]])  
    attr(RET, "hatmatrix") <- op[[1]]
    attr(RET, "trace") <- op[[2]]
    RET
}

AIC.mboost <- function(object, method = c("corrected", "classical", "gMDL"),
                       df = c("trace", "actset"), ..., k = 2) {

    df <- match.arg(df)
    if (df == "trace") {
        hatval <- hatvalues(object)
        RET <- AICboost3(object, method = method,
                         df = attr(hatval, "trace"), k = k)
    } 
    if (df == "actset") {
        ### compute active set: number of non-zero coefficients
        ### for each boosting iteration
        xs <- object$xselect(cw = TRUE)
        xu <- sort(sapply(unique(xs), function(i) which(xs == i)[1]))
        xu <- c(xu, mstop(object) + 1)
        df <- rep(1:(length(xu) - 1), diff(xu))
        ### <FIXME>: offset = 0 may mean hat(offset) = 0 or
        ### no offset computed at all!
        if (object$offset != 0) df <- df + 1
        ### </FIXME>
        RET <- AICboost3(object, method = method, 
                         df = df, k = k)
    }
    return(RET)
}

AICboost3 <- function(object, method = c("corrected", "classical", "gMDL"), df, k = 2) {

    if (object$control$risk != "inbag")
        return(NA)
    method <- match.arg(method)

    if (checkL2(object) && method == "classical")
        stop("classical AIC method not implemented for Gaussian family")
    if (!checkL2(object) && method == "corrected")
        stop("corrected AIC method not implemented for non-Gaussian family")

    sumw <- sum(model.weights(object))
    if (method == "corrected") 
        AIC <- log(object$risk() / sumw) +
               (1 + df/sumw) / (1 - (df + 2)/sumw)

    ### loss-function is to be MINIMIZED, take -2 * logLik == 2 * risk
    if (method == "classical")
        AIC <- 2 * object$risk() + k * df
    if (method == "gMDL"){
        s <- object$risk()/(sumw - df)
        AIC <- log(s) + df/sumw * log((sum(object$response^2) - object$risk())
                      /(df * s))
        }
    mstop <- which.min(AIC)
    RET <- AIC[mstop]

    attr(RET, "mstop") <- which.min(AIC)
    attr(RET, "df") <- df  
    attr(RET, "AIC") <- AIC
    attr(RET, "corrected") <- method == "corrected"

    class(RET) <- "gbAIC"
    return(RET)
}


### compute fitted values
fitted.mboost <- function(object, ...) {
    args <- list(...)
    if (length(args) == 0)
        return(object$fitted())
    object$predict(...)
}

### residuals (the current negative gradient)
resid.mboost <- function(object, ...) 
    object$resid()

logLik.mboost <- function(object, ...)
    object$logLik()

### restrict or enhance models to less/more
### boosting iterations.
### ATTENTION: x gets CHANGED!
"[.mboost" <- function(x, i, return = TRUE, ...) {
    stopifnot(length(i) == 1 && i > 0)
    x$subset(i)
    if (return) return(x)
    return(NULL)
}

mstop.mboost <- function(object, ...) object$mstop()

update.mboost <- function(object, weights, ...)
    object$update(weights)

model.frame.mboost <- function(formula, ...)
    formula$model.frame(...)

response.mboost <- function(object, ...)
    object$response

predict.Glmboost <- function(object, newdata = NULL, ...) {

    if (!is.null(newdata)) {
        newdata <- object$newX(newdata)
    }
    object$predict(newdata = newdata, ...)
}

coef.Glmboost <- function(object, ...) {
    cf <- object$coef(...)
    off <- attr(cf, "offset")
    cf <- cf[[1]]
    names(cf) <- variable.names(object)
    attr(cf, "offset") <- off
    cf
}

hatvalues.Glmboost <- function(model, ...) {

    if (!checkL2(model)) return(hatvalues.mboost(model))
    Xf <- t(model$basemodel[[1]]$MPinv()) * model$control$nu
    X <- model$baselearner[[1]]$get_data()
    op <- .Call("R_trace_glmboost", X, Xf,
                as.integer(model$xselect(cw = TRUE)),
                PACKAGE = "mboost")
    RET <- diag(op[[1]])
    attr(RET, "hatmatrix") <- op[[1]]  
    attr(RET, "trace") <- op[[2]] 
    RET
}

### methods: print
print.mboost <- function(x, ...) {

    cat("\n")
    cat("\t Model-based Boosting\n")
    cat("\n")
    if (!is.null(x$call))
    cat("Call:\n", deparse(x$call), "\n\n", sep = "")
    show(x$family)
    cat("\n")
    cat("Number of boosting iterations: mstop =", mstop(x), "\n")
    cat("Step size: ", x$control$nu, "\n")
    cat("Offset: ", x$offset, "\n")
    cat("Baselearner(s): \n")
    print(names(variable.names(x)))
    cat("\n")
    invisible(x)
}

### methods: print
print.Glmboost <- function(x, ...) {

    cat("\n")
    cat("\t Generalized Linear Models Fitted via Gradient Boosting\n")
    cat("\n")
    if (!is.null(x$call))
    cat("Call:\n", deparse(x$call), "\n\n", sep = "")
    show(x$family)
    cat("\n")
    cat("Number of boosting iterations: mstop =", mstop(x), "\n")
    cat("Step size: ", x$control$nu, "\n")
    cat("Offset: ", x$offset, "\n")
    cat("\n")
    cat("Coefficients: \n")
    cf <- coef(x)
    attr(x, "offset") <- NULL
    print(cf)
    cat("\n")
    invisible(x)
}

variable.names.mboost <- function(object, ...) {
    ret <- sapply(object$baselearner, function(x) x$get_names())
    if (is.matrix(ret)) ret <- ret[, , drop = TRUE]
    ret
}
