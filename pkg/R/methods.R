
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
                    as.integer(model$xselect(cw = TRUE)), PACKAGE = "mboost")
    } else {
        fitm <- predict(model, aggregate = "cumsum")
        op <- bhatmat(n, H, model$xselect(cw = TRUE), fitm, model$family@fW)
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
        RET <- AICboost(object, method = method,
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
        RET <- AICboost(object, method = method, 
                         df = df, k = k)
    }
    return(RET)
}

AICboost <- function(object, method = c("corrected", "classical", "gMDL"), df, k = 2) {

    if (object$control$risk != "inbag")
        return(NA)
    method <- match.arg(method)

    if (checkL2(object) && method == "classical")
        stop("classical AIC method not implemented for Gaussian family")
    if (!checkL2(object) && method == "corrected")
        stop("corrected AIC method not implemented for non-Gaussian family")

    sumw <- sum(model.weights(object)[!is.na(fitted(object))])
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

predict.glmboost <- function(object, newdata = NULL, type = c("lp", "response"), ...) {

    if (!is.null(newdata)) {
        newdata <- object$newX(newdata)
    }
    pr <- object$predict(newdata = newdata, ...)
    type <- match.arg(type)
    if (is.factor(y <- object$response) && type == "response")
        return(factor(levels(y)[(pr > 0) + 1], levels = levels(y)))
    return(pr)
}

coef.glmboost <- function(object, ...) {
    cf <- object$coef(...)
    off <- attr(cf, "offset")
    cf <- cf[[1]]
    names(cf) <- variable.names(object)
    attr(cf, "offset") <- off
    cf
}

hatvalues.glmboost <- function(model, ...) {

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
    cat("\n")
    cat("Baselearner(s): \n")
    print(names(variable.names(x)))
    cat("\n")
    invisible(x)
}

### methods: print
print.glmboost <- function(x, ...) {

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
    ### <FIXME> is pasting what we want?
    ret <- sapply(object$baselearner, function(x) 
                  paste(x$get_names(), collapse = ", "))
    ### </FIXME>
    if (is.matrix(ret)) ret <- ret[, , drop = TRUE]
    ret
}

variable.names.glmboost <- function(object, ...) {
    ret <- object$baselearner[[1]]$get_names()
    names(ret) <- ret
    ret
}


selected <- function(object) UseMethod("selected", object)

selected.mboost <- function(object)
    object$xselect()

selected.glmboost <- function(object)
   object$xselect(cw = TRUE)

summary.mboost <- function(object, ...) {

    ret <- list(object = object, selprob = NULL)
    xs <- selected(object)
    nm <- variable.names(object)
    selprob <- tabulate(xs, nbin = length(nm)) / length(xs)
    names(selprob) <- names(nm)
    selprob <- sort(selprob, decreasing = TRUE)
    ret$selprob <- selprob[selprob > 0]
    class(ret) <- "summary.mboost"
    return(ret)
}

print.summary.mboost <- function(x, ...) {

    print(x$object)
    cat("Selection frequencies:\n")
    print(x$selprob)
    cat("\n")
}
