### compute predictions
### <FIXME>:
###          add link argument (family needs to be touched)
### </FIXME>
.predictmboost <- function(y, pr, type, nm, family) {
    if (!isMATRIX(pr)) {
        pr <- as.vector(pr)
        names(pr) <- nm
    }
    if (is.list(pr)) {
        rownames(pr) <- nm
        if (type != "link")
            warning("argument", sQuote("type"), " is ignored")
        return(pr)
    }
    if (type == "link") ret <- pr
    if (type == "response") ret <- family@response(pr)
    if (type == "class") {
        if (!is.factor(y)) stop("response is not a factor")
            ret <- factor(levels(y)[family@rclass(pr)], levels = levels(y))
        }
    names(ret) <- names(pr)
    return(ret)
}


predict.mboost <- function(object, newdata = NULL,
    type = c("link", "response", "class"), which = NULL,
    aggregate = c("sum", "cumsum", "none"), ...) {

    type <- match.arg(type)
    aggregate <- match.arg(aggregate)
    if (aggregate != "sum")
        stopifnot(type == "link")
    pr <- object$predict(newdata = newdata,
                         which = which, aggregate = aggregate)
    nm <- rownames(newdata)
    if (is.null(newdata)) nm <- object$rownames
    if (is.list(pr)){
        RET <- lapply(pr, .predictmboost, y = object$response,
                      type = type, nm = nm, family = object$family)
        if (type == "link") attr(RET, "offset") <- attr(pr, "offset")
        return(RET)
    }
    RET <- .predictmboost(object$response, pr, type = type, nm = nm,
                          family = object$family)
    if (type != "link") attr(RET, "offset") <- NULL
    return(RET)
}

### extract coefficients
coef.mboost <- function(object, which = NULL,
    aggregate = c("sum", "cumsum", "none"), ...) {

    args <- list(...)
    if (length(args) > 0)
        warning("Arguments ", paste(names(args), sep = ", "), " unknown")
    if (grepl("Negative Binomial Likelihood", object$family@name))
        message("\nNOTE: Coefficients from a Binomial model are half the size of ",
                "coefficients\n from a model fitted via ",
                "glm(... , family = 'binomial').\n",
                "See Warning section in ?coef.mboost\n")
    object$coef(which = which, aggregate = aggregate)
}

### compute boosting hat matrix and its trace
hatvalues.gamboost <- function(model, ...) {
    H <- model$hatvalues(...)
    n <- length(model$response)

    ### <FIXME> better checks
    L2 <- FALSE
    if (!extends(class(model$family), "boost_family_glm")) {
        warning("AIC might not be reasonable for family ", model$family@name)
        L2 <- TRUE
    } else {
        fW <- model$family@fW
    }
    ### </FIXME>
    if (checkL2(model) || L2) {
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

hatvalues.mboost <- function(model, ...) {
    stop("hatvalues are not implemented for models fitted via function ",
         sQuote("mboost"), ".\n Please use functions ", sQuote("gamboost"),
         " or ", sQuote("glmboost"), " instead.")
}

AIC.mboost <- function(object, method = c("corrected", "classical", "gMDL"),
                       df = c("trace", "actset"), ..., k = 2) {

    df <- match.arg(df)
    if (df == "actset" && !inherits(object, "glmboost")) {
        df <- "trace"
        warning("df = ", dQuote("actset"), " can only be used with ",
                sQuote("glmboost"), " models. df = ", dQuote("trace"),
                " is used instead.")
    }
    if (df == "trace") {
        hatval <- hatvalues(object)
        RET <- AICboost(object, method = method,
                        df = attr(hatval, "trace"), k = k)
    }
    if (df == "actset") {
        ### compute active set: number of non-zero coefficients
        ### for each boosting iteration
        xs <- object$xselect()
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
    if (length(args) == 0) {
        ret <- object$fitted()
        names(ret) <- object$rownames
    } else {
        ret <- predict(object, newdata=NULL, ...)
        #if (NROW(ret) == length(ret))
        #    rownames(ret) <- object$rownames
    }
    ret
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
    invisible(NULL)
}

"mstop<-" <- function(x, value) {
    return(x[value, return = TRUE])
}

mstop.mboost <- function(object, ...) object$mstop()

update.mboost <- function(object, weights, oobweights = NULL,
                          risk = NULL, trace = NULL, ...)
    object$update(weights = weights, oobweights = oobweights,
                  risk = risk, trace = trace)

model.frame.mboost <- function(formula, ...)
    formula$model.frame(...)

response.mboost <- function(object, ...)
    object$response

predict.glmboost <- function(object, newdata = NULL,
    type = c("link", "response", "class"), which = NULL,
    aggregate = c("sum", "cumsum", "none"), ...) {

    args <- list(...)
    if (length(args) > 0)
        warning("Arguments ", paste(names(args), sep = ", "), " unknown")

    aggregate <- match.arg(aggregate)

    ### X <- object$baselearner[[1]]$get_data()
    if (!is.null(newdata))
        newdata <- object$newX(newdata)

    ### <FIXME> implement predictions here </FIXME>
    ### which <- object$which(which, usedonly = nw <- is.null(which))
    ### cf <- coef(object, which = which, aggregate = aggregate)
    ### pr <- X[, which, drop = FALSE] %*% cf

    pr <- object$predict(newdata = newdata, which = which,
                         aggregate = aggregate)
    type <- match.arg(type)
    nm <- rownames(newdata)
    if (is.null(newdata)) nm <- object$rownames
    if (is.list(pr))
        return(lapply(pr, .predictmboost, y = object$response,
                      type = type, nm = nm, family = object$family))
    .predictmboost(object$response, pr, type = type, nm = nm,
                   family = object$family)
}

coef.glmboost <- function(object, which = NULL,
    aggregate = c("sum", "cumsum", "none"), off2int = FALSE, ...) {

    args <- list(...)
    if (length(args) > 0)
        warning("Arguments ", paste(names(args), sep = ", "), " unknown")

    if (grepl("Negative Binomial Likelihood", object$family@name))
        message("\nNOTE: Coefficients from a Binomial model are half the size of ",
                "coefficients\n from a model fitted via ",
                "glm(... , family = 'binomial').\n",
                "See Warning section in ?coef.mboost\n")

    aggregate <- match.arg(aggregate)
    cf <- object$coef(which = which, aggregate = aggregate)
    offset <- attr(cf, "offset")

    ### intercept = hat(beta[1]) - bar(x) %*% hat(beta[-1])
    assign <- object$assign
    cm <- object$center
    INTERCEPT <- sum(assign == 0) == 1
    if (!INTERCEPT && off2int)
        warning(sQuote("object"), " does not contain an intercept, ",
                sQuote("off2int = TRUE"), " ignored.")
    if (INTERCEPT && any(cm != 0)) {
        intercept <- which(assign == 0)
        if (is.null(which)) {
            which <- object$which(usedonly = TRUE)
            if (!intercept %in% which) {
                which <- c(intercept, which)
                cf <- object$coef(which = which, aggregate = aggregate)
            }
        }
        which <- object$which(which)
        if (intercept %in% which) {
            if (all(which %in% object$which(usedonly = TRUE))) {
                cm <- cm[which]
                scf <- sapply(1:length(cf), function(i) cf[[i]] * cm[i])
                if (!is.matrix(scf)) scf <- matrix(scf, nrow = 1)
                cf[[intercept]] <- cf[[intercept]] - rowSums(scf)
            } else {
                cm <- cm[object$which(usedonly = TRUE)]
                tmp <- object$coef(which = object$which(usedonly = TRUE),
                                   aggregate = aggregate)
                scf <- sapply(1:length(tmp), function(i) tmp[[i]] * cm[i])
                if (!is.matrix(scf)) scf <- matrix(scf, nrow = 1)
                cf[[intercept]] <- cf[[intercept]] - rowSums(scf)
            }
        }
    }

    if (aggregate == "sum") cf <- unlist(cf)
    if (aggregate == "none") {
        attr(cf, "offset") <- offset
        if (off2int)
            warning(sQuote("off2int = TRUE"), " ignored for ",
                    sQuote("aggregate = \"none\""))
    } else {
        if (off2int & length(offset) == 1) {
            cf[[1]] <- cf[[1]] + offset
        } else {
            if (off2int)
                warning(sQuote("off2int = TRUE"),
                        " ignored for non-scalar offset")
            attr(cf, "offset") <- offset
        }
    }
    cf
}

hatvalues.glmboost <- function(model, ...) {

    args <- list(...)
    if (length(args) > 0)
        warning("Arguments ", paste(names(args), sep = ", "), " unknown")

    if (!checkL2(model)) return(hatvalues.gamboost(model))
    Xf <- t(model$basemodel[[1]]$MPinv()) * model$control$nu
    X <- model$baselearner[[1]]$get_data()
    op <- .Call("R_trace_glmboost", as(X, "matrix"), as(Xf, "matrix"),
                as.integer(model$xselect()),
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
    if (!is.null(x$call)) {
        if(length(deparse(x$call$data)) > 20)
            x$call$data <- deparse(x$call$data, nlines = 1)
        cat("Call:\n", deparse(x$call), "\n\n", sep = "")
    }
    show(x$family)
    cat("\n")
    cat("Number of boosting iterations: mstop =", mstop(x), "\n")
    cat("Step size: ", x$control$nu, "\n")
    cat("Offset: ", x$offset, "\n")
    cat("Number of baselearners: ", length(variable.names(x)), "\n")
    cat("\n")
    invisible(x)
}

### methods: print
print.glmboost <- function(x, ...) {

    cat("\n")
    cat("\t Generalized Linear Models Fitted via Gradient Boosting\n")
    cat("\n")
    if (!is.null(x$call)) {
        if(length(deparse(x$call$data)) > 20)
            x$call$data <- deparse(x$call$data, nlines = 1)
        cat("Call:\n", deparse(x$call), "\n\n", sep = "")
    }
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

variable.names.mboost <- function(object, which = NULL, usedonly = FALSE, ...) {

    which <- object$which(which, usedonly = usedonly)

    args <- list(...)
    if (length(args) > 0)
        warning("Arguments ", paste(names(args), sep = ", "), " unknown")

    ### <FIXME> is pasting what we want?
    ret <- sapply(object$baselearner, function(x)
                  paste(x$get_names(), collapse = ", "))
    ### </FIXME>
    if (is.matrix(ret)) ret <- ret[, , drop = TRUE]
    ret[which]
}

variable.names.glmboost <- function(object, which = NULL, usedonly = FALSE, ...) {

    if (usedonly) {
        which <- object$which(usedonly = TRUE)
        ## if center = TRUE for model fitting intercept is implicitly selected
        center <- get("center", envir = environment(object$newX))
        if (center){
            intercept <- which(object$assign == 0)
            INTERCEPT <- sum(object$assign == 0) == 1
            if (INTERCEPT && !intercept %in% which)
                which <- c(intercept, which)
        }
    } else {
        which <- object$which(which, usedonly = usedonly)
    }

    args <- list(...)
    if (length(args) > 0)
        warning("Arguments ", paste(names(args), sep = ", "), " unknown")

    ret <- object$baselearner[[1]]$get_names()
    names(ret) <- ret
    ret[which]
}

selected.mboost <- function(object, ...)
    object$xselect()

summary.mboost <- function(object, ...) {

    ret <- list(object = object, selprob = NULL)
    xs <- selected(object)
    nm <- variable.names(object)
    selprob <- tabulate(xs, nbins = length(nm)) / length(xs)
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

nuisance <- function(object)
    UseMethod("nuisance")

nuisance.mboost <- function(object)
    object$nuisance()[[mstop(object)]]



extract <- function(object, ...)
    UseMethod("extract")

extract.mboost <- function(object, what = c("design", "penalty", "lambda", "df",
                                   "coefficients", "residuals",
                                   "variable.names", "bnames", "offset",
                                   "nuisance", "weights", "index", "control"),
                           which = NULL, ...){
    what <- match.arg(what)
    which <- object$which(which, usedonly = is.null(which))
    if (what %in% c("design", "penalty", "lambda", "df", "index")){
        fun <- function(which)
            extract(object$basemodel[[which]], what = what, ...)
        ret <- lapply(which, fun)
        names(ret) <- extract(object, what = "bnames", which = which)
        return(ret)
    }
    switch(what,
           "coefficients" = return(coef(object, which = which)),
           "residuals" = return(residuals(object)),
           "variable.names" = return(variable.names(object, which)),
           "bnames" = return(get("bnames", envir = environment(object$update))[which]),
           "offset" = return(object$offset),
           "nuisance" = return(nuisance(object)),
           "weights" = return(model.weights(object)),
           "control" = return(object$control))
}

extract.glmboost <- function(object, what = c("design", "coefficients", "residuals",
                                     "variable.names", "bnames", "offset",
                                     "nuisance", "weights", "control"),
                             which = NULL, asmatrix = FALSE, ...){
    what <- match.arg(what)
    center <- get("center", envir = environment(object$newX))
    if (is.null(which)) {
        which <- object$which(usedonly = TRUE)
        ## if center = TRUE for model fitting intercept is implicitly selected
        if (center){
            intercept <- which(object$assign == 0)
            INTERCEPT <- sum(object$assign == 0) == 1
            if (INTERCEPT && !intercept %in% which)
                which <- c(intercept, which)
        }
    } else {
        which <- object$which(which)
    }

    if (what == "design"){
        mat <- object$baselearner[[1]]$get_data()[,which]
        if (asmatrix)
            mat <- as.matrix(mat)
        return(mat)
    }

    switch(what,
           "coefficients" = return(coef(object, which = which)),
           "residuals" = return(residuals(object)),
           "variable.names" = return(variable.names(object, which)),
           "bnames" = return(get("bnames", envir = environment(object$update))[which]),
           "offset" = return(object$offset),
           "nuisance" = return(nuisance(object)),
           "weights" = return(model.weights(object)),
           "control" = return(object$control))
}

extract.blackboost <- function(object, ...)
    stop("function not yet implemented")

extract.blg <- function(object, what = c("design", "penalty", "index"),
                        asmatrix = FALSE, expand = FALSE, ...){
    what <- match.arg(what)
    #object <- object$dpp(rep(1, NROW(object$model.frame())))
    # return(extract(object, what = what,
    #               asmatrix = asmatrix, expand = expand))
    if (what == "design")
        mat <- get("X", envir = environment(object$dpp))
    if (what == "penalty")
        mat <- get("K", envir = environment(object$dpp))
    if (what == "index")
        return(get("index", envir = environment(object$dpp)))
    ## only applicable for design and penalty matrix
    if (asmatrix){
        mat <- as.matrix(mat)
    }
    if (expand && !is.null(indx <- extract(object, what = "index"))){
        a <- attributes(mat)
        mat <- mat[indx,]
        a[c("dim", "dimnames")] <- attributes(mat)[c("dim", "dimnames")]
        attributes(mat) <- a
    }
    return(mat)
}

extract.bl_lin <- function(object, what = c("design", "penalty", "lambda", "df",
                                   "weights", "index"),
                           asmatrix = FALSE, expand = FALSE,  ...){
    what <- match.arg(what)
    if (what == "design")
        mat <- get("X", envir = environment(object$fit))
    if (what == "penalty")
        mat <- get("K", envir = environment(object$fit))
    if (what %in% c("lambda", "df"))
        return(object$df()[what])
    if (what %in% c("weights", "index"))
        return(get(what, envir = environment(object$fit)))
    ## only applicable for design and penalty matrix
    if (asmatrix){
        mat <- as.matrix(mat)
    }
    if (expand && !is.null(indx <- extract(object, what = "index"))){
        a <- attributes(mat)
        mat <- mat[indx,]
        a[c("dim", "dimnames")] <- attributes(mat)[c("dim", "dimnames")]
        attributes(mat) <- a
    }
    return(mat)
}

extract.bl_tree <- function(object, what = c("design", "penalty", "lambda", "df",
                                    "weights", "index"),
                            ...){
    what <- match.arg(what)
    if (what == "weights"){
        return(get("weights", envir = environment(object$fit)))
    } else {
        warning(paste("model matrix, penalty matrix, lambda and index",
                      "do not exist for tree base-learners"))
        invisible(NULL)
    }
}

residuals.mboost <- function(object, ...){
    if(object$family@name == "Squared Error (Regression)")
        return(object$resid())
    stop(sQuote("residuals()"), " only implemented for ",
         sQuote("family = Gaussian()"))
}

risk <- function(object, ...)
    UseMethod("risk")

risk.mboost <- function(object, ...) {
    object$risk()
}
