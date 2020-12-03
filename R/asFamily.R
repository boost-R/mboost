
as.Family <- function(object, ...)
    UseMethod("as.Family")

.copy_variables <- function(call, envir) {
    nm <- names(call)
    nm <- nm[nm != ""]
    vars <- do.call("c", sapply(call[nm], all.vars))
    ret <- new.env()
    for(n in vars) {
        if (n %in% names(envir))
            assign(n, get(n, envir), ret)
    }
    ret
}

as.Family.lm <- function(object) {

    cf <- coef(object)
    X <- model.matrix(object)
    Y <- model.response(model.frame(object))
    N <- nrow(model.frame(object))

    risk <- function(y, f, w = 1) {
        if (length(f) == 1) f <- rep(f, N)
        if (length(w) == 1) w <- rep(w, N)

        sum(lm.wfit(x = X, y = Y, w = w, offset = c(f))$residuals^2)
    }

    ngradient <- function(y, f, w = 1) {
        if (length(f) == 1) f <- rep(f, N)
        if (length(w) == 1) w <- rep(w, N)

        ### this is not sufficient, data used to fit "object" needs to be 
        ### the same as mboost(..., data); na.action can and will mess this
        ### up!!!
        stopifnot(NROW(f) == N)

        model <<- lm.wfit(x = X, y = Y, w = w, offset = c(f))
        cf <<- coef(model)

        ### residual wrt a constant for _all_ observations
        model$residuals
    }

    Family(ngradient = ngradient, risk = risk,
           check_y = function(y) rep(TRUE, N),
           offset = function(y, w) 0,
           nuisance = function() return(list(coefficients = cf)),
           name = "lm",
           response = function(f) f)
}

as.Family.glm <- function(object) {

    if (!require("sandwich"))
        stop("Package sandwich is required for this family")

    model <- object
    cf <- coef(object)
    X <- model.matrix(object)
    stopifnot("(Intercept)" %in% colnames(X))
    Y <- model.response(model.frame(object))
    N <- nrow(model.frame(object))

    fam <- family(object)$family
    p <- object$rank
    if (fam %in% c("gaussian", "Gamma", "inverse.gaussian")) 
        p <- p + 1

    risk <- function(y, f, w = 1) {
        if (length(f) == 1) f <- rep(f, N)
        if (length(w) == 1) w <- rep(w, N)

        tmp <- glm.fit(x = X, y = Y, w = w, offset = c(f), family = object$family)
        return(-(p - tmp$aic / 2))
    }

    ngradient <- function(y, f, w = 1) {
        if (length(f) == 1) f <- rep(f, N)
        if (length(w) == 1) w <- rep(w, N)

        ### this is not sufficient, data used to fit "object" needs to be 
        ### the same as mboost(..., data); na.action can and will mess this
        ### up!!!
        stopifnot(NROW(f) == N)

        tmp <- glm.fit(x = X, y = Y, w = w, offset = c(f), family = object$family)
        if ((p - tmp$aic / 2) < logLik(model))
            warning("risk increase; decrease stepsize nu")

        model[names(tmp)] <<- tmp
        cf <<- coef(model)

        ### residual wrt a constant for _all_ observations
        sandwich::estfun(model)[,"(Intercept)"]
    }

    Family(ngradient = ngradient, risk = risk,
           check_y = function(y) rep(TRUE, N),
           offset = function(y, w) 0,
           nuisance = function() return(list(coefficients = cf)),
           name = "glm",
           response = function(f) object$family$linkinv(f))
}

model.matrix.merFamily <- function(object, ...)
    matrix(1, nrow = N, ncol = 1)

as.Family.merMod <- function(object) {

    if (!require("lme4"))
        stop("package lme4 is required for this family")
    if (!require("sandwich"))
        stop("Package sandwich is required for this family")

    model <- object
    rf <- lme4::ranef(model)
    ff <- lme4::fixef(model)
    cl <- object@call
    ffm <- formula(object)
    env <- .copy_variables(cl, environment(ffm))
    environment(model@call$formula) <- env
    environment(attr(model@frame, "formula")) <- env

    N <- nrow(model.frame(object))

    ffm <- as.formula(paste(deparse(ffm[[2]]), "~ 1 - 1"))
    cl$formula <- ffm
    fun <- gsub("^lme4::", "", deparse(cl[[1L]]))
    stopifnot(fun %in% c("lmer", "glmer"))
    fun <- gsub("er$", "", fun)
    cl[[1L]] <- as.name(fun)
    environment(cl$formula) <- env
    m0 <- eval(cl, envir = env)

    risk <- function(y, f, w = 1) {
        if (length(f) == 1) f <- rep(f, nrow(model.frame(model)))
        if (length(w) == 1) w <- rep(w, nrow(model.frame(model)))

        assign("offset_.", f, env)
        assign("weights_.", w, env)

        tmp <- eval(update(model, offset = offset_., weights = weights_.), env)
        lp <- predict(tmp)
        assign("lp_.", lp, env)

        -logLik(eval(update(m0, offset = lp_., weights = weights_., evaluate = FALSE), env))
    }

    ngradient <- function(y, f, w = 1) {
        if (length(f) == 1) f <- rep(f, N)
        if (length(w) == 1) w <- rep(w, N)

        ### this is not sufficient, data used to fit "object" needs to be 
        ### the same as mboost(..., data); na.action can and will mess this
        ### up!!!
        stopifnot(NROW(f) == N)

        assign("offset_.", f, env)
        assign("weights_.", w, env)
        assign("start_.", getME(model, c("theta", "fixef")), env)

        model <<- eval(update(model, offset = offset_., weights = weights_., start = start_.), env)
        rf <<- lme4::ranef(model)
        ff <<- lme4::fixef(model)
        lp <- predict(model)

        assign("lp_.", lp, env)

        tmp <- eval(update(m0, offset = lp_., weights = weights_., evaluate = FALSE), env)
        if (logLik(tmp) < logLik(m0))
            warning("risk increase; decrease stepsize nu")
        m0 <<- tmp

        ### note: inherits(m0, "lm") will also be true for glms
        if (class(m0)[1] == "lm")
            return(m0$residuals)

        ### for glms fake an intercept = 0 model
        class(tmp) <- c("merFamily", class(tmp))
        tmp$coefficients <- 0

        ### residual wrt a constant for _all_ observations
        sandwich::estfun(tmp)
    }

    Family(ngradient = ngradient, risk = risk,
           check_y = function(y) rep(TRUE, N),
           offset = function(y, w) 0,
           nuisance = function() return(list(ff = ff, rf = rf)),
           name = "glmer",
           response = function(f) object@resp$family$linkinv(f))
}

as.Family.coxph <- function(object) {

    model <- object
    cf <- coef(model)
    ffm <- formula(object)
    cl <- object$call
    env <- .copy_variables(cl, environment(ffm))
    environment(model$formula) <- env
    N <- nrow(model.frame(object))

    fm <- . ~ . + offset(offset_.)
    m0 <- update(model, formula = fm, weights = weights_., evaluate = FALSE)
    environment(m0$formula) <- env

    risk <- function(y, f, w = 1) {
        if (length(f) == 1) f <- rep(f, nrow(model.frame(model)))
        if (length(w) == 1) w <- rep(w, nrow(model.frame(model)))

        assign("offset_.", f, env)
        assign("weights_.", w, env)

        tmp <- eval(m0, env)
        -logLik(tmp)
    }

    ngradient <- function(y, f, w = 1) {
        if (length(f) == 1) f <- rep(f, N)
        if (length(w) == 1) w <- rep(w, N)

        ### this is not sufficient, data used to fit "object" needs to be 
        ### the same as mboost(..., data); na.action can and will mess this
        ### up!!!
        stopifnot(NROW(f) == N)

        assign("offset_.", f, env)
        assign("weights_.", w, env)

        tmp <<- eval(m0, env)
        if (logLik(tmp) < logLik(model))
            warning("risk increase; decrease stepsize nu")
        model <<- tmp
        cf <<- coef(model)

        ### residual wrt a constant for _all_ observations
        residuals(model, type = "martingale")
    }

    Family(ngradient = ngradient, risk = risk,
           check_y = function(y) rep(TRUE, N),
           offset = function(y, w) 0,
           nuisance = function() return(list(coefficients = cf)),
           name = "coxph",
           ### conditional survivor function???
           response = function(f) NA)
}
