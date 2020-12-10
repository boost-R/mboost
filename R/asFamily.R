
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

as.Family.lm <- function(object, ...) {

    model <- object
    cf <- coef(object)
    X <- model.matrix(object)
    Y <- model.response(model.frame(object))
    N <- nrow(model.frame(object))

    risk <- function(y, f, w = 1)
        sum(w * (Y - (c(f) + X %*% cf))^2)

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

as.Family.glm <- function(object, ...) {

    model <- object
    cf <- coef(object)
    X <- model.matrix(object)
    stopifnot("(Intercept)" %in% colnames(X))
    Y <- model.response(model.frame(object))
    N <- nrow(model.frame(object))

    fam <- family(object)
    p <- object$rank
    if (fam$family %in% c("gaussian", "Gamma", "inverse.gaussian")) 
        p <- p + 1

    risk <- function(y, f, w = 1) {
        mu <- fam$linkin(c(f) + X %*% cf)
        dev <- sum(fam$dev.resids(object$y, mu, w))
        aic <- fam$aic(y = object$y, n = 1, mu = mu, wt = w)
        ### see stats:::glm.wfit and stats:::logLik.glm
        return(aic / 2)
    }

    ngradient <- function(y, f, w = 1) {
        if (length(f) == 1) f <- rep(f, N)
        if (length(w) == 1) w <- rep(w, N)

        ### this is not sufficient, data used to fit "object" needs to be 
        ### the same as mboost(..., data); na.action can and will mess this
        ### up!!!
        stopifnot(NROW(f) == N)

        tmp <- glm.fit(x = X, y = Y, weights = w, offset = c(f), 
                       family = fam)
        #### see stats:::glm.wfit and stats:::logLik.glm
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
           response = function(f) fam$linkinv(f))
}

model.matrix.merFamily <- function(object, ...)
    matrix(1, nrow = object$N, ncol = 1)

as.Family.merMod <- function(object, ...) {

    model <- object
    rf <- lme4::ranef(model)
    ff <- lme4::fixef(model)
    theta <- lme4::getME(model, "theta")
    cl <- object@call
    ffm <- formula(object)
    env <- .copy_variables(cl, environment(ffm))
    environment(model@call$formula) <- env
    environment(attr(model@frame, "formula")) <- env

    N <- nrow(mf <- model.frame(object))
    Y <- model.response(mf)

    ffm <- as.formula(paste(deparse(ffm[[2]]), "~ 1 - 1"))
    cl$formula <- ffm
    fun <- gsub("^lme4::", "", deparse(cl[[1L]]))
    stopifnot(fun %in% c("lmer", "glmer"))
    fun <- gsub("er$", "", fun)
    cl[[1L]] <- as.name(fun)
    environment(cl$formula) <- env
    m0 <- eval(cl, envir = env)

    offset_. <- weights_. <- start_. <- lp_. <- NA

    fam <- family(object)

    risk <- function(y, f, w = 1) {
        ### newdata removes offset from predict
        lp <- predict(model, newdata = mf) + c(f)
        if (inherits(object, "lmerMod"))
            return(sum(w * (Y - lp)^2))
        mu <- fam$linkinv(lp)
        dev <- sum(fam$dev.resids(m0$y, mu, w))
        aic <- fam$aic(y = m0$y, n = 1, mu = mu, wt = w, dev = dev)
        ### see as.Family.glm@risk
        return(aic / 2)
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
        assign("start_.", lme4::getME(model, c("theta", "fixef")), env)

        model <<- eval(update(model, offset = offset_., weights = weights_., start = start_.), env)
        rf <<- lme4::ranef(model)
        ff <<- lme4::fixef(model)
        theta <<- lme4::getME(model, "theta")
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
        tmp$N <- N
        tmp$coefficients <- 0

        ### residual wrt a constant for _all_ observations
        sandwich::estfun(tmp)
    }

    Family(ngradient = ngradient, risk = risk,
           check_y = function(y) rep(TRUE, N),
           offset = function(y, w) 0,
           nuisance = function() return(list(ff = ff, rf = rf, theta = theta)),
           name = "glmer",
           response = function(f) fam$linkinv(f))
}

as.Family.coxph <- function(object, ...) {

    model <- object
    cf <- coef(model)
    ffm <- formula(object)
    cl <- object$call
    env <- .copy_variables(cl, environment(ffm))
    environment(model$formula) <- env
    N <- nrow(model.frame(object))

    fm <- . ~ . + offset(offset_.)
    m0 <- update(model, formula = fm, 
                 subset = subset_., weights = weights_., evaluate = FALSE)
    environment(m0$formula) <- env

    offset_. <- weights_. <- subset_. <- NA

    risk <- function(y, f, w = 1) {

        assign("offset_.", f, env)
        assign("weights_.", w, env)
        assign("subset_.", which(w > 0), env)
        assign("cf_.", cf, env)

        ### do not re-estimate parameters in cf but rather evaluate
        ### partial likelihood with cf and f
        tmp <- eval(update(model, formula = fm, 
                           subset = subset_., 
                           weights = weights_.,
                           init = cf_., 
                           control = coxph.control(iter.max = 0), 
                           evaluate = FALSE), env)
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
        assign("subset_.", which(w > 0), env)

        tmp <- eval(m0, env)
        if (logLik(tmp) < logLik(model))
            warning("risk increase; decrease stepsize nu")
        model <<- tmp
        cf <<- coef(model)

        ### residual wrt a constant for _all_ observations
        ret <- numeric(N)
        ret[w > 0] <- residuals(model, type = "martingale")
        ret
    }

    Family(ngradient = ngradient, risk = risk,
           check_y = function(y) rep(TRUE, N),
           offset = function(y, w) 0,
           nuisance = function() return(list(coefficients = cf)),
           name = "coxph",
           ### conditional survivor function???
           response = function(f) NA)
}
