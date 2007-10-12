
### 
### Experimental version of gradient boosting with conditional trees 
### as base learner
### 


### Fitting function
blackboost_fit <- function(object, 
                           tree_controls = ctree_control(teststat = "max",
                               testtype = "Teststatistic",
                               mincriterion = 0,
                               maxdepth = 2),
                           fitmem = ctree_memory(object, TRUE), 
                           family = GaussReg(), control = boost_control(), 
                           weights = NULL) {

    ### number of observations in the learning sample
    ### make sure this gets _copied_
    obj <- .Call("copymem",  object, package = "mboost")
    y <- .Call("copymem", party:::get_variables(obj@responses)[[1]], 
               package = "mboost")
    check_y_family(y, family)
    if (is.factor(y)) {
        y <- (2 * (as.numeric(y) - 1)) - 1
        obj@responses <- party:::initVariableFrame(data.frame(y = y), NULL, response = TRUE)
    }
    if (is.null(weights)) weights <- obj@weights
    storage.mode(weights) <- "double"
    oobweights <- as.numeric(weights == 0)


    ### hyper parameters
    mstop <- control$mstop
    risk <- control$risk
    constraint <- control$constraint
    nu <- control$nu
    trace <- control$trace
    tracestep <- options("width")$width / 2

    if (control$center)
        warning("inputs are not centered in ", sQuote("gamboost"))

    ### the ensemble, essentially a list of trees
    ens <- vector(mode = "list", length = mstop)

    ### vector of empirical risks for all boosting iterations
    ### (either in-bag or out-of-bag)
    mrisk <- numeric(mstop)
    mrisk[1:mstop] <- NA

    ### extract negative gradient function
    ngradient <- family@ngradient
    riskfct <- family@risk
    if (!family@weights && any(max(abs(weights - 1))))
        stop(sQuote("family"), " is not able to deal with weights")

    fit <- offset <- family@offset(y, weights)
    u <- ustart <- ngradient(y, fit, weights)

    where <- rep(1, obj@nobs)
    storage.mode(where) <- "integer"

    ### start boosting iteration
    for (m in 1:mstop) {
  
        ### fit tree to residuals
        .Call("R_modify_response", as.double(u), obj@responses, 
              PACKAGE = "party")
        ens[[m]] <- .Call("R_TreeGrow", obj, weights, fitmem, tree_controls,
                          where, PACKAGE = "party")

        ### check if first node is terminal, i.e., if at least 
        ### one split was performed
        if (ens[[m]][[4]])
            warning("could not split root node in iteration ", m, 
                    ", with mincriterion ", sQuote("mincriterion"))

        ### update step
        if (risk == "oobag")
            where <- .Call("R_get_nodeID", ens[[m]], obj@inputs, 0.0,
                           PACKAGE = "party")
        fit <- fit + nu * unlist(.Call("R_getpredictions", ens[[m]], where, 
                                       PACKAGE = "party"))

        ### L2 boost with constraints (binary classification)
        if (constraint)
            fit <- sign(fit) * pmin(abs(fit), 1)

        ### negative gradient vector, the new `residuals'
        u <- ngradient(y, fit, weights)

        ### evaluate risk, either for the learning sample (inbag)
        ### or the test sample (oobag)
        if (risk == "inbag") mrisk[m] <- riskfct(y, fit, weights)
        if (risk == "oobag") mrisk[m] <- riskfct(y, fit, oobweights)

        ### print status information
        if (trace) 
            do_trace(m, risk = mrisk, step = tracestep, width = mstop)
    }

    updatefun <- function(object, control, weights)
        blackboost_fit(object, family = family, tree_controls = tree_controls,
                       fitmem = fitmem, control = control, weights = weights)

    RET <- list(ensemble = ens, 	### list of trees
                fit = fit,              ### vector of fitted values   
                offset = offset,        ### offset
                ustart = ustart,        ### first negative gradients
                risk = mrisk,           ### empirical risks for m = 1, ..., mstop
                control = control,      ### control parameters   
                family = family,        ### family object
                response = y,           ### the response variable
                weights = weights,      ### weights used for fitting
                update = updatefun,     ### a function for fitting with new weights
                tree_controls = tree_controls
    )

    ### save learning sample
    if (control$savedata) RET$data <- object

    ### prediction function (linear predictor only)
    RET$predict <- function(newdata = NULL, mstop = mstop, ...) {

        newinp <- party:::newinputs(object, newdata)

        p <- offset
        for (m in 1:mstop) {
            wh <- .Call("R_get_nodeID", RET$ensemble[[m]], newinp, 0.0, 
                        PACKAGE = "party")
            p <- p + nu * unlist(.Call("R_getpredictions", 
                 RET$ensemble[[m]], wh, PACKAGE = "party"))
        }
        if (constraint) p <- sign(p) * pmin(abs(p), 1)
        return(p)
    }

    class(RET) <- "blackboost"
    return(RET)
}

### methods: subset
"[.blackboost" <- function(x, i, ...) { 
    mstop <- mstop(x)
    if (i == mstop) return(x)
    if (length(i) != 1)
        stop("not a positive integer")
    if (i < 1 || i > mstop)
        warning("invalid number of boosting iterations")
    indx <- 1:min(max(i, 1), mstop)
    x$ensemble <- x$ensemble[indx]
    x$risk <- x$risk[indx]
    x$fit <- x$predict(mstop = max(indx))
    x
}

### methods: prediction
predict.blackboost <- function(object, newdata = NULL, 
                              type = c("lp", "response"), ...) {
    y <- party:::get_variables(object$data@responses)[[1]]
    type <- match.arg(type)
    lp <- object$predict(newdata = newdata, mstop = mstop(object), ...)
    if (type == "response" && is.factor(y))
        return(factor(levels(y)[(lp > 0) + 1], levels = levels(y)))
    return(lp)
}

blackboost <- function(x, ...) UseMethod("blackboost")

blackboost.formula <- function(formula, data = list(), weights = NULL, ...) {

    ### construct design matrix etc.
    object <- party:::ctreedpp(formula, data, ...)
    fitmem <- ctree_memory(object, TRUE)

    ### fit the ensemble
    RET <- blackboost_fit(object, fitmem = fitmem, weights = weights, ...)

    RET$call <- match.call()

    return(RET)
}

blackboost.matrix <- function(x, y, weights = NULL, ... ) {

    object <- party:::LearningSample(object = x, response = y)
    fitmem <- ctree_memory(object, TRUE)

    ### fit the ensemble
    RET <- blackboost_fit(object, fitmem = fitmem, weights = weights, ...)

    RET$call <- match.call()

    return(RET)
}

print.blackboost <- function(x, ...) {

    cat("\n")
    cat("\t Tree-Based Gradient Boosting\n")
    cat("\n")
    if (!is.null(x$call))
    cat("Call:\n", deparse(x$call), "\n\n", sep = "")
    show(x$family)
    cat("\n")
    cat("Number of boosting iterations: mstop =", mstop(x), "\n")
    cat("Step size: ", x$control$nu, "\n")
    cat("Offset: ", x$offset, "\n")
    cat("\n")
    invisible(x)
}
