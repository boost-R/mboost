## stabsel method for mboost; requires stabs
stabsel.mboost <- function(x, cutoff, q, PFER,
                    folds = subsample(model.weights(x), B = B),
                    B = ifelse(sampling.type == "MB", 100, 50),
                    assumption = c("unimodal", "r-concave", "none"),
                    sampling.type = c("SS", "MB"),
                    papply = mclapply, verbose = TRUE, FWER, eval = TRUE, ...) {

    cll <- match.call()
    p <- length(variable.names(x))
    ibase <- 1:p

    sampling.type <- match.arg(sampling.type)
    if (sampling.type == "MB")
        assumption <- "none"
    else
        assumption <- match.arg(assumption)

    B <- ncol(folds)

    pars <- stabsel_parameters(p = p, cutoff = cutoff, q = q,
                               PFER = PFER, B = B,
                               verbose = verbose, sampling.type = sampling.type,
                               assumption = assumption, FWER = FWER)
    ## return parameter combination only if eval == FALSE
    if (!eval)
        return(pars)

    cutoff <- pars$cutoff
    q <- pars$q
    PFER <- pars$PFER

    fun <- function(model) {
        xs <- selected(model)
        qq <- sapply(1:length(xs), function(x) length(unique(xs[1:x])))
        xs[qq > q] <- xs[1]
        xs
    }
    if (sampling.type == "SS") {
        ## use complementary pairs
        folds <- cbind(folds, model.weights(x) - folds)
    }
    ss <- cvrisk(x, fun = fun,
                 folds = folds,
                 papply = papply, ...)

    if (verbose){
        qq <- sapply(ss, function(x) length(unique(x)))
        sum_of_violations <- sum(qq < q)
        if (sum_of_violations > 0)
            warning(sQuote("mstop"), " too small in ",
                    sum_of_violations, " of the ", ncol(folds),
                    " subsampling replicates to select ", sQuote("q"),
                    " base-learners; Increase ", sQuote("mstop"),
                    " bevor applying ", sQuote("stabsel"))
    }


    ## if grid specified in '...'
    if (length(list(...)) >= 1 && "grid" %in% names(list(...))) {
        m <- max(list(...)$grid)
    } else {
        m <- mstop(x)
    }
    ret <- matrix(0, nrow = length(ibase), ncol = m)
    for (i in 1:length(ss)) {
        tmp <- sapply(ibase, function(x)
            ifelse(x %in% ss[[i]], which(ss[[i]] == x)[1], m + 1))
        ret <- ret + t(sapply(tmp, function(x) c(rep(0, x - 1), rep(1, m - x + 1))))
    }

    phat <- ret / length(ss)
    rownames(phat) <- names(variable.names(x))
    if (extends(class(x), "glmboost"))
        rownames(phat) <- variable.names(x)
    ret <- list(phat = phat, selected = which((mm <- apply(phat, 1, max)) >= cutoff),
                max = mm, cutoff = cutoff, q = q, PFER = PFER, p = p, B = B, 
                sampling.type = sampling.type, assumption = assumption,
                call = cll)
    ret$call[[1]] <- as.name("stabsel")
    class(ret) <- c("stabsel_mboost", "stabsel")
    ret
}

stabsel_parameters.mboost <- function(p, ...) {
    stabsel(p, ..., eval = FALSE)
}
