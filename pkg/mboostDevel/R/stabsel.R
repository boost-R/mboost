## stabsel method for mboost; requires stabs
stabsel.mboost <- function(x, cutoff, q, PFER,
                           mstop = NULL,
                           folds = subsample(model.weights(x), B = B),
                           B = ifelse(sampling.type == "MB", 100, 50),
                           assumption = c("unimodal", "r-concave", "none"),
                           sampling.type = c("SS", "MB"),
                           papply = mclapply, verbose = TRUE, FWER, eval = TRUE, ...) {

    cll <- match.call()
    p <- length(variable.names(x))
    ibase <- 1:p
    n <- sum(model.weights(x))

    sampling.type <- match.arg(sampling.type)
    if (sampling.type == "MB") {
        assumption <- "none"
    } else {
        assumption <- match.arg(assumption)
    }

    ## check mstop
    if (is.null(mstop)) {
        ## if grid specified in '...'
        if (length(list(...)) >= 1 && "grid" %in% names(list(...))) {
            mstop <- max(list(...)$grid)
        } else {
            mstop <- mstop(x)
        }
    } else {
        if (length(list(...)) >= 1 && "grid" %in% names(list(...)))
            warning(sQuote("grid"), " is ignored if ", sQuote("mstop"), " is specified")
    }

    ## get names for results
    nms <- variable.names(x)
    if (!extends(class(x), "glmboost"))
        nms <- names(nms)

    ## if verbose, count violations (i.e., mstop too small)
    violations <- rep(FALSE, ifelse(sampling.type == "MB", B, 2 * B))

    if (verbose)
        cat("Run stabsel ")

    ## define the fitting function (args.fitfun is not used but needed for
    ## compatibility with run_stabsel
    fit_model <- function(i, folds, q, args.fitfun) {
        if (verbose)
            cat(".")
        ## start by fitting q steps (which are needed at least to obtain q
        ## selected base-learners)
        mod <- update(x, weights = folds[, i])
        mstop(mod) <- q
        ## make sure dispatch works correctly
        class(mod) <- class(x)
        xs <- selected(mod)
        nsel <- length(unique(xs))
        ## now update model until we obtain q different base-learners altogether
        for (m in (q + 1):mstop) {
            if (nsel >= q)
                break
            mstop(mod) <- m
            nsel <- length(unique(selected(mod)))
        }

        ## selected base-learners
        xs <- selected(mod)

        ## logical vector that indicates if the base-learner was selected
        selected <-  logical(p)
        names(selected) <- nms
        selected[unique(xs)] <- TRUE

        if (verbose && sum(selected) < q)
            violations[i] <<- TRUE

        ## compute selection paths
        sel_paths <- matrix(FALSE, nrow = length(nms), ncol = mstop)
        rownames(sel_paths) <- nms
        for (j in 1:length(xs))
            sel_paths[xs[j], j:mstop] <- TRUE

        return(list(selected = selected, path = sel_paths))
    }

    ret <- run_stabsel(fitter = fit_model, args.fitter = list(),
                 n = n, p = p, cutoff = cutoff, q = q,
                 PFER = PFER, folds = folds, B = B, assumption = assumption,
                 sampling.type = sampling.type, papply = papply,
                 verbose = verbose, FWER = FWER, eval = eval,
                 names = unlist(nms), ...)

    if (verbose)
        cat("\n")

    if (!eval)
        return(ret)

    if (any(violations))
        warning(sQuote("mstop"), " too small in ",
                sum(violations), " of the ", length(violations),
                " subsampling replicates to select ", sQuote("q"),
                " base-learners; Increase ", sQuote("mstop"),
                " bevor applying ", sQuote("stabsel"))

    ret$call <- cll
    ret$call[[1]] <- as.name("stabsel")
    class(ret) <- c("stabsel_mboost", "stabsel")
    ret
}

stabsel_parameters.mboost <- function(p, ...) {
    stabsel(p, ..., eval = FALSE)
}
