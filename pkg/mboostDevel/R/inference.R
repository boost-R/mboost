
stabsel <- function(object, cutoff, q, PFER,
                    folds = cv(model.weights(object), type = "subsampling",
                               B = ifelse(error.bound = "MB", 100, 50)),
                    papply = mclapply, verbose = TRUE, FWER,
                    error.bound = c("MB", "SS"), ...) {

    p <- length(variable.names(object))
    ibase <- 1:p

    error.bound <- match.args(error.bound)

    ## only two of the four arguments can be specified
    if ((nmiss <- sum(missing(PFER), missing(cutoff),
                      missing(q), missing(FWER))) != 2) {
        if (nmiss > 2)
            stop("Two of the three argumnets ",
                 sQuote("PFER"), ", ", sQuote("cutoff"), " and ", sQuote("q"),
                 " must be specifed")
        if (nmiss < 2)
            stop("Only two of the three argumnets ",
                 sQuote("PFER"), ", ", sQuote("cutoff"), " and ", sQuote("q"),
                 " can be specifed at the same time")
    }

    if (!missing(FWER)) {
        if (!missing(PFER))
            stop(sQuote("FWER"), " and ", sQuote("PFER"),
                 " cannot be spefified at the same time")
        PFER <- FWER
        warning(sQuote("FWER"), " is deprecated. Use ", sQuote("PFER"),
                " instead.")
    }

    if ((!missing(PFER) || !missing(FWER)) && PFER < 0)
        stop(sQuote("PFER"), " must be greater 0")

    if (!missing(cutoff) && (cutoff < 0.5 | cutoff > 1))
        stop(sQuote("cutoff"), " must be between 0.5 and 1")

    if (!missing(q)) {
        if (p < q)
            stop("Average number of selected base-learners ", sQuote("q"),
                 " must be smaller \n  than the number of base-learners",
                 " specified in the model ", sQuote("object"))
        if (q < 0)
            stop("Average number of selected base-learners ", sQuote("q"),
                 " must be greater 0")
    }

    if (missing(cutoff)) {
        if (error.bound == "MB") {
            cutoff <- min(0.9, tmp <- (q^2 / (PFER * p) + 1) / 2)
            upperbound <- q^2 / p / (2 * cutoff - 1)
        } else {
            objective <- function(cutoff) {
                PFER / p - minD(q, p, cutoff, B)
            }
            root <- uniroot(objective, lower = 0.5, upper = 0.9)$root
            cutoff <- floor(root * 2 * B) / (2* B)
            upperbound <- minD(q, p, cutoff, B) * p
        }
        upperbound <- signif(upperbound, 3)
        if (verbose && tmp > 0.9 && upperbound - PFER > PFER/2) {
            warning("Upper bound for PFER > ", PFER,
                    " for the given value of ", sQuote("q"),
                    " (true upper bound = ", round(upperbound, 2), ")")
        }
    }

    if (missing(q)) {
        if (error.bound == "MB") {
            q <- ceiling(sqrt(PFER * (2 * cutoff - 1) * p))
            upperbound <- q^2 / p / (2 * cutoff - 1)
        } else {
            objective <- function(q) {
                PFER / p - minD(q, p, cutoff, B)
            }
            root <- uniroot(objective, lower = 1,
                            upper = min(sqrt((B - 1) / (2 * B) * p^2),
                            (B - 1) / (2 * B) * p))$root
            q <- ceiling(root)
            upperbound <- minD(q, p, cutoff, B) * p
        }
        upperbound <- signif(upperbound, 3)
        if (verbose && upperbound - PFER > PFER/2)
            warning("Upper bound for PFER > ", PFER,
                    " for the given value of ", sQuote("cutoff"),
                    " (true upper bound = ", upperbound, ")")
    }

    if (missing(PFER)) {
        if (error.bound == "MB") {
            upperbound <- PFER <- q^2 / p / (2 * cutoff - 1)
        } else {
            upperbound <- PFER <- minD(q, p, cutoff, B) * p
        }
        upperbound <- signif(upperbound, 3)
    }
    if (verbose && PFER >= p)
        warning("Upper bound for PFER larger than the number of base-learners.")

    fun <- function(model) {
        xs <- selected(model)
        qq <- sapply(1:length(xs), function(x) length(unique(xs[1:x])))
        xs[qq > q] <- xs[1]
        xs
    }
    if (error.bound == "SS") {
        ## use complementary pairs
        folds <- cbind(folds, model.weights(object) - folds)
    }
    ss <- cvrisk(object, fun = fun,
                 folds = folds,
                 papply = papply, ....)

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
        m <- mstop(object)
    }
    ret <- matrix(0, nrow = length(ibase), ncol = m)
    for (i in 1:length(ss)) {
        tmp <- sapply(ibase, function(x)
            ifelse(x %in% ss[[i]], which(ss[[i]] == x)[1], m + 1))
        ret <- ret + t(sapply(tmp, function(x) c(rep(0, x - 1), rep(1, m - x + 1))))
    }

    phat <- ret / length(ss)
    rownames(phat) <- names(variable.names(object))
    if (extends(class(object), "glmboost"))
        rownames(phat) <- variable.names(object)
    ret <- list(phat = phat, selected = which((mm <- apply(phat, 1, max)) >= cutoff),
                max = mm, cutoff = cutoff, q = q, PFER = upperbound)
    class(ret) <- "stabsel"
    ret
}

print.stabsel <- function(x, decreasing = FALSE, ...) {

    cat("\tStability Selection\n")
    if (length(x$selected) > 0) {
        cat("\nSelected base-learners:\n")
        print(x$selected)
    } else {
        cat("\nNo base-learner selected\n")
    }
    cat("\nSelection probabilities:\n")
    print(sort(x$max[x$max > 0], decreasing = decreasing))
    cat("\nCutoff: ", x$cutoff, "; ", sep = "")
    cat("q: ", x$q, "; ", sep = "")
    cat("PFER: ", x$PFER, "\n\n")
    invisible(x)
}

plot.stabsel <- function(x, main = deparse(x$call), col = NULL, ...) {

    h <- x$phat
    h <- h[rowSums(h) > 0, , drop = FALSE]
    if (is.null(col))
        col <- hcl(h = 40, l = 50, c = h[,ncol(h)] / max(h) * 490)
    matplot(t(h), type = "l", lty = 1, xlab = "Number of boosting iterations",
            ylab = "Selection Probability", main = main, col = col, ylim = c(0, 1), ...)
    abline(h = x$cutoff, lty = 1, col = "lightgray")
    axis(4, at = x$phat[rowSums(x$phat) > 0, ncol(x$phat)],
         labels = rownames(x$phat)[rowSums(x$phat) > 0], las = 1)
}


fitsel <- function(object, newdata = NULL, which = NULL, ...) {
    fun <- function(model) {
        tmp <- predict(model, newdata = newdata,
                       which = which, agg = "cumsum")
        ret <- c()
        for (i in 1:length(tmp))
            ret <- rbind(ret, tmp[[i]])
        ret
    }
    ss <- cvrisk(object, fun = fun, ...)
    ret <- matrix(0, nrow = nrow(ss[[1]]), ncol = ncol(ss[[1]]))
    for (i in 1:length(ss))
        ret <- ret + sign(ss[[i]])
    ret <- abs(ret) / length(ss)
    ret
}


### Modified version of the code accompanying the paper:
###   Shah, R. D. and Samworth, R. J. (2013), Variable selection with error
###   control: Another look at Stability Selection, J. Roy. Statist. Soc., Ser.
###   B, 75, 55-80. DOI: 10.1111/j.1467-9868.2011.01034.x
###
### Original code available from
###   http://www.statslab.cam.ac.uk/~rds37/papers/r_concave_tail.R
### or
###   http://www.statslab.cam.ac.uk/~rjs57/r_concave_tail.R

D <- function(theta, which, B, r) {
    ## If pi = ceil{ B * 2 * eta} / B + 1/B,..., 1 return the tail probability.
    ## If pi < ceil{ B * 2 * eta} / B return 1

    if (which <= 0)
        return(1)
    ### pi muss ein vielfaches von 1/(2 * B) sein, oder?

    s <- 1/r
    thetaB <- theta * B
    k_start <- (ceiling(2 * thetaB) + 1)
    if(k_start > B)
        stop("theta to large")

    Find.a <- function(prev_a)
        uniroot(Calc.a, lower = 0.0001, upper = prev_a,
                tol = .Machine$double.eps^0.75)$root

    Calc.a <- function(a) {
        denom <- sum((a + 0:k)^s)
        num <- sum((0:k) * (a + 0:k)^s)
        num / denom - thetaB
    }

    OptimInt <- function(a) {
        num <- (k + 1 - thetaB) * sum((a + 0:(t-1))^s)
        denom <- sum((k + 1 - (0:k)) * (a + 0:k)^s)
        1 - num / denom
    }

    ## initialize a
    a_vec <- rep(100000, B)

    ## compute a values
    for(k in k_start:B)
        a_vec[k] <- Find.a(a_vec[k-1])

    t <- which
    cur_optim <- rep(0, B)
    for (k in k_start:(B-1))
        cur_optim[k] <- optimize(f=OptimInt, lower = a_vec[k+1],
                                 upper = a_vec[k], maximum  = TRUE)$objective
    return(max(cur_optim))
}

minD <- function(q, p, pi, B, r = c(-1/2, -1/4)) {
    which <- ceiling(signif(pi / (1/(2* B)), 10))
    #if ((which) %% 1 != 0)
    #    stop(sQuote("pi"), " must be a multiple of 1/(2 * B)")

    maxQ <- min(sqrt((B - 1) / (2 * B) * p^2),
                (B - 1) / (2 * B) * p)
    if (q > maxQ)
        stop(sQuote("q"), " must be <= ", maxQ)
    min(c(1, D(q^2 / p^2, which - B, B, r[1]), D(q / p, which , 2*B, r[2])))
}
