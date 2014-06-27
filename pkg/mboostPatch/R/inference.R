
stabsel <- function(object, cutoff, q, PFER,
                    folds = cv(model.weights(object), type = "subsampling",
                               B = ifelse(sampling.type == "MB", 100, 50)),
                    assumption = c("unimodal", "r-concave", "none"),
                    sampling.type = c("SS", "MB"),
                    papply = mclapply, verbose = TRUE, FWER, eval = TRUE, ...) {

    call <- match.call()
    p <- length(variable.names(object))
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
                               assumption = assumption)
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
        folds <- cbind(folds, model.weights(object) - folds)
    }
    ss <- cvrisk(object, fun = fun,
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
                max = mm, cutoff = cutoff, q = q, PFER = PFER,
                sampling.type = sampling.type, assumption = assumption,
                call = call)
    class(ret) <- "stabsel"
    ret
}

stabsel_parameters <- function(p, cutoff, q, PFER,
                               B = ifelse(sampling.type == "MB", 100, 50),
                               assumption = c("unimodal", "r-concave", "none"),
                               sampling.type = c("SS", "MB"),
                               verbose = FALSE, FWER) {

    sampling.type <- match.arg(sampling.type)
    if (sampling.type == "MB")
        assumption <- "none"
    else
        assumption <- match.arg(assumption)


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
        if (assumption == "none") {
            cutoff <- min(1, tmp <- (q^2 / (PFER * p) + 1) / 2)
            upperbound <- q^2 / p / (2 * cutoff - 1)
        } else {
            if (assumption == "unimodal") {
                cutoff <- tmp <- optimal_cutoff(p, q, PFER, B,
                                                assumption = assumption)
                upperbound <- q^2 / p / um_const(cutoff, B, theta = q/p)
            } else {
                cutoff <- tmp <- optimal_cutoff(p, q, PFER, B,
                                                assumption = assumption)
                upperbound <- minD(q, p, cutoff, B) * p
            }
        }
        upperbound <- signif(upperbound, 3)
        if (verbose && tmp > 0.9 && upperbound - PFER > PFER/2) {
            warning("Upper bound for PFER > ", PFER,
                    " for the given value of ", sQuote("q"),
                    " (true upper bound = ", round(upperbound, 2), ")")
        }
    }

    if (missing(q)) {
        if (assumption == "none") {
            q <- floor(sqrt(PFER * (2 * cutoff - 1) * p))
            upperbound <- q^2 / p / (2 * cutoff - 1)
        } else {
            if (assumption == "unimodal") {
                q <- optimal_q(p, cutoff, PFER, B, assumption = assumption)
                upperbound <- q^2 / p / um_const(cutoff, B, theta = q/p)
            } else {
                q <- optimal_q(p, cutoff, PFER, B, assumption = assumption)
                upperbound <- minD(q, p, cutoff, B) * p
            }
        }
        upperbound <- signif(upperbound, 3)
        if (verbose && upperbound - PFER > PFER/2)
            warning("Upper bound for PFER > ", PFER,
                    " for the given value of ", sQuote("cutoff"),
                    " (true upper bound = ", upperbound, ")")
    }

    if (missing(PFER)) {
        if (assumption == "none") {
            upperbound <- PFER <- q^2 / p / (2 * cutoff - 1)
        } else {
            if (assumption == "unimodal") {
                upperbound <- PFER <- q^2 / p / um_const(cutoff, B, theta = q/p)
            } else {
                upperbound <- PFER <- minD(q, p, cutoff, B) * p
            }
        }
        upperbound <- signif(upperbound, 3)
    }

    if (verbose && PFER >= p)
        warning("Upper bound for PFER larger than the number of base-learners.")

    res <- list(cutoff = cutoff, q = q, PFER = upperbound,
                sampling.type = sampling.type, assumption = assumption)
    class(res) <- "stabsel_parameters"
    res
}

print.stabsel <- function(x, decreasing = FALSE, ...) {

    cat("\tStability Selection")
    if (x$assumption == "none")
        cat(" without further assumptions\n")
    if (x$assumption == "unimodal")
        cat(" with unimodality assumption\n")
    if (x$assumption == "r-concave")
        cat(" with r-concavity assumption\n")
    if (length(x$selected) > 0) {
        cat("\nSelected base-learners:\n")
        print(x$selected)
    } else {
        cat("\nNo base-learner selected\n")
    }
    cat("\nSelection probabilities:\n")
    print(sort(x$max[x$max > 0], decreasing = decreasing))
    cat("\n")
    print.stabsel_parameters(x, heading = FALSE)
    cat("\n")
    invisible(x)
}

print.stabsel_parameters <- function(x, heading = TRUE, ...) {
    if (heading) {
        cat("Stability Selection")
        if (x$assumption == "none")
            cat(" without further assumptions\n")
        if (x$assumption == "unimodal")
            cat(" with unimodality assumption\n")
        if (x$assumption == "r-concave")
            cat(" with r-concavity assumption\n")
    }
    cat("Cutoff: ", x$cutoff, "; ", sep = "")
    cat("q: ", x$q, "; ", sep = "")
    if (x$sampling.type == "MB")
        cat("PFER: ", x$PFER, "\n")
    else
        cat("PFER(*): ", x$PFER,
            "\n   (*) or expected number of low selection probability variables\n")
    invisible(x)
}

plot.stabsel <- function(x, main = deparse(x$call), type = c("paths", "maxsel"),
                         col = NULL, ymargin = 10, np = sum(x$max > 0),
                         labels = NULL, ...) {

    type <- match.arg(type)

    if (is.null(col))
        col <- hcl(h = 40, l = 50, c = x$max / max(x$max) * 490)

    if (type == "paths") {
        ## if par(mar) not set by user ahead of plotting
        if (all(par()[["mar"]] == c(5, 4, 4, 2) + 0.1))
            ..old.par <- par(mar = c(5, 4, 4, ymargin) + 0.1)
        h <- x$phat
        h <- h[rowSums(h) > 0, , drop = FALSE]
        matplot(t(h), type = "l", lty = 1,
                xlab = "Number of boosting iterations",
                ylab = "Selection probability",
                main = main, col = col[x$max > 0], ylim = c(0, 1), ...)
        abline(h = x$cutoff, lty = 1, col = "lightgray")
        if (is.null(labels))
            labels <- rownames(x$phat)
        axis(4, at = x$phat[rowSums(x$phat) > 0, ncol(x$phat)],
             labels = labels[rowSums(x$phat) > 0], las = 1)
    } else {
        ## if par(mar) not set by user ahead of plotting
        if (all(par()[["mar"]] == c(5, 4, 4, 2) + 0.1))
            ..old.par <- par(mar = c(5, ymargin, 4, 2) + 0.1)
        if (np > length(x$max))
            stop(sQuote("np"), "is set too large")
        inc_freq <- x$max  ## inclusion frequency
        plot(tail(sort(inc_freq), np), 1:np,
             type = "n", yaxt = "n", xlim = c(0, 1),
             ylab = "", xlab = expression(hat(pi)),
             main = main, ...)
        abline(h = 1:np, lty = "dotted", col = "grey")
        points(tail(sort(inc_freq), np), 1:np, pch = 19,
               col = col[tail(order(inc_freq), np)])
        if (is.null(labels))
            labels <- names(x$max)
        axis(2, at = 1:np, labels[tail(order(inc_freq), np)], las = 2)
        ## add cutoff
        abline(v = x$cutoff, col = "grey")
    }
    if (exists("..old.par"))
        par(..old.par) # reset plotting settings
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
    ## compute upper tail of r-concave distribution function
    ## If q = ceil{ B * 2 * theta} / B + 1/B,..., 1 return the tail probability.
    ## If q < ceil{ B * 2 * theta} / B return 1

    s <- 1/r
    thetaB <- theta * B
    k_start <- (ceiling(2 * thetaB) + 1)

    if (which < k_start)
        return(1)

    if(k_start > B)
        stop("theta to large")

    Find.a <- function(prev_a)
        uniroot(Calc.a, lower = 0.00001, upper = prev_a,
                tol = .Machine$double.eps^0.75)$root

    Calc.a <- function(a) {
        denom <- sum((a + 0:k)^s)
        num <- sum((0:k) * (a + 0:k)^s)
        num / denom - thetaB
    }

    OptimInt <- function(a, t, k, thetaB, s) {
        num <- (k + 1 - thetaB) * sum((a + 0:(t-1))^s)
        denom <- sum((k + 1 - (0:k)) * (a + 0:k)^s)
        1 - num / denom
    }

    ## initialize a
    a_vec <- rep(100000, B)

    ## compute a values
    for(k in k_start:B)
        a_vec[k] <- Find.a(a_vec[k-1])

    cur_optim <- rep(0, B)
    for (k in k_start:(B-1))
        cur_optim[k] <- optimize(f=OptimInt, lower = a_vec[k+1],
                                 upper = a_vec[k],
                                 t = which, k = k, thetaB = thetaB, s = s,
                                 maximum  = TRUE)$objective
    return(max(cur_optim))
}

## minD function for error bound in case of r-concavity
minD <- function(q, p, pi, B, r = c(-1/2, -1/4)) {
    ## get the integer valued multiplier W of
    ##   pi = W * 1/(2 * B)
    which <- ceiling(signif(pi / (1/(2* B)), 10))
    maxQ <- maxQ(p, B)
    if (q > maxQ)
        stop(sQuote("q"), " must be <= ", maxQ)
    min(c(1, D(q^2 / p^2, which - B, B, r[1]), D(q / p, which , 2*B, r[2])))
}

## function to find optimal cutoff in stabsel (when sampling.type = "SS")
optimal_cutoff <- function(p, q, PFER, B, assumption = "unimodal") {
    if (assumption == "unimodal") {
        ## cutoff values can only be multiples of 1/(2B)
        cutoffgrid <- 1/2 + (2:B)/(2*B)
        c_min <- min(0.5 + (q/p)^2, 0.5 + 1/(2*B) + 0.75 * (q/p)^2)
        cutoffgrid <- cutoffgrid[cutoffgrid > c_min]
        upperbound <- rep(NA, length(cutoffgrid))
        for (i in 1:length(cutoffgrid))
            upperbound[i] <- q^2 / p / um_const(cutoffgrid[i], B, theta = q/p)
        cutoff <- cutoffgrid[upperbound < PFER][1]
        return(cutoff)
    } else {
        ## cutoff values can only be multiples of 1/(2B)
        cutoff <- (2*B):1/(2*B)
        cutoff <- cutoff[cutoff >= 0.5]
        for (i in 1:length(cutoff)) {
            if (minD(q, p, cutoff[i], B) * p > PFER) {
                if (i == 1)
                    cutoff <- cutoff[i]
                else
                    cutoff <- cutoff[i - 1]
                break
            }
        }
        return(tail(cutoff, 1))
    }
}

## function to find optimal q in stabsel (when sampling.type = "SS")
optimal_q <- function(p, cutoff, PFER, B, assumption = "unimodal") {
    if (assumption == "unimodal") {
        if (cutoff <= 0.75) {
            upper_q <- max(p * sqrt(cutoff - 0.5),
                           p * sqrt(4/3 * (cutoff - 0.5 - 1/(2*B))))
            ## q must be an integer < upper_q
            upper_q <- ceiling(upper_q - 1)
        } else {
            upper_q <- p
        }
        q <- uniroot(function(q)
                     q^2 / p / um_const(cutoff, B, theta = q/p) - PFER,
                     lower = 1, upper = upper_q)$root
        return(floor(q))
    } else {
        for (q in 1:maxQ(p, B)) {
            if (minD(q, p, cutoff, B) * p > PFER) {
                q <- q - 1
                break
            }
        }
        return(max(1, q))
    }
}

## obtain maximal value possible for q
maxQ <- function(p, B) {
    if(B <= 1)
        stop("B must be at least 2")

    fact_1 <- 4 * B / p
    tmpfct <- function(q)
        ceiling(q * fact_1) + 1 - 2 * B

    res <- tmpfct(1:p)
    length(res[res < 0])
}

## obtain constant for unimodal bound
um_const <- function(cutoff, B, theta) {
    if (cutoff <= 3/4) {
        if (cutoff < 1/2 + min(theta^2, 1 / (2*B) + 3/4 * theta^2))
            stop ("cutoff out of bounds")
        return( 2 * (2 * cutoff - 1 - 1/(2*B)) )
    } else {
        if (cutoff > 1)
            stop ("cutoff out of bounds")
        return( (1 + 1/B)/(4 * (1 - cutoff + 1 / (2*B))) )
    }
}
