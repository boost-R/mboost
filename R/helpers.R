
### try to find duplicated entries / observations with missing values
### for more efficient memory handling
get_index <- function(x) {

    if (isMATRIX(x)) {
        ### handle missing values only
        cc <- Complete.cases(x)
        nd <- which(cc)
        index <- match(1:nrow(x), nd)
    } else {
        ### handle single variables (factors / numerics) factors
        if (length(x) == 1) {
            x <- x[[1]]
            nd <- which(!duplicated(x))
            nd <- nd[complete.cases(x[nd])]
            index <- match(x, x[nd])
        ### go for data.frames with >= 2 variables
        } else {
            tmp <- do.call("paste", x)
            nd <- which(!duplicated(tmp))
            nd <- nd[complete.cases(x[nd,])]
            index <- match(tmp, tmp[nd])
        }
    }
    return(list(nd, index))
}

rescale_weights <- function(w) {
    if (max(abs(w - floor(w))) < sqrt(.Machine$double.eps))
        return(w)
    return(w / sum(w) * sum(w > 0))
}

### check measurement scale of response for some losses
check_y_family <- function(y, family)
    family@check_y(y)

### check for negative gradient corresponding to L2 loss
### <FIXME> better check? </FIXME>
checkL2 <- function(object)
    isTRUE(all.equal(attributes(object$family)[1],
                     attributes(Gaussian())[1]))

### check for classical or Matrix matrices
isMATRIX <- function(x)
    is.matrix(x) || is(x, "Matrix")

### rows without missings in Matrices, matrices and data.frames
Complete.cases <- function(x) {
    if (isMATRIX(x))
        return(rowSums(is.na(x)) == 0)
    complete.cases(x)
}

### fractional polynomials transformation
### all powers `p' of `x', all powers `p' of `x' times `log(x)' and `log(x)'
### see Sauerbrei & Royston (1999), JRSS A (162), 71--94
FP <- function(x, p = c(-2, -1, -0.5, 0.5, 1, 2, 3), scaling = TRUE) {
    xname <- deparse(substitute(x))
    if (is.logical(scaling)) {
        shift <- 0
        scale <- 1
        if (scaling) {
            if (min(x) <= 0) {
                z <- diff(sort(x))
                shift <- min(z[z > 0]) - min(x)
                shift <- ceiling(shift * 10)/10
            }
            range <- mean(x + shift)
            scale <- 10^(sign(log10(range)) * round(abs(log10(range))))
        }
        scaling <- c(shift, scale)
    } else {
        stopifnot(is.numeric(scaling) && length(scaling) == 2)
    }
    x <- (x + scaling[1]) / scaling[2]
    stopifnot(all(x > 0))

    Xp <- sapply(p, function(p) x^p)
    Xp <- cbind(Xp, log(x))
    Xp <- cbind(Xp, Xp * log(x))
    colnames(Xp) <- c(paste(xname, "^", p, sep = ""),
                      paste("log(", xname, ")", sep = ""),
                      paste("log(", xname, ")", xname, "^", p, sep = ""),
                      paste("log(", xname, ")^2", sep = ""))
    attr(Xp, "p") <- p
    attr(Xp, "scaling") <- scaling
    Xp
}

### trace boosting iterations
do_trace <- function(m, mstop, risk, step = options("width")$width / 2,
                     width = 1000) {
    m <- m - mstop
    if (m != width) {
        if ((m - 1) %/% step == (m - 1) / step) {
            mchr <- formatC(m+mstop, format = "d", width = nchar(width) + 1,
                            big.mark = "'")
            cat(paste("[", mchr, "] ",sep = ""))
        } else {
            if ((m %/% step != m / step)) {
                cat(".")
            } else {
                cat(" -- risk:", risk[m+mstop], "\n")
            }
        }
    } else {
        cat("\nFinal risk:", risk[m+mstop], "\n")
    }
}

bhatmat <- function(n, H, xselect, fitm, fW) {

    B <- matrix(0, nrow = n, ncol = n)
    I <- diag(n)
    tr <- numeric(length(xselect))

    for (m in 1:length(xselect)) {
        B <- B + (H[[xselect[m]]] * fW(fitm[,m])) %*% (I - B)
        tr[m] <- sum(diag(B))
    }
    list(hatmatrix = B, trace = tr)
}

nnls1D <- function(XtX, X, y) {
    my <- switch(attr(X, "Ts_constraint"),  "increasing" = {
        ### first column is intercept
        stopifnot(max(abs(X[,1,drop = TRUE] - 1)) < sqrt(.Machine$double.eps))
        min(y)
    }, "decreasing" = {
        stopifnot(max(abs(X[,1,drop = TRUE] + 1)) < sqrt(.Machine$double.eps))
        max(y)
    })
    y <- y - my
    cf <- nnls(XtX, crossprod(X, y))$x
    cf[1] <- cf[1] + switch(attr(X, "Ts_constraint"),  "increasing" = my,
                                                    "decreasing" = -my)
    cf
}

nnls2D <- function(X, XtX, Y) {

    constr <- which(c(!is.null(attr(X$X1, "Ts_constraint")),
                      !is.null(attr(X$X2, "Ts_constraint"))))
    if (length(constr) == 2)
            stop("only one dimension may be subject to constraints")
    Xc <- paste("X", constr, sep = "")
    my <- switch(attr(X[[Xc]], "Ts_constraint"),  "increasing" = {
        ### first column is intercept
        stopifnot(max(abs(X[[Xc]][,1,drop = TRUE] - 1)) < sqrt(.Machine$double.eps))
        min(Y)
    }, "decreasing" = {
        stopifnot(max(abs(X[[Xc]][,1,drop = TRUE] + 1)) < sqrt(.Machine$double.eps))
        max(Y)
    })
    Y <- Y - my

    XWY <- as.vector(crossprod(X$X1, Y) %*% X$X2)
    cf <- nnls(XtX, matrix(as(XWY, "matrix"), ncol = 1))$x
    cf <- matrix(cf, nrow = ncol(X$X1))
    if (constr == 1) cf[1,] <- cf[1,] + switch(attr(X[[Xc]], "Ts_constraint"),
                                               "increasing" = my,
                                               "decreasing" = -my)
    if (constr == 2) cf[,1] <- cf[,1] + switch(attr(X[[Xc]], "Ts_constraint"),
                                               "increasing" = my,
                                               "decreasing" = -my)
    cf
}

## function returns either differences or difference matrix
## needed in bmono
differences <- function(constraint, ncol = NULL) {
    if (length(constraint) == 1) {
        diff <- switch(constraint,
                       none = NULL,
                       positive = 0,
                       negative = 0,
                       increasing = 1,
                       decreasing = -1,
                       convex = 2,
                       concave = -2)
    } else {
        diff <- lapply(constraint, function(x)
                       switch(x,
                              none = NULL,
                              positive = 0,
                              negative = 0,
                              increasing = 1,
                              decreasing = -1,
                              convex = 2,
                              concave = -2))
    }
    if (is.null(ncol))
        return(diff)
    #### and thus quit function

    make_diffs <- function(diff, ncol) {
        if (is.null(diff))
            return(NULL)
        if (diff != 0) {
            D <- sign(diff) * diff(diag(ncol), differences = abs(diff))
        } else {
            D <- ifelse(constraint == "positive", 1, -1) * diag(ncol)
        }
        return(D)
    }

    if (length(diff) == 1)
        return(make_diffs(diff, ncol))
    #### and thus quit function

    D <- vector(mode = "list", length =2)
    tmpD <- make_diffs(diff[[1]], ncol[[1]])
    if (!is.null(tmpD))
        ## to make Matrix quiet
        D[[1]] <- suppressMessages(kronecker(tmpD, diag(ncol[[2]])))
    tmpD <- make_diffs(diff[[2]], ncol[[2]])
    if (!is.null(tmpD))
        ## to make Matrix quiet
        D[[2]] <- suppressMessages(kronecker(diag(ncol[[1]]), tmpD))
    return(D)
}


## least squares with equalities and inequalities
solveLSEI <- function(XtX, Xty, D = NULL) {
    if (is.null(D) || all(sapply(D, is.null)))
        return(solve(XtX, Xty, LINPACK = FALSE))

    if (!isMATRIX(D)) ## i.e. D is a list with 2 entries
        D <- rbind(D[[1]], D[[2]])
    ## NOTE: Currently both constraints get the same weight
    cf <- solve.QP(Dmat = XtX, dvec = as.vector(Xty), Amat = t(D),
                   bvec = rep(0, nrow(D)))$solution
    cf
}

check_newdata <- function(newdata, blg, mf, to.data.frame = TRUE) {
    nm <- names(blg)
    if (!all(nm %in% names(newdata)))
        stop(sQuote("newdata"),
             " must contain all predictor variables,",
             " which were used to specify the model.")
    if (!class(newdata) %in% c("list", "data.frame"))
        stop(sQuote("newdata"), " must be either a data.frame or a list")
    if (any(duplicated(nm)))  ## removes duplicates
        nm <- unique(nm)
    cl1 <- sapply(newdata[nm], class)
    cl2 <- sapply(mf, class)

    if (!isTRUE(all.equal(cl1, cl2))) {
        # idx <- which(cl1 != cl2)
        idx <- which(!sapply(1:length(cl1),
                             function(i) all(cl1[[i]] == cl2[[i]])))
        ## classes can be different when one is integer and the other is double
        if (!all(sapply(newdata[nm][idx], is.numeric)) ||
            !all(sapply(mf[idx], is.numeric)))
            warning("Some variables in ", sQuote("newdata"),
                    " do not have the same class as in the original data set",
                    call. = FALSE)
    }
    ## subset data
    mf <- newdata[nm]
    if (is.list(mf) && to.data.frame)
        mf <- as.data.frame(mf)
    return(mf)
}
