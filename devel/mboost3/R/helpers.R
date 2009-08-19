
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
            if (is.factor(x)) {
                index <- as.integer(x)
            } else {
                index <- match(x, x[nd])
            }
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
checkL2 <- function(object)
    isTRUE(all.equal(attributes(object$family)[-c(4:8)],
                     attributes(GaussReg())[-c(4:8)]))

### check for classical or Matrix matrices
isMATRIX <- function(x)
    is.matrix(x) || is(x, "Matrix")

### rows without missings in Matrices, matrices and data.frames
Complete.cases <- function(x) {
    if (isMATRIX(x)) return(rowSums(is.na(x)) == 0)
    complete.cases(x)
}

MYapply <- function(X, FUN, parallel = FALSE, ...) {

    myapply <- lapply
    if (parallel && .Platform$OS.type == "unix") {
        if (!multicore:::isChild()) 
            myapply <- mclapply
    }
    myapply(X, FUN, ...)
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
do_trace <- function(m, risk, step = options("width")$width / 2,
                     width = 1000) {

    if (m != width) {
        if ((m - 1) %/% step == (m - 1) / step) {
            mchr <- formatC(m, format = "d", width = nchar(width) + 1,
                            big.mark = "'")
            cat(paste("[", mchr, "] ",sep = ""))
        } else {
            if ((m %/% step != m / step)) {
                cat(".")
            } else {
                cat(" -- risk:", risk[m], "\n")
            }
        }    
    } else {
        cat("\nFinal risk:", risk[m], "\n")
    }
}
