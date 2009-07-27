
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
