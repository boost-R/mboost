
get_index <- function(x) {

#    if (length(x) == 1 && is.data.frame(x)) {
#        nd <- which(!duplicated(x[[1]]))
#        nd <- nd[complete.cases(x[nd,])]
#        ux <- sort(x[[1]][nd], na.last = TRUE)
#        index <- match(x[[1]], ux)
    ### handle single factors separately
    if (length(x) == 1 && is.factor(x[[1]])) {
         nd <- which(!duplicated(x[[1]]))
         nd <- nd[complete.cases(x[nd,])]
         index <- as.integer(x[[1]])
    } else {
        tmp <- do.call("paste", x)           
        nd <- which(!duplicated(tmp))
        nd <- nd[complete.cases(x[nd,])]
        index <- match(tmp, tmp[nd])
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

isMATRIX <- function(x)
    is.matrix(x) || is(x, "Matrix")

Complete.cases <- function(x) {
    if (isMATRIX(x)) return(rowSums(is.na(x)) == 0)
    complete.cases(x)
}
