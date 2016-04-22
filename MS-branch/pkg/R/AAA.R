
### attach package `modeltools'

.onLoad <- function(lib, pkg) {
    if (!require("modeltools"))
        stop("cannot load ", sQuote("modeltools"))
    if (!require("party"))
        stop("cannot load ", sQuote("party"))
    ### stats4 redefines AIC as S4 generic which causes dispatch problems
    if (length(grep("stats4", search())) > 0)
        detach(package:stats4)
    return(TRUE)
}
