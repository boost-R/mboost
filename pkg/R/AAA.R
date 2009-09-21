
.onLoad <- function(lib, pkg) {
    if (!require("Matrix"))
        stop("cannot load ", sQuote("Matrix"))
    ### stats4 redefines AIC as S4 generic which causes dispatch problems
    if (length(grep("stats4", search())) > 0)
        detach(package:stats4)
    return(TRUE)
}

#.onLoad <- function(libname, pkgname)
#{
#    op <- options()
#    op.mboost <-
#        list(lambdamax = 1e+06,
#             bsMatrix = c(500, 50),
#             dfindex = 10000)
#    toset <- !(names(op.mboost) %in% names(op))
#    if(any(toset)) options(op.mboost[toset])   
#}
