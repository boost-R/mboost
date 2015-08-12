.onAttach <- function(libname, pkgname) {

    ## get package version
    vers <- packageDescription("mboost")[["Version"]]

    packageStartupMessage("This is mboost ", vers, ". ", "See ",
                          sQuote("package?mboost"), " and ",
                          sQuote('news(package  = "mboost")'), "\n",
                          "for a complete list of changes.\n",
                          appendLF = TRUE)
    return(TRUE)
}

.onLoad <- function(libname, pkgname) {
    if (options("expressions")[[1]] <= 5000)
        options(expressions = 10000)  ### increase maximum number of expressions
    options(mboost_useMatrix = TRUE, ### allow for Matrix package?
            mboost_indexmin = 10000, ### handle ties for n > 10000
            mboost_dftraceS = FALSE,  ### df = trace(S) or df = trace(2 S - StS)
            mboost_lambdaMax = 1e+15,### maximum value for lambda as used in df2lambda
            mboost_Xmonotone = FALSE,### don't force monotonicity in %X%
            mboost_eps = 10e-10, ### factor for dmat in df2lambda
            mboost_check_df2lambda = TRUE) ### check if max(abs(X)) > 10 in df2lambda
                                 ### as this might take a while experts can skip this check
}

.onUnload <- function(libpath) {
    if (options("expressions")[[1]] == 10000)
        options(expression = 5000)  ### decrease maximum number of expressions again
}
