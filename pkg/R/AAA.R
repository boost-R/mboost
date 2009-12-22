
.onAttach <- function(libname, pkgname) {

    sup <- file.path(libname, pkgname, "startup.txt")
    if (file.exists(sup) || !interactive()) return(TRUE)
    if (!suppressWarnings(file.create(sup))) return(TRUE)
    file.remove(sup)
    version <- packageDescription(pkg = pkgname)$Version
    txt <- c("\n",
             paste("	Welcome to ", sQuote("mboost"), 
                   " version ", version, "!", sep = ""),
             "\n",
             "The user-interface had to be changed in some places.",
             paste("Most important, subsetting an", sQuote("mboost"), 
                   "object changes this object now."),
             "Please read the NEWS file, consult the documentation and have fun!",
             "\n",
             "Would you like to see this message on startup again?")
    writeLines(txt)
    choice <- menu(c("Please, no!", "Yes, please!"))
    if (choice == 2) return(TRUE)
    file.create(sup)
    return(TRUE)
}
