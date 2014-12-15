################################################################################
#
# USE THIS SCRIPT TO REPLACE ALL OCCURENCES OF *.Rout.save WITH NEW *.Rout FILES
#
# Author: Benjamin Hofner, 2012
#
# USAGE:
#   Use either
#     ## To copy test output
#     Rscript copy_Rout_to_Routsave.R "which='mboostDevel'" "vignettes=FALSE"
#     ## To copy vignette output
#     Rscript copy_Rout_to_Routsave.R "which='mboostDevel'" "vignettes=TRUE"
#   or
#     ## To copy test output
#     Rscript copy_Rout_to_Routsave.R "which='mboostPatch'" "vignettes=FALSE"
#     ## To copy vignette output
#     Rscript copy_Rout_to_Routsave.R "which='mboostPatch'" "vignettes=TRUE"
#
#
# ALTERNATIVE USAGE (with R CMD BATCH):
#   Use either
#     ## To copy test output
#     R CMD BATCH "--args which='mboostDevel' vignettes=FALSE" copy_Rout_to_Routsave.R
#     ## To copy vignette output
#     R CMD BATCH "--args which='mboostDevel' vignettes=TRUE" copy_Rout_to_Routsave.R
#   or
#     ## To copy test output
#     R CMD BATCH "--args which='mboostPatch' vignettes=FALSE" copy_Rout_to_Routsave.R
#     ## To copy vignette output
#     R CMD BATCH "--args which='mboostPatch' vignettes=TRUE" copy_Rout_to_Routsave.R
#
################################################################################

## Get command line arguments
args <- commandArgs(TRUE)
if (length(args) > 2)
    stop("specify (at maximum) two arguments (i.e., which and vignettes)")
eval(parse(text=args))
if (length(args) == 0) {
    which <- "mboostDevel"
    vignettes <- FALSE
}
if (is.null(which))
    which <- "mboostDevel"
if (is.null(vignettes))
    vignettes <- FALSE

if (which == "mboostDevel") {
    path <- "pkg/mboostDevel"
    check_path <- "mboostDevel.Rcheck/"
} else {
    if (which == "mboostPatch") {
        path <- "pkg/mboostPatch"
        check_path <- "mboost.Rcheck/"
    } else {
        stop("which is not correctly specified")
    }
}

################################################################################
## Copy output of test files

if (vignettes == FALSE) {

    ## Get relevant file names
    ROUT <- list.files(path = check_path, pattern = ".Rout$", recursive = TRUE)
    ROUT2 <- paste(check_path, ROUT, sep ="")

    ROUT.SAVE <- list.files(path = path, pattern = ".Rout.save$", recursive = TRUE)
    ROUT.SAVE <- paste(path, "/", ROUT.SAVE, sep ="")
    ROUT.SAVE <- ROUT.SAVE[grep("test", ROUT.SAVE)]

    ## sort ROUT.SAVE
    idx <- rep(NA, length(ROUT))
    for (i in 1:length(ROUT))
        idx[i] <- grep(ROUT[i], ROUT.SAVE)
    ROUT.SAVE <- ROUT.SAVE[idx]

    cbind(ROUT2, ROUT.SAVE)

    cat("\n\nCopy *.Rout to *.Rout.save:\n---------------------------\n")

    for (i in 1:length(ROUT))
        print(file.copy(ROUT2[i], ROUT.SAVE[i], overwrite = TRUE))

    cat("#########################################################################",
        "# To revert changes simply use:",
        ifelse(which == "mboostDevel",
               "#   svn revert --recursive pkg/mboostDevel/tests\n# or use\n#   git checkout -- pkg/mboostDevel/vignettes",
               "#   svn revert --recursive pkg/mboostPatch/tests\n# or use\n#   git checkout -- pkg/mboostPatch/vignettes"),
        "#########################################################################",
        sep = "\n")

}

################################################################################
## Copy output of vignettes

if (vignettes == TRUE) {
    vpath <- paste(path, "vignettes", sep ="/")

    ## get vignette output as created by R CMD check:
    vROUT <- list.files(path = check_path, pattern = ".Rnw.log$")
    if (length(vROUT) > 0) {
        vROUT2 <- paste(check_path, vROUT, sep ="")

        vROUT.SAVE <- list.files(path = vpath, pattern = ".Rout.save",
                                 recursive = TRUE)
        vROUT.SAVE <- paste(vpath, vROUT.SAVE, sep = "/")

        ## sort
        filenames <- gsub("(.*)\\.Rnw\\.log", "\\1", vROUT)
        idx <- sapply(filenames, function(fn)
                      res <- grep(paste(fn, "\\.Rout\\.save$", sep=""), vROUT.SAVE))

        vROUT.SAVE <- vROUT.SAVE[idx]

        cbind(vROUT2, vROUT.SAVE)

        cat("\n\nCopy *.Rnw.log to *.Rout.save:\n---------------------------\n")

        for (i in 1:length(vROUT))
            print(file.copy(vROUT2[i], vROUT.SAVE[i], overwrite = TRUE))

        cat("#########################################################################",
            "# To revert changes simply use:",
            ifelse(which == "mboostDevel",
                   "#   svn revert --recursive pkg/mboostDevel/vignettes\n# or use\n#   git checkout -- pkg/mboostDevel/vignettes",
                   "#   svn revert --recursive pkg/mboostPatch/vignettes\n# or use\n#   git checkout -- pkg/mboostPatch/vignettes"),
            "#########################################################################",
            sep = "\n")
    } else {
        cat("\n\nNOTE: No changes in output of vignettes.\n\n")
    }
}
