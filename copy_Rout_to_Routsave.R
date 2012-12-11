################################################################################
#
# USE THIS SCRIPT TO REPLACE ALL OCCURENCES OF *.Rout.save WITH NEW *.Rout FILES
#
# Author: Benjamin Hofner, 2012
#
# USAGE:
#   Use either
#     R CMD BATCH "--args which='mboostDevel'" copy_Rout_to_Routsave.R
#   or
#     R CMD BATCH "--args which='mboostPatch'" copy_Rout_to_Routsave.R
#
################################################################################

## Get command line arguments
args <- commandArgs(TRUE)
if (length(args) > 1)
    stop("specify (at maximum) one argument (i.e., which)")
eval(parse(text=args))
if (length(args) == 0)
    which <- "mboostDevel"

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
           "#   svn revert --recursive pkg/mboostDevel/tests",
           "#   svn revert --recursive pkg/mboostPatch/tests"),
    "#########################################################################",
    sep = "\n")

################################################################################
## Copy output of vignettes

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
                  res <- grep(paste(fn, "\\.Rout\\.save", sep=""), vROUT.SAVE))

    vROUT.SAVE <- vROUT.SAVE[idx]

    cbind(vROUT2, vROUT.SAVE)

    cat("\n\nCopy *.Rnw.log to *.Rout.save:\n---------------------------\n")

    for (i in 1:length(vROUT))
        print(file.copy(vROUT2[i], vROUT.SAVE[i], overwrite = TRUE))

    cat("#########################################################################",
    "# To revert changes simply use:",
    ifelse(which == "mboostDevel",
           "#   svn revert --recursive pkg/mboostDevel/vignettes",
           "#   svn revert --recursive pkg/mboostPatch/vignettes"),
    "#########################################################################",
    sep = "\n")
} else {
    cat("\n\nNOTE: No changes in output of vignettes.\n\n")
}
