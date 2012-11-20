######
# USE THIS SCRIPT TO REPLACE ALL OCCURENCES OF *.Rout.save WITH NEW *.Rout FILES
#
# Author: Benjamin Hofner, 2012


################################################################################
# Get command line arguments specified via:
#   R CMD BATCH "--args which='mboostDevel'" copy_Rout_to_Routsave.R
# or
#   R CMD BATCH "--args which='mboostPatch'" copy_Rout_to_Routsave.R
args <- commandArgs(TRUE)
if (length(args) > 1)
    stop("specify (at maximum) one argument (i.e., which)")
eval(parse(text=args[[i]]))
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
