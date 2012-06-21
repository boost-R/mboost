######
# USE THIS SCRIPT TO FIND AND REPLACE ALL OCCURENCES OF mboostDevel with mboost
#
# Author: Benjamin Hofner, 2012

FILES <- list.files(path = "pkg/mboostDevel", recursive = TRUE)
FILES <- paste("pkg/mboostDevel/", FILES, sep ="")

cat("\n\nChanging mboostDevel to mboost:\n-------------------------------\n")

OUT <- lapply(FILES, function(file) {
    a <- readLines(con = file)
    if (any(grepl("mboostDevel", a))){
        b <- gsub("mboostDevel", "mboost", a)
        cat("Write ", file, "\n")
        writeLines(b, con = file)
    }
})

## ATTENTION: DO NOT COMMIT THE NEW, MODIFIED FILES TO mboostDevel
## (only to mboostPatch and CRAN)

## To revert changes simply use:
## svn revert --recursive pkg/mboostDevel
