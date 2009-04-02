
# Note: Requires Data File IAKR51FL.DTA
# Freely available from www.measuredhs.com for research purposes

# Preprocessing of the raw data
source("india_preproc.R", echo=TRUE)

# Fit quantile regressions
source("india_fit.R", echo=TRUE)
