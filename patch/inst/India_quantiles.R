
# Note: Requires Data File IAKR51FL.DTA 
# (Children`s recode for the 2005/2006 India DHS in Stata format)
# Freely available from www.measuredhs.com for research purposes
# (-> Data -> Available Datasets -> India)

# Preprocessing of the raw data
source("india_preproc.R", echo=TRUE)

# Fit quantile regressions
source("india_fit.R", echo=TRUE)
