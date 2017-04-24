###########################################################
##
## This code replicates the analysis presented in 
## 
##   N. Fenske, T. Kneib, and T. Hothorn (2011), 
##   Identifying risk factors for severe childhood malnutrition by 
##   boosting additive quantile regression.
##   Journal of the American Statistical Association, 106:494-510.
## 
###########################################################
## 
## Note: Requires Data File IAKR51FL.DTA 
## (Children`s recode for the 2005/2006 India DHS in Stata format)
## Freely available from www.measuredhs.com for research purposes
## (-> Data -> Available Datasets -> India)
##
###########################################################

# Preprocessing of the raw data
source("india_preproc.R", echo=TRUE)

# Fit additive quantile regressions
source("india_additive.R", echo=TRUE)

# Fit VCM quantile regressions
source("india_vcm.R", echo=TRUE)

# Fit quantile regressions trees
source("india_blackboost.R", echo=TRUE)

# Fit additive quantile regression stumps
source("india_stumps.R", echo=TRUE)

# Summarise results for plotting
# (all algorithms without rqss)
source("india_summary.R", echo=TRUE)
# and create plots
source("india_plots.R", echo=TRUE)


# Fit additive quantile regression by rqss
source("india_rqss.R", echo=TRUE)

# And plot results for empirical risk results for all algorithms
source("india_rqssResults.R", echo=TRUE)