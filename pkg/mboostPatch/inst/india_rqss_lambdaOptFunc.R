#============================================================================
# Function lambdaOpt()
#============================================================================
# Aim: Finding an optimal lambda for rqss (this function is minimized by optim or optimize)
#
# Input:
# --------
# lambdas     Vector of lambda values
# modeldata   Training data from generateData()
# evaldata    Test data lying in the convex hull of modeldata
# tau         Quantile of interest
# --------
# Output:
# --------
# Empirical risk when rqss model is estimated on modeldata and
# evaluated on evaldata
#

lambdaOpt <- function(lambdas, fitdata, evaldata, tau){

          # Model formula
          form <- buildFormula(lambdas=lambdas)

          # Generate rqss object
          rqssLambda <- rqss(form, data=fitdata, tau=tau)

          # Prediction with test data  
          predLambda <- predict(rqssLambda, newdata=evaldata)

          # Calculate empirical risk
          CVrisk <-loss(evaldata$stunting, predLambda, tau=tau)

          return(CVrisk)

}

#==============================================================================
# Check function
#==============================================================================
rho <- function(u, tau){u*(tau - (u<0))}

#==============================================================================
# Risk
#==============================================================================
loss <- function(yTrue,yHat,tau){sum(rho(yTrue-yHat,tau))/length(yTrue)}
  
#======================================================================
# Function buildFormula()
#======================================================================
# Aim: Building formula for function rqss()

buildFormula <- function(lambdas){

      return(as.formula(paste("stunting ~ qss(cage,lambda =", lambdas[1],") + qss(breastfeeding, lambda =", lambdas[2],
      ") + qss(mbmi, lambda =", lambdas[3],") + qss(mage, lambda =", lambdas[4],") + qss(medu, lambda =", lambdas[5],
      ") + qss(edupartner, lambda =", lambdas[6],") + csex + ctwin + cbirthorder + munemployed + mreligion + mresidence
      + deadchildren + wealth + electricity + radio + television + refrigerator + bicycle + motorcycle + car", sep="")))
}



#============================================================================
# Function convexHull()
#============================================================================
# Aim: Reduce data set such as lying in the convex hull of reference data
#
# Input:
# --------
# dataOrg     Origninal data, here: valdata
# dataRef     Reference data, here: traindata
# --------
# Output:
# --------
# reduced data set where each covariate lies in the range of the
# covariates from modeldata
#

convexHull <- function(dataOrg, dataRef){

        conNames <- c("cage", "breastfeeding", "mbmi", "mage", "medu", "edupartner")
        indMat <- matrix(0, nrow=nrow(dataOrg), ncol=6)

        for(i in 1:length(conNames)){
              indMat[,i] <- as.numeric(dataOrg[,conNames[i]] <= min(dataRef[,conNames[i]]) | dataOrg[,conNames[i]] >= max(dataRef[,conNames[i]]))
        }
        indVec <- apply(indMat, MARGIN=1, sum)

        # Observations with indVec>0 are removed from the dataOrg
        evaldata <- dataOrg[indVec==0,]

        return(evaldata=evaldata)

}


convexHullInd <- function(dataOrg, dataRef){

        conNames <- c("cage", "breastfeeding", "mbmi", "mage", "medu", "edupartner")
        indMat <- matrix(0, nrow=nrow(dataOrg), ncol=6)

        for(i in 1:length(conNames)){
              indMat[,i] <- as.numeric(dataOrg[,conNames[i]] <= min(dataRef[,conNames[i]]) | dataOrg[,conNames[i]] >= max(dataRef[,conNames[i]]))
        }
        indVec <- apply(indMat, MARGIN=1, sum)

        # Observations with indVec>0 are removed from the dataOrg
        #evaldata <- dataOrg[indVec==0,]

        return(indVec=indVec)

}



