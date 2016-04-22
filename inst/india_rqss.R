#-----------------------------------------------------------
# Path, Packages, source functions
#-----------------------------------------------------------

library("quantreg")
load("india.Rdata")
source("india_rqss_lambdaOptFunc.R")

#-----------------------------------------------------------
# Quantiles and Indices
#-----------------------------------------------------------

quantiles <- c(0.05, 0.1, 0.5)
inds <- paste("w", 1:50, sep="")


for(i in 1:length(inds)){

  indiaTrain <- india[india[,inds[i]]==1,]
  
  indiaTuning <- india[india[,inds[i]]==2,] 
  indiaTuning <- convexHull(dataOrg=indiaTuning, dataRef=indiaTrain)
  
  indiaTest <- india[india[,inds[i]]==3,] 
  indiaTest <- convexHull(dataOrg=indiaTest, dataRef=indiaTrain)
  
    
  for(j in 1:length(quantiles)){
  
  
       # Starting parameters for lambda equal to 100:
       #optimLambda <- optim(rep(100,6), fn=lambdaOpt, fitdata=indiaTrain, 
       #                     evaldata=indiaTuning, tau=quantiles[j], lower=0)                     
       
       # Next try: Starting parameters equal to 10:
       optimLambda <- optim(rep(10,6), fn=lambdaOpt, fitdata=indiaTrain, 
                            evaldata=indiaTuning, tau=quantiles[j], lower=0) 
       lambdaMin <- optimLambda$par
       
       formule <- buildFormula(lambdas=lambdaMin)
       model <- rqss(formule, data=indiaTrain, tau=quantiles[j])

       # Results with starting parameters for lambda equal to 100:
       save(model, file=paste("rqssObjects/model",i,"_",quantiles[j],".Rdata",sep=""))
       rm(model)
       gc()

  }
}



