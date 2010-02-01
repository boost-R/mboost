\name{birds}
\alias{birds}
\docType{data}
\title{ Habitat Suitability for Breeding Bird Communities }
\description{
  Environmental variables and bird counts for identifying
  suitable bird habitats
}
\usage{data("birds")}
\format{
  A data frame with 258 observations on the following 10 variables.
  \describe{
    \item{\code{GST}}{Growing stock per grid}
    \item{\code{DBH}}{Mean diameter of the largest three trees}
    \item{\code{AOT}}{Age of oldest tree}
    \item{\code{AFS}}{Age of forest stand}
    \item{\code{DWC}}{Amount of dead wood of conifers}
    \item{\code{LOG}}{Amount of logs per grid}
    \item{\code{x_gk}}{grid location, x coordinate}
    \item{\code{y_gk}}{grid location, y coordinate}
    \item{\code{SG4}}{observed number of birds from structural gild 4: 
        Requirement of regeneration (Phylloscopus trochilus, Aegithalos caudatus)}
    \item{\code{SG5}}{observed number of birds from structural gild 5:
        Requirement of regeneration combined with planted conifers (Phylloscopus collybita, Turdus
        merula, Sylvia atricapilla).}
  }
}
\details{

Counts of breeding bird
communities collected at 258 observation plots in a northern
Bavarian forest district are the response variable of interest.  
Along with the number of birds in two structural gilds, 
6 covariates are given here and
one is interested in quantifying their impact on habitat
suitability.

}
\source{

 Joerg Mueller (2005).
 Forest structures as key factor for beetle and bird communities 
 in beech forests. PhD thesis, Munich University of Technology.
 \url{http://mediatum.ub.tum.de}
 
}
\references{

  Thomas Kneib and Joerg Mueller and Torsten Hothorn (2008),
  Spatial smoothing techniques for the assessment of habitat suitability,
  \emph{Environmental and Ecological Statistics}, 
  \bold{15}(3), 343--364.

}
\keyword{datasets}
