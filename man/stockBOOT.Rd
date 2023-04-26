\name{stockBOOT}
\alias{stockBOOT}
\title{Calculates bootstrap estimates of stockSTRUCTURE model}
\description{For a given fitted stockSTRUCTURE model determine uncertainty in parameter values (and outputs) using Bayesian Bootstrap.}
\usage{
  stockBOOT( fm, B=100, mc.cores=NULL)
}
\arguments{
  \item{ fm}{A stockSTRUCTURE.object, typically obtained from a call to \code{\link{stockSTRUCTURE}}. This object MUST NOT be run with the small object option. So, in the initial model fit make sure that small.object=FALSE in the control list (this is the default).}
  \item{ B}{The number of (Bayesian) bootstrap resamples to perform. Default is 100 but for serious applications more will be needed.}
  \item{ mc.cores}{The number of parallel processes to perform at once. Default is NULL, in which case the function will determine the number of cores available and use them all.}
}
\value{
  An object of class stockBOOT.object. This is a list with the following components
  \item{ condMeans}{a three dimensional array with each slice corresponds to the mean frequencies for each allele in each stock for each bootstrap resample.}
  \item{ pis}{a three dimensional array where each slice corresponds to the proportion, in the data, of each of the stocks for each bootstrap resample. Small numbers imply smaller stocks. Note that sum(pi)=1 within each resample.}
  \item{ postProbs}{a three dimensional array where each slice is the soft-classification of sample.grps to stocks for each bootstrap resample. Be aware that the order will be according to levels( sample.grps) and not the \code{natural} ordering that you may have imagined. Blame this on R{}'s internal representation of factors. Consider using mixedsort() from the gtools package.}
  \item{ margLogl}{a vector giving the marginal log-likelihood obtained at the parameter estimates for each bootstrap sample.}
  }
\details{
Bootstrapping is done here using the Bayesian Bootstrap of Rubin (1981) as described in Foster et al (2018). Briefly, the sampling distribution of the model's parameters is obtained by taking weighted re-samples of the original data and then re-estimating the model. The variability between the re-estimations is the same as the sampling variability.
}
\author{Scott D. Foster}
\references{
Foster, S.D., Feutry, P., Grewe, P.M., Berry, O, Hui, F.K.C. and Davies, C.R. (in press) Reliably Discriminating Stock Structure with Genetic Markers: Mixture Models with Robust and Fast Computation. Molecular Ecology Resources

Rubin, D.B. (1981) The Bayesian Bootstrap. The Annals of Statistics \emph{9}:130--134.
}

\seealso{\code{\link{stockSTRUCTURE}}}

\examples{
\donttest{
  ##This example will take a little while to run.
  #This should be very challenging as stock differentiation is non-existant (K=1).
  tmpDat1 <- sim.stock.data( nAnimals=100, nSNP=5000, nSampleGrps=100, K=1, ninform=5000, 
                                    sds=c(alpha=1.6,beta.inform=0.75,beta.noise=0.0005))
  #EM estimation from Kmeans starting values
  tmp <- stockSTRUCTURE( tmpDat1, sample.grps=attr(tmpDat1,"sampleGrps"), K=3, start.grps=NULL)
  #in general, you'll want to use as many cores as possible, or close to.
  #mc.cores=1 is used here to please the CRAN submission checks
  tmpBOOT <- stockBOOT( tmp, B=100, mc.cores=1) 
  print( round( apply( tmpBOOT$postProbs, FUN=quantile, MARGIN=1:2, probs=c(0.025,0.975))), 5)
  #Note that, in this case, the posterior probabilities are not very informative; they could
  #be anywhere between 0 and 1. There are likely to be a few individuals, of course, where 
  #they have a very low chance of belonging to a particular stock (and this is chance). There
  #may even some individuals that get assigned to a group with almost certainty.
  }
}

