\name{stockSTRUCTURE}
\alias{stockSTRUCTURE}
\title{Finds stock structure for sampling groups of animals}
\description{For a given set of markers, scored for each animal (in sampling groups), determine likely stock structure using mixture models.
}
\usage{
stockSTRUCTURE(SNPdata = NULL, sample.grps = as.factor(1:ncol(SNPdata)),  K = 3, 
        weights = rep(1, nlevels( sample.grps)), start.grps = NULL,  control = NULL)
}
\arguments{
\item{ SNPdata}{a numeric matrix of dimension (number of SNPs X number of individuals). As the dimension suggests, each row corresponds to an SNP marker and each column an individual. The entries are the number of copies of a marker present in the animal. So, 0 for no copies (aa), 1 for one copy (Aa), and 2 for two copies (AA). Columns that are constant for all animals are removed.}
\item{ sample.grps}{a factor (or something that can be converted into a factor) of length ncol( SNPdata). This gives the sampling the allocation of animals (columns of SNPdata) to sampling groups. The animals must be ordered the same in sample.grps and in SNPdata. The default sample.grps is to have each animal in its own group (so that no knowledge of breeding groups is assumed). Note that if sample.grps is specified, then the function's output will have (internal) ordering which is based on the alphabetical listing of sample.grps.}
\item{ K}{an integer giving the number of stocks to partition the data into}
\item{ weights}{a numeric vector of length ncol( SNPdata). Gives the weighting for each animal to the analysis. Default is to weight all animals equally. You should have a really good reason not to use equal weights.}
\item{ start.grps}{a numeric vector of length ncol( SNPdata) or NULL. If NULL the (EM-)optimisation algorithm will be started from groups found using the K-means algorithm. If not NULL then these groups will be used to start the EM optimiser. For most purposes this argument should be set as NULL. If specified, then the model's likelihood is less likely to be maximised and the major signals in the data may not be found. This is not a prior, not in the sense of a Bayesian analysis in any case.  Please do not supply the groups you think are in the data as starting values -- you are likely to just be confirming your own beliefs.}
\item{ control}{a list of control parameters for initialisation and optimisation. See below for details.}
}
\value{
An object of class stockSTRUCTURE.object. This is a list with the following components
\item{ condMeans}{the mean frequencies for each allele in each stock}
\item{ pis}{the proportion, in the data, of each of the stocks. Small numbers imply smaller stocks. Note that sum(pi)=1.}
\item{ margLogl}{the marginal log-likelihood obtained at the parameter estimates.}
\item{ postProbs}{the soft-classification of sample.grps to stocks. Be aware that the order will be according to levels( sample.grps) and not the \code{natural} ordering that you may have imagined. Blame this on R{}'s internal representation of factors. Consider using mixedsort() from the gtools package.}
\item{ K}{the number of stocks (specified by user).}
\item{ data}{a list containing all the bits and pieces of data used in the analysis.}
\item{ init.grps}{the initial stock structure for starting the final optimisation (via the EM algorithm).}
\item{ constantSNPs}{a boolean vector giving the locations of the SNP markers with no variation.}
\item{ margLoglTrace}{a trace of the marginal log-likelihood throughout the final optimisation (EM algorithm).}
\item{ control}{the control parameters used in estimation. See below for details.}
Note that only condMeans, pis, margLogl, postProbs, and K are returned if control$small.object=TRUE is specified. See below for details.
Note that if sample.grps is specified, then the value of the function's output is orientated for the alphabetical ordering of sample.grps. This is done to help keep track of the groupped results when the number of sample groups is less than the number of individuals.
}
\details{
Control arguments, and especially their defaults, will vary depending on the type of optimisation used (only a handful are common). The arguments are

\bold{Common to all methods}

\describe{
\item{small.object}{a boolean indicating whether the return object should be made small (ie no returning of data etc). Default is FALSE, so that the return object is normal size.}
\item{quiet}{a boolean indicating whether reporting should be performed throughout the estimation process. Default is FALSE (not quiet) so that reporting is performed. Users could also use \code{suppressMessages} if they prefer.}
\item{method}{a character string. One of "Kmeans.EM" (default), "SA.EM", "EM", or "DA.EM". These specify how the optimisation is performed. "SA.EM" implements the initialisation-then-optimisation strategy in Foster et al (in prep). "EM" performs the EM algorithm from random starts (unless user specifies the start location), this is similar to Chen et al (2006). "DA.EM" implements the deterministic annealing algorithm proposed in Zhou and Lange (2010) from random starts, or user specified start.}
}

*For "Kmeans.EM" (Default and recommended) method*

\describe{
\item{nKmeanStart}{the number of K-menas groupings to perform to obtain starting values. Default (NULL) corresponds to 25 starts.}
\item{nKPCA}{the SNP data is rotated, prior to initial K-means clustering, using PCA. This argument defines the number of PCAs to use. The default is nKPCA=min( nFish, nSNP, 100). The argument is always checked and the value of min( nFish, nSNP, nKPCA) is used.}
\item{EM.eps}{the absolute convergence tolerance for the EM-algorithm. Default is 1e-5, that is successive differences in parameters have to be very small before converged is reached.}
\item{EM.maxit}{the maximum number of iterations for the EM-algorithm. Default is 100.}
\item{EM.minit}{the minimum number of iteraction for the EM-algorithm. Default is 1. EM.minit is really only important for when EM.tau.step < 1.}
\item{EM.tau.step}{the step-size of the initial update for the posterior probabilities. Default is 1. If less than 1 (say 0.5), then the step size is less (half say) than the size that the EM-algorithm suggests. EM.tau.step increases linearly to 1 after EM.minit steps.}
\item{tau.init.max}{After the initial hard-clustering, the groupings are made soft by specifying posterior probabilities that are less than 1 and greater than 0. The default is given in Foster et al (in prep), but makes sure that the hard clustering recieves twice the probability mass than other stocks.}
}

*For "SA.EM" method*

\describe{
\item{SANN.temp}{the initial temperature for simulated annealing. Default is half the number of sampling groups. This means that half of the sampling groups can be swapped to new groups in the first iteration.}
\item{SANN.maxit}{the maximum number of iterations for the simulated annealing. Default is 5000, which is probably overkill for many problems.}
\item{SANN.tmin}{the minimum number of swaps per iteration (ie once the annealing has run for a while). The default is max( nSampGrps \%/\% 20, 1) where nSampGrps is the number of sampling groups.}
\item{SANN.nreport}{the number of iterations to do before printing out report (if printing at all). Default is 100.}
\item{SANN.trace}{a boolean indicating if any trace information should be given. See \code{\link{optim}}.}
\item{EM.eps}{the absolute convergence tolerance for the EM-algorithm. Default is 1e-5, that is successive differences in parameters have to be very small before converged is reached.}
\item{EM.maxit}{the maximum number of iterations for the EM-algorithm. Default is 100.}
\item{EM.minit}{the minimum number of iteraction for the EM-algorithm. Default is 1. EM.minit is really only important for when EM.tau.step < 1.}
\item{EM.tau.step}{the step-size of the initial update for the posterior probabilities. Default is 1. If less than 1 (say 0.5), then the step size is less (half say) than the size that the EM-algorithm suggests. EM.tau.step increases linearly to 1 after EM.minit steps.}
\item{tau.init.max}{After the initial hard-clustering, the groupings are made soft by specifying posterior probabilities that are less than 1 and greater than 0. The default is given in Foster et al (in prep), but makes sure that the hard clustering recieves twice the probability mass than other stocks.}
}

*For "EM" method*

 As per last five entries for "SA.EM" (an having a random start rather than an initial clustering).

*For "DA.EM" method*

\describe{
\item{EM.eps}{the absolute convergence tolerance for the EM-algorithm. Default is 1e-5, that is successive differences in parameters have to be very small before converged is reached.}
\item{EM.maxit}{the maximum number of iterations for the EM-algorithm. Default is 100.}
\item{EM.minit}{the minimum number of iteraction for the EM-algorithm. Default is 25. EM.minit is the number of steps required to reach the non-cooled log-likelihood surface.}
\item{EM.tau.step}{the step-size of the initial update for the posterior probabilities. Default is 1. If less than 1 (say 0.5), then the step size is less (half say) than the size that the EM-algorithm suggests. EM.tau.step increases linearly to 1 after EM.minit steps.}
\item{tau.init.max}{After the initial hard-clustering, the groupings are made soft by specifying posterior probabilities that are less than 1 and greater than 0. The default is to make the hard clusterings ever-so-slightly more likely than others.}
\item{DA.nu.init}{the initial cooling parameter. Default is 1e-4. Algorithm seems particularly sensitive to this choice. The cooling parameter nu is increased gradually until it reaches 1 (max) at EM.tau.step iterations. See Zhou and Lange (2010) for details.}
}
}
\seealso{\code{\link{sim.stock.data}}, \code{\link{stockBOOT}}}
\author{Scott D. Foster}
\references{
Foster, S.D., Feutry, P., Grewe, P.M., Berry, O, Hui, F.K.C. and Davies, C.R. (in press) Reliably Discriminating Stock Structure with Genetic Markers: Mixture Models with Robust and Fast Computation. Molecular Ecology Resources
}
\examples{
\donttest{
set.seed(727)
tmpDat1 <- sim.stock.data( nAnimals=100, nSNP=5000, nSampleGrps=100, K=3, ninform=5000, 
                                    sds=c(alpha=1.6,beta.inform=0.75,beta.noise=0.0005))
#This should not be too challenging as stock differentiation is quite large.
print( calcFst( tmpDat1, attributes( tmpDat1)$grps))
#EM estimation from Kmeans starting values
tmp <- stockSTRUCTURE( tmpDat1, sample.grps=attr(tmpDat1,"sampleGrps"), K=3, start.grps=NULL)
#an easy gold-standard for simulations (only know starting values as this is simulated data)
tmp1 <- stockSTRUCTURE( tmpDat1, sample.grps=attr(tmpDat1,"sampleGrps"), K=3, 
                      start.grps=attr( tmpDat1,"grps"), control=list( method="EM"))
#an easily misled estimation method
tmp2 <- stockSTRUCTURE( tmpDat1, sample.grps=attr(tmpDat1,"sampleGrps"), K=3, 
                                              start.grps=NULL, control=list( method="EM"))
#combine into a single object to investigate results
grpings <- cbind( attributes( tmpDat1)$grps, apply( tmp$postProbs, 1, which.max), 
              apply( tmp2$postProbs, 1, which.max), apply( tmp1$postProbs, 1, which.max))
colnames( grpings) <- c("True","Kmeans.EM Estimated","EM Estimated","Estimated From TRUE Start")
print( "How did we go with the stock identification?")
print( grpings)
#up to label switching this looks good, except for the EM from random start method (a few 
#confused individuals for the second and thrid stocks)
}
}
