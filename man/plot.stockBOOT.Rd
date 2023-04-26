\name{plot.stockBOOT.object}
\alias{plot.stockBOOT.object}
\title{Plots results of a stockR analysis (a barplot-like representation)}
\description{For a given fitted stockBOOT object, visualise the membership of each individual to each stock.}
\section{method}{
  plot.stockBOOT.object( x, locations=NULL, plotTitle=NULL, CI.width=0.95, region.lwd=3.5, ...)
}
\arguments{
  \item{x}{A stockBOOT.object, typically obtained from a call to \code{\link{stockBOOT}}. Note that the accuracy of the results will depend on the number of bootstrap resamples used to create the stockBOOT.object.}
  \item{locations}{A data.frame with two columns and the number of rows defined by the number of individuals in the stockBOOT.object.  Important to note that the ordering of individuals must be the same in the locations data.frame as in the stockBOOT.oject (otherwise garbage will be given). The first column in the locations data.frame is the 'region' from which an individual was sampled.  The second column gives information about how to plot the 'regions' -- in particular their plotting order.  Assign a low number to regions that you want plotted on the left side of the page and a high number for those on the right. Default argument value (for locations) is NULL, in which case all individuals are assumed to come from a single region (plot is for a single block).}
  \item{plotTitle}{The main title (text) to give the plot.}
  \item{CI.width}{The width of the confidence intervals to take transparency from. See details.}
  \item{region.lwd}{The line width of the box around the regions (groups of individuals).}
  \item{...}{Other parameters to be passed through to plotting functions.}
}
\value{
  NULL
}
\details{
These plots give the probability of an individual to be assigned to each group (stock) identified by the mixture model, through stockSTRUCTURE and stockBOOT. This is an assignment probability, not a admixture proportion that is obtained from the STRUCTURE program (for example). 

The intensity of the colour is given by the amount of uncertainty in the estimate of the probabilities -- more solid colours are less uncertain (or more certain if you like to avoid double negatives). Very faint colours have a 100*CI.wdith\% confidence interval width that is essentially 1 (so nothing is known about the probabilities).
  }
  \author{Scott D. Foster}
  \references{
    Foster, S.D., Feutry, P., Grewe, P.M., Berry, O, Hui, F.K.C. and Davies, C.R. (in press) Reliably Discriminating Stock Structure with Genetic Markers: Mixture Models with Robust and Fast Computation. Molecular Ecology Resources
  }
  \seealso{\code{\link{stockSTRUCTURE}}, \code{\link{stockBOOT}}}
  
  \examples{
    \donttest{
      ##This example will take a little while to run.
      #This should be challenging as there are actually 2 stocks and we fit a model with 3.
      tmpDat1 <- sim.stock.data( nAnimals=100, nSNP=5000, nSampleGrps=100, K=2, ninform=5000, 
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
      #Let's visualise it.
      plot( tmpBOOT, locations=NULL, plotTitle="Data contains 2 groups, model fits 3")
      #You can try it with 2 groups.
      
      #Let's now pretend that there are sampling regions
      plot( tmpBOOT, locations=data.frame( locations=rep(1:4, each=25), 
	  order=rep( c(4,3,1,2), each=25)), plotTitle="Plot with grouping")
    }
  }
  
  
