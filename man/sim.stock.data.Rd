\name{sim.stock.data}
\alias{sim.stock.data}
\title{Simulates SNP data according to stockSTRUCTURE model}
\description{For given sample characteristics (number of sampling groups, number of animals, number of SNPs and number/size of stocks) simulates some SNP data.
}
\usage{
sim.stock.data(nAnimals, nSNPs, nSampleGrps, K, ninform = nSNPs, 
  means=c(alpha=-1.9,beta=0), sds=c(alpha=1.6,beta.inform=0.85,
  beta.noise=0.0005), minStockSize=1)
}
\arguments{
\item{ nAnimals}{integer giving the number of animals to simulate. No default, must be specified.}
\item{ nSNPs}{integer giving the number of SNPs to generate. No default, must be specified.}
\item{ nSampleGrps}{integer giving the number of sampling groups. Must be less than nAnimals. No default, must be specified.}
\item{ K}{integer giving the number of stocks}
\item{ ninform}{integer giving the number of informative SNPs. That is the first ninform SNPs will discriminate the stocks and the remainder will not.}
\item{ means}{a named two-element numeric vector consisting of the mean of the intercepts (named "alpha") and the mean of the deviations (named "beta"). The latter should always be zero (the betas will get standardised in any case).}
\item{ sds}{a names three-element numeric vector consisting of the standard deviations of the intercepts (named "alpha"), the sds of the betas that are differentiating between stocks (named "beta.inform"), and the sds of the SNP effects that are not differentiating between stocks (named "beta.noise").}
\item{ minStockSize}{the size of the smallest stock group allowed. That is the size within the sample, not the actual population size.}
}
\value{A numeric matrix containing SNP allele frequencies (in rows) for each animal (in columns). The matrix has attributes, "grps" for the stock groups and "sampleGrps" for the sampling groups.
}
\details{The proportions of each stock in the sample is taken as a single draw from a flat Dirichlet distribution (which is flat over the simplex). The expected frequencies, in each stock, are generated on the logit scale. The mean frequency, for each SNP, is drawn from a random normal with mean means["alpha"] and sd sds["alpha"]. Deviations from the mean, for the first ninform SNPs, are generated from a normal with mean means["beta"] and sd=sds["beta.inform"]. For other SNPs these deviations are drawn from a normal with mean means["alpha"] and sd sds["beta.noise"] (typically really small). These expectations are then inverse logit trasnsformed and SNP data generated using binomial sampling with 2 trials. There was a good reason for choosing these means and standard deviations for generating SNP expectations, but I can't remember what it was...
}
\seealso{\code{\link{stockSTRUCTURE}}}
\author{Scott D. Foster}
\examples{
#a data set with 100 individuals, from 100 different sampling groups (so no individual
#can be assumed a priori to belong to the same stock as any other individual), with 
#5000 SNP markers. There are 3 stocks and only 200 of the 500 SNP markers are 
#differentiated between stocks.
myData <- sim.stock.data( nAnimal=100, nSNP=5000, nSampleGrps=100, K=3, ninform=200)
print( dim( myData))  #should be 5000 x 100
}
