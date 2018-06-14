\name{calcFst}
\alias{calcFst}
\title{Calculate Fst statistics for frequency data}
\description{Calculates Weir and Cockerham's Fst for frequency data
}
\usage{
calcFst(data, grps)
}
\arguments{
You can put normal text in \bold{Arguments}, too, like this. Remember to indent all arguments, as below.
\item{ data}{a nMarker by nAnimal matrix of allele frequencies. That is for animal i and codominant marker j data[j,i] is the number of copies of the SNP allele (0, 1, or 2).}
\item{ grps}{a numeric vector giving the populations of the nAnimals}
}
\value{A numeric matrix containing the pairwise Fst values between all populations.
}
\details{This function is really just a wrapper for the \code{Fstat()} function, taken from the (now unsupported) Geneland package.
}
\section{Authors}{Arnaud Estoup for original code in Turbo Pascal. Translation in Fortran and interface with R{} by Gilles Guillot (for Geneland package). Scott Foster for latest wrapper for stockR package.
}
\examples{
set.seed(717)
data <- sim.stock.data( nAnimals=50, nSNPs=500, nSampleGrps=50, K=5, ninform=50)
calcFst( data, attributes( data)$grps)
}
\keyword{misc}
\author{Scott D. Foster}
