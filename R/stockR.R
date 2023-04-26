# This is package stockR 

#this line of code is stop R CMD check complaining about no visible binding for global variable...
#I don't like global variables, but it makes simulated annealing easier (for now).
if(getRversion() >= "2.15.1")  utils::globalVariables(c( "pkg.globals", ".SANN.iter", ".SANN.maxit",".SANN.temp",".SANN.tmin"))

#This bit of code is for the stupid global variables that seemed to make life easier
#when implementing SANN.EM...
pkg.globals <- new.env()

pkg.globals$data_path <- "data"

#set_data_path <- function(path) {
#  pkg.globals$data_path <- path
#}


"bootFun" <-
function( x, fm1) 
{
  wts <- fm1$data$nSampGrps * my.rdirichlet( n=1, dim=fm1$data$nSampGrps)
  fm1$control$small.object <- TRUE
#  fm1$control$quiet <- TRUE
  fm.booty <- stockSTRUCTURE( SNPdata=fm1$data$SNPdata, sample.grps=fm1$data$sample.grps, K=fm1$K, weights=wts, start.grps=apply( fm1$postProbs, 1, which.max), control=fm1$control)
  fm.booty <- fm.booty[c("condMeans","pis","postProbs","margLogL")]
  
  return( fm.booty)  

}


"calcCompleteLogl_SA" <-
function( grps, all.data, K) 
{
  #grps is a nSampGrps lengthed vector
  taus <- model.matrix( ~-1+as.factor( grps))
  pis <- colMeans( taus)

  #get the cond means
  means <- calcCondMeans( all.data, taus, K)

  #the conditional logl contributions
  condLogls <- calcCondProbs_SA( all.data, grps, means, K)
  #combine into sampleGrp contributions
  condLoglsSampGrps <- tapply( condLogls, all.data$sample.grps, sum)

  #the contributions from logpi
  logPis <- log( pis)
  logPis <- logPis[grps]

  #the complete logl now
  unpenCLogl <- sum( all.data$weights*( condLoglsSampGrps + logPis))

  return( unpenCLogl)
}


"calcCondMeans" <-
function (all.data, taus, K) 
{
    taus.fish <- taus[all.data$sample.grps.num, , drop = FALSE]
    weights.fish <- all.data$weights[all.data$sample.grps.num]
    wted.taus.fish <- taus.fish * rep(weights.fish, times = K)
    #get the cond means
    means <- totWtsPerSNP <- matrix( NA, nrow=nrow( all.data$SNPdata), ncol=K)
    FishPerSNP <- apply( all.data$SNPdata, 1, function(x) which (!is.na( x)))
    if( is( FishPerSNP, "matrix"))
      FishPerSNP <- as.list( as.data.frame( FishPerSNP))
    for( kk in 1:K){
      tmp <- all.data$SNPdata * rep( wted.taus.fish[,kk], each=nrow( all.data$SNPdata))
      means[,kk] <- rowSums( tmp, na.rm=TRUE)
      totWtsPerSNP[,kk] <- sapply( FishPerSNP, function(x) sum( wted.taus.fish[x,kk]))
    }
    totWtsPerSNP[totWtsPerSNP==0] <- .Machine$double.xmin	#guard against zero wts (when group is empty)
    means <- means/ totWtsPerSNP
    means <- means/2
    means[means > 1] <- 1
    means[means < 0] <- 0
    return(means)
}


"calcCondProbs" <-
function( all.data, means, K) 
{
  condProbs <- matrix( NA, nrow=all.data$nFish, ncol=K)
  for( kk in 1:K){
    tmp.probs <- means[,kk]
    #A hack to stop alleles with exactly zero or one instances through 
    #replace exact zeros by a small value and exact ones by a small value.
    #shouldn't be a problem for standard data analysis, just cross-validation/Bootstrap
    #where fish are removed
    abit <- (1/all.data$nFish)/10
    tmp.probs[tmp.probs<abit] <- abit
    tmp.probs[tmp.probs>1-abit] <- 1-abit
    #now get the probs 'safely'
    tmp <- dbinom( as.matrix( all.data$SNPdata)[,drop=FALSE], size=2, prob=rep( tmp.probs, all.data$nFish), log=TRUE)
    condProbs[,kk] <- colSums( tmp, na.rm=TRUE)
  }
  return( condProbs)
}


"calcCondProbs_SA" <-
function( all.data, grps, means, K) 
{
  means[means==0] <- 1e-4 # a bit of a hack to prevent the boundary case
  means[means==1] <- 1-(1e-4)

  condProbs <- rep( NA, all.data$nFish)
  fishOrderGrps <- grps[all.data$sample.grps.num]
  for( kk in 1:K){
    tmp.probs <- means[,kk,drop=FALSE]
    fishGrps <- fishOrderGrps==kk
    tmp <- dbinom( as.matrix( all.data$SNPdata)[,fishGrps,drop=FALSE], size=2, prob=matrix( rep( tmp.probs, times=sum( fishGrps)), nrow=all.data$nSNP), log=TRUE)
    condProbs[fishGrps] <- colSums( tmp, na.rm=TRUE)
  }
  return( condProbs)
}


"calcFst" <-
function( data, grps){
  
  newDat <- matrix( NA, nrow=ncol( data), ncol=2*nrow( data))
  
  for( ii in 1:ncol( data)){
    for( jj in 1:nrow( data)){
      if( is.na( data[jj,ii]))
          newDat[ii,2*(jj-1)+1:2] <- -1
      else{
        if( data[jj,ii]==0)
          newDat[ii,2*(jj-1)+1:2] <- 0
        if( data[jj,ii]==2)
          newDat[ii,2*(jj-1)+1:2] <- 1
        if( data[jj,ii]==1){
          newDat[ii,2*(jj-1)+1] <- 1
          newDat[ii,2*(jj-1)+2] <- 0
        }
      }
    }
  }
  #Geneland is no longer supported for R, so I have taken code from Geneland and dropped it into stockR
  tmp <- Fstat( newDat, length( unique( grps)), pop.mbrship=grps)#attributes( tmpDat1)$grps)

  return( tmp$Fst)
}


"calcMargLogl" <-
function( new.parms, all.data, condProbs.fish, K) 
{
  condProbs.sampGrps <- matrix( NA, ncol=K, nrow=all.data$nSampGrps)
  for( jj in levels( all.data$sample.grps))
    condProbs.sampGrps[jj==levels( all.data$sample.grps),] <- colSums( condProbs.fish[all.data$sample.grps==jj,,drop=FALSE])
  logPis <- log( new.parms$pi)
  tmp <- condProbs.sampGrps + rep( logPis, each=all.data$nSampGrps)

  #find denominator for taus
  maxesIDs <- apply( tmp, 1, which.max)
  maxes <- apply( tmp, 1, max)
  tmp1 <- maxes + log( rowSums( exp( tmp - rep( maxes, times=length( logPis)))))
  tmp1 <- all.data$weights * tmp1
  
  return( sum( tmp1))
}


"calcMargPenLogl_SA4" <-
function( grps, all.data, K) 
{
  cat( grps, " ")
  #grps is a nSampGrps lengthed vector
  taus <- model.matrix( ~-1+as.factor( grps))
  pis <- colMeans( taus)

  #get the cond means
  means <- calcCondMeans( all.data, taus, K)

  #the conditional logl contributions
  condLogls <- calcCondProbs_SA( all.data, grps, means, K)

  #the contributions from logpi
  logPis <- log( pis)
  logPis <- logPis[grps]

  #the complete logl now
  unpenCLogl <- sum( condLogls) + sum( logPis)

  return( unpenCLogl)
}


"calcPis" <-
function( taus) 
{
  #numerator
  numer <- colSums( taus)
  #denomerator
  denom <- sum( taus)
  
  pis <- numer / denom

  return( pis)
}


"calcTaus" <-
function( pi, condProbs.fish, all.data, old.taus, update.prop=0.1, K, nu, method) 
{
  condProbs.sampGrps <- matrix( NA, ncol=K, nrow=all.data$nSampGrps)
  for( jj in levels( all.data$sample.grps))
    condProbs.sampGrps[jj==levels( all.data$sample.grps),] <- colSums( condProbs.fish[all.data$sample.grps==jj,,drop=FALSE])
  logPis <- log( pi)
  tmp <- condProbs.sampGrps + rep( logPis, each=all.data$nSampGrps)
  tmp <- tmp * rep( all.data$weights, times=K)
  if( method == "DA.EM")
    tmp <- tmp * nu

  #find denominator for taus
  maxesIDs <- apply( tmp, 1, which.max)
  maxes <- apply( tmp, 1, max)
  tmp1 <- maxes + log( rowSums( exp( tmp - rep( maxes, times=length( logPis)))))
  #numerator and denominator
  new.taus <- exp( tmp - rep( tmp1, times=length( logPis)))

  taus <- old.taus + update.prop*(new.taus-old.taus)

  return( taus)
}


"converged" <-
function( all, eps=1e-3) 
{
  tmp <- list()
  tmp[["pi"]] <- abs( all$new.parms$pi - all$old.parms$pi)
  tmp[["means"]] <- abs( all$new.parms$means - all$old.parms$means)
  my.max <- function(x){
    if( all( is.na( x)))
      return( Inf)
    tmpWithin <- max( x, na.rm=TRUE)
    return( tmpWithin)
  }
  tmp1 <- sapply( tmp, my.max)
  tmp2 <- all( tmp1 < eps)
  
  return( list( conv=tmp2, crit=tmp1))
}


"CVfun" <-
function( x, fm1, fold=5) 
{
  nholdout <- fm1$data$nFish %/% fold
  holdout <- sample( 1:fm1$data$nFish, nholdout, replace=FALSE)
  wts <- fm1$data$weights
  wts[holdout] <- 0 #remove these fish from the sample
  fm1$control$small.object <- FALSE
  fm1$control$method <- "EM"
  fm1$control$quiet <- TRUE
  fm.cv <- stockSTRUCTURE( SNPdata=fm1$data$SNPdata, sample.grps=fm1$data$sample.grps, K=fm1$K, weights=wts, start.grps=apply( fm1$postProbs, 1, which.max), control=fm1$control)
  condProbs.fish <- calcCondProbs( fm.cv$data, fm.cv$condMeans, fm.cv$K)
  taus <- calcTaus( fm.cv$pis, condProbs.fish, fm1$data, fm1$postProbs, 1, fm.cv$K, nu=1, method="EM")

  tmp.tab <- t( fm1$postProbs[holdout,]) %*% taus[holdout,]
  agree <- sum( diag( tmp.tab))
  tot <- sum( tmp.tab)
  
  prop <- agree / tot
   
  return( prop)  
}


"deleteGlobVars" <-
function() 
{ 
  rm( .SANN.maxit, .SANN.iter, .SANN.temp, .SANN.tmin, pos = pkg.globals)#".GlobalEnv")

  return( NULL)
}


"find.eta.for.DA" <-
function( nu0, m.steps) 
{
  #the eta required to get to optimisation of the full logl in m.steps
  eta <- exp( (log( 1)-log( nu0)) / m.steps) - 1
  
  return( eta)
}


"FormatGenotypes" <-
function (genotypes, ploidy) 
{
    formatted <- genotypes
    nall <- numeric(ncol(genotypes)/2)
    if (ploidy == 2) {
        for (ic in 1:(ncol(genotypes)/2)) {
            ens <- sort(unique(c(genotypes[, 2 * ic - 1], genotypes[, 
                2 * ic])))
            max.ens <- length(ens)
            for (il in 1:(dim(genotypes)[1])) {
                formatted[il, 2 * ic - 1] <- ifelse(is.na(genotypes[il, 
                  2 * ic - 1]), NA, (1:max.ens)[genotypes[il, 
                  2 * ic - 1] == ens])
                formatted[il, 2 * ic] <- ifelse(is.na(genotypes[il, 
                  2 * ic]), NA, (1:max.ens)[genotypes[il, 2 * 
                  ic] == ens])
            }
            nall[ic] <- max.ens
        }
    }
    if (ploidy == 1) {
        for (ic in 1:(ncol(genotypes))) {
            ens <- sort(unique(genotypes[, ic]))
            max.ens <- length(ens)
            for (il in 1:(dim(genotypes)[1])) {
                formatted[il, ic] <- ifelse(is.na(genotypes[il, 
                  ic]), NA, (1:max.ens)[genotypes[il, ic] == 
                  ens])
            }
            nall[ic] <- max.ens
        }
    }
    formatted <- as.matrix(formatted)
    list(genotypes = formatted, allele.numbers = nall)
}


"Fstat" <-
function (genotypes, npop, pop.mbrship, ploidy = 2) 
#function taken directly from unsupported Geneland package - Scott 15-Jan-2018
#minimal changes, just package name in calls.
{
    if (ploidy == 1) 
        stop("Fstat not implemented for haploid data")
    format <- FormatGenotypes(genotypes, ploidy = ploidy)
    genotypes <- format$genotypes
    allele.numbers <- format$allele.numbers
    if (sum(is.na(genotypes)) != 0) {
        warning("Genotypes contain missing values which might bias computations")
        genotypes[is.na(genotypes)] <- -999
    }
    Fis = rep(-999, npop)
    Fst = matrix(nrow = npop, ncol = npop, -999)
    if (npop > 1) {
        for (iclass1 in 1:(npop - 1)) {
            for (iclass2 in (iclass1 + 1):npop) {
                sub1 <- pop.mbrship == iclass1
                sub2 <- pop.mbrship == iclass2
                if ((sum(sub1) != 0) & (sum(sub2) != 0)) {
                  ztmp <- genotypes[sub1 | sub2, ]
                  nindivtmp <- nrow(ztmp)
                  pop.mbrshiptmp <- pop.mbrship[sub1 | sub2]
                  pop.mbrshiptmp[pop.mbrshiptmp == iclass1] <- 1
                  pop.mbrshiptmp[pop.mbrshiptmp == iclass2] <- 2
                  tabindiv <- matrix(nrow = nindivtmp, ncol = 2, 
                    data = -999)
                  kk <- numeric(2)
                  effcl <- table(pop.mbrshiptmp)
                  nloc <- length(allele.numbers)
                  nloc2 <- 2 * nloc
                  Fistmp <- Fsttmp <- Fittmp <- -999
                  res <- .Fortran("ggfst", PACKAGE = "stockR", 
                    as.integer(nindivtmp), as.integer(nloc), 
                    as.integer(nloc2), as.integer(allele.numbers), 
                    as.integer(2), as.integer(effcl), as.integer(ztmp), 
                    as.integer(pop.mbrshiptmp), as.integer(tabindiv), 
                    as.integer(kk), as.double(Fistmp), as.double(Fsttmp), 
                    as.double(Fittmp))
                  Fst[iclass1, iclass2] <- res[[12]][1]
                }
            }
        }
    }
    for (iclass1 in 1:npop) {
        sub1 <- pop.mbrship == iclass1
        if ((sum(sub1) != 0)) {
            ztmp <- rbind(genotypes[sub1, ], genotypes[sub1, 
                ])
            nindivtmp <- nrow(ztmp)
            pop.mbrshiptmp <- c(rep(1, sum(sub1)), rep(2, sum(sub1)))
            tabindiv <- matrix(nrow = nindivtmp, ncol = 2, data = -999)
            kk <- numeric(2)
            effcl <- table(pop.mbrshiptmp)
            nloc <- length(allele.numbers)
            nloc2 <- 2 * nloc
            Fistmp <- Fsttmp <- Fittmp <- -999
            res <- .Fortran("ggfst", PACKAGE = "stockR", as.integer(nindivtmp), 
                as.integer(nloc), as.integer(nloc2), as.integer(allele.numbers), 
                as.integer(2), as.integer(effcl), as.integer(ztmp), 
                as.integer(pop.mbrshiptmp), as.integer(tabindiv), 
                as.integer(kk), as.double(Fistmp), as.double(Fsttmp), 
                as.double(Fittmp))
            Fis[iclass1] <- res[[11]][1]
        }
    }
    Fst[lower.tri(Fst, diag = TRUE)] <- 0
    Fst <- Fst + t(Fst)
    Fis[Fis == -999] <- NA
    Fst[Fst == -999] <- NA
    list(Fis = Fis, Fst = Fst)
}


"get.init.taus" <-
function( grps, magic.num=0.8, K) 
{
  taus <- model.matrix( ~-1+factor( grps, levels=1:K))
  taus <- t( apply( taus, 1, function(x) ifelse( x==1, magic.num, (1-magic.num)/(ncol( taus)-1))))

  return( taus)
}


"get.start.grps" <-
function (control, K, all.data) 
{
    if (control$method == "Kmeans.EM") {
      if (!control$quiet) 
        message("Finding starting values using K-means on rotated data (using specified number of PCAs)")
      if( any( is.na( all.data$SNPdata))){
        tmpData <- all.data$SNPdata
        for( ii in 1:nrow( tmpData))
          tmpData[ii,is.na( tmpData[ii,])] <- mean( tmpData[ii,], na.rm=TRUE)
      }
      else
        tmpData <- all.data$SNPdata
      rotat <- prcomp( t( tmpData), scale.=FALSE, center=FALSE)$rotation
      PCAdat <- t( tmpData) %*% rotat[,1:min(control$nKPCA, ncol( rotat))]#min min needed to avoid situations where there are fewer markers than nKPCA
#      PCAdat <- t( tmpData)    #removed 20-5-2019
      if( sum( !duplicated( t( all.data$SNPdata))) <= K)  #added 3/6/2019 kmeans will through an error in any case.
        stop( "K is larger than the number of distinct individuals. Grouping is stupid?")
      cl <- kmeans( PCAdat, K, nstart=control$nKmeanStart)$cluster          
      tab <- table(cl, all.data$sample.grps)
      probs <- sweep(tab, 2, colSums(tab), FUN = "/")
      start.grps <- apply(probs, 2, function(x) sample(1:nrow(tab), 1, prob = x))
      tmpFun <- NULL
      while (is.null(start.grps) | length(unique(start.grps)) < K) {
        if (is.null(tmpFun)) {
          tmpFun <- function(x) {
            tmp <- 0.9 * (x - (1/K))
            tmp <- tmp + (1/K)
            tmp <- tmp/sum(tmp)
            return(tmp)
          }
        }
        probs <- apply(probs, 2, tmpFun)
        start.grps <- apply(probs, 2, function(x) sample(1:nrow(tab), 1, prob = x))
    }
    return(start.grps)
  }
  start.grps <- NULL
  while (is.null(start.grps) | length(unique(start.grps)) < K) {
    pis <- my.rdirichlet(1, K)
    start.grps <- sample(1:K, all.data$nSampGrps, replace = TRUE, prob = pis)
  }
  if (control$method != "SA.EM") {
    if (!control$quiet) 
      message("Finding starting values using random assignment")
    return(start.grps)
  }
  if (!control$quiet) 
    message("Finding starting values using SANN using the (complete) log-likelihood")
#  tmp1 <- optim(par = start.grps, fn = function(x) -calcCompleteLogl_SA(grps = x, 
#        all.data = all.data, K = K), gr = function(x) proposalFun2(x, 
#        all.data, K), method = "SANN", control = list(maxit = .SANN.maxit, 
#        tmax = 1, trace = control$SANN.trace, REPORT = control$SANN.nreport), 
#        hessian = FALSE)
  tmp1 <- optim(par = start.grps, fn = function(x) -calcCompleteLogl_SA(grps = x, 
        all.data = all.data, K = K), gr = function(x) proposalFun2(x, 
        all.data, K), method = "SANN", control = list(maxit = pkg.globals$.SANN.maxit, 
        tmax = 1, trace = control$SANN.trace, REPORT = control$SANN.nreport), 
        hessian = FALSE)
  return(tmp1$par)
}


"getIC" <-
function( margLogl, all.parms, K, all.data) 
{
  nparms <- K-1 + prod( dim( all.parms$new.parms$means))
  BIC <- -2*margLogl + nparms * log( all.data$nFish)
  AIC <- -2*margLogl + 2*nparms
  
  return( list( BIC=BIC, AIC=AIC))
}


"getMisMatch" <-
function( tab)
{
#  library( gtools)
  
  perms <- permutations( nrow( tab), nrow( tab), 1:nrow( tab))
  kount.wrong <- rep( NA, nrow( perms))
  for( ii in 1:nrow( perms)){
  kount.wrong[ii] <- sum( tab) - sum( tab[cbind( 1:nrow( tab), perms[ii,])])
  }
  
  return( min( kount.wrong))
}


"getReturnObject" <-
function( condMeans, all.parms, margLogl, taus, K, all.data, start.grps, allLogls, control) 
{
  ret <- list( condMeans=condMeans, pis=all.parms$new.parms$pi, margLogL=margLogl, postProbs=taus, K=K)
  if( control$small.object == FALSE){
    ret$data <- all.data
    ret$init.grps=start.grps
    ret$constantSNPs=as.numeric( which( all.data$constantSNPs))
    ret$margLogLTrace=allLogls
    ret$control=control
  }
  rownames( ret$condMeans) <- rownames( all.data$SNPdata)
  rownames( ret$postProbs) <- levels( all.data$sample.grps)
  colnames( ret$postProbs) <- colnames( ret$condMeans) <- names( ret$pis) <- paste( "Pop", 1:K, sep="_")  

  tmp <- getIC( margLogl, all.parms, K, all.data)
  ret$BIC <- tmp$BIC
  ret$AIC <- tmp$AIC
  
  class( ret) <- "stockSTRUCTURE.object"

  return( ret)
}


"getStockGroups" <-
function( K, nSampleGrps, minSize=1) 
{
  #exp( - runif(1,0,1) - 0.5*(0:(K-1)))
  #pis <- pis / sum( pis)
  grps <- NULL
  kount <- 0
  maxKount <- 1000
  while( ( is.null( grps) | any( table( grps)<minSize) | length( unique( grps))<K) & kount < maxKount){#length( unique( grps)) != K){
    pis <- my.rdirichlet( n=1, dim=K, alpha=1.5)#symmetric Dirichlet
    grps <- sample( 1:K, nSampleGrps, prob=pis, replace=TRUE)
    kount <- kount + 1
  }
  if( kount >= maxKount)
    stop("No group found of desired minimum sizes. You probably have the minimum size being quite large in relation to the number of groups")
  grps <- sort( grps) #the stocks of the sample groups

  return( grps)

}


"initialiseParms" <-
function( pis, means, all.data, K) 
{
  old.parms <- list( pi=rep( Inf, K), means=matrix( NA, nrow=all.data$nSNP, ncol=K))
  new.parms <- list( pi=pis, means=means)
    
  names( old.parms$pi) <- names( new.parms$pi) <- paste( "pi", 1:K, sep="_")
  rownames( old.parms$means) <- rownames( new.parms$means) <- paste( "SNP", 1:all.data$nSNP, sep="_")
  colnames( old.parms$means) <- colnames( new.parms$means) <- paste( "pop", 1:K, sep="_")

  return( list( old.parms=old.parms, new.parms=new.parms))

}


"iter.print" <-
function( control, kount, margLogl, nu=NA, all.parms, conv, last) 
{
  #any printing at all?
  if( control$quiet)
    return( invisible( NULL))
  #printing for EM and SA.EM
  if( control$method != "DA.EM"){
    if( last){  #convergence message?
      if( kount < control$EM.maxit)
        message( "Converged | Marg. Logl: ",round( margLogl,3), " | pis: ", paste( round( all.parms$new.parms$pi, 3), collapse=" "), " | Max change: ", round( max( conv$crit['pi'], conv$crit["means"]), 7))
      else
        message( "Not Converged")
      return(NULL)
    }
    message( "Iter: ", kount, " | Marg. Logl: ", round( margLogl, 3), " | pis: ", paste( round( all.parms$new.parms$pi, 3), collapse=" "), " | Max change: ", round( max( conv$crit['pi'], conv$crit["means"]), 7))
    return(NULL)
  }
  #printing for DA.EM
  if( control$method == "DA.EM"){
    if( last){  #convergence message?
      if( kount < control$EM.maxit)
        message( "Converged | Marg. Logl: ",round( margLogl,3), " | nu: ", format( round( nu, 6), nsmall=6), " | pis: ", paste( round( all.parms$new.parms$pi, 3), collapse=" "), " | Max change: ", round( max( conv$crit['pi'], conv$crit["means"]), 7))
      else
        message( "Not Converged")
      return(NULL)
    }
    message( "Iter: ", kount, " | Marg. Logl: ", round( margLogl, 3), " | nu: ", format( round( nu, 6), nsmall=6), " | pis: ", paste( round( all.parms$new.parms$pi, 3), collapse=" "), " | Max change: ", round( max( conv$crit['pi'], conv$crit["means"]), 7))
    return( NULL)
  }  
  
  stop( "Something has gone wrong in iter.print()")
  return(NULL)
}


"K1special" <-
function( all.data, control) 
{
  #just taking means of data etc...  K=1 is largely used as a base-case.
  weights.fish <- all.data$weights[all.data$sample.grps.num]
  condMeans <- rep( weights.fish, each=all.data$nSNP) * all.data$SNPdata
  FishPerSNP <- apply( condMeans, 1, function(x) which( !is.na( x)))
  if( is( FishPerSNP, "matrix"))
    FishPerSNP <- as.list( as.data.frame( FishPerSNP))
  totWtsPerSNP <- sapply( FishPerSNP, function(x) sum( weights.fish[x]))
#  condMeans <- rowSums( condMeans, na.rm=TRUE) / sum( weights.fish) # totWtsPerSNP#the weighted mean
  condMeans <- rowSums( condMeans, na.rm=TRUE) / totWtsPerSNP#the weighted mean
  condMeans <- matrix( condMeans / 2, ncol=1)
  
  tmp.probs <- condMeans
  #A hack to stop alleles with exactly zero or one instances through 
  #replace exact zeros by a small value and exact ones by a small value.
  #shouldn't be a problem for standard data analysis, just cross-validation/Bootstrap
  #where fish are removed
  tmp.probs[tmp.probs==0] <- (1/all.data$nFish)/2
  tmp.probs[tmp.probs==1] <- 1-(1/all.data$nFish)/2
  #now get the probs 'safely'
  tmp <- dbinom( as.matrix( all.data$SNPdata)[,drop=FALSE], size=2, prob=rep( tmp.probs, all.data$nFish), log=TRUE)
  condProbs <- matrix( colSums( tmp, na.rm=TRUE), ncol=1)  #these aren't really conditional probs as there is only 1 group...
  margLogl <- sum( condProbs)
  
  ret <- list( condMeans=condMeans, pis=1, margLogL=margLogl, postProbs=matrix( 1, nrow=all.data$nSampGrps, ncol=1), K=1)
  if( control$small.object == FALSE){
    ret$data <- all.data
    ret$init.grps=rep( 1, all.data$nSampGrps)
    ret$constantSNPs=as.numeric( all.data$constantSNPs)
    ret$margLogLTrace=NA
    ret$control=control
  }
  rownames( ret$condMeans) <- rownames( all.data$SNPdata)
  rownames( ret$postProbs) <- levels( all.data$sample.grps)
  colnames( ret$postProbs) <- colnames( ret$condMeans) <- names( ret$pis) <- "Pop_1"  

  tmp.parms <- list()
  tmp.parms$new.parms <- list()
  tmp.parms$new.parms$means <- condMeans

  tmp <- getIC( margLogl, all.parms=tmp.parms, K=1, all.data)
  ret$BIC <- tmp$BIC
  ret$AIC <- tmp$AIC

  return( ret)
}


"my.rdirichlet" <-
function( n, dim, alpha=1) 
{
  alpha <- rep( alpha, dim)
  gamma.rvs <- rgamma( n*dim, scale=1, shape=rep( alpha, each=n))
  dir.rvs <- matrix( gamma.rvs, nrow=n, ncol=dim)
  dir.rvs <- sweep( dir.rvs, 1, STATS=rowSums( dir.rvs), FUN="/")
  
  return( dir.rvs)
}



"plot.stockBOOT.object" <- 
function (x, locations = NULL, plotTitle = NULL, CI.width=0.95, region.lwd=3.5, ...) 
{
  if (is.null(x) | !is(x, "stockBOOT.object"))
    stop("Need stockBOOT.object from stockBOOT function.  Uncertainty cannot be displayed without this object")

  #to stop new par settings being delivered to user AFTER function exit.  Thanks CRAN maintainer Victoria Wimmer!
  oldpar <- par(no.readonly = TRUE)    # code line i
  on.exit(par(oldpar))            # code line i + 1

  nFish <- dim(x$postProbs)[1]
  nGrps <- dim(x$postProbs)[2]
  if (is.null(locations)) {
    message("No locations supplied. Plotting as a single location.")
    locations <- data.frame(region = rep("", nFish), regionOrder = rep(1, nFish))
  }
  if ( is(locations, "matrix")) 
    locations <- as.data.frame(locations)
  myOrd <- order(locations[, 2])
  #reorder the important things
  posty <- x$postProbs[myOrd,,,drop=FALSE]
  locations <- locations[myOrd,]
  #get average posterior
  vals <- apply(posty, 1:2, mean)
  #get posterior range (of 95% ci)
  tmpQuant <- apply(posty, 1:2, quantile, probs = c( (1-CI.width) / 2, CI.width + (1-CI.width) / 2))#c(0.025, 0.975))
  tmpRange <- matrix( tmpQuant[2,,] - tmpQuant[1,,], nrow=dim( tmpQuant)[2], ncol=dim( tmpQuant)[3])#ncol=dim( tmpQuant)[1]-1)
  
  #padding out the matrix with white space...
  tmp <- !duplicated(locations[, 1])
  tmp1 <- tmp1.expand <- c(which(tmp), nFish + 1)
  vals1 <- matrix(NA, nrow = nFish + 2 * (length(tmp1) - 2), ncol = ncol(vals))
  curr.pos <- 1
  for (jj in 1:(length(tmp1) - 1)) {
    new.ids <- tmp1[jj]:(tmp1[jj + 1] - 1)
    vals1[curr.pos:(curr.pos + length(new.ids) - 1), ] <- vals[new.ids, ]
    curr.pos <- curr.pos + length(new.ids) + 2
  }
  #correct the colour for width of CI
  if (nGrps < 3) 
    myCols <- RColorBrewer::brewer.pal(3, "Set1")[1:nGrps]
  else myCols <- RColorBrewer::brewer.pal(nGrps, "Set1")
  myCols1 <- matrix(rep(myCols, each = nrow(vals1)), nrow = nrow(vals1), ncol = ncol(vals1))
  present.ids <- !apply(vals1, 1, function(x) all(is.na(x)))
  if( ncol( myCols1) > 1)
    for (zz in 1:nrow(tmpRange)) 
      myCols1[present.ids,][zz, ] <- adjustcolor(myCols1[present.ids,][zz, ], alpha.f = min(1, (1 - tmpRange[zz, ]) + 0.00))
  rownames(vals1) <- NULL
  par(mfrow = c(1, 1), oma = c(1, 0, 3, 0), mar = c(4, 0, 2, 0))
  bp <- barplot(t(vals1), col = rep(adjustcolor("white", alpha.f = 0), 
                                    nGrps), border = NA, space = 0, yaxt = "n", main = "", 
                xlab = "", cex.main = 2, cex.names = 1, xpd = TRUE, plot = TRUE, ...)
  for (nn in 1:nrow(vals1)) {
    bot <- 0
    for (mm in 1:ncol(vals1)) {
      topp <- bot + vals1[nn, mm]
      rect(xleft = bp[nn] - 0.5, ybottom = bot, xright = bp[nn] + 0.5, ytop = topp, col = myCols1[nn, mm], border = grey(0.75), lwd = 0.1)
      bot <- topp
    }
  }
  if( any( is.na( vals1))){
    tmppy <- which( is.na( vals1[,1]))
    lefts <- c( -1, tmppy[seq( from=1, to=length( tmppy), by=2)])
    rights <- c( tmppy[seq( from=2, to=length( tmppy), by=2)], nrow( vals1)+2)
    for( ii in 1:length( lefts))
      rect( xleft=lefts[ii]+1, xright=rights[ii]-2, ybottom=0, ytop=1, lwd=region.lwd)
    n.to.print <- tapply(locations[, 2], locations[, 2], length)
    names.to.print <- unique(as.character(locations[, 1]))
    start.pos <- 1
    for (jj in seq(from = 1, to = length(tmp1) - 1, by = 1)) {
      end.pos <- start.pos + diff(tmp1)[jj] - 1
      loc <- floor((bp[start.pos] + bp[end.pos])/2)
      text(loc, -0.02, names.to.print[jj], adj = c(0.5, 1), xpd = TRUE)
      text(loc, -0.13, paste("(n=", n.to.print[jj], ")", sep = ""), xpd = TRUE, adj = c(0.5, 1))
      bit.to.add <- tmp1[jj + 1] - tmp1[jj]
      start.pos <- start.pos + bit.to.add + 2
    }
  }
  else
    rect( xleft=0, xright=nrow( vals1), ybottom=0, ytop=1, lwd=region.lwd)
  mtext(plotTitle, side = 3, line = 0, outer = TRUE, cex = 1.75)
  invisible(NULL)
}


"proposalFun2" <-
function( grps, all.data, K, ...){

  #getting the current annealing temperature
#  assign( ".SANN.iter", .SANN.iter+1, pkg.globals)#.GlobalEnv)#.SANN.iter <<- .SANN.iter + 1
  pkg.globals$.SANN.iter <- pkg.globals$.SANN.iter + 1
  inner.tmp <- pkg.globals$.SANN.temp - ( pkg.globals$.SANN.temp / pkg.globals$.SANN.maxit) * pkg.globals$.SANN.iter
  if( inner.tmp < pkg.globals$.SANN.tmin)
    inner.tmp <- pkg.globals$.SANN.tmin

  #choosing which sampling groups to change cluster
  randNum <- rtriangular( n=1, a=0, b=inner.tmp, mode=inner.tmp/3)
numToExchange <- ceiling( randNum)
toChange <- sample( levels( all.data$sample.grps), numToExchange, replace=FALSE)
toKeep <- setdiff( levels( all.data$sample.grps), toChange)

  newGrps <- grps
for( ii in toChange)
  newGrps[levels( all.data$sample.grps)==ii] <- sample( setdiff( 1:K, grps[levels( all.data$sample.grps)==ii]), 1)
#  newGrps[levels( all.data$sample.grps) %in% toChange] <- sample( 1:K, numToExchange)
  
  #which sample groups not retained?
  grpsToAdd <- which( ! 1:K %in% unique( newGrps))
  while( length( grpsToAdd) != 0){
#    grpsToAdd <- which( ! 1:K %in% unique( newGrps))
    newGrps[sample(1:all.data$nSampGrps, 1)] <- grpsToAdd[1]
    grpsToAdd <- which( ! 1:K %in% unique( newGrps))
  }

  return( newGrps) 
}


"rtriangular" <-
function( n=1, a=0, b=1, mode=1/3)
{
  u <- runif( n)
  Fmode <- (mode-a) / (b-a)
  x <- ifelse( u < Fmode, a + sqrt( u*(b-a)*(mode-a)), b - sqrt( (1-u)*(b-a)*(b-mode)))
    
  return( x)
}


"set.control" <-
function( contr, nSampGrps, K, nFish, nSNP) 
{
  ret <- list() #return control object
  ret$small.object <- contr$small.object
  if( !"small.object" %in% names( contr))
    ret$small.object <- FALSE
  ret$quiet <- contr$quiet
  if( !"quiet" %in% names( contr))
    ret$quiet <- FALSE
  ret$method <- contr$method
  if( !"method" %in% names( contr))
    ret$method <- "Kmeans.EM" #default optimisation method

  if( ret$method == "Kmeans.EM"){
    if( !"nKmeanStart" %in% names( contr))
      ret$nKmeanStart <- 25
    else
      ret$nKmeanStart <- contr$nKmeanStart
    if( !"nKPCA" %in% names( contr))
      ret$nKPCA <- min( nFish, nSNP, 100)
    else
      ret$nKPCA <- min( nFish, nSNP, contr$nKPCA)
    ret$EM.eps <- contr$EM.eps
    if( !"EM.eps" %in% names( contr))
      ret$EM.eps <- 1e-5
    ret$EM.maxit <- contr$EM.maxit
    if( !"EM.maxit" %in% names( contr))
     ret$EM.maxit <- 100
    ret$EM.minit <- contr$EM.minit
    if( !"EM.minit" %in% names( contr))
      ret$EM.minit <- 5
    ret$EM.tau.step <- contr$EM.tau.step
    if( !"EM.tau.step" %in% names( contr))
      ret$EM.tau.step <- 0.5
    ret$EM.tau.init.max <- contr$EM.tau.init.max
    if( !"EM.tau.init.max" %in% names( contr))
      ret$EM.tau.init.max <- 2 / (K+1) #so that there is twice the mass in the most likely group.
 
    return( ret)
  }
  if( ret$method == "EM"){
    #EM particulars
    ret$EM.eps <- contr$EM.eps
    if( !"EM.eps" %in% names( contr))
      ret$EM.eps <- 1e-5
    ret$EM.maxit <- contr$EM.maxit
    if( !"EM.maxit" %in% names( contr))
     ret$EM.maxit <- 100
    ret$EM.minit <- contr$EM.minit
    if( !"EM.minit" %in% names( contr))
      ret$EM.minit <- 1
    ret$EM.tau.step <- contr$EM.tau.step
    if( !"EM.tau.step" %in% names( contr))
      ret$EM.tau.step <- 1
    ret$EM.tau.init.max <- contr$EM.tau.init.max
    if( !"EM.tau.init.max" %in% names( contr))
      ret$EM.tau.init.max <- 2 / (K+1) #so that there is twice the mass in the most likely group.
 
    return( ret)
  }
  if( ret$method=="SA.EM"){
    #global variables to make optim (method="SANN") work
    assign(".SANN.iter", 0, pkg.globals)#environment(stockSTRUCTURE))# .GlobalEnv)#.SANN.iter <<- 0
    if( !"SANN.temp" %in% names( contr))
      assign( ".SANN.temp", nSampGrps/2, pkg.globals)#environment(stockSTRUCTURE))#.GlobalEnv)#.SANN.temp <<- nSampGrps / 2
    else
      assign( ".SANN.temp", contr$SANN.temp, pkg.globals)#, environment(stockSTRUCTURE))#.GlobalEnv)#      .SANN.temp <<- contr$SANN.temp
    if( !"SANN.maxit" %in% names( contr))
      assign( ".SANN.maxit", 5000, pkg.globals)#, environment(stockSTRUCTURE))#.GlobalEnv)#.SANN.maxit <<- 5000
    else
      assign( ".SANN.maxit", contr$SANN.maxit, pkg.globals)#, environment(stockSTRUCTURE))#.GlobalEnv)#.SANN.maxit <<- contr$SANN.maxit
    if( !"SANN.tmin" %in% names( contr))
      assign( ".SANN.tmin", 1, pkg.globals)#, environment(stockSTRUCTURE))#.GlobalEnv)#.SANN.tmin <<- 1#max( nSampGrps %/% 20, 1)
    else
      assign( ".SANN.tmin", contr$SANN.temp, pkg.globals)#, environment(stockSTRUCTURE))#.GlobalEnv)#.SANN.tmin <<- contr$SANN.tmin
    ret$SANN.nreport <- contr$SANN.nreport
    if( !"SANN.nreport" %in% names( contr))
      ret$SANN.nreport <- 100
    ret$SANN.trace <- contr$SANN.trace
    if( !"SANN.trace" %in% names( contr))
      ret$SANN.trace <- TRUE
    #EM particulars for SA.EM
    ret$EM.eps <- contr$EM.eps
    if( !"EM.eps" %in% names( contr))
      ret$EM.eps <- 1e-5
    ret$EM.maxit <- contr$EM.maxit
    if( !"EM.maxit" %in% names( contr))
     ret$EM.maxit <- 100
    ret$EM.minit <- contr$EM.minit
    if( !"EM.minit" %in% names( contr))
      ret$EM.minit <- 1
    ret$EM.tau.step <- contr$EM.tau.step
    if( !"EM.tau.step" %in% names( contr))
      ret$EM.tau.step <- 1
    ret$EM.tau.init.max <- contr$EM.tau.init.max
    if( !"EM.tau.init.max" %in% names( contr))
      ret$EM.tau.init.max <- 2 / (K+1) #so that there is twice the mass in the most likely group.

    return( ret)
  }
  if( ret$method=="DA.EM"){
    #EM particulars for DA.EM
    ret$EM.eps <- contr$EM.eps
    if( !"EM.eps" %in% names( contr))
      ret$EM.eps <- 1e-5
    ret$EM.maxit <- contr$EM.maxit
    if( !"EM.maxit" %in% names( contr))
     ret$EM.maxit <- 100
    ret$EM.minit <- contr$EM.minit
    if( !"EM.minit" %in% names( contr))
      ret$EM.minit <- 25
    ret$EM.tau.step <- contr$EM.tau.step
    if( !"EM.tau.step" %in% names( contr))
      ret$EM.tau.step <- 1
    ret$EM.tau.init.max <- contr$EM.tau.init.max
    if( !"EM.tau.init.max" %in% names( contr))
      ret$EM.tau.init.max <- (1+0.1) / K #so that there is a tiny bit more in the most likely group.
    ret$DA.nu.init <- contr$DA.nu.init
    if( !"DA.nu.init" %in% names( contr)){
      ret$DA.nu.init <- 1e-4
#      N_here <- max( table( all.data$sample.grps.num) * all.data$nSNP)
#      ret$NA.nu.init <- exp( log( 0.1) / N_here)
    }
    ret$DA.eta <- find.eta.for.DA( ret$DA.nu.init, ret$EM.minit)

    return( ret)
  }
  
#  if( !"mc.cores" %in% names( contr))
#    ret$mc.cores <- 4
}


"setData" <-
function( SNPdata, sample.grps, weights, control) 
{
  ret <- list()
  ret$constantSNPs <- apply( SNPdata, 1, function(x) length( unique( x, na.rm=TRUE))==1)
  if( !control$quiet)
    message( "Removing any SNP with no variation amongst *any* fish (there are ", sum( ret$constantSNPs)," of them)")
  ret$SNPdata <- SNPdata[!ret$constantSNPs,]
  ret$weights <- as.numeric( weights)
  ret$nSNP <- nrow( ret$SNPdata)
  ret$nFish <- ncol( ret$SNPdata)
  ret$sample.grps <- as.factor( sample.grps)
  ret$sample.grps <- droplevels( ret$sample.grps)
  ret$nSampGrps <- nlevels( ret$sample.grps)
  ret$sample.grps.num <- as.numeric( ret$sample.grps)
  
  return( ret)
}


"sim.stock.data" <-
function( nAnimals, nSNPs, nSampleGrps, K, ninform=nSNPs, means=c(alpha=-1.9,beta=0), sds=c(alpha=1.6,beta.inform=0.85,beta.noise=0.0005), minStockSize=1) 
{
  grps <- getStockGroups( K, nSampleGrps, minStockSize)
  sampleGrps <- rep( paste( "samp",1:nSampleGrps, sep=""), (nAnimals %/% nSampleGrps) + 1)[1:nAnimals]  #the animal to the sample groups
  sampleGrps <- as.factor( sort( sampleGrps)) #ordering doesn't matter so it might as well be neat
      
  #getting the SNP expected frequencies
  intercepts <- rnorm( nSNPs, mean=means["alpha"], sd=sds["alpha"])
  betas <- matrix( NA, nrow=K, ncol=nSNPs)
  for( kk in 1:K){
    betas[kk,1:min( ninform, nSNPs)] <- rnorm( min( ninform, nSNPs), mean=means['beta'], sd=sds['beta.inform']) #only the first ninform have reall deviation
    if( nSNPs > ninform)
      betas[kk,( ninform+1):nSNPs] <- rnorm( nSNPs-ninform, mean=means['beta'], sd=sds['beta.noise'])
  }
  betas <- sweep( betas, 2, STATS=colMeans( betas), FUN='-')
  lps <- sweep( betas, 2, STATS=intercepts, FUN='+')
  means <- exp( lps) / (1+exp( lps))

  #now get the SNP data
  SNPs <- matrix( NA, nrow=nAnimals, ncol=nSNPs)
  for( ii in 1:nAnimals)
    SNPs[ii,] <- rbinom( nSNPs, size=2, prob=means[grps[sampleGrps[ii]==levels( sampleGrps)],])#Effing awful expression that gets the right mean frequecy for this animal
  SNPs <- t( SNPs)
  
  attr( SNPs, "grps") <- grps
  attr( SNPs, "sampleGrps") <- sampleGrps
  
  return( SNPs)
}


"sim.stock.data.old" <-
function( nAnimals, nSNPs, nSampleGrps, K, ninform=nSNPs, means=c(alpha=-1.9,beta=0), sds=c(alpha=1.6,beta.inform=0.85,beta.noise=0.0005), minStockSize=1) 
{
  grps <- getStockGroups( K, nSampleGrps, minStockSize)
  sampleGrps <- rep( paste( "samp",1:nSampleGrps, sep=""), (nAnimals %/% nSampleGrps) + 1)[1:nAnimals]  #the animal to the sample groups
  sampleGrps <- as.factor( sort( sampleGrps)) #ordering doesn't matter so it might as well be neat
      
  #getting the SNP expected frequencies
  intercepts <- rnorm( nSNPs, mean=means["alpha"], sd=sds["alpha"])
  betas <- matrix( NA, nrow=K, ncol=nSNPs)
  for( kk in 1:K){
    betas[kk,1:min( ninform, nSNPs)] <- rnorm( min( ninform, nSNPs), mean=means['beta'], sd=sds['beta.inform']) #only the first ninform have reall deviation
    if( nSNPs > ninform)
      betas[kk,( ninform+1):nSNPs] <- rnorm( nSNPs-ninform, mean=means['beta'], sd=sds['beta.noise'])
  }
  betas <- sweep( betas, 2, STATS=colMeans( betas), FUN='-')
  lps <- sweep( betas, 2, STATS=intercepts, FUN='+')
  means <- exp( lps) / (1+exp( lps))

  #now get the SNP data
  SNPs <- matrix( NA, nrow=nAnimals, ncol=nSNPs)
  for( ii in 1:nAnimals)
    SNPs[ii,] <- rbinom( nSNPs, size=2, prob=means[grps[sampleGrps[ii]==levels( sampleGrps)],])#Effing awful expression that gets the right mean frequecy for this animal
  SNPs <- t( SNPs)
  
  attr( SNPs, "grps") <- grps
  attr( SNPs, "sampleGrps") <- sampleGrps
  
  return( SNPs)
}


"stockBOOT" <-
function( fm, B=100, mc.cores=NULL) 
{
  if( ! "data" %in% names( fm) | ! "constantSNPs" %in% names( fm))
    stop( "Need data not contained in model object.  Please supply this information (try running stockSTRUCTURE again with (default) argument control=list( small.object=FALSE))")
  
  if( is.null( mc.cores))
    mc.cores <- max( 1, parallel::detectCores())
  cl <- parallel::makeCluster( mc.cores)
  parallel::clusterExport(cl, list( "fm"), envir=environment())
  #parallel::clusterExport( cl, ls( envir=environment( stockR:::bootFun)), envir=environment( stockR:::bootFun))
  parallel::clusterExport( cl, ls( envir=environment( bootFun)), envir=environment( bootFun))
  tmp <- parallel::parLapply( cl, 1:B, bootFun, fm1=fm)
  parallel::stopCluster(cl)

  condMeans <- array( NA, dim=c( nrow(tmp[[1]]$condMeans), ncol( tmp[[1]]$condMeans),B), dimnames=list( rownames( tmp[[1]]$condMeans), colnames( tmp[[1]]$condMeans, paste("bootSamp",1:B,sep=""))))
  pis <- matrix( NA, nrow=B, ncol=length( tmp[[1]]$pis), dimnames=list( paste("bootSamp",1:B,sep=""), names( tmp[[1]]$pis)))
  postProbs <- array( NA, dim=c( nrow( tmp[[1]]$postProbs), ncol(tmp[[1]]$postProbs), B), dimnames=list( rownames( tmp[[1]]$postProbs), colnames( tmp[[1]]$postProbs), paste("bootSamp",1:B,sep="")))
  margLogL <- rep( NA, B); names( margLogL) <- paste( "bootSamp", 1:B,sep="")
  for( bb in  1:B){
    condMeans[,,bb] <- tmp[[bb]]$condMeans
    pis[bb,] <- tmp[[bb]]$pis
    postProbs[,,bb] <- tmp[[bb]]$postProbs
    margLogL[bb] <- tmp[[bb]]$margLogL    
  }

  ret <- list( condMeans=condMeans, pis=pis, postProbs=postProbs, margLogL=margLogL, B=B)
  class( ret) <- "stockBOOT.object"
  return( ret)  
}


"stockCV" <-
function(fm, fold=5, B=25, mc.cores=NULL) 
{
  if( ! "data" %in% names( fm) | ! "constantSNPs" %in% names( fm))
    stop( "Need data not contained in model object.  Please supply this information (try running stockSTRUCTURE again with (default) argument control=list( small.object=FALSE))")
  
  if( is.null( mc.cores))
    mc.cores <- max( 1, parallel::detectCores())
#  cl <- parallel::makeCluster( mc.cores)
##  parallel::clusterExport(cl, list( "fm"), envir=environment())
#  parallel::clusterEvalQ( cl, library( stockR))
##  parallel::clusterExport( cl, list( "stockR:::CVfun", "stockSTRUCTURE", "stockR:::setData", "stockR:::set.control", "stockR:::get.init.taus", "stockR:::calcCondMeans", "stockR:::initialiseParms", "stockR:::calcCondProbs", "stockR:::calcMargLogl", "stockR:::converged", "stockR:::updateOldParms", "stockR:::updateNewParms", "stockR:::calcTaus", "stockR:::deleteGlobVars", "stockR:::calcPis", "stockR:::iter.print", "stockR:::getReturnObject", "stockR:::getIC","stockR:::K1special"))
#  tmp <- parallel::parLapply( cl, 1:B, CVfun, fm1=fm, fold=fold)
#  parallel::stopCluster(cl)
  tmp <- parallel::mclapply( 1:B, CVfun, fm1=fm, fold=fold, mc.cores=mc.cores)

  tmp1 <- as.numeric( do.call( "cbind", tmp))
  attr( tmp1, "pis") <- fm$pis
  
  return( tmp1)
}


"stockSTRUCTURE" <-
function(SNPdata=NULL, sample.grps=as.factor( 1:ncol( SNPdata)), K=3, weights=rep( 1, nlevels( sample.grps)), start.grps=NULL, control=NULL)
{
  #sanity checking  
  if( is.null( SNPdata))
    stop( "You have not given the function any data! Please correct this and try again")
  if( K<1)
    stop( "You have specified less than 1 stock. This function aims to find two or more stocks")

  #set the control bits 'n' bobs (both for SANN and EM and DA.EM parts)
  control <- set.control( control, length( unique( as.character( sample.grps))), K, ncol( SNPdata), nrow( SNPdata))
  
  #convenience object for data storage
  all.data <- setData( SNPdata, sample.grps, weights, control)

  #more sanity checking
  if( length( all.data$weights) != all.data$nSampGrps)
    stop( "Legnths of sample.grps and weights are different! Please correct this and try again")
  if( !is.null( all.data$start.grps) & length( unique( all.data$start.grps)) < K)
    stop( "The starting structure supplied does not contain K groups.  Please correct this and try again")
  if( !is.null( start.grps) & length( start.grps) != all.data$nSampGrps)
    stop( "You must provide the same number of start.grps as there are sample.grps.  Please correct and try again")
  if( all.data$nSampGrps < K)
    stop( "There are fewer things (sample.grps) to cluster than K! There can be no clustering. Please correct and try again")

  #the special case of K=1
  if( K==1){
    ret <- K1special( all.data, control)
    message( "Done! (Simply by taking means etc of data)")
    return( ret)
  }

  if( is.null( start.grps))
    start.grps <- get.start.grps( control, K, all.data)
  
  #initial postprobs
  taus <- get.init.taus( start.grps, magic.num=control$EM.tau.init.max, K=K)
  #marginal probs
  pis <- calcPis( taus)
  #initial means
  condMeans <- calcCondMeans( all.data, taus, K)
  #initial parameter list
  all.parms <- initialiseParms( pis, condMeans, all.data, K)

  #now do EM optimisation for 'real' params.
  if( !control$quiet){
    if( control$method != "DA.EM")
      message( "Maximising the log-likelihood to obtain parameter estimates")
    else
      message( "Maximising the DA log-likelihood to obtain parameter estimates")
  }
  kount <- 0
  allLogls <- rep( NA, 1+control$EM.maxit)
  #This is done just to get the logl for the first iteration
  condProbs.fish <- calcCondProbs( all.data, condMeans, K)
  allLogls[kount+1] <- margLogl <- calcMargLogl( all.parms$new.parms, all.data, condProbs.fish, K)
  conv <- converged( all.parms, control$EM.eps) 
  if( control$method == "DA.EM")    
    nu <- control$DA.nu.init
  #now do an M step, then an E step, until convergence
  repeat{
    iter.print( control, kount, margLogl, nu, all.parms, conv, last=FALSE)
    kount <- kount + 1
    #set new parms to old parms
    all.parms <- updateOldParms( all.parms)
    #update the mean parameters
    condMeans <- calcCondMeans( all.data, taus, K)
    all.parms <- updateNewParms( all.parms, condMeans)
    #get logl contributions, marg logl, and taus
    condProbs.fish <- calcCondProbs( all.data, all.parms$new.parms$means, K)
    #tempering the importance of starting groups -- only in control$EM.tau.step < 1 and when kount < control$EM.minit
    tau.step_iterated <- min( control$EM.tau.step + (1-control$EM.tau.step) * ((kount-1)/max( control$EM.minit,1)), 1)
    #update taus (possibly tempered)
    taus <- calcTaus( all.parms$new.parms$pi, condProbs.fish, all.data, taus, tau.step_iterated, K, nu=nu, method=control$method)
    #update pis
    all.parms$new.parms$pi <- calcPis( taus)#colMeans( taus)
    #assess convergence
    conv <- converged( all.parms, control$EM.eps)
    allLogls[kount+1] <- margLogl <- calcMargLogl( all.parms$new.parms, all.data, condProbs.fish, K)
    #update
    if( control$method=="DA.EM")
      nu <- min( 1, nu+control$DA.eta)#(1+control$DA.eta)*nu)
    if( (conv$conv | kount > control$EM.maxit) & kount > control$EM.minit)
      break
  }
  iter.print( control, kount, margLogl, nu, all.parms, conv, last=TRUE)
  #final update of logls
  if( kount < control$EM.maxit)
    allLogls <- allLogls[-((kount):control$EM.maxit+1)]
  
  #bundling up for return
  ret <- getReturnObject( condMeans, all.parms, margLogl, taus, K, all.data, start.grps, allLogls, control)
  
  #cleaning up
  if( control$method=="SA.EM")
    deleteGlobVars()
  if( !control$quiet)
    message( "Done!")  

  return( ret)
}


"updateNewParms" <-
function( all.parms, condMeans) 
{
  all.parms$new.parms$means <- condMeans
  return( all.parms)
}


"updateOldParms" <-
function( parms) 
{
  parms$old.parms <- parms$new.parms
  return( parms)
}

