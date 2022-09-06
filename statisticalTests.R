#statisticalTests.R
#Charley Wu, 2019
library(BayesFactor)
library(lsr)
####################################################################################
# Statistical Tests
####################################################################################
#Standard t-test
ttestPretty <-function(x,y=NULL,mu=0, var.equal=T, paired=F, maxBF = 100){
  if (!is.null(y)){ #two sample
    freq <- t.test(x,y, var.equal=var.equal, paired = paired)
    d <- cohensD(x = x, y = y)
  }else{ #single sample
    freq <- t.test(x,mu = mu, var.equal=var.equal)
    d <- cohensD(x = x,mu = mu)
  }
  dof <- freq$parameter
  t <- sprintf("%.1f",freq$statistic) 
  p <- pformat(freq$p.value)
  if (!is.null(y)){ #single sample
    BF <-bfformat(extractBF(ttestBF(x,y, paired=paired))$bf, maxValue = maxBF)
  }else{#two sample
    BF <-bfformat(extractBF(ttestBF(x,mu = mu, paired=paired))$bf, maxValue = maxBF)
  }
  return(paste0('$t(',dof,')=',t,'$, $p',p, '$, $d',dFormat(d), '$, $BF',BF,'$'))
}

#Non-parametric tests
ranktestPretty <-function(x,y, oneSample = F, paired=F, maxBF = 100){
  if (oneSample==T){ #one sample wilcoxon signed rank test
    w <- wilcox.test(x, mu=0)
    z <- sprintf("%.1f",qnorm(w$p.value)) #z-statistic
    p <- pformat(w$p.value) #p-value
    r <- corformat(qnorm(w$p.value)/sqrt(length(x))) #effect size; r value
    #Bayes factor
    outsim<-signRankGibbsSampler(x, progBar = T) #sign rank gibbs sampler 
    dense<- density(outsim$deltaSamples)
    ddense <- with(dense, approxfun(x, y, rule=1))
    ifelse(is.na(ddense(0)), denominator <- .Machine$double.xmin, denominator <- ddense(0)) #if no density at 0, then default to the smallest number
    BF<-dcauchy(0, location = 0, scale = 1/sqrt(2), log = FALSE)/denominator
    out <- paste0('$Z=', z , '$, $p',p, '$, $r=',r,'$, $BF', bfformat(BF),'$')
  }else{
    if (paired==T){ #wilcoxon signed rank test (for paired comparisons) 
      w <- wilcox.test(x,y , paired=T)
      z <- sprintf("%.1f",qnorm(w$p.value)) #z-statistic
      p <- pformat(w$p.value) #p-value
      r <- corformat(qnorm(w$p.value)/sqrt(length(x))) #effect size; r value
      #Bayes factor
      outsim<-signRankGibbsSampler(x, y, progBar = T) #sign rank gibbs sampler 
      dense<- density(outsim$deltaSamples)
      ddense <- with(dense, approxfun(x, y, rule=1))
      ifelse(is.na(ddense(0)), denominator <- .Machine$double.xmin, denominator <- ddense(0)) #if no density at 0, then default to the smallest number
      BF<-dcauchy(0, location = 0, scale = 1/sqrt(2), log = FALSE)/denominator
      out <- paste0('$Z=', z , '$, $p',p, '$, $r=',r,'$, $BF', bfformat(BF),'$')
    }else{#Mann whitney U test (for non-paired comparisons)
      w <- wilcox.test(x,y)
      U <- sprintf('%.0f',w$statistic) #extract test statistic (no significant digits)
      p <- pformat(w$p.value) #p-value
      r <- corformat(cor(c(rep(1,length(x)), rep(0,length(y))), c(x,y), method='kendall')) #rank correlation as effect size, comparing data gainst binary identification vector
      #Bayes factor
      outsim<-rankSumGibbsSampler(x, y, progBar = T) #rank sum gibbs sampler 
      dense<- density(outsim$deltaSamples)
      ddense <- with(dense, approxfun(x, y, rule=1))
      ifelse(is.na(ddense(0)), denominator <- .Machine$double.xmin, denominator <- ddense(0)) #if no density at 0, then default to the smallest number
      BF<-dcauchy(0, location = 0, scale = 1/sqrt(2), log = FALSE)/denominator
      out <- paste0('$U=', U , '$, $p',p, '$, $r_{\tau}=',r,'$, $BF', bfformat(BF),'$')
    }
  }
  
  return(out)
}

#Correlation test, for either pearson or kendall's tau
corTestPretty <-function(x,y, maxBF = 100, method = 'pearson'){
  freq <- cor.test(x,y, method = method)
  r <- corformat(freq$estimate) 
  p <- pformat(freq$p.value)
  if (method == 'pearson'){
    BF <-bfformat(extractBF(correlationBF(x,y))$bf, maxValue = maxBF)  
    output <- paste0('$r=',r,'$, $p',p,'$, $BF',BF,'$')
  }else if (method=='kendall'){
    BF <- bfformat(bfCorrieKernelKendallTau(tau = freq$estimate, n = length(x))$bf10, maxValue = maxBF)
    output <- paste0('$r_{\tau}=',r,'$, $p',p,'$, $BF',BF,'$')
  }
  return(output)
}


#Mixed effects regression: reports the standardized fixed effect size, t statistic, dof, and p value. For BF use bridge sampling (not yet implemented)
lmPretty <- function(mod) { #returns normalized regression coefficient
  b <- fixef(mod)[-1] #extract regression coefficient
  sd.x <- apply(matrix(getME(mod,"X")[,-1]),2,sd)
  sd.y <- sd(getME(mod,"y"))
  normalizedBeta <- betaformat(b*sd.x/sd.y) #normalized beta
  t <- sprintf("%.1f", coef(summary(mod))[,"t value"][-1]) #t statistic
  p <- pformat(coef(summary(mod))[,"Pr(>|t|)"][-1]) #pvalue
  dof <-  sprintf("%.0f",coef(summary(mod))[,"df"][-1])
  return(paste0('$\beta=',normalizedBeta, '$, $t(',dof,')=',t,'$, $p',p,'$')) #No simple way to include bayes factor here. Use bridge sampling instead
}


####################################################################################
# Formatting for different test outputs
####################################################################################

#Pearson's R, Kendall's Tau, or Spearman's Rho: Trims leading zero and displays 2 significant digits
corformat <- function(val) { 
  out <- sub("^(-?)0.", "\\1.", sprintf("%.2f", val)) 
  return(out)
}

#p-values: 3 decimal places of precision; <.001 if the case
pformat <- function(val) { 
  if(val <.001){
    out <- '<.001'
  }else{
    out<-paste0('=',sub("^(-?)0.", "\\1.", sprintf("%.3f", val)))
  }
  return(out)
}


#Cohen's D to a single decimal
dFormat <- function(val) { 
  if (val>=.1){
    paste0('=',sprintf("%.1f", val))
  }else if (val>=.01){
    paste0('=',sprintf("%.2f", val))
  }
  else{
    paste0('=',sprintf("%.2f", val))
  }
}


#Bayes Factor: upper ceiling of 100 as >100, if greater than 1 then display as int, if less than zero show 2 decimal places
bfformat <- function(val, maxValue = 100){
  if (is.na(val)){
    out<- "=NA"
  }else if (val>maxValue){ #If larger than max value, just go with BF>maxValue
    out <- paste0('>',maxValue)
  }else if(val>10){ #if larger than 10 then don't add any decimals
    out<-paste0('=',sprintf("%.0f", val))
  }else if(val>1 & val<10){ #if between 1 and 10 add decimals
    out<-paste0('=',sprintf("%.1f", val))
  }else{ #less than 1, add 2 decimals
    out<-paste0('=',sub("^(-?)0.", "\\1.", sprintf("%.2f", val)))
  }
  return(out)
}

#Beta coefficients:  Trims leading zero and displays 2 significant digits
betaformat <- function(val) { 
  out <- sub("^(-?)0.", "\\1.", sprintf("%.2f", val)) 
  return(out)
}


#################################################################################################
#####   Bayes factor for Kendall's tau, adapted from:                                     
#####   van Doorn, J.B., Ly, A., Marsman, M. & Wagenmakers, E.-J. (2018). Bayesian Inference for Kendall’s Rank Correlation Coefficient. The American Statistician, 72:4, 303-308, DOI: 10.1080/00031305.2016.1264998
#################################################################################################
# Prior specification Kendall's Tau
scaledBetaTau <- function(tau, alpha=1, beta=1){
  result <-   ((pi*2^(-2*alpha))/beta(alpha,alpha))  * cos((pi*tau)/2)^(2*alpha-1)
  return(result)
}

priorTau <- function(tau, kappa){
  scaledBetaTau(tau, alpha = (1/kappa), beta = (1/kappa))
}

priorTauPlus <- function(tau, kappa=1) {
  non.negative.index <- tau >=0
  less.than.one.index <- tau <=1
  value.index <- as.logical(non.negative.index*less.than.one.index)
  result <- tau*0
  result[value.index] <- 2*priorTau(tau[value.index], kappa)
  return(result)
}

priorTauMin <- function(tau, kappa=1) {
  negative.index <- tau <=0
  greater.than.min.one.index <- tau >= -1
  value.index <- as.logical(negative.index*greater.than.min.one.index)
  result <- tau*0
  result[value.index] <- 2*priorTau(tau[value.index], kappa)
  return(result)
}


# Posterior specification Kendall's Tau
postDensKendallTau <- function(delta,Tstar,n,kappa=1,var=var,test="two-sided"){ 
  if(test == "two-sided"){priorDens <- priorTau(delta,kappa)
  } else if(test == "positive"){priorDens <- priorTauPlus(delta,kappa)
  } else if(test == "negative"){priorDens <- priorTauMin(delta,kappa)}
  priorDens <- priorTau(delta,kappa)
  dens <- dnorm(Tstar,(1.5*delta*sqrt(n)),sd=sqrt(var))* priorDens
  return(dens)
}
posteriorTau <- function(delta,kentau,n,kappa=1,var=1,test="two-sided"){
  Tstar <- (kentau * ((n*(n-1))/2))/sqrt(n*(n-1)*(2*n+5)/18)
  var <- min(1,var)
  if(test == "two-sided"){lims <- c(-1,1)
  } else if(test == "positive"){lims <- c(0,1)
  } else if(test == "negative"){lims <- c(-1,0)}
  logicalCensor <- (delta >= lims[1] & delta <= lims[2])
  dens <- logicalCensor*postDensKendallTau(delta,Tstar,n,kappa,var,test=test)/
    integrate(function(delta){postDensKendallTau(delta,Tstar,n,kappa,var,test=test)},lims[1],lims[2])$value
} 

# Bayes factor computation Kendall's Tau
bfCorrieKernelKendallTau <- function(tau, n, kappa=1, var=1, ciValue=0.95){ 
  tempList <- list(vector())
  output <- list(n=n, r=tau, bf10=NA, bfPlus0=NA, bfMin0=NA)
  output$bf10 <- priorTau(0,kappa)/posteriorTau(0,tau,n,kappa=kappa,var=var,test="two-sided")
  output$bfPlus0 <- priorTauPlus(0,kappa)/posteriorTau(0,tau,n,kappa=kappa,var=var,test="positive")
  output$bfMin0 <- priorTauMin(0,kappa)/posteriorTau(0,tau,n,kappa=kappa,var=var,test="negative")
  return(output)
}

# Compute credible intervals kendalls tau
credibleIntervalKendallTau <- function(kentau,n,kappa=1,var=1, test="two-sided", ciValue = 0.95){
  nSeqs <- 1000
  lowCI <- (1-ciValue)/2
  upCI <- (1+ciValue)/2
  taus <- seq(-1,1,length.out = (nSeqs-1))
  densVals <- posteriorTau(taus, kentau, n, kappa = kappa, var = var, test = test)
  densVals <- cumsum((densVals[1:(nSeqs-1)]+densVals[2:nSeqs])*0.5*(taus[2]-taus[1]))
  lowerCI <- taus[which(densVals>=lowCI)[1]]
  upperCI <- taus[which(densVals>=upCI)[1]]
  median <- taus[which(densVals>=0.5)[1]]
  return(list(lowerCI = lowerCI, median = median, upperCI = upperCI))
}

sampleTausA <- function(myTau,myN,nSamples = 3e3, var = 1){
  nSeqs <- 1000
  tauSamples <- NULL
  taus <- seq(-1,1,length.out = nSeqs)
  densVals <- posteriorTau(taus, myTau, myN, var = var)
  ceiling <- max(densVals)
  lowerB <- taus[which(round(densVals,digits=6) != 0 )][1]
  upperB <- rev(taus[which(round(densVals,digits=6) != 0 )])[1]
  
  while(length(tauSamples) < nSamples){
    prop <- runif(1,lowerB,upperB)
    propDens <- posteriorTau(prop, myTau, myN, var = var)
    if(propDens > runif(1,0,ceiling)){tauSamples <- c(tauSamples,prop)}
  }
  return(tauSamples)
}


##############################################################################################################################################################################################################################
#Bayes factor for rank tests
#van Doorn, J., Ly, A., Marsman, M., & Wagenmakers, E. J. (2017). Bayesian Latent-Normal Inference for the Rank Sum Test, the Signed Rank Test, and Spearman's $\rho$. arXiv preprint arXiv:1712.06941.
##############################################################################################################################################################################################################################
rankSumGibbsSampler <- function(xVals,yVals, nSamples = 10000, progBar = TRUE, cauchyPriorParameter = 1/sqrt(2), 
                                nGibbsIterations = 10){
  
  if (progBar) {
    myBar <- txtProgressBar(min = 1, max = nSamples, initial = 1, char = "*",style=3,width=50)
  }
  
  n1 <- length(xVals)
  n2 <- length(yVals)
  
  allRanks <- rank(c(xVals,yVals))
  xRanks <- allRanks[1:n1]
  yRanks <- allRanks[(n1+1):(n1+n2)]
  
  currentVals <- sort(rnorm((n1+n2)))[allRanks] # initial values
  
  deltaSamples <- gSamples <- numeric(nSamples)
  sampledX <- matrix(nrow = nSamples, ncol = n1)
  sampledY <- matrix(nrow = nSamples, ncol = n2)
  
  oldDeltaProp <- 0
  
  for (j in 1:nSamples) {
    
    for (i in sample(1:(n1+n2))) {
      
      currentRank <- allRanks[i]
      
      currentBounds <- upperLowerTruncation(ranks=allRanks, values=currentVals, currentRank=currentRank)
      if (i <= n1) {
        oldDeltaProp <- -0.5*oldDeltaProp
      } else if (i > n1) {
        oldDeltaProp <- 0.5*oldDeltaProp
      }
      
      currentVals[i] <- myTruncNormSim(currentBounds[["under"]], currentBounds[["upper"]], mu=oldDeltaProp, sd=1)
      
    }
    
    xVals <- currentVals[1:n1]
    yVals <- currentVals[(n1+1):(n1+n2)]
    
    gibbsResult <- sampleGibbsTwoSample(x = xVals, y = yVals, n1 = n1, n2 = n2, nIter = nGibbsIterations,
                                        rscale = cauchyPriorParameter)
    
    deltaSamples[j] <- oldDeltaProp <- gibbsResult[1]
    gSamples[j] <- gibbsResult[2]
    
    if (progBar) setTxtProgressBar(myBar,j) 
    sampledX[j,] <- xVals
    sampledY[j,] <- yVals
  }
  resultsList <- list(deltaSamples = deltaSamples, gSamples = gSamples,
                      sampledX = sampledX, sampledY = sampledY)
  return(resultsList)
}

sampleGibbsTwoSample <- function(x, y, n1, n2, nIter = 10, rscale = 1/sqrt(2)) {
  
  meanx <- mean(x)
  meany <- mean(y)
  n1 <- length(x)
  n2 <- length(y)
  sigmaSq <- 1 # Arbitrary number for sigma
  g <- 1
  
  for(i in 1:nIter){   
    #sample mu
    varMu <- (4 * g * sigmaSq) / ( 4 + g * (n1 + n2) )
    meanMu <- (2 * g * (n2 * meany - n1 * meanx)) / ((g * (n1 + n2) + 4))
    mu <- rnorm(1, meanMu, sqrt(varMu)) 
    # sample g
    betaG <- (mu^2 + sigmaSq * rscale^2) / (2*sigmaSq)
    g <- 1/rgamma(1, 1, betaG)
    # convert to delta
    delta <- mu / sqrt(sigmaSq)
  }
  
  return(c(delta, mu, g))
  
}


upperLowerTruncation <- function(ranks, values, currentRank, n, ranksAreIndices = FALSE) {
  
  if (currentRank == min(ranks)) {
    under <- -Inf
  } else {
    under <- max(values[ranks < currentRank])
  }
  
  if (currentRank == max(ranks)) {
    upper <- Inf
  } else {
    upper <- min(values[ranks > currentRank])
  }
  
  return(list(under=under, upper=upper))
}

myTruncNormSim <- function(lBound = -Inf, uBound = Inf, mu = 0, sd = 1){
  
  lBoundUni <- pnorm(lBound, mean = mu, sd = sd)
  uBoundUni <- pnorm(uBound, mean = mu, sd = sd)  
  mySample <- qnorm(runif(1, lBoundUni, uBoundUni), mean = mu, sd = sd)
  
  return(mySample)
}

signRankGibbsSampler <- function(xVals, yVals = NULL, testValue = 0, nSamples = 10000, progBar = FALSE, cauchyPriorParameter = 1/sqrt(2),
                                 nGibbsIterations = 5, varyTies = TRUE) {
  
  if (progBar) {
    myBar <- txtProgressBar(min = 1, max = nSamples, initial = 1, char = "*",style=3,width=50)
  }
  
  n <- length(xVals)
  nDF <- n-1
  
  if (!is.null(yVals)) { 
    differenceScores <- xVals - yVals
  } else {rankSumGibbsSampler <- function(xVals,yVals, nSamples = 10000, progBar = TRUE, cauchyPriorParameter = 1/sqrt(2), 
                                          nGibbsIterations = 10){
    
    if (progBar) {
      myBar <- txtProgressBar(min = 1, max = nSamples, initial = 1, char = "*",style=3,width=50)
    }
    
    n1 <- length(xVals)
    n2 <- length(yVals)
    
    allRanks <- rank(c(xVals,yVals))
    xRanks <- allRanks[1:n1]
    yRanks <- allRanks[(n1+1):(n1+n2)]
    
    currentVals <- sort(rnorm((n1+n2)))[allRanks] # initial values
    
    deltaSamples <- gSamples <- numeric(nSamples)
    sampledX <- matrix(nrow = nSamples, ncol = n1)
    sampledY <- matrix(nrow = nSamples, ncol = n2)
    
    oldDeltaProp <- 0
    
    for (j in 1:nSamples) {
      
      for (i in sample(1:(n1+n2))) {
        
        currentRank <- allRanks[i]
        
        currentBounds <- upperLowerTruncation(ranks=allRanks, values=currentVals, currentRank=currentRank)
        if (i <= n1) {
          oldDeltaProp <- -0.5*oldDeltaProp
        } else if (i > n1) {
          oldDeltaProp <- 0.5*oldDeltaProp
        }
        
        currentVals[i] <- myTruncNormSim(currentBounds[["under"]], currentBounds[["upper"]], mu=oldDeltaProp, sd=1)
        
      }
      
      xVals <- currentVals[1:n1]
      yVals <- currentVals[(n1+1):(n1+n2)]
      
      gibbsResult <- sampleGibbsTwoSample(x = xVals, y = yVals, n1 = n1, n2 = n2, nIter = nGibbsIterations,
                                          rscale = cauchyPriorParameter)
      
      deltaSamples[j] <- oldDeltaProp <- gibbsResult[1]
      gSamples[j] <- gibbsResult[2]
      
      if (progBar) setTxtProgressBar(myBar,j) 
      sampledX[j,] <- xVals
      sampledY[j,] <- yVals
    }
    resultsList <- list(deltaSamples = deltaSamples, gSamples = gSamples,
                        sampledX = sampledX, sampledY = sampledY)
    return(resultsList)
  }
  
  sampleGibbsTwoSample <- function(x, y, n1, n2, nIter = 10, rscale = 1/sqrt(2)) {
    
    meanx <- mean(x)
    meany <- mean(y)
    n1 <- length(x)
    n2 <- length(y)
    sigmaSq <- 1 # Arbitrary number for sigma
    g <- 1
    
    for(i in 1:nIter){   
      #sample mu
      varMu <- (4 * g * sigmaSq) / ( 4 + g * (n1 + n2) )
      meanMu <- (2 * g * (n2 * meany - n1 * meanx)) / ((g * (n1 + n2) + 4))
      mu <- rnorm(1, meanMu, sqrt(varMu)) 
      # sample g
      betaG <- (mu^2 + sigmaSq * rscale^2) / (2*sigmaSq)
      g <- 1/rgamma(1, 1, betaG)
      # convert to delta
      delta <- mu / sqrt(sigmaSq)
    }
    
    return(c(delta, mu, g))
    
  }
  
  
  upperLowerTruncation <- function(ranks, values, currentRank, n, ranksAreIndices = FALSE) {
    
    if (currentRank == min(ranks)) {
      under <- -Inf
    } else {
      under <- max(values[ranks < currentRank])
    }
    
    if (currentRank == max(ranks)) {
      upper <- Inf
    } else {
      upper <- min(values[ranks > currentRank])
    }
    
    return(list(under=under, upper=upper))
  }
  
  myTruncNormSim <- function(lBound = -Inf, uBound = Inf, mu = 0, sd = 1){
    
    lBoundUni <- pnorm(lBound, mean = mu, sd = sd)
    uBoundUni <- pnorm(uBound, mean = mu, sd = sd)  
    mySample <- qnorm(runif(1, lBoundUni, uBoundUni), mean = mu, sd = sd)
    
    return(mySample)
  }
  
  
  differenceScores <- xVals - testValue
  }
  
  differenceSigns <- (sign(differenceScores))
  absDifferenceRanked <- rank(abs(differenceScores))
  prodSignAbsRank <- differenceSigns * absDifferenceRanked
  
  initDiffSamples <- sort(abs(rnorm(n)))[absDifferenceRanked]
  sampledDiffsAbs <- abs(initDiffSamples)
  diffSamples <- numeric(n)
  
  deltaSamples <- gSamples <- numeric(nSamples)
  diffSamplesMatrix <- matrix(nrow=nSamples,ncol=n)
  
  oldDeltaProp <- 0
  
  for (j in 1:nSamples) {
    for (i in sample(n)) {
      
      currentRank <- absDifferenceRanked[i]
      
      currentBounds <- upperLowerTruncation(ranks=absDifferenceRanked, values=sampledDiffsAbs, currentRank=currentRank)
      if (is.infinite(currentBounds[["under"]])) {currentBounds[["under"]] <- 0}
      
      sampledDiffsAbs[i] <- myTruncNormSim(currentBounds[["under"]], currentBounds[["upper"]], 
                                           mu = abs(oldDeltaProp), sd=1)
    }
    
    diffSamples <- sampledDiffsAbs * differenceSigns
    
    if (any(differenceSigns == 0) && varyTies) {
      nullSamples <- sampledDiffsAbs[differenceSigns == 0] * sample(c(-1,1), size = sum(differenceSigns == 0), replace = T)
      diffSamples[which(differenceSigns == 0)] <- nullSamples
    }
    
    sampledDiffsAbs <- abs(diffSamples)
    diffSamples <- diffSamples
    
    gibbsOutput <- sampleGibbsOneSample(diffScores = diffSamples, nIter = nGibbsIterations, rscale = cauchyPriorParameter)
    
    deltaSamples[j] <- oldDeltaProp <- gibbsOutput[1]
    gSamples[j] <- gibbsOutput[2]
    
    diffSamplesMatrix[j,] <- diffSamples
    
    if(progBar)setTxtProgressBar(myBar,j) # update the progressbar
  }
  
  resultsList <- list(deltaSamples = deltaSamples, gSamples = gSamples, diffSamples = diffSamplesMatrix)
  return(resultsList)
}

sampleGibbsOneSample <- function(diffScores, nIter = 10, rscale = 1/sqrt(2)){
  
  ybar <- mean(diffScores)
  n <- length(diffScores)
  sigmaSq <- 1
  mu <- ybar
  g <- ybar^2 / sigmaSq + 1
  
  for(i in 1:nIter){   
    #sample mu
    varMu  = sigmaSq / (n + (1 / g))
    meanMu = (n * ybar) / (n + (1 / g))
    mu <- rnorm(1, meanMu, sqrt(varMu) )
    
    # sample g
    scaleg = (mu^2 + sigmaSq * rscale^2) / (2*sigmaSq)
    g = 1 / rgamma(1, 1, scaleg )
    
    delta <- mu / sqrt(sigmaSq)
  }
  return(c(delta, g))
}



upperLowerTruncation <- function(ranks, values, currentRank, n, ranksAreIndices = FALSE) {
  
  if (currentRank == min(ranks)) {
    under <- -Inf
  } else {
    under <- max(values[ranks < currentRank])
  }
  
  if (currentRank == max(ranks)) {
    upper <- Inf
  } else {
    upper <- min(values[ranks > currentRank])
  }
  
  return(list(under=under, upper=upper))
}

myTruncNormSim <- function(lBound = -Inf, uBound = Inf, mu = 0, sd = 1){
  
  lBoundUni <- pnorm(lBound, mean = mu, sd = sd)
  uBoundUni <- pnorm(uBound, mean = mu, sd = sd)  
  mySample <- qnorm(runif(1, lBoundUni, uBoundUni), mean = mu, sd = sd)
  
  return(mySample)
}
