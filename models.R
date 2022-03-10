#Models and sampling strategies
#Charley Wu Feb 2022
library(nnet)
##############################################################################################################
#GP KERNEL
##############################################################################################################

#Radial Basis Kernel
rbf <- function(X1,X2,theta){
  #transfer to matrices
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  #check dimensions
  if(ncol(X1) != ncol(X2)){
    stop("X1 and X2 must contain input values of the same dimension.")
  } else if(!all(theta>=0)){
    stop("All parameters must be >= 0.")
  }
  #get dimensions
  N1 <- nrow(X1)
  N2 <- nrow(X2)
  d <- ncol(X1)
  #initialize sigma
  sigma <-  matrix(rep(0, N1*N2),nrow=N1)
  #observational variance
  sf <- theta[d+1]
  #noise variance
  sn <- theta[d+2]
  #loop through
  for(i in 1:d){
    #length scale
    l <- theta[i] #Note: assumes a unique length scale for each dimension
    #x-diff
    xdiff <- (outer(X1[,i],X2[,i],function(x,y) x - y)/l)^2
    sigma <- sigma + xdiff
  }
  #RBF function
  if(identical(X1,X2)){
    id <- diag(rep(1,N1))
    sigma.final <- sf*exp(-0.5*sigma) + sn*id
  } else {
    sigma.final <- sf*exp(-0.5*sigma)
  }
  #return final covariance matrix
  return(sigma.final)
}
class(rbf)<- c(class(rbf), "GP") #identify the rbf kernel as a gp model


##############################################################################################################
#MATRIX INVERSION
##############################################################################################################

#calculate inverse of the cov-function using sigular value decomposition
cov.inverse.svd <- function(X, tol = sqrt(.Machine$double.eps)){
  # Generalized Inverse of a Matrix
  dnx <- dimnames(X)
  if(is.null(dnx)) dnx <- vector("list", 2)
  #singular value decomposition
  s <- svd(X)
  nz <- s$d > tol * s$d[1]
  #inverse
  K.inv <- structure(
    if(any(nz)) s$v[, nz] %*% (t(s$u[, nz])/s$d[nz]) else X,
    dimnames = dnx[2:1])
  #logarithm of determinant.
  log.K.det <- sum(log(s$d))
  #return inverse plus log determinant
  return(list(Inv = K.inv,lDet = log.K.det))
}

#calculate inverse of the cov-function using Cholesky
cov.inverse.chol <- function(X){
  #cholseky decomposition
  R <- chol(X)
  #complex conjugate
  Rt <- Conj(t(R))
  #invert
  R.inv <- solve(R)
  #invert
  Rt.inv <- solve(Rt)
  #multiply matrices
  X.inv <- R.inv %*% Rt.inv
  #log determinant
  log.X.det <- 2*sum(log(diag(R))) 
  #return both
  return(list(Inv = X.inv, lDet = log.X.det))
}

##############################################################################################################
#GAUSSIAN PROCESS
##############################################################################################################

#Gaussian Process function
#X.test: matrix for predcitions
#theta: vector of hyper-parameters (lambda, Sf, Sn)
#X; matrix of observations
#y: vector of observed outcomes
#kernel: used kernel function, can be "rbf", "oru", or "mat"
gpr <- function(X.test, theta, X,Y, k){
  #make it a matrix
  Xstar <- as.matrix(X.test)
  #dimensions
  d <- ncol(as.matrix(X))
  #calculate capital K
  K <- k(X,X,theta) 
  #Check if matrix is positive semi-definite
  #browser()
  if (is.positive.definite(K)){
    #KK <- cov.inverse.chol(K) #use Cholesky
    KK.inv <- chol2inv(chol(K)) #MASS implementation of Cholesky
  } else {
    KK.inv <- cov.inverse.svd(K) #use SVD
  }
  #times y
  Ky <- KK.inv %*% Y
  #apply the kernel
  result <- apply(Xstar, 1, function(x){
    XX <- matrix(x,nrow=1) 
    Kstar <- k(X, XX, theta)
    Kstarstar <- k(XX,XX,theta)
    #get mean vector
    mu <- t(Kstar) %*% Ky
    #get covariance
    cv <- Kstarstar - (t(Kstar) %*% KK.inv %*% Kstar) #BUG: sometimes cv<0, leading to NaN when we return sqrt(cv)
    #DEBUG
    if (cv<0){ cv <- abs(cv)} #TEMPORARY SOLUTION: MANUALLY SET CV TO BE POSITIVE IF NEGATIVE
    #return means and variance
    return(c(mu, cv))
  })
  #as a data frame with names mu and sig
  prediction <- as.data.frame(t(result))
  prediction[is.na(prediction)] <- 0.01 #remove NaN items with the noise variance 
  colnames(prediction) <- c("mu", "sig")
  return(prediction)
}



##############################################################################################################
##Mean Tracker
##############################################################################################################
bayesianMeanTracker <- function(x, y, theta=c(1), prevPost=NULL){ #Updates the previous posterior based on a single observation
  #parameters
  mu0 <- 0 #prior mean
  var0 <- 5 #prior variance
  vare <- theta[1] #error varriance
  if (is.null(prevPost)){#if no posterior prior, assume it is the first observation
    predictions <- data.frame(mu=rep(mu0,64), sig=rep(var0,64))
  }else{#if previous posterior is provided, update
    predictions <- prevPost
  }
  #Which of the 121 options were chosen at time?
  allopts<-expand.grid(0:7, 0:7)
  chosen <- which(allopts$Var1==x[1] & allopts$Var2==x[2])
  #Kalman gain
  kGain <- predictions$sig[chosen] / (predictions$sig[chosen] + vare^2)
  #update mean
  predictions$mu[chosen] <- predictions$mu[chosen] + (kGain * (y-predictions$mu[chosen]))
  #update variance for observed arm
  predictions$sig[chosen] <- predictions$sig[chosen] * (1 - kGain)
  #return output
  return(predictions)
}
class(bayesianMeanTracker)<- c(class(bayesianMeanTracker), "KalmanFilter")


##############################################################################################################
#NULL MODEL (used with heuristics)
##############################################################################################################
#Returns a mean of 1 and variance of 0 for all models
nullModel <- function(..){ #
  predictions <- data.frame(mu=rep(1,64), sig=rep(0,64)) #defaults all means to zero
  return(predictions)
}
class(nullModel)<- c(class(nullModel), "Null")


##############################################################################################################
#ACQUISITION FUNCTIONS
##############################################################################################################

#Upper Confidence Bound Sampling
ucb<-function(out, pars, refactor=F){
  if (refactor==TRUE){
    gamma <- pars[1]
    beta_star<-pars[2]
    #calulate all the upper confidence bounds
    outtotal<-(gamma*out$mu)+(beta_star*sqrt(out$sig)) #refactored parameters in combination with softmax tau, where gamma = 1/tau and beta_star = beta/tau
    #avoid borderline cases
    outtotal[outtotal<=0]<-0.0001
    outtotal[outtotal>100]<-100
    outtotal<-matrix(outtotal, ncol=nrow(out)/64, byrow=TRUE)
  }else{
    beta <- pars[1]
    #calulate all the upper confidence bounds
    outtotal<-out$mu+(beta*sqrt(out$sig)) #refactored parameters in combination with softmax tau, where gamma = 1/tau and beta_star = beta/tau
    #avoid borderline cases
    outtotal[outtotal<=0]<-0.0001
    outtotal[outtotal>50]<-50
    outtotal<-matrix(outtotal, ncol=nrow(out)/64, byrow=TRUE)
  }
  #return them
  return(outtotal)
}
#add "UCB" to the class of the ucb function, so that modelFit can recognize that it has a longer parameter array
class(ucb)<- c(class(ucb), "UCB")


#The exploration count model weights each option by the inverse of the visitation counts, such that unvisited options are weighted higher
exploreCounts<-function(out, prev, pars){
  beta <- pars[1]
  visits <- as.data.frame(table(prev)) #calculate the number of visits to each optoin
  out$sig <- 1 #replace with ones as the default
  out[as.numeric(as.character(visits$prev)),'sig'] <- 1 /(visits$Freq+1) #replace the uncertainty estimates with the inverse of visit frequency
  outtotal<-out$mu+(beta*(out$sig)) #refactored parameters in combination with softmax tau, where gamma = 1/tau and beta_star = beta/tau
  #avoid borderline cases
  outtotal[outtotal<=0]<-0.0001
  outtotal[outtotal>50]<-50
  return(outtotal)
}
class(exploreCounts) <- c(class(exploreCounts), "exploreCounts")

#Greedy Mean
greedyMean <- function(out){
  outtotal<-out$mu #the value of each choice is solely based on the expectation of reward
  #avoid borderline cases
  outtotal[outtotal<0]<-0.0001
  outtotal[outtotal>50]<-50
  outtotal<-matrix(outtotal, nrow=nrow(out)/64, byrow=TRUE)
  return(t(outtotal))
}
class(greedyMean)<- c(class(greedyMean), "greedyMean")

#Greedy Variance
greedyVar <- function(out){
  outtotal<- sqrt(out$sig) #the value of each choice is solely based on the expected uncertainty
  #avoid borderline cases
  outtotal[outtotal<0]<-0.0001
  outtotal[outtotal>50]<-50
  outtotal<-matrix(outtotal, nrow=nrow(out)/64, byrow=TRUE)
  return(t(outtotal))
}
class(greedyVar)<- c(class(greedyVar), "greedyVar")

#Epsilon greedy
epsilonGreedy <- function(out, beta = 0.5, epsilon=.1){
    n <- length(out$mu)
    p <- rep(1/n*epsilon, n)
    uppercb <- out$mu +(beta *sqrt(out$sig))
    p[which.is.max(uppercb)] <- (1-epsilon) + (1/n*epsilon)
    return(p)
}

class(epsilonGreedy)<- c(class(epsilonGreedy), "epsilonGreedy")

#Add softmax class to all acquisition functions that don't directly output probabilities
class(ucb)<- c(class(ucb), "softmax")
class(exploreCounts)<- c(class(exploreCounts), "softmax")
class(greedyMean)<- c(class(greedyMean), "softmax")
class(greedyVar)<- c(class(greedyVar), "softmax")


##############################################################################################################
#INERTIA MODEL
##############################################################################################################
#The inertial model only uses inverse manhattan distance from the RANDOM INITIALIZATION POINT in order to model behavior
#Here we use it like a posterior, where the means are the inverse euclidean distance, and the variance is always zero
#Combined with the UCB model, these posteriors will yield a behavior akin to inertia
#since the variance is always zero, it will be a greedy mean strategy

readBlockDistance <- function(participantid ){#read precomputed .csv files of manhattan block distance for a specific participant
  blocks <- readRDS(paste0('blockDistance/', participantid, '.RDS'))
  return(blocks)
}

inertia<-function(blocks, trialnumber){
  distance <- blocks[,trialnumber]
  predictions <- data.frame(mu=distance, sig = rep(0,64))
  colnames(predictions) <- c("mu", "sig")
  return(predictions)
}
class(inertia) <- c(class(inertia), "inertia")


##############################################################################################################
#WIN STAY LOSE Shift MODEL
##############################################################################################################
#win := reward_t >= reward* ; else lose
#stay := repeat or cardinal neighbors
#go := any unrevealed tile with equal probability
#WSLG takes a vector of previously observed reweards Yt at inputs Xt (matrix), where the last entry in both is the previous action
WSLS <- function(Yt, Xt){
  mu <- rep(0,121) #initialize matrix storing the value of each tile for t+1
  allopts<-expand.grid(0:10, 0:10) #expand matrix of bivariate function to a vector of 121 unique locations
  y_t <- Yt[length(Yt)] #most recent observation
  x_t <- Xt[length(Yt),]
  y_star <- max(Yt[1:(length(Yt)-1)]) #previous best reward
  if (y_t >= y_star){ #WIN
    #STAY
    #loop through all neighboring tiles as well as current tile, and set value to 1
    for(i in 1:-1){
      for(j in 1:-1){
        mu[which(x_t[1]+i==allopts$Var1 & x_t[2]+j==allopts$Var2)] <- 1 
      }
    }
  }else{ #LOSE
    #GO
    revealed <- apply(Xt,1, FUN=function(x) which(x[1]==allopts$Var1 & x[2]==allopts$Var2))
    mu[-unlist(revealed)] <- 1 #assign all unrevealed tiles a value of 1
  }
  predictions <- data.frame(mu=mu, sig = rep(0,121))
  colnames(predictions) <- c("mu", "sig") #WSLG returns the same kind of data frame as GP regression, but the variance is always set to 0
  return(predictions)
}

class(WSLS) <- c(class(WSLS), "WSLS")
