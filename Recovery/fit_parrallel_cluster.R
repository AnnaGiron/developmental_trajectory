
###
###This script contains modified versions of the modelsfitting functions to run on the cluster
###

modelFit<-function(par, subjD, acquisition, k,  horizonLength, rounds){
  #Extract and process parameters
  if (inherits(acquisition, "epsilonGreedy")){
    epsilon<- 1/(1+exp(-(par[length(par)]))) #transform back from unbounded space; epsilon is the last parameter for epsilon greedy
  }
  par<-exp(par) #exponentiate parameters to make a non-negative and convex optimization surface
  #last parameter for all other models is always inverse temperature for softmax
  tau<-par[length(par)]
  #Which posterior function to use; therefore, which parameters to use
  if (inherits(k, "KalmanFilter")){ #null kernel indicates kalman filter model
    kNoise <- par[1]
    parVec <- c(kNoise) #Vector of parameters to send to the KF posterior function
  }else if(inherits(k, "GP")){ #lambda
    lambda <- par[1]
    parVec <- c(lambda, lambda, 1, .0001) # Vector of parameters to send to the GP posterior vector, where sF and sN are fixed
  }
  #Additional acquisition function dependent parameters
  if (inherits(acquisition, "UCB")| inherits(acquisition, 'exploreCounts')| inherits(acquisition, 'epsilonGreedy')){ #check if UCB is used
    beta <- par[length(par)-1] #If UCB, beta is always 2nd last
    #refactor beta and tau into gamma and beta_star, where gamma = 1/tau and beta_star = beta/tau
  }
  #which rounds to consider?
  trainingSet <- subset(subjD, round %in% rounds)
  #Vector to store negative log likelihods
  nLL <- rep(0,length(rounds))
  for (r in unique(trainingSet$round)){ #Begin looping through each round
    #subset of data for round r
    roundD <- subset(subjD, round==r)
    horizon <- nrow(roundD)
    #Observations of subject choice behavior
    chosen <- roundD$chosen
    chosen <- chosen[2:length(chosen)] # trim first observation, since it wasn't a choice but a randomly revealed tile
    y  <- roundD$z[0:(horizon-1)] #trim off the last observation, because it was not used to inform a choice (round already over)
    x1 <- roundD$x[0:(horizon-1)]
    x2 <- roundD$y[0:(horizon-1)]
    #create observation matrix
    X<-as.matrix(cbind(x1,x2))
    #make sure X is a matrix
    X<-as.matrix(X)
    Xnew<-as.matrix(Xnew)
    #Utilties of each choice
    utilities <- NULL
    prevPost <- NULL #set the previous posterior computation to NULL for the kalman filter
    pMat <- NULL
    #loop through observations
    for (i in 1:(horizon-1)){ #skip the last observation, because no choice was made based on that information
      #new observation
      X1<-matrix(X[1:i,], ncol=2)
      y1<-y[1:i]
      #Which posterior function to use
      if (inherits(k, "KalmanFilter")){# kalman filter model
        out<- bayesianMeanTracker(x = X1[i,], y=y[i], prevPost = prevPost, theta = parVec)
        #update prevPost for the next round
        prevPost <- out
      }else if (inherits(k, 'GP')){# GP with length-scale parameterized kernel
        out <- gpr(X.test=Xnew, theta=parVec, X=X1, Y=y1, k=k) #Mu and Sigma predictions for each of the arms; either GP or Kalman filter
      }else if (inherits(k, 'Null')){ #null model
        out <- nullModel() #Mu and Sigma predictions for each of the arms; either GP or Kalman filter
      }
      #Slightly different function calls for each acquisition function
      if (inherits(acquisition, "UCB")){ #UCB takes a beta parameter
        utilityVec<-acquisition(out, c(beta))
      } else if (inherits(acquisition, 'exploreCounts')){ #count-based exploration
        utilityVec <- exploreCounts(out, roundD$chosen[1:i], c(beta))
      }else if (inherits(acquisition, "epsilonGreedy")){
        p <- epsilonGreedy(out, beta, epsilon)
        pMat <- rbind(pMat, t(p))
      }else{ #any other
        utilityVec <- acquisition(out)
      }
      if (inherits(acquisition, "softmax")){
        utilityVec <- utilityVec - max(utilityVec) #avoid overflow
        utilities <- rbind(utilities, t(utilityVec)) # build horizon_length x options matrix, where each row holds the utilities of each choice at each decision time in the search horizon
      }
    }
    #print(utilities)
    if (inherits(acquisition, "softmax")){
      #Softmax rule
      p <- exp(utilities/tau)
      p <- p/rowSums(p)
      #avoid underflow by setting a floor and a ceiling
      p <- (pmax(p, 0.00001))
      p <- (pmin(p, 0.99999))
      pMat<- p
    }
    #Calculate Negative log likelihood
    nLL[which(unique(trainingSet$round)==r)] <- -sum(log(pMat[cbind(c(1:(horizon-1)),chosen)]))
  }
  #end loop through rounds
  return(sum(nLL))  #Return negative log likelihoods of all observations
}




####
####
####


modelSimulateChoice<-function(par, rounds,k, acquisition,env,envSelector){
  #enselect<-sample(1:40, 1)
  environ<-subset(env, en==envSelector[1])
  #setup stuff i need for later
  nTrials = 25
  Xstar<-cbind(environ$x1[1:64], environ$x2[1:64])#expand.grid(0:7,0:7)#
  #cbind(environ$x1[1:64], environ$x2[1:64])# just sets up x,y positions in a kinda complicated way
  dparams= NULL
  browser()

  par<-exp(par) #exponentiate parameters to make a non-negative and convex optimization surface->SC: seems wierd to me.
  #last parameter is always inverse temperature for softmax
  tau<-par[length(par)]
  beta=NA
  lambda=NA
  epsilon=NA
  
  if (inherits(acquisition, "epsilonGreedy")){
    epsilon<- 1/(1+exp(-(par[length(par)]))) #transform back from unbounded space; epsilon is the last parameter for epsilon greedy
  }
  #Which posterior function to use; therefore, which parameters to use
  if (inherits(k, "KalmanFilter")){ #null kernel indicates kalman filter model
    kNoise <- par[1]
    parVec <- c(kNoise) #Vector of parameters to send to the KF posterior function
  }else if(inherits(k, "GP")){ #lambda
    lambda <- par[1]
    parVec <- c(lambda, lambda, 1, .0001) # Vector of parameters to send to the GP posterior vector, where sF and sN are fixed
  }
  #Additional acquisition function dependent parameters
  if (inherits(acquisition, "UCB")| inherits(acquisition, 'exploreCounts') | inherits(acquisition, "epsilonGreedy") ){ #check if UCB is used
    beta <- par[length(par)-1] #If UCB, beta is always 2nd last
    #refactor beta and tau into gamma and beta_star, where gamma = 1/tau and beta_star = beta/tau
  }
  # browser()
  for (round in rounds){
    # enselect<-sample(1:40, 1)
    prevPost=NULL# for kalman filter
    pMat <- NULL
    utilities <- NULL
    
    environ<-subset(env, en==envSelector[round])
    ind<-sample(1:64,1)
    #X matrix
    X<-cbind(environ$x1[ind], environ$x2[ind])
    #y matrix
    y<-as.matrix(environ$y[ind])
    choice<-ind
    #loop through trials
    #browser()
    for (trial in 1:nTrials){
      X1<-matrix(X[1:trial,], ncol=2)
      y1<-y[1:trial]
      #output by GP with particular parameter settings
      #don't forget mean centering and standardization
      if (inherits(k, "KalmanFilter")){# kalman filter model
        out<- bayesianMeanTracker(x = X1[trial,], y=(y[trial]-25)/50, prevPost = prevPost, theta = parVec)
        #update prevPost for the next round
        prevPost <- out
      }else if (inherits(k, 'GP')){# GP with length-scale parameterized kernel
        out <- gpr(X.test=Xstar, theta=parVec, X=X1, Y=(y1-25)/50, k=k) #Mu and Sigma predictions for each of the arms; either GP or Kalman filter
      }else if (inherits(k, 'Null')){ #null model
        out <- nullModel() #Mu and Sigma predictions for each of the arms; either GP or Kalman filter
      }
      #browser()
      
      #Slightly different function calls for each acquisition function
      if (inherits(acquisition, "UCB")){ #UCB takes a beta parameter
        utilityVec<-acquisition(out, c(beta))
      } else if (inherits(acquisition, 'exploreCounts')){ #count-based exploration
        utilityVec <- exploreCounts(out, roundD$chosen[1:i], c(beta))
      }else if (inherits(acquisition, "epsilonGreedy")){
        p <- epsilonGreedy(out, beta, epsilon)
        pMat <- rbind(pMat, t(p))
      }else{ #any other
        utilityVec <- acquisition(out)
      }
      if (inherits(acquisition, "softmax")){
        utilityVec <- utilityVec - max(utilityVec) #avoid overflow
        utilities <- rbind(utilities, t(utilityVec)) # build horizon_length x options matrix, where each row holds the utilities of each choice at each decision time in the search horizon
      }
      #print(utilities)
      
      if (inherits(acquisition, "softmax")){
        #Softmax rule
        p <- exp(utilities/tau)
        p <- p/rowSums(p)
        #avoid underflow by setting a floor and a ceiling
        p <- (pmax(p, 0.00001))
        p <- (pmin(p, 0.99999))
        pMat<- p
      }
      #browser()
      #index is sampled proportionally to softmaxed utility vector
      ind<-sample(1:64, 1, prob=pMat[trial,])
      #bind X-observations
      X<-rbind(X, cbind(environ$x1[ind], environ$x2[ind]))
      #bind y-observations
      y<-rbind(y, as.matrix(environ$y[ind]))
      choice<-rbind(choice,ind)
    }
    # mu<-c(mu, mean(y))
    cur = data.frame(round=round, 
                     lambda=rep(lambda, nTrials+1), 
                     beta=rep(beta, nTrials+1), 
                     tau=rep(tau, nTrials+1),
                     epsilon=rep(epsilon, nTrials+1),
                     mu=y, 
                     x1=X[,1], 
                     x2=X[,2],
                     chosen=choice)
    dparams = rbind(dparams, cur)
  }
  # dparams$mu[i] = mean(mu)
  dparams
}
#rn negative log likelihoods of all observations



FitParamCV<-function(d1, kernelFun, acquisition, leaveoutindex){
  library(DEoptim)
  
  #subselect participant, horizon and rounds not left out
  #d1<-subset(data, id==uid)
  #training set
  rounds <- 2:9
  #rounds <- count(d1$round)[count(d1$round)$freq>=2,]$x #only rounds where at least 2 clicks have been made
  trainingSet <- rounds[! rounds==leaveoutindex] #remove round specified by leaveoutindex
  #test set
  testSet <- leaveoutindex
  nParams <- 1
  if (inherits(acquisition, 'UCB')| inherits(acquisition, 'exploreCounts')| inherits(acquisition, 'epsilonGreedy')){
    nParams <- nParams + 1 #add beta parameter
  }
  if (inherits(kernelFun, 'GP') | inherits(kernelFun, 'KalmanFilter')){
    nParams <- nParams + 1 #add lambda or error variance
  }
  
  
  #Set upper and lower bounds based on nParams
  lbound <- rep(-5, nParams)
  ubound <- rep(4, nParams)
  if (inherits(acquisition, 'epsilonGreedy') ){ #separate range for epsilon greedy, which after inverse logic transform is bounded between [0,1]
    lbound[length(lbound)] <- -10
    ubound[length(ubound)] <- 10
  }
  #Begin cross validation routine
  if (nParams>=2){#if 2 or more parameters
    #TRAINING SET
    fit<-DEoptim(modelFit, lower=lbound, upper=ubound, subjD=d1, k=kernelFun, rounds = trainingSet, acquisition=acquisition, DEoptim.control(itermax=100))
    paramEstimates <- fit$optim$bestmem #MODEL DEPENDENT PARAMETER ESTIMATES
    #TEST SET
    predict <- modelFit(par=paramEstimates, subjD=d1, acquisition=acquisition, k=kernelFun, rounds=testSet)
    output <- c(leaveoutindex, predict, fit$optim$bestmem) # leaveoutindex, nLL, parameters....
  } else{
    #TRAINING SET
    fit<-optimize(modelFit, lower=lbound, upper=ubound, subjD=d1, k=kernelFun, rounds = trainingSet, acquisition=acquisition)
    paramEstimates <- fit$minimum #MODEL DEPENDENT PARAMETER ESTIMATES
    #TEST SET
    predict <- modelFit(par=paramEstimates, subjD=d1, acquisition=acquisition, k=kernelFun, rounds=testSet)
    output <- c(leaveoutindex, predict, fit$minimum) #leaveoutindex, nLL, parameters....
  }
  return(output) #return optimized value
}
