# Script to run rational models in comparison to human behavior
# Charley Wu, Anna Giron 2022

#############################################################################################################################
# IMPORT BEHAVIORAL DATA AND ENVIRONMENTS
#############################################################################################################################
#house keeping
rm(list=ls())

#load packages
packages <- c('fgpt', 'plyr', 'jsonlite', 'gridExtra', 'reshape2', 'stargazer', "grid", 'matrixcalc', 'parallel')
#load them
invisible(lapply(packages, require, character.only = TRUE))

source("models.R")
source("dataProcessing.R")
#Participant data
d = read.csv('data/behavioralData.csv')

d$agegroup = factor(d$agegroup,
                    levels=c('[25,55]', '[18,25)', '[14,18)', '[11,14)', '[9,11)', '[7,9)', '[5,7)'),
                    labels=c('25-55', '18-24', '14-17', '11-13', '9-10', '7-8', '5-6'))


reps <- 10000 #replications
# reps=2 # for testing
cores <- 2
trials <- 26

#Environments
#load environments from json, unlist, transform into numeric, convert into matrix, and name dimensions
environments <- lapply(fromJSON("data/smoothKernel.json"), FUN=function(x) matrix(as.numeric(unlist(x)), ncol=3, byrow=TRUE, dimnames=list(seq(1,64),  c('x2', 'y', 'x1'))))

#############################################################################################################################
# MODEL FITTING RESTULTS
#############################################################################################################################
# # import modeling results
# params = paramsImport()
# import preprocessed modeling results
params = read.csv('data/modelFit.csv')

params$agegroup = factor(params$agegroup,
                         levels=c('[25,55]', '[18,25)', '[14,18)', '[11,14)', '[9,11)', '[7,9)', '[5,7)'),
                         labels=c('25-55', '18-24', '14-17', '11-13', '9-10', '7-8', '5-6'))

#############################################################################################################################
# Cluster configuration: (model x agegroup combination) per CPU
#############################################################################################################################

#create list of all kernel functions
# kernellist<-list(rbf, bayesianMeanTracker)
kernellist<-list(rbf)
#names of all kernel functions
kernelnames<-c("RBF", "BMT")

#create list of all acquisition functions
acqlist = list(ucb, greedyMean, epsilonGreedy)
#names of all acquisition functions
acqnames = c("UCB", "GM", "EG")

# all combinations of RBF kernel and acqusition functions will be needed
modelcombs = expand.grid(1, 1:length(acqlist))
# additionally, BMT-UCB is needed (lambda lesion)
modelcombs = rbind(modelcombs, expand.grid(2,1))

# list of agegroups
agegroups = levels(d$agegroup)

# learning curves should be simulated for all age groups
combs<-expand.grid(1:nrow(modelcombs), 1:length(agegroups))

#Cluster id from qsub
clusterid <- as.integer(commandArgs(TRUE)[1]) #Cluster id, corresponds to an integer used to indicate which combination of kernel and acquisition function to simulate
if(is.na(clusterid)){clusterid<-sample(1:nrow(combs), 1)} #sample a random number if not provided=

modelId <- combs[clusterid,1] #used to identify unique model
kernelId = modelcombs[modelId, 1] # identify kernel
acqId = modelcombs[modelId, 2] # identify acquisition function

agegroupId <- combs[clusterid,2] #used to identify a unique agegroup

#############################################################################################################################
# RANDOM MODEL
# to read the previous simulation from disk:
# randomDF <- read.csv("rationalModels/random.csv")
#############################################################################################################################

randomModel <- function(replications, outputfile){
  reward <- mcmapply(1:replications,
                     FUN=function(x) sample(environments[[sample(1:40,1)]][,'y'],
                                            trials, replace = TRUE)*50, # rewards between 0-50
                     mc.preschedule=TRUE)#,
                     #mc.cores=cores)
  #put into dataFrame
  randomDF <- data.frame(trial=seq(1:trials)-1,
                         meanReward = rowMeans(reward), #reward averaged over replications
                         meanSE = apply(reward, 1, FUN = function(x) sd(x)/sqrt(length(x))))
  #add model label
  randomDF$model <- "Random"
  #write output
  if(!is.null(outputfile)){
    write.csv(randomDF, outputfile) 
  }
  return(randomDF)
}

#randomDF <- randomModel(reps, "rationalModels/random.csv")

#############################################################################################################################
# BMT-UCB Model
# to read the previous simulation from disk:
# bmtDF <- read.csv("rationalModels/BMTUCB.csv")
#############################################################################################################################

bmtRationalModel <- function(replications, parameters, acq=ucb, cores){
  total <- list()
  
  #choice matrix
  choices <- expand.grid(0:7, 0:7) #build choice matrix
  names(choices)<-c("x1","x2")

  reward <- mcmapply(1:replications, FUN=function(x){
    #sample parameters
    params <- parameters[sample(1:nrow(parameters), 1),]
    kError <- params$kError
    beta <- params$beta
    tau <- params$tau
    #randomly choose environment
    envNum <- sample(1:40,1) 
    #1st trial is random
    location <- sample(1:64,1)#each location as an integer from 1:25; first value is random, so treat both environments as replications of the same actions
    #Observations
    X1 <- choices[location,'x1']
    X2 <- choices[location,'x2']
    #reward
    reward<- c()
    reward[1] <- Y <- (environments[[envNum]][location,"y"])*50 # rewards between 0-50
    #posterior
    prevPost <- NULL  #initialize previous posterior for sequential updating of BMT posterior
    for (j in 2:trials){ #after that, loop through remaining trials and make decisions based on BMT-UCB preditions
      #update posterior predictions
      post <- bayesianMeanTracker(x = as.matrix(cbind(X1,X2))[j-1,], y = (reward[j-1]-25)/50, prevPost = prevPost, theta = c(kError))
      prevPost <- post  #save new posterior as prevPost for next round
      #compute acquisition function evaluation
      utilityVec <- acq(post, pars = c(beta))
      #to prevent overflow, subtract max of q(x) vector 
      utilityVec <- utilityVec - max(utilityVec)
      #compute softmax choice probabilities
      p <- exp(utilityVec/tau)
      p <- p/sum(p)
      #Sample next choice
      location <- sample(1:64,1, prob=p, replace=TRUE)
      #update reward, X1, X2, and Y for both smooth and rough
      reward[j] <- (environments[[envNum]][location,"y"]) * 50
      X1 <- c(X1, choices[location, 'x1'])
      X2 <- c(X2, choices[location, 'x2'])
      Y <- c(Y,  reward[j])}
    reward},
    mc.preschedule = TRUE, mc.cores=cores)
  
  total[[length(total)+1]] <- reward

  total <- do.call("rbind", total)
  #put into dataFrame
  bmtDF <-data.frame(trial=seq(1:trials)-1,
                     meanReward = rowMeans(total), #reward averaged over replications
                     meanSE = apply(total, 1, FUN = function(x) sd(x)/sqrt(length(x))))
  
  #add model label and agegroup
  bmtDF$model <- modelName
  bmtDF$agegroup <- agegroup
  #write to csv
  if(!is.null(outputfile)){
    write.csv(bmtDF, outputfile)
  }
  return(bmtDF)
}

#############################################################################################################################
# GP Model
# to read the previous simulation from disk:
# gpDF <- read.csv("rationalModels/GPUCB.csv")
#############################################################################################################################

gpRationalModel <- function(replications, parameters, acq, kernel=rbf, cores){
  total <- list()
  #choice matrix
  choices <- expand.grid(0:7, 0:7) #build choice matrix
  names(choices)<-c("x1","x2")
  #run for smooth environments
  reward<- mcmapply(1:replications, FUN=function(x){
    #sample parameters
    params <- parameters[sample(1:nrow(parameters), 1),]
    lambda <- params$lambda
    beta <- params$beta
    tau <- params$tau
    epsilon <- params$epsilon
    #randomly choose environment
    envNum <- sample(1:40,1) 
    #1st trial is random
    location <- sample(1:64,1)#each location as an integer from 1:25; first value is random, so treat both environments as replications of the same actions
    #Observations
    X1 <- choices[location,'x1']
    X2 <- choices[location,'x2']
    #reward
    reward<- c()
    reward[1] <- Y <- (environments[[envNum]][location,"y"])*50
    for (j in 2:trials){ #after that, loop through remaining trials and make decisions based on GP preditions
      #compute posterior predictions
      post <- gpr(X.test = choices, theta = c(lambda, lambda, 1, 0.0001), X = cbind(X1,X2), Y = (Y-25)/50, k = kernel) #scale observed Y to zero mean, variance of 1
      
      #Slightly different function calls for each acquisition function
      if (inherits(acq, "UCB")){ #UCB takes a beta parameter
        utilityVec <- acq(post, pars = c(beta))
      }else if (inherits(acq, "epsilonGreedy")){
        p <- acq(post, beta, epsilon)
      }else{ #any other
        utilityVec <- acq(post)
      }
      if (inherits(acq, "softmax")){
        utilityVec <- utilityVec - max(utilityVec) #avoid overflow
        #compute softmax choice probabilities
        p <- exp(utilityVec/tau)
        p <- p/sum(p)
        #avoid underflow by setting a floor and a ceiling
        p <- (pmax(p, 0.00001))
        p <- (pmin(p, 0.99999))
      }
      
      #Sample next choice
      location <- sample(1:64,1, prob=p, replace=TRUE)
      #update reward, X1, X2, and Y
      reward[j] <- (environments[[envNum]][location,"y"]) * 50
      X1 <- c(X1, choices[location, 'x1'])
      X2 <- c(X2, choices[location, 'x2'])
      Y <- c(Y,  reward[j])}
    reward},
    mc.preschedule = TRUE, mc.cores=cores)
  
  total[[length(total)+1]] <- reward

  #put it all together
  total <- do.call("rbind", total)
  gpDF <-data.frame(trial=seq(1:trials)-1,
                    meanReward = rowMeans(total),
                    meanSE = apply(total, 1, FUN = function(x) sd(x)/sqrt(length(x))))
  
  return(gpDF)
}

#############################################################################################################################
# run models
#############################################################################################################################
# subset parameters to use to simulate the learning curves
subparams = subset(params, agegroup==agegroups[agegroupId] & kernel==kernelnames[kernelId] & acq==acqnames[acqId])

# simulate choices
if (kernelId==1) { #gp model
  # we use different acqusition functions for the gp kernel
  outDF <- gpRationalModel(reps, parameters = subparams, acq = acqlist[[acqId]], cores = cores)
} else { #bmt model
  # but we only use the ucb acquisition function for the bmt kernel
  outDF <- bmtRationalModel(reps, parameters = subparams, cores = cores)
}

# add model name and age group
kernelname = ifelse(kernelnames[kernelId]=='RBF', 'GP', 'BMT') # change RBF to GP
modelname = paste0(kernelname, "-", acqnames[acqId])

outDF$model <- modelname
outDF$agegroup <- agegroups[agegroupId]

# save output
# name of the output file
outputfile = paste0("rationalModels/", kernelname, acqnames[acqId], "_", agegroups[agegroupId], ".csv")
#write to csv
write.csv(outDF, outputfile)