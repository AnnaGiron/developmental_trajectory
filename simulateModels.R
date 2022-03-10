#Simulation of models
# Anna Giron, 2021

#housekeeping
rm(list=ls())

#source of modeling code
source("models.R")

#packages
packages <- c('plyr', 'dplyr', 'jsonlite', 'lsr', 'BayesFactor', 'matrixcalc')
#invisible(lapply(packages, install.packages, character.only = TRUE))
invisible(lapply(packages, require, character.only = TRUE)) #loads packages

#extract environments
environments <- lapply(fromJSON("data/smoothKernel.json"), FUN=function(x) matrix(as.numeric(unlist(x)), ncol=3, byrow=TRUE, dimnames=list(seq(1,64),  c('x2', 'y', 'x1'))))

env <- as.data.frame(environments[[1]])
for (i in 2:40){
  env<-rbind(env,as.data.frame(environments[[i]]))
}
env$en<-rep(1:40, each=64)

#Scaling reward range
env$y<-env$y*50

# simulation rounds for each parameter combination
replications = 100
# replications = 2


#############################################################################################################################
# Import fitted parameters
#############################################################################################################################
# get parameter estimates
modelFit = read.csv('data/modelFit.csv')

modelFit = modelFit %>%
  mutate(kernel=factor(kernel, levels=c('RBF', 'BMT'), labels=c('GP', 'BMT'))) %>%
  mutate(agegroup=factor(agegroup,
                         levels=rev(c('[25,55]', '[18,25)', '[14,18)', '[11,14)', '[9,11)', '[7,9)', '[5,7)')),
                         labels=rev(c('25-55', '18-24', '14-17', '11-13', '9-10', '7-8', '5-6')))) %>%
  mutate(experiment=factor(experiment, levels=c('Meder (2021)', 'Schulz (2019)', 'Adolescents')))

# only use GP-UCB model
modelFit = subset(modelFit, kernel=='GP' & acq=='UCB')


#############################################################################################################################
# Simulating performance for different parameter values
#############################################################################################################################
# Tukey's fence to compute upper and lower bound for each parameter
tukeysFence <- function(x, k = 1.5, na.rm = TRUE) {
  quar <- quantile(x, probs = c(0.25, 0.75), na.rm = na.rm)
  iqr <- diff(quar)
  
  lower = quar[1] - k * iqr
  upper = quar[2] + k * iqr
  
  return(c(lower, upper))
}

# parameter range based on Tukey's fence
lambda = tukeysFence(log(modelFit$lambda))
beta = tukeysFence(log(modelFit$beta))
tau = tukeysFence(log(modelFit$tau))

# lower bound for tau below range defined for model fitting - set to lower bound
tau[1] = -5

# parameters to simulate
params = expand.grid(lambda=exp(seq(lambda[1], lambda[2], len=100)), beta=exp(seq(beta[1], beta[2], len=100)),
                     tau=exp(seq(tau[1], tau[2], len=100)))

#############################################################################################################################
# Simulation
#############################################################################################################################
# Cluster id from qsub
# per job, run 1000 simulations with different parameter combinations
# ids 1 - 1000
# clusterid <- sample(1:1000, 1) # sample random cluster id for testing
clusterid = as.integer(commandArgs(TRUE)[1]) # Cluster id, corresponds to an integer used to indicate which combination of kernel and acquisition function to simulate

params$clusterid = rep(1:1000, each=1000)

# parameters for simulation based on cluster id
lambdas = params[params$clusterid==clusterid,1]
betas = params[params$clusterid==clusterid,2]
taus = params[params$clusterid==clusterid,3]

dparams = data.frame(lambda=numeric(), beta=numeric(), tau=numeric(), mu=numeric(), replication=numeric())

Xstar<-cbind(env$x1[1:64], env$x2[1:64])
k<-rbf

for (i in 1:length(lambdas)) {
# for (i in 1:1) {
  lambda = lambdas[i]
  beta = betas[i]
  tau = taus[i]
  parVec <- c(lambda, lambda, 1, .0001) 
  
  mu = numeric()
  for (round in 1:replications){
    enselect<-sample(1:40, 1)
    environ<-subset(env, en==enselect)
    ind<-sample(1:64,1)
    #X matrix
    X<-cbind(environ$x1[ind], environ$x2[ind])
    #y matrix
    y<-as.matrix(environ$y[ind])
    #loop through trials
    for (trial in 1:25){
      #output by GP with particular parameter settings
      #don't forget mean centering and standardization
      out<-gpr(X.test=Xstar, theta=parVec, X=X, Y=(y-25)/50, k=k)
      #utility vector by UCB
      utilityVec<-ucb(out, beta)
      #avoid overflow
      utilities <- utilityVec - max(utilityVec)
      #softmaximization
      p <- exp(utilities/tau)
      #probabilities
      p <- p/colSums(p)
      #numerical overflow
      p <- (pmax(p, 0.00001))
      p <- (pmin(p, 0.99999))
      #index is sampled proportionally to softmaxed utility vector
      ind<-sample(1:64, 1, prob=p)
      #bind X-observations
      X<-rbind(X, cbind(environ$x1[ind], environ$x2[ind]))
      #bind y-observations
      y<-rbind(y, as.matrix(environ$y[ind]))
    }
    # remove first randomly revealed reward value
    y = y[-1]
    
    cur = data.frame(lambda=lambda, beta=beta, tau=tau, mu=mean(y), replication=round)
    dparams = rbind(dparams, cur)
  }
}

filename = paste0("modelResults/simulatedModels/batch1/simulatedModels_", clusterid, ".csv")
write.csv(dparams, filename)
