# Run optimization algorithms in the parameter space
# Anna Giron, 2021

# house keeping
rm(list=ls())

packages <- c('plyr', 'dplyr')
#invisible(lapply(packages, install.packages, character.only = TRUE))
invisible(lapply(packages, require, character.only = TRUE)) #loads packages

source('dataProcessing.R')

################################################################################################
# Data import
################################################################################################
# get parameter estimates
modelFit<-read.csv('data/paramsYoungestAgegroup.csv')

# import mean rewards
sim = read.csv('modelResults/simulatedModels/meanRewards.csv')


# normalize reward
sim$meanReward = sim$meanReward / 50

sim = sim[,!names(sim) %in% 'X']


################################################################################################
# Hill climbing algorithm
################################################################################################
hillClimbing = function(init, tempFunction, imax, shc) {

  theta = data.frame(lambda=init$lambda, beta=init$beta, tau=init$tau, reward=init$meanReward, i=0)
  theta$rewardScaled = NA
  theta$maxRewardDiff = NA
  theta$prob = NA
  
  i = 1
  
  # continue until max change over the last 10 iterations is below a specific threshold
  # while((dim(theta)[1] < 10 | (max(tail(theta$reward, n=10)) - min(tail(theta$reward, n=10))) > .005) &
  #       i <= imax) {
  while(i < imax) {
    # determine possible parameter changes
    prevLambda = tail(theta$lambda, n=1)
    prevLambdaPos = which(lambdas==prevLambda)
    lambdaHigher = lambdas[prevLambdaPos+1]
    lambdaLower = lambdas[prevLambdaPos-1]
    
    prevBeta = tail(theta$beta, n=1)
    prevBetaPos = which(betas==prevBeta)
    betaHigher = betas[prevBetaPos+1]
    betaLower = betas[prevBetaPos-1]
    
    prevTau = tail(theta$tau, n=1)
    prevTauPos = which(taus==prevTau)
    tauHigher = taus[prevTauPos+1]
    tauLower = taus[prevTauPos-1]
    
    prevReward = tail(theta$reward, n=1)
    
    # all possible steps
    steps = expand.grid(lambda=na.omit(c(lambdaLower, prevLambda, lambdaHigher)),
                        beta=na.omit(c(betaLower, prevBeta, betaHigher)),
                        tau=na.omit(c(tauLower, prevTau, tauHigher)))
    
    # expected reward for each step
    rewards = lapply(seq(1, dim(steps)[1]), FUN=function(x) {
      subset(sim, lambda==steps$lambda[x] & beta==steps$beta[x] & tau==steps$tau[x])$meanReward
    })
    steps$reward = unlist(rewards)
    
    if(shc) { # stochastic hill climbing
      
      # temperature
      temp = eval(call(tempFunction, i, imax))
      
      steps$maxRewardDiff = max(steps$reward) - min(steps$reward)
      
      # probability of being selected for each option (softmax)
      steps$rewardScaled = steps$reward - max(steps$reward)
      steps$prob = exp(steps$rewardScaled / temp)
      
      # in last round of linear cooling: temperature is 0
      # change prob of best option from NA to 1
      if(tempFunction=='linearCooling' & i==imax-1) {
        steps$prob = ifelse(is.na(steps$prob), 1, steps$prob)
      }
      steps$prob = steps$prob / sum(steps$prob)
      
      # select next step
      step = steps[sample(1:dim(steps)[1], size=1, prob=steps$prob),]
      
      #step = step[, !names(step) %in% c('prob')]
      step$i = i
      
      theta = rbind(theta, step)
      
      
    } else { # simulated annealing
      
      # randomly pick next step
      step = steps[sample(nrow(steps), 1), ]
      step$i = i
      step$rewardScaled = NA
      step$maxRewardDiff = NA
      step$prob = NA
      
      # temperature
      temp = eval(call(tempFunction, i, imax))
      
      newReward = step$reward
      prevReward = tail(theta$reward, n=1)
      
      if(newReward >= prevReward) { # always accept if new option is better than previous option
        theta = rbind(theta, step)
      } else if(exp((newReward - prevReward) / temp) >= runif(1, min=0, max=1)) { # accept worse option with some probability
        theta = rbind(theta, step)
      } else { # keep previous option
        step = theta[dim(theta)[1], ]
        step$i = i
        theta = rbind(theta, step)
      }
    }
    
    i = i+1
  }
  
  return(theta)
}


################################################################################################
# Temperature Functions
################################################################################################

fastCooling = function(iteration, imax) {
  return(1 / (1+iteration))
}

expCooling = function(iteration, imax) {
  return(exp(-iteration^(1/3)))
}

linearCooling = function(iteration, imax) {
  return(1 - (iteration+1)/imax)
}


################################################################################################
# Run the algorithm
################################################################################################
coolingFuncs = c('fastCooling', 'expCooling', 'linearCooling')

methods = c('SHC', 'SA')

ids = unique(modelFit$id)

trajectories = expand.grid(temp=coolingFuncs, method=methods, id=ids, leaveoutindex=c(2:5))
trajectories$temp = as.character(trajectories$temp)

clusterid = as.integer(commandArgs(TRUE)[1]) # Cluster id, corresponds to an integer used to indicate which combination of kernel and acquisition function to simulate
# clusterid = sample(1:nrow(trajectories),1) #sample random cluster id for testing

# initialization
# use median parameter estimates of youngest age group
i = trajectories[clusterid, "id"]
idx = trajectories[clusterid, "leaveoutindex"]

lambdas = sort(unique(sim$lambda))
closestLambda = which.min(abs(lambdas - subset(modelFit, id==i & leaveoutindex==idx)$lambda))

betas = sort(unique(sim$beta))
closestBeta = which.min(abs(betas - subset(modelFit, id==i & leaveoutindex==idx)$beta))

taus = sort(unique(sim$tau))
closestTau = which.min(abs(taus - subset(modelFit, id==i & leaveoutindex==idx)$tau))

init = subset(sim, lambda==lambdas[closestLambda] & beta==betas[closestBeta] & tau==taus[closestTau])

# SHC or SA?
shc = ifelse(trajectories[clusterid,]$method == 'SHC', TRUE, FALSE)

# SA needs more trajectories than SHC
nIter=1500

# run simulated annealing algorithm
t = hillClimbing(init=init, tempFunction=trajectories[clusterid,]$temp, imax=nIter, shc=shc)
t$coolingFunc = trajectories[clusterid,]$temp
t$method = trajectories[clusterid,]$method
t$id = i
t$leaveoutindex = idx

# save trajectory
filename = paste0('hillClimbingAlgorithm/batch2/trajecotry_', clusterid, '.csv')
write.csv(t, filename)

