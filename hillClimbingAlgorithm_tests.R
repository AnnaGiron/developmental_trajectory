# Compare hill climbing trajectories
# Anna Giron, Charley Wu, 2022

rm(list=ls())

packages <- c('plyr', 'dplyr', 'tidyr', 'ggplot2', 'cowplot', 'viridis', 'colorspace', 'devtools', 'ggnewscale', 'rmisc')
invisible(lapply(packages, require, character.only = TRUE)) #loads packages

source('dataProcessing.R')
source('statisticalTests.R')

modelPal <- c('black', '#6c7195', '#ea9d67', '#7ec3aa')
paramPal = c("#FFEE67", '#27AD88', "#D1495B")
################################################################################################
# Data import
################################################################################################
path = 'hillClimbingAlgorithm/batch2/'
filenames = list.files(path=path, pattern='*.csv')
filenames = paste0(path, filenames)
trajectories = ldply(filenames, read.csv)
trajectories$Algorithm = paste(trajectories$method, trajectories$coolingFunc, sep='-')
trajectories$trajectory = rep(1:720, each=1500)


# model results
params = read.csv('data/modelFit.csv')

params = params %>%
  mutate(kernel=factor(kernel, levels=c('RBF', 'BMT'), labels=c('GP', 'BMT'))) %>%
  mutate(agegroup=factor(agegroup,
                         levels=rev(c('[25,55]', '[18,25)', '[14,18)', '[11,14)', '[9,11)', '[7,9)', '[5,7)')),
                         labels=rev(c('25-55', '18-24', '14-17', '11-13', '9-10', '7-8', '5-6')))) %>%
  mutate(experiment=factor(experiment, levels=c('Meder (2021)', 'Schulz (2019)', 'Adolescent')))

params$ModelName = paste(params$kernel, params$acq, sep="-")
params$ModelName = factor(params$ModelName)

# only G-UCB model
params = subset(params, ModelName=='GP-UCB')

# import previously saved data frame
behavior = read.csv('data/behavioralData.csv')

# only smooth condition
behavior = subset(behavior, condition=='Smooth')

behavior = behavior %>%
  mutate(type_choice=factor(type_choice, levels=c('Repeat', 'Near', 'Far'))) %>%
  mutate(experiment=factor(experiment, levels=c('Meder (2021)', 'Schulz (2019)', 'Adolescent'))) %>%
  mutate(agegroup=factor(agegroup,
                         levels=c('[25,55]', '[18,25)', '[14,18)', '[11,14)', '[9,11)', '[7,9)', '[5,7)'),
                         labels=c('25-55', '18-24', '14-17', '11-13', '9-10', '7-8', '5-6')))

# append mean reward to model results
meanRew = ddply(behavior, ~id, plyr::summarize, meanReward=mean(z))
params$meanReward = params$meanReward / 50


# simulation results
sim = read.csv('modelResults/simulatedModels/meanRewards.csv')

# normalize reward
sim$meanReward = sim$meanReward / 50
sim = sim[,!names(sim) %in% 'X']


################################################################################################
# Compare human and algorithm trajectories
################################################################################################

#Test if the SGD performed better than humans
humanAdults <-  behavior %>% filter(agegroup=='25-55') %>% group_by(id) %>% summarize(reward = mean(z)/50) 
bestAlgo <- trajectories %>% filter(method ==  'SGD' & i == 1499 & coolingFunc=='fastCooling') %>% group_by(id) %>% summarize(reward = mean(reward))

ttestPretty(bestAlgo$reward, humanAdults$reward) #average across simulated participants
ttestPretty(subset(trajectories,method ==  'SGD' & i == 1499 & coolingFunc=='fastCooling')$reward, humanAdults$reward) #each simulated trajectory
ttestPretty(subset(trajectories,method ==  'SA' & i == 1499 & coolingFunc=='fastCooling')$reward, humanAdults$reward) #each simulated trajectory

ttestPretty(subset(trajectories,method ==  'SGD' & i == 1499 & coolingFunc=='fastCooling')$reward, subset(trajectories,method ==  'SA' & i == 1499 & coolingFunc=='fastCooling')$reward)


adultPerformance <- behavior %>% filter(agegroup=='25-55') %>% group_by(id) %>% summarize(reward = mean(z)/50) %>% pull(reward)
adultRange <- t.test(adultPerformance, conf.level = 0.95)$conf.int
