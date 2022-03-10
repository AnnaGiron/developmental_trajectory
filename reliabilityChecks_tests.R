# Reliability checks: Compare behavioral data and model results from different experiments
# Anna Giron, 2021

# house keeping
rm(list=ls())

# load packages
packages <- c('plyr', 'jsonlite', 'gridExtra', 'reshape2', 'stargazer', 'coefplot', 'cowplot',
              "grid", 'corrplot', 'ggbeeswarm', 'tidyverse', 'viridis', 'colorspace', 'ggrepel')
lapply(packages, require, character.only = TRUE)

source("dataProcessing.R")
source('statisticalTests.R')

theme_set(theme_bw(base_size=16)) # use the b&w theme

#############################################################################################################################
# Import data
#############################################################################################################################
# import fitted parameters
# params = paramsImport()

# import previously saved data frame
params = read.csv('data/modelFit.csv')

params = params %>%
  mutate(kernel=factor(kernel, levels=c('RBF', 'BMT'), labels=c('GP', 'BMT'))) %>%
  mutate(agegroup=factor(agegroup,
                         levels=c('[25,55]', '[18,25)', '[14,18)', '[11,14)', '[9,11)', '[7,9)', '[5,7)'),
                         labels=c('25-55', '18-24', '14-17', '11-13', '9-10', '7-8', '5-6'))) %>%
  mutate(experiment=factor(experiment, levels=c('Meder (2021)', 'Schulz (2019)', 'Adolescents')))

params$ModelName = paste(params$kernel, params$acq, sep="-")
params$ModelName = factor(params$ModelName)


# import behavioral data
# behavior = dataImport()

# import previously saved data frame
behavior = read.csv('data/behavioralData.csv')

# only smooth condition
behavior = subset(behavior, condition=='Smooth')

behavior = behavior %>%
  mutate(type_choice=factor(type_choice, levels=c('Repeat', 'Near', 'Far'))) %>%
  mutate(experiment=factor(experiment,
                           levels=c('Meder (2021)', 'Schulz (2019)', 'Adolescent'),
                           labels=c('Meder (2021)', 'Schulz (2019)', 'Adolescents'))) %>%
  mutate(agegroup=factor(agegroup,
                         levels=c('[25,55]', '[18,25)', '[14,18)', '[11,14)', '[9,11)', '[7,9)', '[5,7)'),
                         labels=c('25-55', '18-24', '14-17', '11-13', '9-10', '7-8', '5-6')))

# append mean reward to model results
meanRew = ddply(behavior, ~id, plyr::summarize, meanReward=mean(z))
params = merge(params, meanRew, by='id')

params = subset(params, acq=='UCB')


#############################################################################################################################
# Behavior
#############################################################################################################################
# compare performance of participants from different experiments for overlapping age ranges
dReward = ddply(behavior, ~id+agegroup+age_years+experiment, plyr::summarize, meanReward=mean(z))
# Schulz (2019) and Meder (2021) data
# age 7-9
ttestPretty(subset(dReward, experiment=='Schulz (2019)' & age_years>=7 & age_years <= 9)$meanReward,
            subset(dReward, experiment=='Meder (2021)' & age_years>=7)$meanReward, paired=FALSE)
# adolescent data set and Schulz (2019) data
# age > 20
ttestPretty(subset(dReward, experiment=='Adolescents' & age_years>20)$meanReward, 
            subset(dReward, experiment=='Schulz (2019)' & age_years>20)$meanReward, paired=FALSE)


#############################################################################################################################
# GP-UCB Model Results
#############################################################################################################################
# compare r2 values of data from different experiments for overlapping age ranges
# Schulz (2019) and Meder (2021) data
# age 7-9
ttestPretty(subset(params, experiment=='Schulz (2019)' & kernel=='GP' & acq=='UCB' & age_years>=7 & age_years <= 9)$R2,
            subset(params, experiment=='Meder (2021)' & kernel=='GP' & acq=='UCB' & age_years>=7)$R2, paired=FALSE)
# adolescent data set and Schulz (2019) data
# age > 20
ttestPretty(subset(params, experiment=='Adolescents' & kernel=='GP' & acq=='UCB' & age_years>20)$R2, 
            subset(params, experiment=='Schulz (2019)' & kernel=='GP' & acq=='UCB' & age_years>20)$R2, paired=FALSE)


# compare parameter estimates for data from different experiments for overlapping age ranges
# lambda
# Schulz (2019) and Meder (2021) data
# age 7-9
shapiro.test(subset(params, experiment=='Schulz (2019)' & kernel=='GP' & acq=='UCB' & age_years>=7 & age_years <= 9)$lambda)
shapiro.test(subset(params, experiment=='Meder (2021)' & kernel=='GP' & acq=='UCB' & age_years>=7)$lambda)
# p<.05 -> not normally distributed
ranktestPretty(subset(params, experiment=='Schulz (2019)' & kernel=='GP' & acq=='UCB' & age_years>=7 & age_years <= 9)$lambda,
               subset(params, experiment=='Meder (2021)' & kernel=='GP' & acq=='UCB' & age_years>=7)$lambda, paired=FALSE)


# adolescent data set and Schulz (2019) data
# age > 20
shapiro.test(subset(params, experiment=='Adolescents' & kernel=='GP' & acq=='UCB' & age_years>20)$lambda)
shapiro.test(subset(params, experiment=='Schulz (2019)' & kernel=='GP' & acq=='UCB' & age_years>20)$lambda)
# p<.05 -> not normally distributed
ranktestPretty(subset(params, experiment=='Adolescents' & kernel=='GP' & acq=='UCB' & age_years>20)$lambda, 
               subset(params, experiment=='Schulz (2019)' & kernel=='GP' & acq=='UCB' & age_years>20)$lambda, paired=FALSE)


# beta
# Schulz (2019) and Meder (2021) data
# age 7-9
shapiro.test(subset(params, experiment=='Schulz (2019)' & kernel=='GP' & acq=='UCB' & age_years>=7 & age_years <= 9)$beta)
shapiro.test(subset(params, experiment=='Meder (2021)' & kernel=='GP' & acq=='UCB' & age_years>=7)$beta)
# p<.05 -> not normally distributed
ranktestPretty(subset(params, experiment=='Schulz (2019)' & kernel=='GP' & acq=='UCB' & age_years>=7 & age_years <= 9)$beta,
               subset(params, experiment=='Meder (2021)' & kernel=='GP' & acq=='UCB' & age_years>=7)$beta, paired=FALSE)

# adolescent data set and Schulz (2019) data
# age > 20
shapiro.test(subset(params, experiment=='Adolescents' & kernel=='GP' & acq=='UCB' & age_years>20)$beta)
shapiro.test(subset(params, experiment=='Schulz (2019)' & kernel=='GP' & acq=='UCB' & age_years>20)$beta)
# p<.05 -> not normally distributed
ranktestPretty(subset(params, experiment=='Adolescents' & kernel=='GP' & acq=='UCB' & age_years>20)$beta, 
               subset(params, experiment=='Schulz (2019)' & kernel=='GP' & acq=='UCB' & age_years>20)$beta, paired=FALSE)


# tau
# Schulz (2019) and Meder (2021) data
# age 7-9
shapiro.test(subset(params, experiment=='Schulz (2019)' & kernel=='GP' & acq=='UCB' & age_years>=7 & age_years <= 9)$tau)
shapiro.test(subset(params, experiment=='Meder (2021)' & kernel=='GP' & acq=='UCB' & age_years>=7)$tau)
# p<.05 -> not normally distributed
ranktestPretty(subset(params, experiment=='Schulz (2019)' & kernel=='GP' & acq=='UCB' & age_years>=7 & age_years <= 9)$tau,
               subset(params, experiment=='Meder (2021)' & kernel=='GP' & acq=='UCB' & age_years>=7)$tau, paired=FALSE)

# adolescent data set and Schulz (2019) data
# age > 20
shapiro.test(subset(params, experiment=='Adolescents' & kernel=='GP' & acq=='UCB' & age_years>20)$tau)
shapiro.test(subset(params, experiment=='Schulz (2019)' & kernel=='GP' & acq=='UCB' & age_years>20)$tau)
# p<.05 -> not normally distributed
ranktestPretty(subset(params, experiment=='Adolescents' & kernel=='GP' & acq=='UCB' & age_years>20)$tau, 
               subset(params, experiment=='Schulz (2019)' & kernel=='GP' & acq=='UCB' & age_years>20)$tau, paired=FALSE)
