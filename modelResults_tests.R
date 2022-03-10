# Analyze model results
# Anna Giron, Charley Wu, 2022

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
# IMPORT DATA 
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
  mutate(experiment=factor(experiment, levels=c('Meder (2021)', 'Schulz (2019)', 'Adolescent'))) %>%
  mutate(agegroup=factor(agegroup,
                         levels=c('[25,55]', '[18,25)', '[14,18)', '[11,14)', '[9,11)', '[7,9)', '[5,7)'),
                         labels=c('25-55', '18-24', '14-17', '11-13', '9-10', '7-8', '5-6')))

# append mean reward to model results
meanRew = ddply(behavior, ~id, plyr::summarize, meanReward=mean(z))
params = merge(params, meanRew, by='id')

# normalize reward
params$meanReward = params$meanReward / 50


# import protected exceedance probability
pxpAll = read.csv('modelResults/pxp/PXP.csv', header=F)
pxp1 = read.csv('modelResults/pxp/PXP5-6.csv', header=F)
pxp2 = read.csv('modelResults/pxp/PXP7-8.csv', header=F)
pxp3 = read.csv('modelResults/pxp/PXP9-10.csv', header=F)
pxp4 = read.csv('modelResults/pxp/PXP11-13.csv', header=F)
pxp5 = read.csv('modelResults/pxp/PXP14-17.csv', header=F)
pxp6 = read.csv('modelResults/pxp/PXP18-24.csv', header=F)
pxp7 = read.csv('modelResults/pxp/PXP25-55.csv', header=F)

# combine
pxp = rbind(pxpAll, pxp1, pxp2, pxp3, pxp4, pxp5, pxp6, pxp7)
colnames(pxp) = c('GP-UCB', 'BMT-UCB', 'GP-GM', 'GP-EG')
pxp$agegroup = rev(c('25-55', '18-24', '14-17', '11-13', '9-10', '7-8', '5-6', 'Overall'))
pxp$agegroup = factor(pxp$agegroup, levels=rev(c('25-55', '18-24', '14-17', '11-13', '9-10', '7-8', '5-6', 'Overall')))
pxp = gather(pxp, ModelName, pxp, `GP-UCB`:`GP-EG`)

# only GP-UCB, BMT-UCB, GP-GM, GP-EG
params = subset(params, ModelName %in% c('GP-UCB', 'BMT-UCB', 'GP-GM', 'GP-EG'))

#############################################################################################################################
# Inspect Parameters
#############################################################################################################################
# how many estimates from Schulz et al. (2019) are beyond exp(4) (upper limit in the other experiments)
params %>% 
  filter(experiment=='Schulz (2019)') %>%
  group_by(ModelName) %>%
  summarise(`lambda > exp(4)` = length(which(lambda > exp(4))),
            `kError > exp(4)` = length(which(kError > exp(4))),
            `beta > exp(4)` = length(which(beta > exp(4))),
            `tau > exp(4)` = length(which(tau > exp(4))))


#############################################################################################################################
# Summary of Modeling Results
#############################################################################################################################
# Number of participants best described
models <-  rbind(rep(0, length(levels(params$ModelName))), rep(0, length(levels(params$ModelName))),
                 rep(0, length(levels(params$ModelName))), rep(0, length(levels(params$ModelName))),
                 rep(0, length(levels(params$ModelName))), rep(0, length(levels(params$ModelName))),
                 rep(0, length(levels(params$ModelName))))
colnames(models) <- levels(params$ModelName)
rownames(models) <- levels(params$agegroup)
twoMax = c()
for (pid in unique(params$id)){
  subDF <- subset(params, id==pid)
  best <- subDF$ModelName[which(subDF$R2==max(subDF$R2))]
  if(length(best)>1) {twoMax = c(twoMax, pid)}
  models[subDF$agegroup,best] <- models[subDF$agegroup,best] + 1
}
# # add best described to modelFit 
# params$bestDescribed <- NA
# for (age in rownames(models)){
#   for (mod in colnames(models)){
#     params[params$ModelName == mod & params$agegroup==age,]$bestDescribed <- models[age,mod]
#   }
# }

# overall
View(params %>% 
  group_by(ModelName) %>% 
  summarise(R2 = round(mean(R2), digits=2),
            nLL = round(mean(nLL), digits=2),
            lambda = round(median(lambda, na.rm=TRUE), digits=2),
            kError = round(median(kError, na.rm=TRUE), digits=2),
            beta = round(median(beta, na.rm=TRUE), digits=2),
            tau = round(median(tau, na.rm=TRUE), digits=2),
            epsilon = round(median(epsilon, na.rm=TRUE), digits=2)))

colSums(models)

# separated by agegroup
View(params %>% 
  group_by(agegroup, ModelName) %>% 
  summarise(R2 = round(mean(R2), digits=2),
            nLL = round(mean(nLL), digits=2),
            lambda = round(median(lambda, na.rm=TRUE), digits=2),
            kError = round(median(kError, na.rm=TRUE), digits=2),
            beta = round(median(beta, na.rm=TRUE), digits=2),
            tau = round(median(tau, na.rm=TRUE), digits=2),
            epsilon = round(median(epsilon, na.rm=TRUE), digits=2)))
models


#############################################################################################################################
# Protected exceedance probability
#############################################################################################################################
# save nLL to compute protected exceedance probability
# pxp across all age groups
dnLL = params[, (names(params) %in% c('id', 'ModelName', 'nLL'))]
dnLL = spread(dnLL, ModelName, nLL)
dnLL = dnLL[, !(names(dnLL) %in% 'id')]
write.table(dnLL, file="modelResults/pxp/nLL.csv", sep=",", row.names = F, col.names = F)
# compute pxp with PXP.ipynb

# pxp for each age group individually
# 5-6
dnLL1 = subset(params, agegroup=='5-6')[, (names(params) %in% c('id', 'ModelName', 'nLL'))]
dnLL1 = spread(dnLL1, ModelName, nLL)
dnLL1 = dnLL1[, !(names(dnLL1) %in% 'id')]
write.table(dnLL1, file="modelResults/pxp/nLL5-6.csv", sep=",", row.names = F, col.names = F)

# 7-8
dnLL2 = subset(params, agegroup=='7-8')[, (names(params) %in% c('id', 'ModelName', 'nLL'))]
dnLL2 = spread(dnLL2, ModelName, nLL)
dnLL2 = dnLL2[, !(names(dnLL2) %in% 'id')]
write.table(dnLL2, file="modelResults/pxp/nLL7-8.csv", sep=",", row.names = F, col.names = F)

# 9-10
dnLL3 = subset(params, agegroup=='9-10')[, (names(params) %in% c('id', 'ModelName', 'nLL'))]
dnLL3 = spread(dnLL3, ModelName, nLL)
dnLL3 = dnLL3[, !(names(dnLL3) %in% 'id')]
write.table(dnLL3, file="modelResults/pxp/nLL9-10.csv", sep=",", row.names = F, col.names = F)

# 11-13
dnLL4 = subset(params, agegroup=='11-13')[, (names(params) %in% c('id', 'ModelName', 'nLL'))]
dnLL4 = spread(dnLL4, ModelName, nLL)
dnLL4 = dnLL4[, !(names(dnLL4) %in% 'id')]
write.table(dnLL4, file="modelResults/pxp/nLL11-13.csv", sep=",", row.names = F, col.names = F)

# 14-17
dnLL5 = subset(params, agegroup=='14-17')[, (names(params) %in% c('id', 'ModelName', 'nLL'))]
dnLL5 = spread(dnLL5, ModelName, nLL)
dnLL5 = dnLL5[, !(names(dnLL5) %in% 'id')]
write.table(dnLL5, file="modelResults/pxp/nLL14-17.csv", sep=",", row.names = F, col.names = F)

# 18-25
dnLL6 = subset(params, agegroup=='18-24')[, (names(params) %in% c('id', 'ModelName', 'nLL'))]
dnLL6 = spread(dnLL6, ModelName, nLL)
dnLL6 = dnLL6[, !(names(dnLL6) %in% 'id')]
write.table(dnLL6, file="modelResults/pxp/nLL18-24.csv", sep=",", row.names = F, col.names = F)

# 25-55
dnLL7 = subset(params, agegroup=='25-55')[, (names(params) %in% c('id', 'ModelName', 'nLL'))]
dnLL7 = spread(dnLL7, ModelName, nLL)
dnLL7 = dnLL7[, !(names(dnLL7) %in% 'id')]
write.table(dnLL7, file="modelResults/pxp/nLL25-55.csv", sep=",", row.names = F, col.names = F)


View(pxp %>% 
  group_by(ModelName, agegroup) %>%
  summarise(`pxp > .01` = length(which(pxp > .01))))

