# Plot behavioral results
# Anna Giron, 2021

# house keeping
rm(list=ls())

packages <- c('plyr', 'dplyr', 'tidyr', 'ggplot2', 'cowplot', 'viridis', 'entropy', 'lme4', 'sjPlot', 'brms', 'withr', 'tidyr', 'ggbeeswarm')
#invisible(lapply(packages, install.packages, character.only = TRUE))
invisible(lapply(packages, require, character.only = TRUE)) #loads packages

source('dataProcessing.R')
source('statisticalTests.R')

################################################################################################
# Data import
################################################################################################
# import data
# d = dataImport()

# import previously saved data frame
d = read.csv('data/behavioralData.csv')

d = d %>%
  mutate(type_choice=factor(type_choice, levels=c('Repeat', 'Near', 'Far'))) %>%
  mutate(experiment=factor(experiment, levels=c('Meder (2021)', 'Schulz (2019)', 'Adolescent'))) %>%
  mutate(agegroup=factor(agegroup,
                         levels=c('[25,55]', '[18,25)', '[14,18)', '[11,14)', '[9,11)', '[7,9)', '[5,7)'),
                         labels=c('25-55', '18-24', '14-17', '11-13', '9-10', '7-8', '5-6')))

# normalize reward and previous reward
d$z = d$z / 50
d$previous_reward = d$previous_reward / 50

# Wrapper function for brm models such that it saves the full model the first time it is run, otherwise it loads it from disk
modelpath = 'brmsModels/'

run_model <- function(expr, modelName, path=modelpath, reuse = TRUE) {
  path <- paste0(path,'/', modelName, ".brm")
  if (reuse) {
    fit <- suppressWarnings(try(readRDS(path), silent = TRUE))
  }
  if (is(fit, "try-error")) {
    fit <- eval(expr)
    saveRDS(fit, file = path)
  }
  fit
}


# extract environments
environments <- lapply(fromJSON("data/smoothKernel.json"), FUN=function(x) matrix(as.numeric(unlist(x)), ncol=3, byrow=TRUE, dimnames=list(seq(1,64),  c('x2', 'y', 'x1'))))

env <- as.data.frame(environments[[1]])
for (i in 2:40){
  env<-rbind(env,as.data.frame(environments[[i]]))
}
env$en<-rep(1:40, each=64)


################################################################################################
# Demographic information
################################################################################################
# participants per experiment
dAge = ddply(d, ~id+agegroup+experiment, plyr::summarize, age = mean(age_years))
pAge = ggplot(dAge, aes(x=age, fill=experiment)) +
  geom_histogram(alpha=0.7, color='black', binwidth=1) +
  facet_grid(experiment~.) +
  xlab("Age") +
  ylab("Count") +
  scale_color_viridis(discrete=TRUE, direction=-1) +
  scale_fill_viridis(discrete=TRUE, direction=-1) +
  theme_classic() +
  theme(text = element_text(size=15,  family="sans"), legend.position='none',
        strip.background=element_blank(), panel.grid.minor = element_line(),
        strip.text.y = element_blank())

pAge
# ggsave(filename = "plots/histogramAge.pdf", plot = pAge, height=4, width=5, units = "in")


################################################################################################
# Performance
################################################################################################
# Average rewards
# exclude first randomly revealed tile
dReward = ddply(subset(d, trial!=0), ~id+agegroup, plyr::summarize, meanReward=mean(z))
dReward$agegroup = factor(dReward$agegroup,
                          levels=rev(c('25-55', '18-24', '14-17', '11-13', '9-10', '7-8', '5-6')),
                          labels=rev(c('25-55', '18-24', '14-17', '11-13', '9-10', '7-8', '5-6')))
pReward = ggplot(dReward, aes(x=agegroup, y=meanReward, fill=agegroup, color=agegroup)) +
  geom_boxplot(color='black', outlier.shape = NA, width=.5) +
  geom_hline(yintercept=.5, linetype='dashed', color='red') + # random choice model
  geom_quasirandom(position=position_jitter(0.2), alpha=0.3, color='black', size=.5) +
  stat_summary(fun=mean, geom="point", shape=23, fill="white", color='black', size=1.) + 
  scale_y_continuous(limits = c(.4, 1.), breaks = seq(.4, 1., by = .1)) +
  scale_x_discrete(guide = guide_axis(n.dodge=2))+
  xlab('Age Group') +
  ylab('Normalized Mean Reward') +
  scale_color_viridis(discrete=TRUE, direction=-1) +
  scale_fill_viridis(discrete=TRUE, direction=-1) +
  ggtitle('Performance') +
  theme_classic() +
  theme(legend.position='none',
        strip.background=element_blank())

pReward


# Learning curves - mean reward
dLearningCurves <- ddply(d, ~id+trial+agegroup, plyr::summarize, meanReward = mean(z))
pLearningCurves <- ggplot(dLearningCurves, aes(x=trial, y=meanReward, color=agegroup, fill=agegroup, group=agegroup))+
  stat_summary(fun.y = mean, geom='line', size=.5)+
  stat_summary(fun.data = mean_cl_boot, geom='ribbon', alpha=0.3, color=NA)+
  geom_hline(yintercept=0.5, linetype='dashed', color='red') + # random choice model
  coord_cartesian(ylim = c(.45, 1.)) +
  xlab('Trial')+
  ylab('Normalized Reward ± 95% CI')+
  labs(color='Age Group') +
  labs(fill='Age Group') +
  scale_color_viridis(discrete=TRUE, direction=1) +
  scale_fill_viridis(discrete=TRUE, direction=1) +
  ggtitle('Mean Reward') +
  theme_classic() +
  theme(strip.background=element_blank(), legend.position='none')

pLearningCurves


# Learning curves - max reward
dLearningCurvesMax = data.frame()
for (i in unique(d$id)) {
  for (r in unique(subset(d, id==i)$round)) {
    # compute maximum reward found so far
    maxRews = lapply(seq(0, max(subset(d, id==i)$trial)), FUN=function(t) {
      max(subset(d, id==i & round==r & trial<=t)$z)
    })
  
    trials = length(maxRews)
    subd = subset(d, id==i & round==r)
    
    dat = data.frame(id=rep(i, trials), round=rep(r, trials), trial=seq(0, trials-1), age_months=subd$age_months,
                     age_years=subd$age_years, agegroup=subd$agegroup,
                     maxReward=unlist(maxRews))
    dLearningCurvesMax = rbind(dLearningCurvesMax, dat)
  }
}

dLearningCurvesMax = ddply(dLearningCurvesMax, ~id+trial+agegroup, plyr::summarize, meanMaxRew = mean(maxReward))


# random baseline
dRandom = data.frame(iteration=numeric(), env=numeric(), reward=numeric(), trial=numeric())

for (i in 1:100){
  env_num = sample(1:40, 1) # sample an environment
  
  maxRewards = c()
  
  for (t in 1:26){
    r = sample(subset(env, en==env_num)$y, 1) # sample choice
      
    currentMax = as.double(max(maxRewards, r))
    maxRewards = c(maxRewards, currentMax)
  }
  
  cur = data.frame(iteration=rep(i, 26),
                   env=rep(env_num, 26), 
                   reward=maxRewards, 
                   trial=seq(0,25))
  
  dRandom = rbind(dRandom, cur)
}

dRandom = ddply(dRandom, ~trial, plyr::summarize, meanReward=mean(reward))


pLearningCurvesMax = ggplot()+
  scale_y_continuous(limits = c(0.0, 1.), breaks = seq(0.0, 1., by = .1)) +
  coord_cartesian(ylim=c(.4, 1.)) +
  stat_summary(dLearningCurvesMax, mapping=aes(x=trial, y=meanMaxRew, color=agegroup, group=agegroup), fun = mean, geom='line', size=.5)+
  stat_summary(dLearningCurvesMax, mapping=aes(x=trial, y=meanMaxRew, color=agegroup, fill=agegroup, group=agegroup), fun.data = mean_cl_boot, geom='ribbon', alpha=0.3, color=NA)+
  geom_line(dRandom, mapping=aes(x=trial, y=meanReward), linetype='dashed', color='red', size=.5) + # random choice model
  xlab('Trial') +
  ylab('Normalized Reward ± 95% CI') +
  labs(color='Age Group') +
  labs(fill='Age Group') +
  coord_cartesian(ylim=c(0.45, 1))+
  scale_color_viridis(discrete=TRUE, direction=1) +
  scale_fill_viridis(discrete=TRUE, direction=1) +
  ggtitle('Max Reward') +
  theme_classic() +
  theme(strip.background=element_blank(), legend.position='none')

pLearningCurvesMa


################################################################################################
# Search Decisions
################################################################################################
# Regression model: modulation of distance as a function of reward
DistancePrevReward = run_model(brm(distance ~ previous_reward * agegroup + (previous_reward + agegroup | id),
                                   data = d, cores = 4, iter = 4000, warmup = 1000,
                                   control = list(adapt_delta = 0.99)), modelName = 'DistancePrevReward')

# generate predictions
prevReward = seq(0,50) / 50 # normalized reward
agegroup = levels(d$agegroup)
newdat = expand.grid(previous_reward=prevReward, agegroup=agegroup)
# predict distance based on previous reward
preds = fitted(DistancePrevReward, re_formula=NA, newdata=newdat, probs=c(.025, .975))
predsDF = data.frame(previous_reward=rep(prevReward, 7),
                     agegroup=rep(levels(d$agegroup), each=length(prevReward)),
                     distance=preds[,1],
                     lower=preds[,3],
                     upper=preds[,4])

# average distance
grid = expand.grid(x1=0:7, x2=0:7, y1=0:7, y2=0:7)
grid$distance = NA

for(i in 1:dim(grid)[1]){
  grid$distance[i] <- dist(rbind(c(grid$x1[i], grid$x2[i]), c(grid$y1[i], grid$y2[i])), method = "manhattan")
}

meanDist = mean(grid$distance)

# plot predictions
pDistRewardRegression <- ggplot()+
  stat_summary(d, mapping=aes(x=previous_reward, y=distance, color=agegroup, fill=agegroup), fun.y=mean, geom='point', alpha=0.7, size=.5)+
  geom_line(predsDF, mapping=aes(x=previous_reward, y=distance, color=agegroup), size=.3) +
  geom_ribbon(predsDF, mapping=aes(x=previous_reward, y=distance, ymin=lower, ymax=upper, fill=agegroup), alpha=.3) +
  geom_hline(yintercept=meanDist, linetype='dashed', color='red') + # mean distance
  xlab('Normalized Previous Reward')+
  ylab('Distance to Next Option')+
  scale_color_viridis(discrete=TRUE, name='Age Group', direction=1) +
  scale_fill_viridis(discrete=TRUE, name='Age Group', direction=1) +
  ggtitle('Search Distance ~ Reward') +
  theme_classic() +
  theme(strip.background=element_blank(), legend.position='none')

pDistRewardRegression


# # plot model coefficients
# set_theme(geom.label.color = 'black', base = theme_classic())
# 
# pCoeffs = sjPlot::plot_model(DistancePrevReward) +
#   font_size(title=2) +
#   theme_sjplot() 
# 
# pCoeffs

################################################################################################
# Search decision types

dDecision = ddply(subset(d, !is.na(type_choice)), ~agegroup+id, plyr::summarize,
                  prob=prop.table(table(type_choice)),
                  type_choice=names(table(type_choice)))

dDecision = dDecision %>%
  mutate(type_choice = factor(type_choice, levels=c('Repeat', 'Near', 'Far')))

randBaseline = data.frame(decision = c('Repeat', 'Near', 'Far'),
                          prob = c(1/64, 8/64, 55/64),
                          group = rep('Random\nChoice\nModel', 3))

pDecisionType = ggplot() +
  stat_summary(dDecision, mapping=aes(x=type_choice, y=prob, color=agegroup, group=agegroup),
               fun.y=mean, geom='line') +
  stat_summary(dDecision, mapping=aes(x=type_choice, y=prob, color=agegroup, group=agegroup),
               fun.data=mean_cl_boot, show.legend=F) +
  geom_point(randBaseline, mapping=aes(x=decision, y=prob, group=group, fill=group), color='red', show.legend=FALSE) +
  geom_line(randBaseline, mapping=aes(x=decision, y=prob, group=group), linetype='dashed', color='red',  show.legend=FALSE) +
  scale_y_continuous() +
  xlab('Search Distance') +
  ylab('P(Search Distance)') +
  labs(color='Age Group', fill='') +
  scale_color_viridis(discrete=TRUE, direction=1) +
  ggtitle('Search Distance') +
  theme_classic() +
  theme(strip.background=element_blank(), legend.position='none')

pDecisionType

# variance over options
# don't count first randomly revealed tile as unique option
dUniqueOpts = data.frame()
for (i in unique(d$id)) {
  for (r in unique(subset(d, id==i)$round)) {
    subd = subset(d, id==i & round==r)
    
    uniqueOpts = unique(subset(subd, trial!=0)$chosen)
    firstTile = subset(subd, trial==0)$chosen
    nUniqueOpts = ifelse(firstTile %in% uniqueOpts, length(uniqueOpts)-1, length(uniqueOpts))
    
    dat = data.frame(id=i, round=r, age_months=unique(subd$age_months),
                     age_years=unique(subd$age_years), agegroup=unique(subd$agegroup),
                     nUniqueOpt=nUniqueOpts)
    dUniqueOpts = rbind(dUniqueOpts, dat)
  }
}

dMeanUniqueOpts = ddply(dUniqueOpts, ~id+age_months+agegroup+age_years, plyr::summarize,
                        mUniqueOpt=mean(nUniqueOpt))
randomUnique<- mean(sapply(1:10000, FUN=function(i) length(unique(sample(1:64, 25, replace = T)))))

ttestPretty(subset(dMeanUniqueOpts, agegroup=='5-6')$mUniqueOpt, mu = randomUnique)
ttestPretty(subset(dMeanUniqueOpts, agegroup=='5-6')$mUniqueOpt, mu = 25)

pUniqueOpts <- ggplot(dMeanUniqueOpts, aes(x=age_months/12, y=mUniqueOpt, color = agegroup)) +
  geom_point(aes(group=id), alpha=0.5, size =0.5) +
  geom_smooth(color = '#377eb8', fill = '#377eb8', size=.5) +
  geom_hline(yintercept=randomUnique, linetype='dashed', color='red') +
  coord_cartesian(ylim=c(0,25))+
  scale_x_continuous(trans="log10") +
  scale_color_viridis(discrete=TRUE, direction=1) +
  xlab('Age (Years) [logscale]') +
  ylab('Mean Unique Options') +
  ggtitle('Unique Options per Round') +
  theme_classic() +
  theme(strip.background=element_blank(), axis.line = element_line(),
        axis.ticks = element_line(), legend.position='none')

pUniqueOpts


# combine plots
plots = cowplot::plot_grid(pReward, pLearningCurves, pLearningCurvesMax, pUniqueOpts, pDecisionType,
                           pDistRewardRegression, ncol=2, labels='auto')
plots
# ggsave(filename = "plots/Fig2.pdf", plot=plots, height=8, width=6, units = "in")


################################################################################################
# Learning over rounds
################################################################################################
# hierarchical regression model
avgRewardRound = run_model(brm(meanReward ~ round * agegroup + (round + agegroup | id),
                               data = dLearningCurvesRounds, cores = 4, iter = 4000, warmup = 1000,
                               control = list(adapt_delta = 0.99)), modelName = 'AvgRewardRound')

# generate predictions
round = seq(2,9) # normalized reward
agegroup = levels(d$agegroup)[1:6]
d1 = expand.grid(round=round, agegroup=agegroup)
# agegroup 5-6 only played 4 rounds
round = seq(2,5)
agegroup = levels(d$agegroup)[7]
d2 = expand.grid(round=round, agegroup=agegroup)
# combine
newdat = rbind(d1, d2)
# predict distance based on previous resward
preds = fitted(avgRewardRound, re_formula=NA, newdata=newdat, probs=c(.025, .975))
predsDF = data.frame(round=newdat$round,
                     agegroup=newdat$agegroup,
                     reward=preds[,1],
                     lower=preds[,3],
                     upper=preds[,4])


dLearningCurvesRounds$agegroup = factor(dLearningCurvesRounds$agegroup, levels=c("25-55", "18-24", "14-17", "11-13", "9-10", "7-8", "5-6"))
predsDF$agegroup = factor(predsDF$agegroup, levels=c("25-55", "18-24", "14-17", "11-13", "9-10", "7-8", "5-6"))

# plot predictions
pAvgRewardRound <- ggplot()+
  # geom_point(dLearningCurvesRounds, mapping=aes(x=round, y=meanReward, color=agegroup, group=id), alpha=0.7, size=.5)+
  stat_summary(dLearningCurvesRounds, mapping=aes(x=round, y=meanReward, color=agegroup, fill=agegroup), fun=mean, geom='point', alpha=0.7, size=.5)+
  geom_line(predsDF, mapping=aes(x=round, y=reward, color=agegroup), size=.8) +
  geom_ribbon(predsDF, mapping=aes(x=round, y=reward, ymin=lower, ymax=upper, fill=agegroup), alpha=.3) +
  geom_hline(yintercept=0.5, linetype='dashed', color='red') + # random choice model
  xlab('Round')+
  ylab('Normalized Reward')+
  # scale_x_continuous(limits=c(1.5,9.5), breaks=seq(2,9)) +
  scale_color_viridis(discrete=TRUE, name='Age Group') +
  scale_fill_viridis(discrete=TRUE, name='Age Group') +
  ggtitle('Reward ~ Round') +
  #facet_wrap(~agegroup, ncol=4) +
  theme_classic() +
  theme(strip.background=element_blank())

pAvgRewardRound
