# Analyze behavioral results
# Anna Giron, 2021

rm(list=ls())

packages <- c('plyr', 'dplyr', 'tidyr', 'ggplot2', 'cowplot', 'viridis', 'entropy', 'sjPlot')
invisible(lapply(packages, require, character.only = TRUE)) #loads packages
#invisible(lapply(packages, install.packages, character.only = TRUE))

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

################################################################################################
# Demographic information
################################################################################################
# inspect bins
d %>% 
  filter(round == 2 & trial == 0) %>% 
  group_by(agegroup) %>% 
  summarise(n = n(),
            female = sum(gender == "Female"),
            mean_age = mean(age_years),
            sd_age = sd(age_years),
            mean_duration = mean(duration)) 

d %>% 
  filter(round == 2 & trial == 0) %>% 
  group_by(experiment) %>% 
  summarise(n = n(),
            female = sum(gender == "Female"),
            mean_age = mean(age_years),
            sd_age = sd(age_years),
            mean_duration = mean(duration)) 


################################################################################################
# Performance
################################################################################################
# youngest age group better than chance level?
# age as discrete variable
dReward = ddply(d, ~id+agegroup+age_months, plyr::summarize, meanReward=mean(z))

ttestPretty(subset(dReward, agegroup=='5-6')$meanReward, mu=.5)


# reward as a function of age
corTestPretty(dReward$age_months, dReward$meanReward)


# Learning curves - highest reward
dLearningCurvesMax = data.frame()
for (i in unique(d$id)) {
  for (r in unique(subset(d, id==i)$round)) {
    # compute maximum reward found so far
    maxRews = lapply(seq(0, max(subset(d, id==i)$trial)), FUN=function(t) {
      max(subset(d, id==i & round==r & trial<=t)$z)
    })
    
    trials = length(maxRews)
    subd = subset(d, id==i & round==r)
    
    dat = data.frame(id=rep(i, trials), round=rep(r, trials), trial=seq(1, trials), age_months=subd$age_months,
                     age_years=subd$age_years, agegroup=subd$agegroup,
                     maxReward=unlist(maxRews))
    dLearningCurvesMax = rbind(dLearningCurvesMax, dat)
  }
}

# test if rewards are normally distributed
shapiro.test(subset(dLearningCurvesMax, agegroup=='5-6' & trial==25)$maxReward)
shapiro.test(subset(dLearningCurvesMax, agegroup=='7-8' & trial==25)$maxReward)
shapiro.test(subset(dLearningCurvesMax, agegroup=='9-10' & trial==25)$maxReward)
shapiro.test(subset(dLearningCurvesMax, agegroup=='11-13' & trial==25)$maxReward)
shapiro.test(subset(dLearningCurvesMax, agegroup=='14-17' & trial==25)$maxReward)
shapiro.test(subset(dLearningCurvesMax, agegroup=='18-24' & trial==25)$maxReward)
shapiro.test(subset(dLearningCurvesMax, agegroup=='25-55' & trial==25)$maxReward)


ranktestPretty(subset(dLearningCurvesMax, agegroup=='5-6' & trial==25)$maxReward,
               subset(dLearningCurvesMax, agegroup=='7-8' & trial==25)$maxReward, paired=FALSE)
# "$U=12290$, $p<.001$, $r_{\tau}=-.21$, $BF>100$"

ranktestPretty(subset(dLearningCurvesMax, agegroup=='7-8' & trial==25)$maxReward,
               subset(dLearningCurvesMax, agegroup=='9-10' & trial==25)$maxReward, paired=FALSE)
# "$U=38340$, $p=.031$, $r_{\tau}=-.08$, $BF=.96$"

ranktestPretty(subset(dLearningCurvesMax, agegroup=='9-10' & trial==25)$maxReward,
               subset(dLearningCurvesMax, agegroup=='11-13' & trial==25)$maxReward, paired=FALSE)
# "$U=52501$, $p=.019$, $r_{\tau}=-.08$, $BF=1.2$"

ranktestPretty(subset(dLearningCurvesMax, agegroup=='11-13' & trial==25)$maxReward,
               subset(dLearningCurvesMax, agegroup=='14-17' & trial==25)$maxReward, paired=FALSE)
# "$U=63761$, $p=.218$, $r_{\tau}=-.04$, $BF=.12$"

ranktestPretty(subset(dLearningCurvesMax, agegroup=='14-17' & trial==25)$maxReward,
               subset(dLearningCurvesMax, agegroup=='18-24' & trial==25)$maxReward, paired=FALSE)
# "$U=50773$, $p<.001$, $r_{\tau}=-.11$, $BF=7.7$"

ranktestPretty(subset(dLearningCurvesMax, agegroup=='18-24' & trial==25)$maxReward,
               subset(dLearningCurvesMax, agegroup=='25-55' & trial==25)$maxReward, paired=FALSE)
# "$U=42668$, $p=.625$, $r_{\tau}=-.02$, $BF=.11$"


# age as continuous variable
corTestPretty(subset(dLearningCurvesMax, trial==25)$age_months,
              subset(dLearningCurvesMax, trial==25)$maxReward,
              method='kendall')
# "$r_{\tau}=.23$, $p<.001$, $BF>100$"


################################################################################################
# Search Decision
################################################################################################
# distribution of search decision types
dDecision = ddply(subset(d, !is.na(type_choice)), ~agegroup+id+age_months, plyr::summarize,
                  prob=as.double(prop.table(table(type_choice))),
                  type_choice=names(table(type_choice)))

# repeat clicks of youngest participants different from random choice model?
ttestPretty(subset(dDecision, agegroup=='5-6' & type_choice=='Repeat')$prob, mu=1/64)
# near choices of youngest participants different from random choice model?
ttestPretty(subset(dDecision, agegroup=='5-6' & type_choice=='Near')$prob, mu=8/64)

# proportion of far choices over the lifespan
corTestPretty(subset(dDecision, type_choice=='Far')$age_months, subset(dDecision, type_choice=='Far')$prob, method='kendall')

# when does the rate of repeat and near choices reach parity?
ttestPretty(subset(dDecision, agegroup=='5-6' & type_choice=='Repeat')$prob,
            subset(dDecision, agegroup=='5-6' & type_choice=='Near')$prob, paired=TRUE)

ttestPretty(subset(dDecision, agegroup=='7-8' & type_choice=='Repeat')$prob,
            subset(dDecision, agegroup=='7-8' & type_choice=='Near')$prob, paired=TRUE)

ttestPretty(subset(dDecision, agegroup=='9-10' & type_choice=='Repeat')$prob,
            subset(dDecision, agegroup=='9-10' & type_choice=='Near')$prob, paired=TRUE)

ttestPretty(subset(dDecision, agegroup=='11-13' & type_choice=='Repeat')$prob,
            subset(dDecision, agegroup=='11-13' & type_choice=='Near')$prob, paired=TRUE)

ttestPretty(subset(dDecision, agegroup=='14-17' & type_choice=='Repeat')$prob,
            subset(dDecision, agegroup=='14-17' & type_choice=='Near')$prob, paired=TRUE)

ttestPretty(subset(dDecision, agegroup=='18-24' & type_choice=='Repeat')$prob,
            subset(dDecision, agegroup=='18-24' & type_choice=='Near')$prob, paired=TRUE)

ttestPretty(subset(dDecision, agegroup=='25-55' & type_choice=='Repeat')$prob,
            subset(dDecision, agegroup=='25-55' & type_choice=='Near')$prob, paired=TRUE)



# unique options
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

corTestPretty(dMeanUniqueOpts$age_months, dMeanUniqueOpts$mUniqueOpt)


# entropy
dEntropyRounds = data.frame()

for (i in unique(d$id)) {
  for (r in unique(subset(d, id==i)$round)) {
    subd = subset(d, id==i & round==r)
    entropy = entropy(table(factor(subd$chosen, levels=1:64)))
    dat = data.frame(id=i, age_years=unique(subd$age_years), age_months=unique(subd$age_months), agegroup=unique(subd$agegroup),
                     round=r, entropy=entropy)
    dEntropyRounds = rbind(dEntropyRounds, dat)
  }
}

dEntropy = ddply(dEntropyRounds, ~id+age_months+agegroup+age_years, plyr::summarize,
                 meanEntropy=mean(entropy))

corTestPretty(dEntropy$age_months, dEntropy$meanEntropy)


# Regression model: modulation of distance as a function of reward
DistancePrevReward = run_model(brm(distance ~ previous_reward * agegroup + (previous_reward + agegroup | id),
                                   data = d, cores = 4, iter = 4000, warmup = 1000,
                                   control = list(adapt_delta = 0.99), save_all_pars = T), modelName = 'DistancePrevReward')


DistancePrevRewardNull <- run_model(brm(distance ~ 1 + (previous_reward + agegroup | id),
                                        data = d,  cores = 4,  iter = 4000, warmup = 1000,
                                        control = list(adapt_delta = 0.99), save_all_pars = T), modelName = 'DistancePrevRewardNull')

# # Model comparison
# bayes_factor(DistancePrevReward, DistancePrevRewardNull)
# Estimated Bayes factor in favor of DistancePrevReward over DistancePrevRewardNull:     NA
# Warnmeldungen:
#   1: Infinite value in iterative scheme, returning NA.
# Try rerunning with more samples.
# 2: logml could not be estimated within maxiter, rerunning with adjusted starting value.
# Estimate might be more variable than usual.

# bayes_R2(DistancePrevReward)
# Estimate   Est.Error      Q2.5     Q97.5
# R2 0.403919 0.002604367 0.3987298 0.4090274

# Model results
fixef(DistancePrevReward)
tab_model(DistancePrevReward)
bayes_R2(DistancePrevReward)
# Estimate   Est.Error      Q2.5     Q97.5
# R2 0.403919 0.002604367 0.3987298 0.4090274
