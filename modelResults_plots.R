# Plot model results
# Anna Giron, Simon Ciranka, Charley Wu, 2022

# house keeping
rm(list=ls())


#setwd("./Projects/adolescent_grids/")####THIS IS ONLY BC SIMON HAS HIS PROJECT IN HIS FOLDER 
# load packages
packages <- c('plyr', 'jsonlite', 'gridExtra', 'reshape2', 'stargazer', 'coefplot', 'cowplot', 'lme4', 'sjPlot',
              "grid", 'corrplot', 'ggbeeswarm', 'tidyverse', 'viridis', 'colorspace', 'ggrepel', 'pacman', 'tidybayes','brms','philentropy', 'latex2exp')
lapply(packages, require, character.only = TRUE)

source("dataProcessing.R")
source('statisticalTests.R')

dropLeadingZero <- function(l){
  str_replace(l, '0(?=.)', '')
}


##################################################################################################################
modelPal <- c('black', '#6c7195', '#ea9d67', '#7ec3aa')
paramPal = c("#FFEE67", '#27AD88', "#D1495B")
#############################################################################################################################
# Import & Rsq model comparison.
#############################################################################################################################

params = read.csv('data/modelFit.csv') #all mode results compiled together by running paramImport() from dataProcessing.R

params = params %>%
  mutate(kernel=factor(kernel, levels=c('RBF', 'BMT'), labels=c('GP', 'BMT'))) %>%
  mutate(agegroup=factor(agegroup,
                         levels=rev(c('[25,55]', '[18,25)', '[14,18)', '[11,14)', '[9,11)', '[7,9)', '[5,7)')),
                         labels=rev(c('25-55', '18-24', '14-17', '11-13', '9-10', '7-8', '5-6')))) %>%
  mutate(experiment=factor(experiment, levels=c('Meder (2021)', 'Schulz (2019)', 'Adolescent')))
params$ModelName = paste(params$kernel, params$acq, sep="-")
params$ModelName = factor(params$ModelName, levels = c('GP-UCB', 'BMT-UCB', 'GP-GM', 'BMT-GM', 'GP-GV', 'BMT-GV', 'GP-EG', 'BMT-EG'))
params$acq <- factor(params$acq, levels = c('UCB', 'GM', 'GV', 'EG'))

#Only include key comparisons
params <- subset(params, ModelName %in% c('GP-UCB', 'GP-GM', 'BMT-UCB', 'GP-EG'))

#Two line name for models
params$shortname <- factor(params$ModelName, levels = c('GP-UCB','BMT-UCB', 'GP-GM', 'GP-EG'))
levels(params$shortname) <- c('GP\nUCB', 'lambda\nlesion', 'beta\nlesion', 'tau\nlesion')


params%>%ggplot(aes(x=shortname, y=R2, fill=NA,color=shortname)) +
  #geom_line(aes(group=id), color = 'grey', alpha  = 0.3)+
  geom_quasirandom( size = 0.5)+
  geom_boxplot(width = 0.4, color ='black', outlier.shape=NA, fill = NA)+
  stat_summary(fun.y = mean, geom='point', shape = 23, color = 'black', fill = 'white')+
  xlab('') +
  ylab(expression(R^2)) +
  #ylab(expression(italic(pxp))) +
  scale_color_manual(values=modelPal, name = 'Model', labels = expression('GP-UCB', lambda*' lesion', beta* ' lesion', tau*' lesion')) +
  scale_fill_manual(values=modelPal, name = 'Model', labels = expression('GP-UCB', lambda*' lesion', beta* ' lesion', tau*' lesion')) +
  facet_wrap(~agegroup, nrow = 2)+
  ggtitle('Supplemental Model Comparison') +
  theme_classic() +
  theme(strip.background=element_blank(),
        legend.position = c(.9,0.1), legend.justification = c(1,0), axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())->p_R2_comp 

p_R2_comp
ggsave(filename = 'plots/S2.pdf',p_R2_comp )


#############################################################################################################################
# Plot model fit by age
#############################################################################################################################

# age as continuous variable
p_r2_age = ggplot(params, mapping=aes(x=age_months/12, y=R2, color=shortname, fill=shortname)) +
  geom_hline(yintercept = 0, color = 'red', linetype = 'dashed')+
  geom_point(alpha = 0.5, size=.5) +
  geom_smooth(size=.7) +
  xlab("Age (Years) [logscale]") +
  ylab(expression(R^2)) +
  #facet_grid(experiment~., scales = 'free')+
  scale_x_continuous(trans="log10") +
  annotation_logticks(sides = "b")   +
  scale_fill_manual(values = modelPal, name='', labels = expression('GP-UCB', lambda*' lesion', beta* ' lesion', tau*' lesion')) +
  scale_color_manual(values =modelPal, name='', labels = expression('GP-UCB', lambda*' lesion', beta* ' lesion', tau*' lesion')) +
  ggtitle(expression(paste(R^2, " by Age"))) +
  theme_classic() +
  theme(legend.position=c(0,1.1), legend.justification=c(0,1),
        legend.key=element_blank(), legend.background=element_blank(),strip.background=element_blank(), legend.direction="horizontal") 

p_r2_age


#############################################################################################################################
# PXP #computed in PXP.ipynb
#############################################################################################################################
#Save model nLLs for running pxp
# #all data
# nLLs <- params %>% arrange(shortname) %>% pull(nLL) %>% matrix(ncol = 4)
# write.table(nLLs, file = 'modelResults/pxp/nLL.csv', sep=',', row.names=F,col.names = F)
# #each age group separately
# nLLs <- subset(params, agegroup == "5-6") %>% arrange(shortname) %>% pull(nLL) %>% matrix(ncol = 4)
# write.table(nLLs, file = 'modelResults/pxp/nLL5-6.csv', sep=',', row.names=F,col.names = F)
# 
# nLLs <- subset(params, agegroup == "7-8") %>% arrange(shortname) %>% pull(nLL) %>% matrix(ncol = 4)
# write.table(nLLs, file = 'modelResults/pxp/nLL7-8.csv', sep=',', row.names=F,col.names = F)
# 
# nLLs <- subset(params, agegroup == "9-10") %>% arrange(shortname) %>% pull(nLL) %>% matrix(ncol = 4)
# write.table(nLLs, file = 'modelResults/pxp/nLL9-10.csv', sep=',', row.names=F,col.names = F)
# 
# nLLs <- subset(params, agegroup == "11-13") %>% arrange(shortname) %>% pull(nLL) %>% matrix(ncol = 4)
# write.table(nLLs, file = 'modelResults/pxp/nLL11-13.csv', sep=',', row.names=F,col.names = F)
# 
# nLLs <- subset(params, agegroup == "14-17") %>% arrange(shortname) %>% pull(nLL) %>% matrix(ncol = 4)
# write.table(nLLs, file = 'modelResults/pxp/nLL14-17.csv', sep=',', row.names=F,col.names = F)
# 
# nLLs <- subset(params, agegroup == "18-24") %>% arrange(shortname) %>% pull(nLL) %>% matrix(ncol = 4)
# write.table(nLLs, file = 'modelResults/pxp/nLL18-24.csv', sep=',', row.names=F,col.names = F)
# 
# nLLs <- subset(params, agegroup == "25-55") %>% arrange(shortname) %>% pull(nLL) %>% matrix(ncol = 4)
# write.table(nLLs, file = 'modelResults/pxp/nLL25-55.csv', sep=',', row.names=F,col.names = F)


# protected exceedance probability 
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
colnames(pxp) = c('GP-\nUCB', 'lambda\nlesion', 'beta\nlesion', 'tau\nlesion')
pxp$agegroup = rev(c('25-55', '18-24', '14-17', '11-13', '9-10', '7-8', '5-6', 'Overall'))
pxp$agegroup = factor(pxp$agegroup, levels=rev(c('Overall', '25-55', '18-24', '14-17', '11-13', '9-10', '7-8', '5-6')))
pxp = gather(pxp, ModelName, pxp, `GP-\nUCB`:`tau\nlesion`)
pxp$ModelName <- factor(pxp$ModelName, levels =  c('GP-\nUCB', 'lambda\nlesion', 'beta\nlesion', 'tau\nlesion') )

xlabels <- expression('GP-UCB', lambda * " lesion", beta * " lesion", tau * " lesion")
p_pxp_all = ggplot(pxp, aes(x=ModelName, y=pxp, fill=agegroup)) +
  geom_bar(stat='identity', position="dodge2", color = 'black') +
  xlab('') +
  ylab(expression(italic(pxp))) +
  scale_fill_manual(values=rev(c('black',"#440154FF", "#443A83FF", "#31688EFF" ,"#21908CFF", "#35B779FF", "#8FD744FF", "#FDE725FF")), name = 'Age Group') +
  scale_x_discrete(labels =xlabels)+
  ggtitle("Bayesian Model Selection") +
  theme_classic() +
  theme(legend.position= 'right',strip.background=element_blank(), 
  )

p_pxp_all


#############################################################################################################################
# Learning curves
#############################################################################################################################

# Simulated learning curves
path = 'rationalModels/'
filenames = list.files(path=path, pattern='*.csv')
filenames = paste0(path, filenames)

rationalDF = ldply(filenames, read.csv)
rationalDF = rationalDF[,!names(rationalDF) %in% 'X']

# normalize reward
rationalDF$meanReward = rationalDF$meanReward / 50
rationalDF$meanSE = rationalDF$meanSE / 50

nagegroups = length(levels(params$agegroup))
ntrials = 26

#add human data
behavior = read.csv('data/behavioralData.csv')

# only smooth condition
behavior = subset(behavior, condition=='Smooth')

behavior = behavior %>%
  mutate(type_choice=factor(type_choice, levels=c('Repeat', 'Near', 'Far'))) %>%
  mutate(experiment=factor(experiment, levels=c('Meder (2021)', 'Schulz (2019)', 'Adolescent'))) %>%
  mutate(agegroup=factor(agegroup,
                         levels=rev(c('[25,55]', '[18,25)', '[14,18)', '[11,14)', '[9,11)', '[7,9)', '[5,7)')),
                         labels=rev(c('25-55', '18-24', '14-17', '11-13', '9-10', '7-8', '5-6'))))

# append mean reward to model results
meanRew = ddply(behavior, ~id, plyr::summarize, meanReward=mean(z))


# random model performance should be displayed in all facets
random = subset(rationalDF, model=='Random')

randomModel = data.frame(trial=rep(random$trial, times=nagegroups),
                         meanReward=rep(random$meanReward, times=nagegroups),
                         meanSE=rep(random$meanSE, times=nagegroups),
                         agegroup=rep(levels(behavior$agegroup), each=ntrials),
                         model=rep('Random', times=nagegroups*ntrials))

rationalDF = subset(rationalDF, model!='Random')


# include human data
# normalize reward
behavior$z = behavior$z / 50
dplot = ddply(behavior, ~agegroup+trial, plyr::summarize, meanReward=mean(z), meanSE=sd(z)/sqrt(length(z)))
dplot$model = 'Human'

# combine datasets
rationalDF = rbind(rationalDF, randomModel, dplot)

# reorder factor levels, random model has no agegroup
rationalDF$agegroup = factor(rationalDF$agegroup, levels=c('5-6', '7-8', '9-10', '11-13', '14-17', '18-24', '25-55'),
                             labels=c('5-6', '7-8', '9-10', '11-13', '14-17', '18-24', '25-55'))
rationalDF$model = factor(rationalDF$model, levels=c('Human', 'GP-UCB', 'BMT-UCB',"GP-EG","GP-GM", 'Random'))

# plot
LCmodelPal <- c('grey','black', '#6c7195', '#ea9d67', '#7ec3aa','red')


p_LearningCurves = ggplot(rationalDF, aes(x=trial, y=meanReward, color=model, fill=model)) +
  geom_line() +
  geom_ribbon(aes(ymin=meanReward-meanSE, ymax=meanReward+meanSE), alpha=.3, color=NA) +
  facet_wrap(~agegroup, ncol=4) +
  xlab('Trial')+
  ylab('Normalized Mean Reward ± SEM')+
  ggtitle('Learning Curves') +
  labs(color='') +
  labs(fill='') +
  scale_fill_manual(values = c("grey",modelPal,"red"), breaks=c('Human', 'GP-UCB', 'BMT-UCB',"GP-GM","GP-EG", 'Random'), name='', labels = expression('Human','GP-UCB', lambda*' lesion', beta* ' lesion', tau*' lesion','random')) +
  scale_color_manual(values =c("grey",modelPal,"red"), breaks=c('Human', 'GP-UCB', 'BMT-UCB',"GP-GM","GP-EG", 'Random'),name='', labels = expression('Human','GP-UCB', lambda*' lesion', beta* ' lesion', tau*' lesion','random')) +
  theme_classic() +
  theme( strip.background=element_blank(),legend.position=c(.87,.2))

p_LearningCurves



#########################################################
#combine plots
#########################################################
# plots = cowplot::plot_grid(p_r2_ucb, p_pxp_all, p_LearningCurves, ncol=3)
# ggsave(filename = "plots/modelResults.pdf", plot=plots, height=6, width=14, units = "in")


# insetted<-ggdraw(p_pxp_main+ theme(axis.text.x=element_text(angle=0, hjust=0.5)) + ggtitle('Model comparison')) +
#   draw_plot(p_R2_comp+
#               ggtitle('')+
#               theme(text = element_text(size = 10),
#                     plot.background = element_rect(fill = "transparent",colour = NA)), .35, .2, .65, .8)
# insetted

p_ab <- cowplot::plot_grid(p_pxp_all, p_r2_age,labels = 'auto', ncol = 1)
plotsTop = cowplot::plot_grid(p_ab,p_LearningCurves, ncol=2, rel_widths = c(1, 1.2), labels = c('', 'c'))
plotsTop



##########################################
#Parameter plots
##########################################

#Everything after this is only GP-UCB
params = subset(params, kernel=='GP' & acq=='UCB')%>%mutate(
  agegroup = fct_relevel(agegroup,rev(c('25-55', '18-24', '14-17', '11-13', '9-10', '7-8', '5-6')))
)

params_old<-params


#Pauls function to save & reuse fits
run_model <- function(expr, path, reuse = TRUE) {
  path <- paste0(path, ".Rds")
  if (reuse) {
    fit <- suppressWarnings(try(readRDS(path), silent = TRUE))
  }
  if (is(fit, "try-error")) {
    fit <- eval(expr)
    saveRDS(fit, file = path)
  }
  fit
}


bprior <- c(
  set_prior("normal(0, 2)", resp = c("lambda","beta","tau"),nlpar = "b0"),
  # prior(exponential(0.1), nlpar = "cpunc")+
  set_prior("normal(0, 1)", resp = c("lambda","beta","tau"),nlpar = "b1"),
  set_prior("normal(0, 1)", resp = c("lambda","beta","tau"),nlpar = "b2"),
  set_prior("normal(0, 1)", resp = c("lambda","beta","tau"),nlpar = "alpha")
)# not entirely sure what a good prior is here

# inv logit looses the step function a little
bform <- bf(
  mvbind(lambda,beta,tau) ~ b0 + b1 * (age_years - omega) * step((omega - age_years)) + 
    b2 * (age_years - omega) * step((age_years - omega)),
  b0 + b1 + b2 + alpha ~ 1,
  # to keep omega within the age range of 5 to 25
  nlf(omega ~ 5+inv_logit(alpha) * 20),
  nl = TRUE
)

# rescore here because brms does something wierd to the names that i dont understand.
cpmodelparams<-params%>%mutate(
  lambda=log(lambda),
  beta=log(beta),
  tau=log(tau)
)

multiChange<-run_model(brm(bform, data = cpmodelparams,prior = bprior,cores=4,iter=4000,control=list(adapt_delta=0.99)),path = "brms_modelfits/multi_Change")


#######################################
# Plot change models
#######################################

lambdaParams <- params%>%filter(ModelName=='GP-UCB')%>%cbind(t(posterior_predict(multiChange,newdata=params,nsamples = 400)[,,1]))%>%
  pivot_longer(-colnames(params_old))

#scnd dim is beta
betaParams<-  params%>% filter(ModelName=='GP-UCB')%>%cbind(t(posterior_predict(multiChange,newdata=params,nsamples = 400)[,,2]))%>%
  pivot_longer(-colnames(params_old))

# betaChange%>%spread_draws(b_alpha_Intercept)%>%rowwise()%>%
#   mutate(hm=inv_logit_scaled(b_alpha_Intercept,lb = 0,ub=20))

#third is tau
tauParams<- params%>% filter(ModelName=='GP-UCB')%>%cbind(t(posterior_predict(multiChange,newdata=params,nsamples = 400)[,,3]))%>%
  pivot_longer(-colnames(params_old))

#SC: Here the actual model predictions have been thrown away before but we need to keep them if we want to plot them
allParams<-lambdaParams%>%magrittr::set_colnames(value = c(colnames(lambdaParams)[1:length(colnames(lambdaParams))-1],"lambda_pred"))
allParams$beta_pred <- betaParams$value
allParams$tau_pred <- tauParams$value

# for some reason the factor doesnt work otherwise
allParamsPreds<-lambdaParams
allParamsPreds$lambda<-lambdaParams$value
allParamsPreds$beta <- betaParams$value
allParamsPreds$tau <- tauParams$value
#combine all into a single dataframe

posteriorEstimates <-allParamsPreds %>% select(lambda, beta, tau,age_years, name) %>% 
  pivot_longer(cols = c(lambda, beta, tau), names_to = 'param_pred')


# posteriorEstimates$param_pred <- factor(posteriorEstimates$param_pred, levels=c('lambda', 'beta', 'tau'))
posteriorEstimates$param_pred <- factor(posteriorEstimates$param_pred, levels=c('lambda', 'beta', 'tau'),
                                        labels=c(expression(paste('Generalization ', lambda)),
                                                 expression(paste('Exploration ', beta)),
                                                 expression(paste('Temperature ', tau))))

#Posterior draws
#posteriorParamDraws <- posteriorEstimates %>% group_by(name, age_years, param) %>% summarize(value = mean(value))

# paramPal = c("#FFEE67", '#27AD88', "#D1495B", "#FFEE67", '#27AD88', "#D1495B")
#GP-UCB params as raw data
gpucbParams <- subset(params, ModelName == 'GP-UCB')
gpucbParams <- gpucbParams %>% pivot_longer(cols =c(lambda, beta, tau), names_to='param_pred')
# gpucbParams$param_pred <- factor(gpucbParams$param_pred, levels = c('lambda', 'beta', 'tau'))
gpucbParams$param_pred <- factor(gpucbParams$param_pred, levels=c('lambda', 'beta', 'tau'),
                                 labels=c(expression(paste('Generalization ', lambda)),
                                          expression(paste('Exploration ', beta)),
                                          expression(paste('Temperature ', tau))))

posteriorParamPlots<-posteriorEstimates%>%#mutate(log(value))
  ggplot()+
  stat_summary(aes(x=age_years,y=exp(value)),fun.min = function(z) { quantile(z,0.25) }, fun.max = function(z) { quantile(z,0.75) },
               fun = mean,geom="ribbon",fill="grey",alpha=0.5)+
  geom_point(data=gpucbParams,aes(y=value,x=age_months/12,color=param_pred),alpha=0.3,size=0.3)+
  ##stat_summary(data=allParamsLong,aes(y=log(value),x=age_years,color=param_pred),geom="point",size=1,shape=)+
  #geom_smooth(aes(y=exp(value),group=name,x=age_years,color=param_pred),size=0.04,alpha=0.0001,fill=NA)+
  #geom_smooth(aes(y=exp(value),x=age_years,color=param_pred),size=1,fill=NA)+
  stat_summary(fun = mean, aes(y=exp(value),x=age_years,color=param_pred),geom="line",size=1)+
  coord_cartesian(xlim=c(min(params$age_years),max(params$age_years)))+
  #scale_y_log10(breaks = c(.01, .1, 1.0, 10), labels = c(".01", '0.1', '1.0', '10'))+
  scale_y_log10(labels = dropLeadingZero)+
  scale_x_log10()+
  facet_grid(~param_pred, labeller=label_parsed, scales='free_y')+
  scale_color_manual(values = paramPal)+
  #xlab("Age (Years) [logscale]") +
  xlab("") +
  ggtitle('GP-UCB Parameters')+
  ylab('Estimate')+
  theme_classic() +
  theme(strip.background=element_blank(),
        legend.position='none')

posteriorParamPlots      


#######Histogram of change point

lambdaHist <- multiChange%>%spread_draws(b_lambda_alpha_Intercept)%>%rowwise()%>%
  mutate(hm=inv_logit_scaled(b_lambda_alpha_Intercept,lb = 5,ub=25))%>%
  select(-b_lambda_alpha_Intercept)

lambdaHist$param <- 'lambda'
#when changepoint? uncertainty?
lambdaHist%>%ungroup()%>%
  summarize(m=mean(hm),
            upper_CI=quantile(hm,probs = 0.95),
            lower_CI=quantile(hm,probs = 0.05))#%>%

betaHist <- multiChange%>%spread_draws(b_beta_alpha_Intercept)%>%rowwise()%>%
  mutate(hm=inv_logit_scaled(b_beta_alpha_Intercept,lb = 5,ub=25))%>%
  select(-b_beta_alpha_Intercept)

betaHist$param <- 'beta'

betaHist%>%ungroup()%>%summarize(m=mean(hm),
                                 upper_CI=quantile(hm,probs = 0.95),
                                 lower_CI=quantile(hm,probs = 0.05))

tauHist <- multiChange%>%spread_draws(b_tau_alpha_Intercept)%>%rowwise()%>%
  mutate(hm=inv_logit_scaled(b_tau_alpha_Intercept,lb = 5,ub=25))%>%
  select(-b_tau_alpha_Intercept)

tauHist$param <- 'tau'

tauHist%>%ungroup()%>%summarize(m=mean(hm),
                                upper_CI=quantile(hm,probs = 0.95),
                                lower_CI=quantile(hm,probs = 0.05))

allHist <- rbind(lambdaHist, betaHist,tauHist )
allHist$param <- factor(allHist$param, levels = c('lambda', 'beta', 'tau'))

changeHist <- ggplot(allHist, aes(x=hm, fill = param))+
  geom_histogram(bins=30, color = 'black',  aes(y = stat(width*density)), boundary = 0.5)+
  #geom_density()+
  coord_cartesian(xlim=c(min(params$age_years),max(params$age_years)))+
  scale_x_log10(name="Age (Years) [logscale]")+
  #scale_y_continuous(name=expression(paste(omega, '\n(posterior density)')))+
  scale_y_continuous(name=expression(omega))+
  scale_fill_manual(values = paramPal)+
  facet_wrap(~param, labeller = label_parsed,  scales='free_y')+
  theme_classic() +
  theme( strip.background=element_blank(), legend.position='none', 
         strip.text.x = element_blank())
changeHist


ParamPlot <- cowplot::plot_grid(posteriorParamPlots+ theme(plot.margin = unit(c(7, 0, 7, 7), "pt")), 
                                changeHist+ theme(plot.margin = unit(c(-10, 7, 7, 7), "pt")), ncol=1, rel_heights =c(1, .5))
ParamPlot
#ggsave(ParamPlot,filename = "CP_Model_revision.png",width = 7,height=4)




#######################################################################################
# Parameter correlations
# distance matrix for each age group
########################################################################################


#Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
#map(unique(.$age_years), ~{x})
ageGroupList<-list()

# parse parametercorrs in different agegroups
for (i in 1:length(unique(params$age_years))){
  ageGroupList[[i]]<-params%>%filter(kernel=="GP",acq=="UCB")%>%
    filter(age_years==unique(params$age_years)[i])%>%
    select(tau,lambda,beta)%>%as_tibble()
}

#Confidence interval functions
lower_ci <- function(mean, se, n, conf_level = 0.95){
  lower_ci <- mean - qt(1 - ((1 - conf_level) / 2), n - 1) * se
}

upper_ci <- function(mean, se, n, conf_level = 0.95){
  upper_ci <- mean + qt(1 - ((1 - conf_level) / 2), n - 1) * se
}

########################
#Compute pairwise distance
#Very SLOW
#########################
params<-params%>%arrange(age_years)#important to sort the age groups first.

#Takes a little time
# dissimilarity <- sapply(1:length(params$age_years), FUN=function(x){
#   sapply(1:length(params$age_years), FUN=function(y){
#     vec1<-c(params$tau[x],params$beta[x],params$lambda[x])
#     vec2<-c(params$tau[y],params$beta[y],params$lambda[y])
#     1-cor(vec1,vec2,method="kendall")
#   })
# })
#saveRDS(object = dissimilarity, file = 'data/dissimMatrix.rds')
dissimilarity <- readRDS('data/dissimMatrix.rds') #Pre-saved computations from the commented out block above

comparisonGrid <- expand.grid(from=1:length(params$age_years), to = 1:length(params$age_years))
dissimilarityDF <- data.frame(ageX = params$agegroup[comparisonGrid$from], ageY=params$agegroup[comparisonGrid$to], dissim = as.numeric(dissimilarity))

#prepare correlation plot 
dissimilarityDF %>% group_by(ageX, ageY)%>%
  summarise(mdiss=mean(dissim), sd = sd(dissim), n = length(dissim)) %>%
  mutate(se = sd / sqrt(n),
         lower_ci = lower_ci(mdiss, se, n),
         upper_ci = upper_ci(mdiss, se, n))->getTria

getTria$metric <- 'kendall'
#keep agegroups after melting the correlation matrix with the reshape package. just a lookup vector
agegroupLookup<-params%>%select(age_years,agegroup)%>%unique()%>%arrange(age_years)


dissimData <- lapply(split(getTria,getTria$metric),#turns them into lists and reshapes the list elements individually
                     function(x){dcast(x,ageX~ageY, value.var = 'mdiss')%>%select(-ageX)%>%as.matrix()%>%
                         get_lower_tri()%>%melt()%>%mutate(X3=unique(agegroupLookup$agegroup)[Var1])%>%
                         mutate(
                           X3 = fct_relevel(X3,rev(c('25-55', '18-24', '14-17', '11-13', '9-10', '7-8', '5-6'))),
                           X2 = fct_relevel(Var2,rev(c('25-55', '18-24', '14-17', '11-13', '9-10', '7-8', '5-6')))
                         )
                     }
)

p_SimilarityMetrics <- lapply(seq_along(dissimData),
                              function(i){#why not a forloop you ask? I DONT KNOW
                                ggplot(dissimData[[i]], aes(x = X3,y=X2,fill=value))+
                                  geom_tile()+
                                  scale_y_discrete(name="Age group")+
                                  scale_x_discrete(name="Age group")+
                                  scale_fill_viridis(option="viridis",na.value="white", name = "distance")+
                                  theme_classic() +
                                  ggtitle(names(dissimData)[[i]])+
                                  theme( axis.text.x=element_text(angle = 45,hjust=1),
                                         legend.position = 'right', strip.background=element_blank(),
                                         legend.key=element_blank(),legend.background=element_blank())})

###############
#create plot 
###############

p_Similarity<-ggplot(dissimData[[1]], aes(x = X3,y=X2,fill=1-value))+
  geom_tile()+
  scale_y_discrete(name="Age group")+
  scale_x_discrete(name="Age group")+
  scale_fill_viridis(option="plasma",na.value="white", name = expression(r[tau]))+
  theme_classic() +
  ggtitle("Parameter Similarity")+
  theme( axis.text.x=element_text(angle = 45,hjust=1),
         legend.position = 'right', strip.background=element_blank(),
         legend.key=element_blank(),legend.background=element_blank())
p_Similarity




getTria %>% filter(ageX==ageY)%>%mutate(mdiss=1-mdiss,lower_ci=1-lower_ci, upper_ci=1-upper_ci)%>%
  ggplot(aes(x=ageX, y = mdiss))+
  #geom_line(aes(group=1), color  = 'black')+
  #geom_smooth(method='lm', color = 'black', se = FALSE)+
  geom_errorbar(aes(ymin=lower_ci, ymax = upper_ci, group=1), color  = 'black', width = 0.2)+
  geom_point(aes(group=1))+
  theme_classic()+
  scale_y_continuous(labels = dropLeadingZero)+
  xlab('Age group')+
  ylab(expression(r[tau]))+
  ggtitle('Diagonals')+
  theme(legend.position='none') ->p_selfSimilarity

p_selfSimilarity

#Combine together using an inset
insetSimilarity<-ggdraw(p_Similarity) +
  draw_plot(p_selfSimilarity+
              xlab('')+
              ggtitle('')+
              annotate("text",label= "Diagonals", x= 2.5, y = .73,size = 3 )+
              theme(text = element_text(size = 10),
                    plot.title=element_text(margin=margin(b=-12), hjust=0.1),
                    axis.text.x=element_text(angle = 45,hjust=1),
                    plot.background = element_rect(fill = "transparent",colour = NA)), 0.16, .53, .38, .4)
insetSimilarity

#############################
#Create final plot
##########################

bottomRow <- cowplot::plot_grid(ParamPlot, insetSimilarity ,nrow = 1, rel_widths = c( 1,.7), labels=c('d','e'))

fullPlot <- cowplot::plot_grid(plotsTop, bottomRow, ncol = 1, rel_heights = c(1,.7))
fullPlot

#ggsave('plots/models_27_06_2022.png',fullPlot, width = 10, height = 7.5, units = 'in' )
#ggsave('plots/models_28_06_2022.pdf',fullPlot, width = 12, height = 9, units = 'in' )

