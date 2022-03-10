# Plot simulation of models
# Anna Giron, 2021

rm(list=ls())

packages <- c('plyr', 'dplyr', 'tidyr', 'ggplot2', 'cowplot', 'viridis', 'colorspace')
invisible(lapply(packages, require, character.only = TRUE)) #loads packages
#invisible(lapply(packages, install.packages, character.only = TRUE))

source('dataProcessing.R')

################################################################################################
# Data import
################################################################################################
# model results
params = read.csv('data/modelFit.csv')
# only use parameter estimates from GP-UCB model
params = subset(params, kernel=='RBF' & acq=='UCB')

params = params %>%
  mutate(kernel=factor(kernel, levels=c('RBF', 'BMT'), labels=c('GP', 'BMT'))) %>%
  mutate(agegroup=factor(agegroup,
                         levels=rev(c('[25,55]', '[18,25)', '[14,18)', '[11,14)', '[9,11)', '[7,9)', '[5,7)')),
                         labels=rev(c('25-55', '18-24', '14-17', '11-13', '9-10', '7-8', '5-6'))))

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
params$meanReward = params$meanReward / 50


# simulation results
# path = 'modelResults/simulatedModels/batch1/'
# filenames = list.files(path=path, pattern='*.csv', full.names=TRUE)
#
# sim = lapply(1:length(filenames), FUN=function(x) {
#   cur = read.csv(filenames[x])
#   cur = ddply(cur, ~lambda+beta+tau, plyr::summarize, meanReward=mean(mu))
# })
# sim = as.data.frame(do.call(rbind, sim))
# 
# # save mean rewards
# write.csv(sim, 'modelResults/simulatedModels/meanRewards.csv')

sim = read.csv('modelResults/simulatedModels/meanRewards.csv')

# normalize reward
sim$meanReward = sim$meanReward / 50
sim = sim[,!names(sim) %in% 'X']


################################################################################################
# Plots - all parameter combinations
################################################################################################
################################################################################################
# faceted by tau
taus = sort(unique(sim$tau))
dTau = sim
dTau$tau = factor(dTau$tau, levels=taus, labels=round(taus, digits=4))

# reward continuous
pTau = ggplot(dTau, aes(x=lambda, y=beta, fill=meanReward)) +
  geom_raster(interpolate=TRUE) +
  scale_x_continuous(trans="log10", expand=c(0,0)) +
  scale_y_continuous(trans="log10", expand=c(0,0)) +
  facet_wrap(tau~., nrow=10, labeller=label_parsed) +
  xlab(expression(paste('Generalization ', lambda, ' [logscale]'))) +
  ylab(expression(paste('Exploration ', beta, ' [logscale]'))) +
  scale_fill_gradientn(colours=sequential_hcl(500, 'Inferno'), name='Normalized\nReward') +
  ggtitle(expression(paste('Simulated Reward (Faceted by Temperature ', tau, ')'))) +
  theme_classic() +
  theme(legend.position='right', strip.background=element_blank(), axis.text=element_text(size=6))

pTau

# ggsave(filename = "plots/S5.pdf", plot=pTau, height=8, width=7, units = "in")


################################################################################################
# faceted by lambda
lambdas = sort(unique(sim$lambda))
dLambda = sim
dLambda$lambda = factor(dLambda$lambda, levels=lambdas, labels=round(lambdas, digits=4))

# reward continuous
pLambda = ggplot(dLambda, aes(x=beta, y=tau, fill=meanReward)) +
  geom_raster(interpolate=TRUE) +
  scale_x_continuous(trans="log10", expand=c(0,0)) +
  scale_y_continuous(trans="log10", expand=c(0,0)) +
  facet_wrap(lambda~., nrow=10, labeller=label_parsed) +
  xlab(expression(paste('Exploration ', beta, ' [logscale]'))) +
  ylab(expression(paste('Temperature ', tau, ' [logscale]'))) +
  scale_fill_gradientn(colours=sequential_hcl(500, 'Inferno'), name='Normalized\nReward') +
  ggtitle(expression(paste('Simulated Reward (Faceted by Generlization ', lambda, ')'))) +
  theme_classic() +
  theme(legend.position='right', strip.background=element_blank(), axis.text=element_text(size=6))

pLambda

# ggsave(filename = "plots/S5.1.pdf", plot=pLambda, height=8, width=7, units = "in")


################################################################################################
# faceted by beta
betas = sort(unique(sim$beta))
dBeta = sim
dBeta$beta = factor(dBeta$beta, levels=betas, labels=round(betas, digits=4))


# reward continuous
pBeta = ggplot(dBeta, aes(x=lambda, y=tau, fill=meanReward)) +
  geom_raster(interpolate=TRUE) +
  scale_x_continuous(trans="log10", expand=c(0,0)) +
  scale_y_continuous(trans="log10", expand=c(0,0)) +
  facet_wrap(beta~., nrow=10, labeller=label_parsed) +
  xlab(expression(paste('Generalization ', lambda, ' [logscale]'))) +
  ylab(expression(paste('Temperature ', tau, ' [logscale]'))) +
  scale_fill_gradientn(colours=sequential_hcl(500, 'Inferno'), name='Normalized\nReward') +
  ggtitle(expression(paste('Simulated Reward (Faceted by Exploration ', beta, ')'))) +
  theme_classic() +
  theme(legend.position='right', strip.background=element_blank(), axis.text=element_text(size=6))

pBeta

# ggsave(filename = "plots/S5.2.pdf", plot=pBeta, height=8, width=7, units = "in")
