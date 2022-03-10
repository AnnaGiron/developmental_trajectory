# Plot trajectories of stochastic hill climbing model
# Anna Giron, Charley Wu, Simon Ciranka, 2022

rm(list=ls())

packages <- c('plyr', 'dplyr', 'tidyr', 'ggplot2', 'cowplot', 'viridis', 'colorspace', 'devtools', 'ggnewscale', 'rmisc', 'ggrepel')
invisible(lapply(packages, require, character.only = TRUE)) #loads packages
# invisible(lapply(packages, install.packages, character.only = TRUE))
# devtools::install_github("AckerDWM/gg3D")

source('dataProcessing.R')

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
# Compare algorithms and temperature functions
################################################################################################
dTrajectories <- trajectories %>% group_by(id, coolingFunc, method, i) %>%
  summarize(lambda = mean(lambda), beta = mean(beta), tau = mean(tau), reward = mean(reward))

mTrajectories <- dTrajectories  %>% group_by(coolingFunc, method, i) %>%
  summarize(lambda = median(lambda), beta = median(beta), tau = median(tau), reward = mean(reward)) %>%
  mutate(method=factor(method, levels=c('SA', 'SGD'), labels=c('Simulated Annealing (SA)', 'Stochastic Hill\nClimbing (SHC)'))) %>%
  mutate(coolingFunc=factor(coolingFunc, levels=c('fastCooling', 'expCooling', 'linearCooling'), labels=c('fast', 'exponential', 'linear')))

adultPerformance <- behavior %>% filter(agegroup=='25-55') %>% group_by(id) %>% summarize(reward = mean(z)/50) %>% pull(reward)
adultRange <- t.test(adultPerformance, conf.level = 0.95)$conf.int

# plot reward
pTraj = ggplot(mTrajectories, aes(x=i, y=reward, color=coolingFunc)) +
  # geom_path(aes(group = id), alpha = 0.1) + #replace with 
  geom_ribbon(aes(ymin = adultRange[1], ymax = adultRange[2]), fill = viridis_pal()(7)[1], color =NA, alpha  = 0.2 )+
  #geom_hline(yintercept = adultRange[1], linetype = 'dashed', color = viridis_pal()(7)[1] )+
  #geom_hline(yintercept = adultRange[2], linetype = 'dashed', color = viridis_pal()(7)[1] )+
  geom_path() +
  #annotate("segment", x = .15,xend = .15, y = adultRange[1], yend =adultRange[2],arrow = arrow(ends = "both", angle = 30, length = unit(.2,"cm")))+
  #annotate('text', x=.3, y= adultRange[2]- adultRange[1], label = '25-55 yr olds 95\% CI')+
  #coord_cartesian(ylim = c(.6, .95)) +
  facet_grid(~method) +
  scale_color_manual(values=c('#537fbb', "#FFC300", "#900C3F"), name='Cooling Schedule') +
  xlab("Iteration") +
  ylab("Reward") +
  ggtitle('Algorithms and Cooling Schedules') +
  theme_classic() +
  theme(legend.position='bottom', strip.background=element_blank(), legend.direction='horizontal',
        legend.margin=margin(-5,-5,-5,-5),
        legend.box.margin=margin(0,0,0,0))
        # legend.margin=margin(t=0, unit='cm'), legend.box.margin=margin(t=0, unit='cm'))
pTraj


# trajectories in parameter space
# faceted by tau
# select facet by selecting tau that is closest to the group median
taus = unique(sim$tau)
closestTau = which.min(abs(taus - mean(mTrajectories$tau)))
sim$tau[closestTau]

pTrajParam = ggplot() +
  geom_raster(subset(sim, tau==taus[closestTau]), mapping=aes(x=lambda, y=beta, fill=meanReward)) +
  geom_path(mTrajectories, mapping=aes(x=lambda, y=beta), color='#537fbb', size=.5) +
  scale_x_continuous(trans="log10", expand=c(0,0)) +
  scale_y_continuous(trans="log10", expand=c(0,0)) +
  xlab(expression(paste('Generalization ', lambda, ' [logscale]'))) +
  ylab(expression(paste('Exploration ', beta, ' [logscale]'))) +
  scale_fill_viridis(option='inferno', name='Reward') +
  ggtitle('Optimization Trajectories') +
  facet_grid(method~coolingFunc) +
  theme_classic() +
  theme(legend.position='right', axis.line=element_line(colour='white'),
        strip.background=element_blank(), panel.grid.minor=element_line())

pTrajParam

# human trajectory
dHuman = ddply(params, ~agegroup, plyr::summarize, lambda = median(lambda),
               beta = median(beta), tau = median(tau), mAge = mean(age_years))

pTrajAllAlgos <- pTrajParam +  geom_line(data = dHuman, aes(x=lambda, y=beta), color = 'black', alpha = 1, size = 1) +
  geom_point(data = dHuman, aes(x=lambda, y=beta, color = mAge), alpha = 0.7, size=2) +
  geom_text_repel(data = dHuman, aes(x=lambda, y=beta, label = agegroup, color = mAge), size = 5,
                  alpha = 0.9, segment.alpha = 0.5, min.segment.length = 0, arrow = arrow(length = unit(0.015, "npc")), 
                  seed = 42, box.padding = 1,  max.overlaps = Inf, nudge_x = -.2)+
  scale_color_viridis(direction = -1, name='Age') 

pTrajAllAlgos


################################################################################################
# Simplified trajectory plot 
################################################################################################
closestTau = which.min(abs(taus - median(dHuman$tau))) #For aggregate plot use median across all age groups

#SA and SHC
dModelSASHC <- trajectories %>% 
  subset(Algorithm %in% c('SGD-fastCooling', "SA-fastCooling") ) %>%
  ddply(~id+coolingFunc+method+i+Algorithm, plyr::summarize,
        lambda=mean(lambda), beta=mean(beta), tau=mean(tau), reward=mean(reward)) %>%
  ddply(~coolingFunc+method+i+Algorithm, plyr::summarize,
        lambda=median(lambda), beta=median(beta), tau=median(tau), Reward=mean(reward))

dModelSASHC$Algorithm <- factor(dModelSASHC$Algorithm, labels = c('SA fast cooling', 'SHC fast cooling'))

pDualTraj <- ggplot() +
  facet_grid(~Algorithm)+
  geom_raster(subset(sim, tau==taus[closestTau]), mapping=aes(x=lambda, y=beta, fill=meanReward)) +
  #geom_path(data = dModelPerSim, mapping = aes(x = lambda, y = beta), color='#02a6d4', size=.1, alpha = 0.1)+ #each individual trajectory
  geom_line(data = dHuman, aes(x=lambda, y=beta), color = 'black', alpha = 1, size = 1) +
  geom_point(data = dHuman, aes(x=lambda, y=beta, color = mAge), alpha = 0.7, size=2) +
  geom_text_repel(data = dHuman, aes(x=lambda, y=beta, label = agegroup, color = mAge), size = 6,
                  alpha = 0.9, segment.alpha = 0.5, min.segment.length = 0, arrow = arrow(length = unit(0.015, "npc")), 
                  seed = 42, box.padding = 1,  max.overlaps = Inf, nudge_x = -.2)+
  geom_path(dModelSASHC, mapping=aes(x=lambda, y=beta), color='#537fbb', size=.6) +
  scale_x_continuous(trans="log10", expand=c(0,0)) +
  scale_y_continuous(trans="log10", expand=c(0,0)) +
  scale_alpha_manual(name=NULL, values=c(1,1), breaks=c('Current\nAge Group', 'Algorithm')) +
  scale_fill_viridis(option='inferno', name='Reward') +
  scale_color_viridis(direction = -1, name='Age', guide = 'none') +
  xlab(expression(paste('Generalization ', lambda, ' [logscale]'))) +
  ylab(expression(paste('Exploration ', beta, ' [logscale]'))) +
  #annotate(geom='text', x=.1, y=1.5, label=bquote(~ bar(tau) == .(signif(taus[closestTau],digits=2))),color='white', size=2.7, hjust=0, vjust=0) +
  ggtitle("Human       vs. Algorithm") +
  theme_classic() +
  # guides(fill = guide_legend(override.aes = list(size = 2))) +
  theme( axis.line=element_line(colour='white'), strip.background=element_blank(),
         plot.margin=unit(c(0.1,0.1,0.1,0.1), 'cm'), legend.margin=margin(t=0,unit="cm"),
         legend.key.size = unit(0.4, 'cm'))
pDualTraj


################################################################################################
# Comparison human and algorithm parameters
################################################################################################
colors_human = c("#E69F00", "#009E73", "#D1495B")
colors_algo = c('#537fbb', '#537fbb', '#537fbb')

paramsLong = gather(params, variable, value, c('lambda', 'beta', 'tau'))

mDFPlot <- ddply(paramsLong, ~id+variable+age_months+age_years+experiment+agegroup, summarize,
                 d_ymin = max(min(value), quantile(value, 0.25) - 1.5 * IQR(value)),
                 d_ymax = min(max(value), quantile(value, 0.75) + 1.5 * IQR(value)),
                 d_lower = quantile(value, 0.25), 
                 d_middle = median(value),
                 d_upper = quantile(value, 0.75),
                 mu = mean(value))


mDFPlot$variable = factor(mDFPlot$variable, levels=c('lambda', 'beta', 'tau'),
                          ordered=TRUE, labels=c(expression(paste('Generalization ', lambda)),
                                                 expression(paste('Exploration ', beta)),
                                                 expression(paste('Temperature ', tau))))

dModelShort = dModel %>%
  subset(Algorithm=='SGD-fastCooling' & i<500) %>%
  mutate(age = seq(5,55, len=500)) %>%
  gather(variable, value, c('mLambda', 'mBeta', 'mTau'))

dModelShort$variable = factor(dModelShort$variable, 
                              levels=c('mLambda', 'mBeta', 'mTau'), ordered=TRUE,
                              labels=c(expression(paste('Generalization ', lambda)),
                                       expression(paste('Exploration ', beta)),
                                       expression(paste('Temperature ', tau))))
paramLabels = data.frame('variable' = levels(mDFPlot$variable))

p_params = ggplot(mDFPlot, aes(x=age_months/12, y=mu, color=variable, fill=variable, alpha='Human')) +
  geom_point(mapping=aes(alpha='Human'), alpha = 0.5, size=0.5) +
  geom_smooth(mapping=aes(alpha='Human'), size=0.7) +
  xlab("Human: Age (Years) [logscale]") +
  ylab("Estimate [logscale]")+
  scale_x_continuous(trans="log10", sec.axis = sec_axis(~ (. *10)-5, name = "Algorithm: iterations [logscale]")) +
  scale_y_continuous(trans="log10", breaks = c(.01, .1, 1, 10, 100), labels = c(.01, .1, 1, 10, 100)) +
  annotation_logticks(sides = 'tb')+
  ggtitle('Parameter Comparison') +
  scale_fill_manual(values = paramPal, guide='none') +
  scale_color_manual(values = paramPal, guide='none') +
  #annotate("text", label=rep(unique(mDFPlot$variable),3), x = 15, y = 10, parse=T,  size = 4 ,color = 'black', fill=NA, alpha = 1)+
  new_scale_color() +
  new_scale_fill() +
  geom_smooth(dModelShort, mapping=aes(x=age, y=value, color=variable, alpha='Algorithm'), se=FALSE) +
  theme_classic() +
  theme(strip.background=element_blank(), legend.position=c(0.93, 0.8),
        legend.spacing.y=unit(.1, 'cm'),  strip.placement = "outside") +
  facet_grid(~variable, labeller=label_parsed) +
  scale_color_manual(values = colors_algo, guide='none') +
  scale_alpha_manual(name=NULL, values=c(0.4,0), breaks=c('Human', 'Algorithm'),
                     guide=guide_legend(override.aes=list(linetype=c(1,1), shape=c(16, NA), color=c('black', '#537fbb'))))

p_params


dModelPerSim <- trajectories %>% 
  ddply(~id+coolingFunc+method+i+Algorithm, plyr::summarize,
        lambda=mean(lambda), beta=mean(beta), tau=mean(tau), reward=mean(reward))

dModel <- dModelPerSim %>%
  ddply(~coolingFunc+method+i+Algorithm, plyr::summarize,
        lambda=median(lambda), beta=median(beta), tau=median(tau), Reward=mean(reward))

allAlgos = dModel %>%
  group_by(method, coolingFunc) %>%
  mutate(age = seq(5,55, len=1500)) %>%
  gather(variable, value, c('lambda', 'beta', 'tau'))

allAlgos$variable <- factor(allAlgos$variable, levels = c('lambda', 'beta', 'tau'))
allAlgos$method <- factor(allAlgos$method, levels = c('SA', 'SGD'), labels = c('SA', 'SHC'))
allAlgos$coolingFunc <- factor(allAlgos$coolingFunc, levels = c('fastCooling', 'expCooling', 'linearCooling'), labels = c('fast', 'exponential', 'linear'))
levels(mDFPlot$variable) <- c('lambda', 'beta', 'tau')

p_paramsAll <- ggplot(allAlgos) +
  geom_point(data = mDFPlot, mapping= aes(x=age_months/12, y=mu,alpha='Human'), alpha = 0.2, size=0.5) +
  geom_smooth(data = mDFPlot, mapping= aes(x=age_months/12, y=mu,alpha='Human'), size=0.7, color = 'black') +
  xlab("Human: Age (Years) [logscale]") +
  ylab("Estimate [logscale]")+
  scale_x_continuous(trans="log10", sec.axis = sec_axis(~ (. *30)-5, name = "Algorithm: iterations [logscale]")) +
  scale_y_continuous(trans="log10", breaks = c(.01, .1, 1, 10, 100), labels = c(.01, .1, 1, 10, 100)) +
  #annotation_logticks(sides = 'tb')+
  #ggtitle('Parameter Comparison') +
  #scale_fill_manual(values = paramPal, guide='none') +
  #scale_color_manual(values = paramPal, guide='none') +
  #annotate("text", label=rep(unique(mDFPlot$variable),3), x = 15, y = 10, parse=T,  size = 4 ,color = 'black', fill=NA, alpha = 1)+
  new_scale_color() +
  new_scale_fill() +
  geom_smooth(allAlgos, mapping=aes(x=age, y=value, color=coolingFunc), alpha = 0.7, se=FALSE) +
  theme_classic() +
  theme(strip.background=element_blank(), legend.position='right',strip.placement = "outside") +
  facet_grid(method~variable, labeller=label_parsed) +
  scale_color_manual(values = c('#537fbb', "#FFC300", "#900C3F"), name = 'Cooling Schedule') +
  scale_alpha_manual(name=NULL, values=c(0.4), breaks=c('Human'),
                     guide=guide_legend(override.aes=list(linetype=c(1), shape=c(16), color=c('black'))))

p_paramsAll




########################
#Save plot
########################
#Main plot
combined = cowplot::plot_grid(pTraj, pDualTraj, rel_widths = c(.6, 1), labels='auto')
ggsave(filename = "plots/Fig4.pdf", plot = combined, height=4, width=11, units = "in")



#SI plot
pHillClimbing = cowplot::plot_grid(pTrajAllAlgos, p_paramsAll, rel_heights = c(1,1), labels='auto', ncol=1)
ggsave(filename='plots/S6.pdf', plot=pHillClimbing, height=10, width=11, units = "in")


