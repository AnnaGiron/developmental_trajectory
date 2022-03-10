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
  mutate(experiment=factor(experiment,
                           levels=c('Meder (2021)', 'Schulz (2019)', 'Adolescent'),
                           labels=c('Meder et al. (2021)', 'Schulz et al. (2019)', 'Adolescent Data')))

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
                           labels=c('Meder et al. (2021)', 'Schulz et al. (2019)', 'Adolescent Data'))) %>%
  mutate(agegroup=factor(agegroup,
                         levels=c('[25,55]', '[18,25)', '[14,18)', '[11,14)', '[9,11)', '[7,9)', '[5,7)'),
                         labels=c('25-55', '18-24', '14-17', '11-13', '9-10', '7-8', '5-6')))

behavior$z = behavior$z / 50

# append mean reward to model results
meanRew = ddply(behavior, ~id, plyr::summarize, meanReward=mean(z))
params = merge(params, meanRew, by='id')

params = subset(params, acq=='UCB')


#############################################################################################################################
# Behavior
#############################################################################################################################
# Performance by Age - separated by experiment
dAge <- ddply(subset(behavior, trial!=0), ~id+age_months+experiment, plyr::summarize, meanReward = mean(z))

pPerfAge <- ggplot(dAge, aes(x=age_months/12, y = meanReward, color=experiment, fill=experiment))+
  geom_point(aes(group=id), alpha=0.3, size=.5)+
  geom_smooth(size=.5)+
  geom_hline(yintercept=.5, linetype='dashed', color='red') + # random choice model
  scale_x_continuous(trans="log10") +
  xlab('Age (Years) [logscale]')+
  ylab('Normalized Mean Reward ± 95% CI')+
  ggtitle('a Performance as a Function of Age') +
  scale_fill_viridis(discrete=TRUE, direction=-1, name='Experiment') +
  scale_color_viridis(discrete=TRUE, direction=-1, name='Experiment') +
  theme_classic() +
  theme(text = element_text(size=11,  family="sans"), strip.background=element_blank(),
        legend.position='none', plot.title = element_text(face="bold"))

pPerfAge


#############################################################################################################################
# Model Results
#############################################################################################################################
# R2
pR2 = ggplot(subset(params, kernel=='GP'), mapping=aes(x=age_months/12, y=R2, color=experiment, fill=experiment)) +
  geom_point(alpha = 0.5, size=.7) +
  geom_smooth(size=.7) +
  xlab("Age (Years) [logscale]") +
  ylab(expression(R^2)) +
  guides(color=guide_legend("Data Set"), fill=guide_legend('Data Set')) +
  ggtitle('b Predictive Accuracy: GP-UCB') +
  scale_x_continuous(trans="log10") +
  scale_fill_viridis(discrete=TRUE, direction=-1, name='Experiment') +
  scale_color_viridis(discrete=TRUE, direction=-1, name='Experiment') +
  theme_classic() +
  theme(text = element_text(size=11,  family="sans"), strip.background=element_blank(),
        plot.title = element_text(face="bold")) 

pR2

# GP-UCB parameter estimates
paramsLong = gather(subset(params, kernel=='GP'), variable, value, c('lambda', 'beta', 'tau'))

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

pParams = ggplot(mDFPlot, aes(x=age_months/12, y=mu, color=experiment, fill=experiment)) +
  geom_point(alpha = 0.5, size=0.5)+
  geom_smooth(size=0.7)+
  xlab("Age (Years) [logscale]") +
  ylab("Estimate [logscale]")+
  scale_x_continuous(trans="log10") +
  scale_y_continuous(trans="log10", breaks = c(.01, .1, 1, 10, 100), labels = c(.01, .1, 1, 10, 100)) +
  ggtitle("c Parameter Estimates: GP-UCB") +
  scale_fill_viridis(discrete=TRUE, direction=-1, name='Experiment') +
  scale_color_viridis(discrete=TRUE, direction=-1, name='Experiment') +
  theme_classic() +
  theme(text = element_text(size=11,  family="sans"), strip.background=element_blank(),
        legend.position='none', plot.title = element_text(face="bold")) +
  facet_grid(~variable, labeller=label_parsed)

pParams

pLegend = get_legend(pR2)
pR2 = pR2 + theme(legend.position='none')

pExp = cowplot::plot_grid(cowplot::plot_grid(pPerfAge, pR2, pLegend, ncol=3, rel_widths=c(3,3,1)), pParams, ncol=1)
ggsave(filename = "plots/S1.pdf", plot = pExp, height=6, width=10, units = "in")

