# Compare algorithms to human trajectories
# Simon Ciranka and Charley Wu, 2022
rm(list=ls())

packages <- c('plyr', 'dplyr', 'tidyr', 'ggplot2', 'cowplot', 'viridis', 'entropy', 'sjPlot', 'tidyverse', 'posterior', 
              'Rmisc',"doParallel", "fmsb", "furrr", "zoo", 'forecats', 'purrr', 'future', 'brms', 'tibble', 'ggstance')
invisible(lapply(packages, require, character.only = TRUE)) #loads packages

paramPal = c("#FFEE67", '#27AD88', "#D1495B") #parameter colors


###### Load human data
params = read.csv('data/modelFit.csv')

params = params %>%
  mutate(kernel=factor(kernel, levels=c('RBF', 'BMT'), labels=c('GP', 'BMT'))) %>%
  mutate(agegroup=factor(agegroup,
                         levels=rev(c('[25,55]', '[18,25)', '[14,18)', '[11,14)', '[9,11)', '[7,9)', '[5,7)')),
                         labels=rev(c('25-55', '18-24', '14-17', '11-13', '9-10', '7-8', '5-6')))) %>%
  mutate(experiment=factor(experiment, levels=c('Meder (2021)', 'Schulz (2019)', 'Adolescents')))
params$ModelName = paste(params$kernel, params$acq, sep="-")
params$ModelName = factor(params$ModelName)

params = subset(params, kernel=='GP' & acq=='UCB') 


###### Load algorithm data

# read algorithm trajectory
path = '../hillClimbingAlgorithm/batch2/'
filenames = list.files(path=path, pattern='*.csv')
filenames = paste0(path, filenames)
algotraj = ldply(filenames, read.csv)
algotraj$Algorithm = paste(algotraj$method, algotraj$coolingFunc, sep='-')
algotraj$trajectory = rep(1:720, each=1500)


algotraj_Fast<-algotraj%>%filter(coolingFunc=="fastCooling" & method=="SGD")
#subset algo at human rewards
algotraj_Fast%>%#group_by(id)%>%
  mutate(individual_traj=paste0(id,"_",leaveoutindex))%>%
  ggplot(aes(x=maxRewardDiff))+
  #geom_point(aes(color=as.factor(aaah)))+
  geom_histogram()+
  facet_wrap(.~individual_traj)

#subset each by finding the first time the second derivative or reward is 0, 100 times in a row. 
subetIDX_df<-algotraj_Fast%>%#group_by(id)%>%
  mutate(individual_traj=paste0(id,"_",leaveoutindex))%>%
  group_by(individual_traj)%>%
  mutate(der_2=c(NA,abs(diff(maxRewardDiff))))%>%
  mutate(rollavg_der_2=zoo::rollmean(der_2,k=20,na.pad = T))%>%# rolling average
  filter(rollavg_der_2 == min(rollavg_der_2,na.rm = T))%>% 
  slice(1)%>%mutate(i_subset=i)%>%select(i_subset,individual_traj)

#join with old algo trajectory to be able to filter 
algotraj_Fast<-algotraj_Fast%>%#group_by(id)%>%
  mutate(individual_traj=paste0(id,"_",leaveoutindex))

algotraj_until_converge<-left_join(subetIDX_df,algotraj_Fast,by="individual_traj")%>%
  group_by(individual_traj)%>%
  filter(i<i_subset)#%>%

# make symmetric trajectories based on sampling.

# rank normalize humans  
norm_rank_human<-params%>%mutate(
  norm_rank=(age_years-1)/(max(age_years)-1)
)%>%dplyr::group_by(
  norm_rank
)%>%dplyr::mutate(
  n=n()
)%>%ungroup()%>%
  mutate(norm_rank=ifelse(norm_rank==min(norm_rank),0,norm_rank))

# rank normalize algorithm THIS TAKES ALL ALGORITHM ITERATIONS, change algotraj fast to subetIDX_df for subsetted iterations
algotraj_until_converge<-
  algotraj_Fast%>%group_by(individual_traj)%>%dplyr::mutate(
    norm_rank_a=(i-1)/(max(i)-1)
  )%>%mutate(norm_rank_a=ifelse(norm_rank_a==min(norm_rank_a),0,norm_rank_a))# set lowest to 0

norm_rank_algo_matched<-algotraj_until_converge%>%dplyr::filter(
  round(norm_rank_a,3) %in% unique(round(norm_rank_human$norm_rank,3))
)

# sample from algo so that each rank has the same amount of data
# do this 1000 times and then fit your model on this 1000 times.
# then, build reference distribution

n_samples=100
resampled_dfs_list<-list()

for(i in 1:n_samples){
  # make placeholders
  norm_rank_human$sampled_algo_lambda<-NA
  norm_rank_human$sampled_algo_beta<-NA
  norm_rank_human$sampled_algo_tau<-NA
  for (rank in unique(norm_rank_human$norm_rank)){
    #define vectors to sample from
    lambda_vec=norm_rank_algo_matched[round(norm_rank_algo_matched$norm_rank_a,3)==round(rank,3),]$lambda
    beta_vec=norm_rank_algo_matched[round(norm_rank_algo_matched$norm_rank_a,3)==round(rank,3),]$beta
    tau_vec=norm_rank_algo_matched[round(norm_rank_algo_matched$norm_rank_a,3)==round(rank,3),]$tau
    #how often (define per agegroup)
    n_sample=max(norm_rank_human%>%filter(norm_rank==rank)%>%.$n)
    #allocate samples
    norm_rank_human[norm_rank_human$norm_rank==rank,]$sampled_algo_lambda=sample(lambda_vec,n_sample)
    norm_rank_human[norm_rank_human$norm_rank==rank,]$sampled_algo_beta=sample(beta_vec,n_sample)
    norm_rank_human[norm_rank_human$norm_rank==rank,]$sampled_algo_tau=sample(tau_vec,n_sample)
  }
  # concat different syntehtic
  resampled_dfs_list[[i]]<-norm_rank_human
}
cp_resampled_df<-resampled_dfs_list%>%map_dfr(., .f = cbind, .id = 'index')
#norm_rank_human%>%arrange(norm_rank)
# 
#   ggplot(aes(y=mLambda,x=norm_rank))+
#   stat_summary(aes(group=id),geom="point")+
#   stat_summary(data=norm_rankParams,aes(x=norm_rank,y=lambda),color="red")+theme_bw()



## cp model prototype!
bprior <- c(
  brms::set_prior("normal(0, 2)", resp = c("lambda","beta","tau"),nlpar = "b0"),
  # prior(exponential(0.1), nlpar = "cpunc")+
  brms::set_prior("normal(0, 1)", resp = c("lambda","beta","tau"),nlpar = "b1"),
  brms::set_prior("normal(0, 1)", resp = c("lambda","beta","tau"),nlpar = "b2"),
  brms::set_prior("normal(0, 1)", resp = c("lambda","beta","tau"),nlpar = "alpha")
)# not entirely sure what a good prior is here

# inv logit looses the step function a little
bform_resampled <-  brms::bf(
  brms::mvbind(lambda,beta,tau) ~ b0 + b1 * (norm_rank - omega) * step((omega - norm_rank)) + 
    b2 * (norm_rank - omega) * step((norm_rank - omega)),
  b0 + b1 + b2 + alpha ~ 1+(1|index),#+(1|id),
  nlf(omega ~ inv_logit(alpha)),# bc normrank is already on the 0-1 range we do not need to resample
  nl = TRUE
)
#make dummy model that you update later to avoid recompiling. 

cp_resampled_df_mod<-cp_resampled_df%>%mutate(
  lambda=log(sampled_algo_lambda),
  beta=log(sampled_algo_beta),
  tau=log(sampled_algo_tau)
)#%>%brms::brm(data=.,formula = bform,prior = bprior,chains=3)

#resampled_model<-brms::brm(data=cp_resampled_df,formula = bform_resampled,prior = bprior,chains=3,cores=3)
#saveRDS(resampled_model,file = "resampled_cp_model.rds")
resampled_model<-readRDS("resampled_cp_model.rds")


# coefficients

betas_df_algo<-resampled_model%>%summary()%>%#get credible intervals
  .$fixed%>%.[c(2,3,6,7,10,11),]%>%# fixed effects
  select(`l-95% CI`,`u-95% CI`,Estimate)%>%#subset
  rownames_to_column()%>%
  mutate(slope=case_when(
    grepl("b1",rowname,fixed=T)~"b1",
    grepl("b2",rowname,fixed=T)~"b2"
  ),
  parameter=case_when(
    grepl("lambda",rowname,fixed=T)~"1",
    grepl("beta",rowname,fixed=T)~"2",
    grepl("tau",rowname,fixed=T)~"3"
  ),
  trajectory="algorithm"
  )%>%select(-rowname)



#human coefficients
multi_Change<-readRDS("../brms_modelfits/multi_Change.Rds")

betas_df_human<-multi_Change%>%summary()%>%#get credible intervals
  .$fixed%>%.[c(2,3,6,7,10,11),]%>%# fixed effects
  select(`l-95% CI`,`u-95% CI`,Estimate)%>%#subset
  rownames_to_column()%>%
  mutate(slope=case_when(
    grepl("b1",rowname,fixed=T)~"b1",
    grepl("b2",rowname,fixed=T)~"b2"
  ),
  parameter=case_when(
    grepl("lambda",rowname,fixed=T)~"1",
    grepl("beta",rowname,fixed=T)~"2",
    grepl("tau",rowname,fixed=T)~"3"
  ),
  trajectory="human"
  )%>%select(-rowname)

betas_df<-rbind(betas_df_human,betas_df_algo)
betas_df$trajectory<-factor(betas_df$trajectory,levels=c("human","algorithm"),ordered=T)
levels(betas_df$trajectory) <- c('Human', 'SHC-Fast')

betas_plot<-betas_df%>%ggplot(aes(x=Estimate,color=parameter,y=parameter,shape=slope))+
  geom_vline(aes(xintercept=0), linetype = 'dashed')+
  #geom_errorbarh(aes(xmin=`l-95% CI`,xmax=`u-95% CI`),height =0.1,position = ggstance::position_dodgev(height = 0.5))+
  geom_pointrangeh(aes(xmin=`l-95% CI`,xmax=`u-95% CI`),position = ggstance::position_dodgev(height = 0.5))+
  #scale_alpha_discrete(position_nudge()+
  scale_shape(name = 'Coef')+
  scale_color_manual(breaks=c("1","2","3"),
                     labels=c(expression(lambda),expression(beta),expression(tau)),
                     values=c(paramPal[1],paramPal[2],paramPal[3]), name = 'Param')+
  scale_y_discrete(breaks=c("1","2","3"),
                   labels=c(expression(lambda),expression(beta),expression(tau)), name = 'Parameter',limits=rev )+
  scale_x_continuous(name="Posterior Est. ± 95% CI")+
  facet_wrap(~trajectory,scales="free_x")+
  ggtitle("Changepoint regression coefficients")+
  guides(color="none")+
  theme_classic()+
  theme(strip.background=element_blank(), legend.position= 'right')

betas_plot

####### Histogram of change point


lambda_hist_algo <- resampled_model%>%tidybayes::spread_draws(b_lambda_alpha_Intercept)%>%rowwise()%>%
  mutate(hm=inv_logit_scaled(b_lambda_alpha_Intercept,lb = 0,ub=1500))%>%
  select(-b_lambda_alpha_Intercept)


lambda_hist_algo$param <- 'lambda'

lambdaConvergence <- lambda_hist_algo%>%ungroup()%>%
  summarize(upper_CI=quantile(hm,probs = 0.95))

beta_hist_algo <- resampled_model%>%tidybayes::spread_draws(b_beta_alpha_Intercept)%>%rowwise()%>%
  mutate(hm=inv_logit_scaled(b_beta_alpha_Intercept,lb = 0,ub=1500))%>%
  select(-b_beta_alpha_Intercept)

beta_hist_algo$param <- 'beta'

betaConvergence <- beta_hist_algo%>%ungroup()%>%summarize(upper_CI=quantile(hm,probs = 0.95))

tau_hist_algo <- resampled_model%>%tidybayes::spread_draws(b_tau_alpha_Intercept)%>%rowwise()%>%
  mutate(hm=inv_logit_scaled(b_tau_alpha_Intercept,lb = 0,ub=1500))%>%
  select(-b_tau_alpha_Intercept)

tau_hist_algo$param <- 'tau'

tauConvergence <- tau_hist_algo%>%ungroup()%>%summarize(upper_CI=quantile(hm,probs = 0.95))

#Combine all algorithm omega estimates together
all_hist_algo <- rbind(lambda_hist_algo, beta_hist_algo,tau_hist_algo )
all_hist_algo$param <- factor(all_hist_algo$param, levels = c('lambda', 'beta', 'tau'))

#combine convergence together
all_Algo_convergence <- rbind(lambdaConvergence, betaConvergence, tauConvergence)
all_Algo_convergence$param <- c('lambda', 'beta', 'tau')
all_Algo_convergence$param <- factor(all_Algo_convergence$param, levels = c('lambda', 'beta', 'tau'))


change_hist_algo <- ggplot(all_hist_algo, aes(x=hm, fill = param))+
  geom_histogram(bins=30, color = 'black',  aes(y = stat(width*density)), boundary = 0.5)+
  geom_vline(data = all_Algo_convergence, aes(xintercept = upper_CI ), linetype = 'dashed')+
  scale_x_log10(name="Iteration [logscale]")+
  #scale_y_continuous(name=expression(paste(omega, '\n(posterior density)')))+
  scale_y_continuous(name=expression(omega))+
  scale_fill_manual(breaks=c("lambda","beta","tau"),
                    labels=c(expression(lambda),expression(beta),expression(tau)),
                    values=c(paramPal[1],paramPal[2],paramPal[3]))+
  facet_grid(param ~ ., labeller = label_parsed,  scales='free_y')+
  ggtitle("Changepoint")+
  theme_classic() +
  theme( strip.background=element_blank(), legend.position='none', 
         strip.text.x = element_blank())
change_hist_algo


#### convergence plots

#Human converged parameters

multiChange<-readRDS("../brms_modelfits/multi_Change.Rds") #load change point regression

##Define convergence point based on upper 95% CI
### Lambda
lambdaHist <- multiChange%>%tidybayes::spread_draws(b_lambda_alpha_Intercept)%>%rowwise()%>%
  mutate(hm=inv_logit_scaled(b_lambda_alpha_Intercept,lb = 5,ub=25))%>%
  select(-b_lambda_alpha_Intercept)

lambdaHist$param <- 'lambda'

lambdaConvergence <- lambdaHist%>%ungroup()%>%
  summarize(upper_CI=quantile(hm,probs = 0.95))

### beta
betaHist <- multiChange%>%tidybayes::spread_draws(b_beta_alpha_Intercept)%>%rowwise()%>%
  mutate(hm=inv_logit_scaled(b_beta_alpha_Intercept,lb = 5,ub=25))%>%
  select(-b_beta_alpha_Intercept)

betaHist$param <- 'beta'

betaConvergence <- betaHist%>%ungroup()%>%summarize(upper_CI=quantile(hm,probs = 0.95))

### Tau
tauHist <- multiChange%>%tidybayes::spread_draws(b_tau_alpha_Intercept)%>%rowwise()%>%
  mutate(hm=inv_logit_scaled(b_tau_alpha_Intercept,lb = 5,ub=25))%>%
  select(-b_tau_alpha_Intercept)

tauHist$param <- 'tau'

tauConvergence <- tauHist%>%ungroup()%>%summarize(upper_CI=quantile(hm,probs = 0.95))

###Combine together
allConvergence <- rbind(lambdaConvergence, betaConvergence,tauConvergence )
allConvergence$param <- c('lambda', 'beta', 'tau')

##Now filter human data
humanDF <- data.frame()
for (p in c('lambda', 'beta', 'tau')){
  convergedPs <- params[params$age_years>as.numeric(allConvergence[allConvergence$param==p, 'upper_CI']),p]
  dummyDF <- data.frame(type='Human', param = p, value = convergedPs)
  humanDF <- rbind(humanDF, dummyDF)
}

##Now filter algorithm data
algoPars <- algotraj_Fast %>% group_by(id, i) %>% summarize(lambda = mean(lambda), beta = mean(beta), tau = mean(tau))
algoDF <- data.frame()
for (p in c('lambda', 'beta', 'tau')){
  convergedPs <- algoPars[algoPars$i>as.numeric(all_Algo_convergence[all_Algo_convergence$param==p, 'upper_CI']),p]
  dummyDF <- data.frame(type='SHC-Fast', param = p, value = unlist(convergedPs))
  algoDF <- rbind(algoDF, dummyDF)
}

combinedConvergenceDF <- rbind(humanDF, algoDF)
combinedConvergenceDF$type <- factor(combinedConvergenceDF$type)
combinedConvergenceDF$param <- factor(combinedConvergenceDF$param, levels = c('lambda', 'beta', 'tau'))

#Stats
source('../statisticalTests.R')

ranktestPretty(subset(combinedConvergenceDF, type =='Human' & param == 'lambda')$value
               - median(subset(combinedConvergenceDF, type =='SHC-Fast' & param == 'lambda')$value), oneSample=T ) #$Z=-6.1$, $p<.001$, $r=-.73$, $BF>100$

ranktestPretty(subset(combinedConvergenceDF, type =='Human' & param == 'beta')$value
               - median(subset(combinedConvergenceDF, type =='SHC-Fast' & param == 'beta')$value), oneSample=T) #$Z=0.2$, $p=.599$, $r=.02$, $BF=.11$

ranktestPretty(subset(combinedConvergenceDF, type =='Human' & param == 'tau')$value
               - median(subset(combinedConvergenceDF, type =='SHC-Fast' & param == 'tau')$value), oneSample=T) #$Z=-9.1$, $p<.001$, $r=-.64$, $BF>100$



convergencePlot <- ggplot(combinedConvergenceDF, aes(x = value, fill = param, color = type))+
  geom_histogram(bins=30, color = 'black',  aes(y = stat(width*density)), boundary = 0.5)+
  #geom_vline(data = meanConvergence, aes(xintercept = value), color = 'black', linetype ='dashed')+
  facet_grid(type~param, labeller = label_parsed,  scales='free')+
  scale_x_log10(name="Parameter value [logscale]")+
  theme_classic()+
  scale_y_continuous(name='Density')+
  scale_fill_manual(values = paramPal)+
  ggtitle('Convergence differences')+
  theme(strip.background=element_blank(), legend.position='none')

convergencePlot


#########################3

top <- cowplot::plot_grid(betas_plot, change_hist_algo, labels = 'auto')
finalPlot <- cowplot::plot_grid(top, convergencePlot, ncol = 1, labels = c('','c'))
finalPlot

ggsave(filename = "../plots/S8.pdf",finalPlot,width = 10,height = 7)


#Summary of model for table
fixef(resampled_model)

#Omega values need to be transformed back into the original range of 0-1500

inv_logit_scaled(fixef(resampled_model)['lambda_alpha_Intercept',], ub=1500,lb=0)
inv_logit_scaled(fixef(resampled_model)['beta_alpha_Intercept',], ub=1500,lb=0)
inv_logit_scaled(fixef(resampled_model)['tau_alpha_Intercept',], ub=1500,lb=0)
