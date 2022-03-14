
#params = read.csv('Data/modelFit.csv')
#params = subset(params, kernel=='RBF' & acq=='UCB')


uid = as.integer(commandArgs(TRUE)[1]) 
paramCL= as.integer(commandArgs(TRUE)[2]) 
#d1<-subset(data, id==uid)

#Source dependencies
source('models.R')
source('dataProcessing.R')
source('fit_parrallel_cluster.R')
pacman::p_load('plyr', 'dplyr', 'jsonlite', 'lsr', 'BayesFactor', 'matrixcalc','doParallel','DEoptim')
#invisible(lapply(packages, install.packages, character.only = TRUE))
#invisible(lapply(packages, require, character.only = TRUE)) #loads packages
nParticipants <- 151
# 150 participants, but id 35 is not available
##############################################################################################################
#Cluster configuration: (1 subject x model x 20simpars) per job
##############################################################################################################

#create list of all kernel functions
kernellist<-list(rbf, bayesianMeanTracker)

#names of all kernel functions}

kernelnames<-c("RBF", "BMT")

#list of all acquisition functions
acqlist<-list(greedyMean, greedyVar, ucb) 

#names of all acquisition functions
acqnames<-c("GM", "GV", 'UCB')

#all combinations of kernels and acquisition functions will be needed
#combs<-expand.grid(1:length(kernellist), 1:length(acqlist))

#create a matrix with combinations of subjectIds and model combinations
#subjectComb <- expand.grid(1:nParticipants, 1:(length(kernellist) * length(acqlist))) #1:? defines the number of unique models to be analyzed
#subjectComb = subset(subjectComb, Var1!=35) # remove id 35 as it is not available

#Cluster id from qsub
#clusterid <- sample(1:nrow(subjectComb),1) #sample random cluster id for testing
#clusterid <- 1#as.integer(commandArgs(TRUE)[1]) #Cluster id, corresponds to an integer used to indicate which combination of kernel and acquisition function to simulate

#subjectId <- subjectComb[clusterid,1] #used to identify unique subjects

set.seed(12345) #set seed as the clusterid

##############################################################################################################
#Compile Experimental Data
##############################################################################################################
#only keep people who have completed the task
#data <- dataImport_Adolescent() #sourced from dataProcessing.R
# import preprocessed data
data = read.csv('Data/AdolescentGrid.csv')
#Normalize data
data$z <- (data$z - 25) / 50

#uid <- unique(data$id)[subjectId] #convert subjectId to uid

#extract environments for simulation
environments <- lapply(fromJSON("Data/smoothKernel.json"), FUN=function(x) matrix(as.numeric(unlist(x)), ncol=3, byrow=TRUE, dimnames=list(seq(1,64),  c('x2', 'y', 'x1'))))
env <- as.data.frame(environments[[1]])
for (i in 2:40){
  env<-rbind(env,as.data.frame(environments[[i]]))
}
env$en<-rep(1:40, each=64)

#Scaling reward range
env$y<-env$y*50


fittedModel = read.csv('Data/modelFit.csv')
fittedModel = fittedModel %>%
  mutate(kernel=factor(kernel, levels=c('RBF', 'BMT'), labels=c('GP', 'BMT'))) %>%
  mutate(agegroup=factor(agegroup,
                         levels=rev(c('[25,55]', '[18,25)', '[14,18)', '[11,14)', '[9,11)', '[7,9)', '[5,7)')),
                         labels=rev(c('25-55', '18-24', '14-17', '11-13', '9-10', '7-8', '5-6')))) %>%
  mutate(experiment=factor(experiment, levels=c('Meder (2021)', 'Schulz (2019)', 'Adolescents')))%>%filter(experiment=='Adolescents')

# only use GP-UCB model
fittedModel = subset(fittedModel, kernel=='GP' & acq=='UCB')

#############################################################################################################################
# Simulating performance for different parameter values

#############################################################################################################################
# Tukey's fence to compute upper and lower bound for each parameter
tukeysFence <- function(x, k = 1.5, na.rm = TRUE) {
  quar <- quantile(x, probs = c(0.25, 0.75), na.rm = na.rm)
  iqr <- diff(quar)
  
  lower = quar[1] - k * iqr
  upper = quar[2] + k * iqr
  
  return(c(lower, upper))
}


# parameter range based on Tukey's fence
lambda = tukeysFence(log(fittedModel$lambda))
beta = tukeysFence(log(fittedModel$beta))
tau = tukeysFence(log(fittedModel$tau))

#lambdaQuant=quantile(log(fittedModel$lambda),probs=c(0.05,0.95))
#betaQuant=quantile(log(fittedModel$beta),probs=c(0.05,0.95))
#tauQuant=quantile(log(fittedModel$tau),probs=c(0.05,0.95))

# parameters to simulate take the exp out. Done later.
paramsSims = tibble(lambda=seq(lambda[1], lambda[2], len=20), 
                         beta=seq(beta[1], beta[2], len=20),
                         tau=seq(-5, tau[2], len=20))


# Get envirionments for specific subject
myjson <- fromJSON("Data/BanditData.json")
#env<-myjson$records$data$envOrder#[[k]]

for (i in 1:nrow(myjson)){
  subd <- myjson[i,]
  new<-fromJSON(subd$data)
  if(new$uid==uid){
    break
  }
}

environmentIndex=new$envOrder+1

#library(doParallel)

#####setup parralel
cl <- makeCluster(8)
registerDoParallel(cl)
####

#
allrecov=NULL
for (i in 1:length(paramsSims$lambda)){
  if (paramCL==1){
    print("Varying Lambda")
    simP<-c(paramsSims$lambda[i],
            fittedModel%>%filter(id==uid)%>%pull(beta)%>%log(),#%>%as_vector()
            fittedModel%>%filter(id==uid)%>%pull(tau)%>%log())
  }else if (paramCL==2){
    print("Varying beta")
    simP<-c(fittedModel%>%filter(id==uid)%>%pull(lambda)%>%log(),
            paramsSims$beta[i],#%>%as_vector()
            fittedModel%>%filter(id==uid)%>%pull(tau)%>%log())#%>%as_vector()
  }else if (paramCL==3){
    print("Varying tau")
    simP<-c(fittedModel%>%filter(id==uid)%>%pull(lambda)%>%log(),
            fittedModel%>%filter(id==uid)%>%pull(beta)%>%log(),#%>%as_vector()
            paramsSims$tau[i])#%>%as_vector()
  }else if (paramCL==4){
    print("fit Subs")
    simP<-c(fittedModel%>%filter(id==uid)%>%pull(lambda)%>%log(),
            fittedModel%>%filter(id==uid)%>%pull(beta)%>%log(),#%>%as_vector()
            fittedModel%>%filter(id==uid)%>%pull(tau)%>%log())
  }

  print("simulate GP")
  
  sim<-modelSimulateChoice(par=simP,rounds=2:9, k=rbf, acquisition=ucb,env=env,envSelector=environmentIndex)
  #add simulations to the dataconstruct
  d2=tibble(
    z=(sim$mu -25)/50,
    y=sim$x1,
    x=sim$x2,
    chosen2=sim$chosen,
    round=sim$round
  )
  print(simP)
  print("Setup_CV")

  # i learned this from looking at erics code. just reproduce xy coords. sample isnt right maybe
  Xnew<-expand.grid(x=0:7,y=0:7) #do this outside the loop for better speed
  for (k in 1:nrow(d2)){
    d2$chosen[k]<-which(Xnew$x==d2$x[k] & Xnew$y == d2$y[k])#### x1 and x2 in envirionment as compared to xNew?
  }
  
  
  #d<-read_csv(file = "./paramterRecovery/paramRecovery_1.csv")
  
  #ah<-d%>%group_by(X1,round)%>%dplyr::mutate(id=group_indices())    
  
  #ah# the number of rounds differ across the experiment
  
  #d2<-d%>%filter(id==1)%>%magrittr::set_colnames(c("rowNr","id","round","lambda","beta","tau","z","x1","x2"))
  ### test: where is the fucking error
 # fit<-DEoptim(modelFit, lower=lbound, upper=ubound, subjD=d2, k=rbf, rounds = 2:9, acquisition=ucb, DEoptim.control(itermax=100,trace = 1))
  
  out<-foreach(leaveoutindex = 2:9,.combine = 'rbind')%dopar% {# do the cv on multiple cores
    #source of modeling code
    source("models.R")
    #packages
    packages <- c('plyr', 'dplyr', 'jsonlite', 'lsr', 'BayesFactor', 'matrixcalc')
    invisible(lapply(packages, require, character.only = TRUE)) #loads packages for parallel computing
    
    Xnew<-as.matrix(expand.grid(0:7,0:7)) #do this outside the loop for better speed
    FitParamCV(d1=d2,rbf, ucb,leaveoutindex)
  }
  print("Done")
  out%<>%as_tibble()%>%
    mutate(lambdaSim=simP[1],
         betaSim=simP[2],
         tauSim=simP[3],
         uid=uid)
  allrecov<-rbind(allrecov,out)
}
saveRDS(d2,paste0("./paramterRecovery/Sub_",uid,"_",paramCL,".rds"))
saveRDS(allrecov,paste0("./paramterRecovery/Run_",uid,"_New_",paramCL,".rds"))







# 
# 
# params%>%filter(experiment=='Adolescents')%>%select(lambda,beta,tau,id,age_years)->ah
# 
# 
# 
# 
# ah%>%pivot_longer(c(lambda,beta,tau))%>%
#   mutate(value=(value))%>%
#   ggplot(aes(x=value,fill=age_years))+
#   geom_histogram()+  
#   facet_wrap(.~name)
# 
# ah###
# ##Sth else here:
# ###
# ##
# #
# 
# one<-env%>%
#   ggplot(aes(x=x1,y=x2,fill=y))+
#   geom_tile()+
#   scale_fill_distiller(palette="Reds")+
#   facet_wrap(.~en)+
#   ggtitle("envs x1 x2")
# ggsave(plot = one,filename = "./X_Figures/allenvs_xy.png",width=20,height=20)
# 
# two<-env%>%
#   ggplot(aes(x=x2,y=x1,fill=y))+
#   geom_tile()+
#   scale_fill_distiller(palette="Reds")+
#   facet_wrap(.~en)+
#   ggtitle("envs x2 x1")
# ggsave(plot = two,filename = "./X_Figures/allenvs_yx.png",width=20,height=20)
