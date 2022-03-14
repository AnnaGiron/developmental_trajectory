
#params = read.csv('Data/modelFit.csv')
#params = subset(params, kernel=='RBF' & acq=='UCB')


uid = 1#as.integer(commandArgs(TRUE)[1]) 
simMod= 2#as.integer(commandArgs(TRUE)[2]) 
Rep= 1#as.integer(commandArgs(TRUE)[3]) 

#d1<-subset(data, id==uid)

#Source dependencies
source('../models.R')
source('../dataProcessing.R')
source('fit_parrallel_cluster.R')
pacman::p_load('plyr', 'dplyr', 'jsonlite', 'lsr', 'BayesFactor', 'matrixcalc','doParallel','DEoptim','tidyverse')
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
acqnames<-c("GM", "EG", 'UCB')


newModels=c("BMT-UCB","RBF-UCB","RBF-GM","RBF-EG")
newModelList=list()

newModelList[[1]]=list(rbf,ucb,"RBF","UCB")
newModelList[[2]]=list(rbf,greedyMean,"RBF","GM")
newModelList[[3]]=list(rbf,epsilonGreedy,"RBF","EG")
newModelList[[4]]=list(bayesianMeanTracker,ucb,"BMT","UCB")


#all combinations of kernels and acquisition functions will be needed
#combs<-expand.grid(1:length(kernellist), 1:length(acqlist))

#create a matrix with combinations of subjectIds and model combinations
#subjectComb <- expand.grid(1:nParticipants, 1:(length(kernellist) * length(acqlist))) #1:? defines the number of unique models to be analyzed
#subjectComb = subset(subjectComb, Var1!=35) # remove id 35 as it is not available

#Cluster id from qsub
#clusterid <- sample(1:nrow(subjectComb),1) #sample random cluster id for testing
#clusterid <- 1#as.integer(commandArgs(TRUE)[1]) #Cluster id, corresponds to an integer used to indicate which combination of kernel and acquisition function to simulate

#subjectId <- subjectComb[clusterid,1] #used to identify unique subjects

set.seed(Rep) #set seed as the clusterid

##############################################################################################################
#Compile Experimental Data
##############################################################################################################
#only keep people who have completed the task
#data <- dataImport_Adolescent() #sourced from dataProcessing.R
# import preprocessed data
data = read.csv('../data/AdolescentGrid.csv')
#Normalize data
data$z <- (data$z - 25) / 50

#uid <- unique(data$id)[subjectId] #convert subjectId to uid

#extract environments for simulation
environments <- lapply(fromJSON("../data/smoothKernel.json"), FUN=function(x) matrix(as.numeric(unlist(x)), ncol=3, byrow=TRUE, dimnames=list(seq(1,64),  c('x2', 'y', 'x1'))))
env <- as.data.frame(environments[[1]])
for (i in 2:40){
  env<-rbind(env,as.data.frame(environments[[i]]))
}
env$en<-rep(1:40, each=64)

#Scaling reward range
env$y<-env$y*50

fittedModel = read.csv('../data/modelFit_OriginalID.csv')
EnvFile=read.csv("../data/ykwg_data.csv")

fittedModel = fittedModel %>%
  #mutate(kernel=factor(kernel, levels=c('RBF', 'BMT'), labels=c('GP', 'BMT'))) %>%
  mutate(agegroup=factor(agegroup,
                         levels=rev(c('[25,55]', '[18,25)', '[14,18)', '[11,14)', '[9,11)', '[7,9)', '[5,7)')),
                         labels=rev(c('25-55', '18-24', '14-17', '11-13', '9-10', '7-8', '5-6')))) %>%
  mutate(experiment=factor(experiment, levels=c('Meder (2021)', 'Schulz (2019)', 'Adolescents')))%>%filter(experiment=='Meder (2021)')

# get back original IDs for Epsilon Greedy.
NormalIds=fittedModel%>%
  filter(acq=="UCB")%>%.$id%>%unique()

EGId=fittedModel%>%
  filter(acq=="EG")%>%.$id%>%unique()

# this is the smartest thing i´ve ever coded.
fittedModel[fittedModel$acq=="EG",]%<>%
  rowwise()%>%
  mutate(id=NormalIds[which(.$id==EGId,arr.ind = T)])

fittedmodelid=unique(fittedModel$id)
uidReal=fittedmodelid[uid]# differnt from input in the commandline
unique(EnvFile$id)

environmentIndex<-EnvFile%>%filter(id==uidReal)%>%pull(env)%>%unique()
environmentIndex=environmentIndex+1# try this, for the adolescent data, the json startet at 0 but env file at 1 you might have to skip it if recovery breaks
fittedModel = subset(fittedModel, kernel==unlist(newModelList[[simMod]][3]) & acq==unlist(newModelList[[simMod]][4]))

#############################################################################################################################
# Simulating performance for different parameter values
#############################################################################################################################
# Tukey's fence to compute upper and lower bound for each parameter



#=new$envOrder+1

#library(doParallel)

#####setup parralel
cl <- makeCluster(length(environmentIndex))
registerDoParallel(cl)
####
Xnew<-as.matrix(expand.grid(0:7,0:7))

logitTransform<-function(x){
  y = log(x/(1 - x))
  return(y)
}

simP=c(fittedModel%>%filter(id==uidReal)%>%pull(lambda)%>%log(),
       fittedModel%>%filter(id==uidReal)%>%pull(beta)%>%log(),
       fittedModel%>%filter(id==uidReal)%>%pull(tau)%>%log(),
       fittedModel%>%filter(id==uidReal)%>%pull(epsilon)%>%logitTransform())# other transformation

simP=simP[!is.na(simP)]
#
allrecov=NULL
print("simulate GP")
acqlist

sim<-modelSimulateChoice(par=simP,rounds=1:length(environmentIndex), k=unlist(newModelList[[simMod]][1])[[1]], acquisition=unlist(newModelList[[simMod]][2])[[1]],env=env,envSelector=environmentIndex)
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

print(d2)

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
print("youre already here!")
for (z in 1:length(newModelList)){
  #if(acqnames[z]=="EG")
  out<-foreach(leaveoutindex = 1:length(environmentIndex),.combine = 'rbind')%dopar% {# do the cv on multiple cores
    #source of modeling code
    source("../models.R")
    #packages
    packages <- c('plyr', 'dplyr', 'jsonlite', 'lsr', 'BayesFactor', 'matrixcalc')
    invisible(lapply(packages, require, character.only = TRUE)) #loads packages for parallel computing
    
    Xnew<-as.matrix(expand.grid(0:7,0:7)) #do this outside the loop for better speed
    FitParamCV(d1=d2,unlist(newModelList[[z]][1])[[1]], acquisition=unlist(newModelList[[z]][2])[[1]],leaveoutindex)
  }
  print("Done")
  out%<>%as_tibble()%>%pivot_longer(-c(V1,V2))%>%
    mutate(kernel=newModelList[[z]][3][[1]],
           aq=newModelList[[z]][4][[1]],
           uid=uid)
  allrecov<-rbind(allrecov,out)
}



saveRDS(d2,paste0("./modelRecovery/Meder_Sub_",uid,"_",newModels[simMod],"_",Rep,"_data.rds"))
saveRDS(allrecov,paste0("./modelRecovery/Meder_Sub_",uid,"_",newModels[simMod],"_",Rep,"_data.rds"))







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
