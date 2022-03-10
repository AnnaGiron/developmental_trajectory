#data cleaning
#Charley Wu, Anna Giron, 2021

library('jsonlite')
library('dplyr')
library('plyr')
library('dplyr')

#############################################################################################################################
# Behavioral Data
#############################################################################################################################
#############################################################################################################################
# Adolescent data
dataImport_Adolescent <- function(dataFile ="data/BanditData.json", demographics="data/Demographics.csv", writecsv = FALSE){
  #read in json
  myjson <- fromJSON(dataFile)
  all_opts = expand.grid(0:7, 0:7)
  dat<-data.frame()
  
  envs = data.frame()
  
  # import demographics.csv
  # replace invalid ages in BanditData with more reasonable values from Demographics
  # for IDs 51, 2, 6, 57
  demo = read.csv(demographics, header=T, sep=',')
  demo$sex = as.character(demo$sex)
  ids_age = c(51, 2, 6, 57)
  
  #loop through json to create data frame
  for (i in 1:nrow(myjson)){
    subd <- myjson[i,]
    subdata <- fromJSON(subd$data)
    #x-y-z
    x<-as.vector(t(subdata$searchHistory$xcollect))
    y<-as.vector(t(subdata$searchHistory$ycollect))
    z<-as.vector(t(subdata$searchHistory$zcollect))
    chosen <- apply(cbind(x,y),MARGIN=1, FUN=function(row) which(row[1]==all_opts$Var1 & row[2]==all_opts$Var2))
    zscaled <- as.vector(t(subdata$searchHistory$zcollectScaled))
    #Color value
    #time
    time<- as.vector(t(subdata$searchHistory$tscollect))
    #condition; no conditions
    #cond<-rep(subd$condition, length(x))
    #id
    id<-rep(as.numeric(subd$uid), length(x))
    #age
    age_years<-rep(as.numeric(subdata$age), length(x))
    #demo_age = demo[demo$subject==id[1],]$age
    #if(length(demo_age)!=0) {
    #  if(age_years[1]!=demo_age) {
    #    print(sprintf('ID: %s, BanditData: %s, Demographics: %s', id[1], age_years[1], demo_age))
    #  }
    #}
    if(id[1] %in% ids_age) {
      age_years = rep(demo[demo$subject==id[1],]$age, length(x))
    }
    #age_months <- ageMonths[ageMonths$id==subd$participant$code, 'age_months']
    #trial number
    trial<-rep(0:25, 10)
    #round number
    round<-rep(1:10, each=26)
    #gender
    gender <- subdata$gender
    gender = rep(gender, length(x))
    #demo_gender = demo[demo$subject==id[1],]$sex
    #if(length(demo_gender)!=0) {
    #  demo_gender = if (demo_gender == 'Maennlich') 'Male' else 'Female'
    #  if(gender[1]!=demo_gender) {
    #    print(sprintf('ID: %s, BanditData: %s, Demographics: %s', id[1], gender, demo_gender))
    #  }
    #}
    
    #grade
    grade <- subdata$grade
    # duration (in s, converted to min)
    #duration <-  subd$duration / 60 #not recorded here
    
    #dummy data frame
    dummy<-data.frame(id=id, age_years=age_years, gender=gender, grade=grade, x=x, y=y, chosen=chosen, z=z, zscaled=zscaled, time=time, trial=trial, round=round)
    #bind them
    dat<-rbind(dat, dummy)
    
    ############################################
    # save environment order
    id_env = rep(as.numeric(subd$uid), 10)
    round_env = 1:10
    envOrder = subdata$envOrder
    
    dummy_env = data.frame(id=id_env, round=round_env, env=envOrder)
    envs = rbind(envs, dummy_env)
    
  }
  
  ##########################################
  # compute distance between choices ------
  ##########################################
  dat <- dat %>% 
    arrange(id, round, trial) %>% # make sure rounds, trials etc are in the correct order
    mutate(distance = NA)
  
  # compute (Manhattan) distance between consecutive choices
  for(i in 1:(dim(dat)[1]-1)){
    dat$distance[i+1] <- dist(rbind(c(dat$x[i], dat$y[i]), c(dat$x[i+1], dat$y[i+1])), method = "manhattan")
  }
  
  # set distance for initial trial (=randomly revealed tile) to NA and classify distances
  dat <- dat %>% 
    mutate(distance = ifelse(trial == 0, NA, distance),
           type_choice = case_when(
             distance == 0 ~ "Repeat",
             distance == 1 ~ "Near",
             distance >1   ~ "Far"
             #is.na() ~ "nope",
           )) %>% 
    mutate(type_choice = factor(type_choice, levels = c('Repeat', 'Near', 'Far')))
  
  # distance as function of reward on previous trial ------------------------
  dat$previous_reward <- NA
  
  # add column with reward obtained on previous step
  for(i in 1:nrow(dat)){
    if(dat$trial[i] != 0) {
      dat$previous_reward[i] <-  dat$z[i-1]
    }
  }
  
  if (writecsv==TRUE){
    write.table(dat, file="data/AdolescentGrid.csv", sep=",", row.names = F)
    write.table(envs, file="data/AdolescentEnvironments.csv", sep=",", row.names = F)  
  }
  return(dat)
}


cohensd.ci <- function(d, n1, n2, ci = 0.95) {
  t <- d * sqrt((n1 * n2)/(n1 + n2))
  capture.output(
    fit <- compute.es::tes(t = t, n.1 = n1, n.2 = n2, level = 100 * ci),
    file = "NUL"
  )
  c(lower.ci = fit$l.d, upper.ci = fit$u.d)
}

#############################################################################################################################
# data from Schulz et al. (2019)
dataImport_Schulz = function(datafile='data/kwg_Schulz.json', writecsv=FALSE) {
  # cf https://github.com/ericschulz/kwg/blob/master/C_Code/C1behavioraltests.R
  #read in data
  myjson<-fromJSON(datafile)
  
  all_opts = expand.grid(0:7, 0:7)
  
  #initialize empty data frame
  dat<-data.frame(id=numeric(), cond=numeric(), age=numeric(), x=numeric(), y=numeric(), 
                  z=numeric(), zscaled=numeric(), chosen=numeric(), time=numeric(),
                  trial=numeric(), round=numeric(), gender=character())
  
  #loop through json
  for (i in 1:myjson$count){
    #x-y-z
    x<-as.vector(t(myjson$records$data$searchHistory$xcollect[[i]][2:9,]))
    y<-as.vector(t(myjson$records$data$searchHistory$ycollect[[i]][2:9,]))
    z<-as.vector(t(myjson$records$data$searchHistory$zcollect[[i]][2:9,]))
    chosen <- apply(cbind(x,y),MARGIN=1, FUN=function(row) which(row[1]==all_opts$Var1 & row[2]==all_opts$Var2))
    zscaled <- as.vector(t(myjson$records$data$searchHistory$zcollectScaled[[i]][2:9,]))
    #time
    time<- as.vector(t(myjson$records$data$searchHistory$tscollect[[i]][2:9,]))
    #condition
    cond<-rep(myjson$records$data$condition[i], length(x))
    #age
    age<-rep(myjson$records$data$age[i], length(x))
    #trial number
    trial<-rep(0:25, 8)
    #round number
    round<-rep(2:9, each=26)
    #id
    id<-rep(i, length(x))
    #gender
    gender<-rep(myjson$records$data$gender[i], length(x))
    #dummy frame
    dummy<-data.frame(id=id, cond=cond, age=age, x=x, y=y, z=z, zscaled=zscaled, chosen=chosen,
                      time=time, trial=trial, round=round, gender=gender)
    #bind them
    dat<-rbind(dat, dummy)
  }
  
  ##########################################
  # compute distance between choices ------
  ##########################################
  # compute (Manhattan) distance between consecutive choices
  for(i in 1:(dim(dat)[1]-1)){
    dat$distance[i+1] <- dist(rbind(c(dat$x[i], dat$y[i]), c(dat$x[i+1], dat$y[i+1])), method = "manhattan")
  }
  
  # set distance for initial trial (=randomly revealed tile) to NA and classify distances
  dat <- dat %>% 
    mutate(distance = ifelse(trial == 0, NA, distance),
           type_choice = case_when(
             distance == 0 ~ "Repeat",
             distance == 1 ~ "Near",
             distance >1   ~ "Far"
             #is.na() ~ "nope",
           )) %>% 
    mutate(type_choice = factor(type_choice, levels = c('Repeat', 'Near', 'Far')))
  
  # distance as function of reward on previous trial ------------------------
  dat$previous_reward <- NA
  
  # add column with reward obtained on previous step
  for(i in 1:nrow(dat)){
    if(dat$trial[i] != 0) {
      dat$previous_reward[i] <-  dat$z[i-1]
    }
  }
  
  
  #first 10 were us!
  dat<-subset(dat, id>10)
  #remove younger than 7
  dat<-subset(dat, age>=7)
  
  #age group
  dat$agegroup<-ifelse(dat$age<9, "7-8", dat$age)
  dat$agegroup<-ifelse(dat$age>=9 & dat$age <12, "9-11", dat$agegroup)
  dat$agegroup<-ifelse(dat$age>18, ">18", dat$agegroup)
  
  #condition
  dat$cond<-ifelse(dat$cond==1, "Rough", "Smooth")
  
  #new ids so they match the model comparison results
  ids<-1:160
  uids<-unique(dat$id)

  for (i in 1:nrow(dat)){
    dat$id[i]<-ids[match(dat$id[i], uids)]
  }
  
  
  if (writecsv==TRUE){
    write.table(dat, file="data/kwgdata_Schulz.csv", sep=",", row.names = F)  
  }
  
  return(dat)
}


#############################################################################################################################
# combine behavioral datasets from all experiments
dataImport = function(writecsv=FALSE, bothConds=FALSE) {
  # adolescent data
  # d <- dataImport_Adolescent()
  # import preprocessed data
  dA = read.csv('data/AdolescentGrid.csv')

  # exclude training round and bonus round from analysis
  dA = subset(dA, !round==1 & !round==10)

  dA$experiment = 'Adolescent'
  dA$age_months = dA$age_years*12
  dA$condition = 'Smooth'

  # total duration
  dDuration = data.frame()
  for (i in unique(dA$id)) {
    subd = subset(dA, id==i)

    start = head(subd$time, n=1)
    end = tail(subd$time, n=1)

    cur = data.frame(id=i, duration=(end-start)/60000) # time in minutes
    dDuration = rbind(dDuration, cur)
  }

  # include duration per round in data frame
  dA = merge(dA, dDuration, by='id')


  # Schulz et al. (2019)
  # training and bonus round are already excluded
  #dS = dataImport_Schulz()
  dS = read.csv('data/kwgdata_Schulz.csv')

  dS$experiment = 'Schulz (2019)'
  dS$age_months = dS$age*12

  # total duration
  dDuration = data.frame()
  for (i in unique(dS$id)) {
    subd = subset(dS, id==i)

    start = head(subd$time, n=1)
    end = tail(subd$time, n=1)

    cur = data.frame(id=i, duration=(end-start)/60000) # time in minutes
    dDuration = rbind(dDuration, cur)
  }

  # include duration per round in data frame
  dS = merge(dS, dDuration, by='id')


  # Meder et al. (in press)
  dM = read.csv('data/ykwgdata_Meder.csv')

  # exclude training round and bonus round from analysis
  # exclude rough environments
  dM = subset(dM, !round==1 & !round==6)

  dM$experiment = 'Meder (2021)'
  dM$time = NA

  # combine data frames
  # use ongoing id to separate participants
  dSId = unique(dS$id)
  dSIdNew = (tail(unique(dA$id),n=1)+1) : (tail(unique(dA$id),n=1)+length(dSId))
  for (i in 1:nrow(dS)) {
    dS$id[i] = dSIdNew[match(dS$id[i], dSId)]
  }

  dMId = unique(dM$id)
  dMIdNew = (tail(unique(dS$id),n=1)+1) : (tail(unique(dS$id),n=1)+length(dMId))
  for (i in 1:nrow(dM)) {
    dM$id[i] = dMIdNew[match(dM$id[i], dMId)]
  }

  # drop unused columns that are not available in all data frames
  dA = dA[, !(names(dA) %in% c('grade'))]
  dS = dS[, !(names(dS) %in% c('agegroup'))]
  dM = dM[, !(names(dM) %in% c('grade', 'cond', 'agegroup', 'agegroup2'))]

  dS = dS %>% rename(age_years = age,
                     condition = cond)
  dM = dM %>% rename(condition = Condition)
  dM$gender = mapvalues(dM$gender, from = c('female', 'male'), to = c('Female', 'Male'))

  d = rbind(dA, dS, dM)
  
  if (bothConds==FALSE) {
    # only smooth condition
    d = subset(d, condition=='Smooth')
  }

  # bin data by age
  # set up cut-off values 
  breaks <- c(5, 7, 9, 11, 14, 18, 25, 55)
  # quantile(d$age_months/12, probs = seq(0, 1, by=1/6))
  
  # bucketing values into bins
  d$agegroup <- cut(d$age_months/12,
                    breaks=breaks, 
                    include.lowest=TRUE, 
                    right=FALSE)


  # save complete data frame
  if (writecsv == TRUE) {
    write.table(d, file="data/behavioralData.csv", sep=",", row.names = F)
  }
  
  return(d)
}

####################################################################################
# combine behavioral data (raw) from all experiments
allGridData = function(writecsv = FALSE) {
  # import all behavioral data (from all experiments)
  dMeder = read.csv('data/ykwgdata_Meder.csv')
  dMeder$experiment = 'Meder (2021)'
  dMeder = dMeder[, !(names(dMeder) %in% c('grade', 'agegroup', 'agegroup2', 'age_months', 'cond'))]
  names(dMeder)[names(dMeder) == 'duration'] = 'time'
  names(dMeder)[names(dMeder) == 'Condition'] = 'cond'
  
  dSchulz = read.csv('data/kwgdata_Schulz.csv')
  dSchulz$experiment = 'Schulz (2019)'
  dSchulz = dSchulz[, !(names(dSchulz) %in% c('agegroup'))]
  names(dSchulz)[names(dSchulz) == 'age'] = 'age_years'
  
  dAdolescent = read.csv('data/AdolescentGrid.csv')
  dAdolescent$experiment = 'Adolescent'
  dAdolescent$cond = 'Smooth'
  dAdolescent = dAdolescent[, !(names(dAdolescent) %in% c('grade'))]
  
  # use ongoing id
  SchulzId = unique(dSchulz$id)
  SchulzIdNew = (tail(unique(dAdolescent$id),n=1)+1) : (tail(unique(dAdolescent$id),n=1)+length(SchulzId))
  for (i in 1:nrow(dSchulz)) {
    dSchulz$id[i] = SchulzIdNew[match(dSchulz$id[i], SchulzId)]
  }
  
  MederId = unique(dMeder$id)
  MederIdNew = (tail(unique(dSchulz$id),n=1)+1) : (tail(unique(dSchulz$id),n=1)+length(MederId))
  for (i in 1:nrow(dMeder)) {
    dMeder$id[i] = MederIdNew[match(dMeder$id[i], MederId)]
  }
  
  # combine data
  dat = rbind(dAdolescent, dSchulz, dMeder)
  
  # filter for smooth environment condition
  dat = subset(dat, cond=='Smooth')
  
  # save complete data frame
  if (writecsv == TRUE) {
    write.table(dat, file="data/kwgdata_all.csv", sep=",", row.names = F)
  }
}

#############################################################################################################################
# Model Results
#############################################################################################################################
# imports and preprocesses model results from adolescent data
importModelResults <- function(dataFolder, kernels, acqFuncs){
  #Participant data
  #data<-dataImport()
  # import preprocessed data
  # data = read.csv('data/AdolescentGrid.csv')
  data = read.csv('data/kwgdata_all.csv')
  
  uids = unique(data$id)
  
  #initialize data frames
  modelFit <- data.frame(id=numeric(), nLL=numeric(), kernel=numeric(), acq=numeric(), R2=numeric(), kError=numeric(), lambda=numeric(), beta=numeric(), tau=numeric(), epsilon=numeric()) 
  paramEstimates <- data.frame(id=numeric(), leaveoutindex=numeric(), nLL=numeric(), R2=numeric(), kernel=numeric(), acq=numeric(), kError=numeric(), lambda=numeric(), beta=numeric(), tau=numeric(), epsilon=numeric(), roundnLL=numeric())
  #loop through data
  for (k in kernels){
    for (a in acqFuncs){
      # for (i in 1:length(uids)){ #subjects
      for (i in uids){ #subjects
        filename <- paste0(dataFolder, k, a, i, ".csv") #read output file
        if (file.exists(filename)){
          dp<-read.csv(filename)
          #print(filename)
          #interpret parameters based on model combination
          if (k==""){#Heuristics
            colnames(dp) <- c("", "leaveoutindex", "nLL", "tau")
          }else if (k=="BMT" | k=="LBMT"){#Bayesian mean tracker
            ifelse(a=='UCB' | a=='Counts', colnames(dp)<- c("", "leaveoutindex", "nLL", "kError", "beta","tau"), colnames(dp)<- c("",  "leaveoutindex", "nLL","kError", "tau"))
          }else if (k=="LIN" | k=="LLIN"){ #linear kernel
            ifelse(a=='UCB', colnames(dp)<- c("",  "leaveoutindex", "nLL", "beta","tau"), colnames(dp)<- c("", "leaveoutindex", "nLL", "tau"))
          }else { #normal GP kernel
            if (a=='UCB' | a=='Counts') {
              colnames(dp)<- c("", "leaveoutindex", "nLL", "lambda", "beta","tau")
            } else if (a=='EG') {
              colnames(dp)<- c("", "leaveoutindex", "nLL", "lambda", "beta", "epsilon")
            } else {
              colnames(dp)<- c("", "leaveoutindex", "nLL", "lambda", "tau")
            }
          }
          
          rounds = length(dp$nLL)
          
          #demographics
          dummy <- subset(data, id==i) #subset participant in the dataframe
          #environment <- dummy$Condition[1]
          id <- dummy$id[1]  #replicate ID
          kernel <- k
          acq <- a
          #Total out of sample nLL
          nLL <- sum(dp$nLL)
          randomModelLL <- -log(1/64)*rounds*25
          R2 <- 1 - (nLL/randomModelLL)
          #blank median parameter estimates
          kErrorMed <- NA
          lambdaMed <- NA
          betaMed <- NA
          tauMed <- NA
          epsilonMed <- NA
          #blank mean parameter estimates
          kErrorMean <- NA
          lambdaMean <- NA
          betaMean <- NA
          tauMean <- NA
          epsilonMean <- NA
          #mean parameter estimates for UCB RBF
          if (a=="UCB"| a=='Counts' | a=='EG'){ #UCB has beta
            betaMed <- median(exp(dp$beta))
            betaMean <- mean(exp(dp$beta))
          }
          if (k=="RBF" | k=="LRBF"){
            lambdaMed <- median(exp(dp$lambda))
            lambdaMean <- mean(exp(dp$lambda))
          }
          if (k=="BMT" | k=="LBMT"){ #BMT
            kErrorMed <- median(exp(dp$kError))
            kErrorMean <- mean(exp(dp$kError))
          }
          if (a!='EG') {
            tauMed <- median(exp(dp$tau))
            tauMean <- mean(exp(dp$tau))
          }
          if (a=='EG') {
            epsilon <- 1/(1+exp(-(dp$epsilon)))
            epsilonMed <- median(epsilon)
            epsilonMean <- mean(epsilon)
          }
          #save modelFit
          dadd <- data.frame(id=id, nLL=nLL, kernel=kernel, acq=acq, R2=R2, kError=kErrorMean, lambda=lambdaMean, beta = betaMean, tau=tauMean, epsilon=epsilonMean)
          names(dadd) <- names(modelFit)
          modelFit <-rbind(modelFit, dadd)
          #loop through leave out index to save all 9 parameter estimates for each subject
          for (loo in 2:(rounds+1)){
            subDP <- subset(dp, leaveoutindex == loo)
            roundnLL <- subDP$nLL
            #exponentiation of all parameters
            kError <- ifelse("kError" %in% colnames(subDP), exp(subDP$kError), NA)
            lambda <- ifelse("lambda" %in% colnames(subDP), exp(subDP$lambda), NA)
            beta <- ifelse("beta" %in% colnames(subDP), exp(subDP$beta),  NA)
            tau <- ifelse("tau" %in% colnames(subDP), exp(subDP$tau),  NA)
            epsilon <- ifelse("epsilon" %in% colnames(subDP), 1/(1+exp(-(dp$epsilon))), NA)
            dadd<-data.frame(id=id, leaveoutindex=loo, nLL=nLL, R2 =R2,  kernel=kernel, acq=acq, kError=kError, lambda=lambda, beta=beta, tau=tau, epsilon=epsilon, roundnLL=roundnLL)
            names(dadd) <- names(paramEstimates)
            paramEstimates <-rbind(paramEstimates, dadd)
          }}}}}
  return(list(modelFit, paramEstimates))
}


# imports and combines parameter estimated from all experiments
paramsImport = function(writecsv=FALSE) {
  # # behavioral data
  # # data<-dataImport()
  # # import preprocessed data
  # data = read.csv('data/kwgdata_all.csv')
  # 
  # # Read results from CSV files
  # modelResults <- importModelResults(dataFolder = 'modelResults/batch1/', kernels = c("BMT", "RBF"), acqFuncs = c("UCB", "GM", "GV", "EG"))
  # 
  # #separate into overall per participant and individual cross validation blocks
  # modelFit <- modelResults[[1]]
  # paramEstimates <- modelResults[[2]]
  # 
  # #reorder acqisition function levels and add "ModelName" to identify unique models
  # modelFit$acq <-factor(modelFit$acq)
  # modelFit$ModelName <- paste(modelFit$kernel,modelFit$acq, sep = "-")
  # modelFit$ModelName <- factor(modelFit$ModelName)
  # 
  # # add age in years to modelFit df
  # subData   <- data[!duplicated(data$id), ]
  # modelFit  <- left_join(modelFit, subData[c('id', 'age_years', 'gender', 'experiment')], by="id")
  # paramEstimates <- left_join(paramEstimates, subData[c('id', 'age_years', 'gender', 'experiment')], by="id")
  # 
  # #SAVE MODEL RESULTS TO DISK
  # write.csv(modelFit,'data/modelFit_new.csv')
  # write.csv(paramEstimates,'modelResults/paramEstimates_new.csv')
  
  # load if model results were saved to disk before
  meanParamsNew = read.csv('data/modelFit_new.csv')
  meanParamsNew$age_months = NA
  meanParamsNew$condition = 'Smooth'
  
  # paramsAdolescents = read.csv('modelResults/paramEstimates_Adolescents.csv')
  # paramsAdolescents = paramsAdolescents[, !names(paramsAdolescents) %in% c('X')]

  
  # data from Schulz et al. (2019)
  # GP-UCB
  paramsSchulzGP = read.csv('data/rbfucb_Schulz.csv')
  paramsSchulzGP = paramsSchulzGP %>% dplyr::rename(lambda = par1, beta = par2, tau = par3, leaveoutindex = X.1, nLL = X.2)
  paramsSchulzGP = paramsSchulzGP[, !(names(paramsSchulzGP) %in% c('X', 'X.3'))]
  paramsSchulzGP$kError = NA
  paramsSchulzGP$kernel = 'RBF'
  paramsSchulzGP$acq = 'UCB'
  
  # GP-GM
  paramsSchulzGPGM = read.csv('data/rbfgreedymean_Schulz.csv')
  paramsSchulzGPGM = paramsSchulzGPGM %>% dplyr::rename(lambda = par1, tau = par3, leaveoutindex = X.1, nLL = X.2)
  paramsSchulzGPGM = paramsSchulzGPGM[, !(names(paramsSchulzGPGM) %in% c('X', 'X.3', 'par2'))]
  paramsSchulzGPGM$kError = NA
  paramsSchulzGPGM$beta = NA
  paramsSchulzGPGM$kernel = 'RBF'
  paramsSchulzGPGM$acq = 'GM'
  
  # GP-GV
  paramsSchulzGPGV = read.csv('data/rbfgreedyvariance_Schulz.csv')
  paramsSchulzGPGV = paramsSchulzGPGV %>% dplyr::rename(lambda = par1, tau = par3, leaveoutindex = X.1, nLL = X.2)
  paramsSchulzGPGV = paramsSchulzGPGV[, !(names(paramsSchulzGPGV) %in% c('X', 'X.3', 'par2'))]
  paramsSchulzGPGV$kError = NA
  paramsSchulzGPGV$beta = NA
  paramsSchulzGPGV$kernel = 'RBF'
  paramsSchulzGPGV$acq = 'GV'

  # BMT-UCB
  paramsSchulzBMT = read.csv('data/bmt_Schulz.csv')
  paramsSchulzBMT$r2 = (1-paramsSchulzBMT$X.2/(-25*log(1/64)))
  paramsSchulzBMT = paramsSchulzBMT %>% dplyr::rename(kError = par1, beta = par2, tau = par3, leaveoutindex = X.1, nLL = X.2)
  paramsSchulzBMT = paramsSchulzBMT[, !(names(paramsSchulzBMT) %in% c('X', 'X.3'))]
  paramsSchulzBMT$lambda = NA
  paramsSchulzBMT$kernel = 'BMT'
  paramsSchulzBMT$acq = 'UCB'
  
  # BMT-GM
  paramsSchulzBMTGM = read.csv('data/mtgreedymean_Schulz.csv')
  paramsSchulzBMTGM = paramsSchulzBMTGM %>% dplyr::rename(kError = par1, tau = par3, leaveoutindex = X.1, nLL = X.2)
  paramsSchulzBMTGM = paramsSchulzBMTGM[, !(names(paramsSchulzBMTGM) %in% c('X', 'X.3', 'par2'))]
  paramsSchulzBMTGM$lambda = NA
  paramsSchulzBMTGM$beta = NA
  paramsSchulzBMTGM$kernel = 'BMT'
  paramsSchulzBMTGM$acq = 'GM'
  
  # BMT-GV
  paramsSchulzBMTGV = read.csv('data/mtgreedyvariance_Schulz.csv')
  paramsSchulzBMTGV = paramsSchulzBMTGV %>% dplyr::rename(kError = par1, tau = par3,leaveoutindex = X.1, nLL = X.2)
  paramsSchulzBMTGV = paramsSchulzBMTGV[, !(names(paramsSchulzBMTGV) %in% c('X', 'X.3', 'par2'))]
  paramsSchulzBMTGV$lambda = NA
  paramsSchulzBMTGV$beta = NA
  paramsSchulzBMTGV$kernel = 'BMT'
  paramsSchulzBMTGV$acq = 'GV'

  # combine datasets
  paramsSchulz = rbind(paramsSchulzGP, paramsSchulzGPGM, paramsSchulzGPGV,
                       paramsSchulzBMT, paramsSchulzBMTGM, paramsSchulzBMTGV)
  
  paramsSchulz$lambda = exp(paramsSchulz$lambda)
  paramsSchulz$kError = exp(paramsSchulz$kError)
  paramsSchulz$beta = exp(paramsSchulz$beta)
  paramsSchulz$tau = exp(paramsSchulz$tau)
  
  # mean parameter estimates
  meanParamsSchulz = ddply(paramsSchulz, ~id+kernel+acq, plyr::summarize, R2=mean(r2), nLL=sum(nLL), lambda=mean(lambda),
                           kError=mean(kError), beta=mean(beta), tau=mean(tau))
  
  # include demographic information
  dataSchulz = read.csv('data/kwgdata_Schulz.csv')
  meanParamsSchulz$condition = 0
  meanParamsSchulz$age_years = 0
  meanParamsSchulz$gender = 0
  for (i in 1:nrow(meanParamsSchulz)){
    meanParamsSchulz$condition[i] = dataSchulz$cond[which(meanParamsSchulz$id[i]==dataSchulz$id)[1]]
    meanParamsSchulz$age_years[i] = dataSchulz$age[which(meanParamsSchulz$id[i]==dataSchulz$id)[1]]
    meanParamsSchulz$gender[i] = dataSchulz$gender[which(meanParamsSchulz$id[i]==dataSchulz$id)[1]]
  }
  meanParamsSchulz$age_months = NA
  meanParamsSchulz$experiment = 'Schulz (2019)'
  meanParamsSchulz$epsilon = NA


  # data from Meder et al. (2021)
  paramsMeder = read.csv('data/modelFit_Meder.csv')
  paramsMeder = paramsMeder %>% dplyr::rename(id = participant, condition = environment)
  paramsMeder = subset(paramsMeder, acq!='Counts')
  paramsMeder$experiment = 'Meder (2021)'
  paramsMeder$epsilon = NA

  
  # use ongoing id
  SchulzId = unique(meanParamsSchulz$id)
  SchulzIdNew = (152 : (tail(unique(meanParamsNew$id),n=1)+length(SchulzId))) #150 adolescent participants (id 1-151, id 31 is missing)
  for (i in 1:nrow(meanParamsSchulz)) {
    meanParamsSchulz$id[i] = SchulzIdNew[match(meanParamsSchulz$id[i], SchulzId)]
  }

  MederId = unique(paramsMeder$id)
  MederIdNew = (tail(unique(meanParamsSchulz$id),n=1)+1) : (tail(unique(meanParamsSchulz$id),n=1)+length(MederId))
  for (i in 1:nrow(paramsMeder)) {
    paramsMeder$id[i] = MederIdNew[match(paramsMeder$id[i], MederId)]
  }


  # combine data
  keep = c('id', 'age_years', 'age_months', 'kernel', 'R2', 'nLL', 'kError', 'lambda', 'beta', 'tau', 'epsilon', 'experiment', 'condition', 'acq')
  meanParamsNew = meanParamsNew[, (names(meanParamsNew) %in% keep)]
  meanParamsSchulz = meanParamsSchulz[, (names(meanParamsSchulz) %in% keep)]
  paramsMeder = paramsMeder[, (names(paramsMeder) %in% keep)]

  params = rbind(meanParamsNew, meanParamsSchulz, paramsMeder)
  
  params = subset(params, condition=='Smooth')
  # add age month if not available in the data
  params$age_months = ifelse(is.na(params$age_months), params$age_years*12, params$age_months)
  
  # bin data by age
  # set up cut-off values 
  breaks <- c(5, 7, 9, 11, 14, 18, 25, 55)
  
  # bucketing values into bins
  params$agegroup <- cut(params$age_months/12,
                         breaks=breaks, 
                         include.lowest=TRUE, 
                         right=FALSE)

  # save complete data frame
  if (writecsv == TRUE) {
    write.table(params, file="data/modelFit.csv", sep=",", row.names = F)
  }
  
  return(params)
}


paramsYoungestAgegroup = function(writecsv = FALSE) {
  params = read.csv('data/paramEstimates_Meder.csv')
  
  paramsAvg = read.csv('data/modelFit_Meder.csv')
  
  # add age
  params = merge(params, paramsAvg[, c('participant', 'age_years', 'age_months')], by='participant')
  
  # change ids, such that they match ids in the other data frames
  ids = data.frame(participant = unique(params$participant),
                   id = c(312:(312+length(unique(params$participant))-1)))
  params = merge(params, ids, by='participant')
  
  # only smooth environments and youngest age group, only GP-UCB model
  params = subset(params, environment=='Smooth' & age_years<=6 & kernel=='RBF' & acq=='UCB')
  
  # remove duplicate columns
  params = params[!duplicated(params), ]
  
  if (writecsv == TRUE) {
    write.table(params, file="data/paramsYoungestAgegroup.csv", sep=",", row.names = F)
  }
}
