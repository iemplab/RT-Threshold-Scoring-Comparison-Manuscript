##############################################################
#
#       RESPONSE-TIME THRESHOLD PROCEDURE COMPARISONS        #
#
##############################################################

#clearing memory
rm(list=ls())
#T1 <- Sys.time()

######### Specifying packages needed for analyses ############

library(dplyr) #combines the simulation results into a condition matrix
library(SimDesign) #calculate RMSE
library(mirt) #estimating IRT item parameters
library(mvtnorm) #mvnorm function
library(MBESS) #cor2cov
library(plotrix) #calculate standard error for matrix-format output


##################### Variables that won't be changed ##############
Unmot.perc<-.28 #percentage of sample that engage in RG
I <- 50 # test length
reps<-100 #number of reps for study
reps.EM.imputed<-10 #number of reps for EM-imputation procedure

#################### creating conditions matrix ####################


# level codings will be specified in if else statements

# Sample size: 
sample.size <- c(500, 2500)

#PERCENTAGE OF RG RESPONSES FOR UNMOTIVATED SIMULEES
#there are three levels: 10%; 20% 
rg.percent <- .1

#MISCLASSIFICATION PERCENTAGE
#there are five levels: 2%, 4%, 6%, 8%, 12%
misclassify.percent <-c(0.02, 0.04, 0.06, 0.08, 0.12)

# correlation coefficients between latent abilities
r.ability.RG <- c(-0.1, -0.5, -0.9)

# We need to get a conditions matrix for a fully crossed design. In order to do this
# efficiently we can use the expand grid function which will give us a 48x5 matrix
# 2x2x4 = 16 levels crossed with the 4 variables: impact, NER relationship with ability, NER type, 
#percent unmotivated simulees in focal group, percent NER for unmotivated simulees in focal group
con <- expand.grid(sample.size, rg.percent, misclassify.percent, r.ability.RG)

#creating matrix for overall results by condition
overall.results<-matrix(-NA,nrow=nrow(con),ncol=75) #creating a matrix to place results in by condition, which will be used for the overall descriptive results

colnames(overall.results)<-c(
  
  "Bias.a.ML","Bias.a.penalized","Bias.a.EM","Bias.a.EM.imputed","Bias.a.HG",#1-5
  "Bias.b.ML","Bias.b.penalized","Bias.b.EM","Bias.b.EM.imputed","Bias.b.HG",#6-10
  
  "MAE.a.ML","MAE.a.penalized","MAE.a.EM","MAE.a.EM.imputed","MAE.a.HG",#11-15
  "MAE.b.ML","MAE.b.penalized","MAE.b.EM","MAE.b.EM.imputed","MAE.b.HG",#16-20
  
  "BIC.ML","BIC.penalized","BIC.EM","BIC.EM.imputed","BIC.HG",#21-25
  
  "Non.Converge.ML","Non.Converge.penalized","Non.Converge.EM","Non.Converge.EM.imputed","Non.Converge.HG",#26-30
  
  "Bias.unmot.theta.ML","Bias.unmot.theta.penalized","Bias.unmot.theta.EM","Bias.unmot.theta.EM.imputed","Bias.unmot.theta.HG",#31-35
  "MAE.unmot.theta.ML","MAE.unmot.theta.penalized","MAE.unmot.theta.EM","MAE.unmot.theta.EM.imputed","MAE.unmot.theta.HG",#36-40
  
  "Bias.mot.theta.ML","Bias.mot.theta.penalized","Bias.mot.theta.EM","Bias.mot.theta.EM.imputed","Bias.mot.theta.HG",#41-45
  "MAE.mot.theta.ML","MAE.mot.theta.penalized","MAE.mot.theta.EM","MAE.mot.theta.EM.imputed","MAE.mot.theta.HG",#46-50
  
  "rBias.a.ML","rBias.a.penalized","rBias.a.EM","rBias.a.EM.imputed","rBias.a.HG",#51-55
  "rBias.b.ML","rBias.b.penalized","rBias.b.EM","rBias.b.EM.imputed","rBias.b.HG",#56-60
  "rBias.unmot.theta.ML","rBias.unmot.theta.penalized","rBias.unmot.theta.EM","rBias.unmot.theta.EM.imputed","rBias.unmot.theta.HG",#61-65
  "rBias.mot.theta.ML","rBias.mot.theta.penalized","rBias.mot.theta.EM","rBias.mot.theta.EM.imputed","rBias.mot.theta.HG",#66-70
  
  
  "SE.theta.EM.imputed","Bias.r.theta.RG.HG","rBias.r.theta.RG.HG", "MAE.r.theta.RG.HG",#71-74
  
  "True.RG.Ability.Correlation" #75
)



#condition loop starts here

for(count.con in 1:nrow(con)) { 
  
  
  results.condition<-matrix(-NA,nrow=reps,ncol=75) #creating a matrix to place results in by condition, which will be used for the overall descriptive results
  
  
  # Variables that won't change across replications within a single condision
  N <- con[count.con, 1] #sample size
  theta.cor <- con[count.con, 4] #theta correlations
  N.unmot<-N*Unmot.perc #size of unmotivated sample
  N.mot<-N*(1-Unmot.perc) #size of motivated sample
  
  ######### REPLICATION LOOP STARTS HERE
  for (r in 1:reps){
    #setting seed
    set.seed(123+(r*count.con)+r) #rep to this counter
    
    ##################### Sampling item parameters ##############
    a.3PL.true<-as.matrix(runif(I, min = 0.4, max = 1.5), ncol=1)
    b.3PL.true<-as.matrix(rnorm(I),ncol=1)
    c.3PL.true<-as.matrix(runif(I, min = 0.1, max = 0.25),ncol=1)
    
    ##################### Sampling latent personal parameters (ability, RG propensity) ##############
    #sampling multivariate thetas
    ave.theta <- 0 #mean theta
    sd.theta <- 1 #standard deviation of theta
    
    # unmotivated subgroup
    cor.matrix.unmot<-matrix(c(1,theta.cor,theta.cor,1),ncol=2,nrow=2,byrow = TRUE) #specifying correlation matrix for thetas
    cov.matrix.unmot<-cor2cov(cor.matrix.unmot,c(sd.theta,sd.theta)) #converting correlation to covariance
    personal.par.unmot<-rmvnorm(N.unmot,mean=c(ave.theta,ave.theta),cov.matrix.unmot) #sampling thetas
    ability.theta.true.unmot <- as.matrix(personal.par.unmot[,1], ncol=1)
    RG.theta.true.unmot <- as.matrix(personal.par.unmot[,2], ncol=1)
    
    # motivated subgroupz: theta correlation = 0
    cor.matrix.mot<-matrix(c(1,0,0,1),ncol=2,nrow=2,byrow = TRUE) #specifying correlation matrix for thetas
    cov.matrix.mot<-cor2cov(cor.matrix.mot,c(sd.theta,sd.theta)) #converting correlation to covariance
    personal.par.mot<-rmvnorm(N.mot,mean=c(ave.theta,ave.theta),cov.matrix.mot) #sampling thetas
    ability.theta.true.mot <- as.matrix(personal.par.mot[,1], ncol=1)
    RG.theta.true.mot <- as.matrix(personal.par.mot[,2], ncol=1)
    
    theta.combined <- rbind(ability.theta.true.unmot, ability.theta.true.mot)
    guess.combined <- rbind(RG.theta.true.unmot, RG.theta.true.mot)
    theta.cor.true <- cor(theta.combined, guess.combined, use="complete.obs")
    
    results.condition[r,75] <- theta.cor.true
    
    ##################### Generating response probabilities ##############
    
    #creating function to compute 3pl (for item responses probability)
    PL3<-function(thetas,a.MC,b.MC,c.MC,N.3pl,I.MC){
      alphas.betas <-matrix(a.MC * b.MC,ncol=1)
      alphas.thetas <-matrix(thetas,ncol=1) %*% matrix(a.MC,nrow=1)
      logit <-alphas.thetas - t(matrix(alphas.betas,I.MC,N.3pl))
      c.par <- t(matrix(c.MC,I.MC,N.3pl))
      probs <- (c.par+((1-c.par)*(1/(1+exp(-logit)))))
      return(probs)}
    
    #creating function to compute Rasch (for RG probability)
    Rasch<-function(theta,b,N,I){
      betas <- t(matrix(b,I,N))
      thetas <- matrix(theta, ncol=1) %*% matrix(1, nrow=1, ncol = I)
      logit <-as.matrix((thetas - betas),ncol=1)
      probs <- (1)/(1+exp(-logit))
      return(probs)
    }
    
    #generating true response probability
    IR.prob.unmot <- PL3(ability.theta.true.unmot, a.3PL.true, b.3PL.true, c.3PL.true, N.unmot, I)
    IR.prob.mot <- PL3(ability.theta.true.mot, a.3PL.true, b.3PL.true, c.3PL.true, N.mot, I)
    
    #generating true RG probability
    RG.prob.unmot <- Rasch(RG.theta.true.unmot, b.3PL.true, N.unmot, I)
    RG.prob.mot <- Rasch(RG.theta.true.mot, b.3PL.true, N.mot, I)
    
    ######################### SAMPLING RANDOM NUMBER OF RGs FOR EACH UNMOTIVATED SIMULEE #######################
    #this is the random number of RGs for condition with 10% RGs
    
    aa<-N.unmot*I #total number of item responses (N x I)
    bb<-aa*con[count.con,2] #number of RG responses to sample
    prob.sample<-1/N.unmot #probability of sampling each simulee
    samples<-1 #number of samples of random numbers to pull 
    
    RG.1<-matrix(rmultinom(samples, size = bb,prob=c(rep(prob.sample,N.unmot))),ncol=samples,nrow=N.unmot)
    #sampling number of RG responses per simulee
    
    ######################### SAMPLING RANDOM NUMBER OF RGs FOR EACH MOTIVATED SIMULEE #######################
    #this is the random number of RGs for condition with 10% RGs
    
    aa.mot<-.3*N.mot*I #total number of item responses (N x I) .3=30% of top motivated simulees sampled for overclassification
    bb.mot<-aa.mot*con[count.con,3] #number of RG responses to sample *this group does not have any RG, so this needs to be the number of overclassification
    prob.sample.mot<-1/(N.unmot*.3) #probability of sampling each simulee
    samples<-1 #number of samples of random numbers to pull 
    
    RG.1.mot<-matrix(rmultinom(samples, size = bb.mot,prob=c(rep(prob.sample.mot,(N.mot*.3)))),ncol=samples,nrow=(N.mot*.3))
    #sampling number of RG responses per simulee
    
    ########################  CREATING RG FOR UNDERCLASSIFICATION CONDITION ##############################
    
    RG.matrix.mot<-matrix(0,nrow=N.mot,ncol=I) #CREATING RG Matrix for Motivated Sample (all motivated)
    
    RG.matrix.unmot<-matrix(0,nrow=N.unmot,ncol=I) #CREATING RG Matrix for Unmotivated Sample 
    
    for(p in 1:N.unmot){
      if(RG.1[p]>0){
      
      ############## SELECTING ITEMS TO RECEIVE RG ##################
      #rank ordering the true RG probability from lowest to highest
      rank.prob.RG.guesser<-order(RG.prob.unmot[p,],decreasing=TRUE) #order in descending order the true probability of success from lowest to highest; ties are randomly ordered 
      item.id<-rank.prob.RG.guesser[c(1:RG.1[p,])]
      
      ######################### UPDATING GUESSING PROBABILITIES & RG MATRIX ############################
      IR.prob.unmot[p,item.id]<-.25 #replacing original probability with .25 for the random item number #
      
      #replacing 1 with 0 to represent RG responses
      RG.matrix.unmot[p,item.id]<-1 #this represents accurate RG classifications
      }else if (RG.1[p]==0){
        
        #if there is no RG for an "unmotivated" simulee, then nothing changes
        IR.prob.unmot[p,]<-IR.prob.unmot[p,] #replacing original probability with .25 for the random item number #
        
        #replacing 1 with 0 to represent RG responses
        RG.matrix.unmot[p,]<-RG.matrix.unmot[p,] #this represents accurate RG classifications
      }
    }
    
    
    ###################### ADDING OVERCLASSIFICATION ############################
    
    #step 1: select RG & IR probabilities belonging to top 30 percentile of ability
    
    #RG probabilities
    RG.prob.mot.misclassify<-RG.prob.mot[which(ability.theta.true.mot>=quantile(ability.theta.true.mot,prob=.7)), ]
    RG.prob.mot.no.misclassify<-RG.prob.mot[which(ability.theta.true.mot<quantile(ability.theta.true.mot,prob=.7)), ]
    
    #IR probabilities
    IR.prob.mot.misclassify<-IR.prob.mot[which(ability.theta.true.mot>=quantile(ability.theta.true.mot,prob=.7)), ]
    IR.prob.mot.no.misclassify<-IR.prob.mot[which(ability.theta.true.mot<quantile(ability.theta.true.mot,prob=.7)), ]
    
    #subsetting RG matrix for motivation group based on misclassified and not misclassified subgroups
    
    RG.matrix.mot.misclassify<-RG.matrix.mot[1:(N.mot*.3),]    
    RG.matrix.mot.no.misclassify<-RG.matrix.mot[(((N.mot*.3)+1):N.mot),] 
    
    #step 2: reclassify true SB as RG based on RG.1.mot variable
    for(p in 1:nrow(RG.prob.mot.misclassify)){
      if(RG.1.mot[p]>0){  
      #sample number of items based on RG.1.mot variable
      
        item.id<-order(-IR.prob.mot.misclassify[p,])[1:RG.1.mot[p]] #here I ordered item difficulty from easiest to hardest; rg probability is not used because the covariance between ability and RG propensity is 0
        
      #step 3: reclassify RG matrix
      
      RG.matrix.mot.misclassify[p,item.id]<-1
      
      } else if(RG.1.mot[p]==0){
      #if there is no misclassifications for a given simulee, then the RG.matrix remains unchanged
          RG.matrix.mot.misclassify[p,]<-RG.matrix.mot.misclassify[p,]
      }
    } #closes p loop
    
    #step 4: create updated NA and Incorrect RG matrix for both motivated and unmotivated subgroups
    
    ############## NEED TO RECODE PROBABILITIES BASED ON UPDATED RG MATRIX
    
    #creating updated matrices for recoding
    
    #this code pertains to the unmotivated subgroup
    IR.prob.unmot.NA<-matrix(-999,nrow=N.unmot,ncol=I)
    IR.prob.unmot.zero<-matrix(-999,nrow=N.unmot,ncol=I)
    
    for(p in 1:N.unmot){
      for (i in 1:I){
        ifelse(RG.matrix.unmot[p,i]==1,IR.prob.unmot.NA[p,i]<-NA,IR.prob.unmot.NA[p,i]<-IR.prob.unmot[p,i])
        ifelse(RG.matrix.unmot[p,i]==1,IR.prob.unmot.zero[p,i]<-0,IR.prob.unmot.zero[p,i]<-IR.prob.unmot[p,i])
        
      }
    }
    
    #this code pertains to the motivated subgroup that is misclassified
    IR.prob.mot.misclassify.NA<-matrix(-999,nrow=(N.mot*.3),ncol=I)
    IR.prob.mot.misclassify.zero<-matrix(-999,nrow=(N.mot*.3),ncol=I)
    
    for(p in 1:(N.mot*.3)){
      for (i in 1:I){
        ifelse(RG.matrix.mot.misclassify[p,i]==1,IR.prob.mot.misclassify.NA[p,i]<-NA,IR.prob.mot.misclassify.NA[p,i]<-IR.prob.mot.misclassify[p,i])
        ifelse(RG.matrix.mot.misclassify[p,i]==1,IR.prob.mot.misclassify.zero[p,i]<-0,IR.prob.mot.misclassify.zero[p,i]<-IR.prob.mot.misclassify[p,i])
        
      }
    }
    
    
    
    ################################################################
    #
    #               COMBINING RG AND IR MATRICES                   #
    #
    ################################################################
    
    
    RG.combined<-rbind(RG.matrix.unmot,RG.matrix.mot.misclassify,RG.matrix.mot.no.misclassify) #combining RG classifications
    
    RG.combined[N,]<-1 #adding in case here to make sure that I can estimate model (need at least two response options for each item)
    
    #combining response probabilities for unmotivated and motivated samples
    probs.combined.ML<-rbind(IR.prob.unmot,IR.prob.mot)
    probs.combined.NA<-rbind(IR.prob.unmot.NA,IR.prob.mot.misclassify.NA,IR.prob.mot.no.misclassify)
    probs.combined.zero<-rbind(IR.prob.unmot.zero,IR.prob.mot.misclassify.zero,IR.prob.mot.no.misclassify)
    
    ########################################################################################################
    # ASSIGNING 1'S AND 0'S USING RANDOM NUMBERS #
    ########################################################################################################
    #getting random probabilities
    
    
    rnd.unif <- matrix(runif(N*I, 0, 1), nrow = N, ncol = I)
    
    #transforming probs to 0/1 for when noneffortful responses are in the data matrix
    IR.combined.ML <- (ifelse(probs.combined.ML > rnd.unif,1,0)) #Coding responses as 0/1 for IRT software#
    IR.combined.NA <- (ifelse(probs.combined.NA > rnd.unif,1,0)) #Coding responses as 0/1 for IRT software#
    IR.combined.zero <- (ifelse(probs.combined.zero > rnd.unif,1,0)) #Coding responses as 0/1 for IRT software#
    
    
    
    #adding column names for use in MIRT Package
    colnames(IR.combined.ML)<-c(1:I)
    colnames(IR.combined.NA)<-c(1:I)
    colnames(IR.combined.zero)<-c(1:I)
    
    
    #############################################################################
    #                                                                           #
    #                   ESTIMATING ITEM PARAMETERS & ABILITIES                  #
    #                                                                           #
    #############################################################################
    
    
    ###### ML scoring (ignoring RG responses) ########
    
    ML.tmp<-mirt(IR.combined.ML, 1, itemtype ='2PL',guess=.25, TOL = .0001, technical = list(NCYCLES = 10000), verbose = FALSE)
    
    if(extract.mirt(ML.tmp,'converged')==TRUE){ 
      ipars.ML.tmp <- coef(ML.tmp,IRTpars=TRUE,as.data.frame=TRUE) 
      a.ML.tmp<-ipars.ML.tmp[c(seq(1,200,4))] 
      b.ML.tmp<-ipars.ML.tmp[c(seq(2,200,4))] 
      a.ML.tmp[a.ML.tmp>3] <- 3
      a.ML.tmp[a.ML.tmp<(-3)] <- (-3)
      b.ML.tmp[b.ML.tmp>3] <- 3
      b.ML.tmp[b.ML.tmp<(-3)] <- (-3)
      modifiedpars<-mirt(IR.combined.ML, 1, pars = "values", itemtype ='2PL',guess=.25, TOL = .0001, technical = list(NCYCLES = 10000), verbose = FALSE)
      modifiedpars[c(seq(1,200,4)), "value"] <- a.ML.tmp
      modifiedpars[c(seq(2,200,4)), "value"] <- -(b.ML.tmp)*(a.ML.tmp)
      modifiedpars[c(seq(1,200,4)), "est"] <- FALSE
      modifiedpars[c(seq(2,200,4)), "est"] <- FALSE
      
      ML<-mirt(IR.combined.ML, 1, pars = modifiedpars, verbose=FALSE)
      ipars.ML<-coef(ML,IRTpars=TRUE,as.data.frame=TRUE) #need to read in item parameters
      results.condition[r,26]<-ifelse(extract.mirt(ML,'converged')==TRUE,0,1) #if model failed to converge, give 1; otherwise, 0
      if(results.condition[r,26]==0){ 
        a.ML<-ipars.ML[c(seq(1,200,4))] 
        b.ML<-ipars.ML[c(seq(2,200,4))] 
        theta.ML<-as.matrix(fscores(ML,method='ML',theta_lim = c(-4,4), verbose = FALSE))
        theta.ML[theta.ML>4 | theta.ML<(-4)] <- NA
        results.condition[r,21]<-anova(ML)$BIC
      }else{
        a.ML<-rep(NA,I) 
        b.ML<-rep(NA,I) 
        theta.ML<-rep(NA,N)
        results.condition[r,21]<-NA
      }
      
    }else{
      a.ML<-rep(NA,I) 
      b.ML<-rep(NA,I) 
      theta.ML<-rep(NA,N)
      results.condition[r,21]<-NA
      results.condition[r,26]<-1
    }
    
    ###### penalized scoring ######
    
    penalized.tmp<-mirt(IR.combined.zero, 1, itemtype ='2PL',guess=.25, TOL = .0001, technical = list(NCYCLES = 10000), verbose = FALSE)
    
    if(extract.mirt(penalized.tmp,'converged')==TRUE){ 
      ipars.penalized.tmp <- coef(penalized.tmp,IRTpars=TRUE,as.data.frame=TRUE) 
      a.penalized.tmp<-ipars.penalized.tmp[c(seq(1,200,4))] 
      b.penalized.tmp<-ipars.penalized.tmp[c(seq(2,200,4))] 
      a.penalized.tmp[a.penalized.tmp>3] <- 3
      a.penalized.tmp[a.penalized.tmp<(-3)] <- (-3)
      b.penalized.tmp[b.penalized.tmp>3] <- 3
      b.penalized.tmp[b.penalized.tmp<(-3)] <- (-3)
      modifiedpars<-mirt(IR.combined.zero, 1, pars = "values", itemtype ='2PL',guess=.25, TOL = .0001, technical = list(NCYCLES = 10000), verbose = FALSE)
      modifiedpars[c(seq(1,200,4)), "value"] <- a.penalized.tmp
      modifiedpars[c(seq(2,200,4)), "value"] <- -(b.penalized.tmp)*(a.penalized.tmp)
      modifiedpars[c(seq(1,200,4)), "est"] <- FALSE
      modifiedpars[c(seq(2,200,4)), "est"] <- FALSE
      
      penalized<-mirt(IR.combined.zero, 1, pars = modifiedpars, verbose=FALSE)
      ipars.penalized<-coef(penalized,IRTpars=TRUE,as.data.frame=TRUE) #need to read in item parameters
      results.condition[r,27]<-ifelse(extract.mirt(penalized,'converged')==TRUE,0,1) #if model failed to converge, give 1; otherwise, 0
      if(results.condition[r,27]==0){ 
        a.penalized<-ipars.penalized[c(seq(1,200,4))] 
        b.penalized<-ipars.penalized[c(seq(2,200,4))] 
        # b.penalized[b.penalized>6 | b.penalized<(-6)] <- NA
        #c.ML<-ipars.ML[c(seq(3,200,4))] 
        theta.penalized<-as.matrix(fscores(penalized,method='ML',theta_lim = c(-4,4), verbose = FALSE))
        theta.penalized[theta.penalized>4 | theta.penalized<(-4)] <- NA
        results.condition[r,22]<-anova(penalized)$BIC
      }else{
        a.penalized<-rep(NA,I) 
        b.penalized<-rep(NA,I) 
        #c.ML<-rep(NA,I) 
        theta.penalized<-rep(NA,N) 
        results.condition[r,22]<-NA
      }
      
    }else{
      a.penalized<-rep(NA,I) 
      b.penalized<-rep(NA,I) 
      #c.ML<-rep(NA,I) 
      theta.penalized<-rep(NA,N) 
      results.condition[r,22]<-NA
      results.condition[r,27]<-1
    }
    
    
    ###### EM scoring ######
    EM.tmp<-mirt(IR.combined.NA, 1, itemtype ='2PL',guess=.25, TOL = .0001, technical = list(NCYCLES = 10000), verbose = FALSE)
    
    if(extract.mirt(EM.tmp,'converged')==TRUE){ 
      ipars.EM.tmp <- coef(EM.tmp,IRTpars=TRUE,as.data.frame=TRUE) 
      a.EM.tmp<-ipars.EM.tmp[c(seq(1,200,4))] 
      b.EM.tmp<-ipars.EM.tmp[c(seq(2,200,4))] 
      a.EM.tmp[a.EM.tmp>3] <- 3
      a.EM.tmp[a.EM.tmp<(-3)] <- (-3)
      b.EM.tmp[b.EM.tmp>3] <- 3
      b.EM.tmp[b.EM.tmp<(-3)] <- (-3)
      modifiedpars<-mirt(IR.combined.NA, 1, pars = "values", itemtype ='2PL',guess=.25, TOL = .0001, technical = list(NCYCLES = 10000), verbose = FALSE)
      modifiedpars[c(seq(1,200,4)), "value"] <- a.EM.tmp
      modifiedpars[c(seq(2,200,4)), "value"] <- -(b.EM.tmp)*(a.EM.tmp)
      modifiedpars[c(seq(1,200,4)), "est"] <- FALSE
      modifiedpars[c(seq(2,200,4)), "est"] <- FALSE
      
      EM<-mirt(IR.combined.NA, 1, pars = modifiedpars, verbose=FALSE)
      ipars.EM<-coef(EM,IRTpars=TRUE,as.data.frame=TRUE) #need to read in item parameters
      results.condition[r,28]<-ifelse(extract.mirt(EM,'converged')==TRUE,0,1) #if model failed to converge, give 1; otherwise, 0
      if(results.condition[r,28]==0){ 
        a.EM<-ipars.EM[c(seq(1,200,4))] 
        b.EM<-ipars.EM[c(seq(2,200,4))] 
        c.EM<-ipars.ML[c(seq(3,200,4))] 
        theta.EM<-as.matrix(fscores(EM,method='ML',theta_lim = c(-4,4), verbose = FALSE))
        theta.EM[theta.EM > 4 | theta.EM < (-4)] <- NA
        results.condition[r,23]<-anova(EM)$BIC
      }else{
        a.EM<-rep(NA,I) 
        b.EM<-rep(NA,I) 
        c.EM<-rep(NA,I) 
        theta.EM<-rep(NA,N) 
        results.condition[r,23]<-NA
      }
      
    }else{
      a.EM<-rep(NA,I) 
      b.EM<-rep(NA,I) 
      c.EM<-rep(NA,I) 
      theta.EM<-rep(NA,N) 
      results.condition[r,23]<-NA
      results.condition[r,28]<-1
    }
    
    
    ##### EM-IMPUTATION scoring #####
    a.EM.imputed<-matrix(-999,nrow=I,ncol=reps.EM.imputed)
    b.EM.imputed<-matrix(-999,nrow=I,ncol=reps.EM.imputed)
    #c.EM.imputed<-matrix(-999,nrow=I,ncol=reps.EM.imputed)
    converged.EM.imputed<-matrix(-999,nrow=1,ncol=reps.EM.imputed)
    theta.EM.imputed<-matrix(-999,nrow=N,ncol=reps.EM.imputed)
    
    
    if(results.condition[r,28]==1){
      
      a.EM.imputed<-matrix(rep(matrix(rep(NA,I),ncol=1),reps.EM.imputed),ncol=reps.EM.imputed) 
      b.EM.imputed<-matrix(rep(matrix(rep(NA,I),ncol=1),reps.EM.imputed),ncol=reps.EM.imputed) 
      #c.ML<-rep(NA,I) 
      theta.EM.imputed<-matrix(rep(matrix(rep(NA,N),ncol=1),reps.EM.imputed),ncol=reps.EM.imputed)
      converged.EM.imputed<-matrix(1,nrow=1,ncol=reps.EM.imputed)
      
    } else if (results.condition[r,28]==0){
      #imputing probabilities
      imputed.probabilities.all<-PL3(theta.EM,a.EM,b.EM,c.EM,N,I) #based on the theta and item parameter estimates from EM scoring
      
      
      for(r.EM.impute in 1:reps.EM.imputed){
        rnd.unif.impute <- matrix(runif(N*I, 0, 1), nrow = N, ncol = I)
        #transforming probs to 0/1 for when noneffortful responses are in the data matrix
        IR.EM.imputed.all <- (ifelse(imputed.probabilities.all > rnd.unif.impute,1,0)) #Coding responses as 0/1 for IRT software#
        
        IR.EM.imputed <- matrix(-999, nrow = N, ncol = I)
        for (i in 1:I){
          for (n in 1:N){
            IR.EM.imputed[n,i] <- ifelse(is.na(IR.combined.NA[n,i])==TRUE, IR.EM.imputed.all[n,i], IR.combined.NA[n,i])
          }
        } 
        
        
        
        colnames(IR.EM.imputed)<-rep(1:I)
        
        
        EM.imputed<-mirt(IR.EM.imputed,  1, itemtype ='2PL',guess=.25, TOL = .0001, technical = list(NCYCLES = 10000), verbose = FALSE)
        ipars.EM.imputed<-coef(EM.imputed,IRTpars=TRUE,as.data.frame=TRUE) #need to read in item parameters
        
        converged.EM.imputed[,r.EM.impute]<-ifelse(extract.mirt(EM.imputed,'converged')==TRUE,0,1) #if model failed to converge, give 1; otherwise, 0
        
        #if model converged, take parameter estimates. if not, impute missing values
        if(extract.mirt(EM.imputed,'converged')==TRUE){ 
          a.EM.imputed[,r.EM.impute]<-ipars.EM.imputed[c(seq(1,200,4))] 
          b.EM.imputed[,r.EM.impute]<-ipars.EM.imputed[c(seq(2,200,4))] 
          #c.ML<-ipars.ML[c(seq(3,200,4))] 
          theta.EM.imputed[,r.EM.impute]<-as.matrix(fscores(EM.imputed,method='ML',theta_lim = c(-4,4), verbose = FALSE))
          results.condition[r,24]<-anova(EM.imputed)$BIC
          
        } else if (extract.mirt(EM.imputed,'converged')==FALSE){
          a.EM.imputed[,r.EM.impute]<-matrix(rep(NA,I),ncol=1) 
          b.EM.imputed[,r.EM.impute]<-rep(NA,I) 
          #c.ML<-rep(NA,I) 
          theta.EM.imputed[,r.EM.impute]<-rep(NA,N) 
          results.condition[r,24]<-NA
          
        }
        print(paste("EM.imputation",r.EM.impute,sep = " "))
        #print(r.EM.impute)
      }# close 10  reps loop
      
      # theta
      theta.EM.imputed[theta.EM.imputed==(-999)|theta.EM.imputed>4 | theta.EM.imputed<(-4)] <- NA
      theta.EM.imputed.1 <- theta.EM.imputed
      theta.EM.imputed.1[is.na(theta.EM.imputed)] <- 0
      
      var.w.EMI <- (std.error(theta.EM.imputed.1, na.rm = TRUE)^2)/(reps.EM.imputed)
      var.b.EMI <- (apply(theta.EM.imputed.1,MARGIN = 2, var, na.rm=TRUE))/(reps.EM.imputed-1)
      var.t.EMI <- var.w.EMI+ var.b.EMI+ (var.b.EMI/reps.EM.imputed)
      
    } #closes if/else statement above
    
    
    a.EM.imputed<-matrix(apply(a.EM.imputed,1,mean,na.rm=TRUE),ncol=1)
    b.EM.imputed<-matrix(apply(b.EM.imputed,1,mean,na.rm=TRUE),ncol=1)
    theta.EM.imputed<-matrix(apply(theta.EM.imputed,1,mean,na.rm=TRUE),ncol=1)
    
    results.condition[r,29]<-matrix(rowMeans(converged.EM.imputed,na.rm=TRUE),ncol=1) #average reps nonconverged
    
    
    
    ##### HG scoring ########      
    
    #COMBINING IR AND RG MATRICES TO ESTIMATE 2-FACTOR CORRELATED-TRAITS MODEL#
    
    
    #combining data for both factors into single dataset
    data.factor<-data.matrix(cbind(IR.combined.NA,RG.combined))
    
    colnames(data.factor)<-c(1:100) #adding in column names to estimate data in mirt r package
    
    #estimating multidimensional 2pl modelwith mirt.model definition
    model <- 'F1 = 1-50
    F2 = 51-100
    COV = F1*F2'
    
    item.type<-c(rep('2PL',50),rep('Rasch',50))
    HG<-mirt(data.factor, model,itemtype =item.type,guess=as.vector(c(rep(.25,50),rep(0,50))), TOL = .0001, technical = list(NCYCLES = 10000), verbose = FALSE)
    ipars.HG<-coef(HG,as.data.frame=TRUE)
    
    results.condition[r,30]<-ifelse(extract.mirt(HG,'converged')==TRUE,0,1) #if model failed to converge, give 1; otherwise, 0
    
    #if model converged, take parameter estimates. if not, impute missing values
    if(results.condition[r,30]==0){ 
      #need to transform parameter estimates onto theta scale
      a.HG<-ipars.HG[c(seq(1,250,5))] 
      b.HG<--(ipars.HG[c(seq(3,250,5))])/a.HG #I am not sure if there is a problem; trying to convert cfa parameters to irt b=threshold/slope
      #c.HG<-ipars.HG[c(seq(4,250,5))] 
      r.HG<-ipars.HG[(504)] #Covariance between theta and RG propensity
      theta.HG<-matrix(fscores(HG,method='ML',theta_lim = c(-4,4), verbose = FALSE)[,1],ncol=1) #setting upper and lower limits of theta at 4 and -4
      theta.HG[theta.HG>4 | theta.HG<(-4)] <- NA
      results.condition[r,25]<-anova(HG)$BIC
      
    } else {
      a.HG<-rep(NA,50) 
      b.HG<-rep(NA,50)  #I am not sure if there is a problem; trying to convert cfa parameters to irt b=threshold/slope
      #c.HG<-ipars.HG[c(seq(4,250,5))] 
      r.HG<-NA #Covariance between theta and RG propensity
      theta.HG<-rep(NA,N) #setting upper and lower limits of theta at 4 and -4
      results.condition[r,25]<-NA
      
    }
    
    
    ###################################
    #                                 #
    # CALCULATE DEPENDENT VARIABLES   #
    #                                 #
    ###################################
    #have to drop missing data to run simdesign functions
    
    #dropping missing cases for a.ML
    a.ML.data<-cbind(a.ML,a.3PL.true)
    a.ML.recoded<-matrix(a.ML.data[complete.cases(a.ML.data), ],ncol=2)
    colnames(a.ML.recoded)<-c("a.ML","a.True")
    
    #dropping missing cases for b.ML
    b.ML.data<-cbind(b.ML,b.3PL.true)
    b.ML.recoded<-matrix(b.ML.data[complete.cases(b.ML.data), ],ncol=2)
    colnames(b.ML.recoded)<-c("b.ML","b.True")
    
    #dropping missing cases for theta.ML
    theta.ML.data<-cbind(theta.ML,theta.combined)
    theta.ML.recoded<-matrix(theta.ML.data[complete.cases(theta.ML.data), ],ncol=2)
    theta.ML.recoded <- theta.ML.recoded[!is.infinite(rowSums(theta.ML.recoded)),] #removing any cases with inifite theta estimates
    colnames(theta.ML.recoded)<-c("theta.ML","theta.True")
    
    
    #dropping missing cases for a.penalized
    a.penalized.data<-cbind(a.penalized,a.3PL.true)
    a.penalized.recoded<-matrix(a.penalized.data[complete.cases(a.penalized.data), ],ncol=2)
    colnames(a.penalized.recoded)<-c("a.penalized","a.True")
    
    #dropping missing cases for b.penalized
    b.penalized.data<-cbind(b.penalized,b.3PL.true)
    b.penalized.recoded<-matrix(b.penalized.data[complete.cases(b.penalized.data), ],ncol=2)
    colnames(b.penalized.recoded)<-c("b.penalized","b.True")
    
    #dropping missing cases for theta.penalized
    theta.penalized.data<-cbind(theta.penalized,theta.combined)
    theta.penalized.recoded<-matrix(theta.penalized.data[complete.cases(theta.penalized.data), ],ncol=2)
    theta.penalized.recoded <- theta.penalized.recoded[!is.infinite(rowSums(theta.penalized.recoded)),] #removing any cases with inifite theta estimates
    colnames(theta.penalized.recoded)<-c("theta.penalized","theta.True")
    
    
    #dropping missing cases for a.EM
    a.EM.data<-cbind(a.EM,a.3PL.true)
    a.EM.recoded<-matrix(a.EM.data[complete.cases(a.EM.data), ],ncol=2)
    colnames(a.EM.recoded)<-c("a.EM","a.True")
    
    #dropping missing cases for b.EM
    b.EM.data<-cbind(b.EM,b.3PL.true)
    b.EM.recoded<-matrix(b.EM.data[complete.cases(b.EM.data), ],ncol=2)
    colnames(b.EM.recoded)<-c("b.EM","b.True")
    
    #dropping missing cases for theta.EM
    theta.EM.data<-cbind(theta.EM,theta.combined)
    theta.EM.recoded<-matrix(theta.EM.data[complete.cases(theta.EM.data), ],ncol=2)
    theta.EM.recoded <- theta.EM.recoded[!is.infinite(rowSums(theta.EM.recoded)),] #removing any cases with inifite theta estimates
    colnames(theta.EM.recoded)<-c("theta.EM","theta.True")
    
    
    #dropping missing cases for a.EM.imputed
    a.EM.imputed.data<-cbind(a.EM.imputed,a.3PL.true)
    a.EM.imputed.recoded<-matrix(a.EM.imputed.data[complete.cases(a.EM.imputed.data), ],ncol=2)
    colnames(a.EM.imputed.recoded)<-c("a.EM.imputed","a.True")
    
    #dropping missing cases for b.EM.imputed
    b.EM.imputed.data<-cbind(b.EM.imputed,b.3PL.true)
    b.EM.imputed.recoded<-matrix(b.EM.imputed.data[complete.cases(b.EM.imputed.data), ],ncol=2)
    colnames(b.EM.imputed.recoded)<-c("b.EM.imputed","b.True")
    
    #dropping missing cases for theta.EM.imputed
    theta.EM.imputed.data<-cbind(theta.EM.imputed,theta.combined)
    theta.EM.imputed.recoded<-matrix(theta.EM.imputed.data[complete.cases(theta.EM.imputed.data), ],ncol=2)
    theta.EM.imputed.recoded <- theta.EM.imputed.recoded[!is.infinite(rowSums(theta.EM.imputed.recoded)),] #removing any cases with inifite theta estimates
    colnames(theta.EM.imputed.recoded)<-c("theta.EM.imputed","theta.True")
    
    
    #dropping missing cases for a.HG
    a.HG.data<-cbind(a.HG,a.3PL.true)
    a.HG.recoded<-matrix(a.HG.data[complete.cases(a.HG.data), ],ncol=2)
    colnames(a.HG.recoded)<-c("a.HG","a.True")
    
    #dropping missing cases for b.HG
    b.HG.data<-cbind(b.HG,b.3PL.true)
    b.HG.recoded<-matrix(b.HG.data[complete.cases(b.HG.data), ],ncol=2)
    colnames(b.HG.recoded)<-c("b.HG","b.True")
    
    #dropping missing cases for theta.HG
    theta.HG.data<-cbind(theta.HG,theta.combined)
    theta.HG.recoded<-matrix(theta.HG.data[complete.cases(theta.HG.data), ],ncol=2)
    theta.HG.recoded <- theta.HG.recoded[!is.infinite(rowSums(theta.HG.recoded)),] #removing any cases with inifite theta estimates
    colnames(theta.HG.recoded)<-c("theta.HG","theta.True")
    
    r.HG.recoded <- matrix(cbind(r.HG, theta.cor.true), ncol=2)
    colnames(r.HG.recoded) <- c("cor.HG", "cor.True")
    
    ################ CALCULATING DEPENDENT VARIABLES ###########
    
    #ML
    if(nrow(a.ML.recoded)>0){
      results.condition[r,1]<-bias(a.ML.recoded[,"a.ML"],a.ML.recoded[,"a.True"])#standard bias
      results.condition[r,6]<-bias(b.ML.recoded[,"b.ML"],b.ML.recoded[,"b.True"])
      results.condition[r,11]<-MAE(a.ML.recoded[,"a.ML"],a.ML.recoded[,"a.True"],type="MAE")#MAE
      results.condition[r,16]<-MAE(b.ML.recoded[,"b.ML"],b.ML.recoded[,"b.True"],type="MAE")
      results.condition[r,51]<-bias(a.ML.recoded[,"a.ML"],a.ML.recoded[,"a.True"],type="relative")#relative bias
      results.condition[r,56]<-bias(b.ML.recoded[,"b.ML"],b.ML.recoded[,"b.True"],type="relative")
      
      results.condition[r,31]<-bias(theta.ML.recoded[1:N.unmot,"theta.ML"],theta.ML.recoded[1:N.unmot,"theta.True"])
      results.condition[r,36]<-MAE(theta.ML.recoded[1:N.unmot,"theta.ML"],theta.ML.recoded[1:N.unmot,"theta.True"],type="MAE")
      results.condition[r,61]<-bias(theta.ML.recoded[1:N.unmot,"theta.ML"],theta.ML.recoded[1:N.unmot,"theta.True"],type="relative")
      
      results.condition[r,41]<-bias(theta.ML.recoded[-(1:N.unmot),"theta.ML"],theta.ML.recoded[-(1:N.unmot),"theta.True"])
      results.condition[r,46]<-MAE(theta.ML.recoded[-(1:N.unmot),"theta.ML"],theta.ML.recoded[-(1:N.unmot),"theta.True"],type="MAE")
      results.condition[r,66]<-bias(theta.ML.recoded[-(1:N.unmot),"theta.ML"],theta.ML.recoded[-(1:N.unmot),"theta.True"],type="relative")
      
    } else {
      results.condition[r,1]<-NA
      results.condition[r,6]<-NA
      results.condition[r,11]<-NA
      results.condition[r,16]<-NA
      results.condition[r,51]<-NA
      results.condition[r,56]<-NA
      
      results.condition[r,31]<-NA
      results.condition[r,36]<-NA
      results.condition[r,61]<-NA
      
      results.condition[r,41]<-NA
      results.condition[r,46]<-NA
      results.condition[r,66]<-NA
    }
    
    #Penalized
    
    if(nrow(a.penalized.recoded)>0){
      results.condition[r,2]<-bias(a.penalized.recoded[,"a.penalized"],a.penalized.recoded[,"a.True"])
      results.condition[r,7]<-bias(b.penalized.recoded[,"b.penalized"],b.penalized.recoded[,"b.True"])
      results.condition[r,12]<-MAE(a.penalized.recoded[,"a.penalized"],a.penalized.recoded[,"a.True"],type="MAE")
      results.condition[r,17]<-MAE(b.penalized.recoded[,"b.penalized"],b.penalized.recoded[,"b.True"],type="MAE")
      results.condition[r,52]<-bias(a.penalized.recoded[,"a.penalized"],a.penalized.recoded[,"a.True"],type="relative")
      results.condition[r,57]<-bias(b.penalized.recoded[,"b.penalized"],b.penalized.recoded[,"b.True"],type="relative")
      
      results.condition[r,32]<-bias(theta.penalized.recoded[1:N.unmot,"theta.penalized"],theta.penalized.recoded[1:N.unmot,"theta.True"])
      results.condition[r,37]<-MAE(theta.penalized.recoded[1:N.unmot,"theta.penalized"],theta.penalized.recoded[1:N.unmot,"theta.True"],type="MAE")
      results.condition[r,62]<-bias(theta.penalized.recoded[1:N.unmot,"theta.penalized"],theta.penalized.recoded[1:N.unmot,"theta.True"],type="relative")
      
      results.condition[r,42]<-bias(theta.penalized.recoded[-(1:N.unmot),"theta.penalized"],theta.penalized.recoded[-(1:N.unmot),"theta.True"])
      results.condition[r,47]<-MAE(theta.penalized.recoded[-(1:N.unmot),"theta.penalized"],theta.penalized.recoded[-(1:N.unmot),"theta.True"],type="MAE")
      results.condition[r,67]<-bias(theta.penalized.recoded[-(1:N.unmot),"theta.penalized"],theta.penalized.recoded[-(1:N.unmot),"theta.True"],type="relative")
      
    } else {
      results.condition[r,2]<-NA
      results.condition[r,7]<-NA
      results.condition[r,12]<-NA
      results.condition[r,17]<-NA
      results.condition[r,52]<-NA
      results.condition[r,57]<-NA
      
      results.condition[r,32]<-NA
      results.condition[r,37]<-NA
      results.condition[r,62]<-NA
      
      results.condition[r,42]<-NA
      results.condition[r,47]<-NA
      results.condition[r,67]<-NA
    }
    
    
    
    #EM
    
    if(nrow(a.EM.recoded)>0){
      results.condition[r,3]<-bias(a.EM.recoded[,"a.EM"],a.EM.recoded[,"a.True"])
      results.condition[r,8]<-bias(b.EM.recoded[,"b.EM"],b.EM.recoded[,"b.True"])
      results.condition[r,13]<-MAE(a.EM.recoded[,"a.EM"],a.EM.recoded[,"a.True"],type="MAE")
      results.condition[r,18]<-MAE(b.EM.recoded[,"b.EM"],b.EM.recoded[,"b.True"],type="MAE")
      results.condition[r,53]<-bias(a.EM.recoded[,"a.EM"],a.EM.recoded[,"a.True"],type="relative")
      results.condition[r,58]<-bias(b.EM.recoded[,"b.EM"],b.EM.recoded[,"b.True"],type="relative")
      
      results.condition[r,33]<-bias(theta.EM.recoded[1:N.unmot,"theta.EM"],theta.EM.recoded[1:N.unmot,"theta.True"])
      results.condition[r,38]<-MAE(theta.EM.recoded[1:N.unmot,"theta.EM"],theta.EM.recoded[1:N.unmot,"theta.True"],type="MAE")
      results.condition[r,63]<-bias(theta.EM.recoded[1:N.unmot,"theta.EM"],theta.EM.recoded[1:N.unmot,"theta.True"],type="relative")
      
      results.condition[r,43]<-bias(theta.EM.recoded[-(1:N.unmot),"theta.EM"],theta.EM.recoded[-(1:N.unmot),"theta.True"])
      results.condition[r,48]<-MAE(theta.EM.recoded[-(1:N.unmot),"theta.EM"],theta.EM.recoded[-(1:N.unmot),"theta.True"],type="MAE")
      results.condition[r,68]<-bias(theta.EM.recoded[-(1:N.unmot),"theta.EM"],theta.EM.recoded[-(1:N.unmot),"theta.True"],type="relative")
      
    } else {
      
      results.condition[r,3]<-NA
      results.condition[r,8]<-NA
      results.condition[r,13]<-NA
      results.condition[r,18]<-NA
      results.condition[r,53]<-NA
      results.condition[r,58]<-NA
      
      results.condition[r,33]<-NA
      results.condition[r,38]<-NA
      results.condition[r,63]<-NA
      
      results.condition[r,43]<-NA
      results.condition[r,48]<-NA
      results.condition[r,68]<-NA
    }
    
    
    #EM Imputation
    
    if(nrow(a.EM.imputed.recoded)>0){
      
      results.condition[r,4]<-bias(a.EM.imputed.recoded[,"a.EM.imputed"],a.EM.imputed.recoded[,"a.True"])
      results.condition[r,9]<-bias(b.EM.imputed.recoded[,"b.EM.imputed"],b.EM.imputed.recoded[,"b.True"])
      results.condition[r,14]<-MAE(a.EM.imputed.recoded[,"a.EM.imputed"],a.EM.imputed.recoded[,"a.True"],type="MAE")
      results.condition[r,19]<-MAE(b.EM.imputed.recoded[,"b.EM.imputed"],b.EM.imputed.recoded[,"b.True"],type="MAE")
      results.condition[r,54]<-bias(a.EM.imputed.recoded[,"a.EM.imputed"],a.EM.imputed.recoded[,"a.True"],type="relative")
      results.condition[r,59]<-bias(b.EM.imputed.recoded[,"b.EM.imputed"],b.EM.imputed.recoded[,"b.True"],type="relative")
      
      results.condition[r,34]<-bias(theta.EM.imputed.recoded[1:N.unmot,"theta.EM.imputed"],theta.EM.imputed.recoded[1:N.unmot,"theta.True"])
      results.condition[r,39]<-MAE(theta.EM.imputed.recoded[1:N.unmot,"theta.EM.imputed"],theta.EM.imputed.recoded[1:N.unmot,"theta.True"],type="MAE")
      results.condition[r,64]<-bias(theta.EM.imputed.recoded[1:N.unmot,"theta.EM.imputed"],theta.EM.imputed.recoded[1:N.unmot,"theta.True"],type="relative")
      
      results.condition[r,44]<-bias(theta.EM.imputed.recoded[-(1:N.unmot),"theta.EM.imputed"],theta.EM.imputed.recoded[-(1:N.unmot),"theta.True"])
      results.condition[r,49]<-MAE(theta.EM.imputed.recoded[-(1:N.unmot),"theta.EM.imputed"],theta.EM.imputed.recoded[-(1:N.unmot),"theta.True"],type="MAE")
      results.condition[r,69]<-bias(theta.EM.imputed.recoded[-(1:N.unmot),"theta.EM.imputed"],theta.EM.imputed.recoded[-(1:N.unmot),"theta.True"],type="relative")
      
      results.condition[r,71] <- mean(var.t.EMI, na.rm = TRUE)#all
      
    } else {
      results.condition[r,4]<-NA
      results.condition[r,9]<-NA
      results.condition[r,14]<-NA
      results.condition[r,19]<-NA
      results.condition[r,54]<-NA
      results.condition[r,59]<-NA
      
      results.condition[r,34]<-NA
      results.condition[r,39]<-NA
      results.condition[r,64]<-NA
      
      results.condition[r,44]<-NA
      results.condition[r,49]<-NA
      results.condition[r,69]<-NA
      
      results.condition[r,71]<-NA
    }
    
    
    
    #HG
    
    if(nrow(a.HG.recoded)>0){
      results.condition[r,5]<-bias(a.HG.recoded[,"a.HG"],a.HG.recoded[,"a.True"])
      results.condition[r,10]<-bias(b.HG.recoded[,"b.HG"],b.HG.recoded[,"b.True"])
      results.condition[r,15]<-MAE(a.HG.recoded[,"a.HG"],a.HG.recoded[,"a.True"],type="MAE")
      results.condition[r,20]<-MAE(b.HG.recoded[,"b.HG"],b.HG.recoded[,"b.True"],type="MAE")
      results.condition[r,55]<-bias(a.HG.recoded[,"a.HG"],a.HG.recoded[,"a.True"],type="relative")
      results.condition[r,60]<-bias(b.HG.recoded[,"b.HG"],b.HG.recoded[,"b.True"],type="relative")
      
      results.condition[r,35]<-bias(theta.HG.recoded[1:N.unmot,"theta.HG"],theta.HG.recoded[1:N.unmot,"theta.True"])
      results.condition[r,40]<-MAE(theta.HG.recoded[1:N.unmot,"theta.HG"],theta.HG.recoded[1:N.unmot,"theta.True"],type="MAE")
      results.condition[r,65]<-bias(theta.HG.recoded[1:N.unmot,"theta.HG"],theta.HG.recoded[1:N.unmot,"theta.True"],type="relative")
      
      results.condition[r,72] <- bias(r.HG.recoded[,"cor.HG"], r.HG.recoded[,"cor.True"])
      results.condition[r,73] <- bias(r.HG.recoded[,"cor.HG"], r.HG.recoded[,"cor.True"], type="relative")
      results.condition[r,74] <- MAE(r.HG.recoded[,"cor.HG"], r.HG.recoded[,"cor.True"], type="MAE")
      
      results.condition[r,45]<-bias(theta.HG.recoded[-(1:N.unmot),"theta.HG"],theta.HG.recoded[-(1:N.unmot),"theta.True"])
      results.condition[r,50]<-MAE(theta.HG.recoded[-(1:N.unmot),"theta.HG"],theta.HG.recoded[-(1:N.unmot),"theta.True"],type="MAE")
      results.condition[r,70]<-bias(theta.HG.recoded[-(1:N.unmot),"theta.HG"],theta.HG.recoded[-(1:N.unmot),"theta.True"],type="relative")
      
    } else {
      results.condition[r,5]<-NA
      results.condition[r,10]<-NA
      results.condition[r,15]<-NA
      results.condition[r,20]<-NA
      results.condition[r,55]<-NA
      results.condition[r,60]<-NA
      
      results.condition[r,35]<-NA
      results.condition[r,40]<-NA
      results.condition[r,65]<-NA
      
      results.condition[r,72]<-NA
      results.condition[r,73]<-NA
      results.condition[r,74]<-NA
      
      results.condition[r,45]<-NA
      results.condition[r,50]<-NA
      results.condition[r,70]<-NA
    }
    
    
    print(paste("rep",r,sep = " "))
    #print(r) 
    
  } #closes the rep loop
  
  #place descriptive results for each condition into an overall matrix that will be printed out
  overall.results[count.con,1:75]<-apply(results.condition[,1:75],2,mean,na.rm = TRUE)
  
  print(paste("con",count.con,sep = " "))
  #print(count.con) #printing number of condition
  #count.con = count.con + 1
  
  
} #closes the condition loop; the parenthesis around the bracket will calculate the run time

##### WRITING OUT RESULTS #####
overall.results.output<-data.frame(cbind(matrix(con[,1],ncol=1),
                                         matrix(con[,2],ncol=1),
                                         matrix(con[,3],ncol=1),
                                         matrix(con[,4],ncol=1),
                                         overall.results))


names(overall.results.output)[names(overall.results.output) == "V1"]<-"Sample Size"                 
names(overall.results.output)[names(overall.results.output) == "V2"]<-"RG Percent"
names(overall.results.output)[names(overall.results.output) == "V3"]<-"Misclassify Percent" 
names(overall.results.output)[names(overall.results.output) == "V4"]<-"Thetas Correlation" 

#T2 <- Sys.time()
#(T2-T1)


write.csv(overall.results.output,"Overall Results - Overclassification.csv",row.names = FALSE)
