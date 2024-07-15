##################################################
# R script to simulate two-level data sets
# Author: Clémence Leyrat
# Adapted by: Andriana Kostouraki
# Last update: 2024-02-16
##################################################
#AK: to print the rho values in the excel results.
#to add a counter for the singularity warnings.

#install required packages if not already installed
if (!requireNamespace("lme4")) install.packages("lme4")
if (!requireNamespace("GLMMadaptive")) install.packages("GLMMadaptive")
if (!requireNamespace("boot")) install.packages("boot")

#load required libraries
library(lme4) #Mixed models
library(GLMMadaptive) #integrate over the random effects distribution; numerical method: Monte Carlo integration.
#(if time) to add the lines where we specify m, k, correlated, etc., outside the generation functions.
library(boot) #for the bootstrap application

#rh0 = 0.8 
#delta0 = 1
#a1 = 0.2
#a2 = -0.5
#a3 = 0.6
#eta0 = 0.1
#b1 = 0.2
#b2 = -1
#b3 = 2
r = 20 # number of bootstrap replicates

#load the bootstrap statistic function
setwd("N:\\ICON\\ICON_all\\Users\\Andriana Kostouraki\\PhD\\Latent unconfoundedness work\\simulation\\code")
source("bootstrap_statistic_draft.R",encoding = "UTF-8")

generation_bin<-function(m,k,correlated,model_trt,model_outcome,rh0=0.8,delta0=1,a1=0.2,a2=-0.5,a3=0.6,eta0=0.1,b1=0.2,b2=-1,b3=2){
  
  #Covariate generation
  #Generate the cluster IDs, cluster-level covariate and random effect (for the outcome and the treatment)
  dat_cl<-data.frame("ID_clus"=1:k,"V"=rnorm(k,0,1),"gamma_Y"=rnorm(k,0,1),"gamma_T"=rnorm(k,0,1),"gamma_cov"=rnorm(k,0,1))
  
  #Generate the patient IDs and individual covariates (1 continuous, 1 binary)
  dat_ind<-data.frame("ID"=rep(1:m,k),"ID_clus"=rep(1:k, each=m),"X1"=rnorm(m*k,0,1),"X2"=rbinom(m*k,1,0.5))
  
  #Merge
  dat<-merge(dat_cl,dat_ind,by="ID_clus")
  dat$X1_2<-rh0*dat$V+sqrt(1-rh0**2)*dat$X1 #e.g.: rh0 = 0.3 or 0.5 or 0.8
  
  if(correlated=="corr"){dat$X1<-dat$X1_2}
  if(correlated=="slope"){dat$X1<-dat$X1+dat$gamma_cov}
  
  
  #Treatment model
  if (model_trt=="none"){
    pT<-exp(delta0+a1*dat$X1+a2*dat$X2)/(1+exp(delta0+a1*dat$X1+a2*dat$X2))
  }
  
  if (model_trt=="cl_cov"){
    pT<-exp(delta0+a1*dat$X1+a2*dat$X2+a3*dat$V)/(1+exp(delta0+a1*dat$X1+a2*dat$X2+a3*dat$V))
  }
  
  if (model_trt=="random"){
    pT<-exp(delta0+a1*dat$X1+a2*dat$X2+dat$gamma_T)/(1+exp(delta0+a1*dat$X1+a2*dat$X2+dat$gamma_T))
  }
  
  if (model_trt=="both"){
    pT<-exp(delta0+a1*dat$X1+a2*dat$X2+a3*dat$V+dat$gamma_T)/(1+exp(delta0+a1*dat$X1+a2*dat$X2+a3*dat$V+dat$gamma_T))
  }
  
  dat$Trt<-rbinom(m*k,1,pT)
  
  
  #Outcome model
  if (model_outcome=="cl_cov"){
    pY<-exp(eta0+b1*dat$X1+b2*dat$X2+b3*dat$V)/(1+exp(eta0+b1*dat$X1+b2*dat$X2+b3*dat$V))
  }
  #observations are generated from the cluster-specific probabilities...
  #does this mean we anticipate the ebe ATEs be closer to the true ATE?
  if (model_outcome=="random"){ 
    pY<-exp(eta0+b1*dat$X1+b2*dat$X2+dat$gamma_Y)/(1+exp(eta0+b1*dat$X1+b2*dat$X2+dat$gamma_Y))
  }
  
  if (model_outcome=="both"){
    pY<-exp(eta0+b1*dat$X1+b2*dat$X2+b3*dat$V+dat$gamma_Y)/(1+exp(eta0+b1*dat$X1+b2*dat$X2+b3*dat$V+dat$gamma_Y))
  }
  dat$Y<-rbinom(m*k,1,pY)
  
  return(dat)
  
}

  analysis_bin<-function(tab, R){
    
    t0<-t1<-tab
    t0$Trt<-0
    t1$Trt<-1
    
    
    #Unadjusted-glm
    m1<-glm(Y~Trt, data=tab, family="binomial")
    res1<-t(summary(m1)$coefficient[2,1])
    
    y0<-predict(m1, newdata = t0,
            type = "response")
    y1<-predict(m1, newdata = t1,
                type = "response")
    p1<-mean(y1)
    p0<-mean(y0)
    ATE_m1 <- p1 - p0
    
    marg_m1<-log((p1/(1-p1))/(p0/(1-p0)))
    
    # integrate over the random intercept distribution: for methods 1,3,5 we do not include random effects into the model
    # hence we apply GLMMadaptive only for methods 2,4,6
    # I. via package GLMMadaptive
    # call mixed_model() and then marginal_coefs(); expit of that to get y1 and y0; then mean(y1) - mean(y0)....
    # II. "manually", by using the formula given in the supplementary material of Pavlou et al., 2015 - formula is specific to random intercept models
    # III. approximating the results by the approximate formula for random intercept logistic models given by Zegger et al., 2015.
    
    #Unadjusted-mixed
    m2<-glmer(Y~Trt+(1|ID_clus), data=tab, family="binomial") #AK: consider applying only the GLMMadaptive - to avoid repetition of commands that are equivalent.
    res2<-t(summary(m2)$coefficient[2,1])
    
    y0<-predict(m2, newdata = t0,
                type = "response")
    y1<-predict(m2, newdata = t1,
                type = "response")
    p1<-mean(y1)
    p0<-mean(y0)
    ATE_m2 <- p1 - p0
    marg_m2<-log((p1/(1-p1))/(p0/(1-p0)))
    
    # integrate over the random intercept distribution: GLMMadaptive
    m2 <- mixed_model(fixed = Y ~ as.factor(Trt), random = ~1|ID_clus, data = tab, family = binomial())
    beta.vec <- marginal_coefs(m2)$betas #we obtain the marginal coefficients by averaging over the random effects distribution
    #we apply Monte Carlo integration for the integral approximation
    X_y1 <- as.matrix(cbind(1,t0=1))
    y1 <- as.vector(1/(1 + exp(-(X_y1 %*% beta.vec))))
    X_y0 <- as.matrix(cbind(1,t0=0))
    y0 <- as.vector(1/(1 + exp(-(X_y0 %*% beta.vec))))
    p1<-mean(y1)
    p0<-mean(y0)
    RD_integ2 <- p1 - p0
    pOR_integ2 <-log((p1/(1-p1))/(p0/(1-p0)))
    
    #cs-RD (for random intercept=0)
    beta.vec <- fixef(m2)
    y1 <- as.vector(1/(1 + exp(-(X_y1 %*% beta.vec))))
    y0 <- as.vector(1/(1 + exp(-(X_y0 %*% beta.vec))))
    p1<-mean(y1)
    p0<-mean(y0)
    RD_cs2 <- p1 - p0
    
    #Fully Adjusted-glm
    m3<-glm(Y~Trt+X1+X2+V, data=tab, family="binomial") #to rename X1, X2 as U1, U2 - to maintain correspondence with the main article (if time)
    
    res3b0<-t(summary(m3)$coefficient[1,1]) #intercept
    res3<-t(summary(m3)$coefficient[2,1]) #coefficient of Trt (i.e., gamma)
    res3X1<-t(summary(m3)$coefficient[3,1]) #coefficient estimate of X1
    res3X2<-t(summary(m3)$coefficient[4,1]) #coefficient estimate of X2
    res3V<-t(summary(m3)$coefficient[5,1]) #coefficient estimate of V
    
    y0<-predict(m3, newdata = t0,
                type = "response")
    y1<-predict(m3, newdata = t1,
                type = "response")
    p1<-mean(y1)
    p0<-mean(y0)
    ATE_m3 <- p1 - p0
    marg_m3<-log((p1/(1-p1))/(p0/(1-p0)))
    
    #Fully Adjusted-mixed
    m4<-glmer(Y~Trt+X1+X2+V+(1|ID_clus), data=tab, family="binomial")
    
    res4b0<-t(summary(m4)$coefficient[1,1]) #intercept
    res4<-t(summary(m4)$coefficient[2,1]) #coefficient of Trt
    res4X1<-t(summary(m4)$coefficient[3,1]) #coefficient estimate of X1
    res4X2<-t(summary(m4)$coefficient[4,1]) #coefficient estimate of X2
    res4V<-t(summary(m4)$coefficient[5,1]) #coefficient estimate of V
  
    
    y0<-predict(m4, newdata = t0,
                type = "response")
    y1<-predict(m4, newdata = t1,
                type = "response")
    p1<-mean(y1)
    p0<-mean(y0)
    ATE_m4 <- p1 - p0
    marg_m4<-log((p1/(1-p1))/(p0/(1-p0)))
    
    # integrate over the random intercept distribution: GLMMadaptive
    m4 <- mixed_model(fixed = Y ~ as.factor(Trt) + X1 + as.factor(X2) + V, random = ~1|ID_clus, data = tab, family = binomial())
    beta.vec <- marginal_coefs(m4)$betas #we obtain the marginal coefficients by averaging over the random effects distribution
    #we apply Monte Carlo integration for the integral approximation
    X_y1 <- as.matrix(cbind(1,t0=1,tab[,c("X1","X2","V")]))
    y1 <- as.vector(1/(1 + exp(-(X_y1 %*% beta.vec))))
    X_y0 <- as.matrix(cbind(1,t0=0,tab[,c("X1","X2","V")]))
    y0 <- as.vector(1/(1 + exp(-(X_y0 %*% beta.vec))))
    p1<-mean(y1)
    p0<-mean(y0)
    RD_integ4 <- p1 - p0
    pOR_integ4 <-log((p1/(1-p1))/(p0/(1-p0)))
    
    #cs-RD (for random intercept=0)
    beta.vec <- fixef(m4)
    y1 <- as.vector(1/(1 + exp(-(X_y1 %*% beta.vec))))
    y0 <- as.vector(1/(1 + exp(-(X_y0 %*% beta.vec))))
    p1<-mean(y1)
    p0<-mean(y0)
    RD_cs4 <- p1 - p0
    
    
    
    #Adjusted individual-glm
    m5<-glm(Y~Trt+X1+X2, data=tab, family="binomial")
    
    res5b0<-t(summary(m5)$coefficient[1,1]) #intercept
    res5<-t(summary(m5)$coefficient[2,1]) #coefficient of Trt
    res5X1<-t(summary(m5)$coefficient[3,1]) #coefficient estimate of X1
    res5X2<-t(summary(m5)$coefficient[4,1]) #coefficient estimate of X2
    
    y0<-predict(m5, newdata = t0,
                type = "response")
    y1<-predict(m5, newdata = t1,
                type = "response")
    p1<-mean(y1)
    p0<-mean(y0)
    ATE_m5 <- p1 - p0
    marg_m5<-log((p1/(1-p1))/(p0/(1-p0)))
    
   
    
    #Adjusted individual-mixed
    m6<-glmer(Y~Trt+X1+X2+(1|ID_clus), data=tab, family="binomial")
    
    res6b0<-t(summary(m6)$coefficient[1,1]) #intercept
    res6<-t(summary(m6)$coefficient[2,1]) #coefficient of Trt
    res6X1<-t(summary(m6)$coefficient[3,1]) #coefficient estimate of X1
    res6X2<-t(summary(m6)$coefficient[4,1]) #coefficient estimate of X2
    
    y0<-predict(m6, newdata = t0,
                type = "response")
    y1<-predict(m6, newdata = t1,
                type = "response")
    p1<-mean(y1)
    p0<-mean(y0)
    ATE_m6 <- p1 - p0
    marg_m6<-log((p1/(1-p1))/(p0/(1-p0)))
    
    # integrate over the random intercept distribution: GLMMadaptive
    m6 <- mixed_model(fixed = Y ~ as.factor(Trt) + X1 + as.factor(X2), random = ~1|ID_clus, data = tab, family = binomial())
    beta.vec <- marginal_coefs(m6)$betas #we obtain the marginal coefficients by averaging over the random effects distribution
    #we apply Monte Carlo integration for the integral approximation
    X_y1 <- as.matrix(cbind(1,t0=1,tab[,c("X1","X2")]))
    y1 <- as.vector(1/(1 + exp(-(X_y1 %*% beta.vec))))
    X_y0 <- as.matrix(cbind(1,t0=0,tab[,c("X1","X2")]))
    y0 <- as.vector(1/(1 + exp(-(X_y0 %*% beta.vec))))
    p1<-mean(y1)
    p0<-mean(y0)
    RD_integ6 <- p1 - p0
    pOR_integ6 <-log((p1/(1-p1))/(p0/(1-p0)))
    
    #cs-RD (for random intercept=0)
    beta.vec <- fixef(m6)
    y1 <- as.vector(1/(1 + exp(-(X_y1 %*% beta.vec))))
    y0 <- as.vector(1/(1 + exp(-(X_y0 %*% beta.vec))))
    p1<-mean(y1)
    p0<-mean(y0)
    RD_cs6 <- p1 - p0
    
    
    res<-c(res1,marg_m1,ATE_m1,res2,marg_m2,pOR_integ2,RD_cs2,ATE_m2,RD_integ2,res3b0 - 0.1,res3,res3X1 - 0.2,res3X2 + 1,res3V - 2,marg_m3,
           ATE_m3,res4b0 - 0.1,res4,res4X1 - 0.2,res4X2 + 1,res4V - 2,marg_m4,pOR_integ4,RD_cs4,
           ATE_m4,RD_integ4,res5b0 - 0.1,res5,res5X1 - 0.2,res5X2 + 1,marg_m5,ATE_m5,
           res6b0 - 0.1,res6,res6X1 - 0.2,res6X2 + 1,marg_m6,pOR_integ6,RD_cs6,ATE_m6,RD_integ6) #we print just the bias; we are interested of comparison of the same parameter across methods
    
    #return(res) #AK: to discard tab
    
    #SEs calculation: non-parametric bootstrap ##(for the se3_gamma_boot, to consider obtaining it directly from the model)
    # Non-parametric bootstrap
    
    # a data set for the cluster ids
    vect_clus <- unique(tab$ID_clus) #we need the vector with the cluster ids of the simulated dataset
    data <- data.frame("clusterid"=vect_clus)
    boot_res <- boot(data = tab, statistic = boot_stat, R = r)
    #sd(boot.res.nprcl$t)
    
      se1_gamma_boot <- sd(boot_res$t[,1])
      se1_mar_or_boot <- sd(boot_res$t[,2])
      se1_mar_ate_boot <- sd(boot_res$t[,3])
      
      se2_mixed_bin_csOR_unadj_boot <- sd(boot_res$t[,4])
      se2_mixed_bin_pORebe_boot <- sd(boot_res$t[,5])
      se2_mixed_bin_pORinteg_unadj_boot <- sd(boot_res$t[,6])
      se2_mixed_bin_csRD_unadj_boot <-sd(boot_res$t[,7])
      se2_mixed_bin_RDebe_unadj_boot <-sd(boot_res$t[,8])
      se2_mixed_bin_RDinteg_unadj_boot <-sd(boot_res$t[,9])
    
      se3_glm_b0_adj_boot <- sd(boot_res$t[,10])
      se3_glm_cOR_adj_boot <- sd(boot_res$t[,11])
      se3_glm_b1_adj_boot <- sd(boot_res$t[,12])
      se3_glm_b2_adj_boot <- sd(boot_res$t[,13])
      se3_glm_b3_adj_boot <- sd(boot_res$t[,14])
      se3_glm_mOR_adj_boot <- sd(boot_res$t[,15])  
      se3_glm_RD_adj_boot <- sd(boot_res$t[,16])
      
      se4_mixed_bin_b0_adj_boot <- sd(boot_res$t[,17])
      se4_mixed_bin_csOR_adj_boot <- sd(boot_res$t[,18])
      se4_mixed_bin_b1_adj_boot <- sd(boot_res$t[,19])
      se4_mixed_bin_b2_adj_boot <- sd(boot_res$t[,20])
      se4_mixed_bin_b3_adj_boot <- sd(boot_res$t[,21])
      se4_mixed_bin_pORebe_adj_boot <- sd(boot_res$t[,22])
      se4_mixed_bin_pORinteg_adj_boot <- sd(boot_res$t[,23])
      se4_mixed_bin_csRD_adj_boot <- sd(boot_res$t[,24])
      se4_mixed_bin_RDebe_adj_boot <- sd(boot_res$t[,25])
      se4_mixed_bin_RDinteg_adj_boot <- sd(boot_res$t[,26])
      
      se5_glm_b0_ind_boot <- sd(boot_res$t[,27])
      se5_glm_cOR_ind_boot <- sd(boot_res$t[,28])
      se5_glm_b1_ind_boot <- sd(boot_res$t[,29])
      se_glm_b2_ind_boot <-  sd(boot_res$t[,30])
      se5_glm_mOR_ind_boot <- sd(boot_res$t[,31])
      se5_glm_RD_ind_boot <- sd(boot_res$t[,32])
      
      se6_mixed_bin_b0_ind_boot <- sd(boot_res$t[,33])
      se6_mixed_bin_csOR_ind_boot <- sd(boot_res$t[,34])
      se6_mixed_bin_b1_ind_boot <- sd(boot_res$t[,35])
      se6_mixed_bin_b2_ind_boot <- sd(boot_res$t[,36])  
      se6_mixed_bin_pORebe_ind_boot <- sd(boot_res$t[,37])
      se6_mixed_bin_pORinteg_ind_boot <- sd(boot_res$t[,38])  
      se6_mixed_bin_csRD_ind_boot <- sd(boot_res$t[,39])
      se6_mixed_bin_RDebe_ind_boot <- sd(boot_res$t[,40])
      se6_mixed_bin_RDinteg_ind_boot <- sd(boot_res$t[,41]) 
      
      res <- c(res, se1_gamma_boot,se1_mar_or_boot,se1_mar_ate_boot,se2_mixed_bin_csOR_unadj_boot,se2_mixed_bin_pORebe_boot,se2_mixed_bin_pORinteg_unadj_boot,
               se2_mixed_bin_csRD_unadj_boot, se2_mixed_bin_RDinteg_unadj_boot, se3_glm_b0_adj_boot, se3_glm_cOR_adj_boot, se3_glm_b1_adj_boot,
               se3_glm_b2_adj_boot, se3_glm_b3_adj_boot,se3_glm_mOR_adj_boot, se3_glm_RD_adj_boot,se4_mixed_bin_b0_adj_boot,
               se4_mixed_bin_csOR_adj_boot, se4_mixed_bin_b1_adj_boot, se4_mixed_bin_b2_adj_boot,se4_mixed_bin_b3_adj_boot,se4_mixed_bin_pORebe_adj_boot,
               se4_mixed_bin_pORinteg_adj_boot,se4_mixed_bin_csRD_adj_boot,se4_mixed_bin_RDebe_adj_boot,se5_glm_b0_ind_boot,se5_glm_cOR_ind_boot,se5_glm_b1_ind_boot, 
               se_glm_b2_ind_boot,se5_glm_mOR_ind_boot,se5_glm_RD_ind_boot, se6_mixed_bin_b0_ind_boot,se6_mixed_bin_csOR_ind_boot,se6_mixed_bin_b1_ind_boot,
               se6_mixed_bin_b2_ind_boot,se6_mixed_bin_pORebe_ind_boot,se6_mixed_bin_pORinteg_ind_boot,se6_mixed_bin_csRD_ind_boot,se6_mixed_bin_RDebe_ind_boot,
               se6_mixed_bin_RDinteg_ind_boot)
      
      names(res)<-c("Est_glm_cOR_unadj","Est_glm_mOR_unadj","Est_glm_RD_unadj","Est_mixed_bin_csOR_unadj","Est_mixed_bin_pORebe_unadj",
                    "Est_mixed_bin_pORinteg_unadj","Est_mixed_bin_csRD_unadj","Est_mixed_bin_RDebe_unadj","Est_mixed_bin_RDinteg_unadj",
                    "Est_glm_b0_adj","Est_glm_cOR_adj","Est_glm_b1_adj","Est_glm_b2_adj","Est_glm_b3_adj","Est_glm_mOR_adj","Est_glm_RD_adj",
                    "Est_mixed_bin_b0_adj","Est_mixed_bin_csOR_adj","Est_mixed_bin_b1_adj","Est_mixed_bin_b2_adj","Est_mixed_bin_b3_adj",
                    "Est_mixed_bin_pORebe_adj","Est_mixed_bin_pORinteg_adj","Est_mixed_bin_csRD_adj",
                    "Est_mixed_bin_RDebe_adj","Est_mixed_bin_RDinteg_adj",
                    "Est_glm_b0_ind","Est_glm_cOR_ind","Est_glm_b1_ind","Est_glm_b2_ind","Est_glm_mOR_ind","Est_glm_RD_ind",
                    "Est_mixed_bin_b0_ind","Est_mixed_bin_csOR_ind","Est_mixed_bin_b1_ind","Est_mixed_bin_b2_ind",
                    "Est_mixed_bin_pORebe_ind","Est_mixed_bin_pORinteg_ind","Est_mixed_bin_csRD_ind",
                    "Est_mixed_bin_RDebe_ind","Est_mixed_bin_RDinteg_ind","se1_gamma_boot","se1_mar_or_boot","se1_mar_ate_boot","se2_mixed_bin_csOR_unadj_boot",
                    "n_pORebe_boot","se2_mixed_bin_pORinteg_unadj_boot",
                    "se2_mixed_bin_csRD_unadj_boot", "se2_mixed_bin_RDinteg_unadj_boot", "se3_glm_b0_adj_boot", "se3_glm_cOR_adj_boot", "se3_glm_b1_adj_boot",
                    "se3_glm_b2_adj_boot", "se3_glm_b3_adj_boot","se3_glm_mOR_adj_boot", "se3_glm_RD_adj_boot","se4_mixed_bin_b0_adj_boot",
                    "se4_mixed_bin_csOR_adj_boot", "se4_mixed_bin_b1_adj_boot", "se4_mixed_bin_b2_adj_boot","se4_mixed_bin_b3_adj_boot","se4_mixed_bin_pORebe_adj_boot",
                    "se4_mixed_bin_pORinteg_adj_boot","se4_mixed_bin_csRD_adj_boot","se4_mixed_bin_RDebe_adj_boot","se5_glm_b0_ind_boot","se5_glm_cOR_ind_boot",
                    "se5_glm_b1_ind_boot", 
                    "se_glm_b2_ind_boot","se5_glm_mOR_ind_boot","se5_glm_RD_ind_boot", 
                    "se6_mixed_bin_b0_ind_boot","se6_mixed_bin_csOR_ind_boot","se6_mixed_bin_b1_ind_boot",
                    "se6_mixed_bin_b2_ind_boot","se6_mixed_bin_pORebe_ind_boot","se6_mixed_bin_pORinteg_ind_boot",
                    "se6_mixed_bin_csRD_ind_boot",
                    "se6_mixed_bin_RDebe_ind_boot",
                    "se6_mixed_bin_RDinteg_ind_boot")
      
        }


simulation<-function(nsim,m,k,correlated,model_trt,model_outcome){
  res<-NULL
  
  for (i in 1:nsim){
    print(i)  
    dgen<-generation_bin(m,k,correlated,model_trt,model_outcome)
    res_sim<-analysis_bin(dgen)
    
    res<- rbind(res,res_sim)
  }
  
  summary<-data.frame(t(colMeans(res)))
  summary$m<-m
  summary$k<-k
  summary$corr<-correlated
  summary$model_trt<-model_trt
  summary$model_outcome<-model_outcome
  summary$nsim<-nsim
  
  return(list(res,summary))
  
}

set.seed(20230523)
system.time(s11<-simulation(50,100,200,"uncorr","none","cl_cov"))
system.time(s12<-simulation(1000,20,50,"uncorr","random","both"))
system.time(s21<-simulation(50,50,100,"uncorr","cl_cov","cl_cov"))
system.time(s22<-simulation(1000,20,50,"uncorr","both","both"))
system.time(s31<-simulation(1000,20,50,"corr","cl_cov","cl_cov"))
system.time(s32<-simulation(1000,20,50,"corr","both","both"))

  s1<-simulation(100,100,50,"uncorr","none","cl_cov")
  s2<-simulation(100,100,50,"uncorr","none","random")
  s3<-simulation(100,100,50,"uncorr","none","both")
  s4<-simulation(100,100,50,"uncorr","cl_cov","cl_cov")
  s5<-simulation(100,100,50,"uncorr","cl_cov","random")
  s6<-simulation(100,100,50,"uncorr","cl_cov","both")
  s7<-simulation(100,100,50,"uncorr","random","cl_cov")
  s8<-simulation(100,100,50,"uncorr","random","random")
  s9<-simulation(100,100,50,"uncorr","random","both")
  s10<-simulation(100,100,50,"uncorr","both","cl_cov")
  s11<-simulation(100,100,50,"uncorr","both","random")
  s12<-simulation(100,100,50,"uncorr","both","both")
  
  t1<-simulation(100,100,50,"corr","none","cl_cov")
  t2<-simulation(100,100,50,"corr","none","random")
  t3<-simulation(100,100,50,"corr","none","both")
  t4<-simulation(100,100,50,"corr","cl_cov","cl_cov")
  t5<-simulation(100,100,50,"corr","cl_cov","random")
  t6<-simulation(100,100,50,"corr","cl_cov","both")
  t7<-simulation(100,100,50,"corr","random","cl_cov")
  t8<-simulation(100,100,50,"corr","random","random")
  t9<-simulation(100,100,50,"corr","random","both")
  t10<-simulation(100,100,50,"corr","both","cl_cov")
  t11<-simulation(100,100,50,"corr","both","random")
  t12<-simulation(100,100,50,"corr","both","both")
  
  pooled_AK <-rbind(#s1[[2]],s2[[2]],s3[[2]],s4[[2]],s5[[2]],s6[[2]],s7[[2]],s8[[2]],s9[[2]],s10[[2]],s11[[2]],s12[[2]],
                t1[[2]],t2[[2]],t3[[2]],t4[[2]],t5[[2]],t6[[2]],t7[[2]],t8[[2]],t9[[2]],t10[[2]],t11[[2]],t12[[2]])
  #pooled
  
  write.csv(pooled, "N:/ICON/ICON_all/Users/Andriana Kostouraki/PhD/Simulation (binary outcome)/pooled_AK.csv", row.names=F)

  ##########################################################################################
  #
  #use rsimsum to calculate the bias estimate and its respective Monte Carlo SE estimate
  #
  ##########################################################################################
  
  #let's keep methods 4 and 6 for the moment (i.e., "fully adjusted - GLMM" and "adjusted - individual GLMM", respectively)
  #1. Convert the data set of results into a tidy object
  s11[[1]] <- as.data.frame(s11[[1]])
  data4 <- cbind("rep"=c(1:s11[[2]]$nsim),"dgm"=c("1.1"),"beta0"=c(s11[[1]]$Est_mixed_bin_b0_adj),"method"=c("fully adjusted - GLMM"))
  data6 <- cbind("rep"=c(1:s11[[2]]$nsim),"dgm"=c("1.1"),"beta0"=c(s11[[1]]$Est_mixed_bin_b0_ind),"method"=c("adjusted - individual GLMM"))
  data <- rbind(data4,dat6) #non anticipated error message; to resolve asap
  #2. Use rsimsum to print the Monte Carlo SEs for the relative bias/bias of the theta estimates 
  #(where, theta is: beta0,beta1,beta2,beta3,gamma,cs-ATE,ATE-ebe,ATE-integ)
  #DGM 1.1; methods: "unadjusted - GLM", "unadjusted - GLMM", "fully adjusted - GLM", "fully adjusted - GLMM", "adjusted - individual GLM"
  #"adjusted - individual GLMM"
  