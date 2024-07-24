###########################
#                         #
#  Variance estimation    #
#                         #
###########################


# Non-parametric bootstrap

# a data set for the cluster ids
#vect_clus <- unique(tab$ID_clus) #we need the vector with the cluster ids of the simulated dataset
#data <- data.frame("clusterid"=vect_clus)

boot_stat <- function(data,indices){
  
  
  # the data set: resample the clusters
  d <- data[indices, , drop = FALSE] #drop = FALSE: to prevent converting the one-column data.frame into a vector
  # create the full data set: multiple obs per cluster
  d <- merge(tab, d, by = "clusterid") #tab: the simulated (full) data set; d: the bootstrap data set
  # order by cluster id & patient id
  d <- d[order(d$clusterid,d$Id),]
  
  t0<-t1<-d
  t0$Trt<-0
  t1$Trt<-1
  
  # method 1 - glm unadjusted
  m1_boot<-glm(Y~Trt, data=d, family="binomial")
  res1_boot<-t(summary(m1_boot)$coefficient[2,1])
  y0_boot<-predict(m1_boot, newdata = t0,
              type = "response")
  y1_boot<-predict(m1_boot, newdata = t1,
              type = "response")
  p1_boot<-mean(y1_boot)
  p0_boot<-mean(y0_boot)
  ATE_m1_boot <- p1_boot - p0_boot
  marg_m1_boot<-log((p1_boot/(1-p1_boot))/(p0_boot/(1-p0_boot)))
  
  # method 2 - glmm unadjusted
  m2_boot<-glmer(Y~Trt+(1|clusterid), data=d, family="binomial")
  res2_boot<-t(summary(m2_boot)$coefficient[2,1])
  y0_boot<-predict(m2_boot, newdata = t0,
              type = "response")
  y1_boot<-predict(m2_boot, newdata = t1,
              type = "response")
  p1_boot<-mean(y1_boot)
  p0_boot<-mean(y0_boot)
  ATE_m2_boot <- p1_boot - p0_boot
  marg_m2_boot <-log((p1_boot/(1-p1_boot))/(p0_boot/(1-p0_boot)))
  
  # integrate over the random intercept distribution: GLMMadaptive
  m2_boot <- mixed_model(fixed = Y ~ as.factor(Trt), random = ~1|clusterid, data = d, family = binomial())
  beta.vec_boot <- marginal_coefs(m2_boot)$betas #we obtain the marginal coefficients by averaging over the random effects distribution
  #we apply Monte Carlo integration for the integral approximation
  X_y1_boot <- as.matrix(cbind(1,t0=1))
  y1_boot <- as.vector(1/(1 + exp(-(X_y1_boot %*% beta.vec_boot))))
  X_y0_boot <- as.matrix(cbind(1,t0=0))
  y0_boot <- as.vector(1/(1 + exp(-(X_y0_boot %*% beta.vec_boot))))
  p1_boot<-mean(y1_boot)
  p0_boot<-mean(y0_boot)
  RD_integ2_boot <- p1_boot - p0_boot
  pOR_integ2_boot <-log((p1_boot/(1-p1_boot))/(p0_boot/(1-p0_boot)))
  
  #cs-RD (for random intercept=0)
  beta.vec_boot <- fixef(m2_boot)
  y1_boot <- as.vector(1/(1 + exp(-(X_y1_boot %*% beta.vec_boot))))
  y0_boot <- as.vector(1/(1 + exp(-(X_y0_boot %*% beta.vec_boot))))
  p1_boot<-mean(y1_boot)
  p0_boot<-mean(y0_boot)
  RD_cs2_boot <- p1_boot - p0_boot
  
  #method 3 - fully adjusted glm
  m3_boot<-glm(Y~Trt+X1+X2+V, data=d, family="binomial") #to rename X1, X2 as U1, U2 - to maintain correspondence with the main article (if time)
  res3b0_boot<-t(summary(m3_boot)$coefficient[1,1]) #intercept
  res3_boot<-t(summary(m3_boot)$coefficient[2,1]) #coefficient of Trt (i.e., gamma)
  res3X1_boot<-t(summary(m3_boot)$coefficient[3,1]) #coefficient estimate of X1
  res3X2_boot<-t(summary(m3_boot)$coefficient[4,1]) #coefficient estimate of X2
  res3V_boot<-t(summary(m3_boot)$coefficient[5,1]) #coefficient estimate of V
  y0_boot<-predict(m3_boot, newdata = t0,
              type = "response")
  y1_boot<-predict(m3_boot, newdata = t1,
              type = "response")
  p1_boot<-mean(y1_boot)
  p0_boot<-mean(y0_boot)
  ATE_m3_boot <- p1_boot - p0_boot
  marg_m3_boot<-log((p1_boot/(1-p1_boot))/(p0_boot/(1-p0_boot)))
  
  #method 4 - fully adjusted glmm
  m4_boot<-glmer(Y~Trt+X1+X2+V+(1|clusterid), data=d, family="binomial")
  res4b0_boot<-t(summary(m4_boot)$coefficient[1,1]) #intercept
  res4_boot<-t(summary(m4_boot)$coefficient[2,1]) #coefficient of Trt
  res4X1_boot<-t(summary(m4_boot)$coefficient[3,1]) #coefficient estimate of X1
  res4X2_boot<-t(summary(m4_boot)$coefficient[4,1]) #coefficient estimate of X2
  res4V_boot<-t(summary(m4_boot)$coefficient[5,1]) #coefficient estimate of V
  y0_boot<-predict(m4_boot, newdata = t0,
              type = "response")
  y1_boot<-predict(m4_boot, newdata = t1,
              type = "response")
  p1_boot<-mean(y1_boot)
  p0_boot<-mean(y0_boot)
  ATE_m4_boot <- p1_boot - p0_boot
  marg_m4_boot <-log((p1_boot/(1-p1_boot))/(p0_boot/(1-p0_boot)))
  
  # integrate over the random intercept distribution: GLMMadaptive
  m4_boot <- mixed_model(fixed = Y ~ as.factor(Trt) + X1 + as.factor(X2) + V, random = ~1|clusterid, data = d, family = binomial())
  beta.vec_boot <- marginal_coefs(m4_boot)$betas #we obtain the marginal coefficients by averaging over the random effects distribution
  #we apply Monte Carlo integration for the integral approximation
  X_y1_boot <- as.matrix(cbind(1,t0=1,tab[,c("X1","X2","V")]))
  y1_boot <- as.vector(1/(1 + exp(-(X_y1_boot %*% beta.vec_boot))))
  X_y0_boot <- as.matrix(cbind(1,t0=0,tab[,c("X1","X2","V")]))
  y0_boot <- as.vector(1/(1 + exp(-(X_y0_boot %*% beta.vec_boot))))
  p1_boot<-mean(y1_boot)
  p0_boot<-mean(y0_boot)
  RD_integ4_boot <- p1_boot - p0_boot
  pOR_integ4_boot <-log((p1_boot/(1-p1_boot))/(p0_boot/(1-p0_boot)))
  
  #cs-RD (for random intercept=0)
  beta.vec_boot <- fixef(m4_boot)
  y1_boot <- as.vector(1/(1 + exp(-(X_y1_boot %*% beta.vec_boot))))
  y0_boot <- as.vector(1/(1 + exp(-(X_y0_boot %*% beta.vec_boot))))
  p1_boot<-mean(y1_boot)
  p0_boot<-mean(y0_boot)
  RD_cs4_boot <- p1_boot - p0_boot
  
  #method 5 - adjusted, individual glm
  m5_boot<-glm(Y~Trt+X1+X2, data=d, family="binomial")
  
  res5b0_boot<-t(summary(m5_boot)$coefficient[1,1]) #intercept
  res5_boot<-t(summary(m5_boot)$coefficient[2,1]) #coefficient of Trt
  res5X1_boot<-t(summary(m5_boot)$coefficient[3,1]) #coefficient estimate of X1
  res5X2_boot<-t(summary(m5_boot)$coefficient[4,1]) #coefficient estimate of X2
  y0_boot<-predict(m5_boot, newdata = t0,
              type = "response")
  y1_boot<-predict(m5_boot, newdata = t1,
              type = "response")
  p1_boot<-mean(y1_boot)
  p0_boot<-mean(y0_boot)
  ATE_m5_boot <- p1_boot - p0_boot
  marg_m5_boot<-log((p1_boot/(1-p1_boot))/(p0_boot/(1-p0_boot)))
  
  #method 6 - adjusted, individual glmm
  m6_boot<-glmer(Y~Trt+X1+X2+(1|clusterid), data=tab, family="binomial")
  res6b0_boot<-t(summary(m6_boot)$coefficient[1,1]) #intercept
  res6_boot<-t(summary(m6_boot)$coefficient[2,1]) #coefficient of Trt
  res6X1_boot<-t(summary(m6_boot)$coefficient[3,1]) #coefficient estimate of X1
  res6X2_boot<-t(summary(m6_boot)$coefficient[4,1]) #coefficient estimate of X2
  y0_boot<-predict(m6_boot, newdata = t0,
              type = "response")
  y1_boot<-predict(m6_boot, newdata = t1,
              type = "response")
  p1_boot<-mean(y1_boot)
  p0_boot<-mean(y0_boot)
  ATE_m6_boot <- p1_boot - p0_boot
  marg_m6_boot<-log((p1_boot/(1-p1_boot))/(p0_boot/(1-p0_boot)))
  # integrate over the random intercept distribution: GLMMadaptive
  m6_boot <- mixed_model(fixed = Y ~ as.factor(Trt) + X1 + as.factor(X2), random = ~1|clusterid, data = d, family = binomial())
  beta.vec_boot <- marginal_coefs(m6_boot)$betas #we obtain the marginal coefficients by averaging over the random effects distribution
  #we apply Monte Carlo integration for the integral approximation
  X_y1_boot <- as.matrix(cbind(1,t0=1,tab[,c("X1","X2")]))
  y1_boot <- as.vector(1/(1 + exp(-(X_y1_boot %*% beta.vec_boot))))
  X_y0_boot <- as.matrix(cbind(1,t0=0,tab[,c("X1","X2")]))
  y0_boot <- as.vector(1/(1 + exp(-(X_y0_boot %*% beta.vec_boot))))
  p1_boot<-mean(y1_boot)
  p0_boot<-mean(y0_boot)
  RD_integ6_boot <- p1_boot - p0_boot
  pOR_integ6_boot <-log((p1_boot/(1-p1_boot))/(p0_boot/(1-p0_boot)))
  #cs-RD (for random intercept=0)
  beta.vec_boot <- fixef(m6_boot)
  y1_boot <- as.vector(1/(1 + exp(-(X_y1_boot %*% beta.vec_boot))))
  y0_boot <- as.vector(1/(1 + exp(-(X_y0_boot %*% beta.vec_boot))))
  p1_boot<-mean(y1_boot)
  p0_boot<-mean(y0_boot)
  RD_cs6_boot <- p1_boot - p0_boot
  
  
  #return all obtained parameter estimates
  res<-c(res1_boot,marg_m1_boot,ATE_m1_boot,res2_boot,marg_m2_boot,pOR_integ2_boot,RD_cs2_boot,ATE_m2_boot,RD_integ2_boot,res3b0_boot,res3_boot,res3X1_boot,
         res3X2_boot,res3V_boot,marg_m3_boot,
         ATE_m3_boot,res4b0_boot,res4_boot,res4X1_boot,res4X2_boot,res4V_boot,marg_m4_boot,pOR_integ4_boot,RD_cs4_boot,
         ATE_m4_boot,RD_integ4_boot,res5b0_boot,res5_boot,res5X1_boot,res5X2_boot,marg_m5_boot,ATE_m5_boot,
         res6b0_boot,res6_boot,res6X1_boot,res6X2_boot,marg_m6_boot,pOR_integ6_boot,RD_cs6_boot,ATE_m6_boot,RD_integ6_boot)
  return(res)
 
  
}
