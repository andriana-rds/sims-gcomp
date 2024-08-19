#Bias plots: print two lolliplots for DGM 2.2 and 3.2 respectively
#setwd("N:/ICON/ICON_all/Users/Andriana Kostouraki/PhD/Latent unconfoundedness work/simulation/results")
setwd("C:/Users/AndrianaKostouraki/Downloads/marginal_estimands_to_resubmission (2)/marginal_estimands_to_resubmission/tentative/results")
n_sim <- 1000 #number of iterations
n_methods <- 6
n_dgm <- 2 

#Target: the cluster-specific coefficient estimates for treatment 
#results <- read.csv("dgm22_1000rep_100_50_all.csv")
#res <- read.csv("dgm32_1000rep_05_100_50_all.csv")

results <- readRDS("dgm_32_1000rep_08_1000_5.rds")
results <- as.data.frame(results[[1]])

res <- readRDS("dgm_22_1000rep_1000_5.rds")
res <- as.data.frame(res[[1]])

results$rep <- res$rep <- 1:n_sim
results$dgm <- 1 #DGM 2.2
res$dgm <- 2 #DGM 3.2
results <- rbind(results,res)

#now, create a method and an estimand column
#results1 <- results2 <- results3 <- results4 <- results5 <- results6 <- NA
#for cs-ATE (RD): coded as estimand 1
#unadjusted-GLM
results1 <- as.data.frame(cbind(estimand=1,estimate=NA,method=1,dgm=results$dgm,rep=results$rep))

#results1$estimand <- 1
#results1$estimate <- "not available"
#results1$method <- 1
#results1$dgm <- results$dgm

#unadjusted-GLMM
results2 <- as.data.frame(cbind(estimand=1, estimate=results$Est_mixed_bin_csRD_unadj,method=2,dgm=results$dgm,rep=results$rep))
#results2$estimand <- 1
#results2$estimate <- results$Est_mixed_bin_csRD_unadj
#results2$method <- 2
#results2$dgm <- results$dgm

#fully adjusted - GLM
results3 <- as.data.frame(cbind(estimand=1,estimate=NA,method=3,dgm=results$dgm,rep=results$rep))
#results3$estimand <- 1
#results3$estimate <- "not available"
#results3$method <- 3
#results3$dgm <- results$dgm

#fully adjusted - GLMM
results4 <- as.data.frame(cbind(estimand=1,estimate=results$Est_mixed_bin_csRD_adj,method=4,dgm=results$dgm,rep=results$rep))
#results4$estimand <- 1
#results4$estimate <- results$Est_mixed_bin_csRD_adj
#results4$method <- 4
#results4$dgm <- results$dgm

#adjusted - individidual GLM
results5 <- as.data.frame(cbind(estimand=1,estimate=NA,method=5,dgm=results$dgm,rep=results$rep))
#results5$estimand <- 1
#results5$estimate <- "not available"
#results5$method <- 5
#results5$dgm <- results$dgm

#adjusted - individual GLMM
results6 <- as.data.frame(cbind(estimand=1,estimate=results$Est_mixed_bin_csRD_ind,method=6,dgm=results$dgm,rep=results$rep))
#results6$estimand <- 1
#results6$estimate <- results$Est_mixed_bin_csRD_ind
#results6$method <- 6
#results6$dgm <- results$dgm
results_cs_ATE <- rbind(results2,results4,results6)

#for ATE-ebe (RD): coded as estimand 2
#unadjusted-GLM
results1 <- as.data.frame(cbind(estimand=2, estimate=results$Est_glm_RD_unadj,method=1,dgm=results$dgm,rep=results$rep))
#results1$estimand <- 2
#results1$estimate <- results$Est_glm_RD_unadj
#results1$method <- 1
#results1$dgm <- results$dgm

#unadjusted-GLMM
results2 <- as.data.frame(cbind(estimand=2, estimate=results$Est_mixed_bin_RDebe_unadj,method=2,dgm=results$dgm,rep=results$rep))

#results2$estimand <- 2
#results2$estimate <- results$Est_mixed_bin_RDebe_unadj
#results2$method <- 2
#results2$dgm <- results$dgm

#fully adjusted - GLM
results3 <- as.data.frame(cbind(estimand=2, estimate=results$Est_glm_RD_adj,method=3,dgm=results$dgm,rep=results$rep))
#results3$estimand <- 2
#results3$estimate <- results$Est_glm_RD_adj
#results3$method <- 3
#results3$dgm <- results$dgm

#fully adjusted - GLMM
results4 <- as.data.frame(cbind(estimand=2, estimate=results$Est_mixed_bin_RDebe_adj,method=4,dgm=results$dgm,rep=results$rep))

#results4$estimand <- 2
#results4$estimate <- results$Est_mixed_bin_RDebe_adj
#results4$method <- 4
#results4$dgm <- results$dgm

#adjusted - individidual GLM
results5 <- as.data.frame(cbind(estimand=2, estimate=results$Est_glm_RD_ind,method=5,dgm=results$dgm,rep=results$rep))
#results5$estimand <- 2
#results5$estimate <- results$Est_glm_RD_ind
#results5$method <- 5
#results5$dgm <- results$dgm

#adjusted - individual GLMM
results6 <- as.data.frame(cbind(estimand=2, estimate=results$Est_mixed_bin_RDebe_ind,method=6,dgm=results$dgm,rep=results$rep))
#results6$estimand <- 2
#results6$estimate <- results$Est_mixed_bin_RDebe_ind
#results6$method <- 6
#results6$dgm <- results$dgm
results_ATE_ebe <- rbind(results1,results2,results3,results4,results5,results6)


#for ATE-integ (RD): coded as estimand 3
#unadjusted-GLM
results1 <- as.data.frame(cbind(estimand=3, estimate=results$Est_glm_RD_unadj,method=1,dgm=results$dgm,rep=results$rep))
#results1$estimand <- 3
#results1$estimate <- results$Est_glm_RD_unadj
#results1$method <- 1
#results1$dgm <- results$dgm
#unadjusted-GLMM
results2 <- as.data.frame(cbind(estimand=3, estimate=results$Est_mixed_bin_RDinteg_unadj,method=2,dgm=results$dgm,rep=results$rep))

#results2$estimand <- 3
#results2$estimate <- results$Est_mixed_bin_RDinteg_unadj
#results2$method <- 2
#results2$dgm <- results$dgm

#fully adjusted - GLM
results3 <- as.data.frame(cbind(estimand=3, estimate=results$Est_glm_RD_adj,method=3,dgm=results$dgm, rep=results$rep))

#results3$estimand <- 3
#results3$estimate <- results$Est_glm_RD_adj
#results3$method <- 3
#results3$dgm <- results$dgm

#fully adjusted - GLMM
results4 <- as.data.frame(cbind(estimand=3, estimate=results$Est_mixed_bin_RDinteg_adj,method=4,dgm=results$dgm, rep=results$rep))

#results4$estimand <- 3
#results4$estimate <- results$Est_mixed_bin_RDinteg_adj
#results4$method <- 4
#results4$dgm <- results$dgm

#adjusted - individidual GLM
results5 <- as.data.frame(cbind(estimand=3, estimate=results$Est_glm_RD_ind,method=5,dgm=results$dgm,rep=results$rep))

#results5$estimand <- 3
#results5$estimate <- results$Est_glm_RD_ind
#results5$method <- 5
#results5$dgm <- results$dgm

#adjusted - individual GLMM
results6 <- as.data.frame(cbind(estimand=3, estimate=results$Est_mixed_bin_RDinteg_ind,method=6,dgm=results$dgm,rep=results$rep))

#results6$estimand <- 3
#results6$estimate <- results$Est_mixed_bin_RDinteg_ind
#results6$method <- 6
#results6$dgm <- results$dgm
results_ATE_integ <- rbind(results1,results2,results3,results4,results5,results6)
results <- rbind(results_cs_ATE,results_ATE_ebe,results_ATE_integ)
