#coverage
#last modified: 4th December 2024
#100 simulations; 100 bootstrap replicates
#caution: there must be a coding mistake in the bootstrap function that
#generated these estimates
#author: Andriana 

#read the .RDS file
dt <- readRDS("C:/Users/Andriana/Downloads/test_2024_12_18_cl_boot_500_SIMS_45prcnt_event_rate_rho_0.5.rds")
dat <- as.data.frame(dt[[1]])

dat$`lower_mar-ATE unadjusted GLM`<- dat$`mar-ATE unadjusted GLM` - 1.96*dat$`mar-ATE unadjusted GLM se`
dat$`upper_mar-ATE unadjusted GLM`<- dat$`mar-ATE unadjusted GLM` + 1.96*dat$`mar-ATE unadjusted GLM se`
dat$`coverage_mar-ATE unadjusted GLM`<- ifelse(dat$`lower_mar-ATE unadjusted GLM`<=0&
                  dat$`upper_mar-ATE unadjusted GLM`>=0 & 
                    !is.na(dat$`lower_mar-ATE unadjusted GLM`) &
                    !is.na(dat$`upper_mar-ATE unadjusted GLM`),1,0)  


dat$`lower_cs-ATE unadjusted GLMM`<- dat$`cs-ATE unadjusted GLMM` - 1.96*dat$`cs-ATE unadjusted GLMM se`
dat$`upper_cs-ATE unadjusted GLMM`<- dat$`cs-ATE unadjusted GLMM` + 1.96*dat$`cs-ATE unadjusted GLMM se`
dat$`coverage_cs-ATE unadjusted GLMM`<- ifelse(dat$`lower_cs-ATE unadjusted GLMM`<=0 &
                                                 dat$`upper_cs-ATE unadjusted GLMM`>=0 & 
                                                 !is.na(dat$`lower_cs-ATE unadjusted GLMM`) &
                                                 !is.na(dat$`upper_cs-ATE unadjusted GLMM`),1,0)  

dat$`lower_p-ATE EBE unadjusted GLMM`<- dat$`p-ATE EBE unadjusted GLMM` - 1.96*dat$`p-ATE EBE unadjusted GLMM se`
dat$`upper_p-ATE EBE unadjusted GLMM`<- dat$`p-ATE EBE unadjusted GLMM` + 1.96*dat$`p-ATE EBE unadjusted GLMM se`
dat$`coverage_p-ATE EBE unadjusted GLMM`<- ifelse(dat$`lower_p-ATE EBE unadjusted GLMM`<=0 &
                                                 dat$`upper_p-ATE EBE unadjusted GLMM`>=0 & 
                                                 !is.na(dat$`lower_p-ATE EBE unadjusted GLMM`) &
                                                 !is.na(dat$`upper_p-ATE EBE unadjusted GLMM`),1,0)


dat$`lower_p-ATE integ unadjusted GLMM`<- dat$`p-ATE integ unadjusted GLMM` - 1.96*dat$`p-ATE integ unadjusted GLMM se`
dat$`upper_p-ATE integ unadjusted GLMM`<- dat$`p-ATE integ unadjusted GLMM` + 1.96*dat$`p-ATE integ unadjusted GLMM se`
dat$`coverage_p-ATE integ unadjusted GLMM`<- ifelse(dat$`lower_p-ATE integ unadjusted GLMM`<=0 &
                                                     dat$`upper_p-ATE integ unadjusted GLMM`>=0 & 
                                                     !is.na(dat$`lower_p-ATE integ unadjusted GLMM`) &
                                                     !is.na(dat$`upper_p-ATE integ unadjusted GLMM`),1,0)  



dat$`lower_mar-ATE fully adjusted GLM`<- dat$`mar-ATE fully adjusted GLM` - 1.96*dat$`mar-ATE fully adjusted GLM se`
dat$`upper_mar-ATE fully adjusted GLM`<- dat$`mar-ATE fully adjusted GLM` + 1.96*dat$`mar-ATE fully adjusted GLM se`
dat$`coverage_mar-ATE fully adjusted GLM`<- ifelse(dat$`lower_mar-ATE fully adjusted GLM`<=0 &
                                                    dat$`upper_mar-ATE fully adjusted GLM`>=0 & 
                                                    !is.na(dat$`lower_mar-ATE fully adjusted GLM`) &
                                                    !is.na(dat$`upper_mar-ATE fully adjusted GLM`),1,0)  

dat$`lower_cs-ATE fully adjusted GLMM`<- dat$`cs-ATE fully adjusted GLMM` - 1.96*dat$`cs-ATE fully adjusted GLMM se`
dat$`upper_cs-ATE fully adjusted GLMM`<- dat$`cs-ATE fully adjusted GLMM` + 1.96*dat$`cs-ATE fully adjusted GLMM se`
dat$`coverage_cs-ATE fully adjusted GLMM`<- ifelse(dat$`lower_cs-ATE fully adjusted GLMM`<=0 &
                                                     dat$`upper_cs-ATE fully adjusted GLMM`>=0 & 
                                                     !is.na(dat$`lower_cs-ATE fully adjusted GLMM`) &
                                                     !is.na(dat$`upper_cs-ATE fully adjusted GLMM`),1,0)  


dat$`lower_p-ATE EBE fully adjusted GLMM`<- dat$`p-ATE EBE fully adjusted GLMM` - 1.96*dat$`p-ATE EBE fully adjusted GLMM se`
dat$`upper_p-ATE EBE fully adjusted GLMM`<- dat$`p-ATE EBE fully adjusted GLMM` + 1.96*dat$`p-ATE EBE fully adjusted GLMM se`
dat$`coverage_p-ATE EBE fully adjusted GLMM`<- ifelse(dat$`lower_p-ATE EBE fully adjusted GLMM`<=0 &
                                                     dat$`upper_p-ATE EBE fully adjusted GLMM`>=0 & 
                                                     !is.na(dat$`lower_p-ATE EBE fully adjusted GLMM`) &
                                                     !is.na(dat$`upper_p-ATE EBE fully adjusted GLMM`),1,0)  

dat$`lower_p-ATE integ fully adjusted GLMM`<- dat$`p-ATE integ fully adjusted GLMM` - 1.96*dat$`p-ATE integ fully adjusted GLMM se`
dat$`upper_p-ATE integ fully adjusted GLMM`<- dat$`p-ATE integ fully adjusted GLMM` + 1.96*dat$`p-ATE integ fully adjusted GLMM se`
dat$`coverage_p-ATE integ fully adjusted GLMM`<- ifelse(dat$`lower_p-ATE integ fully adjusted GLMM`<=0 &
                                                        dat$`upper_p-ATE integ fully adjusted GLMM`>=0 & 
                                                        !is.na(dat$`lower_p-ATE integ fully adjusted GLMM`) &
                                                        !is.na(dat$`upper_p-ATE integ fully adjusted GLMM`),1,0)  


dat$`lower_mar-ATE ind adjusted GLM`<- dat$`mar-ATE ind adjusted GLM` - 1.96*dat$`mar-ATE ind adjusted GLM se`
dat$`upper_mar-ATE ind adjusted GLM`<- dat$`mar-ATE ind adjusted GLM` + 1.96*dat$`mar-ATE ind adjusted GLM se`
dat$`coverage_mar-ATE ind adjusted GLM`<- ifelse(dat$`lower_mar-ATE ind adjusted GLM`<=0 &
                                                          dat$`upper_mar-ATE ind adjusted GLM`>=0 & 
                                                          !is.na(dat$`lower_mar-ATE ind adjusted GLM`) &
                                                          !is.na(dat$`upper_mar-ATE ind adjusted GLM`),1,0)  

dat$`lower_cs-ATE ind adjusted GLMM`<- dat$`cs-ATE ind adjusted GLMM` - 1.96*dat$`cs-ATE ind adjusted GLMM se`
dat$`upper_cs-ATE ind adjusted GLMM`<- dat$`cs-ATE ind adjusted GLMM` + 1.96*dat$`cs-ATE ind adjusted GLMM se`
dat$`coverage_cs-ATE ind adjusted GLMM`<- ifelse(dat$`lower_cs-ATE ind adjusted GLMM`<=0 &
                                                   dat$`upper_cs-ATE ind adjusted GLMM`>=0 & 
                                                   !is.na(dat$`lower_cs-ATE ind adjusted GLMM`) &
                                                   !is.na(dat$`upper_cs-ATE ind adjusted GLMM`),1,0)  


dat$`lower_p-ATE EBE ind adjusted GLMM`<- dat$`p-ATE EBE ind adjusted GLMM` - 1.96*dat$`p-ATE EBE ind adjusted GLMM se`
dat$`upper_p-ATE EBE ind adjusted GLMM`<- dat$`p-ATE EBE ind adjusted GLMM` + 1.96*dat$`p-ATE EBE ind adjusted GLMM se`
dat$`coverage_p-ATE EBE ind adjusted GLMM`<- ifelse(dat$`lower_p-ATE EBE ind adjusted GLMM`<=0 &
                                                   dat$`upper_p-ATE EBE ind adjusted GLMM`>=0 & 
                                                   !is.na(dat$`lower_p-ATE EBE ind adjusted GLMM`) &
                                                   !is.na(dat$`upper_p-ATE EBE ind adjusted GLMM`),1,0)  


dat$`lower_p-ATE integ ind adjusted GLMM`<- dat$`p-ATE integ ind adjusted GLMM` - 1.96*dat$`p-ATE integ ind adjusted GLMM se`
dat$`upper_p-ATE integ ind adjusted GLMM`<- dat$`p-ATE integ ind adjusted GLMM` + 1.96*dat$`p-ATE integ ind adjusted GLMM se`
dat$`coverage_p-ATE integ ind adjusted GLMM`<- ifelse(dat$`lower_p-ATE integ ind adjusted GLMM`<=0 &
                                                      dat$`upper_p-ATE integ ind adjusted GLMM`>=0 & 
                                                      !is.na(dat$`lower_p-ATE integ ind adjusted GLMM`) &
                                                      !is.na(dat$`upper_p-ATE integ ind adjusted GLMM`),1,0)  


prop.table(table(dat$`coverage_mar-ATE unadjusted GLM`)) #0 %; 2%
prop.table(table(dat$`coverage_cs-ATE unadjusted GLMM`)) #58%; 42%
prop.table(table(dat$`coverage_p-ATE EBE unadjusted GLMM`)) #50%; 32%
prop.table(table(dat$`coverage_p-ATE integ unadjusted GLMM`)) #40%; 32%

prop.table(table(dat$`coverage_mar-ATE fully adjusted GLM`)) #94%; 90% 
#(AK: to calculate the width of the CI)
prop.table(table(dat$`coverage_cs-ATE fully adjusted GLMM`)) #98%; 88% (undercoverage this time maybe due to low number 
#of bootstrap replicates)
prop.table(table(dat$`coverage_p-ATE EBE fully adjusted GLMM`)) #98%; 86%
prop.table(table(dat$`coverage_p-ATE integ fully adjusted GLMM`)) #98%; 86%

prop.table(table(dat$`coverage_mar-ATE ind adjusted GLM`)) #4%
#First: create a graph - we expect a result like the 2nd graph
#in what Aurelien sent in the Zoom chat
#0% (AK: to calculate the width of the CI - if similar
#the problem is that we have large variance)

prop.table(table(dat$`coverage_cs-ATE ind adjusted GLMM`)) #98%; 90%
prop.table(table(dat$`coverage_p-ATE EBE ind adjusted GLMM`)) #0% - AK to repeat the checks; 4%
prop.table(table(dat$`coverage_p-ATE integ ind adjusted GLMM`)) #98%; 82%





