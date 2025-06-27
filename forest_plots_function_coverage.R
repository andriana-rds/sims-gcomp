#author: Andriana
#last updated: 27 June 2025

#test forest plot
if (!requireNamespace("forestplot")) install.packages("forestplot")
library(forestplot)
if (!requireNamespace("dplyr")) install.packages("dplyr")
library(dplyr)
if (!requireNamespace("tidyr")) install.packages("tidyr")
library(tidyr)

dat$mean <- dat$`p-ATE EBE fully adjusted GLMM`
dat$lower <- dat$`lower_p-ATE EBE fully adjusted GLMM`
dat$upper <- dat$`upper_p-ATE EBE fully adjusted GLMM`
dat$group <- dat$`coverage_p-ATE EBE fully adjusted GLMM`

dat$mean <- dat$`p-ATE integ fully adjusted GLMM`
dat$lower <- dat$`lower_p-ATE integ fully adjusted GLMM`
dat$upper <- dat$`upper_p-ATE integ fully adjusted GLMM`
dat$group <- dat$`coverage_p-ATE integ fully adjusted GLMM`

dat$mean <- dat$`cs-ATE fully adjusted GLMM`
dat$lower <- dat$`lower_cs-ATE fully adjusted GLMM`
dat$upper <- dat$`upper_cs-ATE fully adjusted GLMM`
dat$group <- dat$`coverage_cs-ATE fully adjusted GLMM`

dat$mean <- dat$`mar-ATE fully adjusted GLM`
dat$lower <- dat$`lower_mar-ATE fully adjusted GLM`
dat$upper <- dat$`upper_mar-ATE fully adjusted GLM`
dat$group <- dat$`coverage_mar-ATE fully adjusted GLM`

dat$mean <- dat$`mar-ATE ind adjusted GLM`
dat$lower <- dat$`lower_mar-ATE ind adjusted GLM`
dat$upper <- dat$`upper_mar-ATE ind adjusted GLM`
dat$group <- dat$`coverage_mar-ATE ind adjusted GLM`
# 
dat$mean <- dat$`cs-ATE ind adjusted GLMM`
dat$lower <- dat$`lower_cs-ATE ind adjusted GLMM`
dat$upper <- dat$`upper_cs-ATE ind adjusted GLMM`
dat$group <- dat$`coverage_cs-ATE ind adjusted GLMM`

dat$mean <- dat$`p-ATE EBE ind adjusted GLMM`
dat$lower <- dat$`lower_p-ATE EBE ind adjusted GLMM`
dat$upper <- dat$`upper_p-ATE EBE ind adjusted GLMM`
dat$group <- dat$`coverage_p-ATE EBE ind adjusted GLMM`

dat$mean <- dat$`p-ATE integ ind adjusted GLMM`
dat$lower <- dat$`lower_p-ATE integ ind adjusted GLMM`
dat$upper <- dat$`upper_p-ATE integ ind adjusted GLMM`
dat$group <- dat$`coverage_p-ATE integ ind adjusted GLMM`

dat$mean <- dat$`mar-ATE unadjusted GLM`
dat$lower <- dat$`lower_mar-ATE unadjusted GLM`
dat$upper <- dat$`upper_mar-ATE unadjusted GLM`
dat$group <- dat$`coverage_mar-ATE unadjusted GLM`

dat$mean <- dat$`cs-ATE unadjusted GLMM`
dat$lower <- dat$`lower_cs-ATE unadjusted GLMM`
dat$upper <- dat$`upper_cs-ATE unadjusted GLMM`
dat$group <- dat$`coverage_cs-ATE unadjusted GLMM`

dat$mean <- dat$`p-ATE EBE unadjusted GLMM`
dat$lower <- dat$`lower_p-ATE EBE unadjusted GLMM`
dat$upper <- dat$`upper_p-ATE EBE unadjusted GLMM`
dat$group <- dat$`coverage_p-ATE EBE unadjusted GLMM`

dat$mean <- dat$`p-ATE integ unadjusted GLMM`
dat$lower <- dat$`lower_p-ATE integ unadjusted GLMM`
dat$upper <- dat$`upper_p-ATE integ unadjusted GLMM`
dat$group <- dat$`coverage_p-ATE integ unadjusted GLMM`
#dat <- as_tibble(dat)
nsims = 500

row_names <- list(rep("",nsims))

dat |>
  forestplot(labeltext = row_names,
             mean = mean,
             lower = lower,
             upper = upper,
             zero = 0,
             cex  = 1,
             lineheight = "auto",
             xlab = "p-ATE integ fully adjusted GLMM and 95% CI: 90% coverage rate") |>
  fp_add_header("") |>
  fp_set_style(lines = gpar(col = "darkblue"))


#alternative - but, will use the snipping tool afterwards, to cut the y-axis label
dat$row_names <- rep(1:nsims)

dat |> group_by(group) |>
  forestplot(legend = c("coverage","no coverage"),
             labeltext = row_names, boxsize = 1.65, 
             
             graphwidth = unit(10, "cm"), #customize graph width
             line.margin = .65, #to avoid crowding
             mean = mean,
             lower = lower,
             upper = upper,
             zero = 0,
             cex  = 1,
             lineheight = "auto",
             xlab = "p-ATE integ unadjusted GLMM and 95% CI: 39% coverage rate") |>
  fp_set_style(lines = c("blue", "red"), box = c("blue", "darkred"))


#note: after having the graph printed, we applied the snipping tool to cut
#the unnecessary numbering of each 95% confidence interval. 
#AK: if time, to find a more efficient way of printing the desired graph, only via the 
#R script.

