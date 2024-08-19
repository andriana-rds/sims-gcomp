# s32: rho=0.3
set.seed(20230523)
system.time(s32<-simulation(1000,5,1000,"corr","both","both"))
setwd("N:/ICON/ICON_all/Users/Andriana Kostouraki/PhD/Latent unconfoundedness work/simulation/results")
write.csv(s32[[2]],"dgm32_1000rep_03_1000_5.csv",row.names = FALSE)