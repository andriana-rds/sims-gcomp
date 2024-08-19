#s22
set.seed(20230523)
system.time(s22<-simulation(1000,5,1000,"uncorr","both","both"))
setwd("N:/ICON/ICON_all/Users/Andriana Kostouraki/PhD/Latent unconfoundedness work/simulation/results")
write.csv(s22[[2]],"dgm22_1000rep_1000_5.csv",row.names = FALSE)
#write.csv(new_m4_2_2,"new_m4_2_2.csv", row.names = FALSE)
#write.csv(new_m6_2_2,"new_m6_2_2.csv", row.names = FALSE)
#write.csv(s22[[1]],"dgm22_1000rep_100_50_all.csv",row.names = FALSE)
