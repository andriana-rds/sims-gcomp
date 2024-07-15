if (!requireNamespace("ggplot2")) install.packages("ggplot2")
if (!requireNamespace("gridExtra")) install.packages("gridExtra")

library(ggplot2)
library(gridExtra)
setwd("N:/ICON/ICON_all/Users/Andriana Kostouraki/PhD/Latent unconfoundedness work/simulation/results")

df<-read.csv("results_2024_02_29.csv", sep=",", header=T)

tab<-data.frame(aggregate(estimate ~ method + estimand + dgm, data = df, mean))

add_row<-data.frame("method"=rep(c(1,3,5),2), estimand=rep(1,6), dgm=rep(c(1,2), each=3), estimate=0)
add_row2<-data.frame("method"=rep(c(1:6),2), estimand=rep(0,12), dgm=rep(c(1,2), each=6), estimate=0)


tab<-rbind(tab,add_row,add_row2)
tab<-tab[order(tab$dgm, tab$method, tab$estimand),]

t1<-tab[tab$dgm==1,]
t1_GLMM<-t1[t1$method %in% c(2,4,6),]
t1_GLMM$b1<-t1_GLMM$estimate
t1_GLM<-t1[t1$method %in% c(1,3,5),]
t1_GLM$b2<--(t1_GLM$estimate)
t1_GLM$method<-t1_GLM$method+1

t1_ok<-merge(t1_GLM,t1_GLMM, by=c("dgm", "method", "estimand"))
t1_ok$name<-12:1

t1_ok$name <- factor(t1_ok$name, levels = t1_ok$name[order(t1_ok$name)])
t1_ok$name

t2<-tab[tab$dgm==2,]
t2_GLMM<-t2[t2$method %in% c(2,4,6),]
t2_GLMM$b1<-t2_GLMM$estimate
t2_GLM<-t2[t2$method %in% c(1,3,5),]
t2_GLM$b2<--(t2_GLM$estimate)
t2_GLM$method<-t2_GLM$method+1

t2_ok<-merge(t2_GLM,t2_GLMM, by=c("dgm", "method", "estimand"))
t2_ok$name<-12:1

t2_ok$name <- factor(t2_ok$name, levels = t2_ok$name[order(t2_ok$name)])
t2_ok$name


p<-ggplot(t1_ok) +
      geom_segment(aes(x=name, xend=name, y=b1, yend=b2), color=c("white","orange","purple","lightblue",
                                                                   "white","orange","purple","lightblue",
                                                                   "white","orange","purple","lightblue"), linewidth=1) +
      geom_point( aes(x=name, y=b1), color=c("white","orange","purple","lightblue",
                                             "white","orange","purple","lightblue",
                                             "white","orange","purple","lightblue"), size=3 ) +
      geom_point( aes(x=name, y=b2), color=c("white","white","purple","lightblue",
                                             "white", "orange","purple","lightblue",
                                             "white", "white","purple","lightblue"), size=3 ) +
     
       theme(panel.background=element_blank(),legend.title=element_text(size=20),
             plot.title = element_text(size = 16, face = "bold"),
            legend.text=element_text(size=14), 
            axis.line.x = element_line(color = "black",
                                                                        linewidth = 1
                                                                        ),
            axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold")) +
      coord_flip()+
      geom_hline(yintercept = 0, linetype="dotted", 
                   color = "black", linewidth=1)+
    
      scale_y_continuous("Absolute value of bias", breaks=c(-0.2,-0.1,0,0.1,0.2), labels=c("0.2","0.1","0","0.1","0.2"), limits=c(-0.20,0.20))+
      scale_x_discrete("", breaks="", labels="")+
      annotate("text", x = c(4,8,12), y = c(-0.14,-0.185,-0.19), label = c("Adjusted on individual-level covariates","Fully adjusted","Unadjusted"),
               fontface="bold", size=3.4)+
    
    annotate("text", x = c(12.4,12.4), y = c(-0.1,0.1), label = c("GLM","GLMM"),
             fontface="bold")+
      
      geom_segment(x = 2, y = 0.1, xend = 2, yend = 0.11, color="orange", linewidth=1)+
      geom_segment(x = 1.5, y = 0.1, xend = 1.5, yend = 0.11, color="purple", linewidth=1)+
      geom_segment(x = 1, y = 0.1, xend = 1, yend = 0.11, color="lightblue", linewidth=1)+
      annotate("text", x = c(2,1.5,1), y = c(0.127,0.134, 0.137), label = c("cs-ATE","p-ATE-EBE", "p-ATE-integ"),
               fontface="bold", color=c("orange","purple","lightblue"))+
      ggtitle("No correlation between V and U1")

p2<-ggplot(t2_ok) +
  geom_segment( aes(x=name, xend=name, y=b1, yend=b2), color=c("white","orange","purple","lightblue",
                                                               "white","orange","purple","lightblue",
                                                               "white","orange","purple","lightblue"), linewidth=1) +
  geom_point( aes(x=name, y=b1), color=c("white","orange","purple","lightblue",
                                         "white","orange","purple","lightblue",
                                         "white","orange","purple","lightblue"), size=3 ) +
  geom_point( aes(x=name, y=b2), color=c("white","white","purple","lightblue",
                                         "white", "orange","purple","lightblue",
                                         "white", "white","purple","lightblue"), size=3 ) +
  
  theme(panel.background=element_blank(),legend.title=element_text(size=20),
        plot.title = element_text(size = 16, face = "bold"),
        legend.text=element_text(size=14), 
        axis.line.x = element_line(color = "black",
                                   linewidth = 1
        ),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  coord_flip()+
  geom_hline(yintercept = 0, linetype="dotted", 
             color = "black", linewidth=1)+
  
  scale_y_continuous("Absolute value of bias", breaks=c(-0.2,-0.1,0,0.1,0.2), labels=c("0.2","0.1","0","0.1","0.2"), limits=c(-0.20,0.20))+
  scale_x_discrete("", breaks="", labels="")+
  annotate("text", x = c(4,8,12), y = c(-0.14,-0.185,-0.19), label = c("Adjusted on individual-level covariates","Fully adjusted","Unadjusted"),
           fontface="bold", size=3.4)+
  
  annotate("text", x = c(12.25,12.25), y = c(-0.1,0.1), label = c("GLM","GLMM"),
           fontface="bold")+
  
  geom_segment(x = 2, y = 0.1, xend = 2, yend = 0.11, color="orange", linewidth=1)+
  geom_segment(x = 1.5, y = 0.1, xend = 1.5, yend = 0.11, color="purple", linewidth=1)+
  geom_segment(x = 1, y = 0.1, xend = 1, yend = 0.11, color="lightblue", linewidth=1)+
  annotate("text", x = c(2,1.5,1), y = c(0.127,0.134, 0.137), label = c("cs-ATE","p-ATE-EBE", "p-ATE-integ"),
           fontface="bold", color=c("orange","purple","lightblue"))+
  ggtitle("Correlation between V and U1")

tiff(filename = "Figure1.tiff", 
      width = 200, height = 260, units = 'mm', res=300)
grid.arrange(p, p2, nrow = 2)
dev.off()
