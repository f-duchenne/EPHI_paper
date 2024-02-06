library(bipartite)
library(plot3D)
library(ggplot2)
library(gridExtra)
library(data.table)
library(dplyr)
library(rgbif)
library(ggforce)
library(ggeffects)
library(ggExtra)
library(viridis)
library(lme4)
library(cowplot)
library(scales)
library(car)
library(glmmTMB)
library(ggnewscale)
library(geiger)
library(ggthemes)
library(ggplotify)
library(randomForest)
library(ggtern)
library(ks)
library(sp)
library(RColorBrewer)

filepath="C:/Users/Duchenne/Documents/scripts_divers"
setwd(dir=filepath)

dat=fread("Table2_review.csv")
tab=melt(dat,id.vars=c("model","type","site"))

ggplot(data=tab,aes(y=value,x=site,color=model,shape=variable))+geom_point(size=3,alpha=0.7)+facet_wrap(~type,ncol=2,scales="free_y")+
theme_bw()+theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
legend.position="right",panel.border = element_blank(),axis.line= element_line(),strip.background=element_blank())+
ylab("KS distance between predicted and empirical overlap distributions")+xlab("Sites")+
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
labs(shape="guild")+scale_color_brewer(palette = "Set1")