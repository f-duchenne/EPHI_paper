source("C:/Users/Duchenne/Documents/EPHI_paper/scripts/function_to_predict_from_bayesian.r")
setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper/data")
library(mcmcOutput)
library(mcmcplots)
library(MCMCvis)
library(ggridges)
library(runjags)
library(pastecs)
library(ggplot2)
library(cowplot)
library(gridExtra)

#Costa Rica
load("data_formodel_Costa-Rica.RData")
dat_cr=dat
#### free prior
load("chain_model_free_Costa-Rica_1.RData")
model1=results1
load("chain_model_free_Costa-Rica_2.RData")
model2=results1
load("chain_model_free_Costa-Rica_3.RData")
model3=results1
obj1=as.mcmc(model1)
obj2=as.mcmc(model2)
obj3=as.mcmc(model3)
plot(mcmc.list(obj1,obj3,obj2)[,"Intercept"])
#combine chains and get a summary:
mco_cr <- mcmcOutput(mcmc.list(obj1,obj2,obj3))
suma_cr=summary(mco_cr)
suma_cr$varia=rownames(suma_cr)


#Ecuador
load("data_formodel_Ecuador.RData")
dat_ec=dat
#### free prior
load("chain_model_free_Ecuador_1.RData")
model1=results1
load("chain_model_free_Ecuador_2.RData")
model2=results1
load("chain_model_free_Ecuador_3.RData")
model3=results1
obj1=as.mcmc(model1)
obj2=as.mcmc(model2)
obj3=as.mcmc(model3)
plot(mcmc.list(obj1,obj3,obj2)[,"Intercept"])
#combine chains and get a summary:
mco_ec <- mcmcOutput(mcmc.list(obj1,obj2,obj3))
suma_ec=summary(mco_ec)
suma_ec$varia=rownames(suma_ec)

#Brazil
load("data_formodel_Brazil.RData")
dat_br=dat
#### free prior
load("chain_model_free_Brazil_1.RData")
model1=results1
load("chain_model_free_Brazil_2.RData")
model2=results1
load("chain_model_free_Brazil_3.RData")
model3=results1
obj1=as.mcmc(model1)
obj2=as.mcmc(model2)
obj3=as.mcmc(model3)
plot(mcmc.list(obj1,obj3,obj2)[,"Intercept"])
#combine chains and get a summary:
mco_br <- mcmcOutput(mcmc.list(obj1,obj2,obj3))
suma_br=summary(mco_br)
suma_br$varia=rownames(suma_br)


couleurs=c("#679436","#0B4F6C","#064789")

pre1=predict_model(data=dat_cr,chains=mco_cr,nsampling=5000,site=NA,bird=NA,plant=NA,trait_plant=NA,mismatch=seq(0,60,length.out=100),
barrier=NA,pheno=NA,abond=NA,type="frequency",random_effects=NA,month=NA,year=NA)
pre1$Country="Costa Rica"
pre2=predict_model(data=dat_br,chains=mco_br,nsampling=5000,site=NA,bird=NA,plant=NA,trait_plant=NA,mismatch=seq(0,60,length.out=100),
barrier=NA,pheno=NA,abond=NA,type="frequency",random_effects=NA,month=NA,year=NA)
pre2$Country="Brazil"

pred_mism=rbind(pre1,pre2)

pl1=ggplot(data=pred_mism,aes(x=-mismatch,y=average_lambda,color=Country,fill=Country))+geom_ribbon(aes(ymin=lwr_lambda,ymax=upr_lambda),alpha=0.2,col=NA)+
geom_line(size=1)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),
plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="bottom")+xlab("Trait complementarity")+coord_cartesian(expand=F)+
scale_y_continuous(name = bquote(atop("Average interaction frequency", (interactions.week^-1))),sec.axis = sec_axis(~., name ="" ))+
scale_color_manual(values=couleurs)+scale_fill_manual(values=couleurs)+
labs(col="",fill="")+ggtitle("a")+guides(y="none")


pre1=predict_model(data=dat_cr,chains=mco_cr,nsampling=5000,site=NA,bird=NA,plant=NA,trait_plant=NA,mismatch=NA,barrier=c(0,1),pheno=NA,abond=NA,type="barrier",random_effects=NA,month=NA,year=NA)
pre1$Country="Costa Rica"
pre2=predict_model(data=dat_br,chains=mco_br,nsampling=5000,site=NA,bird=NA,plant=NA,trait_plant=NA,mismatch=NA,barrier=c(0,1),pheno=NA,abond=NA,type="barrier",random_effects=NA,month=NA,year=NA)
pre2$Country="Brazil"
pred_barr=rbind(pre1,pre2)

pl2=ggplot(data=pred_barr,aes(x=as.factor(barrier),y=average_proba,color=Country))+
geom_pointrange(aes(ymin=lwr_proba,ymax=upr_proba),fatten=3)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
plot.title=element_text(size=14,face="bold",hjust = 0),
legend.position="none")+ylab("Probability of interaction")+scale_x_discrete(labels=c("no barrier","barrier"))+xlab("Exploitation barrier")+
labs(col="",fill="")+ggtitle("b")+scale_color_manual(values=couleurs)




#### L2 Match optimum
tabu=as.data.frame(dat_cr[c("hummingbird_numu","culmen_lengthu")])
tabu=cbind(tabu,suma_cr[grep("match_infer",suma_cr$varia),])
tabu$Country="Costa Rica"

tabu2=as.data.frame(dat_br[c("hummingbird_numu","culmen_lengthu")])
tabu2=cbind(tabu2,suma_br[grep("match_infer",suma_br$varia),])
tabu2$Country="Brazil"

pred_L2=rbind(tabu,tabu2)

ribb=data.frame(culmen_lengthu=seq(0.7*min(tabu$culmen_lengthu),1.05*max(tabu$culmen_lengthu),length.out=100))
ribb$ymax=1.5*ribb$culmen_lengthu
ribb$ymin=0.5*ribb$culmen_lengthu

pl3=ggplot(data=pred_L2,aes(x=culmen_lengthu,y=mean,ymax=l95,ymin=u95,color=Country))+
geom_ribbon(data=ribb,aes(y=culmen_lengthu,ymin=ymin,ymax=ymax),fill="white",col=NA)+
geom_abline(intercept=0,slope=1,col="grey")+geom_pointrange(fatten=2,alpha=0.6,size=0.2)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),
plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",panel.background = element_rect(fill = "lightgrey"))+
ylab(bquote(atop("Optimal corolla length",paste("(",italic("L2"),", in mm)"))))+
xlab("Hummingbird bill length (mm)")+
scale_color_manual(values=couleurs)+scale_fill_manual(values=couleurs)+
labs(col="",fill="")+ggtitle("c")+coord_cartesian(expand=F)


#### L1 Barrier threshold
tabu=as.data.frame(dat_cr[c("hummingbird_numu","culmen_lengthu")])
tabu=cbind(tabu,suma_cr[grep("barrier_infer",suma_cr$varia),])
tabu$Country="Costa Rica"

tabu2=as.data.frame(dat_br[c("hummingbird_numu","culmen_lengthu")])
tabu2=cbind(tabu2,suma_br[grep("barrier_infer",suma_br$varia),])
tabu2$Country="Brazil"

pred_L1=rbind(tabu,tabu2)


ribb=data.frame(culmen_lengthu=seq(0.7*min(tabu$culmen_lengthu),1.05*max(tabu$culmen_lengthu),length.out=100))
ribb$ymax=1*ribb$culmen_lengthu
ribb$ymin=2*ribb$culmen_lengthu

pl4=ggplot(data=pred_L1,aes(x=culmen_lengthu,y=mean,ymax=l95,ymin=u95,color=Country))+
geom_ribbon(data=ribb,aes(y=culmen_lengthu,ymin=ymin,ymax=ymax),fill="white",col=NA)+
geom_abline(intercept=0,slope=1,col="grey")+geom_pointrange(fatten=2,alpha=0.6,size=0.2)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_rect(fill = "lightgrey"),
plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none")+ylab(bquote(atop("Threshold for exploitation barrier",paste("(",italic("L1"),", in mm)"))))+
xlab("Hummingbird bill length (mm)")+
scale_color_manual(values=couleurs)+scale_fill_manual(values=couleurs)+
labs(col="",fill="")+ggtitle("d")+coord_cartesian(expand=F)



leg <- ggpubr::as_ggplot(cowplot::get_legend(pl1))
pl1=pl1+theme(legend.position="none")


obj=plot_grid(pl1,pl2,pl3,pl4,align = "hv", rel_heights = c(3, 3),ncol=2)

setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper/")
pdf("Figure2.pdf",width=7,height=6)
grid.arrange(obj,leg,nrow=2,heights=c(7,1))
dev.off();



load("chain_model_binary_1.RData")
model1=results1
load("chain_model_binary_2.RData")
model2=results1
obj1=as.mcmc(model1)
obj2=as.mcmc(model2)
plot(mcmc.list(obj1,obj2)[,"Intercept"])
#combine chains and get a summary:
mco_free <- mcmcOutput(mcmc.list(obj1,obj2))
suma_free=summary(mco_free)
suma_free$varia=rownames(suma_free)
