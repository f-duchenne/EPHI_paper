###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("randomForest","data.table", "dplyr", "lubridate","lme4","R2jags","mcmcOutput","mcmcplots","MCMCvis",
			"pastecs","ggplot2","cowplot","gridExtra","scales") 


inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)

EPHI_version="2024-03-18"

setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper/data")
source("C:/Users/Duchenne/Documents/EPHI_paper/scripts/additional scripts/function_to_predict_from_bayesian.r")


#Costa Rica
load("data_formodel_Costa-Rica.RData")
dat_cr=dat
#### free prior
load("chain_model_Costa-Rica_1.RData")
model1=results1
load("chain_model_Costa-Rica_2.RData")
model2=results1
load("chain_model_Costa-Rica_3.RData")
model3=results1
obj1=as.mcmc(model1)
obj2=as.mcmc(model2)
obj3=as.mcmc(model3)
plot(mcmc.list(obj1,obj3,obj2)[,"Intercept"])
plot(mcmc.list(obj1,obj3,obj2)[,"traitBarrier"])
plot(mcmc.list(obj1,obj3,obj2)[,"traitMismatch"])
#combine chains and get a summary:
mco_cr <- mcmcOutput(mcmc.list(obj1,obj2,obj3))
suma_cr=summary(mco_cr)
suma_cr$varia=rownames(suma_cr)


#Ecuador
load("data_formodel_Ecuador.RData")
dat_ec=dat
#### free prior
load("chain_model_Ecuador_1.RData")
model1=results1
load("chain_model_Ecuador_2.RData")
model2=results1
load("chain_model_Ecuador_3.RData")
model3=results1
obj1=as.mcmc(model1)
obj2=as.mcmc(model2)
obj3=as.mcmc(model3)
plot(mcmc.list(obj1,obj2,obj3)[,"Intercept"])
#combine chains and get a summary:
mco_ec <- mcmcOutput(mcmc.list(obj1,obj2,obj3))
suma_ec=summary(mco_ec)
suma_ec$varia=rownames(suma_ec)

#Brazil
load("data_formodel_Brazil.RData")
dat_br=dat
#### free prior
load("chain_model_Brazil_1.RData")
model1=results1
load("chain_model_Brazil_2.RData")
model2=results1
load("chain_model_Brazil_3.RData")
model3=results1
obj1=as.mcmc(model1)
obj2=as.mcmc(model2)
obj3=as.mcmc(model3)
plot(mcmc.list(obj1,obj3,obj2)[,"Intercept"])
#combine chains and get a summary:
mco_br <- mcmcOutput(mcmc.list(obj1,obj2,obj3))
suma_br=summary(mco_br)
suma_br$varia=rownames(suma_br)


couleurs=c("#679436","#0B4F6C","deeppink")
##################### FIGURE 2 ##################

##################### Panel a ##################
pre1=predict_model(data=dat_cr,chains=mco_cr,nsampling=5000,site=NA,bird=NA,plant=NA,trait_plant=NA,mismatch=seq(0,60,length.out=100),
barrier=NA,pheno=NA,abond=NA,type="frequency",random_effects=NA,month=NA,year=NA,duration=NA)
pre1$Country="Costa Rica"
pre2=predict_model(data=dat_br,chains=mco_br,nsampling=5000,site=NA,bird=NA,plant=NA,trait_plant=NA,mismatch=seq(0,60,length.out=100),
barrier=NA,pheno=NA,abond=NA,type="frequency",random_effects=NA,month=NA,year=NA,duration=NA)
pre2$Country="Brazil"
pre3=predict_model(data=dat_ec,chains=mco_ec,nsampling=5000,site=NA,bird=NA,plant=NA,trait_plant=NA,mismatch=seq(0,60,length.out=100),
barrier=NA,pheno=NA,abond=NA,type="frequency",random_effects=NA,month=NA,year=NA,duration=NA)
pre3$Country="Ecuador"


pred_mism=rbind(pre1,pre2,pre3)

pl1=ggplot(data=pred_mism,aes(x=-mismatch,y=average_lambda,color=Country,fill=Country))+geom_ribbon(aes(ymin=lwr_lambda,ymax=upr_lambda),alpha=0.2,col=NA)+
geom_line(size=1)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),
plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none")+xlab("Trait complementarity")+coord_cartesian(expand=F)+
scale_y_continuous(name = bquote(atop("Average interaction frequency", (interactions.day^-1))),sec.axis = sec_axis(~., name ="" ))+
scale_color_manual(values=couleurs)+scale_fill_manual(values=couleurs)+
labs(col="",fill="")+ggtitle("a")+guides(y="none")

##################### Panel b ##################
pre1=predict_model(data=dat_cr,chains=mco_cr,nsampling=5000,site=NA,bird=NA,plant=NA,trait_plant=NA,mismatch=NA,barrier=c(0,1),pheno=NA,abond=NA,
type="barrier",random_effects=NA,month=NA,year=NA,duration=NA)
pre1$Country="Costa Rica"
pre2=predict_model(data=dat_br,chains=mco_br,nsampling=5000,site=NA,bird=NA,plant=NA,trait_plant=NA,mismatch=NA,barrier=c(0,1),pheno=NA,abond=NA,
type="barrier",random_effects=NA,month=NA,year=NA,duration=NA)
pre2$Country="Brazil"
pre3=predict_model(data=dat_ec,chains=mco_ec,nsampling=5000,site=NA,bird=NA,plant=NA,trait_plant=NA,mismatch=NA,barrier=c(0,1),pheno=NA,abond=NA,
type="barrier",random_effects=NA,month=NA,year=NA,duration=NA)
pre3$Country="Ecuador"
pred_barr=rbind(pre1,pre2,pre3)

pl2=ggplot(data=pred_barr,aes(x=as.factor(barrier),y=average_proba,color=Country))+
geom_pointrange(aes(ymin=lwr_proba,ymax=upr_proba),fatten=3)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
plot.title=element_text(size=14,face="bold",hjust = 0))+ylab("Probability of interaction")+scale_x_discrete(labels=c("no barrier","barrier"))+xlab("Exploitation barrier")+
labs(col="",fill="")+ggtitle("b")+scale_color_manual(values=couleurs)


##################### Panel c ##################
setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper/data")
res_partition=fread("sparsity_estimates.txt")
forplot=melt(res_partition,id.vars=c("site","min_transect_elev","duration","Country","tot","pzerom"))

couleurs2=c("#477998","#291F1E","#F64740")

forplot=forplot[order(forplot$Country,forplot$min_transect_elev),]
forplot$site=factor(forplot$site,levels=unique(forplot$site))
forplot$variable=as.character(forplot$variable)
forplot$variable[forplot$variable=="barrier"]="exploitation barrier"
forplot$variable=factor(forplot$variable,levels=c("unknown","exploitation barrier"))

summary(subset(res_partition,duration %in% c(2000))$barrier)

pl3=ggplot(data=subset(forplot,duration %in% c(2000) & variable!="complementarity"),aes(y=value,fill=variable,x=site))+geom_bar(stat="identity")+
facet_grid(col=vars(Country),scales = "free_x",shrink=T,drop=TRUE,space ="free_x")+
ylab("% of invariant zeros in networks")+
theme_bw()+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background = element_blank(),
axis.text.x=element_text(color="white"),legend.position="bottom")+
labs(fill="Generated by:")+xlab("Sites")+
scale_color_manual(values=couleurs2)+scale_fill_manual(values=couleurs2)+
scale_y_continuous(labels = scales::percent)+ggtitle("c")+coord_cartesian(expand=F)

##### PLOT
leg <- ggpubr::as_ggplot(cowplot::get_legend(pl2))
leg=leg+theme(legend.margin=margin(c(0,0,0,0)))
pl2=pl2+theme(legend.position="none")

top=plot_grid(pl1,pl2,leg,align = "hv",ncol=3,rel_widths=c(1,1,0.1))
bottom=pl3
plot_grid(top,bottom,align = "hv",ncol=1,rel_heights=c(1,1.5))

setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper/")
pdf("Figure2.pdf",width=7,height=6)
plot_grid(top,bottom,align = "hv",ncol=1,rel_heights=c(1,1.5))
dev.off();


##################### Possible extra panel ##################
setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper/data")
res_partition=fread("sparsity_estimates.txt")

pl4=ggplot(data=res_partition,aes(x=duration,color=Country,y=pzerom))+geom_point()+
ylab("% of zeros (non-interactions)")+
theme_bw()+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background = element_blank(),legend.position="bottom")+
labs(col="")+xlab("Duration of sampling per camera (hours)")+
scale_color_manual(values=couleurs)+scale_fill_manual(values=couleurs)+
scale_y_continuous(labels = scales::percent)+ggtitle("c")


##################### FIGURE S3 ##################
pl=ggplot(data=subset(forplot,duration %in% c(100,1000)),aes(y=value,fill=variable,x=site))+geom_bar(stat="identity")+
facet_grid(duration~Country,scales = "free_x",shrink=T,drop=TRUE,space ="free_x",
labeller = label_bquote(rows=paste("Duration of sampling: ",.(duration),"h")))+
ylab("% of invariant zero generated")+
theme_bw()+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background = element_blank(),
axis.text.x=element_text(angle=45,hjust=1))+
labs(fill="Linkage rule")+xlab("Sites")+
scale_color_manual(values=couleurs2)+scale_fill_manual(values=couleurs2)+
scale_y_continuous(labels = scales::percent)

setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper/")
png("Figure_S3.png",width=1000,height=800,res=120)
pl
dev.off();

####################### FIGURE S2 ##################
#### L2 Match optimum
tabu=as.data.frame(dat_cr[c("hummingbird_numu","culmen_lengthu")])
tabu=cbind(tabu,suma_cr[grep("match_infer",suma_cr$varia),])
tabu$Country="Costa Rica"

tabu2=as.data.frame(dat_br[c("hummingbird_numu","culmen_lengthu")])
tabu2=cbind(tabu2,suma_br[grep("match_infer",suma_br$varia),])
tabu2$Country="Brazil"

tabu3=as.data.frame(dat_ec[c("hummingbird_numu","culmen_lengthu")])
tabu3=cbind(tabu3,suma_ec[grep("match_infer",suma_ec$varia),])
tabu3$Country="Ecuador"

pred_L2=rbind(tabu,tabu2,tabu3)

nrow(subset(pred_L2,l95>culmen_lengthu | u95<culmen_lengthu))/nrow(pred_L2)

ribb=data.frame(culmen_lengthu=seq(0.7*min(pred_L2$culmen_lengthu),1.05*max(pred_L2$culmen_lengthu),length.out=100))
ribb$ymax=1.5*ribb$culmen_lengthu
ribb$ymin=0.5*ribb$culmen_lengthu

pl5=ggplot(data=pred_L2,aes(x=culmen_lengthu,y=mean,ymax=l95,ymin=u95,color=Country))+
geom_ribbon(data=ribb,aes(y=culmen_lengthu,ymin=ymin,ymax=ymax),fill="white",col=NA)+
geom_abline(intercept=0,slope=1,col="black")+geom_pointrange(fatten=3,alpha=0.6,size=0.01)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),
plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",panel.background = element_rect(fill = "lightgrey"))+
ylab(bquote(atop("Optimal corolla length",paste("(",italic("L2"),", in mm)"))))+
xlab("Hummingbird bill length (mm)")+
scale_color_manual(values=couleurs)+scale_fill_manual(values=couleurs)+
labs(col="",fill="")+ggtitle("a")+coord_cartesian(expand=F)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x))) +
     scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)))+annotation_logticks()


#### L1 Barrier threshold
tabu=as.data.frame(dat_cr[c("hummingbird_numu","culmen_lengthu")])
tabu=cbind(tabu,suma_cr[grep("barrier_infer",suma_cr$varia),])
tabu$Country="Costa Rica"

tabu2=as.data.frame(dat_br[c("hummingbird_numu","culmen_lengthu")])
tabu2=cbind(tabu2,suma_br[grep("barrier_infer",suma_br$varia),])
tabu2$Country="Brazil"

tabu3=as.data.frame(dat_ec[c("hummingbird_numu","culmen_lengthu")])
tabu3=cbind(tabu3,suma_ec[grep("barrier_infer",suma_ec$varia),])
tabu3$Country="Ecuador"

pred_L1=rbind(tabu,tabu2,tabu3)

nrow(subset(pred_L1,l95>(culmen_lengthu+0.05*culmen_lengthu)))/nrow(pred_L1)

ribb=data.frame(culmen_lengthu=seq(0.7*min(pred_L1$culmen_lengthu),1.05*max(pred_L1$culmen_lengthu),length.out=100))
ribb$ymax=1*ribb$culmen_lengthu
ribb$ymin=2*ribb$culmen_lengthu

pl6=ggplot(data=pred_L1,aes(x=culmen_lengthu,y=mean,ymax=l95,ymin=u95,color=Country))+
geom_ribbon(data=ribb,aes(y=culmen_lengthu,ymin=ymin,ymax=ymax),fill="white",col=NA)+
geom_abline(intercept=0,slope=1,col="black")+geom_pointrange(fatten=3,alpha=0.6,size=0.01)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_rect(fill = "lightgrey"),
plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none")+
ylab(bquote(atop("Threshold for exploitation barrier",paste("(",italic("L1"),", in mm)"))))+
xlab("Hummingbird bill length (mm)")+
scale_color_manual(values=couleurs)+scale_fill_manual(values=couleurs)+
labs(col="",fill="")+ggtitle("b")+coord_cartesian(expand=F)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x))) +
     scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)))+annotation_logticks()

setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper/")
png("Figure_S2.png",width=2000,height=1000,res=300)
plot_grid(pl5,pl6,align = "hv",ncol=2)
dev.off();

#Export table for latter
names(pred_L1)[3:10]=paste0("barrier_",names(pred_L1)[3:10])
names(pred_L2)[3:10]=paste0("match_",names(pred_L2)[3:10])
tabu=unique(as.data.frame(dat_cr[c("hummingbird_num","hummingbird_species")]))
tabu$Country="Costa Rica"
tabu2=unique(as.data.frame(dat_br[c("hummingbird_num","hummingbird_species")]))
tabu2$Country="Brazil"
tabu3=unique(as.data.frame(dat_ec[c("hummingbird_num","hummingbird_species")]))
tabu3$Country="Ecuador"
table=rbind(tabu,tabu2,tabu3)
names(table)[1]="hummingbird_numu"

humm_table=merge(pred_L1,pred_L2,by=c("hummingbird_numu","Country"))
humm_table=merge(humm_table,table,by=c("hummingbird_numu","Country"))

setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper/data/")
data.table::fwrite(humm_table,"humm_table.txt")


















