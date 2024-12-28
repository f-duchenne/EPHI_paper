#######################################################################################################################
###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
pkgs <- c("randomForest","data.table", "dplyr", "lubridate","lme4","qgraph","igraph","ggeffects","DHARMa",
			"pastecs","ggplot2","cowplot","gridExtra","scales","reshape2","bipartite","stringr","mgcv","ggpubr","emmeans","piecewiseSEM","car","glmmTMB","here") 


# Check if packages are already installed
inst <- pkgs %in% installed.packages()
# Install missing packages
if (any(!inst)) install.packages(pkgs[!inst])
# Load packages
pkg_out <- lapply(pkgs, require, character.only = TRUE)

# Set the working directory for the "here" package
here::i_am("EPHI_paper.Rproj")

#color vector for figures
couleurs=c("#679436","#0B4F6C","deeppink")

# load robustness results
extinctions=fread(here("data_zenodo","robustness_simulations_all.csv"))

###################

#### CALCULATE PERSISTENCE FOR EACH REMAINING FRACTION OF PLANT SPECIES STEP (ATAC CURVE)
b=extinctions %>% group_by(r,site,Country,site2,barrier,rank2,min_transect_elev,nbh,nbp,scenario,essai,prop_forbidden_m,C_m,N_m,Inter_eve_m,Cperso_m) %>% summarise(persistence=mean(hummingbird_pers/nbh))
b$symmetrie=b$nbp/b$nbh
b2=extinctions %>% group_by(r,site,Country,site2,barrier,rank2,min_transect_elev,nbh,nbp,scenario,prop_forbidden_m,C_m,N_m,Inter_eve_m,Cperso_m) %>% summarise(persistence_mean=mean(hummingbird_pers/nbh),lower=quantile(hummingbird_pers/nbh,probs=0.025),
higher=quantile(hummingbird_pers/nbh,probs=0.975))

pays="Brazil"

pl1=ggplot()+
geom_line(data=subset(b,scenario=="generalists first" & r==1 & Country==pays),
aes(x=rank2,y=persistence,col=barrier,group=paste0(essai,r,site,barrier,scenario)),alpha=0.01)+
geom_line(data=subset(b2,scenario=="generalists first" & r==1 & Country==pays),aes(y=persistence_mean,x=rank2,col=barrier),size=1.5)+facet_wrap(~site2)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none")+
labs(col="")+ggtitle("a")+ylab("Persistence")+xlab("Fraction of plant species removed")+
scale_color_manual(values=c("#30343F","red"))+scale_fill_manual(values=c("#30343F","red"))+
scale_x_continuous(n.breaks=2,labels = scales::number_format(accuracy = 1))


### AREA UNDER ATAC CURVE (ROBUSTNESS)
b3=subset(b,rank2!=0) %>% group_by(site,site2,Country,barrier,r,min_transect_elev,nbh,nbp,scenario,prop_forbidden_m,C_m,N_m,Inter_eve_m,Cperso_m,symmetrie,essai) %>%
summarise(robustness=sum(persistence)/length(persistence))

# prepare variables for statistical models
b3$barrierbis=0
b3$barrierbis[b3$barrier=="with forbidden links"]=1
b3$barrierbis=as.factor(b3$barrierbis)
b3$sitebis=as.factor(b3$site)
b3$alpha=as.factor(b3$r)

### FIT STATISTICAL MODELS
#when generalists get extinct first
model=glmmTMB(robustness~barrierbis*sitebis*alpha,data=subset(b3,scenario=="generalists first"),family=beta_family(),control=glmmTMBControl(optCtrl=list(iter.max=1e5,eval.max=1e5)))
#when specialists get extinct first
model2=glmmTMB(robustness~barrierbis*sitebis*alpha,data=subset(b3,scenario=="specialists first"),family=beta_family(),control=glmmTMBControl(optCtrl=list(iter.max=1e5,eval.max=1e5)))

# predict for the scenario when  generalists get extinct first
pred=ggpredict(model,c("barrierbis","sitebis","alpha"))
pred=merge(pred,unique(b3[,c("site","Country","prop_forbidden_m","barrierbis")]),by.x=c("group","x"),by.y=c("site","barrierbis"))

# calculate effect of epxloitation barrier on robustness for each site 
model.emm <- as.data.frame(contrast(model.emm, "trt.vs.ctrl", ref = "barrierbis0"))
#merge with site level predictions
pred=merge(pred,model.emm,by.x=c("group","facet"),by.y=c("sitebis","alpha"))
#calculate significancy of effects
pred$signifi="significant"
pred$signifi[pred$p.value>0.05]="not significant"
pred$signifi=factor(pred$signifi,levels=c("significant","not significant"))
pred$scenario="generalists first"

#focus on effects for exploitation barrier only
pred2=subset(pred,facet==1)
dim(pred2[pred2$p.value<0.05 & pred2$odds.ratio<1,])
dim(pred2[pred2$p.value<0.05 & pred2$odds.ratio>1,])
dim(pred2[pred2$p.value>0.05,])

pl2=ggplot(data=pred2,aes(x=x,y=predicted))+
geom_boxplot(aes(color=x))+
geom_line(aes(linetype=signifi,group=group),show.legend = F)+
geom_pointrange(aes(ymin=conf.low,ymax=conf.high,color=x),alpha=0.7)+
scale_color_manual(values=c("#30343F","red"),labels = c("no forbidden links", "with forbidden links"))+
facet_wrap(~Country)+ylab("Robustness")+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),axis.text.x=element_text(),axis.ticks.x=element_blank(),
panel.border = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_blank(),legend.position="bottom")+labs(col="",linetype="")+xlab("")+scale_x_discrete(labels=c("",""))+ggtitle("b")

pl3=ggplot(data=subset(pred2,x==1),
aes(y=odds.ratio-1,x=prop_forbidden_m,col=Country,shape=Country))+
geom_hline(yintercept=0,linetype="dashed")+
geom_pointrange(aes(ymin=odds.ratio-1-1.96*SE,ymax=odds.ratio-1+1.96*SE))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_blank(),legend.position="bottom")+
labs(col="",fill="",shape="")+ggtitle("c")+ylab("Effect of forbidden links\non robustness (odds ratio)")+xlab("Proportion of forbidden links")+
scale_color_manual(values=couleurs)+scale_y_continuous()

bas=plot_grid(pl2,pl3,align="hv",ncol=2,rel_widths=c(1,1.2))

#grid.arrange(pl1,bas,ncol=1,heights=c(1.2,1))

pdf(paste0(here("figures"),"/Figure_4_per_day.pdf"), width=7,height=8)
grid.arrange(pl1,bas,ncol=1,heights=c(1.2,1))
dev.off();


##### FIGURE 5
#extract predicts and effects for the other scenario
pred3=ggpredict(model2,c("barrierbis","sitebis","alpha"))
pred3=merge(pred3,unique(b3[,c("site","Country","prop_forbidden_m","barrierbis")]),by.x=c("group","x"),by.y=c("site","barrierbis"))
model.emm <- emmeans(model2, ~ barrierbis | sitebis + alpha,typ="response")

model.emm <- as.data.frame(contrast(model.emm, "trt.vs.ctrl", ref = "barrierbis0"))
pred3=merge(pred3,model.emm,by.x=c("group","facet"),by.y=c("sitebis","alpha"))
pred3$signifi="significant"
pred3$signifi[pred3$p.value>0.05]="not significant"
pred3$signifi=factor(pred3$signifi,levels=c("significant","not significant"))
pred3$scenario="specialists first"

predf=rbind(pred,pred3) # put all the predict for both scenario in the same table
predf=subset(predf,x==1)

pl4=ggplot(data=predf,aes(y=odds.ratio-1,x=facet))+
geom_hline(yintercept=0,linetype="dashed")+
geom_boxplot(outlier.shape = NA)+
geom_pointrange(aes(ymin=odds.ratio-1-1.96*SE,ymax=odds.ratio-1+1.96*SE,color=Country),position=position_jitter(width =0.2),alpha=0.4)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_blank(),legend.position="bottom")+
labs(col="",fill="",shape="")+ggtitle("")+ylab("Effect of forbidden links\non robustness (odds ratio)")+xlab(expression(paste("Rewiring constraint ",(alpha))))+
scale_color_manual(values=couleurs)+scale_y_continuous()+facet_wrap(~scenario,ncol=1)

pdf(paste0(here("figures"),"/Figure_5_per_day.pdf"), width=4,height=7)
pl4
dev.off();

##### FIGURE S6
bidon=unique(b3[,c("site","Country","prop_forbidden_m","barrierbis","barrier","Cperso_m","N_m","Inter_eve_m")])

dataplot=dcast(bidon,Country+site~barrier,value.var="Cperso_m")
p1 <- ggpaired(dataplot, cond1 = "no forbidden links", cond2 = "with forbidden links",
          color = "condition",palette = c("#30343F","#89023E"),
		  line.color = "gray")
# Change method
p1=p1 + stat_compare_means(paired = TRUE,method = "wilcox.test")+theme(legend.position="none",strip.background=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.x=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+ylab("Connectance")+facet_wrap(~Country)+ggtitle("a")


dataplot=dcast(bidon,Country+site~barrier,value.var="N_m")
p2 <- ggpaired(dataplot, cond1 = "no forbidden links", cond2 = "with forbidden links",
          color = "condition",palette = c("#30343F","#89023E"),
		  line.color = "gray")
# Change method
p2=p2 + stat_compare_means(paired = TRUE,method = "wilcox.test")+theme(legend.position="none",strip.background=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.x=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+ylab("Nestedness")+facet_wrap(~Country)+ggtitle("b")

plot_grid(p1,p2,align="hv",ncol=1)







