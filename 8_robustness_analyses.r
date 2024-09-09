#######################################################################################################################
###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("randomForest","data.table", "dplyr", "lubridate","lme4","qgraph","igraph","ggeffects","DHARMa",
			"pastecs","ggplot2","cowplot","gridExtra","scales","reshape2","bipartite","stringr","mgcv","ggpubr","emmeans","piecewiseSEM","car","glmmTMB") 


inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)

EPHI_version="2024-03-18"

couleurs=c("#679436","#0B4F6C","deeppink")

### PUT all files together
setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper/data")

extinctions=fread("robustness_simulations_all.csv")

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


### FIT STATISTICAL MODEL
b3$barrierbis=0
b3$barrierbis[b3$barrier=="with forbidden links"]=1
b3$barrierbis=as.factor(b3$barrierbis)
b3$sitebis=as.factor(b3$site)
b3$alpha=as.factor(b3$r)

#when generalists get extinct first
model=glmmTMB(robustness~barrierbis*sitebis*alpha,data=subset(b3,scenario=="generalists first"),family=beta_family(),control=glmmTMBControl(optCtrl=list(iter.max=1e5,eval.max=1e5)))

#when specialists get extinct first
model2=glmmTMB(robustness~barrierbis*sitebis*alpha,data=subset(b3,scenario=="specialists first"),family=beta_family(),control=glmmTMBControl(optCtrl=list(iter.max=1e5,eval.max=1e5)))


##### FIGURE 4, let's consider only the first scenario here
pred=ggpredict(model,c("barrierbis","sitebis","alpha"))
pred=merge(pred,unique(b3[,c("site","Country","prop_forbidden_m","barrierbis")]),by.x=c("group","x"),by.y=c("site","barrierbis"))

model.emm <- emmeans(model, ~ barrierbis | sitebis + alpha,typ="response")

model.emm <- as.data.frame(contrast(model.emm, "trt.vs.ctrl", ref = "barrierbis0"))
pred=merge(pred,model.emm,by.x=c("group","facet"),by.y=c("sitebis","alpha"))
pred$signifi="significant"
pred$signifi[pred$p.value>0.05]="not significant"
pred$signifi=factor(pred$signifi,levels=c("significant","not significant"))
pred$scenario="generalists first"

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

setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper")
pdf("Figure_4_per_day.pdf", width=7,height=8)
grid.arrange(pl1,bas,ncol=1,heights=c(1.2,1))
dev.off();


##### FIGURE 5
#extract predicts and effects for the other scenario
pred3=ggpredict(model2,c("barrierbis","sitebis","alpha"))
pred3=merge(pred3,unique(bidon2[,c("site","Country","prop_forbidden_m","barrierbis")]),by.x=c("group","x"),by.y=c("site","barrierbis"))
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


setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper")
pdf("Figure_5_per_day.pdf", width=4,height=7)
pl4
dev.off();





##### FIGURE S6
bidon=unique(b3[,c("site","Country","prop_forbidden_m","barrierbis","barrier","Cperso_m","N_m","Inter_eve_m")])

dataplot=dcast(deter,Country+site~barrier,value.var="Cperso_m")
p1 <- ggpaired(dataplot, cond1 = "no forbidden links", cond2 = "with forbidden links",
          color = "condition",palette = c("#30343F","#89023E"),
		  line.color = "gray")
# Change method
p1=p1 + stat_compare_means(paired = TRUE,method = "wilcox.test")+theme(legend.position="none",strip.background=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.x=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+ylab("Connectance")+facet_wrap(~Country)+ggtitle("a")


dataplot=dcast(deter,Country+site~barrier,value.var="N_m")
p2 <- ggpaired(dataplot, cond1 = "no forbidden links", cond2 = "with forbidden links",
          color = "condition",palette = c("#30343F","#89023E"),
		  line.color = "gray")
# Change method
p2=p2 + stat_compare_means(paired = TRUE,method = "wilcox.test")+theme(legend.position="none",strip.background=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.x=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+ylab("Nestedness")+facet_wrap(~Country)+ggtitle("b")

plot_grid(p1,p2,align="hv",ncol=1)













tab=expand.grid(c("generalists first","specialists first"),c(0,1,2))
predf=NULL
for(i in 1:nrow(tab)){
bidon=subset(b3,scenario==tab$Var1[i] & r==tab$Var2[i])
bidon$barrierbis=0
bidon$barrierbis[bidon$barrier=="with forbidden links"]=1
bidon$barrierbis=as.factor(bidon$barrierbis)
bidon$sitebis=as.factor(bidon$site)
model=glmmTMB(robustness~barrierbis*sitebis,data=bidon,family=beta_family(),control=glmmTMBControl(optCtrl=list(iter.max=1e5,eval.max=1e5)))
# resi= simulateResiduals(model)
# plot(resi)


pred=ggpredict(model,c("barrierbis","sitebis"))
pred=merge(pred,unique(bidon[,c("site","Country","prop_forbidden_m","barrierbis")]),by.x=c("group","x"),by.y=c("site","barrierbis"))
pred$scenario=tab$Var1[i]
pred$r=tab$Var2[i]

model.emm <- emmeans(model, ~ barrierbis | sitebis,typ="response")
model.emm <- as.data.frame(contrast(model.emm, "trt.vs.ctrl", ref = "barrierbis0"))
pred=merge(pred,model.emm,by.x="group",by.y="sitebis")
pred$signifi="significant"
pred$signifi[pred$p.value>0.05]="not significant"
pred$signifi=factor(pred$signifi,levels=c("significant","not significant"))

predf=rbind(pred,predf)

}

figs1=ggplot(data=subset(predf,scenario=="generalists first"),aes(x=x,y=predicted))+
geom_boxplot(aes(color=x))+
geom_line(aes(linetype=signifi,group=group),show.legend = F)+
geom_pointrange(aes(ymin=conf.low,ymax=conf.high,color=x),alpha=0.7)+
scale_color_manual(values=c("#30343F","red"),labels = c("no forbidden links", "with forbidden links"))+
facet_wrap(~Country)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),axis.text.x=element_text(),axis.ticks.x=element_blank(),
panel.border = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_blank(),legend.position="bottom")+labs(col="",linetype="")+xlab("")+scale_x_discrete(labels=c("",""))+
facet_grid(rows=vars(r),cols=vars(Country),scales="free",labeller = label_bquote(rows=alpha == .(r)))+
ggtitle(label="a",subtitle="generalists first")

figs2=ggplot(data=subset(predf,scenario=="specialists first"),aes(x=x,y=predicted))+
geom_boxplot(aes(color=x))+
geom_line(aes(linetype=signifi,group=group),show.legend = F)+
geom_pointrange(aes(ymin=conf.low,ymax=conf.high,color=x),alpha=0.7)+
scale_color_manual(values=c("#30343F","red"),labels = c("no forbidden links", "with forbidden links"))+
facet_wrap(~Country)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),axis.text.x=element_text(),axis.ticks.x=element_blank(),
panel.border = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_blank(),legend.position="bottom")+labs(col="",linetype="")+xlab("")+scale_x_discrete(labels=c("",""))+
facet_grid(rows=vars(r),cols=vars(Country),scales="free",labeller = label_bquote(rows=alpha == .(r)))+
ggtitle(label="b",subtitle="specialists first")

setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper")
png("Figure_S5_per_day.png", width=1200,height=1900,res=180)
grid.arrange(figs1,figs2)
dev.off();



######################################## DETERMINISM OF Robustness
bidon=subset(b3,scenario=="generalists first")
bidon$barrierbis=0
bidon$barrierbis[bidon$barrier=="with forbidden links"]=1
bidon$barrierbis=as.factor(bidon$barrierbis)
bidon$sitebis=as.factor(bidon$site)
model=glmmTMB(robustness~barrierbis*sitebis*as.factor(r),data=bidon,family=beta_family(),control=glmmTMBControl(optCtrl=list(iter.max=1e5,eval.max=1e5)))
# resi= simulateResiduals(model)
# plot(resi)


pred=ggpredict(model,c("barrierbis","sitebis"))
pred=merge(pred,unique(bidon[,c("site","Country","prop_forbidden_m","barrierbis")]),by.x=c("group","x"),by.y=c("site","barrierbis"))

model.emm <- emmeans(model, ~ barrierbis | sitebis+r,typ="response")
model.emm <- as.data.frame(contrast(model.emm, "trt.vs.ctrl", ref = "barrierbis0"))

ggplot(data=model.emm,aes(y=odds.ratio-1,x=r,color=Country))+
geom_hline(yintercept=0,linetype="dashed")+
geom_pointrange(aes(ymin=odds.ratio-1-1.96*SE,ymax=odds.ratio-1+1.96*SE),position=position_jitter(width =0.1))+
geom_boxplot()

b4=b3 %>% group_by(site,site2,Country,barrier,r,min_transect_elev,nbh,nbp,scenario,C_m,N_m,Cperso_m,Inter_eve_m,
 prop_forbidden_m) %>% summarise(robustness_m=mean(robustness))
#b3=merge(b3,forb,by=c("site","min_transect_elev","Country","nbp","nbh"))
setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper")

pdf("Figure5_per_day.pdf", width=7,height=10)
par(mfrow=c(2,1))
for( re in 0:1){
deter=as.data.frame(subset(b4,scenario=="specialists first" & r==re))
deter$symmetrie=deter$nbp/deter$nbh
deter$dive=deter$nbh+deter$nbp
deter$divelog=log(deter$dive)
deter$Country2=as.numeric(as.factor(deter$Country))

deter$roblogit=logit(deter$robustness_m,adjust=0.01)

cor(deter[,c("dive","prop_forbidden_m","nbh","nbp","Inter_eve_m","N_m")])

model_C=lm(C_m~dive+prop_forbidden_m,data=deter)
model_N=lm(N_m~dive+prop_forbidden_m+C_m,data=deter)

model_rob=lm(roblogit~dive+prop_forbidden_m+N_m+C_m,data=deter)
# resi = simulateResiduals(model_rob)
# plot(resi)
print(vif(model_rob))
obj <- psem(model_C,model_N,model_rob,data=as.data.frame(deter), dive %~~% prop_forbidden_m)

objb=summary(obj)
left=0
right=30
center=(right-left)/2
haut=20
bas=0
cex=1.5
cex.arrows=1.5
echelle_fleche=0.5
l=objb$coefficients[,c("Predictor","Response","Std.Estimate","P.Value")]
l=l[grep("~~",l$Predictor,fixed=T,invert=T),]
l$colo=ifelse(l$Std.Estimate<0,"firebrick4","dodgerblue4")
l$lty=1
l$lty=ifelse(l$P.Value>0.05,2,l$lty)
l$curv=NA
l$curv[l$Predictor=="dive" & l$Response=="roblogit"]=-4.3
l$curv[l$Predictor=="prop_forbidden_m" & l$Response=="roblogit"]=4.3
g <- graph.data.frame(l, directed=T)
g= g %>% set_edge_attr("color", value =l$colo)
g= g %>% set_vertex_attr("name", value =
c("Species Richness","Proportion of\nforbidden links",
paste0("Interaction eveness\n(R2 = ",objb$R2$R.squared[objb$R2$Response=="Inter_eve_m"],")"),
paste0("Nestedness\n(R2 = ",objb$R2$R.squared[objb$R2$Response=="N_m"],")"),
paste0("Robustness\n(R2 = ",objb$R2$R.squared[objb$R2$Response=="roblogit"],")")))
coord=data.frame(label=vertex_attr(g, "name"),
x=c(left,right,right-2,left+2,center),
y=c(haut,haut,bas+11,bas+11,bas),vsize=20)

EL=as_edgelist(g)
EL=cbind(EL,l[,3])
asi=abs(l[,3])*echelle_fleche
asi[asi<5]=5
asi[asi>=15]=15


qgraph(EL,layout=as.matrix(coord[,c("x","y")]),edge.color=l$colo,
border.color="white",label.cex=cex,label.scale=F,lty=l$lty,
edge.label.cex = cex.arrows,edge.label.position=0.65,vsize2=9,vsize=coord$vsize,curve=l$curv,
shape="rectangle",edge.labels=T,fade=F,asize=asi,
mar=c(5,5,5,5),knot.border.color="white",curveShape=-1.2,title=paste0(letters[re+1]," - ",ifelse(re==0,paste0("high rewiring (alpha = ",re," )"),paste0("moderate rewiring (alpha = ",re," )"))))

}
dev.off();


qgraph(EL,layout=as.matrix(coord[,c("x","y")]),edge.color=l$colo,
border.color="white",label.cex=cex,label.scale=F,lty=l$lty,
edge.label.cex = cex.arrows,edge.label.position=0.65,vsize2=9,vsize=coord$vsize,curve=l$curv,
shape="rectangle",edge.labels=T,fade=F,asize=asi,
mar=c(5,5,5,5),knot.border.color="white",curveShape=-1.2)
dev.off();




######################### SUPP
pdf("Figure_S6.pdf", width=14,height=9)
par(mfrow=c(2,3))
for(sc in c("generalists first","specialists first")){
for(alpha in c(0,1,2)){
deter=as.data.frame(subset(b4,scenario==sc & r==alpha))
deter$symmetrie=deter$nbp/deter$nbh
deter$dive=deter$nbh+deter$nbp
deter$divelog=log(deter$dive)
deter$Country2=as.numeric(as.factor(deter$Country))

deter$roblogit=logit(deter$robustness_m,adjust=0.01)

model_C=lm(Inter_eve_m~dive+prop_forbidden_m,data=deter)
model_N=lm(N_m~dive+prop_forbidden_m+Inter_eve_m,data=deter)

model_rob=lm(roblogit~dive+prop_forbidden_m+N_m+Inter_eve_m,data=deter)
# resi = simulateResiduals(model_rob)
# plot(resi)
print(vif(model_rob))
obj <- psem(model_C,model_N,model_rob,data=as.data.frame(deter), dive %~~% prop_forbidden_m)

objb=summary(obj)
left=0
right=30
center=(right-left)/2
haut=20
bas=0
cex=1.5
cex.arrows=1.5
echelle_fleche=0.5
l=objb$coefficients[,c("Predictor","Response","Std.Estimate","P.Value")]
l=l[grep("~~",l$Predictor,fixed=T,invert=T),]
l$colo=ifelse(l$Std.Estimate<0,"firebrick4","dodgerblue4")
l$lty=1
l$lty=ifelse(l$P.Value>0.05,2,l$lty)
l$curv=NA
l$curv[l$Predictor=="dive" & l$Response=="roblogit"]=-4.3
l$curv[l$Predictor=="prop_forbidden_m" & l$Response=="roblogit"]=4.3
g <- graph.data.frame(l, directed=T)
g= g %>% set_edge_attr("color", value =l$colo)
g= g %>% set_vertex_attr("name", value =
c("Species Richness","Proportion of\nforbidden links",
paste0("Interaction eveness\n(R2 = ",objb$R2$R.squared[objb$R2$Response=="Inter_eve_m"],")"),
paste0("Nestedness\n(R2 = ",objb$R2$R.squared[objb$R2$Response=="N_m"],")"),
paste0("Robustness\n(R2 = ",objb$R2$R.squared[objb$R2$Response=="roblogit"],")")))
coord=data.frame(label=vertex_attr(g, "name"),
x=c(left,right,right-2,left+2,center),
y=c(haut,haut,bas+11,bas+11,bas),vsize=20)

EL=as_edgelist(g)
EL=cbind(EL,l[,3])
asi=abs(l[,3])*echelle_fleche
asi[asi<5]=5
asi[asi>=15]=15

qgraph(EL,layout=as.matrix(coord[,c("x","y")]),edge.color=l$colo,
border.color="white",label.cex=cex,label.scale=F,lty=l$lty,
edge.label.cex = cex.arrows,edge.label.position=0.65,vsize2=9,vsize=coord$vsize,curve=l$curv,
shape="rectangle",edge.labels=T,fade=F,asize=asi,
mar=c(5,5,5,5),knot.border.color="white",curveShape=-1,title=paste(sc,alpha,sep=" - alpha = "),title.cex=2)
}
}
dev.off();

