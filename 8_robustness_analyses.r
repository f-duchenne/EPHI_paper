#######################################################################################################################
###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("randomForest","data.table", "dplyr", "lubridate","lme4","R2jags","mcmcOutput","mcmcplots","MCMCvis","qgraph","igraph",
			"pastecs","ggplot2","cowplot","gridExtra","scales","reshape2","bipartite","stringr","ungeviz","lme4","mgcv","ggpubr","emmeans","piecewiseSEM","car") 


inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)

EPHI_version="2024-03-18"

couleurs=c("#679436","#0B4F6C","deeppink")


setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper/data/robustness_per_day")
extinctions=NULL
for(pays in c("Costa-Rica","Ecuador","Brazil")){
for(e in 1:100){
extinctionsc=fread(paste0("robustness_simulations_per_day_",pays,"_",e,"_.txt"))
extinctionsc$Country=pays
extinctions=rbind(extinctions,extinctionsc)
}}

extinctions=extinctions %>% group_by(r,site,Country,barrier,essai,scenario) %>% mutate(nbh=max(hummingbird_pers),nbp=max(rank))
extinctions$persistence=extinctions$hummingbird_pers/extinctions$nbh
extinctions$rank2=extinctions$rank/extinctions$nbp
extinctions=extinctions[order(extinctions$min_transect_elev),]
extinctions$site=factor(extinctions$site,levels=unique(extinctions$site))
extinctions$site2=paste0(extinctions$site," (",extinctions$min_transect_elev,"m)")
extinctions$site2=factor(extinctions$site2,levels=unique(extinctions$site2))

extinctions$Country=gsub("-"," ",extinctions$Country,fixed=T)

#### CALCULATE PERSISTENCE FOR EACH REMAINING FRACTION OF PLANT SPECIES STEP (ATAC CURVE)
b=extinctions %>% group_by(r,site,Country,site2,barrier,rank2,min_transect_elev,nbh,nbp,scenario,essai,M,prop_forbidden,C,N,Cperso,Cperso2) %>% summarise(persistence=mean(hummingbird_pers/nbh))
b$symmetrie=b$nbp/b$nbh
b2=extinctions %>% group_by(r,site,Country,site2,barrier,rank2,min_transect_elev,nbh,nbp,scenario) %>% summarise(persistence_mean=mean(hummingbird_pers/nbh),lower=quantile(hummingbird_pers/nbh,probs=0.025),
higher=quantile(hummingbird_pers/nbh,probs=0.975), prop_forbidden=-1*mean(prop_forbidden),C=mean(C),N=mean(N),M=mean(M),compl=mean(compl),Cperso=mean(Cperso),Cperso2=mean(Cperso2))


pays="Brazil"

pl1=ggplot()+
geom_line(data=subset(b,scenario=="generalists first" & r==1 & Country==pays),
aes(x=rank2,y=persistence,col=barrier,group=paste0(essai,r,site,barrier,scenario)),alpha=0.1)+
geom_line(data=subset(b2,scenario=="generalists first" & r==1 & Country==pays),aes(y=persistence_mean,x=rank2,col=barrier),size=1.5)+facet_wrap(~site2)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none")+
labs(col="")+ggtitle("a")+ylab("Persistence")+xlab("Fraction of plant species removed")+
scale_color_manual(values=c("#30343F","red"))+scale_fill_manual(values=c("#30343F","red"))+
scale_x_continuous(n.breaks=2,labels = scales::number_format(accuracy = 1))


### AREA UNDER ATAC CURVE (ROBUSTNESS)
b3=subset(b,rank2!=0 & rank2<1) %>% group_by(site,site2,Country,barrier,r,min_transect_elev,nbh,nbp,scenario,C,N,M,
prop_forbidden,symmetrie,essai) %>%
summarise(robustness=sum(persistence)/length(persistence))
b3=b3 %>% group_by(site,scenario,r,essai) %>% mutate(pers_eff=(robustness[barrier=="with forbidden links"]-robustness[barrier=="no forbidden links"])/robustness[barrier=="no forbidden links"])

b4=b3 %>% group_by(site,site2,Country,barrier,r,min_transect_elev,nbh,nbp,scenario,C,N,M,prop_forbidden,symmetrie) %>% summarise(robustness_mean=mean(robustness),sde_rob=sd(robustness)/sqrt(length(robustness)),
pers_eff_mean=mean(robustness),sde_eff=sd(pers_eff)/sqrt(length(pers_eff)))





bidon=subset(b2,scenario=="generalists first" & r==1,select=c("barrier","robustness","site","Country"))
pl2=ggpaired(bidon, x = "barrier", y = "robustness",id="site",
   color = "barrier", line.color = "gray", line.size = 0.4,
   palette = "npg")+
 stat_compare_means(paired = TRUE,size = 3,label.x.npc="left",aes(label = paste0("p = ",round(as.numeric(..p.format..),digits=3))))+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_blank(),legend.position="bottom",axis.text.x=element_text())+
labs(col="",fill="")+ggtitle("b")+ylab("Robustness")+xlab("")+scale_x_discrete(labels=c("",""))+
scale_color_manual(values=c("#30343F","red"))+scale_fill_manual(values=c("#30343F","red"))+
facet_wrap(~Country)

pl3=ggplot(data=subset(b2,scenario=="generalists first" & r==1 & barrier=="with forbidden links"),
aes(y=pers_eff,x=prop_forbidden,col=Country,shape=Country))+
geom_point(size=1.5)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_blank(),legend.position="bottom")+
labs(col="",fill="",shape="")+ggtitle("c")+ylab("Effect of forbidden links\non robustness")+xlab("Proportion of forbidden links")+
scale_color_manual(values=couleurs)+scale_y_continuous(labels = scales::percent_format())

bas=plot_grid(pl2,pl3,align="hv",ncol=2,rel_widths=c(1,1.2))

grid.arrange(pl1,bas,ncol=1,heights=c(1.2,1))

setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper")
pdf("Figure_4_per_day.pdf", width=7,height=8)
grid.arrange(pl1,bas,ncol=1,heights=c(1.2,1))
dev.off();


wilcox.test(subset(b2,scenario=="generalists first" & barrier=="no forbidden links" & r==2 & Country=="Brazil")$robustness,
subset(b2,scenario=="generalists first" & barrier=="with forbidden links" & r==2 & Country=="Brazil")$robustness,paired=T)

bidon=subset(b2,scenario=="generalists first")
figs1=ggpaired(bidon, x = "barrier", y = "robustness",id="site",
   color = "barrier", line.color = "gray", line.size = 0.4,
   palette = "npg")+
 stat_compare_means(paired = TRUE,size = 2.5,vjust=0.3)+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_blank(),legend.position="bottom",axis.text.x=element_text())+
labs(col="",fill="")+ggtitle("b")+ylab("Robustness")+xlab("")+scale_x_discrete(labels=c("",""))+
scale_color_manual(values=c("#30343F","red"))+scale_fill_manual(values=c("#30343F","red"))+
facet_grid(rows=vars(r),cols=vars(Country),scales="free",labeller = label_bquote(rows=alpha == .(r)))+
ggtitle(label="a",subtitle="generalists first")

bidon=subset(b2,scenario=="specialists first")
figs2=ggpaired(bidon, x = "barrier", y = "robustness",id="site",
   color = "barrier", line.color = "gray", line.size = 0.4,
   palette = "npg")+
 stat_compare_means(paired = TRUE,size = 2.5,vjust=0.3)+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_blank(),legend.position="bottom",axis.text.x=element_text())+
labs(col="",fill="")+ggtitle("b")+ylab("Robustness")+xlab("")+scale_x_discrete(labels=c("",""))+
scale_color_manual(values=c("#30343F","red"))+scale_fill_manual(values=c("#30343F","red"))+
facet_grid(rows=vars(r),cols=vars(Country),scales="free",labeller = label_bquote(rows=alpha == .(r)))+
ggtitle(label="b",subtitle="specialists first")

setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper")
png("Figure_S5_per_day.png", width=1200,height=1900,res=180)
grid.arrange(figs1,figs2)
dev.off();


######################################## DETERMINISM OF Robustness

b3=subset(b,rank2!=0 & rank2<1) %>% group_by(site,site2,Country,barrier,r,min_transect_elev,nbh,nbp,scenario,C,N,M,
prop_forbidden,compl,Cperso2) %>%
summarise(robustness=sum(persistence)/length(persistence))
#b3=merge(b3,forb,by=c("site","min_transect_elev","Country","nbp","nbh"))
deter=as.data.frame(subset(b3,scenario=="generalists first" & r==1))
deter$symmetrie=deter$nbp/deter$nbh
deter$dive=deter$nbh+deter$nbp
deter$Clog=deter$Cperso2
deter$Mlogit=logit(deter$M,adjust=0.01)
deter$divelog=log(deter$dive)
deter$Country2=as.numeric(as.factor(deter$Country))

deter$roblogit=logit(deter$robustness,adjust=0.01)

dataplot=dcast(deter,Country+site~barrier,value.var="robustness")

p <- ggpaired(dataplot, cond1 = "no forbidden links", cond2 = "with forbidden links",
          color = "condition",palette = c("#30343F","#89023E"),
		  line.color = "gray")
# Change method
p + stat_compare_means(paired = TRUE,method = "wilcox.test")+theme(legend.position="none")+ylab("Robustness")+facet_wrap(~Country)

cor(deter[,c("dive","prop_forbidden","compl","nbh","nbp","M","N")])

model_C=lm(Clog~dive+prop_forbidden,data=deter)
model_N=lm(N~dive+prop_forbidden+Clog,data=deter)
model_M=lm(Mlogit~dive+prop_forbidden+Clog,data=deter)

model_rob=lm(roblogit~dive+prop_forbidden+Clog+N,data=deter)
print(vif(model_rob))
obj <- psem(model_C,model_N,model_M,model_rob,data=as.data.frame(deter), dive %~~% prop_forbidden)

objb=summary(obj)
left=0
right=30
center=(right-left)/2
haut=20
bas=0
cex=1.5
cex.arrows=0.7
echelle_fleche=0.5
l=objb$coefficients[,c("Predictor","Response","Std.Estimate","P.Value")]
l=l[grep("~~",l$Predictor,fixed=T,invert=T),]
l$colo=ifelse(l$Std.Estimate<0,"firebrick4","dodgerblue4")
l$lty=1
l$lty=ifelse(l$P.Value>0.05,2,l$lty)
l$curv=NA
l$curv[l$Predictor=="dive" & l$Response=="roblogit"]=-4.3
l$curv[l$Predictor=="prop_forbidden" & l$Response=="roblogit"]=4.3
g <- graph.data.frame(l, directed=T)
g= g %>% set_edge_attr("color", value =l$colo)
g= g %>% set_vertex_attr("name", value =
c("Species Richness","Proportion of\nforbidden links",
paste0("Connectance\n(R2 = ",objb$R2$R.squared[objb$R2$Response=="Clog"],")"),
paste0("Nestedness\n(R2 = ",objb$R2$R.squared[objb$R2$Response=="N"],")"),
paste0("Modularity\n(R2 = ",objb$R2$R.squared[objb$R2$Response=="Mlogit"],")"),
paste0("Robustness\n(R2 = ",objb$R2$R.squared[objb$R2$Response=="roblogit"],")")))
coord=data.frame(label=vertex_attr(g, "name"),
x=c(left,right,center,left,right,center),
y=c(haut,haut,bas+11,bas+11,
bas+11,bas),vsize=20)

EL=as_edgelist(g)
EL=cbind(EL,l[,3])
asi=abs(l[,3])*echelle_fleche
asi[asi<5]=5
asi[asi>=15]=15


setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper")
pdf("Figure5_per_day.pdf", width=8,height=6)
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
b3=subset(b,rank2!=0 & rank2<1) %>% group_by(site,site2,Country,barrier,r,min_transect_elev,nbh,nbp,scenario,C,N,M,
prop_forbidden,compl,Cperso2) %>%
summarise(robustness=sum(persistence)/length(persistence))
#b3=merge(b3,forb,by=c("site","min_transect_elev","Country","nbp","nbh"))
deter=as.data.frame(subset(b3,scenario==sc & r==alpha))
deter$symmetrie=deter$nbp/deter$nbh
deter$dive=deter$nbh+deter$nbp
deter$Clog=deter$Cperso2
deter$Mlogit=logit(deter$M,adjust=0.01)
deter$divelog=log(deter$dive)
deter$Country2=as.numeric(as.factor(deter$Country))

deter$roblogit=logit(deter$robustness,adjust=0.01)

dataplot=dcast(deter,Country+site~barrier,value.var="robustness")

p <- ggpaired(dataplot, cond1 = "no forbidden links", cond2 = "with forbidden links",
          color = "condition",palette = c("#30343F","#89023E"),
		  line.color = "gray")
# Change method
p + stat_compare_means(paired = TRUE,method = "wilcox.test")+theme(legend.position="none")+ylab("Robustness")+facet_wrap(~Country)

cor(deter[,c("dive","prop_forbidden","compl","nbh","nbp")])

model_C=lm(Clog~dive+prop_forbidden,data=deter)
model_N=lm(N~dive+prop_forbidden+Clog,data=deter)
model_M=lm(Mlogit~dive+prop_forbidden+Clog,data=deter)

model_rob=lm(roblogit~dive+prop_forbidden+Clog+N+Mlogit,data=deter)
print(vif(model_rob))
obj <- psem(model_C,model_N,model_M,model_rob,data=as.data.frame(deter), dive %~~% prop_forbidden)

objb=summary(obj)
left=0
right=30
center=(right-left)/2
haut=20
bas=0
cex=1.5
cex.arrows=0.7
echelle_fleche=0.5
l=objb$coefficients[,c("Predictor","Response","Std.Estimate","P.Value")]
l=l[grep("~~",l$Predictor,fixed=T,invert=T),]
l$colo=ifelse(l$Std.Estimate<0,"firebrick4","dodgerblue4")
l$lty=1
l$lty=ifelse(l$P.Value>0.05,2,l$lty)
l$curv=NA
l$curv[l$Predictor=="dive" & l$Response=="roblogit"]=-4.3
l$curv[l$Predictor=="prop_forbidden" & l$Response=="roblogit"]=4.3
g <- graph.data.frame(l, directed=T)
g= g %>% set_edge_attr("color", value =l$colo)
g= g %>% set_vertex_attr("name", value =
c("Species Richness","Proportion of\nforbidden links",
paste0("Connectance\n(R2 = ",objb$R2$R.squared[objb$R2$Response=="Clog"],")"),
paste0("Nestedness\n(R2 = ",objb$R2$R.squared[objb$R2$Response=="N"],")"),
paste0("Modularity\n(R2 = ",objb$R2$R.squared[objb$R2$Response=="Mlogit"],")"),
paste0("Robustness\n(R2 = ",objb$R2$R.squared[objb$R2$Response=="roblogit"],")")))
coord=data.frame(label=vertex_attr(g, "name"),
x=c(left,right,center,left,right,center),
y=c(haut,haut,bas+11,bas+11,
bas+11,bas),vsize=20)

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

