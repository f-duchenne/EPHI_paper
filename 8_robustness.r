library(dplyr)
library(reshape2)
library(bipartite)
library(data.table)

setwd(dir="/home/duchenne/EPHI_paper/robust_data_res/")

# Collect command arguments
args <- commandArgs(trailingOnly = TRUE)
args_contents <- strsplit(args, ' ')
# Get first argument
jj <- as.numeric(args_contents[[1]])
print(jj)

tab=expand.grid(c("Brazil","Costa-Rica","Ecuador"),1:100)
pays=tab$Var1[jj]
e=tab$Var2[jj]

pre1=fread(paste0("initial_network_",pays,".txt"))
sites=fread(paste0("Site_metadata_",pays,".txt"),na.strings = c("",NA))
sites=subset(sites,habitat!="deforested")

r_vec=0:2
scenarios=c("generalists first","specialists first")
rules=c("with forbidden links","no forbidden links")
sites_vec=unique(pre1$site)

extinctions=NULL
pre1b=pre1
for(r in r_vec){
for(scenario in scenarios){
for(bar in rules){

for(z in sites_vec){
bidon=subset(pre1b,site==z)
if(bar=="no forbidden links"){
bidon$binary=1
bidon$average_proba=1
}
bidon$inter=bidon$freq*bidon$average_proba
bidon=bidon %>% dplyr::group_by(plant_species) %>% mutate(degree=length(unique(hummingbird_species[inter>0])))
if(bar=="generalists first"){bidon=bidon[order(bidon$degree,bidon$plant_species,decreasing=TRUE),]}else{
bidon=bidon[order(bidon$degree,bidon$plant_species,decreasing=FALSE),]}

prop_forbidden=mean(1-bidon$average_proba)
average_freq=mean(bidon$freq)

mat=dcast(bidon,plant_species~hummingbird_species,value.var="inter")
mat=as.matrix(mat[,-1])
C=networklevel(mat,index="weighted connectance",weighted=T)
N=wine(mat, nreps = 1000)$wine
M=computeModules(mat, method="Beckett")@likelihood[1]

plant_list=unique(bidon$plant_species)
extinctions=rbind(extinctions,data.frame(plant_ext=NA,hummingbird_pers=length(unique(bidon$hummingbird_species)),site=z,essai=e,r=r,
rank=0,barrier=bar,prop_forbidden=prop_forbidden,average_freq=average_freq,scenario=scenario,C=C,N=N,M=M))
for(j in 1:length(plant_list)){
ext=plant_list[j]
if(length(unique(bidon$plant_species))>1){
for(i in unique(bidon$hummingbird_species)){
bidon2=subset(bidon,hummingbird_species==i)
distances=sqrt((bidon2$compl[bidon2$plant_species==ext]-bidon2$compl[bidon2$plant_species!=ext])^2)
R=mean(bidon2$average_proba[bidon2$plant_species!=ext]*exp(-r*distances/max(distances)))

Pext=(1-R)*bidon2$inter[bidon2$plant_species==ext]/sum(bidon2$inter)

EXT=rbinom(1,1,Pext)

if(EXT==1){bidon=subset(bidon,hummingbird_species!=i)}

}
}
bidon=subset(bidon,plant_species!=ext)
extinctions=rbind(extinctions,data.frame(plant_ext=ext,hummingbird_pers=length(unique(bidon$hummingbird_species)),site=z,essai=e,r=r,
rank=j,barrier=bar,prop_forbidden=prop_forbidden,average_freq=average_freq,scenario=scenario,C=C,N=N,M=M))


}
}
}
}
}


#LOAD SITE METADATA:
extinctionf=merge(extinctionf,sites,by="site")

setwd(dir="/home/duchenne/EPHI_paper/")
fwrite(extinctionf,paste0("robustness_simulations_",pays,"_",e,"_.txt"))



#######################################################################################################################
source("C:/Users/Duchenne/Documents/EPHI_paper/scripts/function_to_predict_from_bayesian.r")
library(bipartite)
library(Rmpfr)
library(circular)
library(CircStats)
library(plot3D)
library(ggradar)
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
library(DHARMa)
library(glmmTMB)
library(qgraph)
library(igraph)
library(piecewiseSEM)
library(randomForest)
library(FactoMineR)
require(factoextra)
library("ggdendro")
library(dendextend)
library(ape)
library(ggtree)
library(ggnewscale)
library(geiger)
library(diversityForest)
library(caper)
library(phytools)
library(ggthemes)
library(ggplotify)
library(ggtern)
library(ks)
library(sp)
library(mgcv)
library(spatialEco)
library(ggpubr)
library(emmeans)

setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper/data/robustness")
extinctions=NULL
for(pays in c("Costa-Rica","Ecuador","Brazil")){
for(e in 1:100){
extinctionsc=fread(paste0("robustness_simulations_",pays,"_",e,"_.txt"))
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

b=extinctions %>% group_by(r,site,Country,site2,barrier,rank2,min_transect_elev,nbh,nbp,scenario) %>% summarise(persistence=mean(hummingbird_pers/nbh),
pers_sde=sd(hummingbird_pers)/sqrt(length(hummingbird_pers)),prop_forbidden=-1*mean(prop_forbidden),
C=mean(C),N=mean(N),M=mean(M),compl=mean(compl),Cperso=mean(Cperso),Cperso2=mean(Cperso2))
b$symmetrie=b$nbp/b$nbh

ggplot()+
geom_line(data=b,aes(y=persistence,x=rank2,col=barrier,linetype=scenario),size=1)+facet_grid(rows=vars(site),cols=vars(r))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA))+
labs(col="")+ggtitle("a")+ylab("Persistence")+xlab("Fraction of plant species removed")+
scale_color_manual(values=c("#776274","#C8CC92"))+scale_fill_manual(values=c("#776274","#C8CC92"))+
scale_x_continuous(n.breaks=3,labels = scales::number_format(accuracy = 0.1))

pays="Brazil"

pl1=ggplot()+
geom_line(data=subset(extinctions,scenario=="generalists first" & r==1 & Country==pays),
aes(x=rank2,y=persistence,col=barrier,group=paste0(essai,r,site,barrier,scenario)),alpha=0.1)+
geom_line(data=subset(b,scenario=="generalists first" & r==1 & Country==pays),aes(y=persistence,x=rank2,col=barrier),size=1.5)+facet_wrap(~site2)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none")+
labs(col="")+ggtitle("a")+ylab("Persistence")+xlab("Fraction of plant species removed")+
scale_color_manual(values=c("#30343F","red"))+scale_fill_manual(values=c("#30343F","red"))+
scale_x_continuous(n.breaks=2,labels = scales::number_format(accuracy = 1))


###AREA
b2=subset(b,rank2!=0 & rank2<1) %>% group_by(site,site2,Country,barrier,r,min_transect_elev,nbh,nbp,scenario,C,N,M,
prop_forbidden,symmetrie) %>%
summarise(robustness=sum(persistence)/length(persistence))
b2=b2 %>% group_by(site,scenario,r) %>% mutate(pers_eff=(robustness[barrier=="with forbidden links"]-robustness[barrier=="no forbidden links"])/robustness[barrier=="no forbidden links"])
b2=b2 %>% group_by(Country) %>%
mutate(divemax=max(nbp+nbh),low_forbidden=min(prop_forbidden[barrier=="with forbidden links"]),symmetrie_max=max(symmetrie))

# model=glm(robustness~poly(min_transect_elev,2)*barrier*as.factor(Country),data=subset(b2,scenario=="generalists first" & r==1),family=quasibinomial)
# AIC(model)
# pre=ggpredict(model,c("min_transect_elev[10:3410]","barrier","Country"))
# names(pre)[names(pre) %in% c("group","facet")]=c("barrier","Country")
# lims=b2 %>% group_by(Country) %>% summarise(mini=min(min_transect_elev),maxi=max(min_transect_elev))
# pre=merge(pre,lims,by="Country")
# pre[pre$x<pre$mini | pre$x>pre$maxi,c("conf.high","conf.low","predicted")]=NA
# b2$dive=(b2$nbh+b2$nbp)/b2$divemax

# pl2=ggplot()+
# geom_ribbon(data=pre,aes(x=x,y=predicted,ymin=conf.low,ymax=conf.high,fill=barrier),alpha=0.2)+
# geom_line(data=pre,aes(x=x,y=predicted,col=barrier),size=1.3)+
# geom_point(data=subset(b2,scenario=="generalists first" & r==1),
# aes(y=robustness,x=min_transect_elev,col=barrier,fill=barrier),size=1.5)+
# geom_vline(data=subset(b2,scenario=="generalists first" & r==1 & symmetrie==symmetrie_max),
# aes(xintercept=min_transect_elev),col="darkgrey",linetype="dotdash")+
# geom_vline(data=subset(b2,scenario=="generalists first" & r==1 & prop_forbidden==low_forbidden),
# aes(xintercept=min_transect_elev),col="lightpink",linetype="dashed")+
# theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
# panel.grid.minor = element_blank(),
# panel.border = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
# strip.background=element_blank(),legend.position="bottom")+
# labs(col="",fill="")+ggtitle("b")+ylab("Robustness")+xlab("Elevation (meters above sea level)")+
# scale_color_manual(values=c("#30343F","red"))+scale_fill_manual(values=c("#30343F","red"))+
# facet_wrap(~Country)

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


couleurs=c("#679436","#0B4F6C","deeppink")

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
pdf("Figure_4.pdf", width=7,height=8)
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
png("Figure_S5.png", width=1200,height=1900,res=180)
grid.arrange(figs1,figs2)
dev.off();


######################################## DETERMINISM OF Robustness
b3=subset(b,rank2!=0 & rank2<1) %>% group_by(site,site2,Country,barrier,r,min_transect_elev,nbh,nbp,scenario,C,N,M,
prop_forbidden,compl,Cperso2) %>%
summarise(robustness=sum(persistence)/length(persistence))
#b3=merge(b3,forb,by=c("site","min_transect_elev","Country","nbp","nbh"))
deter=as.data.frame(subset(b3,scenario=="specialists first" & r==1))
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

cor(deter[,c("dive","prop_forbidden","compl","nbh","nbp","symmetrie")])

model_C=lm(Clog~dive+prop_forbidden+symmetrie,data=deter)
model_N=lm(N~dive+prop_forbidden+Clog+symmetrie,data=deter)
model_M=lm(Mlogit~dive+prop_forbidden+Clog+symmetrie,data=deter)

model_rob=lm(roblogit~dive+prop_forbidden+Clog+N+Mlogit+symmetrie,data=deter)

obj <- psem(model_C,model_N,model_M,model_rob,data=as.data.frame(deter), dive %~~% symmetrie, dive %~~% prop_forbidden)

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
l$curv[l$Predictor=="dive" & l$Response=="roblogit"]=1.7
l$curv[l$Predictor=="symmetrie" & l$Response=="roblogit"]=-4.3
l$curv[l$Predictor=="prop_forbidden" & l$Response=="roblogit"]=4.3
g <- graph.data.frame(l, directed=T)
g= g %>% set_edge_attr("color", value =l$colo)
g= g %>% set_vertex_attr("name", value =
c("Species Richness","Proportion of\nforbidden links","Np/Nh",
paste0("Connectance\n(R2 = ",objb$R2$R.squared[objb$R2$Response=="Clog"],")"),
paste0("Nestedness\n(R2 = ",objb$R2$R.squared[objb$R2$Response=="N"],")"),
paste0("Modularity\n(R2 = ",objb$R2$R.squared[objb$R2$Response=="Mlogit"],")"),
paste0("Robustness\n(R2 = ",objb$R2$R.squared[objb$R2$Response=="roblogit"],")")))
coord=data.frame(label=vertex_attr(g, "name"),
x=c(center,right,left,center,left,right,center),
y=c(haut,haut,haut,bas+11,bas+11,
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
mar=c(5,5,5,5),knot.border.color="white",curveShape=-1)



setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper")
pdf("Figure5.pdf", width=8,height=6)
qgraph(EL,layout=as.matrix(coord[,c("x","y")]),edge.color=l$colo,
border.color="white",label.cex=cex,label.scale=F,lty=l$lty,
edge.label.cex = cex.arrows,edge.label.position=0.65,vsize2=9,vsize=coord$vsize,curve=l$curv,
shape="rectangle",edge.labels=T,fade=F,asize=asi,
mar=c(5,5,5,5),knot.border.color="white",curveShape=-1.2)
dev.off();








mult=multigroup(obj,"Country2")



#######################
b2b=dcast(b2,site+site2+Country+r+min_transect_elev+nbh+nbp+scenario~barrier,value.var="robustness")

b2b$effet=b2b[,10]-b2b[,9]

model=lm(effet~min_transect_elev+Country,data=subset(b2b,scenario=="generalists first" & r==1))
AIC(model)
pre=ggpredict(model,c("min_transect_elev[10:3410]","Country"))
names(pre)[names(pre) %in% c("group","facet")]=c("Country")
lims=b2 %>% group_by(Country) %>% summarise(mini=min(min_transect_elev),maxi=max(min_transect_elev))
pre=merge(pre,lims,by="Country")
pre[pre$x<pre$mini | pre$x>pre$maxi,c("conf.high","conf.low","predicted")]=NA

ggplot()+
geom_ribbon(data=pre,aes(x=x,y=predicted,ymin=conf.low,ymax=conf.high,fill=Country),alpha=0.2)+
geom_line(data=pre,aes(x=x,y=predicted,col=Country),size=1.3)+
geom_point(data=subset(b2b,scenario=="generalists first" & r==1) ,aes(y=effet,x=min_transect_elev,col=Country,fill=Country),size=1.5)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+
labs(col="",fill="")+ggtitle("b")+ylab("Effect of forbidden links on robustness")+xlab("Elevation (meters above sea level)")

scale_color_manual(values=c("#30343F","#89023E"))+scale_fill_manual(values=c("#30343F","#89023E"))+
facet_wrap(~Country)

pairs(deter[,c("min_transect_elev","prop_forbidden","compl","nbh","nbp","symmetrie")])
cor(deter[,c("min_transect_elev","prop_forbidden","compl","nbh","nbp","symmetrie")])


model=randomForest(robustness~prop_forbidden+compl+nbp+nbh+symmetrie,data=deter,importance=T,ntree=4000)
varImpPlot(model)
partialPlot(model, pred.data=deter,x.var="prop_forbidden")
partialPlot(model, pred.data=deter,x.var="compl")

