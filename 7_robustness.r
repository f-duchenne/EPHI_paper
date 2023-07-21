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
library(dplyr)
library(reshape2)
library(bipartite)
library(stringr)
library(ungeviz)
library(data.table)
library(doParallel)
library(foreach)
library(ggpubr)
library(parallel)
ncpus <- makeCluster(5)
registerDoParallel(ncpus)

for(pays in c("Costa-Rica","Ecuador","Brazil")){
pre1b=fread(paste0("initial_network_",pays,".txt"))
sites=fread(paste0("C:/Users/Duchenne/Documents/EPHI_data_clean/",pays,"_",EPHI_version,"/Site_metadata_",pays,".txt"),na.strings = c("",NA))
sites=subset(sites,habitat!="deforested")

r_vec=0:2
scenarios=c("generalists first","specialists first")
rules=c("with forbidden links","no forbidden links")
sites_vec=unique(pre1b$site)


extinctionf=foreach(e=1:max(pre1$essai),.combine=rbind)%dopar%{
library(dplyr)
library(reshape2)
library(bipartite)
extinctions=NULL
pre1b=subset(pre1,essai==e)
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

prop_forbidden=mean(1-bidon$binary)

C=unique(bidon$C)
N=unique(bidon$N)
M=unique(bidon$M)

plant_list=unique(bidon$plant_species)
extinctions=rbind(extinctions,data.frame(plant_ext=NA,hummingbird_pers=length(unique(bidon$hummingbird_species)),site=z,essai=e,r=r,
rank=0,barrier=bar,prop_forbidden=prop_forbidden,scenario=scenario,C=C,N=N,M=M))
for(j in 1:length(plant_list)){
ext=plant_list[j]
if(length(unique(bidon$plant_species))>1){
for(i in unique(bidon$hummingbird_species)){
bidon2=subset(bidon,hummingbird_species==i)
distances=sqrt((bidon2$comp[bidon2$plant_species==ext]-bidon2$comp[bidon2$plant_species!=ext])^2)
R=mean(bidon2$average_proba[bidon2$plant_species!=ext]*exp(-r*distances/max(distances)))

Pext=(1-R)*bidon2$inter[bidon2$plant_species==ext]/sum(bidon2$inter)

EXT=rbinom(1,1,Pext)

if(EXT==1){bidon=subset(bidon,hummingbird_species!=i)}

}
}
bidon=subset(bidon,plant_species!=ext)
extinctions=rbind(extinctions,data.frame(plant_ext=ext,hummingbird_pers=length(unique(bidon$hummingbird_species)),site=z,essai=e,r=r,
rank=j,barrier=bar,prop_forbidden=prop_forbidden,scenario=scenario,C=C,N=N,M=M))


}
}
}
}
}
return(extinctions)
}

#LOAD SITE METADATA:
extinctionf=merge(extinctionf,sites,by="site")

fwrite(extinctionf,paste0("robustness_simulations_",pays,".txt"))
}

stopCluster(ncpus)
stopImplicitCluster()



#######################################################################################################################
extinctions=extinctions %>% group_by(r,site,barrier,essai,scenario) %>% mutate(nbh=max(hummingbird_pers),nbp=max(rank))
extinctions$persistence=extinctions$hummingbird_pers/extinctions$nbh
extinctions$rank2=extinctions$rank/extinctions$nbp
extinctions=extinctions[order(extinctions$min_transect_elev),]
extinctions$site=factor(extinctions$site,levels=unique(extinctions$site))
extinctions$site2=paste0(extinctions$site," (",extinctions$min_transect_elev,"m)")
extinctions$site2=factor(extinctions$site2,levels=unique(extinctions$site2))

b=extinctions %>% group_by(r,site,site2,barrier,rank2,min_transect_elev,nbh,nbp,scenario) %>% summarise(persistence=mean(hummingbird_pers/nbh),
pers_sde=sd(hummingbird_pers)/sqrt(length(hummingbird_pers)),prop_forbidden=mean(prop_forbidden))


ggplot()+
geom_line(data=b,aes(y=persistence,x=rank2,col=barrier,linetype=scenario),size=1)+facet_grid(rows=vars(site),cols=vars(r))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA))+
labs(col="")+ggtitle("a")+ylab("Persistence")+xlab("Fraction of plant species removed")+
scale_color_manual(values=c("#776274","#C8CC92"))+scale_fill_manual(values=c("#776274","#C8CC92"))+
scale_x_continuous(n.breaks=3,labels = scales::number_format(accuracy = 0.1))

pl1=ggplot()+
geom_line(data=subset(extinctions,scenario=="generalists first" & r==1),
aes(x=rank2,y=persistence,col=barrier,group=paste0(essai,r,site,barrier,scenario)),alpha=0.1)+
geom_line(data=subset(b,scenario=="generalists first" & r==1),aes(y=persistence,x=rank2,col=barrier),size=1.5)+facet_wrap(~site2)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none")+
labs(col="")+ggtitle("a")+ylab("Persistence")+xlab("Fraction of plant species removed")+
scale_color_manual(values=c("#30343F","#89023E"))+scale_fill_manual(values=c("#30343F","#89023E"))+
scale_x_continuous(n.breaks=2,labels = scales::number_format(accuracy = 1))


###AREA
b2=subset(b,rank2!=0 & rank2<1) %>% group_by(site,site2,barrier,r,min_transect_elev,nbh,nbp,scenario) %>%
summarise(robustness=sum(persistence)/length(persistence))

model=glm(robustness~poly(min_transect_elev,2)*barrier*as.factor(r)*as.factor(scenario),data=b2,family=quasibinomial)
AIC(model)
pre=ggpredict(model,c("min_transect_elev[605:3110]","barrier","r","scenario"))
names(pre)[names(pre) %in% c("group","facet","panel")]=c("barrier","r","scenario")
pre=subset(pre,r!=0 | barrier!="no forbidden links")

figs=ggplot()+
geom_ribbon(data=pre,aes(x=x,y=predicted,ymin=conf.low,ymax=conf.high,fill=barrier),alpha=0.2)+
geom_line(data=pre,aes(x=x,y=predicted,col=barrier),size=1.3)+
geom_point(data=b2 ,aes(y=robustness,x=min_transect_elev,col=barrier,fill=barrier),size=1.5)+stat_smooth(alpha=0.2)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),strip.background=element_rect(fill=NA,color=NA))+
labs(col="",fill="")+ylab("Robustness")+xlab("Elevation (meters above sea level)")+
scale_color_manual(values=c("#30343F","#89023E"))+scale_fill_manual(values=c("#30343F","#89023E"))+
facet_grid(rows=vars(r),cols=vars(scenario),scales="free",labeller = label_bquote(rows=alpha == .(r)))

setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper")
png("Figure_S2.png", width=1200,height=1000,res=150)
figs
dev.off();

model=glm(robustness~poly(min_transect_elev,2)*barrier,data=subset(b2,scenario=="generalists first" & r==1),family=quasibinomial)
AIC(model)
pre=ggpredict(model,c("min_transect_elev[605:3110]","barrier"))

pl2=ggplot()+
geom_ribbon(data=pre,aes(x=x,y=predicted,ymin=conf.low,ymax=conf.high,fill=group),alpha=0.2)+
geom_line(data=pre,aes(x=x,y=predicted,col=group),size=1.3)+
geom_point(data=subset(b2,scenario=="generalists first" & r==1) ,aes(y=robustness,x=min_transect_elev,col=barrier,fill=barrier),size=1.5)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+
labs(col="",fill="")+ggtitle("b")+ylab("Robustness")+xlab("Elevation (meters above sea level)")+
scale_color_manual(values=c("#30343F","#89023E"))+scale_fill_manual(values=c("#30343F","#89023E"))


grid.arrange(pl1,pl2,ncol=2)
setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper")
pdf("Figure_4.pdf", width=10,height=6)
grid.arrange(pl1,pl2,ncol=2)
dev.off();

######################################## DETERMINISM OF Robustness
deter=merge(b2,forb,by=c("site","min_transect_elev","nbp","nbh"))
deter=subset(deter,scenario=="generalists first" & r==1 & barrier=="with forbidden links")
deter$symmetrie=deter$nbp/deter$nbh
pairs(deter[,c("min_transect_elev","prop_forbidden","compl","nbh","nbp","symmetrie")])
cor(deter[,c("min_transect_elev","prop_forbidden","compl","nbh","nbp","symmetrie")])


model=randomForest(robustness~prop_forbidden+compl+nbp+nbh+symmetrie,data=deter,importance=T,ntree=4000)
varImpPlot(model)
partialPlot(model, pred.data=deter,x.var="prop_forbidden")
partialPlot(model, pred.data=deter,x.var="compl")

