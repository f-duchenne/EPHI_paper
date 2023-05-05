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
library(parallel)

load("data_formodel.RData")

#### conservative prior
load("chain_model_1.RData")
model1=results1
load("chain_model_2.RData")
model2=results1
load("chain_model_3.RData")
model3=results1
obj1=as.mcmc(model1)
obj2=as.mcmc(model2)
obj3=as.mcmc(model3)
#combine chains and get a summary:
mco <- mcmcOutput(mcmc.list(obj1,obj2,obj3))
suma=summary(mco)
suma$varia=rownames(suma)

#### L1 Barrier threshold
tabu=as.data.frame(dat[c("hummingbird_numu","bill_lengthu")])
tabu$barrier_infer=suma$mean[grep("barrier_infer",suma$varia)]
#### L2 Match optimum
tabu$match_infer=suma$mean[grep("match_infer",suma$varia)]
tabu=merge(tabu,unique(as.data.frame(dat[c("hummingbird_species","hummingbird_num")])),by.x="hummingbird_numu",by.y="hummingbird_num")

##### NETWORK FOR ROBUSTESSE
plant_res=fread("plants_per_site_per_month.txt")
plant_res=subset(plant_res,!is.na(Tube_length))
dat$elev2=dat$elev*attr(dat$elev,"scaled:scale")+attr(dat$elev,"scaled:center")
sitebird=unique(as.data.frame(dat[c("site","hummingbird_species","month","phenoh","elev2")]))
tab2=merge(plant_res,sitebird,by=c("site","month"),allow.cartesian=TRUE)
tab2=tab2 %>% group_by(site,plant_species,hummingbird_species,Tube_length,abond_flower_moy) %>% summarise(pheno_predict=sum(phenoh*phenop)/12)

pre1=predict_model(chains=mco,nsampling=5000,site=tab2$site,bird=tab2$hummingbird_species,plant=tab2$plant_species,trait_plant=tab2$Tube_length,mismatch=NA,barrier=NA,
pheno=tab2$pheno_predict,abond=tab2$abond_flower_moy,type="network",
random_effects=c("site","plant"),month=NA,year=NA,nb_net=100)

pre1$inter=pre1$freq*pre1$binary

pre1=merge(pre1,tabu,by="hummingbird_species")
pre1$comp=-abs(pre1$Tube_length-pre1$match_infer)
EPHI_version="2023-02-14"

pre1$C=NA
pre1$N=NA
pre1$M=NA
for(i in unique(pre1$site)){
bidon=subset(pre1,site==i)
for(j in unique(bidon$essai)){
bidonb=subset(bidon,essai==j)
mat=dcast(bidonb,plant_species~hummingbird_species,value.var="inter")
mat=as.matrix(mat[,-1])
pre1$C[pre1$site==i & pre1$essai==j]=networklevel(mat,index="weighted connectance",weighted=T)
pre1$N[pre1$site==i & pre1$essai==j]=wine(mat, nreps = 1000)$wine
pre1$M[pre1$site==i & pre1$essai==j]=computeModules(mat, method="Beckett")@likelihood[1]
}}


sites_cr=fread(paste0("C:/Users/Duchenne/Documents/EPHI_data_clean/Costa-Rica_",EPHI_version,"/Site_metadata_Costa-Rica.txt"),na.strings = c("",NA))
pre1b=merge(pre1,sites_cr,by="site")
fwrite(pre1b,"initial_network.txt")


ncpus <- makeCluster(5)
registerDoParallel(ncpus)

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
if(bar=="no forbidden links"){bidon$binary=1}
bidon$inter=bidon$freq*bidon$binary
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
R=mean(bidon2$binary[bidon2$plant_species!=ext]*exp(-r*distances/max(distances)))

Pext=(1-R)*bidon2$freq[bidon2$plant_species==ext]/sum(bidon2$freq)

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

stopCluster(ncpus)
stopImplicitCluster()

#LOAD SITE METADATA:
extinctionf=merge(extinctionf,sites_cr,by="site")

fwrite(extinctionf,"robustness_simulations.txt")

########################################################################
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
library(emmeans)
setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper/data")

######################## TRAITS
ini_net=fread("initial_network.txt")
extinctions=fread("robustness_simulations.txt")

bp=unique(ini_net[,c("site","plant_species","Tube_length","min_transect_elev")])
bp$type="plants"
names(bp)[2:3]=c("species","trait")
bh=unique(ini_net[,c("site","hummingbird_species","bill_lengthu","min_transect_elev")])
bh$type="birds"
names(bh)[2:3]=c("species","trait")

bf=rbind(bh,bp)
bf$type=as.factor(bf$type)

model=lmer(trait~min_transect_elev*type+(1|site),data=bf)
emtrends(model,specs="type",var="min_transect_elev")

pre=ggpredict(model,c("min_transect_elev","type"))

pl1=ggplot()+
geom_jitter(data=bf,aes(x=min_transect_elev,y=trait,col=type,fill=type),alpha=0.3,width = 20)+
geom_ribbon(data=pre,aes(x=x,y=predicted,ymin=conf.low,ymax=conf.high,fill=group),alpha=0.2)+
geom_line(data=pre,aes(x=x,y=predicted,col=group,linetype=group),size=1.3)+
xlab("Elevation (m)")+ylab("Corolla or Bill length (mm)")+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+
scale_color_manual(values=c("darkorchid4","#C8CC92"))+scale_fill_manual(values=c("darkorchid4","#C8CC92"))+ggtitle("a")+labs(color="",fill="")+
scale_linetype_manual(values=c("dashed","solid"))+guides(linetype="none")

pl1b=ggplot()+
geom_density(data=bf,aes(x=trait,col=min_transect_elev,fill=min_transect_elev,group=min_transect_elev),alpha=0.1)+
xlab("Bill or Corolla length(mm)")+ylab("Density")+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),strip.background = element_blank())+
ggtitle("b")+labs(color="Elev.",fill="Elev.")+
facet_wrap(~type)+scale_color_viridis(option="mako")+scale_fill_viridis(option="mako")+coord_cartesian(expand=F)

####Panel b
forb=as.data.frame(ini_net %>% group_by(site,min_transect_elev) %>% summarise(compl=mean(-1*abs(trait_plant-match_infer)),prop_forbidden=mean(1-average_proba),
nbh=length(unique(hummingbird_species)),nbp=length(unique(plant_species))))

model=lm(compl~min_transect_elev,data=forb)
Anova(model)
pre=ggpredict(model,"min_transect_elev[all]")

pl2=ggplot()+
geom_point(data=forb,aes(y=compl,x=min_transect_elev),size=1.5)+
geom_ribbon(data=pre,aes(x=x,y=predicted,ymin=conf.low,ymax=conf.high),alpha=0.2,fill="lightgrey")+
geom_line(data=pre,aes(x=x,y=predicted),size=1.3,linetype="dashed")+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+
labs(col="",fill="")+ggtitle("b")+ylab("Average trait complementarity")+xlab("Elevation (m)")

bidon=melt(forb,id.vars=c("site","min_transect_elev","compl","prop_forbidden"))
gammodel=gam(value~s(min_transect_elev),data=subset(bidon,variable=="nbh"),family=poisson)
pre1=as.data.frame(predict(gammodel,se=TRUE,type="response"))
gammodel=gam(value~s(min_transect_elev),data=subset(bidon,variable=="nbp"),family=poisson)
pre2=as.data.frame(predict(gammodel,se=TRUE,type="response"))
bidon=cbind(bidon,rbind(pre1,pre2))
bidon$variable=as.character(bidon$variable)
bidon$variable[bidon$variable=="nbp"]="Plants"
bidon$variable[bidon$variable=="nbh"]="H."

pl2b=ggplot(data=bidon,aes(y=value,x=min_transect_elev,col=variable,fill=variable))+
geom_ribbon(aes(ymin=fit-1.96*se.fit,ymax=fit+1.96*se.fit),alpha=0.2,col=NA)+
geom_line(aes(y=fit))+
geom_point(size=1.5)+
xlab("Elevation (m)")+ylab("Species richness")+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+
scale_color_manual(values=c("darkorchid4","#C8CC92"))+scale_fill_manual(values=c("darkorchid4","#C8CC92"))+ggtitle("a")+labs(color="",fill="")


####Panel c
model=glm(prop_forbidden~min_transect_elev,data=forb,family=quasibinomial)

pre=ggpredict(model,c("min_transect_elev [605:3110 by=1]"))

pl3=ggplot()+
geom_point(data=forb,aes(y=prop_forbidden,x=min_transect_elev),size=1.5)+
geom_ribbon(data=pre,aes(x=x,y=predicted,ymin=conf.low,ymax=conf.high),alpha=0.2,fill="lightgrey")+
geom_line(data=pre,aes(x=x,y=predicted),size=1.3)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+
labs(col="",fill="")+ggtitle("c")+ylab("Propotion of forbidden links")+xlab("Elevation (m)")

setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper")
pdf("Figure3.pdf",width=9,height=3)
grid.arrange(pl1,pl2,pl3,ncol=3,widths=c(1.5,1,1))
dev.off();

setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper")
png("FigureS1.png",width=1500,height=800,res=160)
grid.arrange(pl2b,pl1b,ncol=2,widths=c(1,1.5))
dev.off();


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

