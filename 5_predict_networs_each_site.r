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
EPHI_version="2023-07-04"

for(pays in c("Costa-Rica","Ecuador","Brazil")){
load(paste0("data_formodel_",pays,".RData"))
plant_res=fread(paste0("plants_per_site_per_month_",pays,".txt"))
plant_res=subset(plant_res,!is.na(Tube_length))

if(pays=="Ecuador"){
load(paste0("chain_model_",pays,"_1.RData"))
model1=results1
load(paste0("chain_model_",pays,"_2.RData"))
model2=results1
obj1=as.mcmc(model1)
obj2=as.mcmc(model2)
#combine chains and get a summary:
mco <- mcmcOutput(mcmc.list(obj1,obj2))
}else{
#### conservative prior
load(paste0("chain_model_",pays,"_1.RData"))
model1=results1
load(paste0("chain_model_",pays,"_2.RData"))
model2=results1
load(paste0("chain_model_",pays,"_3.RData"))
model3=results1
obj1=as.mcmc(model1)
obj2=as.mcmc(model2)
obj3=as.mcmc(model3)
#combine chains and get a summary:
mco <- mcmcOutput(mcmc.list(obj1,obj2,obj3))
}
suma=summary(mco)
suma$varia=rownames(suma)

#### L1 Barrier threshold
tabu=as.data.frame(dat[c("hummingbird_numu","culmen_lengthu")])
tabu$barrier_infer=suma$mean[grep("barrier_infer",suma$varia)]
#### L2 Match optimum
tabu$match_infer=suma$mean[grep("match_infer",suma$varia)]
tabu=merge(tabu,unique(as.data.frame(dat[c("hummingbird_species","hummingbird_num")])),by.x="hummingbird_numu",by.y="hummingbird_num")

##### NETWORK FOR ROBUSTESSE
dat$elev2=dat$elev*attr(dat$elev,"scaled:scale")+attr(dat$elev,"scaled:center")
sitebird=unique(as.data.frame(dat[c("site","hummingbird_species","month","phenoh","elev2")]))
tab2=merge(plant_res,sitebird,by=c("site","month"),allow.cartesian=TRUE)
tab2=tab2 %>% group_by(site,plant_species,hummingbird_species,Tube_length,abond_flower_moy) %>% summarise(pheno_predict=sum(phenoh*phenop)/12)

pre1=predict_model(chains=mco,nsampling=5000,site=tab2$site,bird=tab2$hummingbird_species,plant=tab2$plant_species,trait_plant=tab2$Tube_length,mismatch=NA,barrier=NA,
pheno=tab2$pheno_predict,abond=tab2$abond_flower_moy,type="network",
random_effects=c("site","plant"),month=NA,year=NA,nb_net=100,duration=NA)

pre1$inter=pre1$freq*pre1$binary

pre1=merge(pre1,tabu,by="hummingbird_species")
pre1$comp=-abs(pre1$Tube_length-pre1$match_infer)


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

sites=fread(paste0("C:/Users/Duchenne/Documents/EPHI_data_clean/",pays,"_",EPHI_version,"/Site_metadata_",pays,".txt"),na.strings = c("",NA))
sites=subset(sites,habitat!="deforested")
pre1b=merge(pre1,sites,by="site")
pre1b$Country=pays
fwrite(pre1b,paste0("initial_network_",pays,".txt"))
}




