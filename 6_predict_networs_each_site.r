###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("randomForest","data.table", "dplyr", "lubridate","lme4","R2jags","mcmcOutput","mcmcplots","MCMCvis",
			"pastecs","ggplot2","cowplot","gridExtra","scales","reshape2","bipartite","stringr","ungeviz") 


inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)

EPHI_version="2024-01-30"

setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper/data")
source("C:/Users/Duchenne/Documents/EPHI_paper/scripts/function_to_predict_from_bayesian.r")


for(pays in c("Costa-Rica","Ecuador","Brazil")){
load(paste0("data_formodel_",pays,".RData"))
plant_res=fread(paste0("plants_per_site_per_month_",pays,".txt"))
plant_res=subset(plant_res,!is.na(Tubelength))

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
tab2=tab2 %>% group_by(site,plant_species,hummingbird_species,Tubelength,abond_flower_moy) %>%
summarise(pheno_predict=sum(phenoh*phenop,na.rm=T)/12)

pre1=predict_model(data=dat,chains=mco,nsampling=5000,site=tab2$site,bird=tab2$hummingbird_species,
plant=tab2$plant_species,trait_plant=tab2$Tubelength,mismatch=NA,barrier=NA,
pheno=tab2$pheno_predict,abond=rep(10,nrow(tab2)),type="network",
random_effects=c("site","plant"),month=NA,year=NA,nb_net=1,duration=rep(1000,nrow(tab2)))
pre1$barrier="with forbidden links"
pre1$inter=pre1$freq*pre1$average_proba

pre2=pre1
pre2$barrier="no forbidden links"
pre2$inter=pre2$freq*pre2$average_proba_without_barrier

pre1=rbind(pre2,pre1)

pre1=merge(pre1,tabu,by="hummingbird_species")
pre1$compl=-abs(pre1$Tube_length-pre1$match_infer)

pre1$C=NA
pre1$N=NA
pre1$M=NA
pre1$Cperso=NA
pre1$Cperso2=NA
for(i in unique(pre1$site)){
bidon=subset(pre1,site==i)
for(j in unique(bidon$barrier)){
bidonb=subset(bidon,barrier==j)
mat=dcast(bidonb,plant_species~hummingbird_species,value.var="inter")
mat=as.matrix(mat[,-1])
pre1$C[pre1$site==i & pre1$barrier==j]=networklevel(mat,index="weighted connectance",weighted=T)
pre1$N[pre1$site==i & pre1$barrier==j]=wine(mat, nreps = 1000)$wine
pre1$M[pre1$site==i & pre1$barrier==j]=computeModules(mat, method="Beckett")@likelihood[1]
pre1$Cperso[pre1$site==i & pre1$barrier==j]=length(mat[mat>0.01*sum(mat)])/length(mat)
vec=sort(c(mat),decreasing=T)
vec1=vec[1:length(cumsum(vec)[cumsum(vec)<0.95*sum(vec)])]
pre1$Cperso2[pre1$site==i & pre1$barrier==j]=length(vec1)/length(mat)
}}

sites=fread(paste0("C:/Users/Duchenne/Documents/EPHI_data_clean/",pays,"_",EPHI_version,"/Site_metadata_",pays,".txt"),na.strings = c("",NA))
sites=subset(sites,habitat!="deforested")
pre1b=merge(pre1,sites,by="site")
pre1b$Country=pays
fwrite(pre1b,paste0("initial_network_",pays,".txt"))
}

