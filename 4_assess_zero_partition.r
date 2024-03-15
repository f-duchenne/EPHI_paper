###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("randomForest","data.table", "dplyr", "lubridate","lme4","R2jags","mcmcOutput","mcmcplots","MCMCvis") 


inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)

EPHI_version="2024-01-30"

setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper/data")

res_partition=NULL
for(pays in c("Costa-Rica","Ecuador","Brazil")){

load(paste0("data_formodel_",pays,".RData"))
plant_res=fread(paste0("plants_per_site_per_month_",pays,".txt"))
plant_res=subset(plant_res,!is.na(Tubelength))
sites=fread(paste0("C:/Users/Duchenne/Documents/EPHI_data_clean/",pays,"_",EPHI_version,"/Site_metadata_",pays,".txt"),na.strings = c("",NA))


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

for (dure in c(10,100,250,500,750,1000,2000)){
pre1=predict_model(data=dat,chains=mco,nsampling=5000,site=tab2$site,bird=tab2$hummingbird_species,
plant=tab2$plant_species,trait_plant=tab2$Tubelength,mismatch=NA,barrier=NA,
pheno=tab2$pheno_predict,abond=rep(10,nrow(tab2)),type="network",
random_effects=c("site","plant","bird_site"),month=NA,year=NA,nb_net=1,duration=rep(dure,nrow(tab2)))

pre1$zero_freq=dnbinom(0,size=unique(pre1$size_param),mu=pre1$freq)

pre1=pre1 %>% group_by(site) %>%
mutate(pzero_infl=1-average_proba_without_barrier,
pzero_barr=(average_proba_without_barrier-average_proba),
pzero_freq=zero_freq*average_proba,
pzero=pzero_freq+(1-average_proba))

b=pre1 %>% group_by(site) %>% summarise(pzerom=mean(pzero),unknown=mean(pzero_infl),
barrier=mean(pzero_barr),complementarity=mean(pzero_freq))

b$tot=apply(b[,3:5],1,sum)
b$duration=dure
b$Country=pays
b=merge(b,sites[,c("site","min_transect_elev")],by="site")
res_partition=rbind(res_partition,b)
}

}

fwrite(res_partition,"sparsity_estimates.txt")

###############