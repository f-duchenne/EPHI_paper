###########################################
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
library(spatialEco)
library(MuMIn)
library(INLA)
library(inlabru)
library("inlaVP")
library(R2jags)
setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper/data")

tab=fread("data_for_analyses.txt")

#### CALCULATE TRAIT MATCHING:
tab$complementarity=tab$Tube_length-tab$bill_length
hist(tab$complementarity)
tab$barrier=0
tab$barrier[tab$complementarity>(0.25*tab$bill_length)]=1
hist(tab$barrier)
tab$mismatch=abs(tab$complementarity) #trait mismatch, scaled for convergence
tab$mismatch2=abs(1-(tab$Tube_length-tab$bill_length)/tab$bill_length) #trait mismatch, scaled for convergence
tab$elev=scale(tab$min_transect_elev)

### CALCULATE PHENOLOGICAL MATCH:
#hummingbird phenologies
phenoh=subset(tab, !is.na(hummingbird_species)) %>% group_by(hummingbird_species,month) %>% summarise(Ys=sum(Yfreq))
phenoh_mat=dcast(phenoh,month~hummingbird_species,value.var="Ys",fill=0)
func=function(x){
t(t(x[,-1])/apply(x[,-1],2,sum))
}
phenoh_mat[,-1]=func(phenoh_mat)
phenoh=melt(phenoh_mat,id.vars="month",variable.name="hummingbird_species")

#plant phenologies
phenop=subset(tab, !is.na(plant_species)) %>% group_by(plant_species,month) %>% summarise(Ys=sum(Yfreq))
phenop_mat=dcast(phenop,month~plant_species,value.var="Ys",fill=0)
phenop_mat[,-1]=func(phenop_mat)
phenop=melt(phenop_mat,id.vars="month",variable.name="plant_species")

#merge both
pheno=merge(phenoh,phenop,by="month")
#calculate match of both
pheno_match=pheno %>% group_by(hummingbird_species,plant_species) %>% summarise(pheno_match=sum(value.x*value.y))
pheno_match$pheno_barrier=0
pheno_match$pheno_barrier[pheno_match$pheno_match==0]=1
#merge pheno_match with data
tab=merge(tab,pheno_match,by=c("hummingbird_species","plant_species"),all.x=T,all.y=F)
tab=merge(tab,phenoh,by=c("hummingbird_species","month"),all.x=T,all.y=F)

#REMOVE DATA WITHOUT TRAITS
tab=subset(tab,!is.na(complementarity))

#PREPARE DATA FOR TEMPORAL AND SPATIAL AUTOCORRELATION:
tab$site <- factor(tab$site)
tab$group <- factor(rep(1, nrow(tab)))
tab$hummingbird_species_site=paste0(tab$hummingbird_species,tab$site)

#EVENTUALLY ROUND INTERACTION FREQUENCIES
tab$Yfreq_r=round(tab$Yfreq)

#convert factors to numerical indices for JAGS:
tab$hummingbird_species_site_num=as.numeric(as.factor(tab$hummingbird_species_site))
tab$hummingbird_num=as.numeric(as.factor(tab$hummingbird_species))
tab$plant_num=as.numeric(as.factor(tab$plant_species))
tab$site_num=as.numeric(as.factor(tab$site))
tab$num_time=as.numeric(as.factor(tab$num_time)) #rescale this factor to start at one

#create a table with bill length for hummingbirds
tabu=unique(tab[,c("hummingbird_num","bill_length")])
names(tabu)=paste0(names(tabu),"u")

sites=unique(tab[,c("site","site_num","midpoint_Longitude","midpoint_Latitude")])
Distance=dist(sites[,c("midpoint_Longitude","midpoint_Latitude")], method="euclidean", diag=TRUE, upper=TRUE)

#assemble data
dat=c(as.list(tab),as.list(tabu),list(N=nrow(tab),
Nbirds=length(unique(tab$hummingbird_num)),Nbird_site=length(unique(tab$hummingbird_species_site_num)),Nsites=length(unique(tab$site_num)),Nplants=length(unique(tab$plant_num)),
Ntemp=length(unique(tab$num_time)),Nsite=length(unique(tab$site_num)),Distance=as.matrix(Distance)))

save(dat,file=paste0("data_formodel.RData"))

######################################################################################################
######################################################################################################
#### JAGS model file written by runjags version 2.2.0-2 on 2021-11-23 15:22:45 
######################################################################################################
######################################################################################################

### Model template as follows - ensure this is syntactically correct before running the model!

model_string="
model{

# Barrier assessment
for(j in 1:Nbirds){
barrier_infer[j] ~ dnorm(bill_lengthu[j],1/(0.2*bill_lengthu[j]*0.2*bill_lengthu[j]))T(bill_lengthu[j],2*bill_lengthu[j])
#barrier_infer[j] ~ dunif(bill_lengthu[j],2*bill_lengthu[j])
match_infer[j] ~ dnorm(bill_lengthu[j],1/(0.2*bill_lengthu[j]*0.2*bill_lengthu[j]))T(0.5*bill_lengthu[j],1.5*bill_lengthu[j])
#match_infer[j] ~ dunif(0.5*bill_lengthu[j],1.5*bill_lengthu[j])
}

# Process model
for(i in 1:N){
	#new variables:
	barrier_var[i] <- ifelse(Tube_length[i] > barrier_infer[hummingbird_num[i]], 1, 0)
	mismatch_var[i] <- abs(Tube_length[i]-match_infer[hummingbird_num[i]])
	
	# Zero inflation:
	logit(pz[i]) <- Interceptpz + traitBarrier * barrier_var[i]
	Z[i] ~ dbern(min(pz[i]+0.00000000000000001,0.9999999999999999))
	
	# Interaction frequency:
	log(lambda[i]) <- Intercept + traitMismatch * mismatch_var[i] + pheno * value[i] +plant_effect[plant_num[i]]+site_effect[site_num[i]]+temp_effect[site_num[i],num_time[i]]+sitebird_effect[hummingbird_species_site_num[i]]
	p[i] <- r/(r+(lambda[i] * Z[i]))
	Yfreq_r[i] ~ dnegbin(min(p[i]+0.00000000000000001,0.9999999999999999),r)
	#Yfreq_r[i] ~ dpois(lambda[i] * Z[i] + 0.00000000000000001)
	
}

# PRIORS FIXED EFFECTS:
Interceptpz ~ dnorm(0, 0.5)
Intercept ~ dnorm(0, 0.01) 
traitMismatch ~ dnorm(0, 0.01)
pheno ~ dnorm(0, 0.01)
traitBarrier ~ dnorm(0, 0.5)T(,0)

# PRIORS RANDOM EFFECTS:
#spatial autocorrelation:
edec ~ dgamma(1, 2)
for(j in 1:Nsite){
for(i in 1:Nsite){
D.covar[i,j] <- exp(-edec*Distance[i,j])
}}

site_effect[1:Nsite] ~ dmnorm.vcov(rep(0,Nsite), spatialvariance * D.covar[1:Nsite,1:Nsite])
spatialvariance=sd.site * sd.site

for(i in 1:Nsite){
#site_effect[i] ~ dnorm(0, tau.site)
temp_effect[i,1] <- 0
for(j in 2:Ntemp){
temp_effect[i,j] ~ dnorm(temp_effect[i,(j-1)], tau.temp)
}
}

for(j in 1:Nbird_site){
sitebird_effect[j] ~ dnorm(0, tau.bird)
}

for(j in 1:Nplants){
plant_effect[j] ~ dnorm(0, tau.plant)
}
tau.plant <- 1/(sd.plant * sd.plant)
sd.plant ~ dt(0, 1, 1)T(0,)

#hyperpriors:
tau.site<- 1/(sd.site * sd.site)
sd.site ~ dt(0, 1, 1)T(0,)
tau.temp<- 1/(sd.temp * sd.temp)
sd.temp ~ dt(0, 1, 1)T(0,)
tau.bird=1/(sd.bird * sd.bird)
sd.bird ~ dt(0, 1, 1)T(0,)


# PRIOR OVERDISPERTION:
r ~ dunif(0,50)

}
"

writeLines(model_string,con="model.txt")

ParsStage <- c("barrier","match","Interceptpz","traitBarrier","Intercept","traitMismatch","sitebird_effect","plant_effect","site_effect","temp_effect","sd.plant","sd.bird","sd.site","sd.temp","r")

Inits <- function(){list()}

t1=Sys.time()
results1 <- jags(model.file="model.txt", parameters.to.save=ParsStage, n.chains=1, data=dat,n.burnin = 50,  n.iter = 100, n.thin = 2,inits =Inits,jags.seed =2)
t2=Sys.time()
t2-t1

save(results1,file=paste0("chain_model_",j,".RData"))
#


setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper/data")
library(mcmcOutput)
library(mcmcplots)
library(MCMCvis)
library(ggridges)
library(runjags)
library(pastecs)

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
bidon=summary(mco)
bidon$varia=rownames(bidon)


tabu$barrier_infer=bidon$mean[grep("barrier_infer",bidon$varia)]
plot(barrier_infer~bill_lengthu,data=tabu)
abline(0,1)

tabu$match_infer=bidon$mean[grep("match_infer",bidon$varia)]
plot(match_infer~bill_lengthu,data=tabu)
abline(0,1)


predict_model=function(nsampling,site,bird,plant,trait_plant,mismatch,barrier,pheno,type,random_effects,month,year){

rm(list=c("barrier_var","mismatch_var"))

####SOME FUNCTIONS:
inv.logit=function(y){exp(y)/(1+exp(y))}

rand_comb=function(y,name){
resi=NULL
if(is.vector(y)){
for(i in 1:length(y)){
resi=rbind(resi,obj[x,paste0(name,"[",y[i],"]")])
}
return(resi)
}else{
for(i in 1:length(y)){
resi=rbind(resi,obj[x,paste0(name,"[",y[[1]][i],",",y[[2]][i],,"]")])
}
return(resi)
}
}

###CONTROL FOR BASIC ERRORS
if(nsampling>10000){stop("nsampling is higher than 10000")}
if(nsampling<100){stop("nsampling is lower than 100")}
if(!(type %in% c('barrier','frequency','total'))){stop("type of predict does correspond to 'barrier', 'frequency' or 'total'")}
if(!(type %in% c("site","plants","bird_site","temp"))){stop("random effects to include in predictions do correspond to 'site','plants','bird_site' or 'temp'")}


#convert discrete variables to numerical index
bidon=unique(dat[,c("hummingbird_num","hummingbird_species")])
if(!is.na(bird)){bird_num=bidon$hummingbird_num[bidon$hummingbird_species==bird]}
bidon=unique(dat[,c("plant_num","plant_species")])
if(!is.na(plant)){plant_num=bidon$plant_num[bidon$plant_species==plant]}
bidon=unique(dat[,c("site_num","site")])
if(!is.na(site)){site_num=bidon$site_num[bidon$site==site]}
bidon=unique(dat[,c("hummingbird_species_site_num","hummingbird_species_site")])
if(!is.na(site) & !is.na(bird)){bird_site_num=bidon$hummingbird_species_site_num[bidon$hummingbird_species_site==paste0(bird,site)]}
bidon=unique(dat[,c("num_time","month","year")])
if(!is.na(month) & !is.na(year)){temp_num=bidon$num_time[bidon$month==month & bidon$year==year]}

if(length(bird)==0){stop("The humming species does fit the list of the model, please check the name")}
if(length(pla)==0){stop("The plant species does fit the list of the model, please check the name or set plant to NA to neglect the random effect")}
if(length(site)==0){stop("The site does fit the list of the model, please check the name or set site to NA to neglect the random effect")}

x=sample(1:nrow(obj),nsampling,replace=T)

if(type %in% c("barrier","total")){
intercept_proba=obj[x,"Interceptpz"]


if(!is.na(barrier)){
barrier_var=t(replicate(nsampling,barrier))
print("Use provided values for barrier instead of latent ones")
}else{
barrier_var <- ifelse(trait_plant > rand_comb(bird_num,"barrier_infer"), 1, 0)
print("No values provided for barrier, use latent ones")
}

proba=inv.logit(intercept_proba+t(replicate(nsampling,barrier_var))*obj[x,paste0("traitBarrier")])
proba_suma=data.frame(site=site,hummingbird_species=bird,plant_species=plant,Tube_length=trait_plant,
average_proba=apply(proba,1,mean),median_proba=apply(proba,1,median),lwr_proba=apply(proba,1,quantile,prob=0.025),upr_proba=apply(proba,1,quantile,prob=0.975),
link_binary=rbinom(1,1,apply(proba,1,mean)))
}




if(type %in% c("frequency","total")){

intercept=obj[x,"Intercept"]

if(!is.na(mismatch)){
mismatch_var=t(replicate(nsampling,mismatch))
print("Use provided values for mismatch instead of latent ones")
}else{
mismatch_var <-  abs(trait_plant-rand_comb(bird_num,"match_infer"))
print("No values provided for mismatch, use latent ones")
}

lambda_log=intercept+obj[x,"traitMismatch"]*t(replicate(nsampling,match_var))+obj[x,"pheno"]*t(replicate(nsampling,pheno))+
			if(is.na(plant) | !("plant" %in% random_effects)){0}else{rand_comb(plant_num,"plant_effect")}+
			if(is.na(site) | !("site" %in% random_effects)){0}else{rand_comb(site_num,"site_effect")}+
			if(is.na(bird_site) | is.na(temp) | !("temp" %in% random_effects)){0}else{rand_comb(list(site_num,num_time),"site_effect")}+
			if(is.na(bird_site) | !("bird_site" %in% random_effects)){0}else{rand_comb(bird_site_num,"sitebird_effect")}


r=obj[x,"r"]

if(type=="total"){Z=matrix(sapply(proba,function(y){return(rbinom(1,1,y))}),ncol=ncol(proba))}else{Z=1}

p <- r/(r+(exp(lambda) * Z))

frequencies=matrix(mapply(function(y,z){rnbinom(1, y, z)},r,p),ncol=ncol(proba))

freq_suma=data.frame(site=site,hummingbird_species=bird,plant_species=plant,Tube_length=trait_plant,
average_freq=apply(frequencies,1,mean),median_freq=apply(frequencies,1,median),lwr_freq=apply(frequencies,1,quantile,prob=0.025),upr_freq=apply(frequencies,1,quantile,prob=0.975),
link_binary=rbinom(1,1,apply(frequencies,1,mean)))
}


if(type=="barrier"){
return(proba_suma)
}else{
return(freq_suma)
}

}

