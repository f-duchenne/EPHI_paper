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

for(pays in c("Costa-Rica","Ecuador")){

tab=fread(paste0("data_for_analyses_",pays,".txt"),na.string=c("",NA))

#### CALCULATE TRAIT MATCHING:
tab$mismatch=abs(tab$Tube_length-tab$bill_length) #trait mismatch, scaled for convergence
tab$elev=scale(tab$min_transect_elev)

#### CALCULATE PHENOLOGICAL MATCHING:
tab$phenomatch=tab$phenoh*tab$phenop

#ABONDANCE FLOWER THAT ARE NA ARE NON DETECTED PLANT ON THE TRANSECT, set a low abundance value:
tab$abond_flower[is.na(tab$abond_flower) | tab$abond_flower==0]=1
tab$abond_flower_log=log(tab$abond_flower)

#REMOVE DATA WITHOUT TRAITS
tab=subset(tab,!is.na(mismatch))

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
tabu=tabu[order(tabu$hummingbird_numu),]

sites=unique(tab[,c("site","site_num","midpoint_Longitude","midpoint_Latitude")])
sites=sites[order(sites$site_num),]
Distance=dist(sites[,c("midpoint_Longitude","midpoint_Latitude")], method="euclidean", diag=TRUE, upper=TRUE)
Distance=as.matrix(Distance)/max(Distance)

#assemble data
dat=c(as.list(tab),as.list(tabu),list(N=nrow(tab),
Nbirds=length(unique(tab$hummingbird_num)),Nbird_site=length(unique(tab$hummingbird_species_site_num)),Nsites=length(unique(tab$site_num)),Nplants=length(unique(tab$plant_num)),
Ntemp=length(unique(tab$num_time)),Nsite=length(unique(tab$site_num)),Distance=Distance))

save(dat,file=paste0("data_formodel_",pays,".RData"))
}

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
match_infer[j] ~ dnorm(bill_lengthu[j],1/(0.2*bill_lengthu[j]*0.2*bill_lengthu[j]))T(0.5*bill_lengthu[j],1.5*bill_lengthu[j])
#match_infer[j] ~ dunif(0.5*bill_lengthu[j],1.5*bill_lengthu[j])
barrier_infer[j] ~ dnorm(bill_lengthu[j],1/(0.3*bill_lengthu[j]*0.3*bill_lengthu[j]))T(max(bill_lengthu[j],match_infer[j]),2*bill_lengthu[j])
#barrier_infer[j] ~ dunif(bill_lengthu[j],2*bill_lengthu[j])
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
	log(lambda[i]) <- Intercept + traitMismatch * mismatch_var[i] + pheno * phenoh[i]+abond *abond_flower_log[i] +plant_effect[plant_num[i]]+site_effect[site_num[i]]+temp_effect[site_num[i],num_time[i]]+sitebird_effect[hummingbird_species_site_num[i]]
	p[i] <- r/(r+(lambda[i] * Z[i]))
	Yfreq_r[i] ~ dnegbin(min(p[i]+0.00000000000000001,0.9999999999999999),r)
	#Yfreq_r[i] ~ dpois(lambda[i] * Z[i] + 0.00000000000000001)
	
}

# PRIORS FIXED EFFECTS:
Interceptpz <- 10
Intercept ~ dnorm(0, 0.01) 
traitMismatch ~ dnorm(0, 0.01)T(,0)
pheno ~ dnorm(0, 0.01)
abond ~ dnorm(0, 0.01)
traitBarrier ~ dunif(-20,0)

# PRIORS RANDOM EFFECTS:
#spatial autocorrelation:
edec ~ dgamma(3, 0.1)
for(j in 1:Nsite){
for(i in 1:Nsite){
D.covar[i,j] <- exp(-edec*Distance[i,j])
}}

site_effect[1:Nsite] ~ dmnorm.vcov(rep(0,Nsite), spatialvariance * D.covar[1:Nsite,1:Nsite])
spatialvariance=sd.site * sd.site

for(i in 1:Nsite){
#site_effect[i] ~ dnorm(0, tau.site)
temp_effect[i,1] ~ dnorm(0, tau.temp)
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
r ~ dnegbin(0.2,4)

}
"

writeLines(model_string,con="model.txt")






ParsStage <- c("barrier_infer","match_infer","Interceptpz","traitBarrier","Intercept","traitMismatch","pheno","abond","plant_effect","site_effect","temp_effect","sitebird_effect","sd.plant","sd.bird","sd.site","sd.temp","r","edec")

Inits <- function(){list()}

t1=Sys.time()
results1 <- jags(model.file="model.txt", parameters.to.save=ParsStage, n.chains=1, data=dat,n.burnin = 50,  n.iter = 100, n.thin = 2,inits =Inits,jags.seed =2)
t2=Sys.time()
t2-t1

save(results1,file=paste0("chain_model_",j,".RData"))
#


