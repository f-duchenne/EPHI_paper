###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
pkgs <- c("randomForest","data.table", "dplyr", "lubridate","lme4","R2jags","here") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)

here::i_am("EPHI_paper.Rproj")

for(pays in c("Costa-Rica","Ecuador","Brazil")){

tab=fread(here("data_zenodo",paste0("data_for_analyses_",pays,".txt")),na.string=c("",NA))

#### CALCULATE TRAIT MATCHING:
tab$mismatch=abs(tab$Tubelength-tab$culmen_length) #trait mismatch
tab$elev=scale(tab$min_transect_elev)

#### CALCULATE PHENOLOGICAL MATCHING:
tab$phenomatch=tab$phenoh*tab$phenop

#ABONDANCE FLOWER THAT ARE NA ARE NON DETECTED PLANT ON THE TRANSECT, set a low abundance value:
tab$abond_flower[is.na(tab$abond_flower) | tab$abond_flower==0]=1
tab$abond_flower_log=log(tab$abond_flower)

unique(tab$site)
#REMOVE DATA WITHOUT TRAITS
tab=subset(tab,!is.na(mismatch))
unique(tab$site)

#PREPARE DATA FOR TEMPORAL AND SPATIAL AUTOCORRELATION:
tab$site <- factor(tab$site)
tab$group <- factor(rep(1, nrow(tab)))
tab$hummingbird_species_site=paste0(tab$hummingbird_species,tab$site)

#convert factors to numerical indices for JAGS:
tab$hummingbird_species_site_num=as.numeric(as.factor(tab$hummingbird_species_site))
tab$hummingbird_num=as.numeric(as.factor(tab$hummingbird_species))
tab$plant_num=as.numeric(as.factor(tab$plant_species))
tab$site_num=as.numeric(as.factor(tab$site))
tab$num_time=as.numeric(as.factor(tab$num_time)) #rescale this factor to start at one

#create a table with bill length for hummingbirds
tabu=unique(tab[,c("hummingbird_num","culmen_length")])
names(tabu)=paste0(names(tabu),"u")
tabu=tabu[order(tabu$hummingbird_numu),]

#create a table with corolla length for plants
# tabp=unique(tab[,c("plant_num","Tubelength")])
# names(tabp)=paste0(names(tabp),"p")
# tabp=tabp[order(tabu$plant_nump),]

sites=unique(tab[,c("site","site_num","midpoint_Longitude","midpoint_Latitude")])
sites=sites[order(sites$site_num),]
Distance=dist(sites[,c("midpoint_Longitude","midpoint_Latitude")], method="euclidean", diag=TRUE, upper=TRUE)
Distance=as.matrix(Distance)/max(Distance)

#assemble data
dat=c(as.list(tab),as.list(tabu),list(N=nrow(tab),
Nbirds=length(unique(tab$hummingbird_num)),Nbird_site=length(unique(tab$hummingbird_species_site_num)),Nsites=length(unique(tab$site_num)),Nplants=length(unique(tab$plant_num)),
Ntemp=length(unique(tab$num_time)),Nsite=length(unique(tab$site_num)),Distance=Distance))

save(dat,file=paste0(here("data_zenodo"),"/data_formodel_",pays,".RData"))
}

######################################################################################################
#####################################################################dat#################################
#### JAGS model file written by runjags version 2.2.0-2 on 2021-11-23 15:22:45 
######################################################################################################
######################################################################################################

model_string="
model{

# Inferring barrier and optimal corolla length for each hummingbird species
for(j in 1:Nbirds){
match_infer[j] ~ dnorm(culmen_lengthu[j],1/(0.2*culmen_lengthu[j]*0.2*culmen_lengthu[j]))T(0.5*culmen_lengthu[j],1.5*culmen_lengthu[j])
#match_infer[j] ~ dunif(0.5*culmen_lengthu[j],1.5*culmen_lengthu[j])
barrier_infer[j] ~ dnorm(culmen_lengthu[j],1/(0.2*culmen_lengthu[j]*0.2*culmen_lengthu[j]))T(max(culmen_lengthu[j],match_infer[j]),2*culmen_lengthu[j])
#barrier_infer[j] ~ dunif(culmen_lengthu[j],2*culmen_lengthu[j])
}

# Model
for(i in 1:N){
	#new variable:
	barrier_var[i] <- ifelse(Tubelength[i] > barrier_infer[hummingbird_num[i]], 1, 0)
	mismatch_var[i] <- abs(Tubelength[i]-match_infer[hummingbird_num[i]])
	# Zero inflation / trait barrier:
	logit(pz[i]) <- Interceptpz + traitBarrier * barrier_var[i]
	Z[i] ~ dbern(min(pz[i]+0.00000000000000001,0.9999999999999999))
	
	# Interaction frequency:
	log(lambda[i]) <- Intercept + traitMismatch * mismatch_var[i] + pheno * phenoh[i]+samp*log(duration_sampling_hours[i])+
	plant_effect[plant_num[i]]+site_effect[site_num[i]]+temp_effect[site_num[i],num_time[i]]+sitebird_effect[hummingbird_species_site_num[i]]
	p[i] <- r/(r+(lambda[i] * Z[i]))
	Y[i] ~ dnegbin(min(p[i]+0.00000000000000001,0.9999999999999999),r)
}

# PRIORS FIXED EFFECTS:
Interceptpz ~ dnorm(0, 0.5)
Intercept ~ dnorm(0, 0.01) 
traitMismatch ~ dnorm(0, 0.01)T(,0)
pheno ~ dnorm(0, 0.01)
samp = 1
traitBarrier ~ dnorm(0, 0.5)T(,0)

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

writeLines(model_string,con=paste0(here("data_zenodo"),"/model.txt"))

############# TABLE S1

plant_list=list()
bird_list=list()
tabf=NULL
for(pays in c("Costa-Rica","Ecuador","Brazil")){

tab=fread(here("data_zenodo",paste0("data_for_analyses_",pays,".txt")),na.string=c("",NA))
tabf=rbind(tabf,tab)
plant_list[[pays]]=unique(tab$plant_species)
bird_list[[pays]]=unique(tab$hummingbird_species)
}
intersect(bird_list[[1]],bird_list[[2]])

site=as.data.frame(unique(tabf[,c("site","midpoint_Longitude","midpoint_Latitude","min_transect_elev","Country")]))
site=site[order(site$Country,site$min_transect_elev),]

fwrite(site,paste0(here("data_zenodo"),"/table_S1.csv"))