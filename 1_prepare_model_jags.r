###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
pkgs <- c("randomForest","data.table", "dplyr", "lubridate","lme4","R2jags","here") 

# Check if packages are already installed
inst <- pkgs %in% installed.packages()
# Install missing packages
if (any(!inst)) install.packages(pkgs[!inst])
# Load packages
pkg_out <- lapply(pkgs, require, character.only = TRUE)

# Set the working directory for the "here" package
here::i_am("EPHI_paper.Rproj")

# Loop over countries to process data for each
for(pays in c("Costa-Rica","Ecuador","Brazil")){

	# Load pre-processed data for the current country
	tab=fread(here("data_zenodo",paste0("data_for_analyses_",pays,".txt")),na.string=c("",NA))

	#### Calculate Trait Matching:
    # Compute mismatch between plant corolla length (Tubelength) and hummingbird bill length (culmen_length)
	tab$mismatch=abs(tab$Tubelength-tab$culmen_length) #trait mismatch
	
	# Scale elevation values for numerical stability
	tab$elev=scale(tab$min_transect_elev)

	 #### Calculate Phenological Matching:
     # Multiply hummingbird phenology (`phenoh`) by plant phenology (`phenop`) to compute phenological overlap
	tab$phenomatch=tab$phenoh*tab$phenop

	#### Handle Missing Flower Abundance Data:
	# For missing or zero flower abundance values, set to a small default value (1)
	tab$abond_flower[is.na(tab$abond_flower) | tab$abond_flower==0]=1
	# Compute the logarithm of flower abundance for the model
	tab$abond_flower_log=log(tab$abond_flower)

	#### Data Cleaning:
	# Remove rows where mismatch data is missing
	tab=subset(tab,!is.na(mismatch))
	unique(tab$site)

	#### Prepare Data for Temporal and Spatial Autocorrelation:
	# Convert site to a factor
	tab$site <- factor(tab$site)
	# Create a placeholder group variable
	tab$group <- factor(rep(1, nrow(tab)))
	# Combine hummingbird species and site into a single variable
	tab$hummingbird_species_site <- paste0(tab$hummingbird_species, tab$site)

	#### Convert Factors to Numeric Indices for JAGS:
	# Numeric indices for hummingbird species by site
	tab$hummingbird_species_site_num <- as.numeric(as.factor(tab$hummingbird_species_site))
	# Numeric indices for hummingbird species
	tab$hummingbird_num <- as.numeric(as.factor(tab$hummingbird_species))
	# Numeric indices for plant species
	tab$plant_num <- as.numeric(as.factor(tab$plant_species))
	# Numeric indices for sites
	tab$site_num <- as.numeric(as.factor(tab$site))
	# Numeric indices for temporal points (starting from 1)
	tab$num_time <- as.numeric(as.factor(tab$num_time))
	
	#### Extract Unique Hummingbird Traits:
	# Create a table with unique bill lengths for hummingbird species
	tabu <- unique(tab[, c("hummingbird_num", "culmen_length")])
	names(tabu) <- paste0(names(tabu), "u")  # Rename columns to avoid confusion with original columns
	tabu <- tabu[order(tabu$hummingbird_numu), ]  # Sort by hummingbird_numu


	#create a table with corolla length for plants
	# tabp=unique(tab[,c("plant_num","Tubelength")])
	# names(tabp)=paste0(names(tabp),"p")
	# tabp=tabp[order(tabu$plant_nump),]

	#### Calculate Geographic Distances Between Sites:
	# Extract unique site information with geographic coordinates
	sites=unique(tab[,c("site","site_num","midpoint_Longitude","midpoint_Latitude")])
	sites=sites[order(sites$site_num),]
	# Compute a distance matrix based on geographic coordinates
	Distance=dist(sites[,c("midpoint_Longitude","midpoint_Latitude")], method="euclidean", diag=TRUE, upper=TRUE)
	 # Normalize distances to range [0, 1]
	Distance=as.matrix(Distance)/max(Distance)

	#### Assemble Data for JAGS Model:
	dat=c(as.list(tab),as.list(tabu),list(N=nrow(tab),
	Nbirds=length(unique(tab$hummingbird_num)),Nbird_site=length(unique(tab$hummingbird_species_site_num)),Nsites=length(unique(tab$site_num)),Nplants=length(unique(tab$plant_num)),
	Ntemp=length(unique(tab$num_time)),Nsite=length(unique(tab$site_num)),Distance=Distance))

	# Save the processed data to a file for modeling
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