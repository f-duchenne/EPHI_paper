###########################################
###########################################
#' Check for packages and install them if necessary
#+ message = FALSE
pkgs <- c("randomForest","data.table", "dplyr", "lubridate","lme4","R2jags",
          "mcmcOutput","mcmcplots","MCMCvis","here") 

# Check if packages are already installed
inst <- pkgs %in% installed.packages()
# Install missing packages
if (any(!inst)) install.packages(pkgs[!inst])
# Load packages
pkg_out <- lapply(pkgs, require, character.only = TRUE)

# Set the working directory for the "here" package
here::i_am("EPHI_paper.Rproj")
# Source custom functions for Bayesian predictions
source(here("scripts/additional scripts", "function_to_predict_from_bayesian.r"))

# Initialize an empty data frame to store results
res_partition = NULL

# Loop through the three countries
for(pays in c("Costa-Rica","Ecuador","Brazil")){ #loop over the three countries

	# Load country-specific data
	load(here("data_zenodo",paste0("data_formodel_",pays,".RData"))) #loading data
	plant_res=fread(here("data_zenodo",paste0("plants_per_site_per_month_",pays,".txt")))
	plant_res=subset(plant_res,!is.na(Tubelength))
	sites=fread(here("data_zenodo",paste0("Site_metadata_",pays,".txt")),na.strings = c("",NA))

	# Load Bayesian model chains
	load(here("data_zenodo",paste0("chain_model_",pays,"_1.RData")))
	model1=results1
	load(here("data_zenodo",paste0("chain_model_",pays,"_2.RData")))
	model2=results1
	load(here("data_zenodo",paste0("chain_model_",pays,"_3.RData")))
	model3=results1
	
	# Convert models to MCMC objects
	obj1=as.mcmc(model1)
	obj2=as.mcmc(model2)
	obj3=as.mcmc(model3)

	# Combine MCMC chains and summarize
	mco <- mcmcOutput(mcmc.list(obj1,obj2,obj3))
	suma=summary(mco)
	suma$varia=rownames(suma)

	#### Extracting model-specific values
	# L1 Barrier threshold
	tabu=as.data.frame(dat[c("hummingbird_numu","culmen_lengthu")])
	tabu$barrier_infer=suma$mean[grep("barrier_infer",suma$varia)]
	# L2 Match optimum
	tabu$match_infer=suma$mean[grep("match_infer",suma$varia)]
	tabu=merge(tabu,unique(as.data.frame(dat[c("hummingbird_species","hummingbird_num")])),by.x="hummingbird_numu",by.y="hummingbird_num")

	#### ASSESS ZERO PARTITION IN THE NETWORKS
	# scale elevation
	dat$elev2=dat$elev*attr(dat$elev,"scaled:scale")+attr(dat$elev,"scaled:center")
	#merge data with site specific information
	sitebird=unique(as.data.frame(dat[c("site","hummingbird_species","month","phenoh","elev2")]))
	tab2=merge(plant_res,sitebird,by=c("site","month"),allow.cartesian=TRUE)
	# calculate phenological overlap between plants and hummingbirds
	tab2=tab2 %>% group_by(site,plant_species,hummingbird_species,Tubelength,abond_flower_moy) %>%
	summarise(pheno_predict=sum(phenoh*phenop,na.rm=T)/12)

	# Loop through different sampling durations to predict
	for (dure in c(10,100,250,500,750,1000,2000)){
		# use the homde made function to predict probability of interactions between each species pair for each site, using the Bayesian model
		pre1=predict_model(data=dat,chains=mco,nsampling=5000,site=tab2$site,bird=tab2$hummingbird_species,
		plant=tab2$plant_species,trait_plant=tab2$Tubelength,mismatch=NA,barrier=NA,
		pheno=tab2$pheno_predict,abond=rep(10,nrow(tab2)),type="network",
		random_effects=c("site","plant","bird_site"),month=NA,year=NA,nb_net=1,duration=rep(dure,nrow(tab2)))

		# Zero-inflation frequency calculations
		pre1$zero_freq=dnbinom(0,size=unique(pre1$size_param),mu=pre1$freq)
		#convert proba of interactions into probability of non-interactions
		pre1=pre1 %>% group_by(site) %>%
		mutate(pzero_infl=1-average_proba_without_barrier,
		pzero_barr=(average_proba_without_barrier-average_proba),
		pzero_freq=zero_freq*average_proba,
		pzero=pzero_freq+(1-average_proba))
		
		# Summarize sparsity estimates
		b=pre1 %>% group_by(site) %>% summarise(pzerom=mean(pzero),unknown=mean(pzero_infl),
		barrier=mean(pzero_barr),complementarity=mean(pzero_freq))

		b$tot=apply(b[,3:5],1,sum) # Calculate total probabilities
		b$duration=dure
		b$Country=pays
		# Add site elevation metadata
		b=merge(b,sites[,c("site","min_transect_elev")],by="site")
		res_partition=rbind(res_partition,b)
	}
}

# Write results to a text file
fwrite(res_partition,paste0(here("data_zenodo"),"/sparsity_estimates.txt"))

###############