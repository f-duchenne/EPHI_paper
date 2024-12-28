###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
pkgs <- c("randomForest","data.table", "dplyr", "lubridate","lme4","R2jags","mcmcOutput","mcmcplots","MCMCvis","pastecs","ggplot2","cowplot","gridExtra","scales","reshape2","bipartite","stringr") 



# Install missing packages
inst <- pkgs %in% installed.packages()
if (any(!inst)) install.packages(pkgs[!inst])

# Load packages
pkg_out <- lapply(pkgs, require, character.only = TRUE)

# Set project path
here::i_am("EPHI_paper.Rproj")

# Source additional functions
source(here("scripts/additional scripts", "function_to_predict_from_bayesian.r"))

# Loop through the three countries
for(pays in c("Costa-Rica","Ecuador","Brazil")){

	# Load data
	load(here("data_zenodo",paste0("data_formodel_",pays,".RData")))
	plant_res=fread(here("data_zenodo",paste0("plants_per_site_per_month_",pays,".txt")))
	plant_res=subset(plant_res,!is.na(Tubelength))

	#### Load MCMC chains
	load(here("data_zenodo",paste0("chain_model_",pays,"_1.RData")))
	model1=results1
	load(here("data_zenodo",paste0("chain_model_",pays,"_2.RData")))
	model2=results1
	load(here("data_zenodo",paste0("chain_model_",pays,"_3.RData")))
	model3=results1
	# Convert models to MCMC objects and summarize
	obj1=as.mcmc(model1)
	obj2=as.mcmc(model2)
	obj3=as.mcmc(model3)
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
	 # Rescale elevation
	dat$elev2=dat$elev*attr(dat$elev,"scaled:scale")+attr(dat$elev,"scaled:center")
	# Calculate phenological overlap
	sitebird=unique(as.data.frame(dat[c("site","hummingbird_species","month","phenoh","elev2")]))
	tab2=merge(plant_res,sitebird,by=c("site","month"),allow.cartesian=TRUE)
	tab2=tab2 %>% group_by(site,plant_species,hummingbird_species,Tubelength,abond_flower_moy) %>%
	summarise(pheno_predict=sum(phenoh*phenop,na.rm=T)/12)

	# Predict network
	pre1=predict_model(data=dat,chains=mco,nsampling=5000,site=tab2$site,bird=tab2$hummingbird_species,
	plant=tab2$plant_species,trait_plant=tab2$Tubelength,mismatch=NA,barrier=NA,
	pheno=tab2$pheno_predict,abond=rep(10,nrow(tab2)),type="network",
	random_effects=c("site","plant"),month=NA,year=NA,nb_net=1,duration=rep(12,nrow(tab2)))

	#merge with hummingbird identity (latin species names)
	pre1=merge(pre1,tabu,by="hummingbird_species")

	sites=fread(here("data_zenodo",paste0("Site_metadata_",pays,".txt")),na.strings = c("",NA))
	sites=subset(sites,habitat!="deforested")
	# Combine predictions with site metadata
	pre1b=merge(pre1,sites,by="site")
	pre1b$Country=pays
	 # Save results
	fwrite(pre1b,paste0(here("data_zenodo"),"/initial_network_per_day",pays,".txt"))
}

