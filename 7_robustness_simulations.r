library(dplyr)
library(reshape2)
library(bipartite)
library(data.table)

setwd(dir="path_to_data_on_cluster") #script running on HPC, path should be adjusted for the new cluster path

# Collect command arguments
args <- commandArgs(trailingOnly = TRUE)
args_contents <- strsplit(args, ' ')
# Get first argument
jj <- as.numeric(args_contents[[1]])
print(jj)

tab=expand.grid(c("Brazil","Costa-Rica","Ecuador"),1:100)
pays=tab$Var1[jj]
e=tab$Var2[jj]

pre1=fread(paste0("initial_network_",pays,".txt"))
sites=fread(paste0("Site_metadata_",pays,".txt"),na.strings = c("",NA))
sites=subset(sites,habitat!="deforested")

r_vec=0:2
scenarios=c("generalists first","specialists first")
rules=c("with forbidden links","no forbidden links")
sites_vec=unique(pre1$site)

extinctions=NULL
pre1b=pre1
for(r in r_vec){
	for(scenario in scenarios){
		for(bar in rules){
			for(z in sites_vec){

				bidon=subset(pre1b,site==z)

				if(bar=="no forbidden links"){
				bidon$proba=bidon$average_proba_without_barrier
				}else{
				bidon$proba=bidon$average_proba
				}

				bidon$inter=bidon$freq*bidon$proba
				bidon$compl=-abs(bidon$Tube_length-bidon$match_infer)

				bidon=bidon %>% dplyr::group_by(plant_species) %>% mutate(degree=length(unique(hummingbird_species[inter>0])))

				if(bar=="generalists first"){bidon=bidon[order(bidon$degree,bidon$plant_species,decreasing=TRUE),]}else{
				bidon=bidon[order(bidon$degree,bidon$plant_species,decreasing=FALSE),]}

				prop_forbidden=mean(1-bidon$average_proba_without_barrier)-mean(1-bidon$proba)
				average_freq=mean(bidon$freq)
				compl=mean(bidon$compl)

				mat=dcast(bidon,plant_species~hummingbird_species,value.var="inter")
				mat=as.matrix(mat[,-1])
				C=networklevel(mat,index="weighted connectance",weighted=T)
				N=wine(mat, nreps = 1000)$wine
				M=computeModules(mat, method="Beckett")@likelihood[1]
				Cperso=mean(mat)
				Inter_eve=networklevel(mat,index="interaction evenness",weighted=T)

				plant_list=unique(bidon$plant_species)
				extinctions=rbind(extinctions,data.frame(plant_ext=NA,hummingbird_pers=length(unique(bidon$hummingbird_species)),site=z,essai=e,r=r,
				rank=0,barrier=bar,prop_forbidden=prop_forbidden,compl=compl,average_freq=average_freq,scenario=scenario,C=C,N=N,M=M,Cperso=Cperso,Inter_eve=Inter_eve))
				for(j in 1:length(plant_list)){
					ext=plant_list[j]
					if(length(unique(bidon$plant_species))>1){
						for(i in unique(bidon$hummingbird_species)){
							bidon2=subset(bidon,hummingbird_species==i)
							distances=sqrt((bidon2$compl[bidon2$plant_species==ext]-bidon2$compl[bidon2$plant_species!=ext])^2)
							R=mean(bidon2$proba[bidon2$plant_species!=ext]*exp(-r*distances/max(distances)))

							Pext=(1-R)*bidon2$inter[bidon2$plant_species==ext]/sum(bidon2$inter)

							EXT=rbinom(1,1,Pext)

							if(EXT==1){bidon=subset(bidon,hummingbird_species!=i)}
						}
					}
					bidon=subset(bidon,plant_species!=ext)
					extinctions=rbind(extinctions,data.frame(plant_ext=ext,hummingbird_pers=length(unique(bidon$hummingbird_species)),site=z,essai=e,r=r,
					rank=j,barrier=bar,prop_forbidden=prop_forbidden,compl=compl,average_freq=average_freq,scenario=scenario,C=C,N=N,M=M,Cperso=Cperso,Inter_eve))
				}
			}
		}
	}
}


#LOAD SITE METADATA:
extinctions=merge(extinctions,sites,by="site")

setwd(dir="path_to_data_on_cluster/robust_data_res/")
fwrite(extinctions,paste0("robustness_simulations_",pays,"_",e,"_.txt"))

#######################################################################################################################
###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
pkgs <- c("data.table", "dplyr","here") 

# Check if packages are already installed
inst <- pkgs %in% installed.packages()
# Install missing packages
if (any(!inst)) install.packages(pkgs[!inst])
# Load packages
pkg_out <- lapply(pkgs, require, character.only = TRUE)

# Set the working directory for the "here" package
here::i_am("EPHI_paper.Rproj")

#running this part requires to have run the first part of Script 7 on a HPC and to have transfer results to the computer

### PUT all files together
# initilaize an empty data frame
extinctions=NULL
#loop over 3 countries and over the 500 replicates
for(pays in c("Costa-Rica","Ecuador","Brazil")){
	for(e in 1:500){
		#load simulations
		extinctionsc=fread(here("data_zenodo",paste0("robustness_simulations_per_day_",pays,"_",e,"_.txt")))
		extinctionsc$Country=pays
		# combine simulations together
		extinctions=rbind(extinctions,extinctionsc)
	}
}

#some useful variables
extinctions=extinctions %>% group_by(r,site,Country,barrier,essai,scenario) %>% mutate(nbh=max(hummingbird_pers),nbp=max(rank))
extinctions$persistence=extinctions$hummingbird_pers/extinctions$nbh
extinctions$rank2=extinctions$rank/extinctions$nbp
extinctions=extinctions[order(extinctions$min_transect_elev),]
extinctions$site=factor(extinctions$site,levels=unique(extinctions$site))
extinctions$site2=paste0(extinctions$site," (",extinctions$min_transect_elev,"m)")
extinctions$site2=factor(extinctions$site2,levels=unique(extinctions$site2))
extinctions= extinctions %>% group_by(site,Country,site2,barrier) %>% mutate(prop_forbidden_m=-1*mean(prop_forbidden),C_m=mean(C),N_m=mean(N),M_m=mean(M),compl_m=mean(compl),Cperso_m=mean(Cperso),
										Inter_eve_m=mean(Inter_eve))

extinctions$Country=gsub("-"," ",extinctions$Country,fixed=T)

#save all simulations results in one single file
fwrite(extinctions,paste0(here("data_zenodo"),"/robustness_simulations_all.csv")
###################