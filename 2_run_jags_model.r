###########################################
library(data.table)
library(dplyr)
library(R2jags)
setwd(dir="path_to_data_on_cluster") #script running on HPC, path should be adjusted for the new cluster path

# Collect command arguments
args <- commandArgs(trailingOnly = TRUE)
args_contents <- strsplit(args, ' ')
# Get first argument
i <- as.numeric(args_contents[[1]]) #iterations happens outside of R, in the slurm scripts used to run the R script on the cluster

set.seed=i+5
nbchains=3
lili=expand.grid(1:nbchains,c("Costa-Rica","Ecuador","Brazil"))

load(paste0("data_formodel_",lili$Var2[i],".RData"))

ParsStage <- c("barrier_infer","match_infer","Interceptpz","traitBarrier","Intercept","traitMismatch","pheno","abond","plant_effect","site_effect","temp_effect","sitebird_effect","sd.plant","sd.bird","sd.site","sd.temp","r","edec","samp")

Inits <- function(){list()}

results1 <- jags(model.file="model.txt", parameters.to.save=ParsStage,  n.chains=1, data=dat,n.burnin = 180000,  n.iter = 200000, n.thin = 6,inits =Inits,jags.seed =i)

save(results1,file=paste0("chain_model_",lili$Var2[i],"_",lili$Var1[i],".RData"))
#