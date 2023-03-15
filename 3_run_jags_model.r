###########################################
library(data.table)
library(dplyr)
library(R2jags)
setwd(dir="/home/duchenne/EPHI_paper")

# Collect command arguments
args <- commandArgs(trailingOnly = TRUE)
args_contents <- strsplit(args, ' ')
# Get first argument
i <- as.numeric(args_contents[[1]])

set.seed=i+5

load("/home/duchenne/EPHI_paper/data_formodel.RData")

ParsStage <- c("barrier_infer","match_infer","Interceptpz","traitBarrier","Intercept","traitMismatch","pheno","abond","plant_effect","site_effect","temp_effect","sitebird_effect","sd.plant","sd.bird","sd.site","sd.temp","r","edec")

Inits <- function(){list()}

results1 <- jags(model.file="model.txt", parameters.to.save=ParsStage, n.chains=1, data=dat,n.burnin = 25000,  n.iter = 30000, n.thin = 3,inits =Inits,jags.seed =i)

save(results1,file=paste0("chain_model_",i,".RData"))
#
