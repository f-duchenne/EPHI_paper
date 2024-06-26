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

print(i)

set.seed=i+5
nbchains=3
lili=expand.grid(1:nbchains,c("Costa-Rica","Ecuador","Brazil"))


load(paste0("/home/duchenne/EPHI_paper/data_formodel_",lili$Var2[i],".RData"))

ParsStage <- c("barrier_infer","match_infer","Interceptpz","traitBarrier","Intercept","traitMismatch","pheno","abond","plant_effect","site_effect","temp_effect","sitebird_effect",
"sd.plant","sd.bird","sd.site","sd.temp","r","edec","samp")

Inits <- function(){list()}
#load(paste0("/home/duchenne/EPHI_paper/initial_conditions_",lili$Var2[i],".RData"))
#Inits =lisf[lili$Var1[i]]


if(lili$Var2[i]=="Ecuador"){
  burnin=150000
  iter=250000
}else{
  burnin=650000
  iter=750000
}

results1 <- jags(model.file="model.txt", parameters.to.save=ParsStage,  n.chains=1, data=dat,n.burnin = burnin,  n.iter = iter, n.thin = 5,inits =Inits,jags.seed =i)

save(results1,file=paste0("chain_model_",lili$Var2[i],"_",lili$Var1[i],".RData"))
#