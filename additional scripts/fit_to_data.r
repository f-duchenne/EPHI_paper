source("C:/Users/Duchenne/Documents/EPHI_paper/scripts/function_to_predict_from_bayesian.r")
setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper/data")
library(mcmcOutput)
library(mcmcplots)
library(MCMCvis)
library(ggridges)
library(runjags)
library(pastecs)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(dplyr)
library(reshape2)
library(bipartite)
library(stringr)
library(ungeviz)
library(data.table)
library(doParallel)
library(foreach)
library(parallel)

load("data_formodel.RData")

#### conservative prior
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
suma=summary(mco)
suma$varia=rownames(suma)

load("chain_model_free_1.RData")
model1=results1
load("chain_model_free_2.RData")
model2=results1
load("chain_model_free_3.RData")
model3=results1
obj1=as.mcmc(model1)
obj2=as.mcmc(model2)
obj3=as.mcmc(model3)
plot(mcmc.list(obj1,obj2,obj3)[,"pheno"])
#combine chains and get a summary:
mco_free <- mcmcOutput(mcmc.list(obj1,obj2,obj3))
suma_free=summary(mco_free)
suma_free$varia=rownames(suma_free)

pre1=predict_model(chains=mco_free,nsampling=5000,site=dat$site,bird=dat$hummingbird_species,plant=dat$plant_species,trait_plant=dat$Tube_length,mismatch=NA,barrier=NA,
pheno=dat$phenoh,abond=exp(dat$abond_flower_log),type="total",
random_effects=c("site","plant","bird_site","temp"),month=dat$month,year=dat$year,nb_net=NA)

pre1$Y=dat$Yfreq_r

cor(pre1[,c("Y","average_lambda","average_freq")])
1-var(pre1$Y-pre1$average_freq)/var(pre1$Y)

plot(Y~average_freq,data=pre1)
abline(0,1)
