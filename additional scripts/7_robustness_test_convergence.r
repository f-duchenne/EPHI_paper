#######################################################################################################################
###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("randomForest","data.table", "dplyr", "lubridate","lme4","R2jags","mcmcOutput","mcmcplots","MCMCvis","qgraph","igraph",
			"pastecs","ggplot2","cowplot","gridExtra","scales","reshape2","bipartite","stringr","ungeviz","lme4","mgcv","ggpubr","emmeans","piecewiseSEM","car") 


inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)

EPHI_version="2024-03-18"

couleurs=c("#679436","#0B4F6C","deeppink")

setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper/data")

extinctions=fread("robustness_simulations_all.csv")

pval=NULL
extinctions=subset(extinctions,r==1 & scenario=="generalists first")

for(j in c(10,50,150,200,300)){
for(i in 1:10){
b=extinctions[extinctions$essai %in% sample(1:500,j),] %>% group_by(r,site,Country,site2,barrier,rank2,min_transect_elev,nbh,nbp,scenario) %>% summarise(persistence=mean(hummingbird_pers/nbh),
pers_sde=sd(hummingbird_pers)/sqrt(length(hummingbird_pers)))
b$symmetrie=b$nbp/b$nbh

pays="Brazil"

###AREA
b2=subset(b,rank2!=0 & rank2<1) %>% group_by(site,site2,Country,barrier,r,min_transect_elev,nbh,nbp,scenario) %>%
summarise(robustness=sum(persistence)/length(persistence))
b2=b2 %>% group_by(site,scenario,r) %>% mutate(pers_eff=(robustness[barrier=="with forbidden links"]-robustness[barrier=="no forbidden links"])/robustness[barrier=="no forbidden links"])


obj=wilcox.test(subset(b2,scenario=="generalists first" & barrier=="no forbidden links" & r==1 & Country=="Costa Rica")$robustness,
subset(b2,scenario=="generalists first" & barrier=="with forbidden links" & r==1 & Country=="Costa Rica")$robustness,paired=T)
pval=rbind(pval,data.frame(pval=obj$p.value,j=j))
}}

boxplot(pval~j,data=pval)

