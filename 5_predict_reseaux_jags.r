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

#### free prior
load("chain_model_free_1.RData")
model1=results1
load("chain_model_free_2.RData")
model2=results1
load("chain_model_free_3.RData")
model3=results1
obj1=as.mcmc(model1)
obj2=as.mcmc(model2)
obj3=as.mcmc(model3)
#combine chains and get a summary:
mco_free <- mcmcOutput(mcmc.list(obj1,obj2,obj3))
suma_free=summary(mco_free)
suma_free$varia=rownames(suma_free)

#### L1 Barrier threshold
tabu=as.data.frame(dat[c("hummingbird_numu","bill_lengthu")])
tabu=tabu[order(tabu$hummingbird_numu),]
tabu$barrier_infer=suma$mean[grep("barrier_infer",suma$varia)]
#### L2 Match optimum
tabu$match_infer=suma$mean[grep("match_infer",suma$varia)]
tabu=merge(tabu,unique(as.data.frame(dat[c("hummingbird_species","hummingbird_num")])),by.x="hummingbird_numu",by.y="hummingbird_num")


#####FIT TO THE DATA

tab=as.data.frame(dat[c("plant_species","hummingbird_species","site","month","year","waypoint","Tube_length","value","abond_flower_log","Yfreq_r")])
tab=tab[sample(1:nrow(tab),500,replace=F),]

pre1=predict_model(chains=mco,nsampling=5000,site=tab$site,bird=tab$hummingbird_species,plant=tab$plant_species,trait_plant=tab$Tube_length,mismatch=NA,barrier=NA,
pheno=tab$value,abond=tab$abond,type="total",
random_effects=c("site","plants","bird_site","temp"),month=tab$month,year=tab$year,nb_net=1)
pre1$obs=tab$Yfreq_r

pre1=merge(pre1,tabu,by=c("hummingbird_species"))
cor(pre1$average_lambda,pre1$obs)

plot(average_lambda~obs,data=pre1)
subset(pre1,obs==0 & average_lambda>40)

##### NETWORK FOR ROBUSTESSE
emp_net=fread("empirical_networks.txt")
plant_res=fread("plants_per_site_per_month.txt")
plant_res=subset(plant_res,!is.na(Tube_length))
sitebird=unique(as.data.frame(dat[c("site","hummingbird_species","month","phenoh")]))
tab2=merge(plant_res,sitebird,by=c("site","month"),allow.cartesian=TRUE)
tab2=tab2 %>% group_by(site,plant_species,hummingbird_species,Tube_length) %>% summarise(pheno_predict=sum(phenoh)/12)

tab2=merge(tab2,tabu,by="hummingbird_species")

tab2$comp=-abs(tab2$Tube_length-tab2$match_infer)

bidon=subset(tab2,site=="QUEB")
bidon=bidon %>% group_by(hummingbird_species,site) %>% mutate(inter_tot=sum(inter))

plant_list=unique(bidon$plant_species)
j=1
ext=plant_list[j]


i="Thalurania colombica"
bidon2=subset(bidon,hummingbird_species==i)
dist_mat=as.matrix(dist(bidon2$comp,upper=T)/max(dist(bidon2$comp)))
diag(dist_mat)=NA






exp(apply(log(dist_mat),1,mean,na.rm=T))

bidon$comp

pre1=predict_model(chains=mco,nsampling=5000,site=tab2$site,bird=tab2$hummingbird_species,plant=tab2$plant_species,trait_plant=tab2$Tube_length,mismatch=NA,barrier=NA,
pheno=NA,abond=NA,type="network",
random_effects=c("site","plants","temp"),month=tab2$month,year=tab2$year,nb_net=100)

networks=pre1 %>% group_by(site,essai,plant_species,hummingbird_species) %>% summarise(binary=mean(binary),freq=sum(freq),inter=sum(inter))

structure_tab=networks %>% group_by(site,essai) %>% summarise(C_bin=length(binary[binary>0])/(length(unique(plant_species))*length(unique(hummingbird_species))),
C_freq=length(freq[freq>0])/(length(unique(plant_species))*length(unique(hummingbird_species))),
C_inter=length(inter[inter>0])/(length(unique(plant_species))*length(unique(hummingbird_species))))

for(i in 1:nrow(structure_tab)){
bidon=subset(networks,site==structure_tab$site[i] & essai==structure_tab$essai[i])

mat=as.data.frame(dcast(bidon,plant_species~hummingbird_species,value.var="binary"))
row.names(mat)=mat[,1]
mat=data.matrix(mat[,-1])
mat[mat>0]=1
structure_tab$NODF_bin[i]=networklevel(mat,index="NODF")[[1]]
structure_tab$Mod_bin[i]=tryCatch({computeModules(as.matrix(mat), method="Beckett")@likelihood},error=function(x){NA})

mat=as.data.frame(dcast(bidon,plant_species~hummingbird_species,value.var="freq"))
row.names(mat)=mat[,1]
mat=data.matrix(mat[,-1])
mat[mat>0]=1
structure_tab$NODF_freq[i]=networklevel(mat,index="NODF")[[1]]
structure_tab$Mod_freq[i]=tryCatch({computeModules(as.matrix(mat), method="Beckett")@likelihood},error=function(x){NA})

mat=as.data.frame(dcast(bidon,plant_species~hummingbird_species,value.var="inter"))
row.names(mat)=mat[,1]
mat=data.matrix(mat[,-1])
mat[mat>0]=1
structure_tab$NODF_inter[i]=networklevel(mat,index="NODF")[[1]]
structure_tab$Mod_inter[i]=tryCatch({computeModules(as.matrix(mat), method="Beckett")@likelihood},error=function(x){NA})
}


structure_tab_obs=emp_net %>% group_by(site) %>% summarise(C_obs=length(Y[Y>0])/(length(unique(plant_species))*length(unique(hummingbird_species))))

for(i in 1:nrow(structure_tab_obs)){
bidon=subset(emp_net,site==structure_tab_obs$site[i])
mat=as.data.frame(dcast(bidon,plant_species~hummingbird_species,value.var="Y"))
row.names(mat)=mat[,1]
mat=data.matrix(mat[,-1])
mat[mat>0]=1
structure_tab_obs$NODF_obs[i]=networklevel(mat,index="NODF")[[1]]
structure_tab_obs$Mod_obs[i]=computeModules(mat, method="Beckett")@likelihood
}


structure_tabf=merge(structure_tab,structure_tab_obs,by=c("site"))

b=melt(structure_tab,id.vars=c("site","essai"))
# Split name column into firstname and last name
b[c('varia', 'type')] <- str_split_fixed(b$variable, '_', 2)

b_obs=melt(structure_tab_obs,id.vars=c("site"))
# Split name column into firstname and last name
b_obs[c('varia', 'type')] <- str_split_fixed(b_obs$variable, '_', 2)

ggplot()+geom_boxplot(data=b,aes(x=site,y=value,color=type))+
geom_hpline(data=b_obs,aes(x=site,y=value),stat = "identity",size=1)+
facet_wrap(~varia,scales="free")+coord_flip()


par(mfrow=c(1,2))
examp=subset(networks,site=="TOLO" & essai==1)
mat=dcast(examp,plant_species~hummingbird_species,value.var="inter")
row.names(mat)=mat[,1]
mat=mat[,-1]
plotweb(mat)

examp=subset(emp_net,site=="TOLO")
mat=as.data.frame(dcast(examp,plant_species~hummingbird_species,value.var="Y"))
row.names(mat)=mat[,1]
mat=data.matrix(mat[,-1])
plotweb(mat)

tab2=unique(as.data.frame(dat[c("site","month","year","num_time")]))
pre1=predict_model(chains=mco,nsampling=500,site=tab2$site,bird=NA,plant=NA,trait_plant=NA,mismatch=rep(0,nrow(tab2)),barrier=rep(0,nrow(tab2)),
pheno=NA,abond=rep(1,nrow(tab2)),type="total",
random_effects=c("site","temp"),month=tab2$month,year=tab2$year,nb_net=NA)
pre1$num_time=tab2$num_time

ggplot(data=pre1,aes(x=num_time,y=average_lambda,color=site))+geom_line()