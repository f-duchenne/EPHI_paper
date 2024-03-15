library(dplyr)
library(R2jags)
library(data.table)
library(mcmcplots)
library(MCMCvis)
library(mcmcOutput)
library(ggridges)
library(runjags)
library(coda)


filepath="C:/Users/Duchenne/Documents/EPHI_paper/data"
setwd(dir=filepath)


set.seed=i+5
nbchains=3
for(i in 1:2){
pays=c("Costa-Rica","Ecuador","Brazil")[i]

load(paste0("chain_model_",pays,"_",1,".RData"))
model1=results1
load(paste0("chain_model_",pays,"_",2,".RData"))
model2=results1
load(paste0("chain_model_",pays,"_",3,".RData"))
model3=results1

obj1=as.mcmc(model1)
obj2=as.mcmc(model2)
obj3=as.mcmc(model3)


mco1 <- mcmcOutput(mcmc.list(obj1,obj2,obj3))
mco=summary(mco1,digits=15)
mco$varia=rownames(mco)
mco$varia2=sapply(strsplit(mco$varia,"[",fixed=T),function(x){x[1]})
mco$index=gsub("]","",sapply(strsplit(mco$varia,"[",fixed=T),function(x){x[2]}),fixed=T)

lisf=list()

for(j in 1:nbchains){
lis=list()
for(i in unique(mco$varia2)){
bidon=subset(mco,varia2==i)
if(length(grep(",",bidon$index,fixed=T))==0){
lis[[paste(i)]]=bidon$mean
if(i=="temp_effect"){lis[[paste(i)]][1]=NA}
}else{
if(i=="lvcoef"){
mat=matrix(bidon$mean,nrow=max(as.numeric(sapply(strsplit(bidon$index,",",fixed=T),function(x){x[1]}))),ncol=max(as.numeric(sapply(strsplit(bidon$index,",",fixed=T),function(x){x[2]}))))
}else{
mat=matrix(bidon$mean,nrow=max(as.numeric(sapply(strsplit(bidon$index,",",fixed=T),function(x){x[1]}))),
ncol=max(as.numeric(sapply(strsplit(bidon$index,",",fixed=T),function(x){x[2]}))))
}
lis[[paste(i)]]=mat
}
}
lis=lis[-grep("deviance",names(lis))]
lisf[[j]]=lis
}


save(lisf,file=paste0("initial_conditions_",pays,".RData"))

}
