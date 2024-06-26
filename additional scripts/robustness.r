library(dplyr)
library(reshape2)
library(bipartite)
library(data.table)

EPHI_version="2023-07-04"

setwd(dir="/home/duchenne/EPHI_paper/robust_data_res/")

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
Cperso=length(mat[mat>=0.05*sum(mat)])/length(mat)
vec=sort(c(mat),decreasing=T)
vec1=vec[1:length(cumsum(vec)[cumsum(vec)<0.95*sum(vec)])]
Cperso2=length(vec1)/length(mat)

plant_list=unique(bidon$plant_species)
extinctions=rbind(extinctions,data.frame(plant_ext=NA,hummingbird_pers=length(unique(bidon$hummingbird_species)),site=z,essai=e,r=r,
rank=0,barrier=bar,prop_forbidden=prop_forbidden,compl=compl,average_freq=average_freq,scenario=scenario,C=C,N=N,M=M,Cperso=Cperso,Cperso2=Cperso2))
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
rank=j,barrier=bar,prop_forbidden=prop_forbidden,compl=compl,average_freq=average_freq,scenario=scenario,C=C,N=N,M=M,Cperso=Cperso,Cperso2=Cperso2))


}
}
}
}
}


#LOAD SITE METADATA:
extinctions=merge(extinctions,sites,by="site")

fwrite(extinctions,paste0("robustness_simulations_",pays,"_",e,"_.txt"))
