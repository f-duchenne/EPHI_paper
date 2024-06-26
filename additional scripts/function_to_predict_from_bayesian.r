predict_model=function(data,chains,nsampling,site,bird,plant,trait_plant,mismatch,barrier,pheno,abond,type,random_effects,month,
year,nb_net,duration){

##### NEED TO HAVE THE DATA FROM THE MODEL LOAD

suppressWarnings(rm(list=c("barrier_var","mismatch_var")))


if(!exists("site")){site=NA}
if(!exists("bird")){bird=NA}
if(!exists("plant")){plant=NA}
if(!exists("trait_plant")){trait_plant=NA}
if(!exists("mismatch")){mismatch=NA}
if(!exists("barrier")){barrier=NA}
if(!exists("pheno")){pheno=NA}
if(!exists("abond")){abond=NA}
if(!exists("random_effects")){random_effects=NA}
if(!exists("month")){month=NA}
if(!exists("year")){year=NA}
if(!exists("type")){stop("please precise the type of prediction you want")}


####SOME FUNCTIONS:
inv.logit=function(y){exp(y)/(1+exp(y))}

rand_comb=function(y,name){
if(!is.list(y)){
resi=chains[x,paste0(name,"[",y,"]")]
return(resi)
}else{
resi=chains[x,paste0(name,"[",y[[1]],",",y[[2]],"]")]
return(resi)
}
}

###CONTROL FOR BASIC ERRORS
if(nsampling>20000){stop("nsampling is higher than 10000, it is too high for your fucking computer")}
if(nsampling<100){stop("nsampling is lower than 100, it is too low to get something trustable")}
if(!(type %in% c('barrier','frequency','total','network'))){stop("type of predict does correspond to 'barrier', 'frequency', 'total' or 'network'")}
if(!is.na(random_effects[1]) & any(!(random_effects %in% c("site","plant","bird_site","temp")))){stop("random effects to include in predictions do not correspond to 'site','plants','bird_site' or 'temp'")}


#convert discrete variables to numerical index
bidon=unique(as.data.frame(data[c("hummingbird_num","hummingbird_species")]))
if(!is.na(bird[1])){
bird_num=bidon$hummingbird_num[match(bird,bidon$hummingbird_species)]
if(length(which(is.na(bird_num)))>0){stop("Some birds were not in the original data")}
}

bidon=unique(as.data.frame(data[c("plant_num","plant_species")]))
if(!is.na(plant[1])){
plant_num=bidon$plant_num[match(plant,bidon$plant_species)]
plant_dontfit=which(is.na(plant_num))
plant_num[plant_dontfit]=1
if(length(plant_dontfit)>0){print("Some plants were not in the original data, random effect set to zero for those species")}
}

bidon=unique(as.data.frame(data[c("site_num","site")]))
if(!is.na(site[1])){
site_num=bidon$site_num[match(site,bidon$site)]
if(length(which(is.na(site_num)))>0){stop("Some sites were not in the original data")}
}

bidon=unique(as.data.frame(data[c("hummingbird_species_site_num","hummingbird_species_site")]))
if(!is.na(site[1]) & !is.na(bird[1])){bird_site_num=bidon$hummingbird_species_site_num[match(paste0(bird,site),bidon$hummingbird_species_site)]}

bidon=unique(as.data.frame(data[c("num_time","month","year")]))
if(!is.na(month[1]) & !is.na(year[1])){
temp_num=bidon$num_time[match(paste(month,year),paste(bidon$month,bidon$year))]
if(length(which(is.na(temp_num)))>0){stop("Some month-year were not in the original data")}
}



x=sample(1:nrow(chains),nsampling,replace=T)

if(type %in% c("barrier","total")){
if(is.na(trait_plant[1]) & is.na(barrier[1])){stop("please provide either values for barrier either for plant trait")}
intercept_proba=chains[x,"Interceptpz"]

if(!is.na(barrier[1])){
barrier_var=t(replicate(nsampling,barrier))
print("Use provided values for barrier instead of latent ones")
}else{
barrier_var <- ifelse(t(replicate(nsampling,trait_plant)) > rand_comb(bird_num,"barrier_infer"), 1, 0)
print("No values provided for barrier, use latent ones")
}

proba=matrix(inv.logit(intercept_proba+barrier_var*chains[x,paste0("traitBarrier")]),nrow=nsampling)
proba_suma=data.frame(site=site,hummingbird_species=bird,plant_species=plant,Tube_length=trait_plant,mismatch=mismatch,barrier=barrier,
average_proba=apply(proba,2,mean),median_proba=apply(proba,2,median),lwr_proba=apply(proba,2,quantile,prob=0.025),upr_proba=apply(proba,2,quantile,prob=0.975),duration=duration)
}


if(type %in% c("frequency","total")){
if(is.na(trait_plant[1]) & is.na(mismatch[1])){stop("please provide either values for mismatch either for plant trait")}
intercept=chains[x,"Intercept"]

if(!is.na(mismatch[1])){
mismatch_var=t(replicate(nsampling,mismatch))
print("Use provided values for mismatch instead of latent ones")
}else{
mismatch_var <-  abs(t(replicate(nsampling,trait_plant))-rand_comb(bird_num,"match_infer"))
print("No values provided for mismatch, use latent ones")
}

if(is.na(pheno[1])){phen=mean(data$phenoh,na.rm=T)*mean(chains[,"pheno"])}else{phen=chains[x,"pheno"]*t(replicate(nsampling,pheno))}
if(is.na(duration[1])){dura=log(12)*mean(chains[,"samp"])}else{dura=chains[x,"samp"]*t(replicate(nsampling,log(duration)))}
#if(is.na(abond[1])){abd=mean(data$abond_flower_log,na.rm=T)*mean(chains[,"abond"])}else{abd=chains[x,"abond"]*log(t(replicate(nsampling,abond)))}
if(is.na(plant[1]) | !("plant" %in% random_effects)){plt=0}else{
plt=rand_comb(plant_num,"plant_effect")
}
if(is.na(site[1]) | !("site" %in% random_effects)){si=0}else{si=rand_comb(site_num,"site_effect")}
if(is.na(site[1]) | is.na(month[1]) | is.na(year[1]) | !("temp" %in% random_effects)){tem=0}else{tem=rand_comb(list(site_num,temp_num),"temp_effect")}
if(is.na(site[1]) | is.na(bird[1]) | !("bird_site" %in% random_effects)){bsi=0}else{bsi=rand_comb(bird_site_num,"sitebird_effect")}


lambda_log=matrix(intercept+chains[x,"traitMismatch"]*mismatch_var+phen+dura+plt+si+tem+bsi,nrow=nsampling)
p=chains[x,"r"]/(chains[x,"r"]+exp(lambda_log))

		
if(type=="total"){Z=matrix(rbinom(ncol(proba)*nrow(proba),1,proba),nrow=nsampling)}else{
Z=1
}

frequencies=matrix(rnbinom(ncol(p)*nrow(p),size=chains[x,"r"],prob=p),nrow=nsampling)

freq_suma=data.frame(site=site,hummingbird_species=bird,plant_species=plant,Tube_length=trait_plant,pheno=pheno,abond=abond,mismatch=mismatch,barrier=barrier,month=month,year=year,
average_lambda=apply(exp(lambda_log)*Z,2,mean),lwr_lambda=apply(exp(lambda_log)*Z,2,quantile,prob=0.025),upr_lambda=apply(exp(lambda_log)*Z,2,quantile,prob=0.975),average_freq=apply(frequencies,2,mean),
median_freq=apply(frequencies,2,median),lwr_freq=apply(frequencies,2,quantile,prob=0.025),upr_freq=apply(frequencies,2,quantile,prob=0.975),duration=duration)
}


if(type=="network"){
print("Since type== network, nsampling will be ignored, all the posterior distribution is used")
###BINARY
if(is.na(trait_plant[1]) & is.na(barrier[1])){stop("please provide either values for barrier either for plant trait")}

bin_link=unique(data.frame(hummingbird_species=bird,bird_num=bird_num,plant_species=plant,trait_plant=trait_plant))
if(!is.na(barrier[1])){
barrier_var=t(replicate(1,barrier))
print("Use provided values for barrier instead of latent ones")
}else{
barrier_var <- ifelse(t(replicate(1,bin_link$trait_plant)) > apply(chains[,paste0("barrier_infer[",bin_link$bird_num,"]")],2,mean), 1, 0)
print("No values provided for barrier, use latent ones")
}

bin_link$average_proba=inv.logit(t(mean(chains[,paste0("Interceptpz")])+barrier_var*mean(chains[,paste0("traitBarrier")])))
bin_link$average_proba_without_barrier=inv.logit(mean(chains[,paste0("Interceptpz")]))

###FREQ
if(is.na(trait_plant[1]) & is.na(mismatch[1])){stop("please provide either values for mismatch either for plant trait")}

if(!is.na(mismatch[1])){
mismatch_var=t(replicate(1,mismatch))
print("Use provided values for mismatch instead of latent ones")
}else{
mismatch_var <-  abs(t(replicate(1,trait_plant))-apply(chains[,paste0("match_infer[",bird_num,"]")],2,mean))
print("No values provided for mismatch, use latent ones")
}

if(is.na(pheno[1])){phen=mean(data$phenoh,na.rm=T)*mean(chains[,"pheno"])}else{phen=mean(chains[,"pheno"])*t(replicate(1,pheno))}
if(is.na(duration[1])){dura=log(12)*mean(chains[,"samp"])}else{dura=mean(chains[,"samp"])*t(replicate(1,log(duration)))}
#if(is.na(abond[1])){abd=mean(data$abond_flower_log,na.rm=T)*mean(chains[,"abond"])}else{abd=mean(chains[,"abond"])*log(t(replicate(1,abond)))}
if(is.na(plant[1]) | !("plant" %in% random_effects)){plt=0}else{
plt=apply(chains[,paste0("plant_effect[",plant_num,"]")],2,mean)
plt[plant_dontfit]=0
}
if(is.na(site[1]) | !("site" %in% random_effects)){si=0}else{si=apply(chains[,paste0("site_effect[",site_num,"]")],2,mean)}
if(is.na(site[1]) | is.na(month[1]) | is.na(year[1]) | !("temp" %in% random_effects)){tem=0}else{tem=apply(chains[,paste0("temp_effect[",site_num,",",temp_num,"]")],2,mean)}
if(is.na(site[1]) | is.na(bird[1]) | !("bird_site" %in% random_effects)){bsi=0}else{bsi=apply(chains[,paste0("sitebird_effect[",bird_site_num,"]")],2,mean)}

lambda_log=t(mean(chains[,paste0("Intercept")])+mean(chains[,"traitMismatch"])*mismatch_var+phen+dura+plt+si+tem+bsi)


networks=NULL
for(j in 1:nb_net){
bin_link$binary=rbinom(nrow(bin_link),1,bin_link$average_proba)
network=data.frame(essai=j,site=site,hummingbird_species=bird,plant_species=plant,Tube_length=trait_plant,pheno=pheno,
abond=abond,mismatch=mismatch,barrier=barrier,
month=month,year=year,duration=duration)
network=merge(network,bin_link,by=c("plant_species","hummingbird_species"))
network$freq=exp(lambda_log)
network$size_param=mean(chains[,"r"])
networks=rbind(networks,network)
}
}


if(type=="barrier"){
return(proba_suma)
}else{
if(type=="network"){
return(networks)
}else{
return(freq_suma)
}}

}


