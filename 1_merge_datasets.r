###########################################
library(bipartite)
library(Rmpfr)
library(circular)
library(CircStats)
library(plot3D)
library(ggradar)
library(ggplot2)
library(gridExtra)
library(data.table)
library(dplyr)
library(rgbif)
library(ggforce)
library(ggeffects)
library(ggExtra)
library(viridis)
library(lme4)
library(cowplot)
library(scales)
library(car)
library(DHARMa)
library(glmmTMB)
library(qgraph)
library(igraph)
library(piecewiseSEM)
library(randomForest)
library(FactoMineR)
require(factoextra)
library("ggdendro")
library(dendextend)
library(ape)
library(ggtree)
library(ggnewscale)
library(geiger)
library(diversityForest)
library(caper)
library(phytools)
library(ggthemes)
library(ggplotify)
library(ggtern)
library(ks)
library(sp)
library(spatialEco)
library(MuMIn)
library(INLA)
library(inlabru)
library("inlaVP")
setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper/data")

EPHI_version="2023-02-14"


#LOAD HUMMINGBIRD DATA FROM COSTA-RICA:
dat_cr=fread(paste0("C:/Users/Duchenne/Documents/EPHI_data_clean/Costa-Rica_",EPHI_version,"/Interactions_data_Costa-Rica.txt"),na.strings = c("",NA))
dat_cr$date=as.IDate(dat_cr$date,"%Y/%m/%d") #be sure date columns is recognize as date
#LOAD CAMERA INFORMATION:
cameras=fread(paste0("C:/Users/Duchenne/Documents/EPHI_data_clean/Costa-Rica_",EPHI_version,"/Cameras_data_Costa-Rica.txt"),na.strings = c("",NA))
cameras$end_date=as.IDate(cameras$end_date,"%Y/%m/%d") #be sure date columns is recognize as date
cameras$start_date=as.IDate(cameras$start_date,"%Y/%m/%d") #be sure date columns is recognize as date
cameras$month=month(cameras$start_date) #extract month from date column
cameras$year=year(cameras$start_date) #extract year from date column
plant_for_res=unique(cameras[,c("plant_species","month","year","site")])
#LOAD TRANSECT DATA:
tr_cr=fread(paste0("C:/Users/Duchenne/Documents/EPHI_data_clean/Costa-Rica_",EPHI_version,"/Transect_data_Costa-Rica.txt"),na.strings = c("",NA))
tr_cr$date=as.IDate(tr_cr$date,"%Y/%m/%d") #be sure date columns is recognize as date
tr_cr$month=month(tr_cr$date) #extract month from date column
tr_cr$year=year(tr_cr$date) #extract year from date column
tr_cr=subset(tr_cr,!is.na(plant_species)) %>% group_by(site,year,month,plant_species) %>% summarise(abond_flower=sum(total_flowers,na.rm=T)) #calculate a total abundance per plant species per month per site
#LOAD SITE METADATA:
sites_cr=fread(paste0("C:/Users/Duchenne/Documents/EPHI_data_clean/Costa-Rica_",EPHI_version,"/Site_metadata_Costa-Rica.txt"),na.strings = c("",NA))


#MERGE THEM:
#interactions and cameras:
dim(dat_cr)
dat_cr=merge(dat_cr,cameras,by=c("waypoint","site"),all.x=T,all.y=T) #we want to keep cameras that did not detect any hummingbird
dim(dat_cr)
dat_cr=subset(dat_cr,!is.na(plant_species)) #remove data with no plant ID
dim(dat_cr)
dat_cr=subset(dat_cr,duration_sampling_hours>2) # remove data from camera with a problem (too short sampling time)
dim(dat_cr)
dim(subset(dat_cr,is.na(hummingbird_species)))
#add transect data:
dim(dat_cr)
dat_cr=merge(dat_cr,tr_cr,by=c("site","month","year","plant_species"),all.x=T,all.y=F)
dim(dat_cr)
#add site data:
dim(dat_cr)
dat_cr=merge(dat_cr,sites_cr,by=c("site"),all.x=T,all.y=F)
dim(dat_cr)

#remove piercing interactions:
dat_cr$piercing[is.na(dat_cr$piercing)]="no"
dim(dat_cr)
dat_cr=subset(dat_cr,piercing!="yes")
dim(dat_cr)

#### CALCULATE INTERACTION FREQUENCIES
tab=dat_cr %>% dplyr::group_by(hummingbird_species,plant_species,waypoint,duration_sampling_hours,site,month,year,abond_flower,midpoint_Longitude,midpoint_Latitude,min_transect_elev) %>% dplyr::summarise(Y=length(time[piercing!="yes"])) #number of interaction detected
tab$Yfreq=tab$Y/tab$duration_sampling_hours*24*7 #divide by sampling time to get interaction frequencies in nb inter./hour
#### INFER ZEROS
mat=dcast(tab,site+year+month+waypoint+plant_species+duration_sampling_hours+abond_flower+midpoint_Longitude+midpoint_Latitude+min_transect_elev~hummingbird_species,value.var="Yfreq",fill=0)
tab=melt(mat,id.vars=c("site","year","month","waypoint","plant_species","duration_sampling_hours","abond_flower","midpoint_Longitude","midpoint_Latitude","min_transect_elev"),variable.name="hummingbird_species",value.name="Yfreq")
# CALCULATE TOTAL ABUNDANCE OF HUMMINGBIRDS PER SITES:
tab=tab %>% dplyr::group_by(hummingbird_species,site) %>% dplyr::mutate(abond_total_h=sum(Yfreq)) #number of interaction detected
#REMOVE HUMMINGBIRDS THAT ARE TOTALLY ABSENT FROM ONE SITE:
tab=subset(tab,abond_total_h>0)
#CREATE A MONTH_YEAR NUMERIC VARIABLES:
bidon=data.frame(year=rep(seq(min(tab$year),max(tab$year)),each=12),month=rep(1:12,length(unique(tab$year))))
bidon$num_time=1:nrow(bidon)
tab=merge(tab,bidon,by=c("year","month"))


getmode <- function(v) {
 uniqv <- unique(na.omit(v))
 uniqv[which.max(tabulate(match(v, uniqv)))]
}

wmean <- function(v,z) {
 uniqv <- v[!is.na(v) & !is.na(z)]
 uniqz <- z[!is.na(v) & !is.na(z)]
 if(length(uniqz>0)){return(sum(uniqv*uniqz)/sum(uniqz))}else{return(NA)}
}

#LOAD HUMMINGBIRD TRAIT DATA AND COMBINE THEM TO HAVE ONE VALUE PER SPECIES
tr=fread("C:/Users/Duchenne/Documents/EPHI_data_clean/hummingbird_traits_2023-01-24/Hummingbird_traits.txt",na.strings = c("",NA))
tr$tail_length=as.numeric(as.character(tr$tail_length))
tr1=tr %>% dplyr::group_by(hummingbird_species) %>% dplyr::summarise(bill_length=wmean(bill_length,N),culmen_length=wmean(culmen_length,N),tail_length=wmean(tail_length,N))

#LOAD PLANT TRAIT DATA AND COMBINE THEM TO HAVE ONE VALUE PER SPECIES
tr=fread("C:/Users/Duchenne/Documents/EPHI_data_clean/plant_traits_2023-01-24/Plant_traits.txt",na.strings = c("",NA))
tr2=subset(tr,!is.na(plant_species))  %>% group_by(plant_species,plant_family,plant_genus) %>% summarise(Tube_length=mean(Tube_length*10,na.rm=T),Anther_length=mean(Anther_length*10,na.rm=T), Stigma_length=mean(Stigma_length*10,na.rm=T),
Opening_corrolla=mean(Opening_corrolla,na.rm=T),Curvature_middle=mean(Curvature_middle,na.rm=T),Type=getmode(Type),sexualsys=getmode(SexualSystem))

#MERGE TRAITS AND DATA
tab=merge(tab,tr1,by=c("hummingbird_species"),all.x=TRUE,all.y=FALSE)
tab=merge(tab,tr2,by=c("plant_species"),all.x=TRUE,all.y=FALSE)

### CALCULATE PHENOLOGICAL INDICES:
#hummingbird phenologies from interaction
phenoh=subset(tab, !is.na(hummingbird_species)) %>% group_by(hummingbird_species,month,site) %>% summarise(Ys=sum(Yfreq))
phenoh_mat=dcast(phenoh,site+month~hummingbird_species,value.var="Ys",fill=0)
phenoh=melt(phenoh_mat,id.vars=c("site","month"),variable.name="hummingbird_species")
phenoh=phenoh %>% group_by(site,hummingbird_species) %>% mutate(maxpheno=max(value))
phenoh=phenoh %>% group_by(site,month,hummingbird_species) %>% summarise(phenoh=value/maxpheno)
phenoh$phenoh[is.na(phenoh$phenoh)]=0

#plant phenologies based on transects
plant_res=merge(plant_for_res,tr_cr,by=c("plant_species","month","year","site"),all=T)
plant_res$abond_flower[is.na(plant_res$abond_flower)]=1
plant_res=merge(plant_res,tr2[,c("plant_species","Tube_length")],by=c("plant_species"),all.x=T)
plant_res=subset(plant_res,!is.na(plant_species))

phenop=plant_res %>% group_by(plant_species,site,month,Tube_length) %>% summarise(value=sum(abond_flower))
bidon=expand.grid(month=1:12,site=unique(phenop$site))
phenop=merge(phenop,bidon,by=c("site","month"))
phenop=phenop %>% group_by(plant_species,site) %>% mutate(sumpheno=sum(value))
phenop=phenop %>% group_by(site,month,plant_species,Tube_length) %>% summarise(phenop=value/sumpheno)
fwrite(phenop,"plants_per_site_per_month.txt")


dim(tab)
tab=merge(tab,phenoh,by=c("site","month","hummingbird_species"),all.x=T,all.y=F)
dim(tab)
tab=merge(tab,phenop[,c("site","month","plant_species","phenop")],by=c("site","month","plant_species"),all.x=T,all.y=F)
dim(tab)

#### EMPIRICAL NETWORKS
emp_net=dat_cr %>% dplyr::group_by(hummingbird_species,plant_species,site) %>% dplyr::summarise(Y=length(time)/sum(duration_sampling_hours[!duplicated(waypoint)])) #number of interaction detected
emp_net$Y=emp_net$Y*24*7
emp_net=subset(emp_net,!is.na(hummingbird_species))
fwrite(emp_net,"empirical_networks.txt")

#EXPORT DATASET
fwrite(tab,"data_for_analyses.txt")


