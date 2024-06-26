###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("randomForest","data.table", "dplyr", "lubridate","lme4") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)

setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper/data")

EPHI_version="2024-03-18"
sitef=NULL
data_descf=NULL
for(pays in c("Costa-Rica","Ecuador","Brazil")){
#LOAD HUMMINGBIRD DATA FROM COSTA-RICA:
dat=fread(paste0("C:/Users/Duchenne/Documents/EPHI_data_clean/",pays,"_",EPHI_version,"/Interactions_data_",pays,".txt"),na.strings = c("",NA))
dat[dat==""]=NA
dat$date=as.IDate(dat$date,"%Y/%m/%d") #be sure date columns is recognize as date
dat$piercing[is.na(dat$piercing)]="no"
#LOAD CAMERA INFORMATION:
cameras=fread(paste0("C:/Users/Duchenne/Documents/EPHI_data_clean/",pays,"_",EPHI_version,"/Cameras_data_",pays,".txt"),na.strings = c("",NA))
cameras$end_date=as.IDate(cameras$end_date,"%Y/%m/%d") #be sure date columns is recognize as date
cameras$start_date=as.IDate(cameras$start_date,"%Y/%m/%d") #be sure date columns is recognize as date
cameras$month=month(cameras$start_date) #extract month from date column
cameras$year=year(cameras$start_date) #extract year from date column
cameras=subset(cameras,year<2023)
plant_for_res=unique(cameras[,c("plant_species","month","year","site")])
#LOAD TRANSECT DATA:
transects=fread(paste0("C:/Users/Duchenne/Documents/EPHI_data_clean/",pays,"_",EPHI_version,"/Transect_data_",pays,".txt"),na.strings = c("",NA))
transects$date=as.IDate(transects$date,"%Y/%m/%d") #be sure date columns is recognize as date
transects$month=month(transects$date) #extract month from date column
transects$year=year(transects$date) #extract year from date column
transects=subset(transects,!is.na(plant_species)) %>% group_by(site,year,month,plant_species) %>% summarise(abond_flower=sum(total_flowers,na.rm=T),
																								abond_plant=length(total_flowers)) #calculate a total abundance per plant species per month per site
transects$abond_flower[transects$abond_flower==0]=1
transects$combi_dates=apply(transects[,c("site","month","year")],1,paste,collapse="-")
plant_for_res$combi_dates=apply(plant_for_res[,c("site","month","year")],1,paste,collapse="-")
plant_for_res$duringEPHI="no"
plant_for_res$duringEPHI[plant_for_res$combi_dates %in% transects$combi_dates]="yes"
transects=merge(transects,plant_for_res,by=c("plant_species","month","year","site","combi_dates"),all=T)
transects$abond_flower[is.na(transects$abond_flower) & transects$duringEPHI=="yes"]=1
transects$abond_plant[is.na(transects$abond_plant) & transects$duringEPHI=="yes"]=1

#LOAD SITE METADATA:
sites=fread(paste0("C:/Users/Duchenne/Documents/EPHI_data_clean/",pays,"_",EPHI_version,"/Site_metadata_",pays,".txt"),na.strings = c("",NA))
sites=subset(sites,habitat!="deforested")

#MERGE THEM:
#interactions and cameras:
dim(dat)
dat=merge(dat,cameras,by=c("waypoint","site"),all.x=F,all.y=T) #we want to keep cameras that did not detect any hummingbird
dim(dat)
dat=subset(dat,!is.na(plant_species)) #remove data with no plant ID
dim(dat)
dat=subset(dat,duration_sampling_hours>2) # remove data from camera with a problem (too short sampling time)
dim(dat)
dim(subset(dat,is.na(hummingbird_species)))
#add transect data:
dim(dat)
dat=merge(dat,transects,by=c("site","month","year","plant_species"),all.x=T,all.y=F)
dim(dat)
#add site data:
dim(dat)
dat=merge(dat,sites,by=c("site"),all=FALSE)
dim(dat)

dat=subset(dat,year(start_date)<=2022)

#Duration from picture when non-available from cameras
dat$duration_sampling_hours[is.na(dat$duration_sampling_hours)]=dat$duration_from_pics[is.na(dat$duration_sampling_hours)]
dat=subset(dat,duration_sampling_hours>=5 & camera_problem!="yes")

#### CALCULATE INTERACTION FREQUENCIES
tab=dat %>% dplyr::group_by(hummingbird_species,plant_species,waypoint,duration_sampling_hours,site,month,year,abond_flower,
midpoint_Longitude,midpoint_Latitude,min_transect_elev,Country) %>% dplyr::summarise(Y=length(time[piercing!="yes"]),Y_pierc=length(time)) #number of interaction detected
total_inter=sum(tab$Y_pierc)
data_desc=data.frame(country=pays,nb_inter=sum(tab$Y),nb_inter_pierc=sum(tab$Y_pierc))
data_descf=rbind(data_descf,data_desc)
#### INFER ZEROS
mat=dcast(tab,site+year+month+waypoint+plant_species+duration_sampling_hours+abond_flower+midpoint_Longitude+
midpoint_Latitude+min_transect_elev+Country~hummingbird_species,value.var="Y",fill=0)
tab=melt(mat,id.vars=c("site","year","month","waypoint","plant_species","duration_sampling_hours",
"abond_flower","midpoint_Longitude","midpoint_Latitude","min_transect_elev","Country"),variable.name="hummingbird_species",value.name="Y")
# CALCULATE TOTAL ABUNDANCE OF HUMMINGBIRDS PER SITES:
tab=tab %>% dplyr::group_by(hummingbird_species,site) %>% dplyr::mutate(abond_total_h=sum(Y)) #number of interaction detected
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
tr=fread(paste0("C:/Users/Duchenne/Documents/EPHI_data_clean/hummingbird_traits_",EPHI_version,"/Hummingbird_traits.txt"),na.strings = c("",NA))
tr$tail_length=as.numeric(as.character(tr$tail_length))
tr1=tr %>% dplyr::group_by(hummingbird_species) %>% dplyr::summarise(nindh=length(culmen_length[!is.na(culmen_length)]),bill_length=wmean(bill_length,N),culmen_length=wmean(culmen_length,N),
tail_length=wmean(tail_length,N))
model=lm(culmen_length~bill_length,data=tr1[!is.na(tr1$bill_length) & !is.na(tr1$culmen_length),])
tr1$culmen_length[is.na(tr1$culmen_length) & !is.na(tr1$bill_length)]=
predict(model,newdata=tr1[is.na(tr1$culmen_length) & !is.na(tr1$bill_length),])

#LOAD PLANT TRAIT DATA AND COMBINE THEM TO HAVE ONE VALUE PER SPECIES
tr=fread(paste0("C:/Users/Duchenne/Documents/EPHI_data_clean/plant_traits_",EPHI_version,"/Plant_traits.txt"),na.strings = c("",NA))
tr2=subset(tr,!is.na(plant_species))  %>% group_by(plant_species,plant_family,plant_genus) %>%
summarise(nindp=length(Tubelength[!is.na(Tubelength)]),Tubelength=mean(Tubelength*10,na.rm=T),StamenLength=mean(StamenLength*10,na.rm=T),
StyleLength=mean(StyleLength*10,na.rm=T),Opening_corolla_lateral=mean(Opening_corolla_lateral,na.rm=T),
CurvatureCentral=mean(CurvatureCentral,na.rm=T),FlowerType=getmode(FlowerType),sexualsys=getmode(SexualSystem))

length(unique(tr2$plant_species[tr2$plant_species %in% dat$plant_species]))
length(unique(dat$plant_species))

#MERGE TRAITS AND DATA
tab=merge(tab,tr1,by=c("hummingbird_species"),all.x=TRUE,all.y=FALSE)
tab=merge(tab,tr2,by=c("plant_species"),all.x=TRUE,all.y=FALSE)

### CALCULATE PHENOLOGICAL INDICES:
#hummingbird phenologies from interaction
phenoh=subset(tab, !is.na(hummingbird_species)) %>% group_by(hummingbird_species,month,site) %>% summarise(Ys=sum(Y/duration_sampling_hours))
phenoh_mat=dcast(phenoh,site+month~hummingbird_species,value.var="Ys",fill=0)
phenoh=melt(phenoh_mat,id.vars=c("site","month"),variable.name="hummingbird_species")
phenoh=phenoh %>% group_by(site,hummingbird_species) %>% mutate(maxpheno=max(value))
phenoh=phenoh %>% group_by(site,month,hummingbird_species) %>% summarise(phenoh=value/maxpheno)
phenoh$phenoh[is.na(phenoh$phenoh)]=0

#plant phenologies based on transects
plant_res=merge(transects,tr2[,c("plant_species","Tubelength")],by=c("plant_species"),all.x=T)
plant_res=subset(plant_res,!is.na(plant_species))
plant_res=plant_res %>% group_by(site) %>% mutate(n_visits=length(unique(paste(month,year))))
plant_res=plant_res %>% group_by(plant_species,site) %>% mutate(abond_flower_moy=sum(abond_flower,na.rm=T)/n_visits,abond_plant_moy=sum(abond_plant,na.rm=T)/n_visits)
subset(plant_res,plant_species=="Thibaudia inflata")
phenop=dcast(plant_res,plant_species+site+Tubelength+abond_plant_moy+abond_flower_moy~month,value.var="abond_flower",fun.aggregate=sum)
phenop=melt(phenop,id.vars=c("plant_species","site","Tubelength","abond_plant_moy","abond_flower_moy"),variable.name = "month")
phenop$value[is.na(phenop$value)]=0
phenop=phenop %>% group_by(plant_species,site) %>% mutate(sumpheno=sum(value,na.rm=T))
phenop=phenop %>% group_by(site,month,plant_species,abond_flower_moy,Tubelength) %>% summarise(phenop=value/sumpheno)

fwrite(phenop,paste0("plants_per_site_per_month_",pays,".txt"))


dim(tab)
tab=merge(tab,phenoh,by=c("site","month","hummingbird_species"),all.x=T,all.y=F)
dim(tab)
tab=merge(tab,phenop[,c("site","month","plant_species","phenop")],by=c("site","month","plant_species"),all.x=T,all.y=F)
dim(tab)
tab$total_inter=total_inter

#### EMPIRICAL NETWORKS
emp_net=tab %>% dplyr::group_by(hummingbird_species,plant_species,site,min_transect_elev,culmen_length,Tubelength) %>%
dplyr::summarise(Y=length(time)/sum(duration_sampling_hours[!duplicated(waypoint)])) #number of interaction detected
emp_net$Y=emp_net$Y*12
emp_net=subset(emp_net,!is.na(hummingbird_species))
fwrite(emp_net,paste0("empirical_networks_",pays,".txt"))

#EXPORT DATASET
fwrite(tab,paste0("data_for_analyses_",pays,".txt"))
}
setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper/data")
fwrite(data_descf,"nb_data.csv")



tabf=NULL
b=0
lili=list()
for(pays in c("Costa-Rica","Ecuador","Brazil")){
b=b+1
tab=fread(paste0("data_for_analyses_",pays,".txt"),na.string=c("",NA))

#### CALCULATE TRAIT MATCHING:
tab$mismatch=abs(tab$Tubelength-tab$culmen_length) #trait mismatch
tab$elev=scale(tab$min_transect_elev)

#### CALCULATE PHENOLOGICAL MATCHING:
tab$phenomatch=tab$phenoh*tab$phenop

#ABONDANCE FLOWER THAT ARE NA ARE NON DETECTED PLANT ON THE TRANSECT, set a low abundance value:
tab$abond_flower[is.na(tab$abond_flower) | tab$abond_flower==0]=1
tab$abond_flower_log=log(tab$abond_flower)

unique(tab$site)
#REMOVE DATA WITHOUT TRAITS
tab=subset(tab,!is.na(mismatch))
lili[[b]]=unique(tab$hummingbird_species)
unique(tab$site)
tabf=rbind(tabf,tab)

}
lili=na.omit(lili)
length(intersect(lili[[1]],lili[[2]]))
length(intersect(lili[[3]],lili[[2]]))
length(intersect(lili[[3]],lili[[1]]))

length(unique(tabf$plant_species))
length(unique(tabf$hummingbird_species))
sum(tabf$Y)
summary(tabf$nindp[!duplicated(tabf$plant_species)])
summary(tabf$nindh[!duplicated(tabf$hummingbird_species)])



