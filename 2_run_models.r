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
library(R2jags)
setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper/data")

tab=fread("data_for_analyses.txt")

#### CALCULATE TRAIT MATCHING:
tab$complementarity=tab$Tube_length-tab$bill_length
hist(tab$complementarity)
tab$barrier=0
tab$barrier[tab$complementarity>(0.25*tab$bill_length)]=1
hist(tab$barrier)
tab$mismatch=scale(abs(tab$complementarity)) #trait mismatch, scaled for convergence
tab$elev=scale(tab$min_transect_elev)

### CALCULATE PHENOLOGICAL MATCH:
phenoh=subset(tab, !is.na(hummingbird_species)) %>% group_by(hummingbird_species,month) %>% summarise(Ys=sum(Yfreq))
phenoh_mat=dcast(phenoh,month~hummingbird_species,value.var="Ys",fill=0)
func=function(x){
t(t(x[,-1])/apply(x[,-1],2,sum))
}

phenoh_mat[,-1]=func(phenoh_mat)
phenoh=melt(phenoh_mat,id.vars="month",variable.name="hummingbird_species")

phenop=subset(tab, !is.na(plant_species)) %>% group_by(plant_species,month) %>% summarise(Ys=sum(Yfreq))
phenop_mat=dcast(phenop,month~plant_species,value.var="Ys",fill=0)
phenop_mat[,-1]=func(phenop_mat)
phenop=melt(phenop_mat,id.vars="month",variable.name="plant_species")

pheno=merge(phenoh,phenop,by="month")

pheno_match=pheno %>% group_by(hummingbird_species,plant_species) %>% summarise(pheno_match=sum(value.x*value.y))
pheno_match$pheno_barrier=0
pheno_match$pheno_barrier[pheno_match$pheno_match==0]=1
#merge pheno_match with data
tab=merge(tab,pheno_match,by=c("hummingbird_species","plant_species"),all.x=T,all.y=F)
tab=merge(tab,phenoh,by=c("hummingbird_species","month"),all.x=T,all.y=F)

#REMOVE DATA WITHOUT TRAITS
tab=subset(tab,!is.na(complementarity))

#PREPARE DATA FOR TEMPORAL AND SPATIAL AUTOCORRELATION:
tab$site <- factor(tab$site)
tab$times <- factor(tab$num_time, levels=sort(unique(tab$num_time)))
tab$pos <- numFactor(tab$midpoint_Longitude, tab$midpoint_Latitude)
tab$group <- factor(rep(1, nrow(tab)))
tab$hummingbird_species_site=paste0(tab$hummingbird_species,tab$site)

#EVENTUALLY ROUND INTERACTION FREQUENCIES
tab$Yfreq_r=round(tab$Yfreq)

# RUN GENERAL MODEL
model=glmmTMB(Yfreq~mismatch+value+ar1(times +0 | pos)+(1|hummingbird_species_site)+(1|plant_species),data=tab,family=nbinom1,ziformula = ~as.factor(barrier),control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS")))

Anova(model)

save(model,file=paste0("model_big_glmmTMB.RData"))


predict_all=predict(model,se.fit=F,type="response",re.form = NULL)
cor(tab$Yfreq,predict_all)
tab2=tab
tab2$mismatch=mean(tab$mismatch)
predict_barrier=predict(model,se.fit=F,type="response",re.form = NULL,newdata=tab2)
cor(tab$Yfreq,predict_barrier)
tab2=tab
tab2$barrier=0
predict_mismatch=predict(model,se.fit=F,type="response",re.form = NULL,newdata=tab2)
cor(tab$Yfreq,predict_mismatch)


#DEFINE MODEL FORMULA:
f.s=formula(paste0("y~-1+Intercept+mismatch+barrier+f(field,model=spde)+f(num_time,model='rw1',group=as.numeric(as.factor(site)))+f(plant_species,model='iid')+f(hummingbird_species_site,model='iid')"))

# RUN MODEL
model <- inla(y~mismatch+value+f(num_time,model='rw1',group=as.numeric(as.factor(site)))+f(plant_species,model='iid')+f(hummingbird_species_site,model='iid'), family = "xpoisson", data =tab)

summary(model)

tab$hummingbird_species_site_num=as.numeric(as.factor(tab$hummingbird_species_site))
tab$hummingbird_num=as.numeric(as.factor(tab$hummingbird_species))
tab$plant_num=as.numeric(as.factor(tab$plant_species))
tab$site_num=as.numeric(as.factor(tab$site))
tab$num_time=as.numeric(as.factor(tab$num_time))

tabu=unique(tab[,c("hummingbird_num","bill_length")])
names(tabu)=paste0(names(tabu),"u")

dat=c(as.list(tab),as.list(tabu),list(N=nrow(tab),
Nbirds=length(unique(tab$hummingbird_num)),Nbird_site=length(unique(tab$hummingbird_species_site_num)),Nsites=length(unique(tab$site_num)),Nplants=length(unique(tab$plant_num)),
Ntemp=length(unique(tab$num_time)),Nsite=length(unique(tab$site_num))))

save(dat,file=paste0("data_formodel.RData"))

######################################################################################################
######################################################################################################
#### JAGS model file written by runjags version 2.2.0-2 on 2021-11-23 15:22:45 
######################################################################################################
######################################################################################################

### Model template as follows - ensure this is syntactically correct before running the model!

model_string="
model{

# Barrier assessment
for(j in 1:Nbirds){
barrier[j] ~ dnorm(bill_lengthu[j],1/(0.2*bill_lengthu[j]*0.2*bill_lengthu[j]))T(bill_lengthu[j],2*bill_lengthu[j])
match[j] ~ dnorm(bill_lengthu[j],1/(0.2*bill_lengthu[j]*0.2*bill_lengthu[j]))T(0.5*bill_lengthu[j],1.5*bill_lengthu[j])
}

# Process model
for(i in 1:N){
	#new variables:
	barrier_var[i] <- ifelse(corolla_length[i] > barrier[hummingbird_num[i]], 1, 0)
	mismatch_var[i] <- abs(corolla_length[i]-match[hummingbird_num[i]])
	
	# Zero inflation:
	logit(pz[i]) <- Interceptpz + traitBarrier * barrier_var[i]
	Z[i] ~ dbern(min(pz[i]+0.00000000000000001,0.9999999999999999))
	
	# Interaction frequency:
	log(lambda[i]) <- Intercept + traitMatch * mismatch_var[i] + pheno * value[i] +plant_num[plant_num[i]]+temp_effect[site_num[i],num_time[i]]+sitebird_effect[hummingbird_species_site_num[i]]
	p[i] <- r/(r+lambda[i])
	p2[i] <- p[i] * Z[i]
	Yfreq_r[i] ~ dnegbin(min(p2[i]+0.00000000000000001,1),r)
	
}

# PRIORS FIXED EFFECTS:
Interceptpz ~ dnorm(0, 0.5)
Intercept ~ dnorm(0, 0.5) 
traitMatch ~ dnorm(0, 0.5)
pheno ~ dnorm(0, 0.5)
traitBarrier ~ dnorm(0, 0.5)

# PRIORS RANDOM EFFECTS:
for(i in 1:Nsite){
temp_effect[i,1] ~ dnorm(0, tau.site)
for(j in 2:Ntemp){
temp_effect[i,j] ~ dnorm(temp_effect[i,(j-1)], tau.temp)
}
}

for(j in 1:Nbird_site){
sitebird_effect[j] ~ dnorm(0, tau.bird)
}

for(j in 1:Nplants){
plant_effect[j] ~ dnorm(0, tau.plant)
}
tau.plant <- 1/(sd.plant * sd.plant)
sd.plant ~ dt(0, 1, 1)T(0,)

#hyperpriors:
tau.site<- 1/(sd.site * sd.site)
sd.site ~ dt(0, 1, 1)T(0,)
tau.temp<- 1/(sd.temp * sd.temp)
sd.temp ~ dt(0, 1, 1)T(0,)
tau.bird=1/(sd.bird * sd.bird)
sd.bird ~ dt(0, 1, 1)T(0,)


# PRIOR OVERDISPERTION:
r ~ dunif(0,50)

}
"

writeLines(model_string,con="model.txt")

ParsStage <- c("barrier","match","Interceptpz","traitbarrier","traitMatch","sitebird_effect","plant_effect","sd.plant","sd.bird","sd.site","temp_effect","r")

Inits <- function(){list()}

results1 <- jags(model.file="model.txt", parameters.to.save=ParsStage, n.chains=1, data=dat,n.burnin = 2500,  n.iter = 3000, n.thin = 2,inits =Inits)

save(results1,file=paste0("chain_model_modelphylo_simpl_",j,".RData"))
#


























































###### LOOP TO ASSESS OPTIMAL TRAIT MATCH PER SPECIES
setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper/data/matches")
resf=NULL
for(i in 1:21){
bidon=fread(paste0("matches_",i,".txt"))
resf=rbind(resf,bidon)
print(i)
}

resf=merge(resf,liste,by=c("hummingbird_species"))

resf=resf %>% group_by(hummingbird_species) %>% mutate(AIC_s=AIC-AIC[1],loglik_s=loglik-loglik[1])

ggplot(data=resf,aes(color=hummingbird_species,x=tongue,y=AIC_s))+geom_line()

b=resf %>% group_by(hummingbird_species,bill_length) %>% summarise(tongue_opt=tongue[which.min(AIC)])

ggplot(data=b,aes(y=tongue_opt,x=bill_length))+geom_point()







model <- inla(round(Yfreq)~mismatch+f(times,model="rw1")+f(site,model="iid"),data=bidon,family="nbinomial")







#MAKE A SPATIAL DATA FRAME FROM THE DATA:
sfpts=sf::st_as_sf(tab, coords = c("midpoint_Longitude","midpoint_Latitude"),crs=CRS("+init=epsg:4326"))
spts=as(sfpts, "Spatial")


###MODEL PREPARATION, SPATIAL AUTOCORRELATION:

# create a simple mesh using the site locations
clm.mesh <- inla.mesh.2d(loc=spts, max.edge = c(1,20)) #max.edge and cut.off can be used to define traingle size in the mesh. Using smaller triangles increases precision but also exponentially increases computing power. 
#PLot the mesh
ggplot() +
  gg(clm.mesh) +
  geom_sf(data=sfpts,col='purple',size=1.7,alpha=0.5) 

# After the mesh has been set up, we need to feed INLA a way to convert this into a model format. This uses an A matrix, which essentially translates spatial locations on the mesh into vectors in the model:
A <- inla.spde.make.A(clm.mesh, loc = as.matrix(tab[,c("midpoint_Longitude","midpoint_Latitude")]))
spde <- inla.spde2.matern(clm.mesh, alpha = 2)
mesh.index <- inla.spde.make.index(name = "field", n.spde = spde$n.spde)

# LIST ALL VARIABLES WE NEED AND GROUP THEM IN A INLA STACK OBJECT:
listevar=list(plant_species = tab$plant_species,num_time=tab$num_time,site=tab$site,mismatch=abs(tab$complementarity),barrier=as.numeric(tab$barrier),hummingbird_species_site=paste0(tab$hummingbird_species,tab$site))
stk.dat <- inla.stack(data = list(y = tab$Yfreq), A = list(A, 1), tag = "est", effects = list(c(mesh.index, list(Intercept = 1)),
listevar))

#DEFINE MODEL FORMULA:
f.s=formula(paste0("y~-1+Intercept+mismatch+barrier+f(field,model=spde)+f(num_time,model='rw1',group=as.numeric(as.factor(site)))+f(plant_species,model='iid')+f(hummingbird_species_site,model='iid')"))

# RUN MODEL
model <- inla(f.s, family = "xpoisson", data = inla.stack.data(stk.dat), verbose = TRUE, 
    control.predictor = list(A = inla.stack.A(stk.dat), compute = TRUE))

summary(model)


#VARIANCE PARTIONING
vp.plot(model, vp.model = "striid", table.type = "mixing", cex.sources = 0.95, title.plot = "", nsim = 1000)
vp.table(model, vp.model = "striid", table.type = "mixing", nsim = 1000)

# save model result
save(model,file=paste0("model_COSTA_RICA.RData"))




tab=tab %>% group_by(plant_species) %>% mutate(nb_sampl=length(unique(waypoint)))








liste=tab %>% group_by(hummingbird_species,nb_sampl) %>% summarise(Ysum=sum(Yfreq))
liste=subset(liste,Ysum>0 & nb_sampl>25)