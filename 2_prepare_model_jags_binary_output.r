setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper/data")
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
barrier_infer[j] ~ dnorm(bill_lengthu[j],1/(0.2*bill_lengthu[j]*0.2*bill_lengthu[j]))T(bill_lengthu[j],2*bill_lengthu[j])
#barrier_infer[j] ~ dunif(bill_lengthu[j],2*bill_lengthu[j])
match_infer[j] ~ dnorm(bill_lengthu[j],1/(0.2*bill_lengthu[j]*0.2*bill_lengthu[j]))T(0.5*bill_lengthu[j],1.5*bill_lengthu[j])
#match_infer[j] ~ dunif(0.5*bill_lengthu[j],1.5*bill_lengthu[j])
}

# Process model
for(i in 1:N){
	#new variables:
	barrier_var[i] <- ifelse(Tube_length[i] > barrier_infer[hummingbird_num[i]], 1, 0)
	mismatch_var[i] <- abs(Tube_length[i]-match_infer[hummingbird_num[i]])
	
	# Zero inflation:
	logit(pz[i]) <- Interceptpz + (-20* traitBarrier * barrier_var[i])
	Z[i] ~ dbern(min(pz[i]+0.00000000000000001,0.9999999999999999))
	
	# Interaction frequency:
	log(lambda[i]) <- Intercept + traitMismatch * mismatch_var[i] + pheno * phenoh[i]+abond *abond_flower_log[i] +plant_effect[plant_num[i]]+site_effect[site_num[i]]+temp_effect[site_num[i],num_time[i]]+sitebird_effect[hummingbird_species_site_num[i]]
	p[i] <- r/(r+(lambda[i] * Z[i]))
	Yfreq_r[i] ~ dnegbin(min(p[i]+0.00000000000000001,0.9999999999999999),r)
	#Yfreq_r[i] ~ dpois(lambda[i] * Z[i] + 0.00000000000000001)
	
}

# PRIORS FIXED EFFECTS:
Interceptpz ~ dnorm(0, 0.5)T(-10,10)
Intercept ~ dnorm(0, 0.01) 
traitMismatch ~ dnorm(0, 0.01)T(,0)
pheno ~ dnorm(0, 0.01)
abond ~ dnorm(0, 0.01)
traitBarrier ~ dbern(0.5)

# PRIORS RANDOM EFFECTS:
#spatial autocorrelation:
edec ~ dgamma(3, 0.1)
for(j in 1:Nsite){
for(i in 1:Nsite){
D.covar[i,j] <- exp(-edec*Distance[i,j])
}}

site_effect[1:Nsite] ~ dmnorm.vcov(rep(0,Nsite), spatialvariance * D.covar[1:Nsite,1:Nsite])
spatialvariance=sd.site * sd.site

for(i in 1:Nsite){
#site_effect[i] ~ dnorm(0, tau.site)
temp_effect[i,1] ~ dnorm(0, tau.temp)
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
r ~ dnegbin(0.2,4)

}
"

writeLines(model_string,con="model_binary.txt")


ParsStage <- c("barrier_infer","match_infer","Interceptpz","traitBarrier","Intercept","traitMismatch","pheno","abond","plant_effect","site_effect","temp_effect","sitebird_effect","sd.plant","sd.bird","sd.site","sd.temp","r","edec")

Inits <- function(){list()}

t1=Sys.time()
results1 <- jags(model.file="model_binary.txt", parameters.to.save=ParsStage, n.chains=1, data=dat,n.burnin = 50,  n.iter = 100, n.thin = 2,inits =Inits,jags.seed =2)
t2=Sys.time()
t2-t1

obj1=as.mcmc(results1)
plot(obj1[,"traitBarrier"])