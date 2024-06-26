###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("randomForest","data.table", "dplyr", "lubridate","lme4","R2jags","mcmcOutput","mcmcplots","MCMCvis",
			"pastecs","ggplot2","cowplot","gridExtra","scales","reshape2","bipartite","stringr","ungeviz","lme4","mgcv","ggpubr","emmeans","ggeffects","viridis") 


inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)

EPHI_version="2024-03-18"

setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper/data")
couleurs=c("#679436","#0B4F6C","deeppink")

######################## TRAITS
humm_table=fread("humm_table.txt")

ini_net=NULL
extinctions=NULL
for(pays in c("Costa-Rica","Ecuador","Brazil")){
ini_netc=fread(paste0("initial_network_",pays,".txt"))
ini_netc$Country=pays
ini_net=rbind(ini_net,ini_netc)
}
ini_net$Country=gsub("-"," ",ini_net$Country,fixed=T)

bp=unique(ini_net[,c("site","Country","plant_species","Tube_length","min_transect_elev")])
bp$type="plants"
names(bp)[3:4]=c("species","trait")
bh=unique(ini_net[,c("site","Country","hummingbird_species","culmen_lengthu","min_transect_elev")])
bh$type="birds"
names(bh)[3:4]=c("species","trait")

bf=rbind(bh,bp)
bf$type=as.factor(bf$type)
bf$facet=bf$Country

model=lmer(trait~min_transect_elev*type*Country+(1|site),data=bf)
obj=as.data.frame(emtrends(model,specs=c("type","Country"),var="min_transect_elev"))
obj$signi="no"
obj$signi[obj$upper.CL<0]="yes"
obj$signi[obj$lower.CL>0]="yes"

pre=ggpredict(model,c("min_transect_elev","type","Country"))
lims=bf %>% group_by(Country) %>% summarise(mini=min(min_transect_elev),maxi=max(min_transect_elev))
pre=merge(pre,lims,by.x="facet",by.y="Country")
pre=merge(pre,obj,by.x=c("group","facet"),by.y=c("type","Country"))
pre[pre$x<pre$mini | pre$x>pre$maxi,c("conf.high","conf.low","predicted")]=NA


tabs21=as.data.frame(summary(model)$coefficients)[,1:2]
tabs21$varia=rownames(tabs21)
tabs21$response="Trait values"

pl1=ggplot()+
geom_jitter(data=bf,aes(x=min_transect_elev,y=trait,col=type,fill=type),alpha=0.3,width = 20,size=0.5)+
geom_ribbon(data=pre,aes(x=x,y=predicted,ymin=conf.low,ymax=conf.high,fill=group),alpha=0.2)+
geom_line(data=pre,aes(x=x,y=predicted,col=group,linetype=signi),size=1.3)+
xlab("Elevation (m)")+ylab("Corolla or Bill length (mm)")+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background = element_blank())+
scale_color_manual(values=c("darkorchid4","chartreuse4"))+scale_fill_manual(values=c("darkorchid4","chartreuse4"))+ggtitle("a")+labs(color="",fill="")+
scale_linetype_manual(values=c("dashed","solid"))+guides(linetype="none")+facet_wrap(~facet)+
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))

pl1b=ggplot()+
geom_density(data=bf,aes(x=trait,col=min_transect_elev,fill=min_transect_elev,group=min_transect_elev),alpha=0.1)+
xlab("Bill or Corolla length(mm)")+ylab("Density")+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),strip.background = element_blank())+
ggtitle("")+labs(color="Elev.",fill="Elev.")+
facet_grid(cols=vars(Country),rows=vars(type))+scale_color_viridis(option="mako")+scale_fill_viridis(option="mako")+
coord_cartesian(expand=F)

setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper")
png("Figure_S4.png",width=1200,height=1000,res=160)
pl1b
dev.off();

####Panel b
forb=as.data.frame(ini_net %>% group_by(site,min_transect_elev,Country) %>%
summarise(compl=mean(-1*abs(trait_plant-match_infer)),prop_forbidden=mean(average_proba_without_barrier-average_proba),
nbh=length(unique(hummingbird_species)),nbp=length(unique(plant_species))))
setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper/data")
fwrite(forb,"network_match_barrier.txt")

model=glm(prop_forbidden~min_transect_elev*Country,data=forb,family=quasibinomial)
obj=as.data.frame(emtrends(model,specs=c("Country"),var="min_transect_elev"))
obj$signi="no"
obj$signi[obj$asymp.UCL<0]="yes"
obj$signi[obj$asymp.LCL>0]="yes"
obj$signi=factor(obj$signi,levels=c("yes","no"))


pre=ggpredict(model,c("min_transect_elev[all]","Country"))
lims=bf %>% group_by(Country) %>% summarise(mini=min(min_transect_elev),maxi=max(min_transect_elev))
pre=merge(pre,lims,by.x="group",by.y="Country")
pre=merge(pre,obj,by.x=c("group"),by.y=c("Country"))
pre[pre$x<pre$mini | pre$x>pre$maxi,c("conf.high","conf.low","predicted")]=NA
pre$group=factor(pre$group,levels=c("Brazil","Costa Rica","Ecuador"))

tabs23=as.data.frame(summary(model)$coefficients)[,1:2]
tabs23$varia=rownames(tabs23)
tabs23$response="Proportion of forbidden links"
fwrite(rbind(tabs21,tabs23),"Table_S2.csv")

pl3=ggplot()+
geom_point(data=forb,aes(y=prop_forbidden,x=min_transect_elev,col=Country),size=1.5)+
geom_ribbon(data=pre,aes(x=x,y=predicted,ymin=conf.low,ymax=conf.high,fill=group),alpha=0.2)+
geom_line(data=pre,aes(x=x,y=predicted,col=group,linetype=signi),size=1.3)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+
labs(col="",fill="")+ggtitle("b")+ylab("Proportion of forbidden links")+xlab("Elevation (m)")+
scale_color_manual(values=couleurs)+scale_fill_manual(values=couleurs)+
scale_linetype_manual(values=c("solid","dashed"))+guides(linetype="none")


setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper")
pdf("Figure3.pdf",width=6,height=6)
plot_grid(pl1,pl3,ncol=1,align="hv")
dev.off();


################FIG S1

bidon=melt(forb,id.vars=c("site","min_transect_elev","compl","prop_forbidden","Country"))
gammodel=gam(value~s(min_transect_elev,by=as.factor(Country)),data=subset(bidon,variable=="nbh"),family=poisson)
newdata=data.frame(min_transect_elev=rep(10:3406,3),Country=rep(c("Brazil","Costa Rica","Ecuador"),each=3397),variable="birds")
pre1=cbind(newdata,as.data.frame(predict(gammodel,se=TRUE,type="response",newdata=newdata)))

gammodel=gam(value~s(min_transect_elev,by=as.factor(Country)),data=subset(bidon,variable=="nbp"),family=poisson)
newdata=data.frame(min_transect_elev=rep(10:3406,3),Country=rep(c("Brazil","Costa Rica","Ecuador"),each=3397),variable="plants")
pre2=cbind(newdata,as.data.frame(predict(gammodel,se=TRUE,type="response",newdata=newdata)))
fit=rbind(pre1,pre2)
fit=merge(fit,lims,by="Country")
fit[fit$min_transect_elev<fit$mini | fit$min_transect_elev>fit$maxi,c("fit","se.fit")]=NA
bidon$variable=as.character(bidon$variable)
bidon$variable[bidon$variable=="nbp"]="plants"
bidon$variable[bidon$variable=="nbh"]="birds"

pl1=ggplot()+
geom_ribbon(data=fit,aes(x=min_transect_elev,col=variable,fill=variable,ymin=fit-1.96*se.fit,ymax=fit+1.96*se.fit),alpha=0.2,col=NA)+
geom_line(data=fit,aes(x=min_transect_elev,col=variable,fill=variable,y=fit))+
geom_point(data=bidon,aes(y=value,x=min_transect_elev,col=variable,fill=variable),size=1.5)+
xlab("Elevation (m)")+ylab("Species richness")+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),strip.background = element_blank())+
scale_color_manual(values=c("darkorchid4","#C8CC92"))+scale_fill_manual(values=c("darkorchid4","#C8CC92"))+ggtitle("a")+labs(color="",fill="")+
facet_wrap(~Country,ncol=1)

forb$symmetrie=forb$nbp/forb$nbh
model=glm(symmetrie~poly(min_transect_elev,2)*as.factor(Country),data=forb,family=quasipoisson)
newdata=data.frame(min_transect_elev=rep(10:3406,3),Country=rep(c("Brazil","Costa Rica","Ecuador"),each=3397))
fit=cbind(newdata,as.data.frame(predict(model,se=TRUE,type="response",newdata=newdata)))
fit=merge(fit,lims,by="Country")
fit[fit$min_transect_elev<fit$mini | fit$min_transect_elev>fit$maxi,c("fit","se.fit")]=NA

pl2=ggplot()+
geom_ribbon(data=fit,aes(x=min_transect_elev,col=Country,fill=Country,ymin=fit-1.96*se.fit,ymax=fit+1.96*se.fit),alpha=0.2,col=NA)+
geom_line(data=fit,aes(x=min_transect_elev,col=Country,fill=Country,y=fit))+
geom_point(data=forb,aes(x=min_transect_elev,y=symmetrie,col=Country,fill=Country),alpha=0.8)+
stat_smooth(method="lm")+
xlab("Elevation (m)")+ylab("Number of plant species / number of hummingbird species")+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),strip.background = element_blank())+
scale_color_manual(values=couleurs)+scale_fill_manual(values=couleurs)+ggtitle("b")+labs(color="",fill="")


setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper")
png("Figure_S1.png",width=1500,height=1000,res=160)
grid.arrange(pl1,pl2,ncol=2,widths=c(1,1))
dev.off();


