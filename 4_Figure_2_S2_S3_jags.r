###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
pkgs <- c("randomForest","data.table", "dplyr", "lubridate","lme4","R2jags","mcmcOutput","mcmcplots","MCMCvis",
			"pastecs","ggplot2","cowplot","gridExtra","scales","here") 


# Check if packages are already installed
inst <- pkgs %in% installed.packages()
# Install missing packages
if (any(!inst)) install.packages(pkgs[!inst])
# Load packages
pkg_out <- lapply(pkgs, require, character.only = TRUE)

# Set the working directory for the "here" package
here::i_am("EPHI_paper.Rproj")
# Source custom functions for Bayesian predictions
source(here("scripts/additional scripts", "function_to_predict_from_bayesian.r"))
# Create figure folder if not existing 
dir.create(paste0(here(),"/figures"), showWarnings = FALSE)

###########################################
# Load Bayesian Models
###########################################

# Function to load model data and chains for a country
load_country_data <- function(country_name, data_file, chain_files) {
  load(here("data_zenodo", data_file))
  dat_country <- dat
  
  chains <- lapply(chain_files, function(file) {
    load(here("data_zenodo", file))
    results1
  })
  
  mcmc_objects <- lapply(chains, as.mcmc)
  mco_country <- mcmcOutput(do.call(mcmc.list, mcmc_objects))
  
  summary_country <- summary(mco_country)
  summary_country$varia <- rownames(summary_country)
  
  list(data = dat_country, mco = mco_country, summary = summary_country)
}

# Load data for each country (it can take several minutes)
costa_rica <- load_country_data("Costa Rica",
                                 "data_formodel_Costa-Rica.RData",
                                 c("chain_model_Costa-Rica_1.RData", "chain_model_Costa-Rica_2.RData", "chain_model_Costa-Rica_3.RData"))

ecuador <- load_country_data("Ecuador",
                             "data_formodel_Ecuador.RData",
                             c("chain_model_Ecuador_1.RData", "chain_model_Ecuador_2.RData", "chain_model_Ecuador_3.RData"))

brazil <- load_country_data("Brazil",
                            "data_formodel_Brazil.RData",
                            c("chain_model_Brazil_1.RData", "chain_model_Brazil_2.RData", "chain_model_Brazil_3.RData"))

## Elements that will be re-use
# Define colors
couleurs <- c("#679436", "#0B4F6C", "deeppink")

# table with all hummingbird species across all countries:
tabsp=unique(rbind(as.data.frame(costa_rica$data[c("hummingbird_num","hummingbird_species","Country")]),as.data.frame(brazil$data[c("hummingbird_num","hummingbird_species","Country")]),as.data.frame(ecuador$data[c("hummingbird_num","hummingbird_species","Country")])))
tabsp$Country=gsub("-"," ",tabsp$Country)

##################### FIGURE 2 ##################

##################### Panel a ##################
pre1=predict_model(data=costa_rica$data, chains=costa_rica$mco,nsampling=5000,site=NA,bird=NA,plant=NA,trait_plant=NA,mismatch=seq(0,60,length.out=100),
barrier=NA,pheno=NA,abond=NA,type="frequency",random_effects=NA,month=NA,year=NA,duration=NA)
pre1$Country="Costa Rica"
pre2=predict_model(data=brazil$data,chains=brazil$mco,nsampling=5000,site=NA,bird=NA,plant=NA,trait_plant=NA,mismatch=seq(0,60,length.out=100),
barrier=NA,pheno=NA,abond=NA,type="frequency",random_effects=NA,month=NA,year=NA,duration=NA)
pre2$Country="Brazil"
pre3=predict_model(data=ecuador$data,chains=ecuador$mco,nsampling=5000,site=NA,bird=NA,plant=NA,trait_plant=NA,mismatch=seq(0,60,length.out=100),
barrier=NA,pheno=NA,abond=NA,type="frequency",random_effects=NA,month=NA,year=NA,duration=NA)
pre3$Country="Ecuador"

# merge predictions
pred_mism=rbind(pre1,pre2,pre3)

pl1=ggplot(data=pred_mism,aes(x=-mismatch,y=average_lambda,color=Country,fill=Country))+geom_ribbon(aes(ymin=lwr_lambda,ymax=upr_lambda),alpha=0.2,col=NA)+
geom_line(size=1)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),
plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none")+xlab("Trait complementarity")+coord_cartesian(expand=F)+
scale_y_continuous(name = bquote(atop("Average interaction frequency", (interactions.day^-1))),sec.axis = sec_axis(~., name ="" ))+
scale_color_manual(values=couleurs)+scale_fill_manual(values=couleurs)+
labs(col="",fill="")+ggtitle("a")+guides(y="none")

##################### Panel b ##################
pre1=predict_model(data=costa_rica$data, chains=costa_rica$mco,nsampling=5000,site=NA,bird=NA,plant=NA,trait_plant=NA,mismatch=NA,barrier=c(0,1),pheno=NA,abond=NA,
type="barrier",random_effects=NA,month=NA,year=NA,duration=NA)
pre1$Country="Costa Rica"
pre2=predict_model(data=brazil$data,chains=brazil$mco,nsampling=5000,site=NA,bird=NA,plant=NA,trait_plant=NA,mismatch=NA,barrier=c(0,1),pheno=NA,abond=NA,
type="barrier",random_effects=NA,month=NA,year=NA,duration=NA)
pre2$Country="Brazil"
pre3=predict_model(data=ecuador$data,chains=ecuador$mco,nsampling=5000,site=NA,bird=NA,plant=NA,trait_plant=NA,mismatch=NA,barrier=c(0,1),pheno=NA,abond=NA,
type="barrier",random_effects=NA,month=NA,year=NA,duration=NA)
pre3$Country="Ecuador"
pred_barr=rbind(pre1,pre2,pre3)

pl2=ggplot(data=pred_barr,aes(x=as.factor(barrier),y=average_proba,color=Country))+
geom_pointrange(aes(ymin=lwr_proba,ymax=upr_proba),fatten=3)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
plot.title=element_text(size=14,face="bold",hjust = 0))+ylab("Probability of interaction")+scale_x_discrete(labels=c("no barrier","barrier"))+xlab("Exploitation barrier")+
labs(col="",fill="")+ggtitle("b")+scale_color_manual(values=couleurs)


##################### Panel c ##################
res_partition=fread(here("data_zenodo","sparsity_estimates.txt"))
forplot=melt(res_partition,id.vars=c("site","min_transect_elev","duration","Country","tot","pzerom"))

couleurs2=c("#477998","#291F1E","#F64740")

forplot=forplot[order(forplot$Country,forplot$min_transect_elev),]
forplot$site=factor(forplot$site,levels=unique(forplot$site))
forplot$variable=as.character(forplot$variable)
forplot$variable[forplot$variable=="barrier"]="exploitation barrier"
forplot$variable=factor(forplot$variable,levels=c("unknown","exploitation barrier"))

summary(subset(res_partition,duration %in% c(2000))$barrier)

pl3=ggplot(data=subset(forplot,duration %in% c(2000) & variable!="complementarity"),aes(y=value,fill=variable,x=site))+geom_bar(stat="identity")+
facet_grid(col=vars(Country),scales = "free_x",shrink=T,drop=TRUE,space ="free_x")+
ylab("% of invariant zeros in networks")+
theme_bw()+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background = element_blank(),
axis.text.x=element_text(color="white"),legend.position="bottom")+
labs(fill="Generated by:")+xlab("Sites")+
scale_color_manual(values=couleurs2)+scale_fill_manual(values=couleurs2)+
scale_y_continuous(labels = scales::percent)+ggtitle("c")+coord_cartesian(expand=F)

##### put together panels
leg <- ggpubr::as_ggplot(cowplot::get_legend(pl2))
leg=leg+theme(legend.margin=margin(c(0,0,0,0)))
pl2=pl2+theme(legend.position="none")

top=plot_grid(pl1,pl2,leg,align = "hv",ncol=3,rel_widths=c(1,1,0.1))
bottom=pl3
plot_grid(top,bottom,align = "hv",ncol=1,rel_heights=c(1,1.5))

# export figuure 2
pdf(paste0(here("figures"),"/Figure2.pdf"),width=7,height=6)
plot_grid(top,bottom,align = "hv",ncol=1,rel_heights=c(1,1.5))
dev.off();


##################### Possible extra panel ##################
pl4=ggplot(data=res_partition,aes(x=duration,color=Country,y=pzerom))+geom_point()+
ylab("% of zeros (non-interactions)")+
theme_bw()+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background = element_blank(),legend.position="bottom")+
labs(col="")+xlab("Duration of sampling per camera (hours)")+
scale_color_manual(values=couleurs)+scale_fill_manual(values=couleurs)+
scale_y_continuous(labels = scales::percent)+ggtitle("c")


##################### FIGURE S3 ##################
pl=ggplot(data=subset(forplot,duration %in% c(100,1000)),aes(y=value,fill=variable,x=site))+geom_bar(stat="identity")+
facet_grid(duration~Country,scales = "free_x",shrink=T,drop=TRUE,space ="free_x",
labeller = label_bquote(rows=paste("Duration of sampling: ",.(duration),"h")))+
ylab("% of invariant zero generated")+
theme_bw()+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background = element_blank(),
axis.text.x=element_text(angle=45,hjust=1))+
labs(fill="Linkage rule")+xlab("Sites")+
scale_color_manual(values=couleurs2)+scale_fill_manual(values=couleurs2)+
scale_y_continuous(labels = scales::percent)

png(paste0(here("figures"),"/Figure_S3.png"),width=1000,height=800,res=120)
pl
dev.off();

####################### FIGURE S2 ##################
#### L2 Match optimum
tabu=as.data.frame(costa_rica$data[c("hummingbird_numu","culmen_lengthu")])
tabu=cbind(tabu,costa_rica$summary[grep("match_infer",costa_rica$summary$varia),])
tabu$Country="Costa Rica"

tabu2=as.data.frame(brazil$data[c("hummingbird_numu","culmen_lengthu")])
tabu2=cbind(tabu2,brazil$summary[grep("match_infer",brazil$summary$varia),])
tabu2$Country="Brazil"

tabu3=as.data.frame(ecuador$data[c("hummingbird_numu","culmen_lengthu")])
tabu3=cbind(tabu3,ecuador$summary[grep("match_infer",ecuador$summary$varia),])
tabu3$Country="Ecuador"

pred_L2=rbind(tabu,tabu2,tabu3)
names(pred_L2)[1]="hummingbird_num"
pred_L2=merge(pred_L2,tabsp,by=c("hummingbird_num","Country"))

pred_bis=pred_L2[pred_L2$hummingbird_species %in% pred_L2$hummingbird_species[duplicated(pred_L2$hummingbird_species)],]
pred_bis=pred_bis[order(pred_bis$hummingbird_species),]
compar1=ggplot(data=pred_bis)+geom_pointrange(aes(x=hummingbird_species,y=mean,ymax=l95,ymin=u95,color=Country))+
geom_point(aes(x=hummingbird_species,y=culmen_lengthu),col="grey")+
scale_color_manual(values=couleurs[-1])+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),
plot.title=element_text(size=14,face="bold",hjust = 0),axis.text.x=element_text(angle=45,hjust=1))+xlab("Hummingbird species")+ylab("Optimal corolla length (mm)")

nrow(subset(pred_L2,l95>culmen_lengthu | u95<culmen_lengthu))/nrow(pred_L2)

# area in which values are possible
ribb=data.frame(culmen_lengthu=seq(0.7*min(pred_L2$culmen_lengthu),1.05*max(pred_L2$culmen_lengthu),length.out=100))
ribb$ymax=1.5*ribb$culmen_lengthu
ribb$ymin=0.5*ribb$culmen_lengthu

pl5=ggplot(data=pred_L2,aes(x=culmen_lengthu,y=mean,ymax=l95,ymin=u95,color=Country))+
geom_ribbon(data=ribb,aes(y=culmen_lengthu,ymin=ymin,ymax=ymax),fill="white",col=NA)+
geom_abline(intercept=0,slope=1,col="black")+geom_pointrange(fatten=3,alpha=0.6,size=0.01)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),
plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="right",panel.background = element_rect(fill = "lightgrey"))+
ylab(bquote(atop("Optimal corolla length",paste("(",italic("L2"),", in mm)"))))+
xlab("Hummingbird bill length (mm)")+
scale_color_manual(values=couleurs)+scale_fill_manual(values=couleurs)+
labs(col="",fill="")+ggtitle("a")+coord_cartesian(expand=F)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x))) +
     scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)))+annotation_logticks()


#### L1 Barrier threshold
tabu=as.data.frame(costa_rica$data[c("hummingbird_numu","culmen_lengthu")])
tabu=cbind(tabu,costa_rica$summary[grep("barrier_infer",costa_rica$summary$varia),])
tabu$Country="Costa Rica"

tabu2=as.data.frame(brazil$data[c("hummingbird_numu","culmen_lengthu")])
tabu2=cbind(tabu2,brazil$summary[grep("barrier_infer",brazil$summary$varia),])
tabu2$Country="Brazil"

tabu3=as.data.frame(ecuador$data[c("hummingbird_numu","culmen_lengthu")])
tabu3=cbind(tabu3,ecuador$summary[grep("barrier_infer",ecuador$summary$varia),])
tabu3$Country="Ecuador"

pred_L1=rbind(tabu,tabu2,tabu3)
names(pred_L1)[1]="hummingbird_num"
pred_L1=merge(pred_L1,tabsp,by=c("hummingbird_num","Country"))

pred_bis=pred_L1[pred_L1$hummingbird_species %in% pred_L1$hummingbird_species[duplicated(pred_L2$hummingbird_species)],]
pred_bis=pred_bis[order(pred_bis$hummingbird_species),]
compar2=ggplot(data=pred_bis)+geom_pointrange(aes(x=hummingbird_species,y=mean,ymax=l95,ymin=u95,color=Country))+
geom_point(aes(x=hummingbird_species,y=culmen_lengthu),col="grey")+
scale_color_manual(values=couleurs[-1])+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),
plot.title=element_text(size=14,face="bold",hjust = 0),axis.text.x=element_text(angle=45,hjust=1))+xlab("Hummingbird species")+ylab("Threshold for exploitation barrier (corolla length, in mm)")

nrow(subset(pred_L1,l95>(culmen_lengthu+0.05*culmen_lengthu)))/nrow(pred_L1)

# area in which values are possible
ribb=data.frame(culmen_lengthu=seq(0.7*min(pred_L1$culmen_lengthu),1.05*max(pred_L1$culmen_lengthu),length.out=100))
ribb$ymax=1*ribb$culmen_lengthu
ribb$ymin=2*ribb$culmen_lengthu

pl6=ggplot(data=pred_L1,aes(x=culmen_lengthu,y=mean,ymax=l95,ymin=u95,color=Country))+
geom_ribbon(data=ribb,aes(y=culmen_lengthu,ymin=ymin,ymax=ymax),fill="white",col=NA)+
geom_abline(intercept=0,slope=1,col="black")+geom_pointrange(fatten=3,alpha=0.6,size=0.01)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_rect(fill = "lightgrey"),
plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none")+
ylab(bquote(atop("Threshold for exploitation barrier",paste("(",italic("L1"),", in mm)"))))+
xlab("Hummingbird bill length (mm)")+
scale_color_manual(values=couleurs)+scale_fill_manual(values=couleurs)+
labs(col="",fill="")+ggtitle("b")+coord_cartesian(expand=F)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x))) +
     scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)))+annotation_logticks()

leg <- ggpubr::as_ggplot(cowplot::get_legend(pl5))
leg=leg+theme(legend.margin=margin(c(0,0,0,0)))
pl5=pl5+theme(legend.position="none")

pl=plot_grid(pl5,pl6,align = "hv",ncol=2,rel_widths=c(1,1))
grid.arrange(pl,leg,ncol=2,widths=c(1,0.2))
png(paste0(here("figures"),"/Figure_S2.png"),width=2300,height=1000,res=300)
grid.arrange(pl,leg,ncol=2,widths=c(1,0.2))
dev.off();


#Export table for later
names(pred_L1)[3:10]=paste0("barrier_",names(pred_L1)[3:10])
names(pred_L2)[3:10]=paste0("match_",names(pred_L2)[3:10])
names(tabsp)[1]="hummingbird_numu"
humm_table=merge(pred_L1,pred_L2,by=c("hummingbird_num","Country"))
humm_table=merge(humm_table,tabsp,by.x=c("hummingbird_num","Country"),by.y=c("hummingbird_numu","Country"))

data.table::fwrite(humm_table,paste0(here("data_zenodo"),"/humm_table.txt"))


















