library(dplyr)
library(reshape2)
library(bipartite)
library(data.table)

EPHI_version="2024-03-18"

BarbersMatrix <- function(MATRIX) {
	return(MATRIX - (cbind(rowSums(MATRIX))%*%rbind(colSums(MATRIX)))/sum(MATRIX))
}

###############################################

WEIGHTEDMODULARITY <- function(BMatrix,Matsum,redlabels,bluelabels) {
	#see equation 8
	holdsum = 0

	for (rr in 1:length(redlabels)) {
	    for (cc in 1:length(bluelabels)) {
	        kroneckerdelta = redlabels[rr] == bluelabels[cc]
	        holdsum = holdsum + BMatrix[rr,cc] * kroneckerdelta
	    }
	}
	return(holdsum/Matsum)
}

TRACE <- function(MATRIX) { return(sum(diag(MATRIX))) }


WEIGHTEDMODULARITY2 <- function(BMatrix,Matsum,redlabels,bluelabels) {
	#see equation 9
	UNIred = unique(redlabels)
	Lred = length(UNIred)
	UNIblu = unique(bluelabels)
	Lblu = length(UNIblu)
	LABELMAT1 = matrix(0,Lred,length(redlabels))
	LABELMAT2 = matrix(0,length(bluelabels),Lblu)

	for(aa in 1:length(redlabels))
		LABELMAT1[which(UNIred == redlabels[aa]),aa] = 1

	for(aa in 1:length(bluelabels))
		LABELMAT2[aa,which(UNIblu == bluelabels[aa])] = 1

	return(TRACE(LABELMAT1 %*% BMatrix %*% LABELMAT2)/Matsum)

}


setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper/data")

Mv=c()
M2v=c()
resf=NULL
for(jj in 1:3){

tab=expand.grid(c("Brazil","Costa-Rica","Ecuador"),1)
pays=tab$Var1[jj]
e=tab$Var2[jj]

pre1=fread(paste0("initial_network_per_day",pays,".txt"))
sites=fread(paste0("Site_metadata_",pays,".txt"),na.strings = c("",NA))
sites=subset(sites,habitat!="deforested")

r_vec=0:2
scenarios=c("generalists first","specialists first")
rules=c("with forbidden links","no forbidden links")
sites_vec=unique(pre1$site)

extinctions=NULL
pre1b=pre1
for(z in sites_vec){

bidon=subset(pre1b,site==z)


bidon$proba=bidon$average_proba

bidon$inter=bidon$freq*bidon$proba
bidon$compl=-abs(bidon$Tube_length-bidon$match_infer)

bidon=bidon %>% dplyr::group_by(plant_species) %>% mutate(degree=length(unique(hummingbird_species[inter>0])))

prop_forbidden=mean(1-bidon$average_proba_without_barrier)-mean(1-bidon$proba)
average_freq=mean(bidon$freq)
compl=mean(bidon$compl)

mat=dcast(bidon,plant_species~hummingbird_species,value.var="inter")
rownames(mat)=mat[,1]
mat=as.matrix(mat[,-1])
C=networklevel(mat,index="weighted connectance",weighted=T)
N=wine(mat, nreps = 1000)$wine
obj=computeModules(mat, method="Beckett")
M=obj@likelihood[1]

###################
mati=obj@moduleWeb
colnames(mati)=colnames(BMatrix)
rownames(mati)=rownames(BMatrix)

Matsum = sum(mati)
col_marginals = colSums(mati)
row_marginals = rowSums(mati)
BMatrix = BarbersMatrix(mati)

#initiliase labels
bluelabels = apply(t(obj@modules[-1,-c(1:2)]),1,function(x){which(x>0)})[(nrow(mati)+1):(nrow(mati)+ncol(mati))]
redlabels = apply(t(obj@modules[-1,-c(1:2)]),1,function(x){which(x>0)})[1:nrow(mati)]

M2=WEIGHTEDMODULARITY(mati,Matsum,redlabels,bluelabels)

################
obj=estimateBipartiteSBM(mat*1000,model = "poisson")
mati=obj$networkData
colnames(mati)=colnames(BMatrix)
rownames(mati)=rownames(BMatrix)

Matsum = sum(mati)
col_marginals = colSums(mati)
row_marginals = rowSums(mati)
BMatrix = BarbersMatrix(mati)

#initiliase labels
bluelabels = obj$memberships$col
redlabels = obj$memberships$row


M3=WEIGHTEDMODULARITY(mati,Matsum,redlabels,bluelabels)

res=data.frame(bipart=M,manual=M2,sbm=M3)


# g=graph_from_incidence_matrix(obj$networkData,weighted=T)

# #M2=prod(obj$nbBlocks)
# M2=modularity(g,c(obj$memberships$row,obj$memberships$col))

# Mv=c(Mv,M)
# M2v=c(M2v,M2)
resf=rbind(resf,res)
}
}


plot(M2v~Mv)
cor(M2v,Mv)

plot(sqrt(M2v)~Mv)
cor(sqrt(M2v),Mv)