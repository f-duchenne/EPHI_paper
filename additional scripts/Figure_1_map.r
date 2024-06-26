library(data.table)
library(rgbif)
library(raster)
library(Hmisc)
library(car)
library(doBy)
library(ggcorrplot)
library(lme4)
library(ape)
library(picante)
library(doBy)
library(ggtree)
library(phytools)
library(phylotools)
library(caper)
library(dplyr)
library(sf)
library(gridExtra)
library(optimx)
library(elevatr)
library(geojsonsf)
library(ggmap)
library(RColorBrewer)
library(viridis)
library(maptools)
library(ggspatial)


#Coast lines
setwd(dir="D:/land use change/continents")
shp2 <- st_read("continents-of-the-world-merged.shp")
shp2=st_crop(shp2, c(xmin=-178.2166, xmax=0, ymin=-55.90223, ymax=60))
# The zoom level, z, impacts on how long it takes to download the imagery
# z ranges from 1 to 14

#Global elevation
setwd(dir="D:/land use change/elevation")
r=raster("alwdgg.tif")
r2=crop(r,extent(shp2))
r3 <- raster::mask(r2, sf:::as_Spatial(shp2))
r4=aggregate(r3, 2, fun=mean, expand=TRUE, na.rm=TRUE)
plot(r4)

b=st_as_sf(as.data.frame(rasterToPoints(r4)),coords=c("x","y"),crs=CRS("EPSG:4326"))
b=subset(b,!is.na(BinValues))

#DATA
setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper/data")
tab1=fread("data_for_analyses_Ecuador.txt",na.string=c("",NA))
tab2=fread("data_for_analyses_Costa-Rica.txt",na.string=c("",NA))
tab3=fread("data_for_analyses_Brazil.txt",na.string=c("",NA))
sites1=unique(tab1[,c("site","midpoint_Longitude","midpoint_Latitude")])
sites2=unique(tab2[,c("site","midpoint_Longitude","midpoint_Latitude")])
sites3=unique(tab3[,c("site","midpoint_Longitude","midpoint_Latitude")])
pts1=st_as_sf(sites1,coords=c("midpoint_Longitude","midpoint_Latitude"),crs=CRS("EPSG:4326"))
stage_bbox1 =st_bbox(pts1)
pts1=st_transform(pts1,st_crs(shp2))
pts2=st_as_sf(sites2,coords=c("midpoint_Longitude","midpoint_Latitude"),crs=CRS("EPSG:4326"))
stage_bbox2 =st_bbox(pts2)
pts2=st_transform(pts2,st_crs(shp2))
pts3=st_as_sf(sites3,coords=c("midpoint_Longitude","midpoint_Latitude"),crs=CRS("EPSG:4326"))
stage_bbox3 =st_bbox(pts3)
pts3=st_transform(pts3,st_crs(shp2))

#Get elevation zoom
ex.df <- data.frame(x= c(stage_bbox2[['xmin']]-0.5, stage_bbox2[['xmax']]+0.5), 
                    y= c(stage_bbox2[['ymin']]-0.5, stage_bbox2[['ymax']]+0.5))
elev_img <- get_elev_raster(ex.df, prj ="EPSG:4326", z = 8, clip = "bbox")
bloc=st_as_sf(as.data.frame(rasterToPoints(elev_img)),coords=c("x","y"),crs="EPSG:4326")
names(bloc)[1]="elev"
bloc=subset(bloc,!is.na(elev))

##CHANGE projection
shp2=st_transform(shp2,CRS("+proj=laea +lat_0=0 +lon_0=-50 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"))
b=st_transform(b,CRS("+proj=laea +lat_0=0 +lon_0=-50 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"))
bloc=st_transform(bloc,CRS("+proj=laea +lat_0=0 +lon_0=-50 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"))
pts1=st_transform(pts1,CRS("+proj=laea +lat_0=0 +lon_0=-50 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"))
pts2=st_transform(pts2,CRS("+proj=laea +lat_0=0 +lon_0=-50 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"))
pts3=st_transform(pts3,CRS("+proj=laea +lat_0=0 +lon_0=-50 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"))

###change elevation zoom to round
bbox_radius <- max(st_distance(pts2))/1.4
bbox_centroid<- data.frame(x = c(mean(st_coordinates(pts1)[,1]),mean(st_coordinates(pts2)[,1]),mean(st_coordinates(pts3)[,1])), 
y = c(mean(st_coordinates(pts1)[,2]),mean(st_coordinates(pts2)[,2]),mean(st_coordinates(pts3)[,2]))) %>%
st_as_sf(coords = c('x','y'), crs = st_crs(shp2))

buffer <- st_buffer(bbox_centroid[2,], dist=bbox_radius)
zoom_dat <- st_intersection(bloc, buffer)


pl1=ggplot()+
geom_sf(data=b,aes(color=BinValues),size=0.3)+
geom_sf(data=shp2,color="black",fill=NA,alpha=0)+
geom_sf(data=bbox_centroid,color="black",fill="white",shape=21,size=2,stroke = 2)+
theme_bw()+
theme(axis.title=element_blank(),axis.text=element_blank(),
axis.ticks=element_blank(),plot.title=element_text(size=14,face="bold"),
panel.grid=element_line(color="darkgrey"),legend.position="left",
axis.line=element_blank(),panel.border=element_blank(),panel.background=element_rect(fill="white"))+
coord_fixed(ratio = 1)+ggtitle("a")+
coord_sf(crs = st_crs(shp2),expand=F)+
scale_color_distiller(palette = "Spectral",limits=c(-50,max(b$BinValues)),n.breaks=4)+xlim(c(-1710854,6242145))+labs(col="Elev. (m)")+
ylim(c(-3062095,9432029))


zoomed_map <- ggplot() +
geom_sf(data = zoom_dat, aes(col=elev), size = 0.35)+
geom_sf(data=pts2,color="black",alpha=0.7,size=2)+
theme_bw()+
theme(axis.title=element_blank(),axis.text=element_blank(),
axis.ticks=element_blank(),plot.title=element_text(size=14,face="bold"),
panel.grid=element_blank(),legend.position="none",
axis.line=element_blank(),panel.border=element_blank(),panel.background=element_rect(fill="white"))+
coord_fixed(ratio = 1)+ggtitle("b")+
coord_sf(crs = st_crs(shp2),expand=F)+
scale_color_distiller(palette = "Spectral",limits=c(-50,max(b$alwdgg)))+
annotation_scale(location="tr")


#map <- grid.arrange(pl1,zoomed_map,ncol=2,widths=c(2,1))

setwd(dir="C:/Users/Duchenne/Documents/EPHI_paper")
png("fig_S1.png",width=600,height=600,res=120)
grid.arrange(pl1,zoomed_map,ncol=2,widths=c(2,1))
dev.off();

pdf("fig_1map.pdf",width=7,height=7)
grid.arrange(pl1,zoomed_map,ncol=2,widths=c(2,1))
dev.off();


zoom_dat <- map_data %>% mutate(colid = factor(row_number())) %>% st_intersection(buffer)