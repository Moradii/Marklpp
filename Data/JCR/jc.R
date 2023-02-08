setwd("~/Documents/Research/Markvariogram.lpp/Data/JCR")
jcdf <- read.csv("jcpd-calls-for-service.csv",sep = ";",header = TRUE)
class(jcdf)
jcdf <- jcdf[!(is.na(jcdf$Y)|is.na(jcdf$X)),]
xy <- jcdf[,c(17,16)]

library(sp)
library(sf)
library(spatstat)
jcdf <- SpatialPointsDataFrame(coords = xy, data = jcdf,
                               proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
jcdf_sf <- st_as_sf(jcdf)
jcdf_sf_transformed <- st_transform(jcdf_sf,3857)

jcdf_ppp <- maptools::as.ppp.SpatialPointsDataFrame(as(jcdf_sf_transformed,"Spatial"))

jcdf_ppp_animal <- jcdf_ppp[jcdf_ppp$marks$Description=="Animal Complaint"]

diff <- difftime(
  strptime(jcdf_ppp_animal$marks$time.arrived, "%Y-%m-%d T %H:%M:%OS"),
  strptime(jcdf_ppp_animal$marks$time.received, "%Y-%m-%d T %H:%M:%OS"),
  units = "mins"
)

JCR_linnet1 <- connected(JCR_linnet,what="component")[[2]]
jcdf_lpp <- as.lpp(jcdf_ppp_animal$x,jcdf_ppp_animal$y,marks = as.numeric(diff),L=JCR_linnet1)
jcdf_lpp_update <- jcdf_lpp[jcdf_lpp$data$marks>0]
plot(jcdf_lpp_update,cols="red",pch=20)

jcdf_corr <- markcorr(as.ppp(jcdf_lpp_update),correction = "none",
                      method = "loess")
jcdf_corr_lpp <- markcorr.lpp(jcdf_lpp_update,ftype = "corr",r=(jcdf_corr$r*1.25))

plot(jcdf_corr,ylim=c(0.6,1.5))
points(jcdf_corr$r,jcdf_corr_lpp,type = "l",col="blue")





jcdf_ppp_Explosion <- jcdf_ppp[jcdf_ppp$marks$Description=="Fire / Explosion"]

diff <- difftime(
  strptime(jcdf_ppp_Explosion$marks$time.arrived, "%Y-%m-%d T %H:%M:%OS"),
  strptime(jcdf_ppp_Explosion$marks$time.received, "%Y-%m-%d T %H:%M:%OS"),
  units = "mins"
)

JCR_linnet1 <- connected(JCR_linnet,what="component")[[2]]
jcdf_lpp <- as.lpp(jcdf_ppp_Explosion$x,jcdf_ppp_Explosion$y,marks = as.numeric(diff),L=JCR_linnet1)
jcdf_lpp_update <- jcdf_lpp[jcdf_lpp$data$marks>0]
plot(jcdf_lpp_update,cols="red",pch=20)

jcdf_corr <- markcorr(as.ppp(jcdf_lpp_update),correction = "none",
                      method = "loess")
jcdf_corr_lpp <- markcorr.lpp(jcdf_lpp_update,ftype = "corr",r=(jcdf_corr$r*1.25))

plot(jcdf_corr,ylim=c(0.6,1.5))
points(jcdf_corr$r,jcdf_corr_lpp,type = "l",col="blue")





jcdf_ppp_Comercial <- jcdf_ppp[jcdf_ppp$marks$Description=="Rob Comercial"]

diff <- difftime(
  strptime(jcdf_ppp_Comercial$marks$time.arrived, "%Y-%m-%d T %H:%M:%OS"),
  strptime(jcdf_ppp_Comercial$marks$time.received, "%Y-%m-%d T %H:%M:%OS"),
  units = "mins"
)

JCR_linnet1 <- connected(JCR_linnet,what="component")[[2]]
jcdf_lpp <- as.lpp(jcdf_ppp_Comercial$x,jcdf_ppp_Comercial$y,marks = as.numeric(diff),L=JCR_linnet1)
jcdf_lpp_update <- jcdf_lpp[jcdf_lpp$data$marks>0]
plot(jcdf_lpp_update,cols="red",pch=20)

jcdf_corr <- markcorr(as.ppp(jcdf_lpp_update),correction = "none",
                      method = "loess")
jcdf_corr_lpp <- markcorr.lpp(jcdf_lpp_update,ftype = "corr",r=(jcdf_corr$r*1.25))

plot(jcdf_corr,ylim=c(0.6,1.5))
points(jcdf_corr$r,jcdf_corr_lpp,type = "l",col="blue")



jcdf_ppp_DPV <- jcdf_ppp[jcdf_ppp$marks$Description=="Disabled Police Vehicle"]

diff <- difftime(
  strptime(jcdf_ppp_DPV$marks$time.arrived, "%Y-%m-%d T %H:%M:%OS"),
  strptime(jcdf_ppp_DPV$marks$time.received, "%Y-%m-%d T %H:%M:%OS"),
  units = "mins"
)

JCR_linnet1 <- connected(JCR_linnet,what="component")[[2]]
jcdf_lpp <- as.lpp(jcdf_ppp_DPV$x,jcdf_ppp_DPV$y,marks = as.numeric(diff),L=JCR_linnet1)
jcdf_lpp_update <- jcdf_lpp[jcdf_lpp$data$marks>0]
plot(jcdf_lpp_update,cols="red",pch=20)

jcdf_corr <- markcorr(as.ppp(jcdf_lpp_update),correction = "none",
                      method = "loess")
jcdf_corr_lpp <- markcorr.lpp(jcdf_lpp_update,ftype = "corr",r=(jcdf_corr$r*1.25))

plot(jcdf_corr,ylim=c(0.2,1.9))
points(jcdf_corr$r,jcdf_corr_lpp,type = "l",col="blue")
