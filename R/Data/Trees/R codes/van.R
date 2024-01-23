rm(list = ls())
setwd("~/Documents")

library(sf)
library(spatstat)
library(parallel)

trees <- read_sf("trees.shp")

trees <- st_transform(trees,3348)
trees_ppp <- maptools::as.ppp.SpatialPointsDataFrame(as(trees,"Spatial"))
trees.diam <- unmark(trees_ppp)
age <- trees_ppp$marks$date_plant

date_1<-as.Date(age)
date_2<-as.Date("2023-05-02")
age <- difftime(date_2,date_1,units="days")


marks(trees.diam) <- data.frame(genus_name=trees_ppp$marks$genus_name,
                                species_name=trees_ppp$marks$species_na,
                                cultivar_name=trees_ppp$marks$cultivar_n,
                                common_name=trees_ppp$marks$common_nam,
                                diameter=trees_ppp$marks$diameter,
                                age=as.vector(age))
  

# plot(unmark(trees.diam),pch=20,cols="red")

setwd("~/Documents")
load("Vancouver_linnet.RData")

Vancouver_linnet <- connected.linnet(Vancouver_linnet,what = "component")[[1]]

trees_lpp <- as.lpp(trees.diam,L=Vancouver_linnet)
# plot(unmark(trees_lpp),pch=20,cols="red",main="")
# marks(trees_lpp)
# head(trees_lpp$data)


species.ppps <- split(trees_lpp,f="species_name")
# genus.ppps <- split(trees_lpp,f="genus_name")
# cultivar.ppps <- split(trees_lpp,f="cultivar_name")
# common.ppps <- split(trees_lpp,f="common_name")

rm(
  trees,
  trees_lpp,
  trees_ppp,
  trees.diam,
  Vancouver_linnet,
  age,
  date_1,
  date_2
)


sort(unlist(lapply(species.ppps, npoints)))
plot(unmark(species.ppps$PROCERA),pch=20,cols=2)

#####################################################
#####################################################
##################################################### INVOLUCRATA
#####################################################
#####################################################
INVOLUCRATA <- species.ppps$INVOLUCRATA
marks(INVOLUCRATA) <- INVOLUCRATA$data$diameter

plot(Smooth(INVOLUCRATA,
            sigma=bw.scott.iso(INVOLUCRATA),
            distance="euclidean",
            positive=T,
            dimyx=512))

library(marklpp)
corrr2 <- markcorr(as.ppp(INVOLUCRATA),correction = "none",
                   method = "density")
corrln <- markcorr.lpp(INVOLUCRATA,ftype = "corr",r=corrr2$r*1.25,method = "density")

plot(corrr2)
points(corrr2$r,as.numeric(corrln),col="red",type="l")


nsim <- 199
nrank <- 5
y1.rand <- rlabel(INVOLUCRATA,nsim = nsim)
r <- corrr2$r

mk_2d_label <- mclapply(X=1:nsim,function(i){
  markcorr(as.ppp(y1.rand[[i]]),
           correction = "none",
           method = "density",
           normalise = TRUE,
           r=r
  )$un
},mc.cores = 8)

mk_2d_label <- do.call(rbind,mk_2d_label)

r_L <- r*1.25
mk_L_label <- mclapply(X=1:nsim,function(i){
  markcorr.lpp(y1.rand[[i]],
               ftype = "corr",
               normalise = TRUE,
               r=r_L,
               method = "density"
  )
},mc.cores = 8)

mk_L_label <- do.call(rbind,mk_L_label)




min_2d_lab <- apply(mk_2d_label, 2, function(x) (sort(x))[nrank])
max_2d_lab <- apply(mk_2d_label, 2, function(x) (sort(x))[nsim-nrank+1])

d_label <- data.frame(r=r,
                      min=min_2d_lab,max=max_2d_lab)






min_L_label <- apply(mk_L_label, 2, function(x) (sort(x))[nrank])
max_L_label <- apply(mk_L_label, 2, function(x) (sort(x))[nsim-nrank+1])

d_L_label <- data.frame(r=r_L,min=min_L_label,max=max_L_label)

png("INVOLUCRATA.png",height = 1200,width = 1200)
par(mar=c(12,12,1,1))
plot(d_L_label$r,d_L_label$min,
     type = "n",
     ylab="", xlab = "",
     ylim = c(0.5,2),
     cex.lab=4,
     cex.axis=4, las=3,
     main=""
)
title(ylab=expression(italic(kappa[mm]^L*(r[L]))),
      line=6,cex.lab=5)
title(xlab = expression(italic(r[L])),
      line=10,cex.lab=5)
points(d_L_label$r,d_L_label$min,type="l",col="white")
points(d_L_label$r,d_L_label$max,type="l",col="white")
polygon(c(d_L_label$r, rev(d_L_label$r)), c(d_L_label$max, rev(d_L_label$min)),
        col = "grey70", border = NA)
points(d_L_label$r,corrln,type="l",lwd=3)
abline(h=1,col="red",lty=2,lwd=3)
dev.off()

png("INVOLUCRATAR2.png",height = 1200,width = 1200)
par(mar=c(12,12,1,1))
plot(d_label$r,d_label$min,
     type = "n",
     ylab="", xlab = "",
     ylim = c(0.5,2),
     cex.lab=4,
     cex.axis=4, las=3,
     main=""
)
title(ylab=expression(italic(kappa[mm]*(r))),
      line=6,cex.lab=5)
title(xlab = expression(italic(r)),
      line=10,cex.lab=5)
points(d_label$r,d_label$min,type="l",col="white")
points(d_label$r,d_label$max,type="l",col="white")
polygon(c(d_label$r, rev(d_label$r)), c(d_label$max, rev(d_label$min)),
        col = "grey70", border = NA)
points(d_label$r,corrr2$un,type="l",lwd=3)
abline(h=1,col="red",lty=2,lwd=3)
dev.off()

INVOLUCRATA <- as.ppp(INVOLUCRATA)
save.image("INVOLUCRATA.RData")


#####################################################
#####################################################
##################################################### BIGNONIOIDES
#####################################################
#####################################################
BIGNONIOIDES <- species.ppps$BIGNONIOIDES
marks(BIGNONIOIDES) <- BIGNONIOIDES$data$diameter

BIGNONIOIDES_smooth <- Smooth(BIGNONIOIDES,
       sigma=bw.scott.iso(BIGNONIOIDES),
       distance="euclidean",
       positive=T,
       dimyx=512)
plot(BIGNONIOIDES,cols="red",
     col="blue")
plot(BIGNONIOIDES_smooth)
plot(as.ppp(unmark(BIGNONIOIDES)),add = T,
     pch=20)

library(marklpp)
corrr2 <- markcorr(as.ppp(BIGNONIOIDES),
                   correction = "none",
                   method = "density")
corrln <- markcorr.lpp(BIGNONIOIDES,
                       ftype = "corr",
                       r=corrr2$r*1.25,
                       method = "density")

plot(corrr2,ylim = c(0,3.8))
points(corrr2$r,as.numeric(corrln),col="red",type="l")

nsim <- 199
nrank <- 5
y1.rand <- rlabel(BIGNONIOIDES,nsim = nsim)
r <- corrr2$r

mk_2d_label <- mclapply(X=1:nsim,function(i){
  markcorr(as.ppp(y1.rand[[i]]),
           correction = "none",
           method = "density",
           normalise = TRUE,
           r=r
  )$un
},mc.cores = 8)

mk_2d_label <- do.call(rbind,mk_2d_label)

r_L <- r*1.25
mk_L_label <- mclapply(X=1:nsim,function(i){
  markcorr.lpp(y1.rand[[i]],
               ftype = "corr",
               normalise = TRUE,
               r=r_L,
               method = "density"
  )
},mc.cores = 6)

mk_L_label <- do.call(rbind,mk_L_label)




min_2d_lab <- apply(mk_2d_label, 2, function(x) (sort(x))[nrank])
max_2d_lab <- apply(mk_2d_label, 2, function(x) (sort(x))[nsim-nrank+1])

d_label <- data.frame(r=r,
                      min=min_2d_lab,max=max_2d_lab)



min_L_label <- apply(mk_L_label, 2, function(x) (sort(x))[nrank])
max_L_label <- apply(mk_L_label, 2, function(x) (sort(x))[nsim-nrank+1])

d_L_label <- data.frame(r=r_L,min=min_L_label,max=max_L_label)


png("BIGNONIOIDES.png",height = 1200,width = 1200)
par(mar=c(12,12,1,1))
plot(d_L_label$r,d_L_label$min,
     type = "n",
     ylab="", xlab = "",
     ylim = c(0.5,2),
     cex.lab=4,
     cex.axis=4, las=3,
     main=""
)
title(ylab=expression(italic(kappa[mm]^L*(r[L]))),
      line=6,cex.lab=5)
title(xlab = expression(italic(r[L])),
      line=10,cex.lab=5)
points(d_L_label$r,d_L_label$min,type="l",col="white")
points(d_L_label$r,d_L_label$max,type="l",col="white")
polygon(c(d_L_label$r, rev(d_L_label$r)), c(d_L_label$max, rev(d_L_label$min)),
        col = "grey70", border = NA)
points(d_L_label$r,corrln,type="l",lwd=3)
abline(h=1,col="red",lty=2,lwd=3)
dev.off()

png("BIGNONIOIDESR2.png",height = 1200,width = 1200)
par(mar=c(12,12,1,1))
plot(d_label$r,d_label$min,
     type = "n",
     ylab="", xlab = "",
     ylim = c(0.5,2),
     cex.lab=4,
     cex.axis=4, las=3,
     main=""
)
title(ylab=expression(italic(kappa[mm]*(r))),
      line=6,cex.lab=5)
title(xlab = expression(italic(r)),
      line=10,cex.lab=5)
points(d_label$r,d_label$min,type="l",col="white")
points(d_label$r,d_label$max,type="l",col="white")
polygon(c(d_label$r, rev(d_label$r)), c(d_label$max, rev(d_label$min)),
        col = "grey70", border = NA)
points(d_label$r,corrr2$un,type="l",lwd=3)
abline(h=1,col="red",lty=2,lwd=3)
dev.off()

y1.rand <- lapply(y1.rand, as.ppp) #reduce space

BIGNONIOIDES <- as.ppp(BIGNONIOIDES)
save.image("BIGNONIOIDES.RData")

#####################################################
#####################################################
##################################################### AQUIFOLIUM
#####################################################
#####################################################


AQUIFOLIUM <- species.ppps$AQUIFOLIUM
marks(AQUIFOLIUM) <- AQUIFOLIUM$data$diameter

plot(Smooth(AQUIFOLIUM,
            sigma=bw.scott.iso(AQUIFOLIUM),
            distance="euclidean",
            positive=T,
            dimyx=512))

library(marklpp)
corrr2 <- markcorr(as.ppp(AQUIFOLIUM),
                   correction = "none",
                   method = "density")
corrln <- markcorr.lpp(AQUIFOLIUM,
                       ftype = "corr",
                       r=corrr2$r*1.25,
                       method = "density")

plot(corrr2,ylim = c(0,3.8))
points(corrr2$r,as.numeric(corrln),col="red",type="l")

nsim <- 199
nrank <- 5
y1.rand <- rlabel(AQUIFOLIUM,nsim = nsim)
r <- corrr2$r

mk_2d_label <- mclapply(X=1:nsim,function(i){
  markcorr(as.ppp(y1.rand[[i]]),
           correction = "none",
           method = "density",
           normalise = TRUE,
           r=r
  )$un
},mc.cores = 8)

mk_2d_label <- do.call(rbind,mk_2d_label)

r_L <- r*1.25
mk_L_label <- mclapply(X=1:nsim,function(i){
  markcorr.lpp(y1.rand[[i]],
               ftype = "corr",
               normalise = TRUE,
               r=r_L,
               method = "density"
  )
},mc.cores = 8)

mk_L_label <- do.call(rbind,mk_L_label)




min_2d_lab <- apply(mk_2d_label, 2, function(x) (sort(x))[nrank])
max_2d_lab <- apply(mk_2d_label, 2, function(x) (sort(x))[nsim-nrank+1])

d_label <- data.frame(r=r,
                      min=min_2d_lab,max=max_2d_lab)






min_L_label <- apply(mk_L_label, 2, function(x) (sort(x))[nrank])
max_L_label <- apply(mk_L_label, 2, function(x) (sort(x))[nsim-nrank+1])

d_L_label <- data.frame(r=r_L,min=min_L_label,max=max_L_label)

png("AQUIFOLIUM.png",height = 1200,width = 1200)
par(mar=c(12,12,1,1))
plot(d_L_label$r,d_L_label$min,
     type = "n",
     ylab="", xlab = "",
     ylim = c(0.5,2),
     cex.lab=4,
     cex.axis=4, las=3,
     main=""
)
title(ylab=expression(italic(kappa[mm]^L*(r[L]))),
      line=6,cex.lab=5)
title(xlab = expression(italic(r[L])),
      line=10,cex.lab=5)
points(d_L_label$r,d_L_label$min,type="l",col="white")
points(d_L_label$r,d_L_label$max,type="l",col="white")
polygon(c(d_L_label$r, rev(d_L_label$r)), c(d_L_label$max, rev(d_L_label$min)),
        col = "grey70", border = NA)
points(d_L_label$r,corrln,type="l",lwd=3)
abline(h=1,col="red",lty=2,lwd=3)
dev.off()

png("AQUIFOLIUMR2.png",height = 1200,width = 1200)
par(mar=c(12,12,1,1))
plot(d_label$r,d_label$min,
     type = "n",
     ylab="", xlab = "",
     ylim = c(0.5,2),
     cex.lab=4,
     cex.axis=4, las=3,
     main=""
)
title(ylab=expression(italic(kappa[mm]*(r))),
      line=6,cex.lab=5)
title(xlab = expression(italic(r)),
      line=10,cex.lab=5)
points(d_label$r,d_label$min,type="l",col="white")
points(d_label$r,d_label$max,type="l",col="white")
polygon(c(d_label$r, rev(d_label$r)), c(d_label$max, rev(d_label$min)),
        col = "grey70", border = NA)
points(d_label$r,corrr2$un,type="l",lwd=3)
abline(h=1,col="red",lty=2,lwd=3)
dev.off()

#reduce spaces
rm(
  trees,
  trees_lpp,
  trees_ppp,
  trees.diam,
  Vancouver_linnet,
  age,
  date_1,
  date_2
)

y1.rand <- lapply(y1.rand, as.ppp) #reduce space
rm(Vancouver_linnet)
AQUIFOLIUM <- as.ppp(AQUIFOLIUM)
save.image("AQUIFOLIUM.RData")

#####################################################
#####################################################
##################################################### POPULUS
#####################################################
#####################################################
POPULUS <- genus.ppps$POPULUS
marks(POPULUS) <- POPULUS$data$diameter

plot(Smooth(POPULUS,
            sigma=bw.scott.iso(POPULUS),
            distance="euclidean",
            positive=T,
            dimyx=512))

library(marklpp)
corrr2 <- markcorr(as.ppp(POPULUS),
                   correction = "none",
                   method = "density")
corrln <- markcorr.lpp(POPULUS,
                       ftype = "corr",
                       r=corrr2$r*1.25,
                       method = "density")

plot(corrr2,ylim = c(0,3.8))
points(corrr2$r,as.numeric(corrln),col="red",type="l")

nsim <- 199
nrank <- 5
y1.rand <- rlabel(POPULUS,nsim = nsim)
r <- corrr2$r

mk_2d_label <- mclapply(X=1:nsim,function(i){
  markcorr(as.ppp(y1.rand[[i]]),
           correction = "none",
           method = "density",
           normalise = TRUE,
           r=r
  )$un
},mc.cores = 8)

mk_2d_label <- do.call(rbind,mk_2d_label)

r_L <- r*1.25
mk_L_label <- mclapply(X=1:nsim,function(i){
  markcorr.lpp(y1.rand[[i]],
               ftype = "corr",
               normalise = TRUE,
               r=r_L,
               method = "density"
  )
},mc.cores = 8)

mk_L_label <- do.call(rbind,mk_L_label)




min_2d_lab <- apply(mk_2d_label, 2, function(x) (sort(x))[nrank])
max_2d_lab <- apply(mk_2d_label, 2, function(x) (sort(x))[nsim-nrank+1])

d_label <- data.frame(r=r,
                      min=min_2d_lab,max=max_2d_lab)






min_L_label <- apply(mk_L_label, 2, function(x) (sort(x))[nrank])
max_L_label <- apply(mk_L_label, 2, function(x) (sort(x))[nsim-nrank+1])

d_L_label <- data.frame(r=r_L,min=min_L_label,max=max_L_label)

png("POPULUS.png",height = 1200,width = 1200)
par(mar=c(12,12,1,1))
plot(d_L_label$r,d_L_label$min,
     type = "n",
     ylab="", xlab = "",
     ylim = c(0.5,2),
     cex.lab=4,
     cex.axis=4, las=3,
     main=""
)
title(ylab=expression(italic(kappa[mm]^L*(r[L]))),
      line=6,cex.lab=5)
title(xlab = expression(italic(r[L])),
      line=10,cex.lab=5)
points(d_L_label$r,d_L_label$min,type="l",col="white")
points(d_L_label$r,d_L_label$max,type="l",col="white")
polygon(c(d_L_label$r, rev(d_L_label$r)), c(d_L_label$max, rev(d_L_label$min)),
        col = "grey70", border = NA)
points(d_L_label$r,corrln,type="l",lwd=3)
abline(h=1,col="red",lty=2,lwd=3)
dev.off()

png("POPULUSR2.png",height = 1200,width = 1200)
par(mar=c(12,12,1,1))
plot(d_label$r,d_label$min,
     type = "n",
     ylab="", xlab = "",
     ylim = c(0.5,2),
     cex.lab=4,
     cex.axis=4, las=3,
     main=""
)
title(ylab=expression(italic(kappa[mm]*(r))),
      line=6,cex.lab=5)
title(xlab = expression(italic(r)),
      line=10,cex.lab=5)
points(d_label$r,d_label$min,type="l",col="white")
points(d_label$r,d_label$max,type="l",col="white")
polygon(c(d_label$r, rev(d_label$r)), c(d_label$max, rev(d_label$min)),
        col = "grey70", border = NA)
points(d_label$r,corrr2$un,type="l",lwd=3)
abline(h=1,col="red",lty=2,lwd=3)
dev.off()

rm(
  # trees,
  # trees_lpp,
  # trees_ppp,
  # trees.diam,
  # Vancouver_linnet,
  # age,
  # date_1,
  # date_2
  genus.ppps
)

y1.rand <- lapply(y1.rand, as.ppp) #reduce space
POPULUS <- as.ppp(POPULUS)
save.image("POPULUS.RData")

#####################################################
#####################################################
##################################################### ARNOLD
#####################################################
#####################################################
ARNOLD <- cultivar.ppps$ARNOLD
marks(ARNOLD) <- ARNOLD$data$diameter

plot(Smooth(ARNOLD,
            sigma=bw.scott.iso(ARNOLD),
            distance="euclidean",
            positive=T,
            dimyx=512))

library(marklpp)
corrr2 <- markcorr(as.ppp(ARNOLD),
                   correction = "none",
                   method = "density")
corrln <- markcorr.lpp(ARNOLD,
                       ftype = "corr",
                       r=corrr2$r*1.25,
                       method = "density")

plot(corrr2,ylim = c(0,3.8))
points(corrr2$r,as.numeric(corrln),col="red",type="l")

nsim <- 199
nrank <- 5
y1.rand <- rlabel(ARNOLD,nsim = nsim)
r <- corrr2$r

mk_2d_label <- mclapply(X=1:nsim,function(i){
  markcorr(as.ppp(y1.rand[[i]]),
           correction = "none",
           method = "density",
           normalise = TRUE,
           r=r
  )$un
},mc.cores = 8)

mk_2d_label <- do.call(rbind,mk_2d_label)

r_L <- r*1.25
mk_L_label <- mclapply(X=1:nsim,function(i){
  markcorr.lpp(y1.rand[[i]],
               ftype = "corr",
               normalise = TRUE,
               r=r_L,
               method = "density"
  )
},mc.cores = 8)

mk_L_label <- do.call(rbind,mk_L_label)

min_2d_lab <- apply(mk_2d_label, 2, function(x) (sort(x))[nrank])
max_2d_lab <- apply(mk_2d_label, 2, function(x) (sort(x))[nsim-nrank+1])

d_label <- data.frame(r=r,
                      min=min_2d_lab,max=max_2d_lab)





min_L_label <- apply(mk_L_label, 2, function(x) (sort(x))[nrank])
max_L_label <- apply(mk_L_label, 2, function(x) (sort(x))[nsim-nrank+1])

d_L_label <- data.frame(r=r_L,min=min_L_label,max=max_L_label)

png("ARNOLD.png",height = 1200,width = 1200)
par(mar=c(12,12,1,1))
plot(d_L_label$r,d_L_label$min,
     type = "n",
     ylab="", xlab = "",
     ylim = c(0.5,2),
     cex.lab=4,
     cex.axis=4, las=3,
     main=""
)
title(ylab=expression(italic(kappa[mm]^L*(r[L]))),
      line=6,cex.lab=5)
title(xlab = expression(italic(r[L])),
      line=10,cex.lab=5)
points(d_L_label$r,d_L_label$min,type="l",col="white")
points(d_L_label$r,d_L_label$max,type="l",col="white")
polygon(c(d_L_label$r, rev(d_L_label$r)), c(d_L_label$max, rev(d_L_label$min)),
        col = "grey70", border = NA)
points(d_L_label$r,corrln,type="l",lwd=3)
abline(h=1,col="red",lty=2,lwd=3)
dev.off()

png("ARNOLDR2.png",height = 1200,width = 1200)
par(mar=c(12,12,1,1))
plot(d_label$r,d_label$min,
     type = "n",
     ylab="", xlab = "",
     ylim = c(0.5,2),
     cex.lab=4,
     cex.axis=4, las=3,
     main=""
)
title(ylab=expression(italic(kappa[mm]*(r))),
      line=6,cex.lab=5)
title(xlab = expression(italic(r)),
      line=10,cex.lab=5)
points(d_label$r,d_label$min,type="l",col="white")
points(d_label$r,d_label$max,type="l",col="white")
polygon(c(d_label$r, rev(d_label$r)), c(d_label$max, rev(d_label$min)),
        col = "grey70", border = NA)
points(d_label$r,corrr2$un,type="l",lwd=3)
abline(h=1,col="red",lty=2,lwd=3)
dev.off()

rm(
  # trees,
  # trees_lpp,
  # trees_ppp,
  # trees.diam,
  # Vancouver_linnet,
  # age,
  # date_1,
  # date_2
  cultivar.ppps
)

y1.rand <- lapply(y1.rand, as.ppp) #reduce space
ARNOLD <- as.ppp(ARNOLD)
save.image("ARNOLD.RData")


################################################## mean dbh
mean(marks(INVOLUCRATA))
mean(marks(ARNOLD))
mean(marks(POPULUS))
mean(marks(BIGNONIOIDES))
mean(marks(AQUIFOLIUM))

hist(marks(INVOLUCRATA))
hist(marks(ARNOLD))
hist(marks(POPULUS))
hist(marks(BIGNONIOIDES))
hist(marks(AQUIFOLIUM))
