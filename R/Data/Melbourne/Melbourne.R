rm(list = ls())
setwd("~/Documents/Research/Journal Papers/Marked point processes on linear networks.lpp/Data/Melbourne/butterfly-biodiversity-survey-2017")
library(sf)
library(spatstat)
library(parallel)
butterfly <- read_sf("butterflypark.shp")

butterfly <- st_transform(butterfly,3857)

butterfly_ppp <- maptools::as.ppp.SpatialPointsDataFrame(as(butterfly,"Spatial"))

# setwd("~/Documents/Research/Journal Papers/Marked point processes on linear networks.lpp/Data/Melbourne/municipal-boundary")
# 
# W <- read_sf("municipal-boundary.shp")
# 
# W <- st_transform(W,crs = 3857)
# W_win <- as.owin(as_Spatial(W))

W_win <- convexhull(butterfly_ppp)

butterfly_ppp <- ppp(butterfly_ppp$x,
                     butterfly_ppp$y,
                     marks = marks(butterfly_ppp)[,20:40],
                     window = W_win)

butterfly_ppp <- as.ppp(butterfly_ppp)
plot(unmark(butterfly_ppp),pch=20,cols="red",main="")

mm <- apply(apply(marks(butterfly_ppp), 2, as.factor), 2, as.numeric)

btpp <- list()
for (i in 1:21) {
btpp[[i]] <- unmark(butterfly_ppp[as.logical(mm[,i])])  
}

colnames(mm)


cabbagewhite <- unmark(butterfly_ppp[as.logical(mm[,"blue"]==0 & mm[,"prap"]==1)])
Littleblue <- unmark(butterfly_ppp[as.logical(mm[,"blue"]==1 & mm[,"prap"]==0)])

plot(cabbagewhite,pch=20,cols="red",main="")
plot(Littleblue,pch=20,cols="blue",main="",add = T)

bw_cabb <- bw.scott(cabbagewhite)
bw_blue <- bw.scott(Littleblue)

sig <- c(exp(mean(log(c(bw_cabb[1],bw_blue[1])))),
         exp(mean(log(c(bw_cabb[2],bw_blue[2])))))

marks(cabbagewhite) <- density.ppp(cabbagewhite,
                      sigma = sig,
                      diggle = TRUE,
                      positive = TRUE,
                      at="points")

marks(Littleblue) <- density.ppp(Littleblue,
                      sigma = sig,
                      diggle = TRUE,
                      positive = TRUE,
                      at="points")


##############################################
##############################################
############################################## Envelope
##############################################
##############################################

###########################       cabbagewhite
m_cabb <- markcorr(cabbagewhite,
                   correction = "none",
                   method = "density",
                   normalise = TRUE)
plot(m_cabb)

nsim <- 199
nrank <- 5
r <- m_cabb$r

cabbagewhite.rand <- rlabel(cabbagewhite,nsim = nsim)


mk_2d_label <- mclapply(X=1:nsim,function(i){
  markcorr(as.ppp(cabbagewhite.rand[[i]]),
           correction = "none",
           method = "density",
           normalise = TRUE,
           r=r
  )$un
},mc.cores = 8)

mk_2d_label <- do.call(rbind,mk_2d_label)

min_2d_lab <- apply(mk_2d_label, 2, function(x) (sort(x))[nrank])
max_2d_lab <- apply(mk_2d_label, 2, function(x) (sort(x))[nsim-nrank+1])

d_label <- data.frame(r=r,min=min_2d_lab,max=max_2d_lab)

png("cabbagewhitecorr.png",height = 1200,width = 1200)
par(mar=c(12,12,1,1))
plot(d_label$r,d_label$min,
     type = "n",
     ylab="", xlab = "",
     ylim = c(0,3),
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
points(d_label$r,m_cabb$un,type="l",lwd=3)
abline(h=1,col="red",lty=2,lwd=3)
dev.off()


###########################       Littleblue
m_blue <- markcorr(Littleblue,
                   correction = "none",
                   method = "density",
                   normalise = TRUE,
                   r=r)
plot(m_blue)


Littleblue.rand <- rlabel(Littleblue,nsim = nsim)


mk_2d_label_blue <- mclapply(X=1:nsim,function(i){
  markcorr(Littleblue.rand[[i]],
           correction = "none",
           method = "density",
           normalise = TRUE,
           r=r
  )$un
},mc.cores = 8)

mk_2d_label_blue <- do.call(rbind,mk_2d_label_blue)

min_2d_lab_blue <- apply(mk_2d_label_blue, 2, function(x) (sort(x))[nrank])
max_2d_lab_blue <- apply(mk_2d_label_blue, 2, function(x) (sort(x))[nsim-nrank+1])

d_label_blue <- data.frame(r=r,min=min_2d_lab_blue,max=max_2d_lab_blue)

png("Littlebluecorr.png",height = 1200,width = 1200)
par(mar=c(12,12,1,1))
plot(d_label_blue$r,d_label_blue$min,
     type = "n",
     ylab="", xlab = "",
     ylim = c(0,3),
     cex.lab=4,
     cex.axis=4, las=3,
     main=""
)
title(ylab=expression(italic(kappa[mm]*(r))),
      line=6,cex.lab=5)
title(xlab = expression(italic(r)),
      line=10,cex.lab=5)
points(d_label_blue$r,d_label_blue$min,type="l",col="white")
points(d_label_blue$r,d_label_blue$max,type="l",col="white")
polygon(c(d_label_blue$r, rev(d_label_blue$r)), c(d_label_blue$max, rev(d_label_blue$min)),
        col = "grey70", border = NA)
points(d_label_blue$r,m_blue$un,type="l",lwd=3)
abline(h=1,col="red",lty=2,lwd=3)
dev.off()

###########################       Cross K and J

marks(Littleblue) <- cbind(int=marks(Littleblue),
                           type=rep(1,npoints(Littleblue)))
marks(cabbagewhite) <- cbind(int=marks(cabbagewhite),
                           type=rep(2,npoints(cabbagewhite)))

X <- superimpose(Littleblue,cabbagewhite)
Y <- X
marks(Y)[marks(Y)==1] <- rep("Little blue",npoints(Littleblue))
marks(Y)[marks(Y)==2] <- rep("Cabbage white",npoints(cabbagewhite))

marks(Y) <- as.factor(marks(Y)[,2])

png("Meltypes.png",height = 3500,width = 3800)
par(mar=c(0,0,0,0))
plot(Y,main="",cols="black",
     lwd=5,leg.side = "left",
     cex=9,
     box=F,
     # symap=symbolmap(size=1),
     leg.args=list(cex.axis=6,
                   cex=9,
                   lwd=10)
)
dev.off()

KK <- Kcross.inhom(Y,
                   i="Little blue",
                   j="Cabbage white",
             lambdaJ = marks(cabbagewhite)[,1],
             lambdaI = marks(Littleblue)[,1],
             correction = "iso")

png("Kcross.png",height = 1200,width = 1200)
par(mar=c(12,12,1,1))
plot(KK,legend=F,
     ylab="", xlab = "",
     cex.lab=4,
     cex.axis=3, las=3,
     main="",
     lwd=3
)
title(ylab=expression(italic(K["Little blue, Cabbage white"]^{inhom}*(r))),
      line=6,cex.lab=4)
title(xlab = expression(italic(r)),
      line=10,cex.lab=4)
dev.off()


JJ <- Jcross.inhom(Y,
                   i="Little blue",
                   j="Cabbage white",
                   lambdaJ = marks(cabbagewhite)[,1],
                   lambdaI = marks(Littleblue)[,1],
                   lambdamin = 2.9e-06
                   )

png("Jcross.png",height = 1200,width = 1200)
par(mar=c(12,12,1,1))
plot(JJ,legend=F,
     ylab="", xlab = "",
     cex.lab=4,
     cex.axis=3, las=3,
     main="",
     lwd=3
)
title(ylab=expression(italic(J["Cabbage white, Little blue"]^{inhom}*(r))),
      line=6,cex.lab=4)
title(xlab = expression(italic(r)),
      line=10,cex.lab=4)
dev.off()


save.image(file = "Melbourne.RData")
plot(localKinhom(Y))
