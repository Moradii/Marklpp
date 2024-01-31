rm(list = ls())

library(marklpp)
library(parallel)
library(spatstat)

################################################################
################################################################
################################################################ load data
################################################################
################################################################

load("~/M1.rda")
load("~/r.rda")

########################################### M1

###############################
############################### corr
nsim <- 199
m1_corr <- mclapply(X=1:nsim,function(i){
  markcorr.lpp(M1[[i]],
               ftype = "corr",
               normalise = TRUE,
               r=r_L,
               method = "density"
  )[,1]
},mc.cores = 8)

m1_corr <- do.call(rbind,m1_corr)

nrank <- 5

min_L_corr <- apply(m1_corr, 2, function(x) (sort(x))[nrank])
max_L_corr <- apply(m1_corr, 2, function(x) (sort(x))[nsim-nrank+1])

d_L_corr <- data.frame(r=r_L,min=min_L_corr,max=max_L_corr)

par(mfrow=c(1,1)) ## put two plots together

png("M1corrsup.png",height = 1200,width = 1200)
par(mar=c(12,12,1,1))
plot(d_L_corr$r,d_L_corr$min,
     type = "n",
     ylab="", xlab = "",
     ylim = c(0.3,3),
     cex.lab=4,
     cex.axis=4, las=3,
     main=""
)
title(ylab=expression(italic(kappa[mm]^L*(r[L]))),
      line=6,cex.lab=5)
title(xlab = expression(italic(r[L])),
      line=10,cex.lab=5)
points(d_L_corr$r,d_L_corr$min,type="l",col="white")
points(d_L_corr$r,d_L_corr$max,type="l",col="white")
polygon(c(d_L_corr$r, rev(d_L_corr$r)), c(d_L_corr$max, rev(d_L_corr$min)),
        col = "grey70", border = NA)
points(d_L_corr$r,colMeans(m1_corr),type="l",lwd=3)
abline(h=1,col="red",lty=2,lwd=3)
dev.off()

###############################
############################### vario

m1_vario <- mclapply(X=1:nsim,function(i){
  markcorr.lpp(M1[[i]],
               ftype = "vario",
               f=function(m1,m2){0.5*(m1-m2)^2},
               normalise = TRUE,
               r=r_L,
               method = "density"
  )[,1]
},mc.cores = 8)

m1_vario <- do.call(rbind,m1_vario)

nrank <- 5

min_L_vario <- apply(m1_vario, 2, function(x) (sort(x))[nrank])
max_L_vario <- apply(m1_vario, 2, function(x) (sort(x))[nsim-nrank+1])

d_L_vario <- data.frame(r=r_L,min=min_L_vario,max=max_L_vario)

par(mfrow=c(1,1)) ## put two plots together

png("M1variosup.png",height = 1200,width = 1200)
par(mar=c(12,12,1,1))
plot(d_L_vario$r,d_L_vario$min,
     type = "n",
     ylab="", xlab = "",
     ylim = c(0,3),
     cex.lab=4,
     cex.axis=4, las=3,
     main=""
)
title(ylab=expression(italic(gamma[mm]^L*(r[L]))),
      line=6,cex.lab=5)
title(xlab = expression(italic(r[L])),
      line=10,cex.lab=5)
points(d_L_vario$r,d_L_vario$min,type="l",col="white")
points(d_L_vario$r,d_L_vario$max,type="l",col="white")
polygon(c(d_L_vario$r, rev(d_L_vario$r)), c(d_L_vario$max, rev(d_L_vario$min)),
        col = "grey70", border = NA)
points(d_L_vario$r,colMeans(m1_vario),type="l",lwd=3)
abline(h=1,col="red",lty=2,lwd=3)
dev.off()

###############################
############################### Shimanti

m1_shimanti <- mclapply(X=1:nsim,function(i){
  markcorr.lpp(M1[[i]],
               ftype = "shimanti",
               normalise = TRUE,
               r=r_L,
               method = "density"
  )[,1]
},mc.cores = 8)

m1_shimanti <- do.call(rbind,m1_shimanti)

nrank <- 5

min_L_shimanti <- apply(m1_shimanti, 2, function(x) (sort(x))[nrank])
max_L_shimanti <- apply(m1_shimanti, 2, function(x) (sort(x))[nsim-nrank+1])

d_L_shimanti <- data.frame(r=r_L,min=min_L_shimanti,max=max_L_shimanti)

par(mfrow=c(1,1)) ## put two plots together

png("M1Shimantisup.png",height = 1200,width = 1200)
par(mar=c(12,12,1,1))
plot(d_L_shimanti$r,d_L_shimanti$min,
     type = "n",
     ylab="", xlab = "",
     ylim = c(-1,1.5),
     cex.lab=4,
     cex.axis=4, las=3,
     main=""
)
title(ylab=expression(italic(I[mm]^paste(Shi,",",L)*(r[L]))),
      line=6,cex.lab=5)
title(xlab = expression(italic(r[L])),
      line=10,cex.lab=5)
points(d_L_shimanti$r,d_L_shimanti$min,type="l",col="white")
points(d_L_shimanti$r,d_L_shimanti$max,type="l",col="white")
polygon(c(d_L_shimanti$r, rev(d_L_shimanti$r)), c(d_L_shimanti$max, rev(d_L_shimanti$min)),
        col = "grey70", border = NA)
points(d_L_shimanti$r,colMeans(m1_shimanti),type="l",lwd=3)
abline(h=0,col="red",lty=2,lwd=3)
dev.off()

###############################
############################### Beisbart

m1_Bei <- mclapply(X=1:nsim,function(i){
  markcorr.lpp(M1[[i]],
               ftype = "Beisbart",
               f=function(m1,m2){m1+m2},
               normalise = TRUE,
               r=r_L,
               method = "density"
  )[,1]
},mc.cores = 8)

m1_Bei <- do.call(rbind,m1_Bei)

nrank <- 5

min_L_Bei <- apply(m1_Bei, 2, function(x) (sort(x))[nrank])
max_L_Bei <- apply(m1_Bei, 2, function(x) (sort(x))[nsim-nrank+1])

d_L_Bei <- data.frame(r=r_L,min=min_L_Bei,max=max_L_Bei)

par(mfrow=c(1,1)) ## put two plots together

png("M1Beisbartsup.png",height = 1200,width = 1200)
par(mar=c(12,12,1,1))
plot(d_L_Bei$r,d_L_Bei$min,
     type = "n",
     ylab="", xlab = "",
     ylim = c(0.5,1.7),
     cex.lab=4,
     cex.axis=4, las=3,
     main=""
)
title(ylab=expression(italic(kappa[mm]^paste(Bei,",",L)*(r[L]))),
      line=6,cex.lab=5)
title(xlab = expression(italic(r[L])),
      line=10,cex.lab=5)
points(d_L_Bei$r,d_L_Bei$min,type="l",col="white")
points(d_L_Bei$r,d_L_Bei$max,type="l",col="white")
polygon(c(d_L_Bei$r, rev(d_L_Bei$r)), c(d_L_Bei$max, rev(d_L_Bei$min)),
        col = "grey70", border = NA)
points(d_L_Bei$r,colMeans(m1_Bei),type="l",lwd=3)
abline(h=1,col="red",lty=2,lwd=3)
dev.off()


save.image("M1_additional.rda")
