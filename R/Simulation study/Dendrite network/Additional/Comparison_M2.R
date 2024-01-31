rm(list = ls())

library(marklpp)
library(parallel)
library(spatstat)

########################################### M2
load("~/M2.rda")
load("~/r.rda")

###############################
############################### corr
nsim <- 199
M2_corr <- mclapply(X=1:nsim,function(i){
  markcorr.lpp(M2[[i]],
               ftype = "corr",
               normalise = TRUE,
               r=r_L,
               method = "density"
  )[,1]
},mc.cores = 8)

M2_corr <- do.call(rbind,M2_corr)

nrank <- 5

min_L_corr_M2 <- apply(M2_corr, 2, function(x) (sort(x))[nrank])
max_L_corr_M2 <- apply(M2_corr, 2, function(x) (sort(x))[nsim-nrank+1])

d_L_corr_M2 <- data.frame(r=r_L,min=min_L_corr_M2,max=max_L_corr_M2)

par(mfrow=c(1,1)) ## put two plots together

png("M2corrsup.png",height = 1200,width = 1200)
par(mar=c(12,12,1,1))
plot(d_L_corr_M2$r,d_L_corr_M2$min,
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
points(d_L_corr_M2$r,d_L_corr_M2$min,type="l",col="white")
points(d_L_corr_M2$r,d_L_corr_M2$max,type="l",col="white")
polygon(c(d_L_corr_M2$r, rev(d_L_corr_M2$r)), c(d_L_corr_M2$max, rev(d_L_corr_M2$min)),
        col = "grey70", border = NA)
points(d_L_corr_M2$r,colMeans(M2_corr),type="l",lwd=3)
abline(h=1,col="red",lty=2,lwd=3)
dev.off()

###############################
############################### vario

M2_vario <- mclapply(X=1:nsim,function(i){
  markcorr.lpp(M2[[i]],
               ftype = "vario",
               f=function(m1,m2){0.5*(m1-m2)^2},
               normalise = TRUE,
               r=r_L,
               method = "density"
  )[,1]
},mc.cores = 8)

M2_vario <- do.call(rbind,M2_vario)

nrank <- 5

min_L_vario_M2 <- apply(M2_vario, 2, function(x) (sort(x))[nrank])
max_L_vario_M2 <- apply(M2_vario, 2, function(x) (sort(x))[nsim-nrank+1])

d_L_vario_M2 <- data.frame(r=r_L,min=min_L_vario_M2,max=max_L_vario_M2)

par(mfrow=c(1,1)) ## put two plots together

png("M2variosup.png",height = 1200,width = 1200)
par(mar=c(12,12,1,1))
plot(d_L_vario_M2$r,d_L_vario_M2$min,
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
points(d_L_vario_M2$r,d_L_vario_M2$min,type="l",col="white")
points(d_L_vario_M2$r,d_L_vario_M2$max,type="l",col="white")
polygon(c(d_L_vario_M2$r, rev(d_L_vario_M2$r)), c(d_L_vario_M2$max, rev(d_L_vario_M2$min)),
        col = "grey70", border = NA)
points(d_L_vario_M2$r,colMeans(M2_vario),type="l",lwd=3)
abline(h=1,col="red",lty=2,lwd=3)
dev.off()

###############################
############################### Shimanti

M2_shimanti <- mclapply(X=1:nsim,function(i){
  markcorr.lpp(M2[[i]],
               ftype = "shimanti",
               normalise = TRUE,
               r=r_L,
               method = "density"
  )[,1]
},mc.cores = 8)

M2_shimanti <- do.call(rbind,M2_shimanti)

nrank <- 5

min_L_shimanti_M2 <- apply(M2_shimanti, 2, function(x) (sort(x))[nrank])
max_L_shimanti_M2 <- apply(M2_shimanti, 2, function(x) (sort(x))[nsim-nrank+1])

d_L_shimanti_M2 <- data.frame(r=r_L,min=min_L_shimanti_M2,max=max_L_shimanti_M2)

par(mfrow=c(1,1)) ## put two plots together

png("M2Shimantisup.png",height = 1200,width = 1200)
par(mar=c(12,12,1,1))
plot(d_L_shimanti_M2$r,d_L_shimanti_M2$min,
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
points(d_L_shimanti_M2$r,d_L_shimanti_M2$min,type="l",col="white")
points(d_L_shimanti_M2$r,d_L_shimanti_M2$max,type="l",col="white")
polygon(c(d_L_shimanti_M2$r, rev(d_L_shimanti_M2$r)), c(d_L_shimanti_M2$max, rev(d_L_shimanti_M2$min)),
        col = "grey70", border = NA)
points(d_L_shimanti_M2$r,colMeans(M2_shimanti),type="l",lwd=3)
abline(h=0,col="red",lty=2,lwd=3)
dev.off()

###############################
############################### Beisbart

M2_Bei <- mclapply(X=1:nsim,function(i){
  markcorr.lpp(M2[[i]],
               ftype = "Beisbart",
               f=function(m1,m2){m1+m2},
               normalise = TRUE,
               r=r_L,
               method = "density"
  )[,1]
},mc.cores = 8)

M2_Bei <- do.call(rbind,M2_Bei)

nrank <- 5

min_L_Bei_M2 <- apply(M2_Bei, 2, function(x) (sort(x))[nrank])
max_L_Bei_M2 <- apply(M2_Bei, 2, function(x) (sort(x))[nsim-nrank+1])

d_L_Bei_M2 <- data.frame(r=r_L,min=min_L_Bei_M2,max=max_L_Bei_M2)

par(mfrow=c(1,1)) ## put two plots together

png("M2Beisbartsup.png",height = 1200,width = 1200)
par(mar=c(12,12,1,1))
plot(d_L_Bei_M2$r,d_L_Bei_M2$min,
     type = "n",
     ylab="", xlab = "",
     ylim = c(0.5,1.7),
     cex.lab=4,
     cex.axis=4, las=3,
     main=""
)
title(ylab=expression(italic(I[mm]^paste(Bei,",",L)*(r[L]))),
      line=6,cex.lab=5)
title(xlab = expression(italic(r[L])),
      line=10,cex.lab=5)
points(d_L_Bei_M2$r,d_L_Bei_M2$min,type="l",col="white")
points(d_L_Bei_M2$r,d_L_Bei_M2$max,type="l",col="white")
polygon(c(d_L_Bei_M2$r, rev(d_L_Bei_M2$r)), c(d_L_Bei_M2$max, rev(d_L_Bei_M2$min)),
        col = "grey70", border = NA)
points(d_L_Bei_M2$r,colMeans(M2_Bei),type="l",lwd=3)
abline(h=1,col="red",lty=2,lwd=3)
dev.off()

save.image("M2_additional.rda")
