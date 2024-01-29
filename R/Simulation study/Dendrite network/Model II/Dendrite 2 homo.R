rm(list = ls())

library(marklpp)
library(parallel)
library(spatstat)
library(LinearJ)

################################################################
################################################################
################################################################ dendrite network
################################################################
################################################################
nsim <- 199
y <- runiflpp(100,L=dendrite$domain,nsim = nsim)
b <- border.lpp(dendrite)
for (i in 1:nsim) {
  marks(y[[i]]) <- apply(crossdist(y[[i]],b), 1, min)
}

library(circlize)
cc <- colorRamp2(
  c(
    quantile(marks(y[[1]]),0.1),
    quantile(marks(y[[1]]),0.2),
    quantile(marks(y[[1]]),0.3),
    quantile(marks(y[[1]]),0.5),
    quantile(marks(y[[1]]),0.75),
    max(marks(y[[1]]))
  ),
  colors=c(1:6)
)

png("model2.png",height = 1200,width = 1200)
par(mar=c(0,0,0,0))
plot(y[[1]],
     main="",
     box=T,
     leg.side = "bottom",lwd=5,
     symap=symbolmap(pch=21,cex=5,
                     lwd=5,bg=cc,
                     range=range(marks(y[[1]])),etch=TRUE,
                     sep=0,size=6,
                     type="continuous"),
     leg.args=list(cex.axis=3,cex=9)
)
dev.off()
########################################### Markcorr
###############################
###############################

r <- seq(0,200,length.out=513)

mk_2d_for_sims <- mclapply(X=1:nsim,function(i){
  markcorr(as.ppp(y[[i]]),
           correction = "none",
           method = "density",
           normalise = TRUE,
           r=r
  )$un
},mc.cores = 8)

mk_2d_for_sims <- do.call(rbind,mk_2d_for_sims)


########################################### Markcorr.lpp
###############################
###############################

r_L <-r*1.25

mk_L_for_sims <- mclapply(X=1:nsim,function(i){
  markcorr.lpp(y[[i]],
               ftype = "corr",
               normalise = TRUE,
               r=r_L,
               method = "density"
  )[,1]
},mc.cores = 8)

mk_L_for_sims <- do.call(rbind,mk_L_for_sims)


nrank <- 5

mk_2d_for_sims <- mk_2d_for_sims

min_2d <- apply(mk_2d_for_sims, 2, function(x) (sort(x))[nrank])
max_2d <- apply(mk_2d_for_sims, 2, function(x) (sort(x))[nsim-nrank+1])

d <- data.frame(r=r,
                min=min_2d,max=max_2d)



mk_L_for_sims <- mk_L_for_sims


min_L <- apply(mk_L_for_sims, 2, function(x) (sort(x))[nrank])
max_L <- apply(mk_L_for_sims, 2, function(x) (sort(x))[nsim-nrank+1])

d_L <- data.frame(r=r_L,min=min_L,max=max_L)

par(mfrow=c(2,2)) ## put two plots together

png("model2corrnet.png",height = 1200,width = 1200)
par(mar=c(12,12,1,1))
plot(d_L$r,d_L$min,
     type = "n",
     ylab="", xlab = "",
     ylim = c(0,2),
     cex.lab=4,
     cex.axis=4, las=3,
     main=""
)
title(ylab=expression(italic(kappa[mm]^L*(r[L]))),
      line=6,cex.lab=5)
title(xlab = expression(italic(r[L])),
      line=10,cex.lab=5)
points(d_L$r,d_L$min,type="l",col="white")
points(d_L$r,d_L$max,type="l",col="white")
polygon(c(d_L$r, rev(d_L$r)), c(d_L$max, rev(d_L$min)),
        col = "grey70", border = NA)
points(d_L$r,colMeans(mk_L_for_sims),type="l",lwd=3)
abline(h=1,col="red",lty=2,lwd=3)
dev.off()

png("model2corr.png",height = 1200,width = 1200)
par(mar=c(12,12,1,1))
plot(d$r,d$min,
     type = "n",
     ylim = c(0,2),
     cex.lab=4,
     cex.axis=4,
     main="", xlab="",ylab="",las=3
)
title(ylab=expression(italic(kappa[mm](r))),
      line=6,cex.lab=5)
title(xlab = expression(italic(r)),
      line=10,cex.lab=5)
points(d$r,d$min,type="l",col="white")
points(d$r,d$max,type="l",col="white")
polygon(c(d$r, rev(d$r)), c(d$max, rev(d$min)),
        col = "grey70", border = NA)
points(d$r,colMeans(mk_2d_for_sims),type="l",lwd=3)
abline(h=1,col="red",lty=2,lwd=3)
dev.off()

png("model2corrcompare.png",height = 1200,width = 1200)
par(mar=c(12,12,1,1))
plot(d$r,colMeans(mk_L_for_sims),type="l",
     ylim = c(0,2),
     xlab = "",
     ylab="",
     main="",
     cex.lab=4,
     cex.axis=4, las=3,lwd=3,xaxt="n")
points(d$r,colMeans(mk_2d_for_sims),type="l",lty=2,lwd=3)
abline(h=1,col="red",lty=2,lwd=3)
legend("topright",lty=c(1,2),lwd=c(2,2),cex=4,
       legend = c(
         expression(italic(bar(kappa)[mm]^L*(r[L]))),
         expression(italic(bar(kappa)[mm](r)))))
dev.off()

# variance comparison
plot(d$r,apply(mk_L_for_sims,2,var),type="l",
     ylim = c(0,.02),
     xlab = "r",
     ylab="",
     main = "variance")
points(d$r,apply(mk_2d_for_sims,2,var),type="l",lty=2)
legend(x=120,y=0.02,lty=c(1,2),
       legend = c("Network","2D"))



################################################################
################################################################
################################################################ rlabel
################################################################
################################################################

y1 <- y[[1]]

y1.rand <- rlabel(y1,nsim = nsim)

########################################### Markcorr
###############################
###############################


mk_2d_label <- mclapply(X=1:nsim,function(i){
  markcorr(as.ppp(y1.rand[[i]]),
           correction = "none",
           method = "density",
           normalise = TRUE,
           r=r
  )$un
},mc.cores = 8)

mk_2d_label <- do.call(rbind,mk_2d_label)


########################################### Markcorr.lpp
###############################
###############################

mk_L_label <- mclapply(X=1:nsim,function(i){
  markcorr.lpp(y1.rand[[i]],
               ftype = "corr",
               normalise = TRUE,
               r=r_L,
               method = "density"
  )[,1]
},mc.cores = 8)

mk_L_label <- do.call(rbind,mk_L_label)


mk_2d_label <- mk_2d_label

min_2d_lab <- apply(mk_2d_label, 2, function(x) (sort(x))[nrank])
max_2d_lab <- apply(mk_2d_label, 2, function(x) (sort(x))[nsim-nrank+1])

d_label <- data.frame(r=r,
                      min=min_2d_lab,max=max_2d_lab)



mk_L_label <- mk_L_label


min_L_label <- apply(mk_L_label, 2, function(x) (sort(x))[nrank])
max_L_label <- apply(mk_L_label, 2, function(x) (sort(x))[nsim-nrank+1])

d_L_label <- data.frame(r=r_L,min=min_L_label,max=max_L_label)


par(mfrow=c(1,2)) ## put two plots together

png("model2singlenet.png",height = 1200,width = 1200)
par(mar=c(12,12,1,1))
plot(d_L_label$r,d_L_label$min,
     type = "n",xlab = "",
     ylab="",
     ylim = c(0,2),
     main="",
     cex.lab=4,
     cex.axis=4, las=3
)
title(ylab=expression(italic(kappa[mm]^L*(r[L]))),
      line=6,cex.lab=5)
title(xlab = expression(italic(r[L])),
      line=10,cex.lab=5)
points(d_L_label$r,d_L_label$min,type="l",col="white")
points(d_L_label$r,d_L_label$max,type="l",col="white")
polygon(c(d_L_label$r, rev(d_L_label$r)), c(d_L_label$max, rev(d_L_label$min)),
        col = "grey70", border = NA)
points(d_L_label$r,mk_L_for_sims[1,],type="l",lwd=3)
abline(h=1,col="red",lty=2,lwd=3)
dev.off()

png("model2singleR2.png",height = 1200,width = 1200)
par(mar=c(12,12,1,1))
plot(d_label$r,d_label$min,
     type = "n",xlab = "",
     ylab="",
     ylim = c(0,2),
     main="",
     cex.lab=4,
     cex.axis=4, las=3
)
title(ylab=expression(italic(kappa[mm](r))),
      line=6,cex.lab=5)
title(xlab = expression(italic(r)),
      line=10,cex.lab=5)
points(d_label$r,d_label$min,type="l",col="white")
points(d_label$r,d_label$max,type="l",col="white")
polygon(c(d_label$r, rev(d_label$r)), c(d_label$max, rev(d_label$min)),
        col = "grey70", border = NA)
points(d_label$r,mk_2d_for_sims[1,],type="l",lwd=3)
abline(h=1,col="red",lty=2,lwd=3)
dev.off()




save.image("dendrite2homo.RData")
