rm(list = ls())
library(spatstat)
library(parallel)


i=37
X <- flu$pattern[[i]]
marks(X) <- as.factor(marks(X))

M2 <- X[X$marks=="M2"]
HA <- X[X$marks=="HA"]

plot(M2,pch=20,cols="red",main="")
plot(HA,pch=20,cols="blue",main="",add = T)

bw_M2 <- bw.scott(M2)
bw_HA <- bw.scott(HA)

sig <- c(exp(mean(log(c(bw_M2[1],bw_HA[1])))),
         exp(mean(log(c(bw_M2[2],bw_HA[2])))))

marks(M2) <- density.ppp(M2,
                         sigma = sig,
                         diggle = TRUE,
                         positive = TRUE,
                         at="points")

marks(HA) <- density.ppp(HA,
                         sigma = sig,
                         diggle = TRUE,
                         positive = TRUE,
                         at="points")


##############################################
##############################################
############################################## Envelope
##############################################
##############################################

###########################       M2
m_M2 <- markcorr(M2,
                   correction = "none",
                   method = "density",
                   normalise = TRUE)
plot(m_M2)

nsim <- 199
nrank <- 5
r <- m_M2$r

M2.rand <- rlabel(M2,nsim = nsim)


mk_2d_label <- mclapply(X=1:nsim,function(i){
  markcorr(as.ppp(M2.rand[[i]]),
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

png("M2corr.png",height = 1200,width = 1200)
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
points(d_label$r,m_M2$un,type="l",lwd=3)
abline(h=1,col="red",lty=2,lwd=3)
dev.off()


###########################       HA
m_HA <- markcorr(HA,
                   correction = "none",
                   method = "density",
                   normalise = TRUE,
                   r=r)
plot(m_HA)


HA.rand <- rlabel(HA,nsim = nsim)


mk_2d_label_HA <- mclapply(X=1:nsim,function(i){
  markcorr(HA.rand[[i]],
           correction = "none",
           method = "density",
           normalise = TRUE,
           r=r
  )$un
},mc.cores = 8)

mk_2d_label_HA <- do.call(rbind,mk_2d_label_HA)

min_2d_lab_HA <- apply(mk_2d_label_HA, 2, function(x) (sort(x))[nrank])
max_2d_lab_HA <- apply(mk_2d_label_HA, 2, function(x) (sort(x))[nsim-nrank+1])

d_label_HA <- data.frame(r=r,min=min_2d_lab_HA,max=max_2d_lab_HA)

png("HAcorr.png",height = 1200,width = 1200)
par(mar=c(12,12,1,1))
plot(d_label_HA$r,d_label_HA$min,
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
points(d_label_HA$r,d_label_HA$min,type="l",col="white")
points(d_label_HA$r,d_label_HA$max,type="l",col="white")
polygon(c(d_label_HA$r, rev(d_label_HA$r)), c(d_label_HA$max, rev(d_label_HA$min)),
        col = "grey70", border = NA)
points(d_label_HA$r,m_HA$un,type="l",lwd=3)
abline(h=1,col="red",lty=2,lwd=3)
dev.off()

###########################       Cross K and J

marks(HA) <- cbind(int=marks(HA),
                           type=rep(1,npoints(HA)))
marks(M2) <- cbind(int=marks(M2),
                             type=rep(2,npoints(M2)))

X <- superimpose(HA,M2)
Y <- X
marks(Y)[marks(Y)==1] <- rep("HA",npoints(HA))
marks(Y)[marks(Y)==2] <- rep("M2",npoints(M2))

marks(Y) <- as.factor(marks(Y)[,2])

png("FLutypes.png",height = 3500,width = 3800)
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
                   i="M2",
                   j="HA",
                   lambdaI = marks(M2)[,1],
                   lambdaJ = marks(HA)[,1],
                   correction = "border")

png("Kcross.png",height = 1200,width = 1200)
par(mar=c(12,12,1,1))
plot(KK,legend=F,
     ylab="", xlab = "",
     cex.lab=4,
     cex.axis=3, las=3,
     main="",
     lwd=3
)
title(ylab=expression(italic(K["HA, M2"]^{inhom}*(r))),
      line=6,cex.lab=4)
title(xlab = expression(italic(r)),
      line=10,cex.lab=4)
dev.off()


JJ <- Jcross.inhom(Y,
                   i="M2",
                   j="HA",
                   lambdaI = marks(M2)[,1],
                   lambdaJ = marks(HA)[,1],
                   lambdamin = min(c(marks(M2)[,1],marks(HA)[,1]))*0.95
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
title(ylab=expression(italic(J["M2, HA"]^{inhom}*(r))),
      line=6,cex.lab=4)
title(xlab = expression(italic(r)),
      line=10,cex.lab=4)
dev.off()


save.image(file = "Flu.RData")
