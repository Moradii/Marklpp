X.all <- superimpose(BIGNONIOIDES,
                     ARNOLD,
                     INVOLUCRATA,
                     AQUIFOLIUM,
                     POPULUS)
plot(X.all,cols="red")

marks(BIGNONIOIDES) <- rep("Bignonioides",npoints(BIGNONIOIDES))
marks(ARNOLD) <- rep("Arnold",npoints(ARNOLD))
marks(INVOLUCRATA) <- rep("Involucrata",npoints(INVOLUCRATA))
marks(AQUIFOLIUM) <- rep("Aquifolium",npoints(AQUIFOLIUM))
marks(POPULUS) <- rep("Populus",npoints(POPULUS))

X <- superimpose(
  BIGNONIOIDES,
  ARNOLD,
  INVOLUCRATA,
  AQUIFOLIUM,
  POPULUS
)

library(circlize)
cc <- colorRamp2(
  c(
  quantile(marks(X.all),0.1),
  quantile(marks(X.all),0.2),
  quantile(marks(X.all),0.3),
  quantile(marks(X.all),0.5),
  quantile(marks(X.all),0.75),
  max(marks(X.all))
  ),
  colors=c(1:6)
  )


png("vantypes.png",height = 3000,width = 3400)
par(mar=c(0,0,0,0))
plot(X,main="",cols="red",
     lwd=5,leg.side = "left",
     col="black",
     cex=5,
     leg.args=list(cex.axis=10,cex=15,lwd=5,las=1))
dev.off()

png("vandiam.png",height = 3000,width = 2800)
par(mar=c(0,0,0,0))
plot(X.all,
     main="",
     leg.side = "bottom",lwd=5,
     symap=symbolmap(pch=21,cex=5,
                     lwd=5, bg=cc,
                     range=c(0,100),etch=TRUE,
                     sep=0,size=6,
                     type="continuous"),
     leg.args=list(cex.axis=5,cex=9)
)
dev.off()


png("van.png",height = 3000,width = 7000)
par(mar=c(0,0,0,0),mfrow=c(1,2))

plot(X,main="",cols="black",
     lwd=5,leg.side = "bottom",
     col="grey",
     cex=5,
     box=T,
     # symap=symbolmap(size=1),
     leg.args=list(cex.axis=6,
                   cex=35,
                   lwd=10))

plot(X.all,
     main="",
     col="grey",
     box=T,
     leg.side = "bottom",lwd=5,
     symap=symbolmap(pch=21,cex=5,
                     lwd=5, bg=cc,
                     range=c(0,100),etch=TRUE),
     leg.args=list(cex.axis=8,cex=35,lwd=10)
)
dev.off()
