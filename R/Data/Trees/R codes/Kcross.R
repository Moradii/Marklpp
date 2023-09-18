
marks(X) <- as.factor(marks(X))

sig <- exp(mean(log(
 c(
   bw.scott.iso(X[X$data$marks=="Bignonioides"]),
   bw.scott.iso(X[X$data$marks=="Arnold"]),
   bw.scott.iso(X[X$data$marks=="Involucrata"]),
   bw.scott.iso(X[X$data$marks=="Aquifolium"]),
   bw.scott.iso(X[X$data$marks=="Populus"])
 ) 
)))

lamimArn <- densityQuick.lpp(X[X$data$marks=="Arnold"],dimyx=512,sigma = sig)
lamimBig <- densityQuick.lpp(X[X$data$marks=="Bignonioides"],dimyx=512,sigma = sig)
plot(lamimArn)
plot(lamimBig)
library(raster)
plot(stack(
  raster(as.im(lamimArn)),
  raster(as.im(lamimBig))
))
plot(lamimBig/lamimArn)

lam <- as.numeric(densityQuick.lpp(X,at="points",sigma = sig))


lamBig <- as.numeric(densityQuick.lpp(X[X$data$marks=="Bignonioides"],at="points",sigma = sig))
lamAr <- as.numeric(densityQuick.lpp(X[X$data$marks=="Arnold"],at="points",sigma = sig))
lamInv <- as.numeric(densityQuick.lpp(X[X$data$marks=="Involucrata"],at="points",sigma = sig))
lamAqu <- as.numeric(densityQuick.lpp(X[X$data$marks=="Aquifolium"],at="points",sigma = sig))
lamPop <- as.numeric(densityQuick.lpp(X[X$data$marks=="Populus"],at="points",sigma = sig))

Kbig <- Knetinhom(X[X$data$marks=="Bignonioides"],lambda = lamBig)
KAr <- Knetinhom(X[X$data$marks=="Arnold"],lambda = lamAr)

K <- Knetinhom(X,lambda = lam)
v <- volume(X$domain)

KBigAr <- (((npoints(X)/v)^2)*K$est - ((npoints(X[X$data$marks=="Bignonioides"])/v)^2)*Kbig$est
- ((npoints(X[X$data$marks=="Arnold"])/v)^2)*KAr$est)/ (
  2*(npoints(X[X$data$marks=="Bignonioides"])/v)*(npoints(X[X$data$marks=="Arnold"])/v)
  )
plot(K$r,KBigAr,type="l")
points(K$r,K$theo,type="l",col="red")
plot(Kbig)
plot(KAr)

K2 <- Kcross.inhom(
  as.ppp(X),
  i= "Bignonioides",
  j= "Arnold",
  lambdaI = lamBig,
  lambdaJ = lamAr
)
plot(K2)
