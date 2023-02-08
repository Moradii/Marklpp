######################################## example



pats <- list()
mvec <- c(5,5,200,5,5)
for (i in 1:5) {
  X <- rpoislpp(0.02,L=chicago$domain)
  Y <- runiflpp(1,L=chicago$domain)
  cd <- crossdist(X,Y)
  pats[[i]] <- X[which(as.vector(cd)<300)]
  marks(pats[[i]]) <- rpois(npoints(pats[[i]]),mvec[i])
}


Xnew <- pats[[1]]
for (i in 2:5) {
  Xnew <- superimpose.lpp(Xnew,pats[[i]])
}

hist(marks(Xnew),breaks = 20)
plot(Xnew,cols="red")

v <- markcorr(as.ppp(Xnew),normalise = T,correction = "none",
              method = "loess",f=function(m1,m2){m1==m2})
r <- v$r*1.25
ch <- markcorr.lpp(as.ppp(Xnew),r=r,ftype = "equ",f=function(m1,m2){m1==m2})


lm <- densityQuick.lpp(Xnew,at="points")
lkspat <- linearKinhom(Xnew,lambda = as.vector(lm),
                       normalise = F,
                       leaveoneout = F)
chinhom <- linearinhommarkk.lpp(Xnew,
                                r=lkspat$r,
                                lambda = lm,
                                normalize = F)
plot(lkspat,ylim = c(0,2200))
points(chinhom$r,chinhom$markKinhom,type="l",
       col="blue")

