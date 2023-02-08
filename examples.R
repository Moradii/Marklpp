library(marklpp)

###########################
###########################
#      linearinhommarkk.lpp
###########################
###########################

########################### Random field
# Random field  is the original mark, but is replaced by estimated intensity

# load JCR_linnet
L <- JCR_linnet
L <- connected(L,what = c( "components"))[[2]]

m <- as.im(function(x, y){5 - 1.5 * (x - 0.5)^2 + 2 * (y - 0.5)^2},
           W=L$window)
X <- rpoislpp(lambda = function(x,y){(0.005*x+20*y)/1000000000000},L)
plot(X,pch=20,cols="red")
marks(X) <- m[X]

plot(X,pch=20,cols="red")

Xnew <- X
lm <- densityQuick.lpp(Xnew,at="points")
marks(Xnew) <- as.vector(lm) # estimated intensity is given as marks

lmim <- densityQuick.lpp(Xnew) # just for visualising intensity map
plot(lmim)

lkspat <- linearKinhom(Xnew,lambda = as.vector(lm),
                       normalise = F,
                       leaveoneout = F,
                       r=seq(0,1500,length.out=513)) # default linearinhomK

chinhom <- linearinhommarkk.lpp(Xnew,
                                r=lkspat$r,
                                lambda = lm,
                                normalize = F) # our NEW K
plot(lkspat,ylim = c(0,2200))
points(chinhom$r,chinhom$markKinhom,type="l",
       col="blue")



########################### intensity field

X <- rpoislpp(lambda = function(x,y){
  ifelse(y<700,0.02,0.005)
},L)

X <- rpoislpp(lambda = function(x,y){(x+5*y)/400000},L)
plot(X,pch=20,cols="red")

Xnew <- X
lm <- densityQuick.lpp(Xnew,at="points")
# lm[Xnew$data$y>600] <- lm[Xnew$data$y>600]
marks(Xnew) <- as.vector(1/lm)

lmim <- densityQuick.lpp(Xnew)
plot(lmim)

lkspat <- linearKinhom(Xnew,lambda = as.vector(lm),
                       normalise = F,
                       leaveoneout = F,
                       r=seq(0,1500,length.out=513))

chinhom <- linearinhommarkk.lpp(Xnew,
                                r=lkspat$r,
                                lambda = lm,
                                normalize = F)

plot(lkspat,ylim = c(0,2200))
points(chinhom$r,chinhom$markKinhom,type="l",
       col="blue")

######################## markcorr.lpp

mcorspat <- markcorr(as.ppp(Xnew),
                     correction = "none",
                      method = "loess")
plot(mcorspat,main="")

mcor <- markcorr.lpp(Xnew,
                     r=(mcorspat$r*1.25),
                     ftype = "corr",
                     f=function(m1,m2){m1 * m2})
plot(mcorspat$r,mcor,type="l",col="blue")
abline(h=1,col="red")

plot(mcorspat,main="",ylim = c(0.9,1.2))
points(mcorspat$r,mcor,type="l",col="blue")

set.seed(1234)
Y <- rlabel(Xnew,nsim = 39)
library(parallel)
env <- mclapply(X=1:39,function(i){
  markcorr.lpp(Y[[i]],r=(mcorspat$r*1.25),ftype = "corr")
},mc.cores = 8)

env.df <- do.call(rbind,env)
minv <- apply(env.df, 2, min)
maxv <- apply(env.df, 2, max)

plot((mcorspat$r*1.25),maxv,type="l",ylim=c(0.5,1.5))
points((mcorspat$r*1.25),minv,type="l")
points((mcorspat$r*1.25),mcor,type="l",col="red")


###################################################
mcor <- markcorr.lpp(Xnew,r=(mcorspat$r*1.25),
                     ftype = "vario",
                     f=function(m1, m2) {0.5*(m1 - m2)^2})
plot((mcorspat$r*1.25),mcor,type="l")
abline(h=1,col="red")
