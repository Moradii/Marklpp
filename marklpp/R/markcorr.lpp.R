#' Mark correlation function for point patterns over a linear network
#'
#' Mark correlation function for point patterns over a linear network
#'
#'
#'
#' @param X an object of class lpp
#' @param r Optional. Numeric vector. The values of the argument r at which the mark correlation function should be evaluated.
#' @param normalise If normalise=FALSE, compute only the numerator of the expression for the mark correlation.
#' @param f  Optional. Test function f used in the definition of the mark correlation function. An R function with at least two arguments. There is a sensible default.
#' @param ftype type of test function used in argument f. Currently any selection of the options "corr","vario","rcorr","schlather","equ","breisgart"
#' @param method type of smoothing, either density or loess.
#' @examples
#'  X <- rpoislpp(10,simplenet)
#'  r <- seq(0,boundingradius(simplenet),length.out=513)
#'  markcorr.lpp(X,r=r,ftype = "equ",f=function(m1,m2){m1==m2})

#' @references Eckardt, M., and Moradi, M. (2023) Marked point processes on linear networks.
#' @return a numeric vector.
#' @author Mehdi Moradi \email{m2.moradi@yahoo.com} and Matthias Eckardt


#' @import spatstat.linnet
#' @import stats
#' @export


markcorr.lpp <- function(X,r,normalise=TRUE,f = function(m1, m2) {m1 * m2},
                         ftype=c("corr","vario","rcorr","schlather","equ","breisgart"),
                         method=c("density","loess"),
                         ...){
  n <- npoints(X)
  d <- pairdist.lpp(X)
  if(missing(r)){
    L <- X$domain
    rmaxdefault <- 0.98 * boundingradius(L)
    W <- Window(L)
    breaks <- handle.r.b.args(r, NULL, W, rmaxdefault = rmaxdefault)
    r <- breaks$r
  }
  rmax <- max(r)
  m <- marks(X)

  df <- cbind(dist=as.vector(d),id.row=rep(c(1:n),each=n),id.col=rep(c(1:n),n))
  df.filter <- df[df[,1]< rmax & df[,1]>0,]
  m1 <- m[df.filter[,2]]
  m2 <- m[df.filter[,3]]

  if(ftype=="schlather"){
    m1 <- m1 - mean(m)
    m2 <- m2 - mean(m)
  }

  dfvario <- data.frame(d=df.filter[,1], ff=(f(m1,m2)))

  if(method=="density"){
    Kf <- unnormdensity(dfvario$d, weights = dfvario$ff,
                        from=min(r), to=max(r), n=length(r),
                        ...)$y
    ## smooth estimate of kappa_1
    K1 <- unnormdensity(dfvario$d, weights=rep(1,nrow(dfvario)),
                        from=min(r), to=max(r), n=length(r),
                        ...)$y
    Eff <- Kf/K1
  }else if(method=="loess"){
    lo <- loess(ff~d,data = dfvario,...)
    Eff <- predict(lo, newdata=data.frame(d=r))
  }else{
    stop("method should currently be either loess or density!!!")
  }


  if(normalise){
    if(ftype=="corr"){
      out <- Eff/(mean(m)^2)
    } else if(ftype=="vario"){
      out <- Eff/var(m)
    }else if(ftype=="rcorr"){
      out <- Eff/mean(m)
    }else if(ftype=="schlather"){
      out <- Eff/var(m)
    }else if(ftype=="equ"){
      tb <- table(m)
      out <- Eff/(sum(tb^2)/n^2)
    }else if(ftype=="breisgart"){
      out <- Eff/(2*mean(m))
    }
    else{
      stop("your ftype is not supported!!")
    }
  }else{
    out <- Eff
  }
  return(out)
}


