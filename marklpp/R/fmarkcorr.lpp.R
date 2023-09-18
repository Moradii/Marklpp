#' Functional mark correlation function for point patterns over a linear network
#'
#' Functional mark correlation function for point patterns over a linear network
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
#' L <- spiders$domain
#' X <- runiflpp(150,L=L)
#' m <- t(replicate(150,runif(513)))
#' marks(X) <- as.data.frame(m)
#' Fcor <- fmarkcorr.linnet(X,r, ftype = "corr", method = "density" ,normalise = TRUE)
#' plot(Fcor$r,Fcor$gw,type = "l")

#' @references Eckardt, M., Mateu, J., and Moradi, M. (2023) Point processes on linear networks with function-valued marks.
#' @return a numeric vector.
#' @author Mehdi Moradi \email{m2.moradi@yahoo.com} and Matthias Eckardt


#' @import spatstat.linnet
#' @import stats
#' @export
fmarkcorr.linnet <- function(X,r,
                             normalise=FALSE,
                             f = function(m1, m2) {m1 * m2},
                             ftype=c("corr","vario","rcorr","schlather","breisgart"),
                             method=c("density","loess"),...
                             ){
  n <- npoints(X)
  d <- pairdist.lpp(unmark(X))

  if(missing(r)){
    L <- X$domain
    rmaxdefault <- 0.98 * boundingradius(L)
    W <- Window(L)
    breaks <- handle.r.b.args(r, NULL, W, rmaxdefault = rmaxdefault)
    r <- breaks$r
  }
  rmax <- max(r)

  m <- as.data.frame(marks(X))
  nf <- dim(m)[1]
  f.len <- dim(m)[2]

  mu <- apply(m, 2, FUN=function(x) mean(x))
  sigma <- apply(m, 2, var)

  mu.twice <- apply(m, 2, FUN=function(x) 2*mean(x))
  res <- data.frame(matrix(NA, nrow=f.len, ncol=f.len))


  df <- cbind(dist=as.vector(d),id.row=rep(c(1:n),each=n),id.col=rep(c(1:n),n))
  df.filter <- df[df[,1]< rmax & df[,1]>0,]
  #id1 <-  df.filter[df.filter[ ,"id.row"] == k,]

  for(h in 1:f.len){

      m1 <- m[df.filter[,2], h]
      m2 <- m[df.filter[,3], h]
            #m1 <- mrx1[, h]
            #m2 <- mrx2[, h]

      if(ftype=="schlather"){
        m1 <- m1 - mu[h]
        m2 <- m2 - mu[h]
      }

      dfvario <- data.frame(d=df.filter[,1], ff=(f(m1,m2)))


      if(method=="density"){
        Kf <- unnormdensity(dfvario$d, weights = dfvario$ff,
                            from=min(r), to=max(r), n=length(r)
                            )$y
        ## smooth estimate of kappa_1
        K1 <- unnormdensity(dfvario$d, weights=rep(1,nrow(dfvario)),
                            from=min(r), to=max(r), n=length(r)
                            )$y
        res[, h] <- Kf/K1
      }else if(method=="loess"){
        lo <- loess(ff~d,data = dfvario,
                    control=loess.control(surface="direct"))
        res[, h] <- predict(lo, newdata=data.frame(d=r))
      }else{
        stop("method should currently be either loess or density!!!")
      }

      # lo <- loess(ff~d,data = dfvario,
      # control=loess.control(surface="direct"))
      # res[, h] <- predict(lo, newdata=data.frame(d=r))
  }

      if(normalise){
           if(ftype=="corr"){
             res <- res/(mu^2)
           } else if(ftype=="vario"){
             res <- res/sigma
           }else if(ftype=="rcorr"){
             res <- res/mu
           }else if(ftype=="breisgart"){
             res <- res/mu.twice
           }else if(ftype=="schlather"){
             res <- res/sigma}
            else{
              stop("your ftype is not supported!!")
               }
      }else{
      res <- res
           }

  # delta <- list()
  # for (i in 1:ncol(res)) {
  #   delta[[i]] <- abs(diff(res[,i]))
  # }
  # delta <- do.call(rbind,delta)

  FcorInt <-  apply(res, 2, mean)
  return(list(pw=res, gw=FcorInt, r=r))
}
