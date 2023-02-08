#' Mark-Weighted inhomogeneous K Function for point patterns over a linear network
#'
#' Mark-Weighted inhomogeneous K Function for point patterns over a linear network
#'
#'
#'
#' @param X an object of class lpp
#' @param r Optional. Numeric vector. The values of the argument r at which the mark correlation function should be evaluated.
#' @param lambda Intensity values at data points.
#' @param normalize Logical.
#'
#' @examples
#'  X <- rpoislpp(10,simplenet)
#'  r <- seq(0,boundingradius(simplenet),length.out=513)
#'  dx <- densityQuick.lpp(X,at = "points")
#'  linearinhommarkk.lpp(X,r=r,lambda=dx)

#' @references Eckardt, M., and Moradi, M. (2023) Marked point processes on linear networks.
#' @return a numeric vector.
#' @author Mehdi Moradi \email{m2.moradi@yahoo.com} and Matthias Eckardt


#' @import spatstat.linnet
#' @import stats
#' @export

linearinhommarkk.lpp <- function(X,
                                 r=r,
                                 lambda=lambda,
                                 normalize=FALSE,
                                 ...){

  if (!inherits(X, "lpp")) stop("X should be from class lpp")


  l <- domain(X)
  tleng <- summary(l)$totlength
  n <- npoints(X)

  m <- marks(X)
  mm <- mean(m)
  mout <- outer(m,m)

  sdist <- pairdist.lpp(X)

  toler <- default.linnet.tolerance(l)
  ml <- matrix(1, n, n)
  for(j in 1:n) {
    ml[ -j, j] <- countends(l, X[-j], sdist[-j,j], toler=toler)
  }


  lamden <- outer(lambda,lambda,FUN = "*")
  diag(lamden) <- 1
  edgetl <- ml*lamden

  maxs <- 0.7*max(sdist[!is.infinite(sdist)])

  if(missing(r)) r <- seq(0,maxs,length.out=513)

  K <- c()

  for (i in 1:length(r)) {
      out <- (sdist<=r[i])
      diag(out) <- 0
      kout <- (out*mout)/(edgetl*mm^2)
      K[i] <- sum(kout[!is.na(kout) & !is.infinite(kout)])

  }

  if(normalize){
    revrho <- outer(1/lambda,1/lambda,FUN = "*")
    appx <- (tleng)/(sum(revrho[lower.tri(revrho, diag = FALSE)])*2)
    K <- K*appx
  }
  else{
    K <- K/(tleng)
  }

  ##
  Kout <- list(markKinhom=K,r=r)
  return(Kout)
}
