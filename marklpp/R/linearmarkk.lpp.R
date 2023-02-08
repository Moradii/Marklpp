#' Mark-Weighted homogeneous K Function for point patterns over a linear network
#'
#' Mark-Weighted homogeneous K Function for point patterns over a linear network
#'
#'
#'
#' @param X an object of class lpp
#' @param r Optional. Numeric vector. The values of the argument r at which the mark correlation function should be evaluated.
#'
#' @examples
#'  X <- rpoislpp(10,simplenet)
#'  r <- seq(0,boundingradius(simplenet),length.out=513)
#'  linearmarkk.lpp(X,r=r)

#' @references Eckardt, M., and Moradi, M. (2023) Marked point processes on linear networks.
#' @return a numeric vector.
#' @author Mehdi Moradi \email{m2.moradi@yahoo.com} and Matthias Eckardt


#' @import spatstat.linnet
#' @import stats
#' @export


linearmarkk.lpp <- function(X,
                            r=r){

  if (!inherits(X, "lpp")) stop("X should be from class lpp")

  l <- domain(X)
  n <- npoints(X)

  m <- marks(X)
  mm <- mean(m)

  tleng <- summary(l)$totlength
  norm <- tleng/(n*(n-1))

  sdist <- pairdist.lpp(X)

  ##
  toler <- default.linnet.tolerance(l)
  ml <- matrix(1, n, n)
  for(j in 1:n) {
    ml[ -j, j] <- countends(l, X[-j], sdist[-j,j], toler=toler)
  }

  edgetl <- ml
  maxs <- 0.7*max(sdist[!is.infinite(sdist)])

  if(imissing(r)) r <- seq(0,maxs,length.out=513)

  K <- c()
  no <- sdist == 0 | sdist==Inf

  mout <- outer(m,m)

  for (i in 1:length(r)) {
    out <- (sdist<=r[i])
    kout <- (out[!no]*mout[!no])/(edgetl[!no]*mm^2)
    K[i] <- sum(kout[!is.na(kout) & !is.infinite(kout)])
  }

  K <- K*norm

  ##
  Kout <- list(markKest=K,r=r)
  return(Kout)
}
