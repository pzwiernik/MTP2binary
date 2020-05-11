#' Numerical procedure to optimize likelihood for binary MTP2 distributions. This procedure sequentially
#' optimizes the likelihood over convex exponential families defined by a single MTP2 inequality.
#'
#' @param phat The sample distribution
#' @param tol Tolerance of the convergence criterion.
#' @return The distribution corresponding to the MLE.
#' @examples
#' phat <- array(0,c(2,2,2,2))
#' phat[1,1,1,1] <- phat[2,1,1,1] <- phat[2,2,1,1] <- phat[2,2,2,1] <- phat[2,2,2,2] <- phat[1,1,1,2] <- phat[1,1,2,2] <- phat[1,2,2,2] <- 1/8
#' mtp2.mle(phat,tol=1e-8)

mtp2.mle <- function(phat,tol = 1e-8) {
  m <- length(dim(phat))
  X <- combinat::hcube(rep(2,m), rep(1,m), rep(0,m))
  pold <- full.independence(phat)
  #old.like <- log.likelihood(pold,phat)
  p <- phat
  while (sum(abs(pold-p))>tol){
  #for (iter in 1:10000){
    pold <- p
 #   old.like <- log.likelihood(pnew,phat)
    for (i in 1:2^m){
      xx <- X[i,]
      Xpaired <- list.paired(xx)
      nX <- nrow(Xpaired)
      if (is.null(nX)==FALSE){
        for (j in 1:nX){
          yy <- Xpaired[j,]
          x <- i
          y <- index(yy)
          xymax <- index(pmax(xx,yy))
          xymin <- index(pmin(xx,yy))
          if (p[x]*p[y]>p[xymax]*p[xymin]){
            lam <- (p[x]*p[y]-p[xymax]*p[xymin])/(p[x]+p[y]+p[xymax]+p[xymin])
            p[x] <- p[x] -lam
            p[y] <- p[y] -lam
            p[xymax] <- p[xymax]+lam
            p[xymin] <- p[xymin]+lam
          }
        }
      }
    }
  }
  return(p)
}
