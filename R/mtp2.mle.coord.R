#' Numerical procedure to optimize likelihood for binary MTP2 distributions. This procedure sequentially
#' optimizes the likelihood over convex exponential families defined by a single MTP2 inequality.
#'
#' @param phat The sample distribution
#' @param tol Tolerance of the convergence criterion.
#' @return The distribution corresponding to the MLE.
#' @examples
#' phat <- array(0,c(2,2,2,2))
#' phat[1,1,1,1] <- phat[2,1,1,1] <- phat[2,2,1,1] <- phat[2,2,2,1] <- phat[2,2,2,2] <- phat[1,1,1,2] <- phat[1,1,2,2] <- phat[1,2,2,2] <- 1/8
#' mtp2.mle.coord(phat,tol=1e-8)

mtp2.mle.coord <- function(phat,tol = 1e-8) {
  m <- length(dim(phat))
  # the rows of X are all vectors xx
  X <- combinat::hcube(rep(2,m), rep(1,m), rep(0,m))
  pold <- phat
  p <- full.independence(phat) # this is the starting point
  while (sum(abs(pold-p))>tol){
    pold <- p
    # cycle through all vectors xx
    atoms <- 1+2^(0:(m-1))
    for (i in setdiff(2:(2^m-1),atoms)){
      xx <- X[i,]
      # for each xx find all yy such that they form an elementary pair (could be improved)
      Xpaired <- list.paired(xx)
      nX <- nrow(Xpaired)
      if (is.null(nX)==FALSE){
        # cycle through all pairs xx, yy
        for (j in 1:nX){
          yy <- Xpaired[j,]
          x <- i
          y <- index(yy)
          xymax <- index(pmax(xx,yy))
          xymin <- index(pmin(xx,yy))
          sT <- phat[x]+phat[y]+phat[xymax]+phat[xymin]
          sS <- p[x]+p[y]+p[xymax]+p[xymin]
          # this is the mu which adjust the four entries
          # now we compute the scaling factor according to the two cases
            s <- (1-sT)/(1-sS)
            mu <- max(c((phat[x]*phat[y]-phat[xymax]*phat[xymin])/sT,0))
            p <- s*p
            p[x] <- phat[x] -mu
            p[y] <- phat[y] -mu
            p[xymax] <- phat[xymax]+mu
            p[xymin] <- phat[xymin]+mu
         print(sum(p))
        }
      }
    }
  }
  return(p)
}
