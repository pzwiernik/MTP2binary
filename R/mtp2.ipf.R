mtp2.ipf <- function(phat,tol = 1e-8) {
  m <- length(dim(phat))
  X <- combinat::hcube(rep(2,m), rep(1,m), rep(0,m))
  pold <- phat
  pnew <- full.independence(phat)
  while (sum(abs(pnew-pold))>tol){
    pold <- pnew
    for (i in 1:2^m){
      xx <- X[i,]
      change <- phat[i]-pold[i]
      #Case 1: The change is positive
      if (change >= 0){
        Xpaired <- list.paired(xx)
        nX <- nrow(Xpaired)
        aux <- 1
        if (is.null(nX)==FALSE){
          aux <- rep(0,nX)
          for (j in 1:nX){
            yy <- Xpaired[j,]
            u <- index(pmax(xx,yy))
            v <- index(pmin(xx,yy))
            aux[j] <- (pold[u]*pold[v]-pold[i]*pold[index(yy)])/(change*pold[index(yy)])
          }
        }
        pnew[i] <- pold[i] + max(aux)*change
        pnew <- pnew/sum(pnew)
      }
      if (change<0){
        Xup <- list.upcover(xx)
        nX <- nrow(Xup)
        aux <- 0
        if (is.null(nX)==FALSE){
          aux <- rep(0,nX)
          for (j in 1:nX){
            yy <- Xup[j,1:m]
            zz <- Xup[j,(m+1):(2*m)]
            u <- index(pmax(yy,zz))
            aux[j] <- -(pold[u]*pold[i]-pold[index(yy)]*pold[index(zz)])/(change*pold[u])
          }
        }
      Xdown <- list.downcover(x)
        nX <- nrow(Xdown)
        aux2 <- 0
        if (is.null(nX)==FALSE){
          aux2 <- rep(0,nX)
          for (j in 1:nX){
            yy <- Xdown[j,1:m]
            zz <- Xdown[j,(m+1):(2*m)]
            u <- index(pmin(yy,zz))
            aux2[j] <- -(pold[u]*pold[i]-pold[index(yy)]*pold[index(zz)])/(change*pold[u])
          }
        }
        pnew[i] <- pold[i] + max(cbind(aux,aux2))*change
        pnew <- pnew/sum(pnew)
      }
    }
    }
  return(pnew)
}
