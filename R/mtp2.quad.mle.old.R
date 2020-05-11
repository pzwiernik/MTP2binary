#' Numerical procedure to optimize likelihood for binary quadratic MTP2 distributions. This procedure sequentially
#' optimizes the likelihood over convex exponential families defined by a single MTP2 inequality.
#'
#' @param phat The sample distribution
#' @param tol Tolerance of the convergence criterion.
#' @return The distribution corresponding to the MLE.
#' @examples
#' phat <- array(0,c(2,2,2,2))
#' phat[1,1,1,1] <- phat[2,1,1,1] <- phat[2,2,1,1] <- phat[2,2,2,1] <- phat[2,2,2,2] <- phat[1,1,1,2] <- phat[1,1,2,2] <- phat[1,2,2,2] <- 1/8
#' mtp2.mle(phat,tol=1e-8)

mtp2.quad.mle.old <- function(phat,edges = NULL,tol = 1e-10,maxit=1000) {
  d <- length(dim(phat)) # number of variables
  # if the set of edges was not provided, all pairs will be used
  if (is.null(edges)){
    edges <- list()
    for (i in 1:(d-1)){
      for (j in (i+1):d){
        edges <- append(edges,list(c(i,j)))
      }
    }
  }
  X <- combinat::hcube(rep(2,d), rep(1,d), rep(0,d))
  mm <- compute.moments.2(phat) # matrix of 2nd order moments in 0/1 coding
  xbar <- diag(mm)
  # compute the second moments matrix in -1/1 coding
  M <- 4*mm-2*outer(rep(1,d), xbar)-2*outer( xbar,rep(1,d))+matrix(1,d,d)
  # compute the mean vector in the -1/1 coding
  xbar <- 2*diag(mm)-1
  # compute the sample covariance matrix
  S <- M-outer(xbar,xbar)
  E <- array(0,c(d,d,2,2))
  for (i in 1:(d-1)){
    for (j in (i+1):d){
      #E[i,j,2,2] <- (1+xbar[i]+xbar[j]+M[i,j])/4
      #E[i,j,2,1] <- (1+xbar[i]-xbar[j]-M[i,j])/4
      #E[i,j,1,2] <- (1-xbar[i]+xbar[j]-M[i,j])/4
      #E[i,j,1,1] <- (1-xbar[i]-xbar[j]+M[i,j])/4
      E[i,j,,] <- apply(phat,c(i,j),sum)
    }
  }
  pold <- phat
  #old.like <- log.likelihood(pold,phat)
  p <- full.independence(phat)
  J <- matrix(0,d,d)
  Shat <- diag(diag(S))
  iters <- 0
  while ((prod(Shat-S>=-tol)==0 || max(abs(J*(Shat-S)))>tol) && (iters<maxit)){
    #while (sum(abs(pold-p))>tol){
    #print(c(p))
    iters <- iters+1
    pold <- p
    for (e in edges){
      i <- e[1]
      j<- e[2]
      pij <- apply(p,c(i,j),sum)
      Q <- matrix(0,2,2)
      for (k in 1:2){
        for (l in 1:2){
          if (pij[k,l]==0) {Q[k,l] <- 0}
          if (pij[k,l]>0) {Q[k,l] <- E[i,j,k,l]/pij[k,l]}
        }
      }
      #if (i==1 && j==2) {print((Q[1,1]*Q[2,2])/(Q[1,2]*Q[2,1]))}
      #Q <- E[i,j,,]/pij
      xx <- rep(1,d)
      # check if updated J_{ij} is nonnegative. If yes, update p.
      if (p[index(xx+diag(d)[i]+diag(d)[j])]*p[index(xx)]*Q[1,1]*Q[2,2]>=p[index(xx+diag(d)[i])]*p[index(xx+diag(d)[j])]*Q[1,2]*Q[2,1]){
        for (x in 1:2^d){
          xx <- X[x,]
          p[index(xx)] <- p[index(xx)]*Q[index(xx[c(i,j)])]
        }
      } else {
        rat <- p[index(xx+diag(d)[i,]+diag(d)[j,])]*p[index(xx)]/(p[index(xx+diag(d)[i,])]*p[index(xx+diag(d)[j,])])
        rat0 <- Q[1,1]*Q[2,2]/(Q[1,2]*Q[2,1])
        t <- 0
        if (rat0!=1){t <- -log(rat)/log(rat0)}
        #print(t)
        #print(c(rat*rat0,log(rat),log(rat0),t))
        if (t!=0){
          print(t)
          for (x in 1:2^d){
            xx <- X[x,]
            p[index(xx)] <- p[index(xx)]*(Q[index(xx[c(i,j)])])^t
          }
          p <- p/sum(p)
        }
      }
      J <- compute.J(p)
      mmhat <- compute.moments.2(p) # matrix of 2nd order moments in 0/1 coding
      xbarhat <- diag(mmhat)
      # compute the second moments matrix in -1/1 coding
      Mhat <- 4*mmhat-2*outer(rep(1,d), xbarhat)-2*outer(xbarhat,rep(1,d))+matrix(1,d,d)
      # compute the mean vector in the -1/1 coding
      xbarhat <- 2*diag(mmhat)-1
      #print(round(xbarhat,2))
      Shat <- Mhat - outer(xbarhat,xbarhat)
      print(sum(phat*log(p)))
    }
  }
  #print('Number of iterations')
  #print(iters-1)
  return(p)
}

