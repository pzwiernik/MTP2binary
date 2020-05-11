#' Numerical procedure to optimize likelihood for binary quadratic MTP2 distributions. This procedure sequentially
#' optimizes the likelihood over convex exponential families defined by a single MTP2 inequality.
#'
#' @param xbar The sample mean.
#' @param M The sample second moment matrix.
#' @param edges An mx2 matrix whose columns are indices of edges.
#' @param tol Tolerance of the convergence criterion.
#' @param maxit The maximum number of iterations of the algortihm.
#' @param print Logical parameter. If FALSE, no output will be printed.
#' @return The distribution corresponding to the MLE.
#' @examples
#' phat <- array(0,c(2,2,2,2))
#' phat[1,1,1,1] <- phat[2,1,1,1] <- phat[2,2,1,1] <- phat[2,2,2,1] <- phat[2,2,2,2] <- phat[1,1,1,2] <- phat[1,1,2,2] <- phat[1,2,2,2] <- 1/8
#' mtp2.mle(phat,tol=1e-8)

mtp2.quad.mle <- function(xbar,M,edges = NULL,tol = 1e-10,maxit=1000,print=TRUE) {
  d <- length(xbar) # number of variables
  S <- M-outer(xbar,xbar)
  # if the set of edges was not provided, all pairs will be used
  if (is.null(edges)){
    edges <- c()
    for (i in 1:(d-1)){
      for (j in (i+1):d){
        if (S[i,j]>0) {edges <- rbind(edges,c(i,j))}
      }
    }
  }
  # note that we do not check what happens when edges remains empty. It is important that
  # we start from the full independence model with the same xbar
  p1 <- cbind(1-xbar,1+xbar)/2
  p <- p1[1,]
  for (i in 2:d){
    p <- outer(p,p1[i,])
  }
  # compute pairise sample distributions
  E <- array(0,c(d,d,2,2))
  for (i in 1:(d-1)){
    for (j in (i+1):d){
      E[i,j,2,2] <- (1+xbar[i]+xbar[j]+M[i,j])/4
      E[i,j,2,1] <- (1+xbar[i]-xbar[j]-M[i,j])/4
      E[i,j,1,2] <- (1-xbar[i]+xbar[j]-M[i,j])/4
      E[i,j,1,1] <- (1-xbar[i]-xbar[j]+M[i,j])/4
      #E[i,j,,] <- apply(phat,c(i,j),sum)
    }
  }
  Mold <- M
  M <- outer(xbar,xbar)
  diag(M) <- 1
  J <- matrix(0,d,d)
  Shat <- diag(diag(S))
  xbarhat <- xbar
  iters <- 0
  X <- combinat::hcube(rep(2,d), rep(1,d), rep(0,d))
  # the convergence criterion relies on the KKT conditions
  while ((prod(Shat-S>=-tol)==0 || max(abs(xbar-xbarhat))>tol ||max(abs(J*(Shat-S)))>tol) && (iters<maxit)){
    iters <- iters+1
    n.edges <- nrow(edges)
    pold <- p
    for (e in 1:n.edges){
      i <- edges[e,1]
      j<- edges[e,2]
      if (S[i,j]>0){
       pij <- apply(p,c(i,j),sum)
       Q <- E[i,j,,]/pij
     #print(c(Q))
       # in case we want to extend the support
       #Q <- matrix(0,2,2)
       #for (k in 1:2){
       #  for (l in 1:2){
        #  if (pij[k,l]==0) {Q[k,l] <- 0}
        #  if (pij[k,l]>0) {Q[k,l] <- E[i,j,k,l]/pij[k,l]}
        #}
       #}
       xx <- rep(1,d)
       # check if updated J_{ij} is nonnegative. If yes, update p.
       if (p[index(xx+diag(d)[i,]+diag(d)[j,])]*p[index(xx)]*Q[1,1]*Q[2,2]>=p[index(xx+diag(d)[i,])]*p[index(xx+diag(d)[j,])]*Q[1,2]*Q[2,1]){
         for (x in 1:2^d){
           xx <- X[x,]
           p[index(xx)] <- p[index(xx)]*Q[index(xx[c(i,j)])]
         }
       } else {
         # if updated J_{ij} would be negative, we need to find lambda* before updating p
         R <- exp(-4*J[i,j])*(pij[1,1]*pij[2,2])/(pij[1,2]*pij[2,1])
         a <- 1-R
         b <- E[i,j,1,1]+E[i,j,2,2]+R*(E[i,j,1,2]+E[i,j,2,1])
         cc <- (E[i,j,1,1]*E[i,j,2,2])-R*(E[i,j,1,2]*E[i,j,2,1])
         #print(c(a,b,cc))
         if (abs(a)<tol) {
           lam <- -cc
         } else {
           Del <- b^2-4*a*cc
           rts <- c(-b-sqrt(Del),-b+sqrt(Del))/(2*a)
           #print(c(rts,cc))
           #lam <- min(rts[which(sign(rts)==-sign(cc))])
           lam <- min(rts[which(rts>=0)])
           #print(lam)
           }
         Q <- (E[i,j,,]+matrix(c(lam,-lam,-lam,lam),2,2))/pij
         #print(lam)
         for (x in 1:2^d){
             xx <- X[x,]
             p[index(xx)] <- p[index(xx)]*(Q[index(xx[c(i,j)])])
         }
       }
       #J <- compute.J(p)
       #print(c(abs(J[1,2])<tol,abs(J[1,3])<tol,abs(J[2,3])<tol))
       #print(sum(phat*log(p)))
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
       #print(sum(phat*log(p)))
    }
    #print(round(c(p),2))
  }
  if (print==TRUE){
    if (iters==maxit) {
      message(sprintf("The maximum number of iterations %d reached. Check the KKT conditions", maxit))
    } else {
      message(cat("The algorithm converged. Number of iterations:", iters))
    }
  }
  return(p)
}

