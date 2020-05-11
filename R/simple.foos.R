#' Compute vector index of an array from the vector of indices. This is inverse to arrayInd.
#'
#' @param x The index vector
#' @return The number
#' @examples
#' index(c(1,2,1,2,1))
index <- function(x){
  m <- length(x)
  return(1+sum((2^(0:(m-1)))*(x-rep(1,m))))
}


#' Compute the full independence distribution with the same one-way margins
#'
#' @param p An binary array with entries of the distribution
#' @return The binary array which is an outer product of one-way margins of \code{p}.
#' @examples
#' p <- array(0,c(2,2,2))
#' p[1,1,1] <- p[1,2,1] <- p[2,1,2] <- 1/3
#' full.independence(p)
full.independence <- function(p){
  m <- length(dim(p))
  q <- apply(p,1,sum)
  for (i in 2:m) {
    q <- outer(q,apply(p,i,sum))
  }
  if (sum(q==0)>0) print('The MLE does not exist for this sample');
  return(q)
}

#' Compute the log-likelihood
#'
#' @param theta An array with parameters.
#' @param phat An array with the sample distribution.
#' @return The likelihood of \code{phat} evaluated at \code{theta}.
#' @examples
#' TBA
log.likelihood <- function(p,phat){
  m <- length(dim(phat))
  like <- 1
  for (i in 1:2^m){
    like <- like * p[i]^(phat[i])
  }
  return(log(like))
}

#' List all vectors that differ from \code{x} only in two coordinates and they are not comparable to \code{x}.
#'
#' @param x A vectors of 1s and 2s.
#' @return A matrix whose rows are the vectors that differ from \code{x} only in two coordinates and they are not comparable to \code{x}.
#' @examples
#' list.paired(c(1,2,1,2,1))
list.paired <- function(xx){
  m <- length(xx)
  A <- which(xx==1)
  B <- setdiff(1:m,A)
  basis <- diag(m)
  ninst <- length(A)*length(B)
  if (ninst==0) return(NULL);
  Xpaired <- matrix(0,ninst,m)
  count <- 1
  for (i in A){
    for (j in B){
      Xpaired[count,] <- xx+basis[i,]-basis[j,]
      count <- count +1
    }
  }
  return(Xpaired)
}

#' List all pairs of vectors that cover \code{x} and their minimum is \code{x}
#'
#' @param x A vector of 1s and 2s.
#' @return A matrix whose rows contain all pairs of vectors that cover \code{x} and their minimum is \code{x}
#' @examples
#' list.upcover(c(1,2,1,2,1))
list.upcover <- function(x){
  m <- length(x)
  A <- which(x==1)
  basis <- diag(m)
  ninst <- length(A)*(length(A)-1)/2
  if (ninst==0) return(NULL);
  Xcoverup <- matrix(0,ninst,2*m)
  count <- 1
  for (i in setdiff(A,max(A))){
    for (j in intersect(A,(i+1):m)){
      Xcoverup[count,] <- cbind(x+basis[i,],x+basis[j,])
      count <- count +1
    }
  }
  return(Xcoverup)
}

#' List all pairs of vectors covered by \code{x} whose maximum is \code{x}
#'
#' @param x A vector of 1s and 2s.
#' @return A matrix whose rows contain all pairs of vectors covered by \code{x} whose maximum is \code{x}
#' @examples
#' list.downcover(c(1,2,1,2,1))
list.downcover <- function(x){
  m <- length(x)
  A <- which(x==2)
  basis <- diag(m)
  ninst <- length(A)*(length(A)-1)/2
  if (ninst==0) return(NULL);
  Xcoverdown <- matrix(0,ninst,2*m)
  count <- 1
  for (i in setdiff(A,max(A))){
    for (j in intersect(A,(i+1):m)){
      Xcoverdown[count,] <- cbind(x-basis[i,],x-basis[j,])
      count <- count +1
    }
  }
  return(Xcoverdown)
}

#' Checks if the given array gives an MTP2 distribution.
#'
#' @param p Binary array
#' @param tol Tolerance with which the inequalities are checked.
#' @return Logical value, where \code{TRUE} means that \code{p} is MTP2.
#' @examples
#' TBA
is.mtp2 <- function(p,tol=1e-8){
  m <- length(dim(p))
  X <- combinat::hcube(rep(2,m), rep(1,m), rep(0,m))
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
        if (p[xymax]*p[xymin]-p[x]*p[y]< -tol){
          print(c(x,y))
          return(FALSE)}
    }
    }
  }
  return(TRUE)
}

#' Computes means for 0/1 coding
#'
#' @param p Binary array
#' @return Vector of means.
#' @examples
#' TBA
compute.moments.1 <- function(p){
  d <- length(dim(p))
  m <- rep(0,d)
  for (i in 1:d){
    m[i] <- apply(p,i,sum)[2]
  }
  return(m)
}

#' Computes second-order moments for 0/1 coding
#'
#' @param p Binary array
#' @return matrix of second moments.
#' @examples
#' TBA
compute.moments.2 <- function(p){
  d <- length(dim(p))
  m <- matrix(0,d,d)
  diag(m) <- compute.moments.1(p)
  for (i in 1:(d-1)){
    for (j in (i+1):d){
      m[i,j] <- m[j,i] <- apply(p,c(i,j),sum)[2,2]
    }
  }
  return(m)
}

ci.check <- function(p,tol=1e-8){
  d <- length(dim(p))
  A <- matrix(1,d,d)-diag(d)
  for (i in 1:(d-1)){
    for (j in (i+1):d){
      a <- rep(1,d)
      if (log(p[index(a+diag(d)[i,]+diag(d)[j,])]*p[index(a)]/(p[index(a+diag(d)[i,])]*p[index(a+diag(d)[j,])]))<tol) A[i,j] <- A[j,i] <- 0
    }
  }
  return(A)
}

