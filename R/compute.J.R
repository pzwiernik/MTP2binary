#' Compute the J matrix for a given distribution in the Ising model
#'
#' @param p A distribution from an Ising model
#' @return The distribution corresponding to the MLE.
#' @examples
#' p <- array(1/16,c(2,2,2,2))
#' compute.J(p)

compute.J <- function(p){
  d <- length(dim(p)) # number of variables
  J <- matrix(0,d,d)
  for (i in 1:(d-1)){
    for (j in (i+1):d){
      a <- rep(1,d)
      J[i,j] <- J[j,i] <- log(p[index(a+diag(d)[i,]+diag(d)[j,])]*p[index(a)]/(p[index(a+diag(d)[i,])]*p[index(a+diag(d)[j,])]))/4
    }
  }
  return(J)
}
