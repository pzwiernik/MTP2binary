#' Compute the S matrix for a given bianry distribution in -1/1 coding
#'
#' @param p A distribution
#' @return S corresponding covariance matrix
#' @examples
#' p <- array(1/16,c(2,2,2,2))
#' compute.S(p)

compute.S <- function(p){
  d <- length(dim(p)) # number of variables
  mm <- compute.moments.2(p) # matrix of 2nd order moments in 0/1 coding
  xbar <- diag(mm)
  # compute the second moments matrix in -1/1 coding
  M <- 4*mm-2*outer(rep(1,d), xbar)-2*outer( xbar,rep(1,d))+matrix(1,d,d)
  # compute the mean vector in the -1/1 coding
  xbar <- 2*diag(mm)-1
  # compute the sample covariance matrix
  return(M-outer(xbar,xbar))
}
