#' Compute the mean vector for a given bianry distribution in -1/1 coding
#'
#' @param p A distribution
#' @return mu corresponding mean vector
#' @examples
#' p <- array(1/16,c(2,2,2,2))
#' compute.S(p)

compute.mu <- function(p){
  d <- length(dim(p)) # number of variables
  mm <- compute.moments.2(p) # matrix of 2nd order moments in 0/1 coding
  xbar <- diag(mm)
  return(2*diag(mm)-1)
}
