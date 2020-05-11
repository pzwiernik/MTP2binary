full.ind <- function(p){
  m <- length(dim(p))
  q <- 1
  for (i in 1:m) {
    q <- outer(q,apply(p,i,sum))
  }
  return(q)
}
