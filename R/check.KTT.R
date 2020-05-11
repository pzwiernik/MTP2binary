#' Check the KKT conditions for the quadratic Ising model
#'
#' @param phat The sampling distribution
#' @param p The MLE candidate
#' @return A report on the KKT conditions
#' @examples
#' p <- array(1/16,c(2,2,2,2))
#' compute.J(p)

check.KKT <- function(xbar,M,p,round=10){
S <- M-outer(xbar,xbar)
mu <- compute.mu(p)
Shat <- compute.S(p)
J <- compute.J(p)
print('The matrix J should be nonnegative')
print(round(J,round))
print('Sigma hat should dominate the sample covariance. Here (Sigma hat - S) is ')
print(round(Shat-S,round))
print('The sample mean should be equal to mu. Here (xbar - mu) is ')
print(round(xbar-mu,round))
print('The complementary slackness holds if all entries below are zero')
print(round(J*(Shat-S),round))
}
