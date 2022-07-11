#' @title get.resid
#' @description Function to get residuals needed for PeDecURe
#' @param X n x p feature matrix
#' @param Y vector specifying outcome variable of interest (e.g., disease status)
#' @param A matrix of nuisance variables
#' @export
#' @return list containing X.tilde (goes in penalty term for PeDecURe) and X.star (goes in maximization term for PeDecURe)

get.resid = function(X, Y, A){
  X = as.matrix(X) # make sure it's a matrix
  A = as.matrix(A)

  X.tilde = c()
  X.star = c()
  for (j in 1:ncol(X)){
    lm.X.j = lm(X[,j] ~ Y + A)
    X.tilde = cbind(X.tilde, (X[,j] - lm.X.j$coefficients["Y"]*Y))
    X.star = cbind(X.star, (X[,j] - A%*%lm.X.j$coefficients[startsWith(names(lm.X.j$coefficients),"A")]))
  }

  return(list(X.tilde = X.tilde,
              X.star = X.star))
}
