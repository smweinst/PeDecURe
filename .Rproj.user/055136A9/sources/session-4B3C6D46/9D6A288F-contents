get.X.resid = function(X, Y, A, Z = NULL){
  X = as.matrix(X) # make sure it's a matrix
  A = as.matrix(A)

  X.tilde = c()
  X.star = c()
  for (j in 1:ncol(X)){
    if (!is.null(Z)){ # including extra covariates (Z) is optional
      lm.X.j = lm(X[,j] ~ Y + A + Z)
    } else{
      lm.X.j = lm(X[,j] ~ Y + A)
    }
    X.tilde = cbind(X.tilde, (X[,j] - lm.X.j$coefficients["Y"]*Y))
    X.star = cbind(X.star, (X[,j] - A%*%lm.X.j$coefficients[startsWith(names(lm.X.j$coefficients),"A")]))
  }

  return(list(X.tilde = X.tilde,
              X.star = X.star))
}
