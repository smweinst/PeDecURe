#' @title pedecureMax
#' @description maximization function for PeDecURe (called within function pedecure()). Note: the code for implementing PeDecURe was adapted from Lin et al. (2016)'s implementation of AC-PCA in R (see https://github.com/linzx06/AC-PCA).
#' @param v PC loading vector that will be maximized via eigs_sym
#' @param args list of arguments (inputted from within pedecure() function)
#' @export
#' @return v which maximizes objective function

pedecureMax = function(v,args){
  X = args$X
  X.penalize = args$X.penalize
  Y = args$Y
  A = args$A
  lambda = args$lambda # tuning parameter lambda
  XTY=crossprod(X,Y)
  K=A%*%t(A)

  out = XTY%*%t(XTY)%*%v - lambda*t(X.penalize)%*%K%*%X.penalize%*%v # based on output from calAv function in acPCA R package (Lin et al., 2016)

  return(out)
}
