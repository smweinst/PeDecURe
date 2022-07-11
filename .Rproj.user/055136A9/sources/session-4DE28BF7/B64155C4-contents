#' @title plsMax
#' @description maximization function for PLS, when no penalty term is specified in pedecure() function
#' @param v PC loading vector that will be maximized via eigs_sym
#' @param args list of arguments (inputted from within pedecure() function)
#' @export
#' @return v which maximizes v'X'YY'Xv

plsMax = function(v,args){
  X = args$X
  Y = args$Y

  out = tcrossprod(crossprod(X,Y))%*%matrix(v,ncol=1)
  return(out)

}
