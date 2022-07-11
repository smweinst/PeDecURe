#' @title PeDecURe
#' @description Function to implement PeDecURe (or PLS or other penalized PLS method, depending on input specifications). Note: the code for implementing PeDecURe was adapted from Lin et al. (2016)'s implementation of AC-PCA in R (see https://github.com/linzx06/AC-PCA).
#' @param X n x p matrix of features to appear in the PeDecURe maximization term. In PeDecURe, this should be defined as X* ("X star"), i.e. the matrix of residuals after subtracting out the effects of A conditional on Y
#' @param X.penalize n x p matrix of features to appear in the penalty term.In PeDecURe, this should be defined as X~ ("X tilde"), i.e. the matrix of residuals after subtracting out the effects of Y conditional on A. If X.penalize is not specified, no penalization will be done.
#' @param A matrix of nuisance variables
#' @param Y vector specifying outcome variable of interest (e.g., disease status)
#' @param lambda tuning parameter (to be determined in other function)
#' @param nPC number of primary components (PCs) to be extracted
#' @param centerX TRUE/FALSE whether to center X by column mean (note: some centering/scaling may have been done before inputting data into pedecure function)
#' @param centerA TRUE/FALSE whether to center A by column mean
#' @param centerY TRUE/FALSE whether to center Y by mean
#' @param scaleX TRUE/FALSE whether to scale each column of X by its standard deviation
#' @param scaleA TRUE/FALSE whether to scale each column of A by its standard deviation
#' @param scaleY TRUE/FALSE whether to scale Y by its standard deviation
#' @export
#' @importFrom RSpectra eigs eigs_sym
#' @return PC loading vectors

pedecure = function(X, X.penalize = NULL, A, Y, lambda, nPC, centerX = F, centerA = F, centerY = F, scaleX = F, scaleA = F, scaleY = F){

  X = scale(X, center = centerX, scale = scaleX)
  Y = scale(Y, center = centerY, scale = scaleY)

  if (!is.null(X.penalize) & lambda > 0){ # if X.penalize is specified, then there will be a penalty term. otherwise, no penalty term
    #cat("\n penalty term supplied.\n penalized decomposition maximizes v'X'YY'Xv - lambda v'X.penalize'A A'X.penalize v")

    X.penalize = scale(X.penalize, center = centerX, scale = scaleX) # center and scale - if TRUE, then this entails centering/scaling *after* residualization
    A = scale(A, center = centerA, scale = scaleA)

    # kernel matrix for penalization (assume linear):
    eig_out = eigs_sym(pedecureMax,
                       k=nPC,
                       which = "LA",
                       n=ncol(X),
                       args = list(X = X,
                                   X.penalize = X.penalize,
                                   Y = Y,A = A,
                                   lambda = lambda))

  } else{
    #cat("\n will only maximize X'YY'X because: ")

    # if (is.null(X.penalize)){
    #   cat("\n ... X.penalize not supplied.")
    # }
    # if (lambda==0){
    #   cat("\n ... lambda = 0.")
    # }

    eig_out = eigs_sym(plsMax,
                       k=nPC,
                       which = "LA",
                       n=ncol(X),
                       args = list(X = X, Y = Y))



  }

}

