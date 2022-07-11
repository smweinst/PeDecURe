pedecure = function(X, X.penalize = NULL, A, Y, lambda, nPC, centerX = F, centerA = F, centerY = F, scaleX = F, scaleA = F, scaleY = F){
  # X = the X that appears in maximization term
  # X.penalize = the X that appears in penalty term
  # A = matrix of nuisance variables
  # Y = matrix of outcome variable / variable of interest
  # lambda = tuning paramter
  # nPC = number of directions of variation to extract
  # centering and scaling - specify if not yet centered/scaled as desired before starting the function

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

pedecureMax = function(v,args){
  X = args$X
  X.penalize = args$X.penalize
  Y = args$Y
  A = args$A
  lambda = args$lambda # tuning parameter lambda
  XTY=crossprod(X,Y)
  K=A%*%t(A)

  out = XTY%*%t(XTY)%*%v - lambda*t(X.penalize)%*%K%*%X.penalize%*%v # based on output from calAv function in acPCA R package (Lin et al., 2016)

  # out = tcrossprod(crossprod(X,Y))%*%matrix(v) - lambda*t(X.penalize)%*% tcrossprod(A) %*%X.penalize%*%matrix(v)

  return(out)
}

plsMax = function(v,args){
  X = args$X
  Y = args$Y

  out = tcrossprod(crossprod(X,Y))%*%matrix(v,ncol=1)
  return(out)

}

pedecure.tune = function(X.orig, X.max, X.penalize, lambdas, A, Y, nPC, centerX = F, scaleX = F, centerY = F, scaleY = F, centerA = F, scaleA = F, n.cores = 1, plot = TRUE, thresh = 0.001){

  sum_abs_diff = vector(mode = "numeric",length = length(lambdas))
  names_A = colnames(A)

  lambdas=sort(lambdas) # sort in case they started out of order (want to consider lambdas in ascending order)

  split_lambdas = matrix(c(lambdas,rep(NA,2-length(lambdas)%%2)), ncol=2,byrow=T)
  best_lambda = best_weighted_sum_abs_diff = vector(mode="numeric",length=nrow(split_lambdas))
  all_weighted_sum_abs_diff = c()
  for (lambda_group in 1:nrow(split_lambdas)){
    weighted_sum_abs_diff = unlist(mclapply(split_lambdas[lambda_group,which(!is.na(split_lambdas[lambda_group,]))], FUN = function(l){
      # penalized decomposition using residuals for a given lambda:
      temp_out = pedecure(X = X.max, X.penalize = X.penalize, A = A, Y = Y, lambda = l, nPC = nPC,
                          centerX = centerX, scaleX = scaleX, centerA = centerA, scaleA = scaleA, centerY = centerY, scaleY=scaleY)

      # in-sample scores for lambda = l:
      temp_scores = X.orig%*%temp_out$vectors # 11/03/2021 - not re-centering and scaling X here because this would have been done already
      #temp_scores = scale(X.orig,center = centerX, scale=scaleX)%*%temp_out$vectors # scores use original data (even though residuals are used in maximization and penalty terms)

      # partial correlations:
      temp_pcor = my.pcor.fun(scores = temp_scores, A = A, Y = Y)$partial$estimates

      # weighted sum of difference of absolute values of partial correlations between Y or A and scores (partial correlations conditional on A or Y)
      ## weighted by proportion of variation explained (eigenvalues / sum of eigenvalues)
      #weighted_sum_abs_diff_pcor = (temp_out$values/sum(temp_out$values))%*%(abs(temp_pcor["Y",]) - abs(temp_pcor["A",]))
      weighted_sum_abs_diff_pcor = sapply(1:ncol(A), FUN = function(a){
        (temp_out$values/sum(temp_out$values))%*%(abs(temp_pcor["Y",]) - abs(temp_pcor[names_A[a],]))
      })

      return(mean(weighted_sum_abs_diff_pcor)) # taking the average of the weighted sum of absolute differences (assuming each confounder is equally prioritized?)

    }, mc.cores = n.cores))

    all_weighted_sum_abs_diff = c(all_weighted_sum_abs_diff,weighted_sum_abs_diff)
    best_weighted_sum_abs_diff[lambda_group] = max(weighted_sum_abs_diff)

    best_lambda[lambda_group] = max(split_lambdas[lambda_group,which.max(weighted_sum_abs_diff)])

    if (lambda_group > 2 & lambda_group < nrow(split_lambdas)-1){
      prop_change = (best_weighted_sum_abs_diff[lambda_group]-best_weighted_sum_abs_diff[lambda_group-1])/best_weighted_sum_abs_diff[lambda_group-1]
      if (prop_change < thresh){
        best_lambda[(lambda_group+1):length(best_lambda)] = NA
        break
      }
    }
  }

  return(list(weighted_sum_abs_diff = all_weighted_sum_abs_diff,
              lambdas = lambdas[1:length(all_weighted_sum_abs_diff)],
              lambda_tune = best_lambda[lambda_group]))

}
