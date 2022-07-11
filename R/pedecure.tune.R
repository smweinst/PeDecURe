#' @title pedecure.tune
#' @description Function to tune lambda for PeDecURe. note: input arguments are similar but not identical to those in pedecure() function
#' @param X.orig n x p matrix of ~original~ features (not residualized). this is the version of the data that's used to obtain PC scores
#' @param X.max n x p matrix of features that goes into the maximization term (X* ("X star") for PeDecURe)
#' @param X.penalize n x p matrix of features to appear in the penalty term.In PeDecURe, this should be defined as X~ ("X tilde"), i.e. the matrix of residuals after subtracting out the effects of Y conditional on A.
#' @param A matrix of nuisance variables
#' @param Y vector specifying outcome variable of interest (e.g., disease status)
#' @param lambdas vector including candidate values of tuning parameter (e.g., lambdas = seq(0,10,by=0.1))
#' @param nPC number of primary components (PCs) to be extracted
#' @param centerX TRUE/FALSE whether to center X by column mean (note: some centering/scaling may have been done before inputting data into pedecure function)
#' @param centerA TRUE/FALSE whether to center A by column mean
#' @param centerY TRUE/FALSE whether to center Y by mean
#' @param scaleX TRUE/FALSE whether to scale each column of X by its standard deviation
#' @param scaleA TRUE/FALSE whether to scale each column of A by its standard deviation
#' @param scaleY TRUE/FALSE whether to scale Y by its standard deviation
#' @param n.cores number of cores for parallelization in parallel::mclapply()
#' @param plot TRUE/FALSE whether to plot the tuning function for different values of lambda
#' @param thresh threshold for deciding when lambda is good enough (unit is proportion change, not percentage)
#' @export
#' @importFrom RSpectra eigs eigs_sym
#' @importFrom parallel mclapply
#' @return best lambda

pedecure.tune = function(X.orig, X.max, X.penalize, lambdas, A, Y, nPC, centerX = F, scaleX = F, centerY = F, scaleY = F, centerA = F, scaleA = F, n.cores = 1, plot = FALSE, thresh = 0.001){

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

      # partial correlations:
      temp_pcor = partial.cor(scores = temp_scores, A = A, Y = Y)$partial$estimates

      # weighted sum of difference of absolute values of partial correlations between Y or A and scores (partial correlations conditional on A or Y)
      ## weighted by proportion of variation explained (eigenvalues / sum of eigenvalues)
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
