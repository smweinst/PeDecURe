my.pcor.fun = function(scores, A, Y){
  scores = as.matrix(scores)
  A = as.matrix(A)
  colnames(scores)  = paste0("PC",1:ncol(scores))
  n = nrow(scores)
  df.partial = n-2-ncol(A) # n - 2 - number of variables conditioning on

  pcor.Xv.j.Y.given.A = vector(mode = "numeric",length = ncol(scores))
  pcor.Xv.j.A.given.Y = rep(list(vector(mode = "numeric", length = ncol(scores))), ncol(A))
  names(pcor.Xv.j.A.given.Y) = colnames(A)#paste0("A",1:ncol(A))
  #pcor.Xv.j.Y.given.A = pcor.Xv.j.A.given.Y = vector(mode = "numeric",length = ncol(scores))

  t.statistics = matrix(nrow=(ncol(A) + 1),ncol=ncol(scores))
  rownames(t.statistics) = c(paste0(colnames(A)),"Y")
  colnames(t.statistics) = paste0("PC",1:ncol(scores))

  for (j in 1:ncol(scores)){
    omega.j = cov(cbind(A,Y,scores[,j]),method = "pearson")
    colnames(omega.j) = rownames(omega.j) = c(rownames(t.statistics),paste0("Xv",j))
    P.j = solve(omega.j)
    pcor.estimate.j = -cov2cor(P.j)
    pcor.Xv.j.Y.given.A[j] = pcor.estimate.j["Y",(ncol(A)+2)]

    for (a in 1:ncol(A)){
      pcor.Xv.j.A.given.Y[[colnames(A)[a]]][j] = pcor.estimate.j[colnames(A)[a],(ncol(A)+2)]
    }

  }



  partial.correlations = rbind(do.call("rbind",pcor.Xv.j.A.given.Y),
                               Y = pcor.Xv.j.Y.given.A)

  colnames(partial.correlations) = paste0("PC",1:ncol(scores))

  t.statistics = partial.correlations*sqrt(df.partial/(1-partial.correlations^2))
  colnames(t.statistics) = paste0("PC",1:ncol(scores))

  p.values = 2*pt(-abs(t.statistics),df.partial)

  out.partial = list(estimates = partial.correlations,
                     statistics = t.statistics,
                     df = df.partial,
                     p.values = p.values)

  # marginal correlations
  cor.Xv.A = matrix(apply(scores, 2, cor, A),nrow = ncol(A)); rownames(cor.Xv.A) = paste0("A",1:ncol(A))
  cor.Xv.Y = apply(scores, 2, cor, Y)
  marginal.correlations = rbind(cor.Xv.A,
                                Y = cor.Xv.Y)
  df.marg = n-2

  marg.cor.statistics = marginal.correlations*sqrt(df.marg/(1-marginal.correlations^2))
  marg.cor.p.values = 2*pt(-abs(marg.cor.statistics),df.marg)

  out.marginal = list(estimates = marginal.correlations,
                      statistics = marg.cor.statistics,
                      df = df.marg,
                      p.values = marg.cor.p.values)



  return(list(partial = out.partial,
              marginal = out.marginal))


}
