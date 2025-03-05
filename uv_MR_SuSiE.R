MR_SuSiE_sel = function(n1, n2, Y2, Z2, corr, gamma_hat){
  res = susie(cbind(Z2,Z2%*%gamma_hat), Y2, L = dim(Z2)[2] + 1, intercept = F)
  Iu = sort(as.vector(unlist(res$sets$cs)))
  Iu = Iu[Iu <= dim(Z2)[2]]
  Ix0 = setdiff(1:dim(corr)[1], Iu)
  return(Ix0)
}

sum_MR_SuSiE_sel = function(corr, betaZX, betaZY, se_betaZY, gamma_hat, n1, n2){
  ZTZ = corr
  p = dim(corr)[1]
  ZTY = matrix(diag(ZTZ) * betaZY, ncol = 1);
  YTY = NULL
  for(SNP in 1:p){
    YTY[SNP] = (n2-1)*ZTZ[SNP, SNP]*(se_betaZY^2)[SNP] + ZTY[SNP]*betaZY[SNP]
  }
  YTY = median(YTY)
  res = susie_suff_stat(XtX = corr, Xty = ZTY, yty = YTY, n = max(n1, n2), L = dim(corr)[2] + 1, estimate_residual_variance = FALSE)
  Iu = sort(as.vector(unlist(res$sets$cs)))
  Iu = Iu[Iu <= p]
  Ix0 = setdiff(1:dim(corr)[1], Iu)
  return(Ix0)
}

