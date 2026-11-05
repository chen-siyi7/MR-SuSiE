fit_comparators <- function(observed,root,seed=1L,native_seed=seed) {
  env <- new.env(parent=globalenv())
  for(file in c('comparators.R','estimators.R','cml_selection.R','cml_estimation.R')) sys.source(file.path(root,'R',file),env)
  capture <- mr_susie_capture
  d <- observed$dat
  k <- d$num_exposures
  out <- list()
  fail <- function(reason) list(result=NULL,error=reason,warnings=character(),elapsed=0)
  rank_ok <- d$m>k && qr(d$betaZX,tol=1e-10)$rank==k
  names <- c('MVMR-IVW','MVMR-Egger','MVMR-Lasso','MVMR-median','MVMR-robust')
  input <- if(requireNamespace('MendelianRandomization',quietly=TRUE)) MendelianRandomization::mr_mvinput(bx=d$betaZX,bxse=d$seZX,by=d$betaZY,byse=d$se_betaZY) else NULL
  for(i in seq_along(names)) {
    name <- names[i]
    out[[name]] <- if(!rank_ok) fail('insufficient_joint_design_rank_or_variants') else capture(function() {
      if(is.null(input)) stop('MendelianRandomization is required')
      set.seed(seed+i)
      z <- switch(name,
        'MVMR-IVW'=MendelianRandomization::mr_mvivw(input,model='random',robust=FALSE,correl=FALSE),
        'MVMR-Egger'=MendelianRandomization::mr_mvegger(input,orientate=1,correl=FALSE),
        'MVMR-Lasso'=MendelianRandomization::mr_mvlasso(input,orientate=1),
        'MVMR-median'=MendelianRandomization::mr_mvmedian(input,iterations=1000,seed=seed+i),
        'MVMR-robust'=MendelianRandomization::mr_mvivw(input,model='random',robust=TRUE,correl=FALSE))
      list(beta=z@Estimate,se=if(name=='MVMR-Egger') z@StdError.Est else z@StdError,
        p=if(name=='MVMR-Egger') z@Pvalue.Est else z@Pvalue)
    })
  }
  functions <- list(
    'MVMR-IVW-LD'=function() env$mvmr_ivw_ld_adjusted(d,integer(),'ref',FALSE),
    'MVMR-IVW-LD-ME'=function() env$mvmr_ivw_ld_adjusted(d,integer(),'ref',TRUE,sigma_exp=d$sigma_exp),
    'MVMR-Egger-LD'=function() env$cmp_mvmr_egger_ld(d,which_R='ref'),
    'MV-IWAS'=function() env$cmp_mv_iwas(d,which_R='ref'),
    'IWAS-Egger-projection'=function() {
      z <- d
      z$betaZX <- cbind(intercept=1,d$betaZX)
      z$num_exposures <- k+1L
      result <- env$cmp_mv_iwas(z,which_R='ref')
      result$intercept <- result$beta[1]
      for(field in c('beta','se','p','q')) result[[field]] <- result[[field]][-1]
      result
    })
  for(name in names(functions)) out[[name]] <- if(rank_ok) capture(functions[[name]]) else fail('insufficient_joint_design_rank_or_variants')
  out[['TWAS-summary-Z']] <- capture(function() {
    B <- MASS::ginv(d$corr_ref)%*%d$betaZX
    variance <- diag(crossprod(B,d$corr_ref%*%B))
    z <- as.numeric(crossprod(B,d$betaZY/d$se_betaZY))/sqrt(variance)
    list(beta=rep(NA_real_,k),se=rep(NA_real_,k),p=2*pnorm(-abs(z)),z=z,status='marginal_summary_Z')
  })
  sys.source(file.path(root,'R/cml_model.R'),environment())
  out[['MVMR-cML-SuSiE']] <- fit_cml_bic(observed,root,native_seed)
  out
}
