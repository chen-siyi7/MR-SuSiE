fit_cml_bic <- function(observed,root,seed=1L) {
  env <- new.env(parent=globalenv())
  for(file in c('cml_selection.R','cml_estimation.R')) sys.source(file.path(root,'R',file),env)
  d <- observed$dat
  k <- d$num_exposures
  mr_susie_capture(function() {
    selected <- d
    selected$betaZY <- observed$selection$beta
    selected$se_betaZY <- observed$selection$se
    selected$n2 <- observed$selection$n_eff
    sel <- env$native_cml(selected,seed)
    if(!isTRUE(sel$completed)) stop(sel$status)
    al <- sel$step3$alpha
    sc <- sel$screened
    ids <- integer()
    if(length(sc)&&!is.null(al)) for(i in seq_len(nrow(al))) {
      eligible <- which(al[i,]>1/length(sc))
      if(length(eligible)) ids <- c(ids,min(sc[eligible[al[i,eligible]==max(al[i,eligible])]]))
    }
    ids <- sort(unique(ids))
    if(!length(ids)) stop('no_native_exposure_selected')
    x <- env$model_data4(d,ids)
    set.seed(seed+1L)
    b <- env$bic4(x)
    result <- list(beta=rep(NA_real_,k),se=rep(NA_real_,k),p=rep(NA_real_,k),
      selection=sel,selected=ids,fit=b,criterion='BIC',data_perturbation=FALSE)
    result$beta[ids] <- b$beta
    result$se[ids] <- b$se
    result$p[ids] <- 2*pnorm(-abs(b$beta/b$se))
    result
  })
}
