fit_true_causal_exposures <- function(model,observed,fit,theta,cohort='',dataset='') {
  env <- new.env(parent=globalenv())
  sys.source(file.path(model$root,'R/load_native.R'),env)
  sys.source(file.path(model$root,'R/causal_exposure_reference.R'),env)
  reference <- env$make_support_reference(model$root,file.path(model$root,'R/point_estimation.R'),file.path(model$work_dir,'reference'),TRUE)
  k <- length(theta)
  causal <- which(theta!=0)
  traits <- observed$dat$traits
  if(length(traits)!=k) traits <- paste0('exposure_',seq_len(k))
  rows <- list()
  results <- list()
  for(panel in names(fit$fits)) {
    base <- fit$fits[[panel]]$joint$result$score
    ii <- fit$fits[[panel]]$retained
    source_ok <- isTRUE(base$selector_converged)&&isTRUE(base$selection$converged)
    invalid <- if(source_ok) ii[base$selection$invalid] else integer()
    a <- if(source_ok) mr_susie_capture(function() reference$fit(observed$dat,causal,invalid,ii,'reference')) else list(result=NULL,error='joint_selection_unavailable',warnings=character())
    z <- a$result
    point <- center <- se <- p <- critical <- rep(NA_real_,k)
    lower <- upper <- rep(NA_real_,k)
    available <- test <- rep(FALSE,k)
    known <- setdiff(seq_len(k),causal)
    point[known] <- 0
    available[known] <- TRUE
    if(isTRUE(z$identified)&&length(causal)) {
      point[causal] <- z$point$beta
      center[causal] <- z$score$beta
      se[causal] <- z$score$se
      p[causal] <- z$score$p
      critical[causal] <- z$score$critical_95
      lower[causal] <- center[causal]-critical[causal]*se[causal]
      upper[causal] <- center[causal]+critical[causal]*se[causal]
      available[causal] <- TRUE
      test[causal] <- TRUE
    }
    rows[[panel]] <- data.frame(cohort=cohort,dataset=dataset,panel=panel,mode='ME',
      method='True causal exposures',estimand='conditional_direct_effect',exposure=seq_len(k),trait=traits,
      beta=point,score_center=center,se=se,p=p,critical=critical,lower=lower,upper=upper,
      point_available=available,test_available=test,status=ifelse(seq_len(k)%in%known,'Given',if(isTRUE(z$identified)) 'estimated' else 'unavailable'),
      warnings=paste(a$warnings,collapse=' | '),error=if(length(a$error)) a$error else '',stringsAsFactors=FALSE)
    results[[panel]] <- a
  }
  list(rows=do.call(rbind,rows),fits=results)
}
