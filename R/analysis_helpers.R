mr_susie_arguments <- function(args) {
  result <- list()
  for(a in args) {
    z <- strsplit(sub('^--','',a),'=',fixed=TRUE)[[1]]
    if(length(z)>=2L) result[[z[1]]] <- paste(z[-1],collapse='=')
  }
  result
}

mr_susie_capture <- function(fun) {
  warnings <- character()
  error <- ''
  start <- proc.time()[3]
  value <- tryCatch(withCallingHandlers(fun(),warning=function(w) {
    warnings <<- c(warnings,conditionMessage(w))
    invokeRestart('muffleWarning')
  }),error=function(e) {
    error <<- conditionMessage(e)
    NULL
  })
  list(result=value,warnings=unique(warnings),error=error,elapsed=unname(proc.time()[3]-start))
}

mr_susie_clump <- function(d,r2=.01) {
  R <- cov2cor(d$corr_ref)
  score <- apply(abs(d$betaZX/d$seZX),1,max)
  active <- rep(TRUE,d$m)
  keep <- integer()
  for(i in order(-score,seq_len(d$m))) if(active[i]) {
    keep <- c(keep,i)
    active[which(R[,i]^2>r2+1e-12)] <- FALSE
  }
  sort(keep)
}

mr_susie_subset_rows <- function(o,ii) {
  d <- o$dat
  d$m <- length(ii)
  for(n in intersect(c('betaZX','seZX'),names(d))) d[[n]] <- d[[n]][ii,,drop=FALSE]
  for(n in intersect(c('betaZY','se_betaZY','rsid','reference_dosage_sd'),names(d))) d[[n]] <- d[[n]][ii]
  for(n in c('corr_in','corr_ref')) d[[n]] <- d[[n]][ii,ii,drop=FALSE]
  o$dat <- d
  o$selection$beta <- o$selection$beta[ii]
  o$selection$se <- o$selection$se[ii]
  o$selection$gamma_hat <- MASS::ginv(d$corr_ref)%*%d$betaZX
  o
}

mr_susie_rows <- function(f,method,estimand,traits,cohort='',dataset='',panel='',mode='',warnings=character(),error='') {
  k <- length(traits)
  val <- function(x) if(length(x)==k) as.numeric(x) else rep(NA_real_,k)
  if(method %in% c('MVMR-SuSiE','UVMR-SuSiE')) {
    s <- f$score
    beta <- val(f$point)
    center <- val(s$beta)
    critical <- if(length(s$critical_95)==1L) s$critical_95 else NA_real_
    available <- rep(isTRUE(f$point_available),k)&is.finite(beta)
    test <- rep(isTRUE(s$identified)&&isTRUE(s$converged)&&isTRUE(s$selector_converged),k)
  } else {
    s <- f
    beta <- center <- val(s$beta)
    critical <- if(length(s$critical_95)==1L) s$critical_95 else qnorm(.975)
    available <- is.finite(beta)&!nzchar(error)
    test <- !nzchar(error)&!any(grepl('did not converge|not converg|non.?converg',warnings,ignore.case=TRUE))
  }
  se <- val(s$se)
  p <- val(s$p)
  test <- test&is.finite(p)
  if(method %in% c('MVMR-SuSiE','UVMR-SuSiE')) test <- test&is.finite(se)&se>0
  status <- if(length(s$status)) as.character(s$status)[1] else if(nzchar(error)) error else 'returned'
  data.frame(cohort=cohort,dataset=dataset,panel=panel,mode=mode,method=method,
    estimand=estimand,exposure=seq_len(k),trait=traits,beta=beta,score_center=center,
    se=se,p=p,critical=critical,lower=center-critical*se,upper=center+critical*se,
    point_available=available,test_available=test,status=status,
    warnings=paste(unique(c(warnings,s$selection_warnings)),collapse=' | '),error=error,
    stringsAsFactors=FALSE)
}

mr_susie_read_input <- function(path,model,n_eff=NULL,application=FALSE) {
  if(grepl('[.]json$',path)) {
    x <- jsonlite::read_json(path,simplifyVector=TRUE)
    k <- ncol(x$B)
    m <- nrow(x$B)
    d <- list(m=as.integer(m),num_exposures=as.integer(k),betaZX=x$B,
      seZX=if(length(x$sx)==1L) matrix(x$sx,m,k) else as.matrix(x$sx),
      betaZY=x$y,se_betaZY=rep(x$sy,length.out=m),corr_in=x$R,corr_ref=x$R,
      sigma_exp=x$Sigma,overlap=0,n1=10000L,n2=10000L)
    return(model$observed(d,x$selection_y,rep(x$sy,m),n_eff,
      if(application) 'binary_score' else 'gaussian',FALSE,n_reference=Inf))
  }
  x <- readRDS(path)
  if(!is.null(x$dat)&&!is.null(x$selection)) return(x)
  d <- if(!is.null(x$betaZX)) x else if(!is.null(x$data)) x$data else x$generated$dat
  stopifnot(is.list(d),!is.null(d$betaZX))
  if(application) d$n2 <- n_eff
  model$observed(d,n_select=if(application) n_eff else d$n2,
    selector_family=if(application) 'binary_score' else 'gaussian',
    outcome_reuse=TRUE,n_reference=if(application) 503 else Inf)
}

mr_susie_analyze <- function(observed,model,controls=list(L=10,prior=.5,coverage=.95,purity=.5),
    context=list(kind='simulation',family='broader'),cohort='',dataset='',mode='ME',comparators=FALSE,seed=1L,me_sensitivity=FALSE,native_seed=seed) {
  keep <- mr_susie_clump(observed$dat)
  panel_inputs <- list(original_full_LD=observed,matched_r2_0.01=if(length(keep)==observed$dat$m) observed else mr_susie_subset_rows(observed,keep))
  rows <- list()
  fits <- list()
  for(panel in names(panel_inputs)) {
    if(panel=='matched_r2_0.01'&&length(keep)==observed$dat$m) {
      fits[[panel]] <- fits[['original_full_LD']]
      old_keys <- names(rows)
      for(key in old_keys) if(startsWith(key,'original_full_LD ')) {
        z <- rows[[key]]
        z$panel <- panel
        rows[[sub('original_full_LD','matched_r2_0.01',key,fixed=TRUE)]] <- z
      }
      next
    }
    o <- panel_inputs[[panel]]
    k <- o$dat$num_exposures
    traits <- o$dat$traits
    if(length(traits)!=k) traits <- paste0('exposure_',seq_len(k))
    set.seed(seed)
    f <- mr_susie_capture(function() model$fit_observed(o,list(baseline=controls,context=context)))
    name <- if(k==1L) 'UVMR-SuSiE' else 'MVMR-SuSiE'
    fits[[panel]] <- list(joint=f,marginal=list(),comparators=list(),retained=if(panel=='original_full_LD') seq_len(observed$dat$m) else keep)
    rows[[paste(panel,name)]] <- mr_susie_rows(f$result,name,'conditional_direct_effect',traits,cohort,dataset,panel,mode,f$warnings,f$error)
    if(k>1L) for(i in seq_len(k)) {
      set.seed(seed+i)
      u <- mr_susie_capture(function() model$uvmr(o,list(baseline=controls,context=context),i))
      fits[[panel]]$marginal[[as.character(i)]] <- u
      z <- mr_susie_rows(u$result,'UVMR-SuSiE','marginal_association',traits[i],cohort,dataset,panel,mode,u$warnings,u$error)
      z$exposure <- i
      rows[[paste(panel,'uvmr',i)]] <- z
    }
    if(me_sensitivity&&identical(context$kind,'application')) {
      me_from <- function(input,base) {
        v <- model$fixed_support_profile(input$dat,base$score$selection,base$score$selection_warnings)
        v$selector_converged <- base$score$selector_converged
        v$selection <- base$score$selection
        model$from_score(v,input$dat$num_exposures)
      }
      me <- mr_susie_capture(function() me_from(o,f$result))
      fits[[panel]]$ME <- me
      rows[[paste(panel,name,'ME')]] <- mr_susie_rows(me$result,name,'conditional_direct_effect',traits,cohort,dataset,panel,'ME',me$warnings,me$error)
      if(k>1L) for(i in seq_len(k)) {
        u <- fits[[panel]]$marginal[[as.character(i)]]
        z <- mr_susie_capture(function() me_from(model$subset_exposure(o,i),u$result))
        fits[[panel]]$marginal_ME[[as.character(i)]] <- z
        rr <- mr_susie_rows(z$result,'UVMR-SuSiE','marginal_association',traits[i],cohort,dataset,panel,'ME',z$warnings,z$error)
        rr$exposure <- i
        rows[[paste(panel,'uvmr',i,'ME')]] <- rr
      }
    }
    if(comparators) {
      cc <- fit_comparators(o,model$root,seed,native_seed)
      fits[[panel]]$comparators <- cc
      for(nm in names(cc)) rows[[paste(panel,nm)]] <- mr_susie_rows(cc[[nm]]$result,nm,
        if(nm=='TWAS-summary-Z') 'marginal_association' else 'conditional_direct_effect',
        traits,cohort,dataset,panel,mode,cc[[nm]]$warnings,cc[[nm]]$error)
    }
  }
  list(rows=do.call(rbind,rows),fits=fits)
}
