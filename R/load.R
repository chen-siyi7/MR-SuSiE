make_mr_susie <- function(root=NULL, work_dir=file.path(tempdir(),'mr_susie')) {
  if (is.null(root)) {
    root <- system.file('engine', package='MRSuSiE')
    if (!nzchar(root)) root <- '.'
  }
  root <- normalizePath(root, mustWork=TRUE)
  stopifnot(file.exists(file.path(root,'R/selector.R')))
  available <- c(
    Rcpp=requireNamespace('Rcpp',quietly=TRUE),
    RcppArmadillo=requireNamespace('RcppArmadillo',quietly=TRUE),
    susieR=requireNamespace('susieR',quietly=TRUE),
    MASS=requireNamespace('MASS',quietly=TRUE),
    jsonlite=requireNamespace('jsonlite',quietly=TRUE),
    digest=requireNamespace('digest',quietly=TRUE)
  )
  missing <- names(available)[!available]
  if(length(missing)) stop('Missing R dependencies: ',paste(missing,collapse=', '))
  stopifnot(as.character(utils::packageVersion('susieR'))=='0.14.2')
  e <- new.env(parent=globalenv())
  for(file in c('load_native.R','method_engine.R','observed_inputs.R')) {
    sys.source(file.path(root,'R',file),envir=e)
  }
  engine <- e$make_mr_susie_engine(root,work_dir,file.path(root,'R/point_estimation.R'))
  format_result <- function(fit,units) {
    score <- fit$score
    k <- fit$M
    value <- function(x) if(length(x)==k) as.numeric(x) else rep(NA_real_,k)
    critical <- if(length(score$critical_95)==1L) score$critical_95 else NA_real_
    center <- value(score$beta)
    se <- value(score$se)
    point <- value(fit$point)
    test_available <- isTRUE(score$identified) && isTRUE(score$converged) && isTRUE(score$selector_converged)
    p <- value(score$p)
    family <- data.frame(raw=p,available=test_available & is.finite(p))
    z <- ifelse(family$available,p,1)
    for(method in c('holm','bonferroni','BH')) {
      family[[method]] <- stats::p.adjust(z,method)
      family[[method]][!family$available] <- NA_real_
    }
    list(method=fit$method,point=point,point_available=isTRUE(fit$point_available),
      coefficient_units=units,
      score_inference=data.frame(center=center,se=se,lower=center-critical*se,upper=center+critical*se,p=p),
      family_tests=family,selection=score$selection,warnings=score$selection_warnings,
      status=score$status,raw_fit=fit)
  }
  fit <- function(B,sX,y,sy,R,Sigma_exp,selection_beta,selection_se,n_select,
      coefficient_units,selector_family='gaussian',outcome_selection_reuse=FALSE,
      n_reference=Inf,controls=list(L=10,prior=.5,coverage=.95,purity=.5),
      context=list(kind='simulation',family='broader')) {
    B <- as.matrix(B)
    sX <- as.matrix(sX)
    k <- ncol(B)
    m <- nrow(B)
    stopifnot(k>=1L,identical(dim(sX),dim(B)),length(y)==m,length(sy)==m,
      all(dim(R)==c(m,m)),all(dim(Sigma_exp)==c(k,k)),
      length(coefficient_units)==k,all(nzchar(coefficient_units)))
    d <- list(m=as.integer(m),num_exposures=as.integer(k),betaZX=B,seZX=sX,
      betaZY=as.numeric(y),se_betaZY=as.numeric(sy),corr_in=R,corr_ref=R,
      sigma_exp=Sigma_exp,overlap=0)
    observed <- e$r2_observed(d,selection_beta,selection_se,n_select,
      selector_family,outcome_selection_reuse,n_reference=n_reference)
    format_result(engine$fit(observed,list(baseline=controls,context=context)),coefficient_units)
  }
  list(fit=fit,fit_observed=engine$fit,uvmr=engine$uvmr,
    subset_exposure=engine$subset_exposure,from_score=engine$from_score,
    fixed_support_profile=engine$fixed_support_profile,observed=e$r2_observed,
    format_result=format_result,warm=function(family) {
      core <- get('engine',environment(engine$fit))
      get('get_engine',environment(core$fit))(family)
      invisible(NULL)
    },root=root,work_dir=work_dir)
}
