make_mr_susie_engine <- function(source_root, work_dir, point_rule_path) {
  e <- new.env(parent=environment())
  sys.source(file.path(source_root,'R/score_engine.R'),e)
  sys.source(point_rule_path,e)
  engine <- e$make_score_engine(source_root,list(runtime=list(work_dir=work_dir)))
  from_score <- function(score,M) {
    point <- e$posterior_mean_estimate(score,'EB',M)
    list(method=if(M==1L)'UVMR-SuSiE'else'MVMR-SuSiE',M=M,
         point=point$beta,point_available=isTRUE(point$identified),
         score=score,point_fit=point,
         estimand=if(M==1L)'single_exposure_model'else'conditional_direct_effect')
  }
  fit <- function(observed,controls) {
    M <- ncol(observed$dat$betaZX)
    stopifnot(M>=1L,observed$dat$num_exposures==M,
              identical(dim(observed$dat$seZX),dim(observed$dat$betaZX)),
              ncol(observed$selection$gamma_hat)==M)
    score <- engine$fit(observed,controls)
    from_score(score,M)
  }
  subset_exposure <- function(observed,k) {
    stopifnot(length(k)==1L,k>=1L,k<=ncol(observed$dat$betaZX))
    z<-observed;d<-z$dat
    d$num_exposures<-1L;d$betaZX<-d$betaZX[,k,drop=FALSE]
    d$seZX<-d$seZX[,k,drop=FALSE];d$sigma_exp<-matrix(1)
    if(length(d$traits))d$traits<-d$traits[k]

    z$selection$gamma_hat<-z$selection$gamma_hat[,k,drop=FALSE]
    z$dat<-d;z
  }
  uvmr <- function(observed,controls,k=1L)fit(subset_exposure(observed,k),controls)

  list(fit=fit,uvmr=uvmr,mvmr=fit,subset_exposure=subset_exposure,
       from_score=from_score,fixed_support_profile=get('profile',environment(engine$fit))$corrected_score_refit)
}
