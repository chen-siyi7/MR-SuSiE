make_score_engine <- function(root,spec) {
  base<-root
  here<-file.path(root,'R')
  local<-new.env(parent=environment())
  sys.source(file.path(here,'localization_engine.R'),local)
  sys.source(file.path(here,'profile_engine.R'),local)
  sys.source(file.path(here,'binary_adapter.R'),local)

  cache<-file.path(spec$runtime$work_dir,'cpp_selector');dir.create(cache,recursive=TRUE,showWarnings=FALSE)
  code<-paste(deparse(local$make_localization_selector),collapse='\n')
  needle<-'file.path(d, "cppcache_replay_ser")'
  stopifnot(grepl(needle,code,fixed=TRUE))
  selector_factory<-eval(parse(text=sub(needle,'cache',code,fixed=TRUE)),environment())
  profile_cache<-file.path(spec$runtime$work_dir,'cpp_profile');dir.create(profile_cache,recursive=TRUE,showWarnings=FALSE)
  code<-paste(deparse(local$make_native_profile_engine),collapse='\n')
  needle<-'file.path(d, "cppcache_profile")'
  stopifnot(grepl(needle,code,fixed=TRUE))
  profile_factory<-eval(parse(text=sub(needle,'profile_cache',code,fixed=TRUE)),environment())
  profile<-profile_factory(root)
  engines<-new.env(parent=emptyenv())
  get_engine<-function(family) {
    key<-if(family=='figure2')'figure2'else if(grepl('^four_',family))'four'else if(family=='applications')'application'else'broader'
    if(exists(key,engines,inherits=FALSE))return(get(key,engines))
    subdir<-switch(key,figure2='figure2/software',four='independent_remaining/software',
      application='applications/software',broader='broader_simulations/vendor')
    v<-new.env(parent=environment());sourcefile<-file.path(root,'R/selector.R')
    sys.source(sourcefile,v)
    n<-selector_factory(root,v$mv_MR_susie_sel_joint,return_posterior=FALSE,native_updates=TRUE)
    b<-local$make_binary_score_selector(v$mv_MR_susie_sel_joint_score,n$namespace$susie_rss,FALSE)
    ans<-list(gaussian=n$select,binary=b);assign(key,ans,engines);ans
  }
  fit<-function(observed,controls) {
    app<-identical(controls$context$kind,'application')
    family<-if(app)'applications'else controls$context$family
    engine<-get_engine(family);ctl<-controls$baseline;d<-observed$dat;s<-observed$selection
    warnings<-character()
    selection<-withCallingHandlers({
      if(s$family=='gaussian') {
        args<-list(corr=d$corr_ref,betaZY=s$beta,se_betaZY=s$se,gamma_hat=s$gamma_hat,
          n2=s$n_eff,L=ctl$L,exposure_prior_mass=ctl$prior,coverage=ctl$coverage,purity=ctl$purity)
        if(family%in%c('figure2','four_baseline','four_high_correlation')) {
          args$se_gamma<-d$seZX;args$sigma_exp<-d$sigma_exp;args$ambiguous<-'unresolved'
        }else args$me_correct<-FALSE
        do.call(engine$gaussian,args)
      } else {
        stopifnot(s$family%in%c('binary_score','binary_score_approximation'))
        engine$binary(corr=d$corr_ref,z=s$beta/s$se,gamma_hat=s$gamma_hat,
          L=ctl$L,exposure_prior_mass=ctl$prior,coverage=ctl$coverage,purity=ctl$purity,
          n_eff=s$n_eff,ambiguous='unresolved')
      }
    },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
    if(app) {
      f<-local$primary_application_profile(observed,selection,profile$corrected_score_refit,unique(warnings))
    }else f<-profile$corrected_score_refit(d,selection,unique(warnings))
    f$selector_converged<-isTRUE(selection$converged)&&!any(grepl('did not converge|not converg|non.?converg',warnings,ignore.case=TRUE))
    if(!f$selector_converged)f$status<-paste(f$status,'selector_nonconvergence',sep=';')
    f$selection<-selection;f$selection_warnings<-unique(warnings)
    f$selection_family<-s$family;f$n_effective_selection<-s$n_eff
    f$reference_policy<-observed$reference;f$measurement_error_correction<-!app

    f$score_covariance<-f$inference$score_covariance
    f$quadratic_term<-f$inference$quadratic_term
    f$fixed_support_covariance<-f$inference$covariance
    f$inference<-NULL
    f$method_scope<-'Working score intervals conditional on localization and supplied LD; empirical coverage is reported separately.'
    f
  }
  list(fit=fit)
}
