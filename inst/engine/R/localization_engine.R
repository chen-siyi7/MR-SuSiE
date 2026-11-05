make_localization_selector <- function(root, original_joint_function, return_posterior=FALSE, native_updates=TRUE) {
  stopifnot(as.character(packageVersion('susieR'))=='0.14.2')
  d<-file.path(root,'src')
  ns<-asNamespace('susieR'); local_ns<-new.env(parent=ns)
  for(n in ls(ns,all.names=TRUE)) {f<-get(n,ns);if(is.function(f)&&identical(environment(f),ns)){environment(f)<-local_ns;assign(n,f,local_ns)}}
  mr_susie_load_native(file.path(d,'selector.cpp'),env=local_ns,cacheDir=file.path(d,'cppcache_replay_ser'),showOutput=FALSE,verbose=FALSE)
  opt<-function(betahat,shat2,prior_weights,V_init,check_null_threshold) {
    context<-native_likelihood_context(betahat,shat2,prior_weights)
    objective<-function(logv)-native_cached_likelihood(context,exp(logv))
    lv<-stats::optim(par=log(max(c(betahat^2-shat2,1),na.rm=TRUE)),fn=objective,
                     method='Brent',lower=-30,upper=15)$par
    if(objective(lv)>objective(log(V_init)))lv<-log(V_init)
    V<-exp(lv)
    if(native_cached_likelihood(context,0)+check_null_threshold>=native_cached_likelihood(context,V))V<-0
    V
  }
  environment(opt)<-local_ns;local_ns$native_optimize_variance<-opt
  local_ns$original_ser<-local_ns$single_effect_regression_ss
  ser<-function(Xty,dXtX,V=1,residual_variance=1,prior_weights=NULL,optimize_V=c('none','optim','uniroot','EM','simple'),check_null_threshold=0) {
    optimize_V<-match.arg(optimize_V)
    if(is.null(prior_weights))prior_weights<-rep(1/length(dXtX),length(dXtX))
    if(!optimize_V%in%c('none','optim')||any(dXtX<=0)||any(!is.finite(c(Xty,dXtX))))return(original_ser(Xty,dXtX,V,residual_variance,prior_weights,optimize_V,check_null_threshold))
    native_single_effect(as.numeric(Xty),dXtX,V,residual_variance,prior_weights,optimize_V,check_null_threshold,native_optimize_variance)
  }
  environment(ser)<-local_ns;local_ns$single_effect_regression_ss<-ser
  if(native_updates) {
    local_ns$original_update<-local_ns$update_each_effect_ss
    update<-function(XtX,Xty,s_init,estimate_prior_variance=FALSE,estimate_prior_method='optim',check_null_threshold=0) {
      if(!estimate_prior_method%in%c('none','optim')||any(attr(XtX,'d')<=0)||any(!is.finite(Xty)))return(original_update(XtX,Xty,s_init,estimate_prior_variance,estimate_prior_method,check_null_threshold))
      native_update_ss(XtX,as.numeric(Xty),s_init,estimate_prior_variance,estimate_prior_method,check_null_threshold,native_optimize_variance)
    }
    environment(update)<-local_ns;local_ns$update_each_effect_ss<-update
  }
  wrapper_env<-new.env(parent=environment(original_joint_function));wrapper_env$fast_susie_suff_stat<-local_ns$susie_suff_stat
  txt<-paste(deparse(original_joint_function),collapse='\n')
  stopifnot(length(gregexpr('susieR::susie_suff_stat',txt,fixed=TRUE)[[1]])==1L)
  txt<-sub('susieR::susie_suff_stat','fast_susie_suff_stat',txt,fixed=TRUE)
  txt<-sub('check_input = TRUE','check_input = TRUE, max_iter = 1000, refine = TRUE',txt,fixed=TRUE)
  f<-eval(parse(text=txt),wrapper_env);environment(f)<-wrapper_env
  if(return_posterior){b<-body(f);last<-length(b);b[[last]]<-substitute({answer<-OLD;answer$susie_fit<-res;answer},list(OLD=b[[last]]));body(f)<-b}
  list(select=f,ser=ser,namespace=local_ns)
}
