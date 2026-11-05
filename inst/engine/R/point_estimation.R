posterior_mean_estimate <- function(v,rule='EB',expected_k=NULL) {
 rule<-match.arg(rule,'EB');k<-if(is.null(expected_k))length(v$beta)else as.integer(expected_k)
 fail<-function(reason)list(beta=rep(NA_real_,k),identified=FALSE,status=reason)
 if(length(v$beta)!=k)return(fail('score_coefficients_unavailable'))
 if(!isTRUE(v$identified)||!isTRUE(v$converged)||!isTRUE(v$selector_converged))return(fail('score_unavailable'))
 V<-v$working_phi_covariance
 if(!identical(dim(V),c(as.integer(k),as.integer(k)))||any(!is.finite(V)))return(fail('full_working_covariance_unavailable'))
 ev<-eigen((V+t(V))/2,symmetric=TRUE);d<-ev$values;Q<-ev$vectors
 if(min(d)<=1e-12*max(d))return(fail('working_covariance_not_positive_definite'))
 u<-as.numeric(crossprod(Q,v$beta));scale<-median(d)
  tau2<-scale*10^seq(-4,4,by=.125)
  denom<-outer(d,tau2,'+');criterion_all<-colSums(log(denom)+u^2/denom)
  i<-which.min(criterion_all);shrink<-tau2[i]/denom[,i];lambda<-1/tau2[i];criterion<-criterion_all[i]
 point<-as.numeric(Q%*%(shrink*u))
 list(beta=point,inference_beta=v$beta,se=v$se,p=v$p,critical_95=v$critical_95,
  lower=v$beta-v$critical_95*v$se,upper=v$beta+v$critical_95*v$se,
  identified=TRUE,converged=TRUE,status='point_available_score_inference_separate',
  rule=rule,lambda=lambda,shrinkage=shrink,criterion=criterion,
  covariance_source='Corrected-score working covariance',
  inference_scope='Confidence intervals and tests target the full coefficient vector using the original corrected score; they are not recentered on the posterior-mean point estimate.',
  selection=v$selection)
}
