make_native_profile_engine <- function(root) {
  d<-root
  e<-new.env(parent=environment())
  sys.source(file.path(d,'R/profile_score.R'),envir=e)
  sys.source(file.path(d,'R/profile_solver.R'),envir=e)
  sys.source(file.path(d,'R/score_refit.R'),envir=e)
  mr_susie_load_native(file.path(d,'src/profile.cpp'),env=e,cacheDir=file.path(d,'cppcache_profile'),showOutput=FALSE,verbose=FALSE)
  e$original_profile_inference<-e$profile_score_inference
  f<-function(B,y,sy,sX,R,Sg,Z=NULL,beta,fixed_W=NULL,bread=c('observed','expected')) {
    bread<-match.arg(bread)
    if(bread=='observed')return(original_profile_inference(B,y,sy,sX,R,Sg,Z,beta,fixed_W,bread))
    k<-ncol(B);st<-profile_score_state(B,y,sy,sX,R,Sg,Z,beta,fixed_W)
    Q<-native_profile_quadratic(st$M,R,sX,Sg,beta)
    S<-crossprod(B,st$M%*%st$Omega%*%st$M%*%B)+Q;S<-(S+t(S))/2
    Ai<-ps_inverse(st$A)
    vv<-if(is.null(Ai))matrix(NA_real_,k,k)else{z<-Ai%*%S%*%t(Ai);(z+t(z))/2}
    ok<-all(is.finite(vv))&&all(diag(vv)>0)&&min(eigen(S,symmetric=TRUE,only.values=TRUE)$values)>-1e-8*max(1,max(abs(S)))
    se<-if(ok)sqrt(diag(vv))else rep(NA_real_,k)
    c(st,list(beta=beta,se=se,p=2*pnorm(-abs(beta/se)),covariance=vv,covariance_expected=vv,
      covariance_observed=matrix(NA_real_,k,k),score_covariance=S,quadratic_term=Q,
      observed_bread=matrix(NA_real_,k,k),observed_bread_computed=FALSE,bread=bread,available=ok,
      score_norm=sqrt(sum(st$score^2)),scope='Expected-bread fixed-support Gaussian inference; equivalent native quadratic contraction. Observed-bread diagnostic not requested.'))
  }
  environment(f)<-e;e$profile_score_inference<-f;e
}
