profile_prepare <- function(dat, invalid_snps=integer(), which_R=c("ref","in")) {
  which_R<-match.arg(which_R);k<-dat$num_exposures
  if(!is.null(dat$overlap)&&is.finite(dat$overlap)&&dat$overlap>0)stop("Overlap unsupported")
  Rfull<-if(which_R=="ref")dat$corr_ref else dat$corr_in
  gsd<-sqrt(diag(Rfull));stopifnot(all(is.finite(gsd)),all(gsd>0))
  R<-cov2cor(Rfull);I<-sort(unique(as.integer(invalid_snps)))
  stopifnot(all(I>=1L&I<=nrow(R)))
  Sg<-dat$sigma_exp;if(is.null(Sg))stop("Exposure error correlation is required")
  list(B=sweep(dat$betaZX,1,gsd,"*"),y=dat$betaZY*gsd,sy=dat$se_betaZY*gsd,
       sX=sweep(dat$seZX,1,gsd,"*"),R=R,Sg=(Sg+t(Sg))/2,
       Z=R[,I,drop=FALSE],invalid=I,k=k,m=nrow(R))
}

profile_fit_unscaled <- function(B,y,sy,sX,R,Sg,Z=NULL,weight_mode=c("iterated","outcome"),
                                 bread=c("expected","observed"),maxit=100L,tol=1e-8) {
  weight_mode<-match.arg(weight_mode);bread<-match.arg(bread);k<-ncol(B)
  W0<-solve(diag(sy)%*%R%*%diag(sy))
  fixed_W<-if(weight_mode=="outcome")W0 else NULL
  state<-function(beta)profile_score_state(B,y,sy,sX,R,Sg,Z,beta,fixed_W)
  initial<-profile_score_state(B,y,sy,sX,R,Sg,Z,rep(0,k),W0)
  S<-crossprod(B,initial$M%*%B)
  if(qr(S,tol=1e-10)$rank<k)stop("Causal design rank deficiency")
  beta<-as.numeric(solve(S,crossprod(B,initial$M%*%y)))
  merit<-function(st)sqrt(sum(st$score^2))/max(1,sqrt(sum(crossprod(B,st$M%*%y)^2)))
  history<-list();converged<-FALSE;status<-"iteration_limit";st<-state(beta)
  for(it in seq_len(maxit)) {
    old_merit<-merit(st)
    if(old_merit<tol){converged<-TRUE;status<-"converged";break}
    ev<-eigen((st$A+t(st$A))/2,symmetric=TRUE,only.values=TRUE)$values
    if(min(ev)<=1e-10*max(1,max(ev))){status<-"nonpositive_corrected_information";break}
    rhs<-as.numeric(crossprod(B,st$M%*%y))
    target<-as.numeric(solve(st$A,rhs));accepted<-FALSE
    for(step in 2^-(0:10)) {
      proposal<-beta+step*(target-beta)
      next_st<-tryCatch(state(proposal),error=function(e)NULL)
      if(!is.null(next_st)&&all(is.finite(next_st$score))&&
         (merit(next_st)<old_merit*(1-1e-4*step)||merit(next_st)<tol)) {
        history[[length(history)+1L]]<-c(iteration=it,step=step,merit=merit(next_st),
          max_beta_change=max(abs(proposal-beta)))
        beta<-proposal;st<-next_st;accepted<-TRUE;break
      }
    }
    if(!accepted){status<-"line_search_stalled";break}
  }
  if(merit(st)<tol){converged<-TRUE;status<-"converged"}
  ev<-eigen((st$A+t(st$A))/2,symmetric=TRUE,only.values=TRUE)$values
  strength_ok<-min(ev)>1e-10*max(1,max(ev))
  inf<-profile_score_inference(B,y,sy,sX,R,Sg,Z,beta,fixed_W,bread)
  identified<-converged&&strength_ok&&isTRUE(inf$available)
  alpha<-if(!is.null(Z)&&ncol(Z))as.numeric(solve(crossprod(Z,st$W%*%Z),crossprod(Z,st$W%*%(y-B%*%beta))))else numeric()
  list(beta=beta,se=if(identified)inf$se else rep(NA_real_,k),
       p=if(identified)inf$p else rep(NA_real_,k),identified=identified,
       converged=converged,status=status,iterations=length(history),
       normalized_score=merit(st),score_norm=inf$score_norm,min_corrected_eigen=min(ev),
       corrected_condition=max(ev)/min(ev),direct=alpha,
       weight_mode=weight_mode,bread=bread,history=history,inference=inf,
       scope="Corrected-score working covariance conditional on selected support and supplied LD")
}

fit_profile_corrected <- function(dat,invalid_snps=integer(),which_R=c("ref","in"),
                                  weight_mode=c("iterated","outcome"),
                                  bread=c("expected","observed"),maxit=100L,tol=1e-8) {
  which_R<-match.arg(which_R);weight_mode<-match.arg(weight_mode);bread<-match.arg(bread)
  ans<-tryCatch({
    a<-profile_prepare(dat,invalid_snps,which_R)
    f<-profile_fit_unscaled(a$B,a$y,a$sy,a$sX,a$R,a$Sg,a$Z,weight_mode,bread,maxit,tol)
    f$invalid_adjusted<-a$invalid;f$which_R<-which_R;f
  },error=function(e)list(beta=rep(NA_real_,dat$num_exposures),se=rep(NA_real_,dat$num_exposures),
      p=rep(NA_real_,dat$num_exposures),identified=FALSE,converged=FALSE,status=conditionMessage(e)))
  ans
}

fit_fun <- function(dat,invalid_snps=integer())
  fit_profile_corrected(dat,invalid_snps,which_R="ref",weight_mode="iterated",bread="expected")
