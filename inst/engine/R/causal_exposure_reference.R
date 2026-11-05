make_support_reference <- function(source_root, point_rule, cache, native=TRUE) {
  base<-source_root
  if(native) {
    env<-new.env(parent=environment())
    sys.source(file.path(base,'R/profile_engine.R'),env)
    factory<-env$make_native_profile_engine
    code<-paste(deparse(factory),collapse='\n')
    needle<-"file.path(d, \"cppcache_profile\")"
    stopifnot(grepl(needle,code,fixed=TRUE))
    local<-new.env(parent=environment(factory));local$oracle_cache<-cache
    f<-eval(parse(text=sub(needle,'oracle_cache',code,fixed=TRUE)),local)
    engine<-f(source_root)
  }else {
    engine<-new.env(parent=environment())
    for(f in c('R/profile_score.R','R/profile_solver.R','R/score_refit.R'))
      sys.source(file.path(base,f),engine)
  }
  sys.source(point_rule,engine)
  fit<-function(dat,causal_exposures,invalid_variants,retained_rows=seq_len(dat$m),which_ld='reference') {
    K<-as.integer(dat$num_exposures);E<-sort(unique(as.integer(causal_exposures)))
    I<-sort(unique(as.integer(invalid_variants)));ii<-as.integer(retained_rows)
    stopifnot(all(E>=1L&E<=K),all(I>=1L&I<=dat$m),all(ii>=1L&ii<=dat$m),!anyDuplicated(ii))
    empty<-function(status)list(identified=FALSE,status=status,causal=E,invalid=I,rows=ii,score=NULL,point=NULL)
    if(!is.null(dat$overlap)&&is.finite(dat$overlap)&&dat$overlap>0)return(empty('overlap_unsupported_by_frozen_score'))
    if(!length(E))return(list(identified=TRUE,status='all_coefficients_known_zero',causal=E,invalid=I,rows=ii,score=NULL,point=NULL))
    Rpop<-if(which_ld=='population')dat$corr_in else dat$corr_ref
    stopifnot(identical(dim(Rpop),c(dat$m,dat$m)),all(is.finite(Rpop)),all(diag(Rpop)>0))
    gsd<-sqrt(diag(Rpop));Rfull<-cov2cor(Rpop);R<-Rfull[ii,ii,drop=FALSE]
    B<-sweep(dat$betaZX[ii,E,drop=FALSE],1,gsd[ii],'*')
    sX<-sweep(dat$seZX[ii,E,drop=FALSE],1,gsd[ii],'*')
    y<-dat$betaZY[ii]*gsd[ii];sy<-dat$se_betaZY[ii]*gsd[ii]
    Sg<-dat$sigma_exp[E,E,drop=FALSE];Sg<-(Sg+t(Sg))/2

    Z<-Rfull[ii,I,drop=FALSE]
    if(ncol(Z)) {
      norms<-sqrt(colSums(Z^2));good<-which(norms>0)
      if(length(good)) {
        zz<-qr(sweep(Z[,good,drop=FALSE],2,norms[good],'/'),tol=1e-10)
        Z<-Z[,good[zz$pivot[seq_len(zz$rank)]],drop=FALSE]
      }else Z<-matrix(0,length(ii),0)
    }
    H<-cbind(B,Z);rank_H<-qr(H,tol=1e-8)$rank
    if(rank_H-ncol(Z)<length(E))return(empty('true_structure_residual_causal_rank_deficient'))
    df<-nrow(H)-rank_H
    if(df<=0)return(empty('true_structure_no_residual_degrees_of_freedom'))
    v<-engine$profile_fit_unscaled(B,y,sy,sX,R,Sg,Z,weight_mode='iterated',bread='expected')
    if(!isTRUE(v$identified))return(c(empty(v$status),list(diagnostics=v)))
    st<-v$inference
    residual<-y-as.numeric(B%*%v$beta)-if(ncol(Z))as.numeric(Z%*%v$direct)else 0
    phi<-max(1,as.numeric(crossprod(residual,st$W%*%residual))/df)
    linear<-crossprod(B,st$M%*%st$Omega%*%st$M%*%B);Ai<-solve(st$A)
    V<-Ai%*%(phi*linear+st$quadratic_term)%*%t(Ai);V<-(V+t(V))/2
    v$working_phi_covariance<-V;v$se<-sqrt(diag(V));v$p<-2*pt(-abs(v$beta/v$se),df)
    v$critical_95<-qt(.975,df);v$residual_df<-df;v$working_phi<-phi
    v$identified<-all(is.finite(v$se))&&all(v$se>0)
    v$selector_converged<-v$selector_usable<-TRUE
    v$selection<-list(invalid=I,exposures=E,converged=TRUE,oracle_structure=TRUE)
    point<-engine$posterior_mean_estimate(v,'EB',length(E))
    list(identified=isTRUE(point$identified),status=if(isTRUE(point$identified))'oracle_true_structure'else point$status,
      causal=E,invalid=I,rows=ii,score=v,point=point,nuisance_rank=ncol(Z),residual_df=df)
  }
  list(fit=fit,profile=engine)
}
