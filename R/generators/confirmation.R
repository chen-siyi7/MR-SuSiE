generate_binary <- function(seed,root) {
 source(file.path(root,'R/load_native.R'),local=TRUE)
 mr_susie_load_native(file.path(root,'src/binary_gwas.cpp'),env=environment(),cacheDir=file.path(root,'runtime/binary_gwas'))
  set.seed(seed);m<-20L;K<-3L;n1<-10000L;n2<-20000L
  Gamma<-matrix(runif(m*K,.04,.12),m,K);theta<-c(.2,.1,0);invalid<-sort(sample.int(m,4));direct<-numeric(m);direct[invalid]<-runif(4,.04,.10)
  Sigma<-matrix(.4,K,K);diag(Sigma)<-1
  G1<-(matrix(rbinom(n1*m,2,.3),n1,m)-.6)/sqrt(.42);G2<-(matrix(rbinom(n2*m,2,.3),n2,m)-.6)/sqrt(.42)
  X1<-G1%*%Gamma+matrix(rnorm(n1*K),n1,K)%*%chol(Sigma)
  X2<-G2%*%Gamma+matrix(rnorm(n2*K),n2,K)%*%chol(Sigma)
  G1c<-scale(G1,center=TRUE,scale=FALSE);X1c<-scale(X1,center=TRUE,scale=FALSE);xx<-colSums(G1c^2)
  B<-crossprod(G1c,X1c)/xx;residss<-matrix(rep(colSums(X1c^2),each=m),m,K)-B^2*xx
  sx<-sqrt(residss/(n1-2)/xx)
  eta<--2+as.numeric(X2%*%theta+G2%*%direct);Y<-rbinom(n2,1,plogis(eta));out<-single_snp_logistic(G2,Y)
  stopifnot(all(out$converged));neff<-4*sum(Y)*sum(1-Y)/n2
  d<-list(m=m,num_exposures=K,betaZX=B,seZX=sx,betaZY=out$beta,se_betaZY=out$se,
    corr_in=diag(m),corr_ref=diag(m),sigma_exp=cor(X1-G1%*%Gamma),overlap=0,n1=n1,n2=n2)
  list(generated=list(dat=d,select=out$beta,theta=theta,truth_invalid=invalid,seed=seed),n_select=neff,n_exposure=n1,n_estimate=n2,n_reference=Inf,
    selector_family='binary_score',truth=theta,invalid=invalid,cases=sum(Y),controls=sum(1-Y),
    scope='Individual-level structural logistic model; summary-data linear MR remains an approximation.')
}
