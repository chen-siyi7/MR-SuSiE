generate_prospective <- function(j) {
  set.seed(j$seed);m<-12L;K<-as.integer(j$M);R<-j$rho^abs(outer(seq_len(m),seq_len(m),'-'));L<-t(chol(R));Sg<-matrix(.4,K,K);diag(Sg)<-1
  Gamma<-matrix(runif(m*K,if(as.logical(j$weak)).005 else .04,if(as.logical(j$weak)).015 else .12),m,K)
  theta<-if(as.logical(j$null))rep(0,K)else if(K==1).2 else c(.2,.1,0)
  invalid<-sort(sample.int(m,3));alpha<-numeric(m);alpha[invalid]<-runif(3,.04,.10)
  B<-R%*%Gamma+.01*L%*%matrix(rnorm(m*K),m,K)%*%chol(Sg)
  mean_y<-as.numeric(R%*%(Gamma%*%theta+alpha));y<-mean_y+.01*as.numeric(L%*%rnorm(m));ysel<-mean_y+.01*as.numeric(L%*%rnorm(m))
  ref_sample<-matrix(rnorm(503*m),503,m)%*%chol(R);Rref<-cor(ref_sample)
  dat<-list(m=m,num_exposures=K,betaZX=B,seZX=matrix(.01,m,K),betaZY=y,se_betaZY=rep(.01,m),corr_in=R,corr_ref=R,sigma_exp=Sg,overlap=0,n1=10000L,n2=10000L)

 list(data=dat,selection_beta=ysel,selection_se=rep(.01,m),reference_503=Rref,truth=theta,invalid=invalid,seed=j$seed)
}
