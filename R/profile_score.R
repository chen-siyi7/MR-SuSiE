ps_trace <- function(x) sum(diag(x))
ps_inverse <- function(x,tol=1e-10) {
 s<-svd(x);if(!length(s$d)||min(s$d)<=tol*max(s$d))return(NULL)
 s$v%*%(t(s$u)/s$d)
}
profile_score_state <- function(B,y,sy,sX,R,Sg,Z,beta,fixed_W=NULL) {
 k<-ncol(B);m<-nrow(B);if(is.null(Z))Z<-matrix(0,m,0)
 Vtheta<-sweep(sX,2,beta,'*');Omega<-diag(sy)%*%R%*%diag(sy)+R*(Vtheta%*%Sg%*%t(Vtheta))
 W<-if(is.null(fixed_W))solve(Omega)else fixed_W
 if(ncol(Z)) {zz<-ps_inverse(crossprod(Z,W%*%Z));if(is.null(zz))stop('Nuisance matrix rank deficiency');M<-W-W%*%Z%*%zz%*%t(Z)%*%W} else M<-W
 M<-(M+t(M))/2;C<-Sg*(t(sX)%*%(M*R)%*%sX);C<-(C+t(C))/2
 A<-crossprod(B,M%*%B)-C
 score<-as.numeric(crossprod(B,M%*%(y-as.numeric(B%*%beta)))+C%*%beta)
 list(W=W,M=M,C=C,A=A,Omega=Omega,score=score)
}
profile_score_inference <- function(B,y,sy,sX,R,Sg,Z=NULL,beta,fixed_W=NULL,bread=c('observed','expected')) {
 bread<-match.arg(bread);k<-ncol(B);m<-nrow(B);st<-profile_score_state(B,y,sy,sX,R,Sg,Z,beta,fixed_W)

 TT<-lapply(seq_len(k),function(j) {
  ans<-matrix(0,m,m)
  for(l in seq_len(k))ans<-ans+beta[l]*Sg[l,j]*(sX[,l]*R)*rep(sX[,j],each=m)
  ans
 })
 S<-crossprod(B,st$M%*%st$Omega%*%st$M%*%B)
 Q<-matrix(0,k,k)
 for(i in seq_len(k))for(j in seq_len(k))Q[i,j]<-ps_trace(st$M%*%TT[[j]]%*%st$M%*%TT[[i]])
 S<-(S+Q+t(S+Q))/2
 J<-matrix(NA_real_,k,k)
 for(j in seq_len(k)) {
  h<-1e-5*max(1,abs(beta[j]));lo<-hi<-beta;lo[j]<-lo[j]-h;hi[j]<-hi[j]+h
  J[,j]<--(profile_score_state(B,y,sy,sX,R,Sg,Z,hi,fixed_W)$score-profile_score_state(B,y,sy,sX,R,Sg,Z,lo,fixed_W)$score)/(2*h)
 }
 makeV<-function(JJ){ii<-ps_inverse(JJ);if(is.null(ii))return(matrix(NA_real_,k,k));v<-ii%*%S%*%t(ii);(v+t(v))/2}
 ve<-makeV(st$A);vo<-makeV(J);vv<-if(bread=='observed')vo else ve
 ok<-all(is.finite(vv))&&all(diag(vv)>0)&&min(eigen(S,symmetric=TRUE,only.values=TRUE)$values)>-1e-8*max(1,max(abs(S)))
 se<-if(ok)sqrt(diag(vv))else rep(NA_real_,k)
 c(st,list(beta=beta,se=se,p=2*pnorm(-abs(beta/se)),covariance=vv,
  covariance_expected=ve,covariance_observed=vo,score_covariance=S,quadratic_term=Q,
  observed_bread=J,bread=bread,available=ok,score_norm=sqrt(sum(st$score^2)),
  scope='Fixed correct nuisance support, known LD, Gaussian sampling errors; no selection or external-LD uncertainty correction'))
}
