


#include <RcppArmadillo.h>

arma::mat native_profile_quadratic(const arma::mat& M,const arma::mat& R,const arma::mat& sx,const arma::mat& Sg,const arma::vec& beta) {
  unsigned int k=sx.n_cols;arma::mat Q(k,k,arma::fill::zeros);
  if(arma::accu(arma::abs(sx))==0 || arma::accu(arma::abs(beta))==0)return Q;
  arma::mat a=sx*arma::diagmat(beta)*Sg;
  std::vector<arma::mat> MT;MT.reserve(k);
  for(unsigned int j=0;j<k;++j){arma::mat T=R;T.each_col()%=a.col(j);T.each_row()%=sx.col(j).t();MT.push_back(M*T);}
  for(unsigned int i=0;i<k;++i)for(unsigned int j=0;j<=i;++j){Q(i,j)=arma::accu(MT[j]%MT[i].t());Q(j,i)=Q(i,j);}
  return Q;
}

RCPP_MODULE(mr_susie_profile) {
  Rcpp::function("native_profile_quadratic", &native_profile_quadratic);
}
