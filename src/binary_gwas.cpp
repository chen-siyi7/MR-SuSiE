
#include <Rcpp.h>
using namespace Rcpp;

List single_snp_logistic(NumericMatrix G, NumericVector y) {
 const int n=G.nrow(),m=G.ncol();NumericVector beta(m),se(m);IntegerVector iterations(m);LogicalVector converged(m);
 double py=mean(y);
 for(int j=0;j<m;++j) {
  double a=std::log(py/(1-py)),b=0;bool ok=false;int it=0;double h00=0,h01=0,h11=0;
  for(it=0;it<50;++it) {
   double u0=0,u1=0;h00=h01=h11=0;
   for(int i=0;i<n;++i) {double eta=a+b*G(i,j);double p=eta>35?1-1e-15:eta< -35?1e-15:1/(1+std::exp(-eta));double w=p*(1-p),r=y[i]-p,x=G(i,j);u0+=r;u1+=r*x;h00+=w;h01+=w*x;h11+=w*x*x;}
   double det=h00*h11-h01*h01;if(!(det>0))break;
   double da=(h11*u0-h01*u1)/det,db=(h00*u1-h01*u0)/det;a+=da;b+=db;
   if(std::max(std::abs(da),std::abs(db))<1e-10){ok=true;break;}
  }
  double det=h00*h11-h01*h01;beta[j]=ok?b:NA_REAL;se[j]=ok?std::sqrt(h00/det):NA_REAL;converged[j]=ok;iterations[j]=it+1;
 }
 return List::create(_["beta"]=beta,_["se"]=se,_["converged"]=converged,_["iterations"]=iterations);
}

RCPP_MODULE(mr_susie_binary_gwas) {
  Rcpp::function("single_snp_logistic", &single_snp_logistic);
}
