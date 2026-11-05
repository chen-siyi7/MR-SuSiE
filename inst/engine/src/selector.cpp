


#pragma STDC FP_CONTRACT OFF
#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;
struct MixtureLikelihood {
  std::vector<double> bh, s2, logprior, base;
  mutable std::vector<double> trace_x,trace_y; bool tracing=false;
  explicit MixtureLikelihood(NumericVector b, NumericVector s, NumericVector p):bh(b.begin(),b.end()),s2(s.begin(),s.end()),logprior(b.size()),base(b.size()) {
    for(size_t j=0;j<bh.size();++j) {logprior[j]=std::log(p[j]+std::sqrt(DBL_EPSILON));base[j]=R::dnorm(bh[j],0,std::sqrt(s2[j]),true);}
  }
  double value(double v) const {
    std::vector<double> lp(bh.size());double mx=R_NegInf;
    for(size_t j=0;j<bh.size();++j) {lp[j]=std::isinf(s2[j])?0:R::dnorm(bh[j],0,std::sqrt(v+s2[j]),true)-base[j]+logprior[j];if(ISNAN(lp[j]))return NA_REAL;mx=std::max(mx,lp[j]);}
    long double su=0;for(double x:lp)su+=std::exp(x-mx);return std::log(static_cast<double>(su))+mx;
  }
  double objective(double logv) const {double ans=-value(std::exp(logv));if(tracing){trace_x.push_back(logv);trace_y.push_back(ans);}return ans;}
};

List native_single_effect(NumericVector Xty,NumericVector dXtX,double V,double residual_variance,NumericVector prior_weights,std::string optimize_V,double check_null_threshold,Function optimize_variance) {
  int p=Xty.size();if(dXtX.size()!=p||prior_weights.size()!=p)stop("length mismatch");
  NumericVector bh(p),s2(p);for(int j=0;j<p;++j){bh[j]=(1/dXtX[j])*Xty[j];s2[j]=residual_variance/dXtX[j];}
  MixtureLikelihood ll(bh,s2,prior_weights);
  if(optimize_V=="optim") {
    V=as<double>(optimize_variance(bh,s2,prior_weights,V,check_null_threshold));
  } else if(optimize_V!="none")stop("Native prototype supports optim and none only");
  NumericVector lbf(p),lpo(p),alpha(p),mu(p),mu2(p);double mx=R_NegInf;
  for(int j=0;j<p;++j){lbf[j]=R::dnorm(bh[j],0,std::sqrt(V+s2[j]),true)-ll.base[j];lpo[j]=lbf[j]+ll.logprior[j];if(std::isinf(s2[j])){lbf[j]=0;lpo[j]=0;}mx=std::max(mx,lpo[j]);}
  long double su=0;for(int j=0;j<p;++j){alpha[j]=std::exp(lpo[j]-mx);su+=alpha[j];}
  double sum=static_cast<double>(su);
  for(int j=0;j<p;++j){alpha[j]/=sum;double var=R_pow(1/V+dXtX[j]/residual_variance,-1.0);mu[j]=(1/residual_variance)*var*Xty[j];mu2[j]=var+mu[j]*mu[j];}
  return List::create(_["alpha"]=alpha,_["mu"]=mu,_["mu2"]=mu2,_["lbf"]=lbf,_["V"]=V,_["lbf_model"]=mx+std::log(sum));
}

List native_update_ss(NumericMatrix XtX,NumericVector Xty,List input,bool estimate_prior_variance,std::string estimate_prior_method,double check_null_threshold,Function optimize_variance) {
  List s=clone(input);
  NumericMatrix alpha=s["alpha"],mu=s["mu"],mu2=s["mu2"],lbfvar=s["lbf_variable"];
  NumericVector V=s["V"],lbf=s["lbf"],KL=s["KL"],pi=s["pi"],dXtX=XtX.attr("d");
  NumericVector raw=s["XtXr"];NumericVector fitted=clone(raw);
  const int L=alpha.nrow(),p=alpha.ncol();double sigma2=as<double>(s["sigma2"]);
  if(XtX.nrow()!=p||XtX.ncol()!=p||Xty.size()!=p)stop("dimension mismatch");
  std::string method=estimate_prior_variance?estimate_prior_method:"none";
  Function matprod("%*%");
  for(int l=0;l<L;++l){
    NumericVector eb(p);for(int j=0;j<p;++j)eb[j]=alpha(l,j)*mu(l,j);
    NumericVector product=matprod(XtX,eb);NumericVector residual(p);
    for(int j=0;j<p;++j){fitted[j]-=product[j];residual[j]=Xty[j]-fitted[j];}
    List post=native_single_effect(residual,dXtX,V[l],sigma2,pi,method,check_null_threshold,optimize_variance);
    NumericVector pa=post["alpha"],pm=post["mu"],pm2=post["mu2"],plbf=post["lbf"];
    V[l]=as<double>(post["V"]);lbf[l]=as<double>(post["lbf_model"]);
    long double cross1=0,cross2=0;
    for(int j=0;j<p;++j){alpha(l,j)=pa[j];mu(l,j)=pm[j];mu2(l,j)=pm2[j];lbfvar(l,j)=plbf[j];eb[j]=pa[j]*pm[j];cross1+=eb[j]*residual[j];cross2+=dXtX[j]*(pa[j]*pm2[j]);}
    double elog=-0.5/sigma2*(-2*static_cast<double>(cross1)+static_cast<double>(cross2));KL[l]=-lbf[l]+elog;
    NumericVector newproduct=matprod(XtX,eb);for(int j=0;j<p;++j)fitted[j]+=newproduct[j];
  }
  fitted.attr("dim")=IntegerVector::create(p,1);fitted.attr("dimnames")=R_NilValue;
  s["alpha"]=alpha;s["mu"]=mu;s["mu2"]=mu2;s["V"]=V;s["lbf"]=lbf;s["lbf_variable"]=lbfvar;s["KL"]=KL;s["XtXr"]=fitted;
  return s;
}


SEXP native_likelihood_context(NumericVector bh, NumericVector s2, NumericVector prior) {
  return XPtr<MixtureLikelihood>(new MixtureLikelihood(bh,s2,prior),true);
}

double native_cached_likelihood(SEXP context,double v) {
  XPtr<MixtureLikelihood> ll(context);return ll->value(v);
}

RCPP_MODULE(mr_susie_selector) {
  Rcpp::function("native_single_effect", &native_single_effect);
  Rcpp::function("native_update_ss", &native_update_ss);
  Rcpp::function("native_likelihood_context", &native_likelihood_context);
  Rcpp::function("native_cached_likelihood", &native_cached_likelihood);
}
