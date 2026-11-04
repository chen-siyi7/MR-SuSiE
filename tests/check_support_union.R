a <- grep('^--file=',commandArgs(),value=TRUE)
root <- normalizePath(file.path(dirname(gsub('~+~',' ',sub('^--file=','',a[1]),fixed=TRUE)),'..'))
source(file.path(root,'R/support_union.R'))
print(support_union_self_check())
stopifnot(ar_chi2_sf(0,3)==1,ar_chi2_sf(Inf,3)==0)
for(df in c(1,2,3,9,12,30,100)) {
  p <- c(.001,.01,.5,.95,.99,.999)
  stopifnot(max(abs(ar_chi2_sf(ar_chi2_ppf(p,df),df)-(1-p)))<1e-12)
}
x <- list(B=matrix(c(.1,.12,.07,.08),ncol=1),y=c(.02,.024,.014,.016),
  R=diag(4),Sigma=matrix(1),smax=0L)
z <- do.call(union_ar,x)
stopifnot(z$supports==1L,ar_contains(z$sets[[1]],.2))
bad <- x
bad$R[1,2] <- .2
stopifnot(inherits(try(do.call(union_ar,bad),silent=TRUE),'try-error'))
bad <- x
bad$smax <- 1.5
stopifnot(inherits(try(do.call(union_ar,bad),silent=TRUE),'try-error'))
m <- 8L
R <- .97^abs(outer(seq_len(m),seq_len(m),'-'))
B <- matrix(seq(.03,.10,length.out=m),ncol=1)
y <- c(10,.003,-.007,.014,-.011,.009,-.005,.018)
conditional <- vapply(seq_len(m),function(t) {
  s <- setdiff(seq_len(m),t)
  e <- y[t]-R[t,s,drop=FALSE]%*%solve(R[s,s],y[s])
  variance <- R[t,t]-R[t,s,drop=FALSE]%*%solve(R[s,s],R[s,t,drop=FALSE])
  as.numeric(e^2/variance)
},numeric(1))
fit <- union_ar(B,y,R,matrix(1),smax=m-1L)
expected <- stats::pchisq(min(conditional)/.01^2,df=1,lower.tail=FALSE)
stopifnot(abs(fit$joint_null_p-expected)<1e-10)
for(scale in c(1e-200,1e-100,1,1e100,1e200)) {
  bounds <- ar_scalar_set(scale,2*scale,scale)
  stopifnot(isTRUE(all.equal(bounds,list(c(2-sqrt(3),2+sqrt(3))),tolerance=1e-12)))
  stopifnot(!length(ar_scalar_set(0,0,scale)))
}
small_root <- ar_scalar_set(1,1e8,1)[[1]]
stopifnot(abs(small_root[1]-5e-9)<1e-23,abs(prod(small_root)-1)<1e-12)
tiny_roots <- ar_scalar_set(1,1e-200,0)[[1]]
stopifnot(tiny_roots[1]==0,abs(tiny_roots[2]/2e-200-1)<1e-14)
tiny_constant <- -1e-320
bounds <- ar_scalar_set(1,0,tiny_constant)[[1]]
stopifnot(max(abs(bounds/c(-sqrt(-tiny_constant),sqrt(-tiny_constant))-1))<1e-14)
stopifnot(isTRUE(all.equal(ar_scalar_set(-1e-12,0,1),list(c(-Inf,-1e6),c(1e6,Inf)),tolerance=1e-12)))
B <- cbind(c(.1,.03,.12,.05,.08,.02),c(.02,.11,.04,.13,.01,.09))
y <- as.numeric(B%*%c(.2,-.1))+c(.005,-.003,.002,0,-.004,.001)
R <- .3^abs(outer(seq_len(6),seq_len(6),'-'))
Sigma <- matrix(c(1,.4,.4,1),2)
reference <- union_ar(B,y,R,Sigma,sx=.005,sy=.01,smax=1L)
stopifnot(all(is.finite(unlist(reference$sets))))
for(scale in c(1e-150,1,1e150)) {
  fit <- union_ar(B*scale,y*scale,R,Sigma,sx=.005*scale,sy=.01*scale,smax=1L)
  stopifnot(isTRUE(all.equal(fit,reference,tolerance=1e-9)))
}
bad <- x
bad$B <- bad$B+1i
stopifnot(inherits(try(do.call(union_ar,bad),silent=TRUE),'try-error'))
cat('Support-union geometry, chi-square calculations and input checks passed\n')
