ar_chi2_sf <- function(x,df) {
  if(length(df)!=1L||!is.finite(df)||df<1||df!=floor(df)) stop('Positive integer df required')
  stats::pchisq(x,df=df,lower.tail=FALSE)
}

ar_chi2_ppf <- function(p,df) {
  if(length(df)!=1L||!is.finite(df)||df<1||df!=floor(df)) stop('Positive integer df required')
  if(any(!is.finite(p)|p<=0|p>=1)) stop('Probabilities must lie strictly between zero and one')
  stats::qchisq(p,df=df)
}

ar_scalar_set <- function(a,b,c,tol=1e-10) {
  coefficients <- c(a,b,c)
  if(length(coefficients)!=3L||!is.numeric(coefficients)||is.complex(coefficients)||any(!is.finite(coefficients))) stop('Finite real quadratic coefficients required')
  if(length(tol)!=1L||!is.finite(tol)||tol<=0) stop('Positive finite tolerance required')
  scale <- max(abs(coefficients))
  if(scale==0) return(list(c(-Inf,Inf)))
  normalized <- coefficients/scale
  if(any(coefficients!=0&normalized==0)) return(list(c(-Inf,Inf)))
  a <- normalized[1]
  b <- normalized[2]
  c <- normalized[3]
  if(a==0) {
    if(b==0) return(if(c<=0) list(c(-Inf,Inf)) else list())
    boundary <- c/(2*b)
    if(!is.finite(boundary)) return(list(c(-Inf,Inf)))
    return(if(b>0) list(c(boundary,Inf)) else list(c(-Inf,boundary)))
  }
  product_root <- sqrt(abs(a))*sqrt(abs(c))
  if(c!=0&&((a<0)!=(c<0))) {
    disc_scale <- max(abs(b),product_root)
    disc_root <- disc_scale*sqrt((b/disc_scale)^2+(product_root/disc_scale)^2)
  } else {
    gap <- abs(b)-product_root
    if(gap<0&&(gap*(abs(b)+product_root)) < -tol) return(if(a<0) list(c(-Inf,Inf)) else list())
    disc_root <- sqrt(max(0,gap))*sqrt(abs(b)+product_root)
  }
  if(disc_root==0) {
    if(a<0) return(list(c(-Inf,Inf)))
    root <- b/a
    return(if(is.finite(root)) list(c(root,root)) else list(c(-Inf,Inf)))
  }
  q <- b+if(b>=0) disc_root else -disc_root
  roots <- if(q==0) c(0,0) else sort(c(q/a,c/q))
  if(any(!is.finite(roots))) return(list(c(-Inf,Inf)))
  if(a>0) list(roots) else list(c(-Inf,roots[1]),c(roots[2],Inf))
}

ar_merge_sets <- function(sets) {
  if(!length(sets)) return(list())
  bounds <- do.call(rbind,sets)
  bounds <- bounds[order(bounds[,1],bounds[,2]),,drop=FALSE]
  out <- list()
  for(i in seq_len(nrow(bounds))) {
    z <- as.numeric(bounds[i,])
    j <- length(out)
    if(!j||z[1]>out[[j]][2]+1e-10) out[[j+1L]] <- z
    else out[[j]][2] <- max(out[[j]][2],z[2])
  }
  out
}

ar_contains <- function(intervals,t) {
  if(!length(intervals)) return(rep(FALSE,length(t)))
  vapply(t,function(x) any(vapply(intervals,function(z) z[1]-1e-10<=x&&x<=z[2]+1e-10,logical(1))),logical(1))
}

ar_project <- function(H,h,c,k) {
  rest <- setdiff(seq_along(h),k)
  if(!length(rest)) return(ar_scalar_set(H[k,k],h[k],c))
  U <- H[rest,rest,drop=FALSE]
  scale <- max(abs(U))
  if(scale==0||min(eigen(U/scale,symmetric=TRUE,only.values=TRUE)$values)<=1e-10) return(list(c(-Inf,Inf)))
  z <- solve(U,H[rest,k],tol=0)
  v <- solve(U,h[rest],tol=0)
  coefficients <- c(H[k,k]-sum(H[k,rest]*z),h[k]-sum(H[k,rest]*v),c-sum(h[rest]*v))
  if(any(!is.finite(coefficients))) return(list(c(-Inf,Inf)))
  do.call(ar_scalar_set,as.list(coefficients))
}

union_ar <- function(B,y,R,Sigma,sx=.01,sy=.01,smax=3L,level=.95) {
  if(any(vapply(list(B,y,R,Sigma,sx,sy),is.complex,logical(1)))) stop('Real summary inputs required')
  B <- as.matrix(B)
  y <- as.numeric(y)
  R <- as.matrix(R)
  Sigma <- as.matrix(Sigma)
  m <- nrow(B)
  K <- ncol(B)
  if(m<1L||K<1L||length(y)!=m||!identical(dim(R),c(m,m))||!identical(dim(Sigma),c(K,K))) stop('Incompatible summary dimensions')
  if(any(!is.finite(c(B,y,R,Sigma)))) stop('Summary inputs must be finite')
  if(length(smax)!=1L||!is.finite(smax)||smax!=floor(smax)||smax<0||smax>=m) stop('Invalid prespecified support bound')
  if(length(level)!=1L||!is.finite(level)||level<=0||level>=1) stop('Invalid confidence level')
  if(length(sx)!=1L||length(sy)!=1L||!is.finite(sx)||!is.finite(sy)||sx<0||sy<=0) stop('Homogeneous valid standard errors required')
  symmetric <- function(x) max(abs(x-t(x)))<=1e-10*max(1,max(abs(x)))
  if(!symmetric(R)||!symmetric(Sigma)) stop('LD and exposure-error covariance must be symmetric')
  if(inherits(try(chol(R),silent=TRUE),'try-error')) stop('Positive-definite supplied LD required')
  if(min(eigen(Sigma,symmetric=TRUE,only.values=TRUE)$values)< -1e-10*max(1,max(abs(Sigma)))) stop('Exposure-error covariance must be positive semidefinite')
  summary_scale <- max(abs(B),abs(y),sx,sy)
  B <- B/summary_scale
  y <- y/summary_scale
  sx <- sx/summary_scale
  sy <- sy/summary_scale
  if(sy*sy==0||(sx>0&&sx*sx==0)) stop('Summary scales exceed numerical range')
  smax <- as.integer(smax)
  supports <- if(smax) combn(seq_len(m),smax,simplify=FALSE) else list(integer())
  Ri <- solve(R,tol=0)
  critical <- ar_chi2_ppf(level,m-smax)
  summaries <- cbind(B,y)
  V <- sx*sx*Sigma
  pieces <- vector('list',K)
  constants <- numeric(length(supports))
  for(j in seq_along(supports)) {
    ids <- supports[[j]]
    residual <- summaries
    if(length(ids)) {
      residual <- residual-R[,ids,drop=FALSE]%*%solve(R[ids,ids,drop=FALSE],summaries[ids,,drop=FALSE],tol=0)
      residual[ids,] <- 0
    }
    gram <- crossprod(residual,Ri%*%residual)
    A <- gram[seq_len(K),seq_len(K),drop=FALSE]
    h <- gram[seq_len(K),K+1L]
    cc <- gram[K+1L,K+1L]
    constants[j] <- cc
    H <- A-critical*V
    constant <- cc-critical*sy*sy
    for(k in seq_len(K)) pieces[[k]] <- c(pieces[[k]],ar_project(H,h,constant,k))
  }
  list(sets=lapply(pieces,ar_merge_sets),
    joint_null_p=as.numeric(ar_chi2_sf(max(0,min(constants))/(sy*sy),m-smax)),
    supports=length(supports),critical=as.numeric(critical))
}

support_union_self_check <- function() {
  grid <- seq(-20,20,length.out=10001)
  cases <- rbind(c(1,0,-1),c(-1,0,1),c(0,1,2),c(0,0,1),c(0,0,-1),c(1,1,1))
  for(i in seq_len(nrow(cases))) {
    z <- cases[i,]
    stopifnot(identical(ar_contains(ar_scalar_set(z[1],z[2],z[3]),grid),z[1]*grid^2-2*z[2]*grid+z[3]<=1e-9))
  }
  B <- matrix(c(.1,.12,.07,.08),ncol=1)
  y <- as.numeric(B*.2)
  y[1] <- y[1]+10
  fit <- union_ar(B,y,diag(4),matrix(1),smax=1L)
  stopifnot(ar_contains(fit$sets[[1]],.2))
  H <- matrix(c(3,1,1,2),2)
  h <- c(1,2)
  ss <- ar_project(H,h,-1,1L)
  for(t in seq(-5,5,length.out=101)) {
    u <- (h[2]-H[2,1]*t)/H[2,2]
    v <- c(t,u)
    stopifnot(ar_contains(ss,t)==as.numeric(crossprod(v,H%*%v)-2*crossprod(h,v)-1<=1e-9))
  }
  stopifnot(length(ar_scalar_set(1,0,1))==0L,
    length(ar_scalar_set(-1,0,1))==2L,
    identical(ar_project(diag(c(1,0)),c(0,0),1,1L),list(c(-Inf,Inf))))
  list(scalar_quadratic_checks=nrow(cases),projection_checks=101L,
    large_invalid_support_check='PASS',empty_disconnected_unbounded_checks='PASS')
}
