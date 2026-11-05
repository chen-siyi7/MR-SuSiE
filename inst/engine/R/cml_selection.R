native_cml <- function(dat, seed, outer_guard=100L, memoize=TRUE) {
  stopifnot(requireNamespace("MVMRcMLSuSiE",quietly=TRUE),
    as.character(packageVersion("MVMRcMLSuSiE"))=="0.1.0")
  suppressPackageStartupMessages(library(MRcML))
  suppressPackageStartupMessages(library(susieR))
  start <- proc.time()[3]
  K <- ncol(dat$betaZX); m <- nrow(dat$betaZX)
  pval <- 2*pnorm(abs(dat$betaZX/dat$seZX),lower.tail=FALSE)

  instruments <- lapply(seq_len(K),function(k)which(pval[,k]<5e-8))
  step1 <- step2 <- step3 <- NULL
  selected <- invalid <- screened <- integer()
  warnings <- character(); error <- ""; stdout <- character(); stage <- "step1"
  finished <- FALSE
  set.seed(seed,kind="Mersenne-Twister",normal.kind="Inversion",sample.kind="Rejection")
  cache <- list(); cache_hits <- 0L; actual_fits <- 0L
  api_version <- list(version=as.character(packageVersion("MRcML")),
    commit=packageDescription("MRcML")$RemoteSha,
    body_sha256=digest::digest(body(MRcML::mr_cML),algo="sha256"))
  cached_mr_cML <- function(...) {
    key <- list(api=api_version,args=list(...))
    found <- which(vapply(cache,function(z)identical(z$args,key),logical(1)))
    if(length(found)) {
      hit <- cache[[found[1]]]; cache_hits <<- cache_hits+1L
      assign(".Random.seed",hit$rng,envir=.GlobalEnv)
      for(w in hit$warnings)warning(w)
      return(hit$result)
    }
    ww <- list(); actual_fits <<- actual_fits+1L
    result <- withCallingHandlers(MRcML::mr_cML(...),
      warning=function(w)ww[[length(ww)+1L]] <<- w)
    cache[[length(cache)+1L]] <<- list(args=key,result=result,
      rng=get(".Random.seed",envir=.GlobalEnv),warnings=ww)
    result
  }
  step1_function <- getExportedValue("MVMRcMLSuSiE","mvmr.cml.susie.step1")
  step2_function <- getExportedValue("MVMRcMLSuSiE","mvmr.cml.susie.step2")
  if(memoize) {
    e1 <- new.env(parent=environment(step1_function)); e1$mr_cML <- cached_mr_cML
    e2 <- new.env(parent=environment(step2_function)); e2$mr_cML <- cached_mr_cML
    environment(step1_function) <- e1; environment(step2_function) <- e2
  }
  body_guard <- getExportedValue("MVMRcMLSuSiE","mvmr.cml.susie.step3")

  bb <- as.list(body(body_guard)); wi <- which(vapply(bb,function(z)
    is.call(z)&&identical(z[[1]],as.name("while")),logical(1)))
  stopifnot(length(wi)==1L)
  original_loop <- as.list(bb[[wi]][[3]])
  guard_expr <- substitute(if(iter>=LIMIT)stop("native_cml_outer_guard_reached",call.=FALSE),
    list(LIMIT=as.integer(outer_guard)))
  bb[[wi]][[3]] <- as.call(c(list(as.name("{"),guard_expr),original_loop[-1]))

  ret <- which(vapply(bb,function(z)is.call(z)&&identical(z[[1]],as.name("return")),logical(1)))
  stopifnot(length(ret)==1L)
  bb <- append(bb,list(quote(res$outer_iterations <- iter),quote(res$outer_diff <- diff)),after=ret-1L)
  body(body_guard) <- as.call(bb)
  tryCatch(withCallingHandlers({
    stdout <- capture.output({
      if(any(lengths(instruments)==0L)) stop("no_significant_instruments_for_step1")
      step1 <- step1_function(
        sample.sizes=c(rep(dat$n1,K),dat$n2),
        beta.exposure.ls=lapply(seq_len(K),function(k)dat$betaZX[instruments[[k]],k]),
        se.exposure.ls=lapply(seq_len(K),function(k)dat$seZX[instruments[[k]],k]),
        beta.outcome.ls=lapply(seq_len(K),function(k)dat$betaZY[instruments[[k]]]),
        se.outcome.ls=lapply(seq_len(K),function(k)dat$se_betaZY[instruments[[k]]]),
        use.openGWAS=FALSE)
      if(length(step1)!=K || any(!is.finite(step1)))stop("nonfinite_step1_screening_pvalues")
      screened <- which(step1 < .05/K)
      if(length(screened)) {
        stage <- "step2"
        step2 <- step2_function(
          sample.sizes.subset=c(rep(dat$n1,length(screened)),dat$n2),
          beta.exposure.mat=dat$betaZX[,screened,drop=FALSE],
          se.exposure.mat=dat$seZX[,screened,drop=FALSE],
          beta.outcome.vec=dat$betaZY,se.outcome.vec=dat$se_betaZY,
          pval.exposure.mat=pval[,screened,drop=FALSE],use.openGWAS=FALSE,cutoff=5e-8)
        invalid <- sort(unique(as.integer(step2$invalid.idx)))
        if(any(!is.finite(step2$theta.vec)))stop("nonfinite_step2_starting_coefficients")
        stage <- "step3"
        rho <- diag(length(screened)+1L)
        rho[seq_along(screened),seq_along(screened)] <- dat$sigma_exp[screened,screened,drop=FALSE]
        step3 <- body_guard(step2$mvdat,step2$invalid.idx,step2$theta.vec,rho,
          S=10,max.iter=200,tol=1e-10)
        if(!isTRUE(step3$converged))stop("native_cml_inner_nonconvergence")
        if(any(!is.finite(step3$alpha)))stop("nonfinite_step3_assignment_probabilities")
        selected <- screened[which(colSums(step3$alpha>1/length(screened))>0)]
      }
      finished <- TRUE; stage <- "complete"
    })
  },warning=function(w){warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}),
  error=function(e)error <<- conditionMessage(e))
  nonconvergence <- any(grepl("did not converge|not converged|failed to converge|unconverged",
    warnings,ignore.case=TRUE))
  usable <- finished && !nonconvergence
  list(exposures=if(usable)selected else integer(),invalid=if(usable)invalid else integer(),
    returned_exposures=selected,returned_invalid=invalid,screened=screened,
    step1=step1,step2=step2,step3=step3,step1_instruments=instruments,
    status=if(usable)"completed" else if(nzchar(error))error else "explicit_nonconvergence_warning",
    stage=stage,completed=usable,returned=finished,warnings=unique(warnings),error=error,
    stdout=stdout,seed=seed,upstream_fixed_random_seed=1L,random_starts=100L,
    memoization_enabled=memoize,cached_duplicate_fits=cache_hits,
    cache_api_version=api_version,
    executed_uvmr_fits=if(memoize)actual_fits else K+length(screened),
    outer_guard=outer_guard,inner_max_iter=200L,outer_tol=1e-10,
    elapsed_seconds=unname(proc.time()[3]-start))
}
