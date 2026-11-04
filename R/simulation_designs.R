make_simulation_generator <- function(root,family,rho=0) {
  e <- new.env(parent=globalenv())
  if(family %in% c('figure2','fresh_independent','four_baseline','four_high_correlation')) {
    e$N <- 50000L
    e$RHO <- rho
    sys.source(file.path(root,'R/generators',if(grepl('four_',family)) 'independent_four.R' else 'independent_six.R'),e)
    return(e$generate)
  }
  if(family=='broader') {
    sys.source(file.path(root,'R/generators/broader.R'),e)
    return(e$generate_data)
  }
  if(family %in% c('central','fresh_central','central_outcome_reuse','central_gaussian')) {
    e$config <- jsonlite::read_json(file.path(root,'config/central_design.json'),simplifyVector=TRUE)
    e$M <- as.integer(e$config$m)
    e$K <- as.integer(e$config$K)
    e$Rtrue <- e$config$ld_rho^abs(outer(seq_len(e$M),seq_len(e$M),'-'))
    e$Ltrue <- t(chol(e$Rtrue))
    e$Sg <- matrix(e$config$genetic_rho,e$K,e$K)
    diag(e$Sg) <- 1
    sys.source(file.path(root,'R/generators/central.R'),e)
    return(function(seed) {
      e$master <- as.integer(seed)-1L
      e$make_data(1L)
    })
  }
  stop('Unknown simulation family: ',family)
}

generate_simulation <- function(root,job,model) {
  seed <- as.integer(job$generation_seed)
  if(length(seed)!=1L||!is.finite(seed)) stop('An explicit generation_seed is required')
  family <- as.character(job$family)
  setting <- as.integer(job$setting)
  grid <- read.csv(file.path(root,'config/independent_settings.csv'))
  metadata <- jsonlite::read_json(file.path(root,'config/setting_metadata.json'),simplifyVector=TRUE)
  if(family %in% c('figure2','fresh_independent','four_baseline','four_high_correlation')) {
    k <- if(grepl('four_',family)) 4L else 6L
    theta <- if(k==6L) list(c(.2,.1,.3,.4,0,0),c(0,-.1,.1,.2,.3,0)) else list(c(.2,.1,.3,.4),c(0,-.1,.1,.2))
    gen <- make_simulation_generator(root,family,if(family=='four_high_correlation') .8 else 0)
    x <- gen(seed,theta[[as.integer(grid$theta_set[setting])]],grid$scenario[setting],grid$prop_invalid[setting],n=50000L)
    o <- model$observed(x$data,n_select=x$data$n2,gamma_mode='supplied',outcome_reuse=TRUE,n_reference=Inf)
    return(list(observed=o,theta=x$truth$theta,invalid=x$truth$true_invalid,raw=x,
      controls=list(L=10,prior=.5,coverage=.95,purity=.5),context=list(kind='simulation',family=if(family=='fresh_independent') 'figure2' else family)))
  }
  if(family=='broader') {
    g <- metadata[[sprintf('broader_%02d',setting)]]
    for(name in c('ref_rho','n_ref')) g[[name]] <- as.numeric(g[[name]])
    if(isTRUE(g$reused_central_data)) {
      gen <- make_simulation_generator(root,'central')
      d <- gen(seed)
      d$theta <- c(.2,.1,.3,.4,0,0)
    } else {
      gen <- make_simulation_generator(root,'broader')
      g$seed_base <- seed-as.integer(job$replicate)
      d <- gen(g,as.integer(job$replicate))
    }
    d$dat$n1 <- g$n_exposure
    d$dat$n2 <- g$n_estimate
    o <- model$observed(d$dat,d$select,rep(1/sqrt(g$n_select),d$dat$m),g$n_select,
      if(g$outcome_family=='binomial') 'binary_score_approximation' else 'gaussian',FALSE,n_reference=g$n_ref)
    return(list(observed=o,theta=d$theta,invalid=d$truth_invalid,raw=d,
      controls=list(L=g$L,prior=g$prior,coverage=.95,purity=g$purity),context=list(kind='simulation',family='broader')))
  }
  if(family=='individual_binary') {
    sys.source(file.path(root,'R/generators/confirmation.R'),envir=environment())
    x <- generate_binary(seed,root)
    d <- x$generated$dat
    o <- model$observed(d,x$generated$select,d$se_betaZY,x$n_select,'binary_score',TRUE,n_reference=Inf)
    return(list(observed=o,theta=x$truth,invalid=x$invalid,raw=x,
      controls=list(L=10,prior=.5,coverage=.95,purity=.5),context=list(kind='simulation',family='broader')))
  }
  gen <- make_simulation_generator(root,family)
  d <- gen(seed)
  d$dat$n1 <- d$dat$n2 <- 50000L
  reused <- family=='central_outcome_reuse'
  o <- model$observed(d$dat,if(reused) d$dat$betaZY else d$select,
    rep(1/sqrt(50000),d$dat$m),50000,'gaussian',reused,n_reference=500)
  list(observed=o,theta=c(.2,.1,.3,.4,0,0),invalid=d$truth_invalid,raw=d,
    controls=list(L=20,prior=.5,coverage=.95,purity=.5),context=list(kind='simulation',family='broader'))
}
