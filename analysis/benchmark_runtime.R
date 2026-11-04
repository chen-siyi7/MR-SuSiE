a <- grep('^--file=',commandArgs(),value=TRUE)
root <- normalizePath(file.path(dirname(gsub('~+~',' ',sub('^--file=','',a[1]),fixed=TRUE)),'..'))
source(file.path(root,'R/load.R'))
source(file.path(root,'R/analysis_helpers.R'))
source(file.path(root,'R/cml_model.R'))
args <- mr_susie_arguments(commandArgs(TRUE))
if(!identical(args$run,'true')) stop('Use --run=true to request the ten paired timing attempts')
Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,VECLIB_MAXIMUM_THREADS=1)
roster <- read.csv(file.path(root,'config/runtime_pairs.csv'),stringsAsFactors=FALSE)
stopifnot(nrow(roster)==10L)
output <- if(is.null(args$out)) file.path(root,'results/runtime_benchmark') else args$out
if(dir.exists(output)&&length(list.files(output))) stop('Output directory must be empty')
dir.create(output,recursive=TRUE,showWarnings=FALSE)
model <- make_mr_susie(root,file.path(root,'runtime/benchmark'))
model$warm('figure2')
model$warm('broader')
for(package in c('MRcML','MVMRcML','MVMRcMLSuSiE')) {
  if(!requireNamespace(package,quietly=TRUE)) stop('Missing dependency: ',package)
}
rows <- list()
for(i in seq_len(nrow(roster))) {
  j <- roster[i,]
  path <- file.path(root,j$prepared_file)
  stopifnot(digest::digest(file=path,algo='sha256',serialize=FALSE)==j$prepared_sha256)
  x <- readRDS(path)
  methods <- if(i%%2L) c('MVMR-SuSiE','MVMR-cML-SuSiE') else c('MVMR-cML-SuSiE','MVMR-SuSiE')
  for(method in methods) {
    seed <- as.integer(if(method=='MVMR-SuSiE') 1100000000+j$extension_index*101 else 1200000000+j$extension_index*101)
    set.seed(seed)
    start <- proc.time()
    z <- if(method=='MVMR-SuSiE') mr_susie_capture(function() model$fit_observed(x$observed,list(baseline=x$controls,context=x$context))) else fit_cml_bic(x$observed,root,seed)
    timing <- proc.time()-start
    k <- x$observed$dat$num_exposures
    n <- if(method=='MVMR-SuSiE') if(isTRUE(z$result$point_available)) sum(is.finite(z$result$point)) else 0L else sum(is.finite(z$result$beta))
    rows[[length(rows)+1L]] <- data.frame(pair=i,job_id=j$job_id,method=method,order=match(method,methods),
      variants=x$observed$dat$m,exposures=k,wall_seconds=unname(timing[3]),cpu_seconds=unname(sum(timing[1:2])),
      available_coefficients=n,error=z$error,warnings=paste(z$warnings,collapse=' | '),seed=seed,
      prepared_sha256=j$prepared_sha256)
    saveRDS(z,file.path(output,paste0('pair_',sprintf('%02d',i),'_',if(method=='MVMR-SuSiE') 'mr' else 'cml','.rds')))
    write.csv(do.call(rbind,rows),file.path(output,'timing.csv'),row.names=FALSE,na='NA')
  }
}
writeLines(capture.output(sessionInfo()),file.path(output,'session_info.txt'))
