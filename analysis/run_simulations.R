a <- grep('^--file=',commandArgs(),value=TRUE)
root <- normalizePath(file.path(dirname(gsub('~+~',' ',sub('^--file=','',a[1]),fixed=TRUE)),'..'))
source(file.path(root,'R/load.R'))
source(file.path(root,'R/analysis_helpers.R'))
source(file.path(root,'R/simulation_designs.R'))
source(file.path(root,'R/comparator_models.R'))
source(file.path(root,'R/reference_models.R'))
args <- mr_susie_arguments(commandArgs(TRUE))
if(is.null(args$plan)) stop('Supply --plan=FILE')
jobs <- read.csv(args$plan,stringsAsFactors=FALSE)
limit <- if(is.null(args$limit)) 1L else if(args$limit=='all') nrow(jobs) else as.integer(args$limit)
stopifnot(limit>=1L)
jobs <- head(jobs,limit)
output <- if(is.null(args$out)) file.path(root,'results/simulations') else args$out
if(dir.exists(output)&&length(list.files(output))) stop('Output directory must be empty')
dir.create(output,recursive=TRUE,showWarnings=FALSE)
model <- make_mr_susie(root,file.path(root,'runtime/simulations'))
results <- list()
for(i in seq_len(nrow(jobs))) {
  job <- as.list(jobs[i,])
  if(!is.null(job$engine)&&job$engine=='confirmation'||!is.null(job$source_group)&&job$source_group=='confirmation') job$family <- job$configuration
  x <- generate_simulation(root,job,model)
  fitting_seed <- if(!is.null(job$fit_seed)&&is.finite(as.numeric(job$fit_seed))) as.integer(job$fit_seed) else as.integer(job$generation_seed)
  native_seed <- if(!is.null(job$native_fit_seed)&&is.finite(as.numeric(job$native_fit_seed))) as.integer(job$native_fit_seed) else fitting_seed
  cohort <- if(!is.null(job$cohort)&&nzchar(job$cohort)) job$cohort else paste0(job$family,'_s',sprintf('%02d',job$setting))
  fit <- mr_susie_analyze(x$observed,model,x$controls,x$context,cohort,job$job_id,
    comparators=identical(args$comparators,'true'),seed=fitting_seed,native_seed=native_seed)
  if(!identical(args$benchmark,'false')) {
    reference <- fit_true_causal_exposures(model,x$observed,fit,x$theta,cohort,job$job_id)
    fit$benchmark <- reference$fits
    fit$rows <- rbind(fit$rows,reference$rows)
  }
  fit$evaluation_truth <- list(theta=x$theta,invalid=x$invalid)
  fit$job <- job
  fit$rows$truth <- ifelse(fit$rows$estimand=='conditional_direct_effect',x$theta[fit$rows$exposure],NA_real_)
  fit$rows$phase <- job$phase
  saveRDS(fit,file.path(output,paste0(job$job_id,'.rds')))
  results[[i]] <- fit$rows
}
write.csv(do.call(rbind,results),file.path(output,'coefficient_results.csv'),row.names=FALSE,na='NA')
writeLines(capture.output(sessionInfo()),file.path(output,'session_info.txt'))
