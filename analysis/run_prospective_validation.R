a <- grep('^--file=',commandArgs(),value=TRUE)
root <- normalizePath(file.path(dirname(gsub('~+~',' ',sub('^--file=','',a[1]),fixed=TRUE)),'..'))
source(file.path(root,'R/load.R'))
source(file.path(root,'R/analysis_helpers.R'))
source(file.path(root,'R/generators/prospective.R'))
args <- mr_susie_arguments(commandArgs(TRUE))
plan <- if(is.null(args$plan)) file.path(root,'config/prospective_validation_plan.csv') else args$plan
jobs <- read.csv(plan,stringsAsFactors=FALSE)
limit <- if(is.null(args$limit)) 1L else if(args$limit=='all') nrow(jobs) else as.integer(args$limit)
output <- if(is.null(args$out)) file.path(root,'results/prospective') else args$out
if(dir.exists(output)&&length(list.files(output))) stop('Output directory must be empty')
dir.create(output,recursive=TRUE,showWarnings=FALSE)
model <- make_mr_susie(root,file.path(root,'runtime/prospective'))
rows <- list()
for(i in seq_len(min(limit,nrow(jobs)))) {
  j <- as.list(jobs[i,])
  generated <- generate_prospective(j)
  modes <- if(j$rho>0) c('known_LD','reference_503') else 'known_LD'
  fits <- list()
  for(mode in modes) {
    d <- generated$data
    if(mode=='reference_503') d$corr_ref <- generated$reference_503
    o <- model$observed(d,generated$selection_beta,generated$selection_se,10000,'gaussian',FALSE,n_reference=if(mode=='known_LD') Inf else 503)
    z <- mr_susie_capture(function() model$fit_observed(o,list(baseline=list(L=10,prior=.5,coverage=.95,purity=.5),context=list(kind='simulation',family='broader'))))
    fits[[mode]] <- z
    method <- if(j$M==1L) 'UVMR-SuSiE' else 'MVMR-SuSiE'
    rr <- mr_susie_rows(z$result,method,'conditional_direct_effect',paste0('exposure_',seq_len(j$M)),j$cohort,j$job_id,mode,'ME',z$warnings,z$error)
    rr$truth <- generated$truth
    rows[[paste(i,mode)]] <- rr
  }
  saveRDS(list(job=j,observed=generated,fits=fits),file.path(output,paste0(j$job_id,'.rds')))
}
write.csv(do.call(rbind,rows),file.path(output,'coefficient_results.csv'),row.names=FALSE,na='NA')
