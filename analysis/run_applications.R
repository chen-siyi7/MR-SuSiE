a <- grep('^--file=',commandArgs(),value=TRUE)
root <- normalizePath(file.path(dirname(gsub('~+~',' ',sub('^--file=','',a[1]),fixed=TRUE)),'..'))
source(file.path(root,'R/load.R'))
source(file.path(root,'R/analysis_helpers.R'))
source(file.path(root,'R/comparator_models.R'))
args <- mr_susie_arguments(commandArgs(TRUE))
if(is.null(args$input)&&is.null(args$manifest)) stop('Supply --input=FILE or --manifest=FILE')
model <- make_mr_susie(root,file.path(root,'runtime/applications'))
output <- if(is.null(args$out)) file.path(root,'results/applications') else args$out
if(dir.exists(output)&&length(list.files(output))) stop('Output directory must be empty')
dir.create(output,recursive=TRUE,showWarnings=FALSE)
if(!is.null(args$input)) {
  manifest <- data.frame(job_id='example',input_file=args$input,n_eff=if(is.null(args$n_eff)) NA_real_ else as.numeric(args$n_eff),configuration='user_input')
} else manifest <- read.csv(args$manifest,stringsAsFactors=FALSE)
limit <- if(is.null(args$limit)) 1L else if(args$limit=='all') nrow(manifest) else as.integer(args$limit)
stopifnot(limit>=1L)
manifest <- head(manifest,limit)
results <- list()
for(i in seq_len(nrow(manifest))) {
  row <- manifest[i,]
  path <- row$input_file
  if(!file.exists(path)) path <- file.path(root,path)
  stopifnot(file.exists(path))
  n_eff <- if(is.na(row$n_eff)) NULL else row$n_eff
  observed <- mr_susie_read_input(path,model,n_eff,TRUE)
  result <- mr_susie_analyze(observed,model,context=list(kind='application',family='application'),
    cohort=row$configuration,dataset=row$job_id,mode='NOME',
    comparators=identical(args$comparators,'true'),seed=601200000L+as.integer(i)*1000L,
    me_sensitivity=identical(args$me,'true'))
  result$input_sha256 <- digest::digest(file=path,algo='sha256',serialize=FALSE)
  saveRDS(result,file.path(output,paste0(row$job_id,'.rds')))
  results[[i]] <- result$rows
}
write.csv(do.call(rbind,results),file.path(output,'coefficient_results.csv'),row.names=FALSE,na='NA')
writeLines(capture.output(sessionInfo()),file.path(output,'session_info.txt'))
