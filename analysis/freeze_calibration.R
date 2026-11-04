a <- grep('^--file=',commandArgs(),value=TRUE)
root <- normalizePath(file.path(dirname(gsub('~+~',' ',sub('^--file=','',a[1]),fixed=TRUE)),'..'))
source(file.path(root,'R/analysis_helpers.R'))
args <- mr_susie_arguments(commandArgs(TRUE))
stopifnot(!is.null(args$input),!is.null(args$out),!file.exists(args$out))
x <- read.csv(args$input,stringsAsFactors=FALSE)
stopifnot(all(x$phase=='calibration'))
x <- x[x$estimand=='conditional_direct_effect'&x$truth==0&x$test_available&is.finite(x$score_center)&is.finite(x$se)&x$se>0,]
if(!nrow(x)) stop('No usable calibration nulls')
key <- interaction(x$cohort,x$panel,x$mode,x$method,drop=TRUE)
result <- lapply(split(x,key),function(d) data.frame(cohort=d$cohort[1],panel=d$panel[1],mode=d$mode[1],method=d$method[1],
  usable_nulls=nrow(d),cutoff=as.numeric(quantile(abs(d$score_center/d$se),.95,type=8)),quantile_type=8))
write.csv(do.call(rbind,result),args$out,row.names=FALSE)
jsonlite::write_json(list(input_sha256=digest::digest(file=args$input,algo='sha256',serialize=FALSE),
  cutoff_sha256=digest::digest(file=args$out,algo='sha256',serialize=FALSE),phase='calibration',
  quantile_type=8,evaluation_used=FALSE),paste0(args$out,'.json'),pretty=TRUE,auto_unbox=TRUE)
