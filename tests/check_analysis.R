a <- grep('^--file=',commandArgs(),value=TRUE)
root <- normalizePath(file.path(dirname(gsub('~+~',' ',sub('^--file=','',a[1]),fixed=TRUE)),'..'))
work <- tempfile('mr_susie_analysis_')
dir.create(work)
rscript <- file.path(R.home('bin'),'Rscript')
run <- function(script,args) system2(rscript,c(shQuote(file.path(root,script)),shQuote(args)),stdout=TRUE,stderr=TRUE)
x <- data.frame(cohort='toy',panel='full',mode='ME',method='MVMR-SuSiE',
  estimand='conditional_direct_effect',dataset=c('a','a','b','b'),exposure=c(1,2,1,2),
  truth=c(.2,0,.2,0),beta=c(.1,.2,NA,.4),score_center=c(.1,.2,NA,.4),se=c(.1,.1,NA,.1),
  p=c(.01,.01,NA,.5),lower=c(0,.1,NA,.3),upper=c(.3,.3,NA,.5),
  point_available=c(TRUE,TRUE,FALSE,TRUE),test_available=c(TRUE,TRUE,FALSE,TRUE),phase='evaluation')
input <- file.path(work,'input.csv')
metrics <- file.path(work,'metrics.csv')
write.csv(x,input,row.names=FALSE,na='NA')
z <- run('analysis/evaluate.R',c(paste0('--input=',input),paste0('--out=',metrics)))
stopifnot(is.null(attr(z,'status')),file.exists(metrics))
y <- read.csv(metrics)
stopifnot(abs(y$mse-.01)<1e-12,y$point_availability==.5,y$coverage==1,
  y$nominal_detection==.5,y$null_rejection==.5,y$coefficient_fdr==.25)
cutoff <- file.path(work,'cutoff.csv')
z <- suppressWarnings(run('analysis/freeze_calibration.R',c(paste0('--input=',input),paste0('--out=',cutoff))))
stopifnot(!is.null(attr(z,'status')),!file.exists(cutoff))
x <- x[rep(2,20),]
x$dataset <- paste0('cal_',seq_len(nrow(x)))
x$phase <- 'calibration'
x$score_center <- 1
x$se <- .2
write.csv(x,input,row.names=FALSE,na='NA')
z <- run('analysis/freeze_calibration.R',c(paste0('--input=',input),paste0('--out=',cutoff)))
stopifnot(is.null(attr(z,'status')),read.csv(cutoff)$cutoff==5)
before <- readLines(cutoff)
z <- suppressWarnings(run('analysis/freeze_calibration.R',c(paste0('--input=',input),paste0('--out=',cutoff))))
stopifnot(!is.null(attr(z,'status')),identical(readLines(cutoff),before))
unlink(work,recursive=TRUE)
cat('Scheduled-denominator, failure, FDR and calibration-freeze checks passed\n')
