a <- grep('^--file=',commandArgs(),value=TRUE)
root <- normalizePath(file.path(dirname(gsub('~+~',' ',sub('^--file=','',a[1]),fixed=TRUE)),'..'))
source(file.path(root,'R/analysis_helpers.R'))
args <- mr_susie_arguments(commandArgs(TRUE))
stopifnot(!is.null(args$input),!is.null(args$out))
x <- read.csv(args$input,stringsAsFactors=FALSE)
x <- x[x$estimand=='conditional_direct_effect'&is.finite(x$truth),]
if(!nrow(x)) stop('No conditional-effect truth is available')
target_class <- if(is.null(args$target)) 'nonzero' else args$target
stopifnot(target_class%in%c('nonzero','null','all'))
keys <- interaction(x$cohort,x$panel,x$mode,x$method,drop=TRUE)
rows <- lapply(split(x,keys),function(d) {
  eligible <- if(target_class=='nonzero') d$truth!=0 else if(target_class=='null') d$truth==0 else rep(TRUE,nrow(d))
  ok <- d$point_available&is.finite(d$beta)&eligible
  interval <- d$test_available&is.finite(d$lower)&is.finite(d$upper)&eligible
  test <- d$test_available&is.finite(d$p)
  means <- function(z) if(length(z)) mean(z) else NA_real_
  discoveries <- test&d$p<.05
  fdp <- vapply(split(seq_len(nrow(d)),d$dataset),function(i) sum(discoveries[i]&d$truth[i]==0)/max(1,sum(discoveries[i])),numeric(1))
  data.frame(cohort=d$cohort[1],panel=d$panel[1],mode=d$mode[1],method=d$method[1],
    datasets=length(unique(d$dataset)),target_class=target_class,scheduled_targets=sum(eligible),scheduled_all_targets=nrow(d),available_points=sum(ok),
    point_availability=if(any(eligible)) mean(ok[eligible]) else NA_real_,mse=means((d$beta[ok]-d$truth[ok])^2),
    rmse=sqrt(means((d$beta[ok]-d$truth[ok])^2)),bias=means(d$beta[ok]-d$truth[ok]),
    available_intervals=sum(interval),coverage=means(d$lower[interval]<=d$truth[interval]&d$upper[interval]>=d$truth[interval]),
    mean_width=means(d$upper[interval]-d$lower[interval]),
    nominal_detection=means(discoveries[d$truth!=0]),null_rejection=if(d$method[1]=='True causal exposures') NA_real_ else means(discoveries[d$truth==0]),
    coefficient_fdr=if(d$method[1]=='True causal exposures') NA_real_ else mean(fdp))
})
write.csv(do.call(rbind,rows),args$out,row.names=FALSE,na='NA')
