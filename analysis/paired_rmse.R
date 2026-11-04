a <- grep('^--file=',commandArgs(),value=TRUE)
root <- normalizePath(file.path(dirname(gsub('~+~',' ',sub('^--file=','',a[1]),fixed=TRUE)),'..'))
source(file.path(root,'R/analysis_helpers.R'))
args <- mr_susie_arguments(commandArgs(TRUE))
stopifnot(!is.null(args$input),!is.null(args$out),!is.null(args$method1),!is.null(args$method2))
x <- read.csv(args$input,stringsAsFactors=FALSE)
x <- x[x$estimand=='conditional_direct_effect'&is.finite(x$truth)&x$truth!=0&x$point_available&is.finite(x$beta),]
keys <- c('cohort','panel','mode','dataset','exposure')
y <- merge(x[x$method==args$method1,c(keys,'beta','truth')],x[x$method==args$method2,c(keys,'beta','truth')],by=keys,suffixes=c('_1','_2'))
stopifnot(nrow(y)>0L,all(y$truth_1==y$truth_2))
y$error_1 <- (y$beta_1-y$truth_1)^2
y$error_2 <- (y$beta_2-y$truth_2)^2
set.seed(20261104)
nboot <- if(is.null(args$bootstrap)) 2000L else as.integer(args$bootstrap)
result <- lapply(split(y,interaction(y$panel,y$mode,drop=TRUE)),function(d) {
  stat <- function(z) {
    s <- aggregate(cbind(error_1,error_2)~cohort,z,mean)
    sqrt(mean(s$error_1))-sqrt(mean(s$error_2))
  }
  draws <- replicate(nboot,{
    z <- lapply(split(d,d$cohort),function(g) {
      ids <- unique(g$dataset)
      chosen <- sample(ids,length(ids),replace=TRUE)
      do.call(rbind,lapply(chosen,function(id) g[g$dataset==id,,drop=FALSE]))
    })
    stat(do.call(rbind,z))
  })
  data.frame(panel=d$panel[1],mode=d$mode[1],settings=length(unique(d$cohort)),
    paired_targets=nrow(d),paired_datasets=length(unique(paste(d$cohort,d$dataset))),
    method1=args$method1,method2=args$method2,rmse_difference=stat(d),
    mc_lower=quantile(draws,.025),mc_upper=quantile(draws,.975),bootstrap=nboot)
})
write.csv(do.call(rbind,result),args$out,row.names=FALSE)
