a <- grep('^--file=',commandArgs(),value=TRUE)
root <- normalizePath(file.path(dirname(gsub('~+~',' ',sub('^--file=','',a[1]),fixed=TRUE)),'..'))
source(file.path(root,'R/analysis_helpers.R'))
args <- mr_susie_arguments(commandArgs(TRUE))
stopifnot(!is.null(args$summary),!is.null(args$ld),!is.null(args$traits),!is.null(args$out),!file.exists(args$out))
traits <- strsplit(args$traits,',',fixed=TRUE)[[1]]
d <- read.csv(args$summary,stringsAsFactors=FALSE,check.names=FALSE)
R <- as.matrix(read.csv(args$ld,row.names=1,check.names=FALSE))
fields <- c('variant','effect_allele','other_allele','reference_effect_allele','reference_other_allele','reference_sd','beta_outcome','se_outcome',paste0('beta_',traits),paste0('se_',traits))
stopifnot(all(fields%in%names(d)),!anyDuplicated(d$variant),identical(rownames(R),colnames(R)))
for(field in c('effect_allele','other_allele','reference_effect_allele','reference_other_allele')) d[[field]] <- toupper(d[[field]])
ambiguous <- paste0(d$effect_allele,d$other_allele)%in%c('AT','TA','CG','GC')
same <- d$effect_allele==d$reference_effect_allele&d$other_allele==d$reference_other_allele
swap <- d$effect_allele==d$reference_other_allele&d$other_allele==d$reference_effect_allele
keep <- !ambiguous&(same|swap)&complete.cases(d[,fields])&d$reference_sd>0&d$se_outcome>0
keep <- keep&apply(d[paste0('se_',traits)],1,function(z) all(z>0))&d$variant%in%rownames(R)
d <- d[keep,,drop=FALSE]
sign <- ifelse(swap[keep],-1,1)
R <- R[d$variant,d$variant,drop=FALSE]
B <- as.matrix(d[paste0('beta_',traits)])*sign*d$reference_sd
sx <- as.matrix(d[paste0('se_',traits)])*d$reference_sd
Sigma <- if(is.null(args$sigma)) diag(length(traits)) else as.matrix(read.csv(args$sigma,row.names=1,check.names=FALSE))
stopifnot(nrow(Sigma)==length(traits),ncol(Sigma)==length(traits),all(diag(R)>0))
value <- list(m=nrow(d),num_exposures=length(traits),betaZX=B,seZX=sx,
  betaZY=d$beta_outcome*sign*d$reference_sd,se_betaZY=d$se_outcome*d$reference_sd,
  corr_in=R,corr_ref=R,sigma_exp=Sigma,overlap=0,rsid=d$variant,traits=traits,
  n1=if(is.null(args$n1)) NA_real_ else as.numeric(args$n1),
  n2=if(is.null(args$n2)) NULL else as.numeric(args$n2),reference_dosage_sd=d$reference_sd)
saveRDS(value,args$out)
