a <- grep('^--file=',commandArgs(),value=TRUE)
root <- normalizePath(file.path(dirname(gsub('~+~',' ',sub('^--file=','',a[1]),fixed=TRUE)),'..'))
source(file.path(root,'R/analysis_helpers.R'))
args <- mr_susie_arguments(commandArgs(TRUE))
stopifnot(!is.null(args$directory),!is.null(args$out))
files <- list.files(args$directory,pattern='[.]rds$',full.names=TRUE)
rows <- list()
for(file in files) {
  x <- readRDS(file)
  if(is.null(x$evaluation_truth)||is.null(x$fits)) next
  for(panel in names(x$fits)) {
    s <- x$fits[[panel]]$joint$result$score$selection
    ii <- x$fits[[panel]]$retained
    usable <- isTRUE(s$converged)
    for(role in c('exposure','invalid')) {
      truth <- if(role=='exposure') which(x$evaluation_truth$theta!=0) else x$evaluation_truth$invalid
      sets <- if(!usable) list() else if(role=='exposure') s$exposure_cs else lapply(s$invalid_cs,function(z) ii[z])
      selected <- sort(unique(unlist(sets)))
      denominator <- length(selected)
      coverage <- if(length(sets)) mean(vapply(sets,function(z) any(z%in%truth),logical(1))) else NA_real_
      rows[[length(rows)+1L]] <- data.frame(dataset=basename(file),panel=panel,role=role,
        selector_available=usable,true_members=length(truth),selected_members=denominator,
        sensitivity=if(length(truth)) length(intersect(selected,truth))/length(truth) else NA_real_,
        false_discovery_proportion=if(denominator) length(setdiff(selected,truth))/denominator else 0,
        credible_sets=length(sets),sets_covered=if(length(sets)) sum(vapply(sets,function(z) any(z%in%truth),logical(1))) else 0,
        credible_set_coverage=coverage,mixed_sets=length(s$unresolved_all))
    }
  }
}
if(!length(rows)) stop('No saved simulation selections with evaluation truth')
write.csv(do.call(rbind,rows),args$out,row.names=FALSE,na='NA')
