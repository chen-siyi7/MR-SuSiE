main <- function() {
  a <- grep('^--file=',commandArgs(),value=TRUE)
  root <- normalizePath(file.path(dirname(gsub('~+~',' ',sub('^--file=','',a[1]),fixed=TRUE)),'..'))
  source(file.path(root,'R/analysis_helpers.R'))
  source(file.path(root,'R/support_union.R'))
  supplied <- commandArgs(TRUE)
  if(any(!grepl('^--(input|out|sx|sy|smax|level)=',supplied))) stop('Use named options of the form --input=FILE --out=FILE.rds')
  args <- mr_susie_arguments(supplied)
  if(is.null(args$input)) {
    if(length(supplied)) stop('Supply --input=FILE or run without options for the self-check')
    print(support_union_self_check())
    return(invisible(NULL))
  }
  if(is.null(args$out)) stop('Supply --out=FILE.rds')
  if(!grepl('[.]rds$',args$out,ignore.case=TRUE)) stop('Output must be an RDS file')
  if(file.exists(args$out)) stop('Output file already exists')
  x <- if(grepl('[.]json$',args$input,ignore.case=TRUE)) jsonlite::read_json(args$input,simplifyVector=TRUE) else readRDS(args$input)
  get <- function(name,default) if(!is.null(args[[name]])) as.numeric(args[[name]]) else if(!is.null(x[[name]])) x[[name]] else default
  result <- union_ar(x$B,x$y,x$R,x$Sigma,get('sx',.01),get('sy',.01),get('smax',3L),get('level',.95))
  dir.create(dirname(args$out),recursive=TRUE,showWarnings=FALSE)
  saveRDS(result,args$out)
  invisible(result)
}
if(sys.nframe()==0L) main()
