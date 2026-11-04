a <- grep('^--file=',commandArgs(),value=TRUE)
root <- normalizePath(file.path(dirname(gsub('~+~',' ',sub('^--file=','',a[1]),fixed=TRUE)),'..'))
source(file.path(root,'R/analysis_helpers.R'))
args <- mr_susie_arguments(commandArgs(TRUE))
input <- if(is.null(args$input)) file.path(root,'data/Table_S1.xlsx') else args$input
output <- if(is.null(args$out)) file.path(root,'results/table_s1') else args$out
if(dir.exists(output)&&length(list.files(output))) stop('Output directory must be empty')
dir.create(output,recursive=TRUE,showWarnings=FALSE)
for(sheet in readxl::excel_sheets(input)) {
  d <- as.data.frame(readxl::read_excel(input,sheet=sheet),check.names=FALSE)
  name <- gsub('[^A-Za-z0-9_]+','_',sheet)
  write.csv(d,file.path(output,paste0(name,'.csv')),row.names=FALSE,na='NA')
}
