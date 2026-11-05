a <- grep('^--file=',commandArgs(),value=TRUE)
root <- normalizePath(file.path(dirname(gsub('~+~',' ',sub('^--file=','',a[1]),fixed=TRUE)),'..'))
source(file.path(root,'R/load.R'))
model <- make_mr_susie(root,file.path(root,'runtime'))
for(name in c('scalar_strong_0001','joint_strong_0001')) {
  x <- jsonlite::read_json(file.path(root,'examples',paste0(name,'.json')),simplifyVector=TRUE)
  k <- ncol(x$B)
  m <- nrow(x$B)
  answer <- model$fit(x$B,matrix(x$sx,m,k),x$y,rep(x$sy,m),x$R,x$Sigma,
    x$selection_y,rep(x$sy,m),10000,rep('Outcome SD per exposure SD',k))
  actual <- data.frame(point=answer$point,answer$score_inference)
  expected <- read.csv(file.path(root,'examples',paste0(name,'_expected.csv')))
  stopifnot(isTRUE(all.equal(as.matrix(actual),as.matrix(expected),tolerance=1e-9,check.attributes=FALSE)))
  cat(name,': PASS\n')
  print(actual)
}
