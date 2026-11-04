make_binary_score_selector <- function(original_score_function,
                                        rss_function=susieR::susie_rss,
                                        return_posterior=FALSE) {
  txt<-paste(deparse(original_score_function),collapse='\n')
  needle<-'res <- do.call(susieR::susie_rss, rss_args)'
  hits<-gregexpr(needle,txt,fixed=TRUE)[[1]]
  stopifnot(length(hits)==1L,hits[1]>0)
  txt<-sub(needle,paste0('rss_args$max_iter <- 1000L\n',
    'rss_args$refine <- TRUE\nres <- do.call(binary_susie_rss, rss_args)'),txt,fixed=TRUE)
  e<-new.env(parent=environment(original_score_function));e$binary_susie_rss<-rss_function
  f<-eval(parse(text=txt),e);environment(f)<-e
  if(return_posterior) {
    b<-body(f);last<-length(b)
    b[[last]]<-substitute({answer<-OLD;answer$susie_fit<-res;answer},list(OLD=b[[last]]))
    body(f)<-b
  }
  f
}

primary_application_profile <- function(observed,selection,refit_function,
                                        selection_warnings=character()) {
  d<-observed$dat;d$seZX[]<-0
  f<-refit_function(d,selection,selection_warnings)
  f$primary_measurement_error_correction<-FALSE
  f$application_scale<-'Outcome log odds per exposure unit; standardized genotype rows already supplied'
  f$information_limits<-'Working conditional intervals; shared outcomes, reference LD, selected support and known source overlap remain unintegrated'
  f
}
