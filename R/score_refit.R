make_joint_selector<-function(original_joint_function) {
 txt<-paste(deparse(original_joint_function),collapse='\n')
 hits<-gregexpr('check_input = TRUE',txt,fixed=TRUE)[[1]]
 if(length(hits)!=1L||hits[1]<0)stop('Expected exactly one pinned SuSiE input-check argument')
 f<-eval(parse(text=sub('check_input = TRUE','check_input = TRUE, max_iter = 1000, refine = TRUE',txt,fixed=TRUE)),envir=environment(original_joint_function))
 environment(f)<-environment(original_joint_function);f
}

corrected_score_refit<-function(dat,selection,selection_warnings=character()) {
 selector_ok<-isTRUE(selection$converged)&&!any(grepl('did not converge|not converg|non.?converg',selection_warnings,ignore.case=TRUE))
 f<-fit_profile_corrected(dat,selection$invalid,which_R='ref',weight_mode='iterated',bread='expected')
 f$selector_usable<-selector_ok
 if(!isTRUE(f$identified)||!selector_ok) {
  f$identified<-FALSE;f$se<-rep(NA_real_,dat$num_exposures);f$p<-rep(NA_real_,dat$num_exposures)
  f$critical_95<-NA_real_;f$working_phi<-NA_real_;f$residual_df<-NA_integer_;return(f)
 }
 a<-profile_prepare(dat,selection$invalid,'ref');st<-f$inference;H<-cbind(a$B,a$Z)
 df<-a$m-qr(H,tol=1e-8)$rank;if(df<=0)stop('No residual degrees of freedom')
 res<-a$y-as.numeric(a$B%*%f$beta)-if(ncol(a$Z))as.numeric(a$Z%*%f$direct)else 0
 phi<-max(1,as.numeric(crossprod(res,st$W%*%res))/df)
 linear<-crossprod(a$B,st$M%*%st$Omega%*%st$M%*%a$B);Ai<-solve(st$A)
 V<-Ai%*%(phi*linear+st$quadratic_term)%*%t(Ai);V<-(V+t(V))/2
 f$gaussian_score_se<-f$se;f$working_phi_covariance<-V;f$se<-sqrt(diag(V));f$p<-2*pt(-abs(f$beta/f$se),df)
 f$critical_95<-qt(.975,df);f$working_phi<-phi;f$residual_df<-df
 f$identified<-all(is.finite(f$se))&&all(f$se>0)
 f$scope<-'Joint selection and corrected-score inference with working residual heterogeneity and Student-t reference'
 f
}
