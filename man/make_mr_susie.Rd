\name{make_mr_susie}
\alias{make_mr_susie}
\title{Construct the common MR-SuSiE model}
\description{Creates the shared fitting interface for a single exposure (UVMR-SuSiE) or jointly modeled exposures (MVMR-SuSiE).}
\usage{
make_mr_susie (root = NULL, work_dir = file.path(tempdir(), "mr_susie"))
}
\arguments{
\item{root}{Optional source-checkout directory. The default uses the engine installed with the package.}
\item{work_dir}{Writable directory for temporary working files. The installed package uses precompiled native helpers.}
}
\value{A list containing \code{fit}, \code{fit_observed}, \code{uvmr}, \code{from_score}, formatting and diagnostic helpers, and the source and working paths.}
\details{The \code{fit} function accepts exposure coefficients \code{B}, their standard errors \code{sX}, outcome coefficients \code{y}, outcome standard errors \code{sy}, signed LD \code{R}, and exposure-error correlation \code{Sigma_exp}. Exposure coefficients and standard errors are matrices with one column per exposure. It also requires \code{selection_beta}, \code{selection_se}, \code{n_select}, and one label per exposure in \code{coefficient_units}.

Optional arguments specify the selection family, reuse of selection outcomes, reference sample size, localization controls and analysis context. The default selection family is Gaussian. Binary score localization uses \code{selector_family="binary_score"}. Setting \code{n_select=NULL} uses the large-sample score convention. The default assumes independent exposure and outcome estimation errors.

The returned \code{point} contains posterior-mean estimates. \code{score_inference} contains separate corrected-score centers, standard errors, interval bounds and p-values. Shrinkage does not recenter those intervals or tests. \code{selection}, \code{warnings}, \code{status} and \code{raw_fit} retain localization and numerical diagnostics. Interval validity depends on localization and the covariance assumptions. Coefficients must use comparable units.}
\examples{
model <- make_mr_susie()
x <- jsonlite::read_json(
  system.file("examples", "scalar_strong_0001.json", package="MRSuSiE"),
  simplifyVector=TRUE
)
m <- nrow(x$B)
fit <- model$fit(
  x$B, matrix(x$sx,m,1), x$y, rep(x$sy,m), x$R, x$Sigma,
  x$selection_y, rep(x$sy,m), 10000, "Outcome SD per exposure SD"
)
fit$point
fit$score_inference
}
