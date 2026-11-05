\name{union_ar}
\alias{union_ar}
\title{Support-union Gaussian sensitivity}
\description{Constructs effect confidence sets by taking a union over candidate invalid-variant supports.}
\usage{
union_ar (B, y, R, Sigma, sx = 0.01, sy = 0.01, smax = 3L, level = 0.95)
}
\arguments{
\item{B}{Matrix of marginal variant-exposure associations.}
\item{y}{Vector of marginal variant-outcome associations.}
\item{R}{Known, positive-definite signed LD matrix.}
\item{Sigma}{Exposure-GWAS estimation-error correlation matrix.}
\item{sx}{Common exposure coefficient standard error.}
\item{sy}{Common outcome coefficient standard error.}
\item{smax}{Prespecified maximum number of invalid variants.}
\item{level}{Requested confidence level, strictly between zero and one.}
}
\value{A list with projected confidence \code{sets}, the union test \code{joint_null_p}, and diagnostic fields. Empty, disconnected and unbounded sets are retained.}
\details{The construction assumes independent Gaussian exposure and outcome errors, homogeneous standard errors and known LD. Its coverage guarantee does not extend automatically to an estimated LD reference. The invalid-variant bound is supplied before examining outcomes.}
\examples{
x <- jsonlite::read_json(
  system.file("examples", "support_union_example.json", package="MRSuSiE"),
  simplifyVector=TRUE
)
fit <- do.call(union_ar, x)
fit$sets
}
