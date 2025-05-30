\name{correlation.test}

\alias{correlation.test}

\title{Continuous Information}

\description{Summarization of the continuous information.}

\usage{
correlation.test (x,
                  y,
                  method = c("pearson", "spearman","MINE"), 
                  name=NA, 
                  perm=100, ...)
}

\arguments{
  \item{x}{a numeric vector.}
  \item{y}{a numeric vector.}
  \item{method}{a character string indicating which correlation method is to be computed. "pearson" (default), "spearman", or "MINE".}
  \item{name}{the name of the feature.}
  \item{perm}{number of permutation needed to estimate the p-value with MINE correlation.}
  \item{...}{further arguments to be passed to or from methods.}
}

\value{
    The function returns a table with the summarized information. 
  }

\author{Stefano Cacciatore}

  \seealso{\code{\link{categorical.test}},\code{\link{continuous.test}}, \code{\link{txtsummary}}}

\examples{
data(prostate)

correlation.test(prostate[,"Age"],prostate[,"BMI"],name="correlation between Age and BMI")

}
\keyword{correlation.test}
