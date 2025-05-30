\name{txtsummary}

\alias{txtsummary}

\title{Median and Coefficient Interval}

\description{Summarization of a numeric vector.}

\usage{
txtsummary(x,
           f = c("median", "mean"),
           digits = 0,
           scientific = FALSE,
           range = c("IQR", "95\%CI", "range", "sd"))
}

\arguments{
  \item{x}{a numeric vector.}
  \item{f}{xxx.}
  \item{digits}{how many significant digits are to be used.}
  \item{scientific}{either a logical specifying whether result should be encoded in scientific format.}
  \item{range}{the range to be visualized.}

}

\value{
    The function returns the median and the range (interquartile or 95\% coefficient interval) of numeric vetor.
  }

\author{Stefano Cacciatore}


  \seealso{\code{\link{categorical.test}},\code{\link{continuous.test}},\code{\link{correlation.test}}, \code{\link{txtsummary}}}

\examples{
data(prostate)

txtsummary(prostate[,"BMI"])


}
\keyword{txtsummary}
