\name{add_analysis}
\alias{add_analysis}
\title{Add analysis results to a clinical object}
\description{
  This function appends the results of an analysis to a clinical object.
}
\usage{
add_analysis(ma, name, x)
}
\arguments{
  \item{ma}{A clinical object to which results will be attached.}
  \item{name}{The name (string) of the analysis.}
  \item{x}{The result object (e.g., data frame, matrix) to attach.}
}
\value{
  Returns the updated clinical object with the new analysis appended.
}
\examples{
data(prostate)


ma=initialization(prostate[,"Hospital"])
ma=add_analysis(ma,"Gender",prostate[,"Gender"])
ma=add_analysis(ma,"Gleason score",prostate[,"Gleason score"])
ma=add_analysis(ma,"BMI",prostate[,"BMI"])
ma=add_analysis(ma,"Age",prostate[,"Age"])
ma
}
\author{Stefano Cacciatore}
