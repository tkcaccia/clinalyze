\name{add_analysis}
\alias{add_analysis}
\title{Add analysis results to a clinical object}
\description{
  This function appends the results of an analysis to a clinical object.
}
\usage{
add_analysis(object, name, results)
}
\arguments{
  \item{object}{A clinical object to which results will be attached.}
  \item{name}{The name (string) of the analysis.}
  \item{results}{The result object (e.g., data frame, matrix) to attach.}
}
\value{
  Returns the updated clinical object with the new analysis appended.
}
\examples{
data(clinical)
clinical <- add_analysis(clinical, "PCA", pca_results)
}
\author{Stefano Cacciatore}
