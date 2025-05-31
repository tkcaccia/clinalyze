\name{add_analysis}
\alias{add_analysis}
\title{Add Statistical Analysis to a Clinical Table}
\description{
Adds results of a statistical test to a \code{"clinical.table"} object. It chooses the appropriate test depending on whether the input variable is numeric or factor.
}
\usage{
add_analysis(ma, name, x)
}
\arguments{
  \item{ma}{An object of class \code{"clinical.table"}.}
  \item{name}{A character string naming the feature being analyzed.}
  \item{x}{A numeric or factor vector representing the feature values to test.}
}
\value{
Returns the updated \code{"clinical.table"} object with test results appended.
}
\examples{
y <- factor(c("A", "B", "A", "B"))
ma <- initialization(y)
x <- c(10, 12, 11, 14)
ma <- add_analysis(ma, "ExampleVar", x)
}
\keyword{utilities}
