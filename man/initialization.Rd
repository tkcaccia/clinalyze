\name{initialization}
\alias{initialization}
\title{Initialize a Clinical Table Object}
\description{
Initializes a new object of class \code{"clinical.table"}, setting up the grouping variable and a flag for whether to include a total column.
}
\usage{
initialization(y, total.column = FALSE)
}
\arguments{
  \item{y}{A numeric or factor vector representing the outcome or grouping variable.}
  \item{total.column}{Logical; whether to include a total column in the results (default is \code{FALSE}).}
}
\value{
Returns an S4 object of class \code{"clinical.table"}.
}
\examples{
y <- factor(c("A", "B", "A", "B"))
ma <- initialization(y)
}
\keyword{utilities}
