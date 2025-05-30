\name{intersect}
\alias{intersect}
\title{Intersection of Multiple Vectors}
\description{
Returns the intersection of all vectors provided as arguments.
}
\usage{
intersect(x, y, ...)
}
\arguments{
  \item{x}{A vector.}
  \item{y}{A vector.}
  \item{...}{Additional vectors to include in the intersection.}
}
\details{
This function computes the intersection of two or more vectors by identifying elements that are present in all of them. Unlike the base R \code{\link{intersect}}, which only compares two vectors, this version can take multiple vectors as input.
}
\value{
A character vector (or \code{NULL} if no intersection exists) containing the elements that are common to all input vectors.
}
\examples{
x <- c("a", "b", "c")
y <- c("b", "c", "d")
z <- c("c", "b", "e")

intersect(x, y, z)  # Returns "b" and "c"
intersect(x, y, c("f", "g"))  # Returns NULL
}
\seealso{\code{\link{intersect}}, \code{\link{union}}, \code{\link{setdiff}}}
\author{Stefano Cacciatore}
\keyword{set operations}
