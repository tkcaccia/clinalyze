\name{multi_analysis}

\alias{multi_analysis}

\title{Continuous Information}

\description{Summarization of the continuous information.}

\usage{
multi_analysis  (data, 
                 y, 
                 FUN=c("continuous.test","correlation.test"), ...)                 
}

\arguments{
  \item{data}{the matrix containing the continuous values. Each row corresponds to a different sample. Each column corresponds to a different variable.}
  \item{y}{the classification of the cohort.}
  \item{FUN}{function to be considered. Choices are "\code{continuous.test}" and "\code{correlation.test}"}
  \item{...}{further arguments to be passed to or from methods.}
}

\value{
    The function returns a table with the summarized information. If the number of group is equal to two, the p-value is computed using the Wilcoxon rank-sum test, Kruskal-Wallis test otherwise.
  }

\author{Stefano Cacciatore}

  \seealso{\code{\link{categorical.test}},\code{\link{continuous.test}},\code{\link{correlation.test}}, \code{\link{txtsummary}}}

\examples{
data(clinical)


multi_analysis(clinical[,c("BMI","Age")],clinical[,"Hospital"],FUN="continuous.test")

}
\keyword{multi_analysis}
