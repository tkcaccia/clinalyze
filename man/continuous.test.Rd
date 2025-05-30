
\name{continuous.test}

\alias{continuous.test}

\title{Summary and Statistical Test for Continuous Variables}

\description{
Summarizes a continuous variable across groups and performs a statistical test to compare distributions between them.
}

\usage{
continuous.test(name,
                x,    
                y,
                center = c("median", "mean"),
                digits = 3,
                scientific = FALSE, 
                range = c("IQR", "95%CI", "range", "sd"),
                logchange = FALSE,
                pos = 2, 
                method = c("non-parametric", "parametric"),
                total.column = FALSE, 
                ...)
}

\arguments{
  \item{name}{A character string specifying the name of the continuous variable.}
  \item{x}{A numeric vector containing the continuous values to summarize.}
  \item{y}{A factor or character vector specifying the groupings or classifications.}
  \item{center}{A character string specifying the measure of central tendency to report: either \code{"median"} or \code{"mean"}.}
  \item{digits}{An integer indicating the number of digits to display.}
  \item{scientific}{Logical value indicating whether to use scientific notation for numbers.}
  \item{range}{A character string specifying the type of dispersion measure to display: \code{"IQR"}, \code{"95\%CI"}, \code{"range"}, or \code{"sd"}.}
  \item{logchange}{Logical value indicating whether to include the log2 fold change between groups (if applicable).}
  \item{pos}{Integer value indicating the layout of the range: \code{1} for displaying in columns, \code{2} for rows.}
  \item{method}{A character string specifying the statistical test to apply: \code{"non-parametric"} (default) or \code{"parametric"}.}
  \item{total.column}{Logical value indicating whether to include a column summarizing the total population (default is \code{FALSE}).}
  \item{...}{Additional arguments passed to internal functions.}
}

\value{
A data frame summarizing the continuous variable for each group with the selected central tendency and range, as well as the corresponding p-value.

For \code{method = "non-parametric"}:
\itemize{
  \item If two groups are present, the Wilcoxon rank-sum test is used.
  \item If more than two groups, the Kruskal-Wallis test is applied.
}

For \code{method = "parametric"}:
\itemize{
  \item If two groups are present, the Student's t-test is used.
  \item If more than two groups, one-way ANOVA is applied.
}
}

\author{Stefano Cacciatore}

\seealso{\code{\link{correlation.test}}, \code{\link{categorical.test}}, \code{\link{txtsummary}}}

\examples{
data(prostate)

hosp <- clinical[,"Hospital"]
gender <- clinical[,"Gender"]
GS <- clinical[,"Gleason score"]
BMI <- clinical[,"BMI"]
age <- clinical[,"Age"]

A <- categorical.test("Gender", gender, hosp)
B <- categorical.test("Gleason score", GS, hosp)

C <- continuous.test("BMI", BMI, hosp, digits = 2)
D <- continuous.test("Age", age, hosp, digits = 1)

rbind(A, B, C, D)
}

\keyword{univariate}
\keyword{summary}
