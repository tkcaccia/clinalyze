\name{categorical.test}

\alias{categorical.test}

\title{Categorical Information}

\description{Summarization of the categorical information.}

\usage{
categorical.test (name,x,y,total.column=FALSE,...)
}

\arguments{
  \item{name}{the name of the feature.}
  \item{x}{the information to summarize.}
  \item{y}{the classification of the cohort.}
  \item{total.column}{option to visualize the total (by default = "\code{FALSE}").}
    \item{...}{further arguments to be passed to the function.}

}

\value{
    The function returns a table with the summarized information and The p-value computated using the Fisher's test.
  }

\author{Stefano Cacciatore}

  \seealso{\code{\link{correlation.test}},\code{\link{continuous.test}}, \code{\link{txtsummary}}}

\examples{

data(prostate)

hosp=clinical[,"Hospital"]
gender=clinical[,"Gender"]
GS=clinical[,"Gleason score"]
BMI=clinical[,"BMI"]
age=clinical[,"Age"]

A=categorical.test("Gender",gender,hosp)
B=categorical.test("Gleason score",GS,hosp)

C=continuous.test("BMI",BMI,hosp,digits=2)
D=continuous.test("Age",age,hosp,digits=1)

rbind(A,B,C,D)


}
\keyword{categorical.test}
