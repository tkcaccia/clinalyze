\name{frequency_matching}

\alias{frequency_matching}

\title{Frequency Matching}

\description{A method to select unbalanced groupd in a cohort.}

\usage{
frequency_matching (data,label,times=5,seed=1234)
}

\arguments{
  \item{data}{a data.frame of data.}
  \item{label}{a classification of the groups.}
  \item{times}{The ratio between the two groups.}
  \item{seed}{a single number for random number generation.}

}



\value{
    The function returns a list with 2 items or 4 items (if a test data set is present):
    \item{data}{the data after the frequency matching.}
    \item{label}{the label after the frequency matching.}
    \item{selection}{the rows selected for the frequency matching.}
  }

\author{Stefano Cacciatore}




\examples{


hosp=prostate[,"Hospital"]
gender=prostate[,"Gender"]
GS=prostate[,"Gleason score"]
BMI=prostate[,"BMI"]
age=prostate[,"Age"]

A=categorical.test("Gender",gender,hosp)
B=categorical.test("Gleason score",GS,hosp)

C=continuous.test("BMI",BMI,hosp,digits=2)
D=continuous.test("Age",age,hosp,digits=1)

# Analysis without matching
rbind(A,B,C,D)



# The order is important. Right is more important than left in the vector
# So, Ethnicity will be more important than Age
var=c("Age","BMI","Gleason score")
data.categorized=prostate[,var]

# Extract the Age vector
x <- data.categorized[["Age"]]

# Compute quantiles (0%, 25%, 50%, 75%, 100%) with NA handling
breaks <- quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)

# Expand the range to ensure all values are included
breaks[1] <- breaks[1] - 1
breaks[length(breaks)] <- breaks[length(breaks)] + 1

# Apply the cut and update the Age column with labeled bins
data.categorized[["Age"]] <- cut(x, breaks = breaks, include.lowest = TRUE)

# Extract the Age vector
x <- data.categorized[["BMI"]]

# Compute quantiles (0%, 25%, 50%, 75%, 100%) with NA handling
breaks <- quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)

# Expand the range to ensure all values are included
breaks[1] <- breaks[1] - 1
breaks[length(breaks)] <- breaks[length(breaks)] + 1

# Apply the cut and update the Age column with labeled bins
data.categorized[["BMI"]] <- cut(x, breaks = breaks, include.lowest = TRUE)

times=c(1,1)
names(times)=c("Hospital A","Hospital B")
t=frequency_matching(data.categorized,prostate[,"Hospital"],times=times)



newdata=prostate[t$selection,]

hosp.new=newdata[,"Hospital"]
gender.new=newdata[,"Gender"]
GS.new=newdata[,"Gleason score"]
BMI.new=newdata[,"BMI"]
age.new=newdata[,"Age"]

A=categorical.test("Gender",gender.new,hosp.new)
B=categorical.test("Gleason score",GS.new,hosp.new)

C=continuous.test("BMI",BMI.new,hosp.new,digits=2)
D=continuous.test("Age",age.new,hosp.new,digits=1)

# Analysis with matching
rbind(A,B,C,D)

}
