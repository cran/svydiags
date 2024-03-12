\name{svyhat}
\alias{svyhat}

\title{
Leverages for models fitted with complex survey data
}
\description{
Compute leverages for fixed effects, linear regression models fitted from complex survey data.
}
\usage{
svyhat(mobj, doplot=FALSE)
}

\arguments{
  \item{mobj}{
 model object produced by \code{svyglm} in the \code{survey} package
 }
  \item{doplot}{
 if \code{TRUE}, plot the standardized residuals vs. their sequence number in data set. A reference line is drawn at 3 times the mean leverage
 }
}
\details{
\code{svyhat} computes the leverages from a model fitted with complex survey data. The model object \code{mobj} must be created by \code{svyglm} in the R \code{survey} package. The output is a vector of the leverages and a scatterplot of them versus the sequence number of the sample element used in fitting the model.  By default, \code{svyglm} uses only complete cases (i.e., ones for which the dependent variable and all independent variables are non-missing) to fit the model. The rows of the data frame used in fitting the model can be retrieved from the \code{svyglm} object via \code{as.numeric(names(mobj$y))}. The data for those rows is in \code{mobj$data}.
}

\value{
 Numeric vector whose names are the rows of the data frame in the \code{svydesign} object that were used in fitting the model.
}

\references{
Belsley, D.A., Kuh, E. and Welsch, R. (1980). \emph{Regression Diagnostics: Identifying Influential Data and Sources of Collinearity}. New York: John Wiley & Sons, Inc.

Li, J., and Valliant, R. (2009). Survey weighted hat matrix and leverages. \emph{Survey Methodology}, 35, 15-24.

Lumley, T. (2010). \emph{Complex Surveys}. New York: John Wiley & Sons.

Lumley, T. (2023). survey: analysis of complex survey samples. R package version 4.2.
}

\author{
Richard Valliant
}
\seealso{
\code{\link{svystdres}}
}
\examples{
require(survey)
data(api)
dstrat <- svydesign(id=~1,strata=~stype, weights=~pw, data=apistrat)
m1 <- svyglm(api00 ~ ell + meals + mobility, design=dstrat)
h <- svyhat(mobj = m1, doplot=TRUE)
100*sum(h > 3*mean(h))/length(h)    # percentage of leverages > 3*mean

require(NHANES)
data(NHANESraw)
dnhanes <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTINT2YR, nest=TRUE, data=NHANESraw)
m1 <- svyglm(BPDiaAve ~ as.factor(Race1) + BMI + AlcoholYear, design = dnhanes)
h <- svyhat(mobj = m1, doplot=TRUE)
}

\keyword{methods}
\keyword{survey} 