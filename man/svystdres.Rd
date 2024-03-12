\name{svystdres}
\alias{svystdres}

\title{
Standardized residuals for models fitted with complex survey data
}
\description{
Compute standardized residuals for fixed effects, linear regression models fitted with data collected from one- and two-stage complex survey designs.
}
\usage{
svystdres(mobj, stvar=NULL, clvar=NULL, doplot=FALSE)
}

\arguments{
  \item{mobj}{
 model object produced by \code{svyglm} in the \code{survey} package
}
  \item{stvar}{
 name of the stratification variable in the \code{svydesign} object used to fit the model
}
  \item{clvar}{
 name of the cluster variable in the \code{svydesign} object used to fit the model}
  \item{doplot}{
 if \code{TRUE}, plot the standardized residuals vs. their sequence number in data set. Reference lines are drawn at +/-3}
}
\details{
\code{svystdres} computes the standardized residuals, i.e., the residuals divided by an estimate of the model standard deviation of the residuals. Residuals are used from a model object created by \code{svyglm} in the R \code{survey} package. The output is a vector of the standardized residuals and a scatterplot of them versus the sequence number of the sample element used in fitting the model.  By default, \code{svyglm} uses only complete cases (i.e., ones for which the dependent variable and all independent variables are non-missing) to fit the model. The rows of the data frame used in fitting the model can be retrieved from the \code{svyglm} object via \code{as.numeric(names(mobj$y))}. The data for those rows is in \code{mobj$data}.
}

\value{
  List object with values:
  \item{stdresids}{Numeric vector whose names are the rows of the data frame in the \code{svydesign} object that were used in fitting the model}
  \item{n}{number of sample clusters}
  \item{mbar}{average number of non-missing, sample elements per cluster}
  \item{rtsighat}{estimate of the square root of the model variance of the residuals, \eqn{\sqrt(\sigma^2)}}
  \item{rhohat}{estimate of the intracluster correlation of the residuals, \eqn{\rho}}
}

\references{
Li, J., and Valliant, R. (2011). Linear regression diagnostics for unclustered survey data. \emph{Journal of Official Statistics}, 27, 99-119.

Li, J., and Valliant, R. (2015). Linear regression diagnostics in cluster samples. \emph{Journal of Official Statistics}, 31, 61-75.

Lumley, T. (2010). \emph{Complex Surveys}. New York: John Wiley & Sons.

Lumley, T. (2023). survey: analysis of complex survey samples. R package version 4.2.
}

\author{
Richard Valliant
}
\seealso{
\code{\link{svyhat}}, \code{\link{svyCooksD}}
}

\examples{
require(survey)
data(api)
    # unstratified design single stage design
d0 <- svydesign(id=~1,strata=NULL, weights=~pw, data=apistrat)
m0 <- svyglm(api00 ~ ell + meals + mobility, design=d0)
svystdres(mobj=m0, stvar=NULL, clvar=NULL)

    # stratified cluster design
require(NHANES)
data(NHANESraw)
dnhanes <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTINT2YR, nest=TRUE, data=NHANESraw)
m1 <- svyglm(BPDiaAve ~ as.factor(Race1) + BMI + AlcoholYear, design = dnhanes)
svystdres(mobj=m1, stvar= "SDMVSTRA", clvar="SDMVPSU")
}

\keyword{methods}
\keyword{survey} 