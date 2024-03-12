\name{svyCooksD}
\alias{svyCooksD}

\title{
Modified Cook's D for models fitted with complex survey data
}
\description{
Compute a modified Cook's D for fixed effects, linear regression models fitted with data collected from one- and two-stage complex survey designs.
}
\usage{
svyCooksD(mobj, stvar=NULL, clvar=NULL, doplot=FALSE)
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
 if \code{TRUE}, plot the modified Cook's D values vs. their sequence number in data set. Reference lines are drawn at 2 and 3}
}
\details{
\code{svyCooksD} computes the modified Cook's D (m-cook; see Atkinson (1982) and Li & Valliant (2011, 2015)) which measures the effect on the vector of parameter estimates of deleting single observations when fitting a fixed effects regression model to complex survey data. The function \code{svystdres} is called for some of the calculations. Values of m-cook are considered large if they are greater than 2 or 3.  The R package \code{MASS} must also be loaded before calling \code{svyCooksD}. The output is a vector of the m-cook values and a scatterplot of them versus the sequence number of the sample element used in fitting the model.  By default, \code{svyglm} uses only complete cases (i.e., ones for which the dependent variable and all independent variables are non-missing) to fit the model. The rows of the data frame used in fitting the model can be retrieved from the \code{svyglm} object via \code{as.numeric(names(mobj$y))}. The data for those rows is in \code{mobj$data}.
}

\value{
Numeric vector whose names are the rows of the data frame in the \code{svydesign} object that were used in fitting the model
}

\references{
Atkinson, A.C. (1982). Regression diagnostics, transformations and constructed variables (with discussion). \emph{Journal of the Royal Statistical Society}, Series B, Methodological, 44, 1-36.

Cook, R.D. (1977). Detection of Influential Observation in Linear Regression. \emph{Technometrics}, 19, 15-18.

Cook, R.D. and Weisberg, S. (1982). \emph{Residuals and Influence in Regression}. London:Chapman & Hall Ltd.

Li, J., and Valliant, R. (2011). Linear regression diagnostics for unclustered survey data. \emph{Journal of Official Statistics}, 27, 99-119.

Li, J., and Valliant, R. (2015). Linear regression diagnostics in cluster samples. \emph{Journal of Official Statistics}, 31, 61-75.

Lumley, T. (2010). \emph{Complex Surveys}. New York: John Wiley & Sons.

Lumley, T. (2023). survey: analysis of complex survey samples. R package version 4.2.
}

\author{
Richard Valliant
}
\seealso{
\code{\link{svydfbetas}}, \code{\link{svydffits}}, \code{\link{svystdres}}
}

\examples{
require(MASS)   # to get ginv
require(survey)
data(api)
    # unstratified design single stage design
d0 <- svydesign(id=~1,strata=NULL, weights=~pw, data=apistrat)
m0 <- svyglm(api00 ~ ell + meals + mobility, design=d0)
mcook <- svyCooksD(m0, doplot=TRUE)

    # stratified clustered design
require(NHANES)
data(NHANESraw)
dnhanes <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTINT2YR, nest=TRUE, data=NHANESraw)
m2 <- svyglm(BPDiaAve ~ as.factor(Race1) + BMI + AlcoholYear, design = dnhanes)
mcook <- svyCooksD(mobj=m2, stvar="SDMVSTRA", clvar="SDMVPSU", doplot=TRUE)
}

\keyword{methods}
\keyword{survey} 