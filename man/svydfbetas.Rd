\name{svydfbetas}
\alias{svydfbetas}

\title{
dfbetas for models fitted with complex survey data
}
\description{
Compute the dfbetas measure of the effect of extreme observations on parameter estimates for fixed effects, linear regression models fitted with data collected from one- and two-stage complex survey designs.
}
\usage{
svydfbetas(mobj, stvar=NULL, clvar=NULL, z=3)
}

\arguments{
  \item{mobj}{
    model object produced by \code{svyglm} in the \code{survey} package
  }
  \item{stvar}{
    name of the stratification variable in the \code{svydesign} object used to fit the model
  }
  \item{clvar}{
    name of the cluster variable in the \code{svydesign} object used to fit the model
  }
  \item{z}{
    numerator of cutoff for measuring whether an observation has an extreme effect on its own predicted value; default is 3 but can be adjusted to control how many observations are flagged for inspection
  }
}
\details{
\code{svydfbetas} computes the values of dfbetas for each observation and parameter estimate, i.e., the amount that a parameter estimate changes when the unit is deleted from the sample. The model object must be created by \code{svyglm} in the R \code{survey} package. The output is a vector of the dfbeta and standardized dfbetas values.  By default, \code{svyglm} uses only complete cases (i.e., ones for which the dependent variable and all independent variables are non-missing) to fit the model. The rows of the data frame used in fitting the model can be retrieved from the \code{svyglm} object via \code{as.numeric(names(mobj$y))}. The data for those rows is in \code{mobj$data}.
}

\value{
  List object with values:
  \item{Dfbeta}{Numeric vector of unstandardized dfbeta values whose names are the rows of the data frame in the \code{svydesign} object that were used in fitting the model}
  \item{Dfbetas}{Numeric vector of standardized dfbetas values whose names are the rows of the data frame in the \code{svydesign} object that were used in fitting the model}
    \item{cutoff}{Value used for gauging whether a value of dffits is large. For a single-stage sample, \code{cutoff}=\eqn{z/\sqrt{n}}; for a 2-stage sample, \code{cutoff}=\eqn{z/\sqrt{n[1+\rho (\bar{m}-1)]}} }
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
\code{\link{svydffits}}, \code{\link{svyCooksD}}
}

\examples{
require(survey)
data(api)
    # unstratified design single stage design
d0 <- svydesign(id=~1,strata=NULL, weights=~pw, data=apistrat)
m0 <- svyglm(api00 ~ ell + meals + mobility, design=d0)
svydfbetas(mobj=m0)

    # stratified cluster
require(NHANES)
data(NHANESraw)
dnhanes <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTINT2YR, nest=TRUE, data=NHANESraw)
m2 <- svyglm(BPDiaAve ~ as.factor(Race1) + BMI + AlcoholYear, design = dnhanes)
yy <- svydfbetas(mobj=m2, stvar= "SDMVSTRA", clvar="SDMVPSU")
apply(abs(yy$Dfbetas) > yy$cutoff,1, sum)
}

\keyword{methods}
\keyword{survey} 