\name{Vmat}
\alias{Vmat}

\title{
Compute covariance matrix of residuals for general linear models fitted with complex survey data
}
\description{
Compute a covariance matrix using residuals from a fixed effects, general linear regression model fitted with data collected from one- and two-stage complex survey designs.
}
\usage{
Vmat(mobj, stvar = NULL, clvar = NULL)
}

\arguments{
  \item{mobj}{
    model object produced by \code{svyglm}
}
  \item{stvar}{
 field in \code{mobj} that contains the stratum variable in the complex sample design; use \code{stvar = NULL} if there are no strata
}
  \item{clvar}{
 field in \code{mobj} that contains the cluster variable in the complex sample design; use \code{clvar = NULL} if there are no clusters
}
}
\details{
\code{Vmat} computes a covariance matrix among the residuals returned from \code{svyglm} in the \code{survey} package. \code{Vmat} is called by \code{svyvif} when computing variance inflation factors. The matrix that is computed by \code{Vmat} is appropriate under these model assumptions: (1) in single-stage, unclustered sampling, units are assumed to be uncorrelated but can have different model variances, (2) in single-stage, stratified sampling, units are assumed to be uncorrelated within strata and between strata but can have different model variances; (3) in unstratified, clustered samples, units in different clusters are assumed to be uncorrelated but units within clusters are correlated; (3) in stratified, clustered samples, units in different strata or clusters are assumed to be uncorrelated but units within clusters are correlated.
}

\value{
 \eqn{n \times n} matrix where \eqn{n} is the number of cases used in the linear regression model
}

\references{
Liao, D, and Valliant, R. (2012). Variance inflation factors in the analysis of complex survey data. \emph{Survey Methodology}, 38, 53-62.

Lumley, T. (2010). \emph{Complex Surveys}. New York: John Wiley & Sons.

Lumley, T. (2021). survey: analysis of complex survey samples. R package version 4.1-1.
}

\author{
Richard Valliant
}
\seealso{
\code{\link{svyvif}}
}

\examples{
require(Matrix)
require(survey)
data(nhanes2007)
black <- nhanes2007$RIDRETH1 == 4
X <- nhanes2007
X <-  cbind(X, black)
X1 <- X[order(X$SDMVSTRA, X$SDMVPSU),]

    # unstratified, unclustered design
nhanes.dsgn <- svydesign(ids = 1:nrow(X1),
                         strata = NULL,
                         weights = ~WTDRD1, data=X1)
m1 <- svyglm(BMXWT ~ RIDAGEYR + as.factor(black) + DR1TKCAL, design=nhanes.dsgn)
summary(m1)

V <- Vmat(mobj = m1,
          stvar = NULL,
          clvar = NULL)

    # stratified, clustered design
nhanes.dsgn <- svydesign(ids = ~SDMVPSU,
                         strata = ~SDMVSTRA,
                         weights = ~WTDRD1, nest=TRUE, data=X1)
m1 <- svyglm(BMXWT ~ RIDAGEYR + as.factor(black) + DR1TKCAL, design=nhanes.dsgn)
summary(m1)
V <- Vmat(mobj = m1,
          stvar = "SDMVSTRA",
          clvar = "SDMVPSU")
}

\keyword{methods}
\keyword{survey} 