\name{svycollinear}
\alias{svycollinear}

\title{
Condition indexes and variance decompositions in linear models fitted with complex survey data
}
\description{
Compute condition indexes and variance decompositions for diagnosing collinearity in fixed effects, linear regression models fitted with data collected from one- and two-stage complex survey designs.
}
\usage{
svycollinear(mod, intcpt=TRUE, w, Vcov, sc=TRUE, svyglm.obj, rnd=3, fuzz=0.3)
}

\arguments{
  \item{mod}{
 Either (i) an \eqn{n \times p} matrix of real-valued covariates used in fitting a linear regression; \eqn{n} = number of observations, \eqn{p} = number of covariates in model, excluding the intercept; the matrix \code{mod} should not contain columns for the strata and cluster identifiers (unless those variables are part of the model). No missing values are allowed. The matrix can be created with \code{model.matrix}; see the example below. The \code{model} list component of an \code{svyglm} object cannot be used as it does not have factor variables coded as zero-ones.

Or, (ii) a model object produced by \code{svyglm} in the \code{survey} package.
}
  \item{intcpt}{
  \code{TRUE} if the model contains an intercept; \code{FALSE} if not.
}
  \item{w}{
 \eqn{n}-vector of survey weights used in fitting the model. No missing values are allowed.
}
  \item{Vcov}{
 Variance-covariance matrix of the estimated slopes in the regression model; component \code{cov.unscaled} in an \code{svyglm} object. This matrix can also be extracted with \code{vcov(m)} where \code{m} is an \code{svyglm} object.
}
  \item{sc}{
    \code{TRUE} if the columns of the weighted model matrix \eqn{\tilde{\mathbf{X}}} (defined in Details) should be scaled for computing condition indexes; \code{FALSE} if not. If \code{TRUE}, each column of  \eqn{\tilde{\mathbf{X}}} is divided by its Euclidean norm, \eqn{\sqrt{\tilde{\mathbf{x}}^T \tilde{\mathbf{x}}}}.
}
  \item{svyglm.obj}{
    Is \code{mod} a \code{svyglm.obj} object? \code{TRUE} or \code{FALSE}. If \code{mod} is a matrix, set \code{svyglm.obj} to \code{FALSE}.
}
  \item{rnd}{
    Round the output to \code{rnd} decimal places.
}
  \item{fuzz}{
    Replace any variance decomposition proportions that are less than \code{fuzz} by '.' in the output.
}
}
\details{
\code{svycollinear} computes condition indexes and variance decomposition proportions to use for diagnosing collinearity in a linear model fitted from complex survey data as discussed in Liao and Valliant (2012). All measures are based on \eqn{\widetilde{\mathbf{X}} = \mathbf{W}^{1/2}\mathbf{X}} where \eqn{\mathbf{W}} is the diagonal matrix of survey weights and \bold{X} is the \eqn{n \times p} matrix of covariates. In a full-rank model with \emph{p} covariates, there are \emph{p} condition indexes, defined as the ratio of the maximum eigenvalue of \eqn{\widetilde{\mathbf{X}}} to its minimum eigenvalue. Before computing condition indexes, as recommended by Belsley (1991), the columns are normalized by their individual Euclidean norms, \eqn{\sqrt{\tilde{\mathbf{x}}^T\tilde{\mathbf{x}}}}, so that each column has unit length. The columns are not centered around their means because that can obscure near-dependencies between the intercept and other covariates (Belsley 1984).

Variance decompositions are for the variance of each estimated regression coefficient and are based on a singular value decomposition of the variance formula. Proportions of the model variance, \eqn{Var_M(\hat{\mathbf{\beta}}_k)}, associated with each column of \eqn{\widetilde{\mathbf{X}}} are displayed in an output matrix described below.
}

\value{
 \eqn{p \times (p+1)} data frame, \eqn{\mathbf{\Pi}}. The first column gives the condition indexes of \eqn{\widetilde{\mathbf{X}}}. Values of 10 or more are usually considered to potentially signal collinearity of two or more columns of \eqn{\widetilde{\mathbf{X}}}. The remaining columns give the proportions (within columns) of variance of each estimated regression coefficient associated with a singular value decomposition into \emph{p} terms. Columns \eqn{2, \ldots, p+1} will each approximately sum to 1. Note that some `proportions' can be negative due to the nature of the variance decomposition. If two proportions in a given row of \eqn{\mathbf{\Pi}} are relatively large and its associated condition index in that row in the first column of \eqn{\mathbf{\Pi}} is also large, then near dependencies between the covariates associated with those elements are influencing the regression coefficient estimates.
}

\references{
Belsley, D.A., Kuh, E. and Welsch, R.E. (1980). \emph{Regression Diagnostics: Identifying Influential Data and Sources of Collinearity}. New York: Wiley-Interscience.

Belsley, D.A. (1984). Demeaning conditioning diagnostics through centering. \emph{The American Statistician}, 38(2), 73-77.

Belsley, D.A. (1991). \emph{Conditioning Diagnostics, Collinearity, and Weak Data in Regression}. New York: John Wiley & Sons, Inc.

Liao, D, and Valliant, R. (2012). Condition indexes and variance decompositions for diagnosing
collinearity in linear model analysis of survey data. \emph{Survey Methodology}, 38, 189-202.

Lumley, T. (2010). \emph{Complex Surveys}. New York: John Wiley & Sons.

Lumley, T. (2023). survey: analysis of complex survey samples. R package version 4.2.
}

\author{
Richard Valliant
}
\seealso{
\code{\link{svyvif}}
}

\examples{
require(survey)
    # example from svyglm help page
data(api)
dstrat <- svydesign(id=~1,strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
m1 <- svyglm(api00~ell+meals+mobility, design=dstrat)
    # send model object from svyglm
CI.out <- svycollinear(mod = m1, w=apistrat$pw, Vcov=vcov(m1), sc=TRUE, svyglm.obj=TRUE,
        rnd=3, fuzz= 0.3)

    # use model.matrix to create matrix of covariates in model
data(nhanes2007)
newPSU <- paste(nhanes2007$SDMVSTRA, nhanes2007$SDMVPSU, sep=".")
nhanes.dsgn <- svydesign(ids = ~newPSU,
                         strata = NULL,
                         weights = ~WTDRD1, data=nhanes2007)
m1 <- svyglm(BMXWT ~ RIDAGEYR + as.factor(RIDRETH1) + DR1TKCAL +
    DR1TTFAT + DR1TMFAT + DR1TSFAT + DR1TPFAT, design=nhanes.dsgn)
X <- model.matrix(~ RIDAGEYR + as.factor(RIDRETH1) + DR1TKCAL + DR1TTFAT + DR1TMFAT
            + DR1TSFAT + DR1TPFAT,
        data = data.frame(nhanes2007))
CI.out <- svycollinear(mod = X, w=nhanes2007$WTDRD1, Vcov=vcov(m1), sc=TRUE, svyglm.obj=FALSE,
        rnd=2, fuzz=0.3)
}

\keyword{methods}
\keyword{survey}
