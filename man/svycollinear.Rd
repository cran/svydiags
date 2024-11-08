\name{svycollinear}
\alias{svycollinear}

\title{
Condition indexes and variance decompositions in general linear models (GLMs) fitted with complex survey data
}
\description{
Compute condition indexes and variance decompositions for diagnosing collinearity in fixed effects, general linear regression models fitted with data collected from one- and two-stage complex survey designs.
}
\usage{
svycollinear(mobj, X, w, sc=TRUE, rnd=3, fuzz=0.05)
}

\arguments{
  \item{mobj}{
   model object produced by \code{svyglm}. The following families of models are allowed: \code{binomial} and \code{quasibinomial} (\code{logit} and \code{probit} links), \code{gaussian} (\code{identity} link), \code{poisson} and \code{quasipoisson} (\code{log} link), \code{Gamma} (\code{inverse} link), and \code{inverse.gaussian} (\code{1/mu^2} link). Other families or links allowed by \code{svyglm} will produce an error in \code{svycollinear}.
}
  \item{X}{
  \eqn{n \times p} matrix of real-valued covariates used in fitting the regression; \eqn{n} = number of observations, \eqn{p} = number of covariates in model, excluding the intercept. A column of 1's for an intercept may be included if the model includes an intercept. \code{X} is most easily produced by the function \code{model.matrix} in the \code{stats} package, which will correctly code factors as 0-1. \code{X} should not contain columns for the strata and cluster identifiers (unless those variables are part of the model). No missing values are allowed.
}
  \item{w}{
  \eqn{n}-vector of survey weights used in fitting the model. No missing values are allowed.
}
  \item{sc}{
    \code{TRUE} if the columns of the weighted model matrix \eqn{\tilde{\mathbf{X}}} (defined in Details) should be scaled for computing condition indexes; \code{FALSE} if not. If \code{TRUE}, each column of  \eqn{\tilde{\mathbf{X}}} is divided by its Euclidean norm, \eqn{\sqrt{\tilde{\mathbf{x}}^T \tilde{\mathbf{x}}}}.
}
  \item{rnd}{
    Round the output to \code{rnd} decimal places.
}
  \item{fuzz}{
    Replace any variance decomposition proportions that are less than \code{fuzz} by `.' in the output.
}
}
\details{
\code{svycollinear} computes condition indexes and variance decomposition proportions to use for diagnosing collinearity in a general linear model fitted from complex survey data as discussed in Liao (2010, ch. 5) and Liao and Valliant (2012). All measures are based on \eqn{\widetilde{\mathbf{X}} = \mathbf{W}^{1/2}\hat{\mathbf{\Gamma}}\mathbf{X}} where \eqn{\mathbf{W}} is the diagonal matrix of survey weights, \eqn{\hat{\mathbf{\Gamma}}} is a diagonal matrix of estimated parameters from the particular type of GLM,  and \bold{X} is the \eqn{n \times p} matrix of covariates. In a full-rank model with \emph{p} covariates, there are \emph{p} condition indexes, defined as the ratio of the maximum eigenvalue of \eqn{\widetilde{\mathbf{X}}} to each of the \emph{p} eigenvalues. If \code{sc=TRUE}, before computing condition indexes, as recommended by Belsley (1991), the columns are normalized by their individual Euclidean norms, \eqn{\sqrt{\tilde{\mathbf{x}}^T\tilde{\mathbf{x}}}}, so that each column has unit length. The columns are not centered around their means because that can obscure near-dependencies between the intercept and other covariates (Belsley 1984).

Variance decompositions are for the variance of each estimated regression coefficient and are based on a singular value decomposition of the variance formula. For linear models, the decomposition is for the sandwich variance estimator, which has both a model-based and design-based interpretation. In the case of nonlinear GLMs (i.e., \code{family} is not \code{gaussian}), the variance is the approximate model variance. Proportions of the model variance, \eqn{Var_M(\hat{\mathbf{\beta}}_k)}, associated with each column of \eqn{\widetilde{\mathbf{X}}} are displayed in an output matrix described below.
}

\value{
 \eqn{p \times (p+1)} data frame, \eqn{\mathbf{\Pi}}. The first column gives the condition indexes of \eqn{\widetilde{\mathbf{X}}}. Values of 10 or more are usually considered to potentially signal collinearity of two or more columns of \eqn{\widetilde{\mathbf{X}}}. The remaining columns give the proportions (within columns) of variance of each estimated regression coefficient associated with a singular value decomposition into \emph{p} terms. Columns \eqn{2, \ldots, p+1} will each approximately sum to 1. When \code{family=gaussian}, some `proportions' can be negative or greater than 1 due to the nature of the variance decomposition (see Liao and Valliant, 2012). For other families the proportions will be in [0,1]. If two proportions in a given row of \eqn{\mathbf{\Pi}} are relatively large and its associated condition index in that row in the first column of \eqn{\mathbf{\Pi}} is also large, then near dependencies between the covariates associated with those elements are influencing the regression coefficient estimates.
}

\references{
Belsley, D.A., Kuh, E. and Welsch, R.E. (1980). \emph{Regression Diagnostics: Identifying Influential Data and Sources of Collinearity}. New York: Wiley-Interscience.

Belsley, D.A. (1984). Demeaning conditioning diagnostics through centering. \emph{The American Statistician}, 38(2), 73-77.

Belsley, D.A. (1991). \emph{Conditioning Diagnostics, Collinearity, and Weak Data in Regression}. New York: John Wiley & Sons, Inc.

Liao, D. (2010). Collinearity Diagnostics for Complex Survey Data. PhD thesis, University of Maryland. \url{ http://hdl.handle.net/1903/10881}.

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
    # linear model
m1 <- svyglm(api00 ~ ell + meals + mobility, design=dstrat)
X.model <- model.matrix(~ ell + meals + mobility, data = apistrat)
    # send model object from svyglm
svycollinear(mobj=m1, X=X.model, w=apistrat$pw, sc=TRUE, rnd=3, fuzz= 0.05)

    # logistic model
data(nhanes2007)
nhanes2007$obese <- nhanes2007$BMXBMI >= 30
nhanes.dsgn <- svydesign(ids = ~SDMVPSU,
                         strata = ~SDMVSTRA,
                         weights = ~WTDRD1, nest=TRUE, data=nhanes2007)
m2 <- svyglm(obese ~ RIDAGEYR + as.factor(RIDRETH1) + DR1TKCAL +
    DR1TTFAT + DR1TMFAT, design=nhanes.dsgn, family=quasibinomial())
X.model <- model.matrix(~ RIDAGEYR + as.factor(RIDRETH1) + DR1TKCAL + DR1TTFAT + DR1TMFAT,
        data = data.frame(nhanes2007))
svycollinear(mobj=m2, X=X.model, w=nhanes2007$WTDRD1, sc=TRUE, rnd=2, fuzz=0.05)
}

\keyword{methods}
\keyword{survey}
