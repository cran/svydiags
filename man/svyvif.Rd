\name{svyvif}
\alias{svyvif}

\title{
Variance inflation factors (VIF) for general linear models fitted with complex survey data
}
\description{
Compute a VIF for fixed effects, general linear regression models fitted with data collected from one- and two-stage complex survey designs.
}
\usage{
svyvif(mobj, X, w, stvar=NULL, clvar=NULL)
}

\arguments{
  \item{mobj}{
 model object produced by \code{svyglm}. The following families of models are allowed: binomial, gaussian, poisson, quasibinomial, and quasipoisson. Other families allowed by \code{svyglm} will produce an error in \code{svyvif}.
}
  \item{X}{
 \eqn{n \times p} matrix of real-valued covariates used in fitting a linear regression; \eqn{n} = number of observations, \eqn{p} = number of covariates in model, excluding the intercept. A column of 1's for an intercept should not be included. \code{X} should not contain columns for the strata and cluster identifiers (unless those variables are part of the model). No missing values are allowed.
}
  \item{w}{
 \eqn{n}-vector of survey weights used in fitting the model. No missing values are allowed.
}
  \item{stvar}{
 field in \code{mobj} that contains the stratum variable in the complex sample design; use \code{stvar = NULL} if there are no strata
}
  \item{clvar}{
 field in \code{mobj} that contains the cluster variable in the complex sample design; use \code{clvar = NULL} if there are no clusters
}
}
\details{
\code{svyvif} computes a variance inflation factor (VIF) appropriate for linear models and some general linear models (GLMs) fitted from complex survey data (see Liao & Valliant 2012). A VIF measures the inflation of a slope estimate caused by nonorthogonality of the predictors over and above what the variance would be with orthogonality (Theil 1971; Belsley, Kuh, and Welsch 1980). The standard VIF equals \eqn{1/(1 - R^2_k)} where \eqn{R_k} is the multiple correlation of the \eqn{k^{th}} column of \code{X} regressed on the remaining columns. The complex sample value of the VIF for a linear model consists of the standard VIF multiplied by two adjustments denoted in the output as \code{zeta} and \code{varrho}. The VIF for a GLM is similar (Liao 2010, chap. 5). There is no widely agreed-upon cutoff value for identifying high values of a VIF, although 10 is a common suggestion.
}

\value{
 \eqn{p \times 5} matrix with columns:
  \item{svy.vif}{complex sample VIF}
  \item{reg.vif}{standard VIF, \eqn{1/(1 - R^2_k)}, that omits the factors, \code{zeta} and \code{varrho}; \eqn{R^2_k} is an R-square from a weighted least squares regression of the \eqn{k^{th}} x on the other x's in the regression}
  \item{zeta}{1st multiplicative adjustment to \code{reg.vif}}
  \item{varrho}{2nd multiplicative adjustment to \code{reg.vif}}
  \item{zeta.x.varrho}{product of the two adjustments to \code{reg.vif}}
  \item{R.square}{R-square in the regression of the \eqn{k^{th}} x on the other x's, including the intercept}
}

\references{
Belsley, D.A., Kuh, E. and Welsch, R.E. (1980). \emph{Regression Diagnostics: Identifying Influential Data and Sources of Collinearity}. New York: Wiley-Interscience.

Liao, D. (2010). Collinearity Diagnostics for Complex Survey Data. PhD thesis, University of Maryland. \url{ http://hdl.handle.net/1903/10881}.

Liao, D, and Valliant, R. (2012). Variance inflation factors in the analysis of complex survey data. \emph{Survey Methodology}, 38, 53-62.

Theil, H. (1971). \emph{Principles of Econometrics}. New York: John Wiley & Sons, Inc.

Lumley, T. (2010). \emph{Complex Surveys}. New York: John Wiley & Sons.

Lumley, T. (2023). survey: analysis of complex survey samples. R package version 4.2.
}

\author{
Richard Valliant
}
\seealso{
\code{\link{Vmat}}
}

\examples{
require(survey)
data(nhanes2007)
X1 <- nhanes2007[order(nhanes2007$SDMVSTRA, nhanes2007$SDMVPSU),]
    # eliminate cases with missing values
delete <- which(complete.cases(X1)==FALSE)
X2 <- X1[-delete,]
nhanes.dsgn <- svydesign(ids = ~SDMVPSU,
                         strata = ~SDMVSTRA,
                         weights = ~WTDRD1, nest=TRUE, data=X2)
    # linear model
m1 <- svyglm(BMXWT ~ RIDAGEYR + as.factor(RIDRETH1) + DR1TKCAL
            + DR1TTFAT + DR1TMFAT, design=nhanes.dsgn)
summary(m1)
    # construct X matrix using model.matrix from stats package
X3 <- model.matrix(~ RIDAGEYR + as.factor(RIDRETH1) + DR1TKCAL + DR1TTFAT + DR1TMFAT,
        data = data.frame(X2))
    # remove col of 1's for intercept with X3[,-1]
svyvif(mobj=m1, X=X3[,-1], w = X2$WTDRD1, stvar=NULL, clvar=NULL)

    # Logistic model
X2$obese <- X2$BMXBMI >= 30
nhanes.dsgn <- svydesign(ids = ~SDMVPSU,
                         strata = ~SDMVSTRA,
                         weights = ~WTDRD1, nest=TRUE, data=X2)
m2 <- svyglm(obese ~ RIDAGEYR + as.factor(RIDRETH1) + DR1TKCAL
             + DR1TTFAT + DR1TMFAT, design=nhanes.dsgn, family="quasibinomial")
summary(m2)
svyvif(mobj=m2, X=X3[,-1], w = X2$WTDRD1, stvar = "SDMVSTRA", clvar = "SDMVPSU")
}

\keyword{methods}
\keyword{survey} 