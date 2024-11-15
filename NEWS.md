---
output:
  html_document: default
  pdf_document: default
---
# Changes and Updates for svydiags package

# svydiags 0.7

*    svycollinear updated to include more families: binomial and quasibinomial (logit and probit links), Gamma (inverse link),
     inverse.guassian ("1/mu^2" link), poisson and quasipoisson (log link). The function has also been redesigned to require 
     the X matrix as an input, which is comparable to the requirements for svyvif.
     The unit variance-covariance matrix is now computed within the function based on information in the input model object (mobj).

*    svyvif updated to include Gamma and inverse.gaussian families

# svydiags 0.6

*    Intercept-excluded version of VIF added for all models covered: binomial, gaussian, poisson, quasibinomial, and quasipoisson. The intercept-excluded version is analogous to the one found in other statistical packages for linear models fitted with non-survey data.

*    Error corrected in calculation of R-squared in intercept-included version of VIF for non-gaussian models (i.e., binomial, poisson, quasibinomial, and quasipoisson). In previous versions SST in the denominator was not corrected for the mean.

# svydiags 0.5

*    svycollinear: Updated to include an error trap on the mod and svyglm.obj parameters. If mod is a model matrix, then svyglm.obj must be FALSE.

*    svycollinear: Check added when mod is not a matrix but is a data frame. If any column is a factor, then an error is issued that the mod must have factors recoded as zeroes and ones using model.matrix.
