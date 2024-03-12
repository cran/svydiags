---
output:
  html_document: default
  pdf_document: default
---
# Changes and Updates for svydiags package

# svydiags 0.5

*    svycollinear: Updated to include an error trap on the mod and svyglm.obj parameters. If mod is a model matrix, then svyglm.obj must be FALSE.

*    svycollinear: Check added when mod is not a matrix but is a data frame. If any column is a factor, then an error is issued that the mod must have factors recoded as zeroes and ones using model.matrix.
