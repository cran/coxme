# Automatically generated from all.nw using noweb
# The objects that do the actual work (not much work)
fixef.coxme <- function(object, ...)
    object$coefficients$fixed

ranef.coxme <- function(object, ...)
    object$coefficients$random
