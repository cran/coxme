# Automatically generated from all.nw using noweb
# The objects that do the actual work (not much work)
fixef.coxme <- function(object, ...)
    object$coefficients

fixef.lmekin <- function(object, ...)
    object$coefficients$fixed

ranef.coxme <- function(object, ...)
    object$frail

ranef.lmekin <- function(object, ...)
    object$coefficients$random

VarCorr.coxme <- function(x, ...) 
    x$vcoef
    
VarCorr.lmekin <- function(x, ...) 
    x$vcoef

vcov.coxme <- function(object, ...) {
    nf <- length(fixef(object))
    indx <- seq(length=nf, to=nrow(object$var))
    as.matrix(object$var[indx, indx])
}

vcov.lmekin <- vcov.coxme    
