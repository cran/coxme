ranef <- function(object, ...) {
    UseMethod("ranef")
    }
random.effects <- function(object,...) {
    UseMethod("ranef")
    }

fixef <- function(object, ...) {
    UseMethod("fixef")
    }
fixed.effects <- function(object, ...) {
    UseMethod("fixef")
    }

fixef.coxme <- function(object, ...)
    object$coefficients$fixed

ranef.coxme <- function(object, ...)
    object$coefficients$random
