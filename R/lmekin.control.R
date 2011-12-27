#
# Gather all of the control parameters for lmekin into one spot
#
lmekin.control <- function(optpar=list(method='BFGS', 
                                      control=list(reltol=1e-8))) {
    list(optpar=optpar)
    }
