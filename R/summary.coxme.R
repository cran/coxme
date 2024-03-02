summary.coxme <- function(object, ...){
    beta <- object$coefficients
    nvar <- length(beta)
    nfrail<- nrow(object$var) - nvar

    loglik <- object$loglik + c(0,0, object$penalty)
    chi1 <- 2*diff(loglik[1:2]) 
    chi2 <- 2*diff(loglik[c(1,3)])
    chi <- rbind(c(round(chi1,2), round(object$df[1],2),
                    signif(1- pchisq(chi1,object$df[1]),5),
                    round(chi1- 2*object$df[1],2),
                    round(chi1- log(object$n[1])*object$df[1],2)),
                  c(round(chi2,2), round(object$df[2],2),
                    signif(1- pchisq(chi2,object$df[2]),5),
                    round(chi2- 2*object$df[2],2),
                    round(chi2- log(object$n[1])*object$df[2],2)))
    dimnames(chi) <- list(c("Integrated loglik", " Penalized loglik"),
                           c("Chisq", "df", "p", "AIC", "BIC"))

    if (length(beta) > 0)  { # Not a ~1 model
        se <- sqrt(diag(object$var)[nfrail+1:nvar])
        coef <- cbind(beta, exp(beta), se, round(beta/se,2),
                     pchisq((beta/ se)^2, 1, lower.tail=FALSE))
        dimnames(coef) <- list(names(beta), c("coef", "exp(coef)",
            "se(coef)", "z", "p"))
        }


    random <- VarCorr(object)
    nrow <-  sapply(random, 
                    function(x) if (is.matrix(x)) nrow(x) else length(x))
    maxcol <-max(sapply(random,
                        function(x) if (is.matrix(x)) 1+ncol(x) else 2))
    temp1 <- matrix(NA, nrow=sum(nrow), ncol=maxcol)
    indx <- 0
    for (term in  random) {
        if (is.matrix(term)) {
            k <- nrow(term)
            nc <- ncol(term)  #assume nc > nr (only cases I know so far)
            for (j in 1:k) {
                temp1[j+indx, 1] <- sqrt(term[j,j])
                temp1[j+indx, 2] <- term[j,j]
                if (nc>j) {
                    indx2 <- (j+1):nc
                    temp1[j+indx, 1+ indx2] <- term[j, indx2]
                    }
                }
            }
        else {
            k <- length(term)
            temp1[1:k + indx,1] <- sqrt(term)
            temp1[1:k + indx,2] <- term
            }
        indx <- indx + k
        }
        
    indx <- cumsum(c(1, nrow))   # starting row of each effect
    temp3 <- rep("", nrow(temp1))
    temp3[indx[-length(indx)]] <- names(random)
    xname <- unlist(lapply(random, 
                  function(x) if (is.matrix(x)) dimnames(x)[[1]] else names(x)))
    
    rmat = data.frame(group= temp3, variable= unname(xname), 
                            sd= temp1[,1], variance= temp1[,2])
    if (maxcol == 3) rmat$corr = temp1[,3]

    output <- list(n=object$n, loglik= object$loglik, chi= chi, 
                   coefficients= coef,  random= rmat,
                   na.action= object$na.action, call= object$call)
    class(output) <- "summary.coxme"
    output
    }

print.summary.coxme <- function(x,  digits=max(1L, getOption("digits") - 3L),
                                signif.stars= FALSE, ...) {
    if (!is.null(x$call)) {
        if (TRUE) { # look like lmer
	    cat("Mixed effects coxme model\n Formula:", deparse(x$call$formula),
                "\n")
	    if (!is.null(x$call$data)) cat("    Data:", 
                                           deparse(x$call$data), "\n")
	    cat("\n")
        } else {
            cat("Call:\n")
            dput(x$call)
            cat("\n")
        }
    }
    
    cat("  events, n = ", x$n[1], ', ', x$n[2], sep='')
    if(length(x$na.action) > 0)
        cat(" (", naprint(x$na.action), ")", sep = "")
    cat("\n")

    cat("\nRandom effects:\n")
    print(x$random)

    print(x$chi, quote=FALSE, digits=digits)
    cat("\nFixed effects:\n")

    printCoefmat(x$coefficients, digits=digits, P.values=TRUE, 
                             has.Pvalue=TRUE,
                             signif.stars = signif.stars, ...)
    invisible(x)
}
