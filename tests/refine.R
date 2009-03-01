library(coxme)
options(na.action='na.exclude', contrasts=c('contr.treatment', 'contr.poly'))
aeq <- function(x,y, ...) all.equal(as.vector(x), as.vector(y), ...)

#
# Same data set as slope1
#
set.seed(56)
n.subject <- seq(180, by=21, length=9) # number of subjects
slope <- sort(-.5 + rnorm(9, sd=.5))         # true treament effects

inst <- rep(1:9, n.subject)
n <- length(inst)
simdata <- data.frame(id=1:n, inst=inst,
                      trt= rep(0:1, length=n),
                      age= runif(n, 40, 70))
#risk goes up 30%/decade of age
simdata$hazard <- .8* exp(simdata$trt * rep(slope, n.subject) +
                          (simdata$age-55) * .03)

rtime <- function(hazard, censor=c(1,2)) {
    stime <- rexp(length(hazard), rate=hazard)
    ctime <- runif(length(hazard), censor[1], censor[2])
    list(time= pmin(stime, ctime), status=1*(stime <=ctime))
    }
temp <- rtime(simdata$hazard)
simdata$time <- temp$time
simdata$status <- temp$status

#
# Test out the refine.n code, using the simdata
#  A simple diagonal variance
#
nsim <- 100
var  <- .3   #sizeable

set.seed(20)
fit1 <- coxme(Surv(time, status) ~ age + trt + (trt|inst) + strata(inst),
              vfixed=.3, simdata, refine.n=nsim)

debug <- fit1$refine.debug  # Optional return item -see end of coxme.fit

nfrail <- length(unlist(fit1$frail))
set.seed(20)
bmat <- matrix(rnorm(nfrail*nsim), nfrail)  # replicate the simulations
if (!is.null(debug)) all.equal(bmat, debug$bmat)

clog <- double(nsim)
Xmat <- scale(as.matrix(simdata[,c('age', 'trt')]), fit1$means, FALSE)
bmat <- bmat * sqrt(.3)  #the actual bmat used or fits, var=.3

# Part 1, loglik for a set of nearby Cox models
fix.lp <- Xmat %*% fixef(fit1)
for (i in 1:nsim) {
    lp <- fix.lp + bmat[simdata$inst,i]*simdata$trt
    tfit <- coxph(Surv(time, status) ~ offset(lp) + strata(inst), simdata)
    clog[i] <- tfit$loglik
    }
if (!is.null(debug)) aeq(clog, debug$loglik)

# Part 2: Taylor series for the PPL
bhat <- unlist(fit1$frail)  
b.sig <- t(bmat-bhat) %*% fit1$hmat[1:nfrail, 1:nfrail]  #b time sqrt(H)
taylor <- rowSums(b.sig^2)/2

temp2 <- cbind(clog, fit1$log[3] + colSums(bmat^2)/.6 - taylor)
m2 <- mean(temp2)
errhat <- exp(temp2[,1]-m2) - exp(temp2[,2]-m2)
escale <- exp(m2-fit1$log[2])

if (!is.null(debug)) {
    aeq(errhat, debug$errhat)
    aeq(escale, debug$escale)
    }

#Interestingly, even though both the above pass with standard tolerance,
#  the one below will fail.  The escale variable is slightly different,
#  since one is based on fit1$log[2] before correction (internal to coxme.fit)
#  and the one here on the value after correction.
aeq(escale * c(mean(errhat), sqrt(var(errhat)/nsim)), fit1$refine,
    tolerance=1e-5)


