# Author: Cong Jiang (c55jiang@uwaterloo.ca)
# Article:

# data setup
n <- 500
# expit function
expit <- function(x) {1/(1+exp(-x))}

x1 <- runif(n,0, 1)
x2 <- round(runif(n, min = 0, max = 3), digits = 0)
x3 <- rbinom(n,1, 0.5)
x4 <- rbinom(n,1, 0.75)
x1r <- runif(n,0, 1)
x2r <- round(runif(n, min = 0, max = 3), digits = 0)
x3r <- rbinom(n,1, 0.5)
x4r <- rbinom(n,1, 0.75)
# marginal propensity
Pmar <-  expit(-1.25 + x1 + 0.05*x2 + 0.25*x3 + 0.6*x4)
Pmarr <-  expit(-1.25 + x1r + 0.05*x2r + 0.25*x3r + 0.6*x4r)
# odds ratio
tao <- exp( -0.25 + 0.5*(x3+ x3r) )
# joint propensity (based on the formula)
fir <- 1 - (1-tao)*(Pmar + Pmarr)
Pir <- (fir - sqrt(fir*fir - 4*tao*(tao - 1)*Pmar*Pmarr))/(2*tao - 2)
Pir <- ifelse(tao == 1, Pmar*Pmarr, Pir)
P10 <- Pmar - Pir
P01 <- Pmarr - Pir
P00 <- 1 - Pmar - Pmarr + Pir
PP <- cbind(Pir, P10, P01, P00)
# test apply(PP, 1, sum)
# generate AA based on the joint propensity
genAA <- function(p){
  A <-sample(c("11", "10", "01", "00"), 1, prob= p, replace=TRUE)
  if (A == '11') {return(c(1,1))}
  else if (A == '10') {return(c(1,0))}
  else if (A == '01') {return(c(0,1))}
  else {return(c(0,0))}
}

AA <- data.frame(t(apply(PP, 1, genAA)))
colnames(AA) <- c('Ai', 'Ar')
hh <- rep(1:n)
## generate dataset
data.frame(hh, x1, x2, x3, x4, x1r, x2r, x3r, x4r, AA)
dati <- data.frame(hh, x1, x2, x3, x4, A = AA$Ai)
datr <- data.frame(hh, x1 = x1r, x2 = x2r, x3 = x3r, x4 = x4r, A = AA$Ar)
dat <- rbind(dati, datr)
dat <- dat[order(dat$hh),]

dat

###########
# estimation
###########
# marginal propensity estimation:
margA <- glm(A ~ x1 + x2 + x3 + x4, data = dat, family = "binomial")
summary(margA)

# Now estimating the OR parameter. OR is constant
# We see a strong dependence with an OR at around 2 that is clearly significant.
estOR <- binomial.twostage(margA, data = dat, var.link = 1,
                           clusters = dat$hh, detail = 0 )
summary(estOR)
# Now estimating the OR parameter. OR depends on covariate
########
zz <- aggregate(dat$x3, list(dat$hh), sum)
dat$z <- rep(zz$x, each=2)
theta.des <- model.matrix( ~ z,data=dat)
estOR2 <- binomial.twostage(margA, data = dat, var.link = 1,
                         clusters = dat$hh,theta.des=theta.des)

summary(estOR2)

######mets
library(mets)
