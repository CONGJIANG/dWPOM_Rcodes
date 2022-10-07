# Author: Cong Jiang (c55jiang@uwaterloo.ca)
# Article:

#######################################################################################
#MODULE: tranAA
#Estimation joint propensity of correlated treatments in a Household.
#
#INPUT
# AS, AR, the vectors of treatments of (S, R) in a household
#
#OUTPUT
# treatment category, ie, (1,1) --> 1; (1,0) --> 2; (0,1) --> 3; (0,0) --> 4; 
#######################################################################################


tranAA <- function(AS, AR){
  tran <- rep(3, length(AS))
  tran[AS == 1] <- 2
  tran[(AS + AR) == 2] <- 1
  tran[(AS + AR) == 0] <- 4
  return(tran)
}


#######################################################################################
#
#MODULE: BEGINEND
#Estimation joint propensity of correlated treatments in a Household.
#
#INPUT
#n: Vector of cluster sample sizes
#
#OUTPUT
#first: Vector with starting row for cluster i
# last: Vector with ending row for cluster i
# require library(mets)
#######################################################################################

# marginal propensity estimation:
margA <- glm(A ~ x1 + x2 + x3 + x4, data = dat, family = "binomial")
summary(margA)
pi_hat <- predict(margA, type = "response")
# Now estimating the OR parameter. OR depends on covariate x3
########
zz <- aggregate(dat$x3, list(dat$hh), sum)
dat$z <- rep(zz$x, each=2)
theta.des <- model.matrix( ~ z,data=dat)
estOR2 <- binomial.twostage(margA, data = dat, var.link = 1,
                            clusters = dat$hh,theta.des=theta.des)

summary(estOR2)
estOR <- exp(estOR2$theta[1] + estOR2$theta[2] * zz$x)
indS <- seq(1, length(pi_hat), by = 2)
indR <- seq(2, length(pi_hat), by = 2)
MarS <- pi_hat[indS]
MarR <- pi_hat[indR]

# joint propensity (based on the formula)
fir <- 1 - (1-estOR)*(MarS + MarR)
Pir <- (fir - sqrt(fir*fir - 4*estOR*(estOR - 1)*MarS*MarR))/(2*estOR - 2)
Pir <- ifelse(estOR == 1, MarS*MarR, Pir)
Pir
P10 <- MarS - Pir
P01 <- MarR - Pir
P00 <- 1 - MarS - MarR + Pir
w.matrix <- cbind(1/Pir, 1/P10, 1/P01, 1/P00)

tranA <- tranAA(dat$A[indS], dat$A[indR])

w <- rep(NA, nrow(w.matrix))
for (i in 1:nrow(w.matrix)) {
  w[i] <- w.matrix[i, tranA[i]]
}


# IPW-based weights
w/apply(w.matrix, 1, sum)






c(1,2) %*% c(1,2)

decision <- function(xi, psi, phi, xx, xs, xh){
  A <- xi %*% xx + psi %*% xs + phi %*% xh
  B <- psi %*% xs + phi %*% xh
  C <- xi %*% xx + phi %*% xh
  if(A > 0 && B > 0 && C >0) return(c(1, 1))
  if(B < 0 && (xi %*% xx > psi %*% xs) && (A - B) > 0) return(c(1, 0))
  if(C < 0 && (xi %*% xx > psi %*% xs) && (A - C) > 0) return(c(0, 1))
  if(A < 0 && (A - B) < 0 && (A - C) < 0) return(c(0, 0))
}


regret <- function(aopt, a, xi, psi, phi, xx, xs, xh){
  (aopt[1] - a[1])*xi %*% xx + (aopt[2] - a[2])* psi %*% xs + (aopt[1]*aopt[2] - a[1]*a[2])* phi %*% xh
}


