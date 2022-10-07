# Author: Cong Jiang (c55jiang@uwaterloo.ca)
# Article:

# data setup
n <- 5000
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
# Assignment mechanism
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
genAA1 <- function(p){
  A <-sample(c("11", "10", "01", "00"), 1, prob= p, replace=TRUE)
  if (A == '11') {return(c(1,1))}
  else if (A == '10') {return(c(1,0))}
  else if (A == '01') {return(c(0,1))}
  else {return(c(0,0))}
}

AA <- data.frame(t(apply(PP, 1, genAA1)))
colnames(AA) <- c('As', 'Ar')
hh <- rep(1:n)
## generate dataset
data.frame(hh, x1, x2, x3, x4, x1r, x2r, x3r, x4r, AA)
dats <- data.frame(hh, x1, x2, x3, x4, A = AA$As)
datr <- data.frame(hh, x1 = x1r, x2 = x2r, x3 = x3r, x4 = x4r, A = AA$Ar)
dat <- rbind(dats, datr)
dat <- dat[order(dat$hh),]

dat

###########
# estimation
###########
# marginal propensity estimation:
margA <- glm(A ~ x1 + x2 + x3 + x4, data = dat, family = "binomial")
summary(margA)
pi_hat <- predict(margA, type = "response")
library(mets)
# Now estimating the OR parameter. OR is constant
# We see a strong dependence with an OR at around 2 that is clearly significant.
estOR <- binomial.twostage(margA, data = dat, var.link = 1,
                           clusters = dat$hh, detail = 0 )
summary(estOR)
# Now estimating the OR parameter. OR depends on covariate x3
########
zz <- aggregate(dat$x3, list(dat$hh), sum)
dat$z <- rep(zz$x, each=2)
theta.des <- model.matrix( ~ z,data=dat)
estOR2 <- binomial.twostage(margA, data = dat, var.link = 1,
                         clusters = dat$hh,theta.des=theta.des)

summary(estOR2)
estOR <- exp(estOR2$theta[1] + estOR2$theta[2] * zz$x)
indS <- seq(1, nrow(dat), by = 2)
indR <- seq(2, nrow(dat), by = 2)
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
tranAA <- function(AS, AR){
  tran <- rep(3, length(AS))
  tran[AS == 1] <- 2
  tran[(AS + AR) == 2] <- 1
  tran[(AS + AR) == 0] <- 4
  return(tran)
}
tranA <- tranAA(dat$A[indS], dat$A[indR])

w <- rep(NA, nrow(w.matrix))
for (i in 1:nrow(w.matrix)) {
  w[i] <- w.matrix[i, tranA[i]]
}
# IPW-based weights
w_ipw <- w/apply(w.matrix, 1, sum)

# Overlap-type weights
w.matrix1 <- cbind(Pir, P10, P01, P00)
w_op <- w * apply(w.matrix1, 1, prod)

length(w_ipw)

nrow(dat)









dat$AAcat <- rep(tranA, each=2) 

dat$ipw <- rep(w_ipw, each = 2)
dat$ow <- rep(w_op, each = 2)
dat.wide <- cbind(dat[indS,],dat[indR,])
subdat11 <- dat.wide[dat.wide$AAcat == 1, ]
subdat10 <- dat.wide[dat.wide$AAcat == 2, ]
subdat01 <- dat.wide[dat.wide$AAcat == 3, ]
subdat00 <- dat.wide[dat.wide$AAcat == 4, ]



cov11 <- apply(subdat11 * subdat11$ipw, 2, sum)/sum(subdat11$ipw)
cov11 <- c(cov11[c(2:5)], cov11[c(12:15)])

cov10 <- apply(subdat10 * subdat10$ipw, 2, sum)/sum(subdat10$ipw)
cov10 <- c(cov10[c(2:5)], cov10[c(12:15)])

cov01 <- apply(subdat01 * subdat01$ipw, 2, sum)/sum(subdat01$ipw)
cov01 <- c(cov01[c(2:5)], cov01[c(12:15)])

cov00 <- apply(subdat00 * subdat00$ipw, 2, sum)/sum(subdat00$ipw)
cov00 <- c(cov00[c(2:5)], cov00[c(12:15)])






matrix(1:12,3,4) * c(2,1,1/3)

c <- c(1,2,3,4)
c[c == 2]


######mets
library(mets)
data("twinstut")
twinstut$binstut <- 1*(twinstut$stutter=="yes")
twinsall <- twinstut
twinstut <- subset(twinstut,zyg%in%c("mz","dz"))
head(twinstut)



margbin <- glm(binstut~factor(sex)+age,data=twinstut,family=binomial())
summary(margbin)
# Now estimating the OR parameter. OR is constant
# We see a strong dependence with an OR at around 8 that is clearly significant.
bina <- binomial.twostage(margbin,data=twinstut,var.link=1,
                          clusters=twinstut$tvparnr,detail=0)
summary(bina)

# Now estimating the OR parameter. OR depends on covariate
theta.des <- model.matrix( ~-1+factor(zyg),data=twinstut)
bin <- binomial.twostage(margbin,data=twinstut,var.link=1,
                         clusters=twinstut$tvparnr,theta.des=theta.des)
summary(bin)


## regular GLM
resp_glm <- glm(A ~ x1 + x2 + x3 + x4, data = dat, family = "binomial")
summary(resp_glm)

## GEE with ind
resp_gee.in <- gee(A ~ x1 + x2 + x3 + x4, data = dat, family = "binomial",
                   id = hh, corstr = "independence",
                   scale.fix = TRUE, scale.value = 1)
summary(resp_gee.in)

## GEE with exchangeable matrix
resp_gee.ex <- gee(A ~ x1 + x2 + x3 + x4, data = dat, family = "binomial",
                   id = hh, corstr = "exchangeable",
                   scale.fix = TRUE, scale.value = 1)
summary(resp_gee.ex)

library(gee)
library(MESS)
library(plm)
data("Grunfeld")
str(Grunfeld)

Grunfeld$firm

qnorm(0.97, mean = -1/0.15, sd = sqrt(1/(3 * 0.15^2)))
2 - 2 * pnorm((log(0.6) + log(0.2) + log(0.9))/3, mean = -1/0.15, sd = sqrt(1/(3 * 0.15^2)))



