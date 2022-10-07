genAA1 <- function(p){
  A <-sample(c("11", "10", "01", "00"), 1, prob= p, replace=TRUE)
  if (A == '11') {return(c(1,1))}
  else if (A == '10') {return(c(1,0))}
  else if (A == '01') {return(c(0,1))}
  else {return(c(0,0))}
}
decision.value <- function(xi, psi, phi, xx, xs, xh){
  A <- xx %*% xi + xs  %*%  psi + xh  %*%  phi
  B <- xx %*% xi 
  C <- xs  %*%  psi
  return(data.frame(A = A, B = B, C = C))
}
decision <- function(value){
  value <- as.vector(value)
  A <- value[1]; B <- value[2]; C <- value[3]; 
  if((A > 0) && (A > B) && (A > C)) return(c(1, 1))
  if((B > A) && (B > C) && (B > 0)) return(c(1, 0))
  if((C > A) && (C > B) && (C > 0)) return(c(0, 1))
  else  return(c(0, 0))
}

tranAA <- function(AS, AR){
  tran <- rep(3, length(AS))
  tran[AS == 1] <- 2
  tran[(AS + AR) == 2] <- 1
  tran[(AS + AR) == 0] <- 4
  return(tran)
}

genTrt <- function(x1, x2, x1r, x2r) {
  # marginal propensity
  Pmar <-  expit(-1 + 1.15*exp(x1) - 0.5*x2)
  Pmarr <-  expit(-1 + 1.15*exp(x1r) - 0.5*x2r)
  # odds ratio
  tao <- exp( -0.15 + 0.25*(x1 + x1r))
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
  AA <- data.frame(t(apply(PP, 1, genAA1)))
  colnames(AA) <- c('As', 'Ar')
  return(AA)
}

n <- 1000
x11 <- rnorm(n,0,1)
x12 <- rbinom(n,1, 0.5)
x11r <- rnorm(n,0,1)
x12r <- rbinom(n,1, 0.5)
AA1 <- genTrt(x11, x12, x11r, x12r)
colnames(AA1) <- c('As1', 'Ar1')
hist(AA1$As1)
hist(AA1$Ar1)

hh <- rep(1:n)
## generate dataset
dat1s <- data.frame(hh, x11, x12, A1 = AA1$As1)
dat1r <- data.frame(hh, x11 = x11r, x12 = x12r, A1 = AA1$Ar1)
dat1 <- rbind(dat1s, dat1r)
datwide1 <- data.frame(hh, x11, x12, x11r, x12r, AA1)
datHH1 <- data.frame(hh, x11, x11r, x12H =  x12 + x12r, AA1)



x21 <- (x11 + rnorm(n,0,1))/2
x22 <- rbinom(n,1, (0.1 + x12/2) )
x21r <- (x11r + rnorm(n,0,1))/2
x22r <- rbinom(n,1, (0.1 + x12r/2) )
AA2 <- genTrt(x21, x22, x21r, x22r)
colnames(AA2) <- c('As2', 'Ar2')
hist(AA2$As2)
hist(AA2$Ar2)

dat2s <- data.frame(hh, x21, x22, A2 = AA2$As2)
dat2r <- data.frame(hh, x21 = x21r, x22 = x22r, A2 = AA2$Ar2)
dat2 <- rbind(dat2s, dat2r)
datwide2 <- data.frame(hh, x21, x22, x21r, x22r, AA2)
datHH2 <- data.frame(hh, x21, x21r, x22H =  x22 + x22r, AA2)


AA1opt.val <- decision.value(c(-0.5, 1),c(-0.5,1), c(1, -0.5), cbind(1, x11), cbind(1,x11r),
                             cbind(1,x12+x12r))
AA1opt <- t(apply(AA1opt.val, 1, decision))
AA2opt.val <- decision.value(c(-0.25, 0.5),c(-0.25,0.5), c(1, -0.5), cbind(1, x21), cbind(1,x21r),
                             cbind(1,x22+x22r))
AA2opt <- t(apply(AA2opt.val, 1, decision))
# regrets



trtfree <-  cos(pi * (datwide1$x11 + datwide1$x11r)) + cos(pi * (datwide2$x21 + datwide2$x21r))
# - (datwide$x1 + datwide$x1r)^3 - log(abs(1/datwide$x1)) - 2*datwide$x1r^2 + (datwide$x3 + datwide$x3r)^3
gammaS1 <- datwide1$As1 *(-0.5 + 1 * datwide1$x11)
gammaR1 <- datwide1$Ar1 *(-0.5 + 1  * datwide1$x11r)
gammaInt1 <- datwide1$As1 * datwide1$Ar1 *(-1 + 0.5 * (datwide1$x12 + datwide1$x12r))

gammaS2 <- datwide2$As2 *(-0.5 + 1 * datwide2$x21)
gammaR2 <- datwide2$Ar2 *(-0.5 + 1  * datwide2$x21r)
gammaInt2 <- datwide2$As2 * datwide2$Ar2 *(-1 + 0.5 * (datwide2$x22 + datwide2$x22r))
y <- trtfree + gammaS1 + gammaR1 + gammaInt1 + gammaS2 + gammaR2 + gammaInt2 
mu <- y



probs <- c(0.65, 0.25, 0.1)

genOrdOut <- function(probs, mu){
  cprop <- cumsum(probs)
  # map cumulative probs to thresholds for reference group
  gamma.c <- qlogis(cprop)
  matlp <- matrix(rep(gamma.c, length(mu)), 
                  ncol = length(cprop), 
                  byrow = TRUE)
  matlpInd <- matlp - mu
  matcump <- 1 / (1 + exp(-matlpInd))
  matcump <- cbind(0, matcump)
  p <- t(t(matcump)[-1,] - t(matcump)[-4,])
  cat <- simstudy:::matMultinom(p)
  catF <- ordered(cat)
  return(catF)
}

datHH2$OrdY <- genOrdOut(probs, mu)
datHH2
###################################

# analysis
# models to be passed to dWGLM
outcome.mod <- OrdY ~ 1
blip.modxi <- list(~ x11, ~x21)
blip.modpsi <- list(~ x11r, ~x21r)
blip.modphi <- list(~ x12H, ~x22H)
treat.mod <- list(A1 ~ x11, A2 ~ x21)
tf.mod <-  list(~ x11 + x11r + x12H, ~ x22)
dataHH <- data.frame(datHH1, datHH2)
dataInd <- data.frame(dat1, dat2)
dataInd <- dataInd[order(dataInd$hh), ]
k <- 2
j = 1
aaHH <- list(As1 ~ Ar1, As2 ~ Ar2)

obj <- list()
# extract original outcome
Y <- model.response(model.frame(outcome.mod,dataHH))

# regress Y on terms in treatment free model and blip model
Hxi <- model.matrix(blip.modxi[[j]], model.frame(blip.modxi[[j]], dataHH))
Hpsi <- model.matrix(blip.modpsi[[j]], model.frame(blip.modpsi[[j]], dataHH))
Hphi <- model.matrix(blip.modphi[[j]], model.frame(blip.modphi[[j]], dataHH))
Hbeta <- model.matrix(tf.mod[[j]], model.frame(tf.mod[[j]], dataHH))
R = 5
Theta <- matrix(NA, R, dim(Hxi)[2] + dim(Hpsi)[2] +dim(Hphi)[2])
As <- model.frame(aaHH[[j]], dataHH)[,1]
Ar <- model.frame(aaHH[[j]], dataHH)[,2]

##########
dataInd
indS <- seq(1, nrow(dataInd), by = 2)
indR <- seq(2, nrow(dataInd), by = 2)
treat.mod <- list(A1 ~ x11 + x12, A2~ x21 + x22)
alpha <- glm(treat.mod[[j]],binomial, dataInd)
pi_hat <- predict(alpha, type = "response")
MarS <- pi_hat[indS]
MarR <- pi_hat[indR]
A <- model.response(model.frame(treat.mod[[j]],dataInd))
# weights
w <- abs(A - fitted(alpha))



zz1 <- aggregate(dataInd$x11, list(dataInd$hh), sum)
zz2 <- aggregate(dataInd$x21, list(dataInd$hh), sum)
dataInd$z1 <- rep(zz1$x, each=2)
dataInd$z2 <- rep(zz2$x, each=2)

theta.des <- model.matrix( ~ z1, data=dataInd)
estOR2 <- binomial.twostage(alpha, data = dataInd, var.link = 1,
                            clusters = dataInd$hh,theta.des=theta.des)
estOR <- exp(estOR2$theta[1] + estOR2$theta[2] * zz1$x)
tranA <- tranAA(dataInd$A1[indS], dataInd$A1[indR])

weights <- function(MarS, MarR, estOR, tranA){
  # joint propensity (based on the formula)
  fir <- 1 - (1-estOR)*(MarS + MarR)
  Pir <- (fir - sqrt(fir*fir - 4*estOR*(estOR - 1)*MarS*MarR))/(2*estOR - 2)
  Pir <- ifelse(estOR == 1, MarS*MarR, Pir)
  P10 <- MarS - Pir
  P01 <- MarR - Pir
  P00 <- 1 - MarS - MarR + Pir
  w.matrix <- cbind(1/Pir, 1/P10, 1/P01, 1/P00)
  w <- rep(NA, nrow(w.matrix))
  for (i in 1:nrow(w.matrix)) {
    w[i] <- w.matrix[i, tranA[i]]
  }
  
  w_b <-w/apply(w.matrix, 1, sum)
  
  # Overlap-type weights
  w.matrix1 <- cbind(Pir, P10, P01, P00)
  w_c <- w * apply(w.matrix1, 1, prod)
  return(data.frame(ipw = w_b, op = w_c))
}

opw <- weights(MarS, MarR, estOR, tranA)$op
##########
#first stage:

clm1 <- clm(Y ~ 0+ Hbeta + I(As*Hxi) +  I(Ar*Hpsi) + I(As*Ar*Hphi), weights = opw)
coeff <- clm1$coefficients[4:length(clm1$coefficients)]
zeta <- clm1$coefficients[1:2]


kappa <- function(Hbeta, Hxi, Hpsi, Hphi, coeff, zeta, opw, As, Ar){
  Hbeta <- Hbeta[, 2:dim(Hbeta)[2]]
  Newdatwide <- data.frame(Hbeta,As*Hxi, Ar*Hpsi, As*Ar*Hphi)
  
  eta1 <- zeta[1] - as.matrix(Newdatwide) %*% coeff 
  eta2 <- zeta[2] - as.matrix(Newdatwide) %*% coeff 
  kappa <- expit(eta2)*(1- expit(eta1))*( expit(eta1) - expit(eta2) + 1)
  
  
  As <- rep(1, nrow(Newdatwide))
  Ar <- rep(1, nrow(Newdatwide))
  Newdatwide <- data.frame(Hbeta,As*Hxi, Ar*Hpsi, As*Ar*Hphi)
  eta1 <- zeta[1] - as.matrix(Newdatwide) %*% coeff 
  eta2 <- zeta[2] - as.matrix(Newdatwide) %*% coeff 
  kappa11 <- expit(eta2)*(1- expit(eta1))*( expit(eta1) - expit(eta2) + 1)
  
  Ar <- rep(0, nrow(Newdatwide))
  Newdatwide <- data.frame(Hbeta,As*Hxi, Ar*Hpsi, As*Ar*Hphi)
  eta1 <- zeta[1] - as.matrix(Newdatwide) %*% coeff 
  eta2 <- zeta[2] - as.matrix(Newdatwide) %*% coeff 
  kappa10 <- expit(eta2)*(1- expit(eta1))*( expit(eta1) - expit(eta2) + 1)
  
  As <- rep(0, nrow(Newdatwide))
  Newdatwide <- data.frame(Hbeta,As*Hxi, Ar*Hpsi, As*Ar*Hphi)
  eta1 <- zeta[1] - as.matrix(Newdatwide) %*% coeff 
  eta2 <- zeta[2] - as.matrix(Newdatwide) %*% coeff 
  kappa00 <- expit(eta2)*(1- expit(eta1))*( expit(eta1) - expit(eta2) + 1)
  
  
  As <- rep(0, nrow(Newdatwide))
  Ar <- rep(1, nrow(Newdatwide))
  Newdatwide <- data.frame(Hbeta,As*Hxi, Ar*Hpsi, As*Ar*Hphi)
  eta1 <- zeta[1] - as.matrix(Newdatwide) %*% coeff 
  eta2 <- zeta[2] - as.matrix(Newdatwide) %*% coeff 
  kappa01 <- expit(eta2)*(1- expit(eta1))*( expit(eta1) - expit(eta2) + 1)
  w_d <- (opw* kappa11 * kappa10 * kappa01 * kappa00)/kappa
  return(w_d/mean(w_d))
}

w_d <- kappa(Hbeta, Hxi, Hpsi, Hphi, coeff, zeta, opw, As, Ar)
#######################


clm2 <- clm(Y ~ 0+ Hbeta + I(As*Hxi) +  I(Ar*Hpsi) + I(As*Ar*Hphi))
par <- clm2$coefficients[-(seq_len(dim(Hxi)[2] + dim(Hpsi)[2] +dim(Hphi)[2]))]


dval <- decision.value(par[1:2],par[3:4], par[5:6], Hxi , Hpsi, Hphi)
aopt <- t(apply(dval, 1, decision))

predCLM <- function(Hbeta, Hxi, Hpsi, Hphi, coeff, zeta, aopt){
  Hbeta <- Hbeta[, 2:dim(Hbeta)[2]]
  As <- aopt[,1]; Ar <- aopt[,2];
  Newdatwide <- data.frame(Hbeta,As*Hxi, Ar*Hpsi, As*Ar*Hphi)
  
  eta1 <- zeta[1] - as.matrix(Newdatwide) %*% coeff 
  eta2 <- zeta[2] - as.matrix(Newdatwide) %*% coeff 
  p1 <- expit(eta1); p2 <- expit(eta2) - expit(eta1); p3 <- 1 - p1 - p2;
  return(data.frame(p1 = p1, p2 = p2, p3 = p3))
}
predP <- predCLM(Hbeta, Hxi, Hpsi, Hphi, coeff, zeta, aopt)
# generate individual level category outcomes based on p
cat <- simstudy:::matMultinom(as.matrix(predP))
catF <- ordered(cat)

