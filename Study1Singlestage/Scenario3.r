# Author: Cong Jiang (c55jiang@uwaterloo.ca)
# Article: Estimating dynamic treatment regimes with household interference for ordinal outcomes 
#           using longitudinal survey data


tranAA <- function(AS, AR){
  tran <- rep(3, length(AS))
  tran[AS == 1] <- 2
  tran[(AS + AR) == 2] <- 1
  tran[(AS + AR) == 0] <- 4
  return(tran)
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

kappa <- function(datwide, coeff, zeta, w_c){
  datwide$Asx1 <- datwide$As*datwide$x1
  datwide$Arx1 <- datwide$Ar*datwide$x1r
  datwide$AsAr <- datwide$As*datwide$Ar
  datwide$AsArx3 <- datwide$As*datwide$Ar*(datwide$x3 + datwide$x3r)
  
  pre_order <- c("x1", "x1r", "x3", "x3r", "As", "Asx1", "Ar", "Arx1","AsAr", "AsArx3" )
  Newdatwide <- datwide[, pre_order]
  eta1 <- zeta[1] - as.matrix(Newdatwide) %*% coeff 
  eta2 <- zeta[2] - as.matrix(Newdatwide) %*% coeff 
  kappa <- expit(eta2)*(1- expit(eta1))*( expit(eta1) - expit(eta2) + 1)
  
  
  Newdatwide$As <- rep(1, nrow(Newdatwide))
  Newdatwide$Ar <- rep(1, nrow(Newdatwide))
  eta1 <- zeta[1] - as.matrix(Newdatwide) %*% coeff 
  eta2 <- zeta[2] - as.matrix(Newdatwide) %*% coeff 
  kappa11 <- expit(eta2)*(1- expit(eta1))*( expit(eta1) - expit(eta2) + 1)
  
  Newdatwide$Ar <- rep(0, nrow(Newdatwide))
  eta1 <- zeta[1] - as.matrix(Newdatwide) %*% coeff 
  eta2 <- zeta[2] - as.matrix(Newdatwide) %*% coeff 
  kappa10 <- expit(eta2)*(1- expit(eta1))*( expit(eta1) - expit(eta2) + 1)
  
  Newdatwide$As <- rep(0, nrow(Newdatwide))
  eta1 <- zeta[1] - as.matrix(Newdatwide) %*% coeff 
  eta2 <- zeta[2] - as.matrix(Newdatwide) %*% coeff 
  kappa00 <- expit(eta2)*(1- expit(eta1))*( expit(eta1) - expit(eta2) + 1)
  
  
  Newdatwide$As <- rep(0, nrow(Newdatwide))
  Newdatwide$Ar <- rep(1, nrow(Newdatwide))
  eta1 <- zeta[1] - as.matrix(Newdatwide) %*% coeff 
  eta2 <- zeta[2] - as.matrix(Newdatwide) %*% coeff 
  kappa01 <- expit(eta2)*(1- expit(eta1))*( expit(eta1) - expit(eta2) + 1)
  w_d <- (w_c* kappa11 * kappa10 * kappa01 * kappa00)/kappa
  return(w_d/mean(w_d))
}


# data setup
n <- 5000
# expit function
expit <- function(x) {1/(1+exp(-x))}

#######################
library(simstudy)
library(ordinal)
library(mets)
#

M = 500
est0 <- matrix(NA, M, 6)
est1 <- matrix(NA, M, 6)
est2 <- matrix(NA, M, 6)
est3 <- matrix(NA, M, 6)
est4 <- matrix(NA, M, 6)
optrat0 <- rep(NA, M)
optrat1 <- rep(NA, M)
optrat2 <- rep(NA, M)
optrat3 <- rep(NA, M)
optrat4 <- rep(NA, M)
optrat0s <- rep(NA, M)
optrat1s <- rep(NA, M)
optrat2s <- rep(NA, M)
optrat3s <- rep(NA, M)
optrat4s <- rep(NA, M)
for (m in 1:M) {
  x1 <- runif(n,0, 1)
  x2 <- rnorm(n,0,1)
  x3 <- rbinom(n,1, 0.5)
  x4 <- rbinom(n,1, 0.75)
  x1r <- runif(n,0, 1)
  x2r <- rnorm(n,0,1)
  x3r <- rbinom(n,1, 0.5)
  x4r <- rbinom(n,1, 0.75)
  
  # marginal propensity
  Pmar <-  expit(-1.15 + 0.5*exp(x1) - 0.25*x2^2 + 0.25*x3 + 0.6*x4)
  Pmarr <-  expit(-1.15 + 0.5*exp(x1r) - 0.25*x2r^2 + 0.25*x3r + 0.6*x4r)
  # odds ratio
  tao <- exp( -0.25 + 0.5*(x3+ x3r) + 0.25*(x1 + x1r))
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
  datwide <- data.frame(hh, x1, x2, x3, x4, x1r, x2r, x3r, x4r, AA)
  dats <- data.frame(hh, x1, x2, x3, x4, A = AA$As)
  datr <- data.frame(hh, x1 = x1r, x2 = x2r, x3 = x3r, x4 = x4r, A = AA$Ar)
  dat <- rbind(dats, datr)
  dat <- dat[order(dat$hh),]
  
  trtfree <-  cos(datwide$x1 + datwide$x1r) - (datwide$x1 + datwide$x1r)^3 - log(abs(1/datwide$x1)) - 2*datwide$x1r^2 + (datwide$x3 + datwide$x3r)^3
  gammaS <- datwide$As *(-0.5 + 1 * datwide$x1)
  gammaR <- datwide$Ar *(-0.5 + 1  * datwide$x1r)
  gammaInt <- datwide$As *datwide$Ar *(-1 + 0.5 * (datwide$x3 + datwide$x3r))
  y <- trtfree + gammaS + gammaR + gammaInt
  datwide$mu <- y
  
  
  probs <- c(0.65, 0.25, 0.1)
  cprop <- cumsum(probs)
  
  # map cumulative probs to thresholds for reference group
  
  (gamma.c <- qlogis(cprop))
  matlp <- matrix(rep(gamma.c, nrow(datwide)), 
                  ncol = length(cprop), 
                  byrow = TRUE)
  
  
  # set individual thresholds based on covariates,
  # which is an additive shift from the reference group
  # based on mu
  
  matlpInd <- matlp - datwide$mu
  
  
  
  # convert log odds to cumulative probabability
  matcump <- 1 / (1 + exp(-matlpInd))
  matcump <- cbind(0, matcump)
  
  
  # convert cumulative probs to category probs:
  # originally, I used a loop to do this, but
  # thought it would be better to vectorize.
  # see 2nd addendum for time comparison - not
  # much difference
  
  p <- t(t(matcump)[-1,] - t(matcump)[-4,])
  
  
  # generate individual level category outcomes based on p
  cat <- simstudy:::matMultinom(p)
  catF <- ordered(cat)
  datwide$cat = catF
  
  
  
  fmm0 <- clm(cat ~  x1 + x1r + x3 + x3r + As + I(As *x1) + Ar + I(Ar *x1r) + I(As *Ar) + I(As *Ar *(x3 + x3r)), data=datwide)
  cv0 <- fmm0$coefficients[7:length(fmm0$coefficients)]
  
  # Pmar <-  expit(-1.5 + exp(x1) - 0.25*x2^2 + 0.25*x3 + 0.6*x4)
  alpha <- glm(A ~ exp(x1) + x2^2 + x3 + x4, family=binomial(link='logit'), dat)
  pi_hat <- predict(alpha, type = "response")
  # weights
  w <- abs(dat$A - pi_hat)
  indS <- seq(1, nrow(dat), by = 2)
  indR <- seq(2, nrow(dat), by = 2)
  w_a <- w[indS]*w[indR]
  fmm1 <- clm(cat ~  x1 + x1r + x3 + x3r  + As + I(As *x1) + Ar + I(Ar *x1r) + I(As *Ar) + I(As *Ar *(x3 + x3r)),
              w = w_a, data=datwide)
  cv1 <- fmm1$coefficients[7:length(fmm1$coefficients)]
  
  # Now estimating the OR parameter. OR depends on covariate x3
  ########
  zz <- aggregate(dat$x3, list(dat$hh), sum)
  dat$z <- rep(zz$x, each=2)
  
  zz1 <- aggregate(dat$x1, list(dat$hh), sum)
  dat$z1 <- rep(zz1$x, each=2)
  
  theta.des <- model.matrix( ~ z + z1, data=dat)
  estOR2 <- binomial.twostage(alpha, data = dat, var.link = 1,
                              clusters = dat$hh,theta.des=theta.des)
  
  #summary(estOR2)
  estOR <- exp(estOR2$theta[1] + estOR2$theta[2] * zz$x + estOR2$theta[3] * zz1$x)
  MarS <- pi_hat[indS]
  MarR <- pi_hat[indR]
  
  # joint propensity (based on the formula)
  fir <- 1 - (1-estOR)*(MarS + MarR)
  Pir <- (fir - sqrt(fir*fir - 4*estOR*(estOR - 1)*MarS*MarR))/(2*estOR - 2)
  Pir <- ifelse(estOR == 1, MarS*MarR, Pir)
  P10 <- MarS - Pir
  P01 <- MarR - Pir
  P00 <- 1 - MarS - MarR + Pir
  w.matrix <- cbind(1/Pir, 1/P10, 1/P01, 1/P00)
  
  tranA <- tranAA(dat$A[indS], dat$A[indR])
  # hist(tranA)
  w <- rep(NA, nrow(w.matrix))
  for (i in 1:nrow(w.matrix)) {
    w[i] <- w.matrix[i, tranA[i]]
  }
  
  w_b <-w/apply(w.matrix, 1, sum)
  fmm2 <- clm(cat ~  x1 + x1r + x3 + x3r  + As + I(As *x1) + Ar + I(Ar *x1r) + I(As *Ar) + I(As *Ar *(x3 + x3r)),
              w = w_b, data=datwide)
  cv2 <- fmm2$coefficients[7:length(fmm2$coefficients)]
  
  
  # Overlap-type weights
  w.matrix1 <- cbind(Pir, P10, P01, P00)
  w_c <- w * apply(w.matrix1, 1, prod)
  w_c <- w_c/mean(w_c)
  fmm3 <- clm(cat ~  x1 + x1r + x3 + x3r  + As + I(As *x1) + Ar + I(Ar *x1r) + I(As *Ar) + I(As *Ar *(x3 + x3r)),
              w = w_c, data=datwide)
  cv3 <- fmm3$coefficients[7:length(fmm3$coefficients)]
  
  coeff <- as.vector(fmm3$coefficients[3:length(fmm3$coefficients)])
  zeta <- as.vector(fmm3$coefficients[1:2])
  
  w_d <- kappa(datwide, coeff, zeta, w_c)
  fmm4 <- clm(cat ~  x1 + x1r + x3 + x3r + As + I(As *x1) + Ar + I(Ar *x1r) + I(As *Ar) + I(As *Ar *(x3 + x3r)),
              w = w_d, data=datwide)
  cv4 <- fmm4$coefficients[7:length(fmm4$coefficients)]
  
  
  est0[m,] <- cv0
  est1[m,] <- cv1
  est2[m,] <- cv2
  est3[m,] <- cv3
  est4[m,] <- cv4
  
  
  dval <- decision.value(cv0[1:2],cv0[3:4], cv0[5:6], cbind(1, datwide$x1), cbind(1,datwide$x1r),
                         cbind(1, datwide$x3+datwide$x3r))
  aopt <- t(apply(dval, 1, decision))
  
  dval.t <- decision.value(c(-0.5, 1),c(-0.5,1), c(-1, 0.5), cbind(1, datwide$x1), cbind(1,datwide$x1r),
                           cbind(1, datwide$x3+datwide$x3r))
  aopt.t <- t(apply(dval.t, 1, decision))
  
  cont <- apply( aopt == aopt.t, 1, sum)
  optrat0[m] <- table(cont)[length(table(cont))]/n
  optrat0s[m] <- 1 - table(cont)[1]/n
  # 1
  dval <- decision.value(cv1[1:2],cv1[3:4], cv1[5:6], cbind(1, datwide$x1), cbind(1,datwide$x1r),
                         cbind(1, datwide$x3+datwide$x3r))
  aopt <- t(apply(dval, 1, decision))
  cont <- apply( aopt == aopt.t, 1, sum)
  optrat1[m] <- table(cont)[length(table(cont))]/n
  optrat1s[m] <- 1 - table(cont)[1]/n
  # 2
  dval <- decision.value(cv2[1:2],cv2[3:4], cv2[5:6], cbind(1, datwide$x1), cbind(1,datwide$x1r),
                         cbind(1, datwide$x3+datwide$x3r))
  aopt <- t(apply(dval, 1, decision))
  cont <- apply( aopt == aopt.t, 1, sum)
  optrat2[m] <- table(cont)[length(table(cont))]/n
  optrat2s[m] <- 1 - table(cont)[1]/n
  # 3
  dval <- decision.value(cv3[1:2],cv3[3:4], cv3[5:6], cbind(1, datwide$x1), cbind(1,datwide$x1r),
                         cbind(1, datwide$x3+datwide$x3r))
  aopt <- t(apply(dval, 1, decision))
  cont <- apply( aopt == aopt.t, 1, sum)
  optrat3[m] <- table(cont)[length(table(cont))]/n
  optrat3s[m] <- 1 - table(cont)[1]/n
  # 4
  dval <- decision.value(cv4[1:2],cv4[3:4], cv4[5:6], cbind(1, datwide$x1), cbind(1,datwide$x1r),
                         cbind(1, datwide$x3+datwide$x3r))
  aopt <- t(apply(dval, 1, decision))
  cont <- apply( aopt == aopt.t, 1, sum)
  optrat4[m] <- table(cont)[length(table(cont))]/n
  optrat4s[m] <- 1 - table(cont)[1]/n
  
}



library(vioplot)
par(mfrow=c(3,2))
vioplot(est0[,1], est1[,1], est2[,1], est3[,1], est4[,1],
        names=c( "M0","M1", "M2", "M3", "M4"),xlab= "H = 3000", ylab = expression(paste( xi[0], " estimates")))
abline(h = -0.5, col = "red")

vioplot(est0[,2], est1[,2], est2[,2], est3[,2],est4[,2],
        names=c( "M0","M1", "M2", "M3", "M4"),xlab= "H = 3000", ylab = expression(paste( xi[1], " estimates")))
abline(h = 1, col = "red")
vioplot(est0[,3], est1[,3], est2[,3], est3[,3],est4[,3],
        names=c( "M0","M1", "M2", "M3", "M4"),xlab= "H = 3000",  ylab = expression(paste( psi[0], " estimates")))
abline(h = -0.5, col = "red") 

vioplot(est0[,4], est1[,4], est2[,4], est3[,4],est4[,4],
        names=c( "M0","M1", "M2", "M3", "M4"), xlab= "H = 3000", ylab = expression(paste( psi[1], " estimates")))
abline(h = 1, col = "red") 

vioplot(est0[,5], est1[,5], est2[,5], est3[,5],est4[,5],
        names=c( "M0","M1", "M2", "M3", "M4"),xlab= "H = 3000",  ylab = expression(paste( phi[0], " estimates")))
abline(h = -1, col = "red") 

vioplot(est0[,6], est1[,6], est2[,6], est3[,6],est4[,6],
        names=c( "M0","M1", "M2", "M3", "M4"), xlab= "H = 3000", ylab = expression(paste( phi[1], " estimates")))
abline(h = 0.5, col = "red") 

mean(optrat0); mean(optrat1); mean(optrat2); mean(optrat3);mean(optrat4)
mean(optrat0s); mean(optrat1s); mean(optrat2s); mean(optrat3s);mean(optrat4s)