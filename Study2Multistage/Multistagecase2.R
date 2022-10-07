library("vioplot")
library("mets")
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

predCLM <- function(datHH, coeff, zeta, aopt){
  datHH$Asx1 <- datHH$As1*datHH$x11
  datHH$Arx1 <- datHH$Ar1*datHH$x11r
  datHH$AsAr1 <- datHH$As1*datHH$Ar1
  datHH$AsArHH1 <- datHH$As1*datHH$Ar1*datHH$x12H
  
  datHH$As2 <- aopt[,1]
  datHH$Ar2 <- aopt[,2]
  datHH$Asx2 <- datHH$As2*datHH$x21
  datHH$Arx2 <- datHH$Ar2*datHH$x21r
  datHH$AsAr2 <- datHH$As2*datHH$Ar2
  datHH$AsArHH2 <- datHH$As2*datHH$Ar2*datHH$x22H
  pre_order1 <- c("x11",  "x11r" , "x12H",  "As1", "Asx1", "Ar1", "Arx1","AsAr1", "AsArHH1")
  pre_order2 <- c("x21",  "x21r" , "x22H",  "As2", "Asx2", "Ar2", "Arx2","AsAr2", "AsArHH2" )
  Newdatwide <- datHH[, c(pre_order1, pre_order2)]
  eta1 <- zeta[1] - as.matrix(Newdatwide) %*% coeff 
  eta2 <- zeta[2] - as.matrix(Newdatwide) %*% coeff 
  p1 <- expit(eta1)
  p2 <- expit(eta2) - expit(eta1)
  p3 <- 1 - p1 - p2
  return(data.frame(p1 = p1, p2 = p2, p3 = p3))
}

genTrt <- function(x1, x2, x1r, x2r) {
  # marginal propensity
  Pmar <-  expit(-1 + 1.15*exp(x1) - 0.5*x2)
  Pmarr <-  expit(-1 + 1.15*exp(x1r) - 0.5*x2r)
  # odds ratio
  tao <- exp( -0.15 + 0.25*(x2 + x2r))
  # joint propensity (based on the formula)
  fir <- 1 - (1-tao)*(Pmar + Pmarr)
  Pir <- (fir - sqrt(fir*fir - 4*tao*(tao - 1)*Pmar*Pmarr))/(2*tao - 2)
  Pir <- ifelse(tao == 1, Pmar*Pmarr, Pir)
  P10 <- Pmar - Pir
  P01 <- Pmarr - Pir
  P00 <- 1 - Pmar - Pmarr + Pir
  PP <- abs(cbind(Pir, P10, P01, P00))
  # test apply(PP, 1, sum)
  # generate AA based on the joint propensity
  AA <- data.frame(t(apply(PP, 1, genAA1)))
  colnames(AA) <- c('As', 'Ar')
  return(AA)
}

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
  return(data.frame(ipw = w_b/mean(w_b), op = abs(w_c/mean(w_c))))
}




kappa1 <- function(datHH, coeff, zeta, w_c){
  datHH$Asx1 <- datHH$As1*datHH$x11
  datHH$Arx1 <- datHH$Ar1*datHH$x11r
  datHH$AsAr1 <- datHH$As1*datHH$Ar1
  datHH$AsArHH1 <- datHH$As1*datHH$Ar1*datHH$x12H
  
  pre_order1 <- c("x11",  "x11r" , "x12H",  "As1", "Asx1", "Ar1", "Arx1","AsAr1", "AsArHH1")
  Newdatwide <- datHH[, pre_order1]
  eta1 <- zeta[1] - as.matrix(Newdatwide) %*% coeff 
  eta2 <- zeta[2] - as.matrix(Newdatwide) %*% coeff 
  kappa <- expit(eta2)*(1- expit(eta1))*( expit(eta1) - expit(eta2) + 1)
  
  
  Newdatwide$As1 <- rep(1, nrow(Newdatwide))
  Newdatwide$Ar1 <- rep(1, nrow(Newdatwide))
  eta1 <- zeta[1] - as.matrix(Newdatwide) %*% coeff 
  eta2 <- zeta[2] - as.matrix(Newdatwide) %*% coeff 
  kappa11 <- expit(eta2)*(1- expit(eta1))*( expit(eta1) - expit(eta2) + 1)
  
  Newdatwide$Ar1 <- rep(0, nrow(Newdatwide))
  eta1 <- zeta[1] - as.matrix(Newdatwide) %*% coeff 
  eta2 <- zeta[2] - as.matrix(Newdatwide) %*% coeff 
  kappa10 <- expit(eta2)*(1- expit(eta1))*( expit(eta1) - expit(eta2) + 1)
  
  Newdatwide$As1 <- rep(0, nrow(Newdatwide))
  eta1 <- zeta[1] - as.matrix(Newdatwide) %*% coeff 
  eta2 <- zeta[2] - as.matrix(Newdatwide) %*% coeff 
  kappa00 <- expit(eta2)*(1- expit(eta1))*( expit(eta1) - expit(eta2) + 1)
  
  
  Newdatwide$As1 <- rep(0, nrow(Newdatwide))
  Newdatwide$Ar1 <- rep(1, nrow(Newdatwide))
  eta1 <- zeta[1] - as.matrix(Newdatwide) %*% coeff 
  eta2 <- zeta[2] - as.matrix(Newdatwide) %*% coeff 
  kappa01 <- expit(eta2)*(1- expit(eta1))*( expit(eta1) - expit(eta2) + 1)
  w_d <- (w_c* kappa11 * kappa10 * kappa01 * kappa00)/kappa
  return(w_d/mean(w_d))
}


kappa2 <- function(datHH, coeff, zeta, w_c){
  datHH$Asx1 <- datHH$As1*datHH$x11
  datHH$Arx1 <- datHH$Ar1*datHH$x11r
  datHH$AsAr1 <- datHH$As1*datHH$Ar1
  datHH$AsArHH1 <- datHH$As1*datHH$Ar1*datHH$x12H
  
  datHH$Asx2 <- datHH$As2*datHH$x21
  datHH$Arx2 <- datHH$Ar2*datHH$x21r
  datHH$AsAr2 <- datHH$As2*datHH$Ar2
  datHH$AsArHH2 <- datHH$As2*datHH$Ar2*datHH$x22H
  
  pre_order1 <- c("x11",  "x11r" , "x12H",  "As1", "Asx1", "Ar1", "Arx1","AsAr1", "AsArHH1")
  pre_order2 <- c("x21",  "x21r" , "x22H",  "As2", "Asx2", "Ar2", "Arx2","AsAr2", "AsArHH2" )
  Newdatwide <- datHH[, c(pre_order1, pre_order2)]
  eta1 <- zeta[1] - as.matrix(Newdatwide) %*% coeff 
  eta2 <- zeta[2] - as.matrix(Newdatwide) %*% coeff 
  kappa <- expit(eta2)*(1- expit(eta1))*( expit(eta1) - expit(eta2) + 1)
  
  
  Newdatwide$As2 <- rep(1, nrow(Newdatwide))
  Newdatwide$Ar2 <- rep(1, nrow(Newdatwide))
  eta1 <- zeta[1] - as.matrix(Newdatwide) %*% coeff 
  eta2 <- zeta[2] - as.matrix(Newdatwide) %*% coeff 
  kappa11 <- expit(eta2)*(1- expit(eta1))*( expit(eta1) - expit(eta2) + 1)
  
  Newdatwide$Ar2 <- rep(0, nrow(Newdatwide))
  eta1 <- zeta[1] - as.matrix(Newdatwide) %*% coeff 
  eta2 <- zeta[2] - as.matrix(Newdatwide) %*% coeff 
  kappa10 <- expit(eta2)*(1- expit(eta1))*( expit(eta1) - expit(eta2) + 1)
  
  Newdatwide$As2 <- rep(0, nrow(Newdatwide))
  eta1 <- zeta[1] - as.matrix(Newdatwide) %*% coeff 
  eta2 <- zeta[2] - as.matrix(Newdatwide) %*% coeff 
  kappa00 <- expit(eta2)*(1- expit(eta1))*( expit(eta1) - expit(eta2) + 1)
  
  
  Newdatwide$As2 <- rep(0, nrow(Newdatwide))
  Newdatwide$Ar2 <- rep(1, nrow(Newdatwide))
  eta1 <- zeta[1] - as.matrix(Newdatwide) %*% coeff 
  eta2 <- zeta[2] - as.matrix(Newdatwide) %*% coeff 
  kappa01 <- expit(eta2)*(1- expit(eta1))*( expit(eta1) - expit(eta2) + 1)
  w_d <- (w_c* kappa11 * kappa10 * kappa01 * kappa00)/kappa
  return(w_d/mean(w_d))
}

n <- 3000
r <- 500
Qstage2 <- matrix(NA, nrow = r, ncol = 6)
colnames(Qstage2) <- c("xi0_hat","xi1_hat", "psi0_hat","psi1_hat", "phi0_hat","phi1_hat")

Qstage1 <- matrix(NA, nrow = r, ncol = 6)
colnames(Qstage1) <- c("xi0_hat","xi1_hat", "psi0_hat","psi1_hat", "phi0_hat","phi1_hat")

stage2 <- matrix(NA, nrow = r, ncol = 6)
colnames(stage2) <- c("xi0_hat","xi1_hat", "psi0_hat","psi1_hat", "phi0_hat","phi1_hat")

stage1 <- matrix(NA, nrow = r, ncol = 6)
colnames(stage1) <- c("xi0_hat","xi1_hat", "psi0_hat","psi1_hat", "phi0_hat","phi1_hat")

stage2 <- matrix(NA, nrow = r, ncol = 6)
colnames(stage2) <- c("xi0_hat","xi1_hat", "psi0_hat","psi1_hat", "phi0_hat","phi1_hat")

stage1 <- matrix(NA, nrow = r, ncol = 6)
colnames(stage1) <- c("xi0_hat","xi1_hat", "psi0_hat","psi1_hat", "phi0_hat","phi1_hat")

for (i in 1:r) {
  x11 <- rnorm(n,0,1)
  x12 <- rbinom(n,1, 0.5)
  x11r <- rnorm(n,0,1)
  x12r <- rbinom(n,1, 0.5)
  AA1 <- genTrt(x11, x12, x11r, x12r)
  colnames(AA1) <- c('As1', 'Ar1')
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
  
  
  dat2s <- data.frame(hh, x21, x22, A2 = AA2$As2)
  dat2r <- data.frame(hh, x21 = x21r, x22 = x22r, A2 = AA2$Ar2)
  dat2 <- rbind(dat2s, dat2r)
  datwide2 <- data.frame(hh, x21, x22, x21r, x22r, AA2)
  datHH2 <- data.frame(hh, x21, x21r, x22H =  (x22 + x22r), AA2)
  
  
#  AA1opt.val <- decision.value(c(-0.5, 1),c(-0.5,1), c(1, -0.5), cbind(1, x11), cbind(1,x11r),
#                               cbind(1,x12+x12r))
 # AA1opt <- t(apply(AA1opt.val, 1, decision))
  
  
  
  
#  AA2opt.val <- decision.value(c(-0.25, 0.5),c(-0.25,0.5), c(1, -0.5), cbind(1, x21), cbind(1,x21r),
#                               cbind(1,x22+x22r))
#  AA2opt <- t(apply(AA2opt.val, 1, decision))
  # regrets
  
  
  # cos(datwide1$x11) + sin(datwide1$x11r)+ 0.5*exp(datwide2$x21 + datwide2$x21r) - 0.5*exp(datwide2$x21 + datwide2$x21r)
  trtfree <-  0.5*cos(pi*(datwide1$x11 + datwide1$x11r)) + 0.5* exp(datwide2$x21 + datwide2$x21r) + 0.2*(datwide1$x12 + datwide1$x12r)^2
  #- (datwide1$x11 + datwide1$x11r)^3 + log(abs(1/datwide2$x21)) - datwide2$x21r^2
  gammaS1 <- datwide1$As1 *(-0.25 + 0.5 * datwide1$x11)
  gammaR1 <- datwide1$Ar1 *(-0.25 + 0.5  * datwide1$x11r)
  gammaInt1 <- datwide1$As1 * datwide1$Ar1 *(-0.5 + 0.25 * (datwide1$x12 + datwide1$x12r))
  
  gammaS2 <- datwide2$As2 *(-0.25 + 0.5 * datwide2$x21)
  gammaR2 <- datwide2$Ar2 *(-0.25 + 0.5 * datwide2$x21r)
  gammaInt2 <- datwide2$As2 * datwide2$Ar2 *(-0.5 + 0.25 * (datwide2$x22 + datwide2$x22r))
  y <- trtfree + gammaS1 + gammaR1 + gammaInt1 + gammaS2 + gammaR2 + gammaInt2 
  mu <- y
  
  
  
  probs <- c(0.6, 0.25, 0.15)
  
  datHH2$OrdY <- genOrdOut(probs, mu)
  ###################################
  
  # analysis
  # models to be passed to dWGLM
  
  datHH <- data.frame(datHH1, datHH2)
  ob.add <- table(datHH2$OrdY)[3]/(table(datHH2$OrdY)[1] + table(datHH2$OrdY)[2])
  dataInd <- data.frame(dat1, dat2)
  dataInd <- dataInd[order(dataInd$hh), ]
  
  
  
  ##########
  indS <- seq(1, nrow(dataInd), by = 2)
  indR <- seq(2, nrow(dataInd), by = 2)
  #   Pmar <-  expit(-1 + 1.15*exp(x1) - 0.5*x2)
  alpha <- glm(A2~ exp(x21) + x22, family = 'binomial', dataInd)
  pi_hat <- predict(alpha, type = "response")
  MarS <- pi_hat[indS]
  MarR <- pi_hat[indR]
  A2 <- dataInd$A2
  # weights
  w <- abs(A2 - fitted(alpha))
  
  zz2 <- aggregate(dataInd$x22, list(dataInd$hh), sum)
  dataInd$z2 <- rep(zz2$x, each=2)
  
  theta.des <- model.matrix( ~ z2, data=dataInd)
  estOR2 <- binomial.twostage(alpha, data = dataInd, var.link = 1,
                              clusters = dataInd$hh,theta.des=theta.des)
  estOR <- exp(estOR2$theta[1] + estOR2$theta[2] * zz2$x)
  tranA <- tranAA(dataInd$A2[indS], dataInd$A2[indR])
  
  
  opw <- weights(MarS, MarR, estOR, tranA)$op
  ##########
  #second stage:
  # Q-learning
  Qclm1 <- clm(OrdY ~ x11 + x11r + x12H + As1 + I(As1*x11) + Ar1 + I(Ar1*x11r) + I(As1 * Ar1) + I(As1 * Ar1 *x12H) +
                x21 + x21r + x22H + As2 + I(As2*x21) + Ar2 + I(Ar2*x21r) + I(As2 * Ar2) + I(As2 * Ar2 *x22H) , data = datHH)
  Qpar <- tail(Qclm1$coefficients,6)
  Qstage2[i,] <- Qpar
  Qdval <- decision.value(Qpar[1:2],Qpar[3:4], Qpar[5:6], cbind(1, datHH$x21), cbind(1,datHH$x21r),
                         cbind(1, datHH$x22H))
  Qaopt <- t(apply(Qdval, 1, decision))
  Qcoeff <-   Qclm1$coefficients[3:length(Qclm1$coefficients)]
  Qzeta <-   Qclm1$coefficients[1:2]
  
  
  clm1 <- clm(OrdY ~ x11 + x11r + x12H + As1 + I(As1*x11) + Ar1 + I(Ar1*x11r) + I(As1 * Ar1) + I(As1 * Ar1 *x12H) +
                x21 + x21r + x22H + As2 + I(As2*x21) + Ar2 + I(Ar2*x21r) + I(As2 * Ar2) + I(As2 * Ar2 *x22H) , weights = opw, data = datHH)
  coeff <- clm1$coefficients[3:length(clm1$coefficients)]
  zeta <- clm1$coefficients[1:2]
  
  w_d <- kappa2(datHH, coeff, zeta, opw)
  #######################
  clm2 <- clm(OrdY ~ x11 + x11r + x12H + As1 + I(As1*x11) + Ar1 + I(Ar1*x11r) + I(As1 * Ar1) + I(As1 * Ar1 *x12H) +
                x21 + x21r + x22H + As2 + I(As2*x21) + Ar2 + I(Ar2*x21r) + I(As2 * Ar2) + I(As2 * Ar2 *x22H) , weights = w_d, data = datHH)
  par <- tail(clm2$coefficients,6)
  stage2[i,] <- par
  coeffS2 <- clm2$coefficients[3:length(clm2$coefficients)]
  zetaS2 <- clm2$coefficients[1:2]
  
  dval <- decision.value(par[1:2],par[3:4], par[5:6], cbind(1, datHH$x21), cbind(1,datHH$x21r),
                         cbind(1, datHH$x22H))
  aoptS2 <- t(apply(dval, 1, decision))
  
 
  # Pmar <-  expit(-1 + 1.15*exp(x1) - 0.5*x2)
  alpha1 <- glm(A1~ exp(x11) + x12, family = 'binomial', dataInd)
  pi_hat1 <- predict(alpha1, type = "response")
  MarS1 <- pi_hat1[indS]
  MarR1 <- pi_hat1[indR]
  A1 <- dataInd$A1
  # weights
  w1 <- abs(A1 - fitted(alpha1))
  
  zz1 <- aggregate(dataInd$x12, list(dataInd$hh), sum)
  dataInd$z1 <- rep(zz1$x, each=2)
  
  
  theta.des <- model.matrix( ~ z1, data=dataInd)
  estOR2 <- binomial.twostage(alpha1, data = dataInd, var.link = 1,
                              clusters = dataInd$hh,theta.des=theta.des)
  estOR <- exp(estOR2$theta[1] + estOR2$theta[2] * zz1$x)
  tranA <- tranAA(dataInd$A1[indS], dataInd$A1[indR])
  
  opwS1 <- weights(MarS1, MarR1, estOR, tranA)$op
  stage11 <- matrix(NA, nrow = 15, ncol = 6)
  Qstage11 <- matrix(NA, nrow = 15, ncol = 6)
  oddQ <- rep(NA, 15)
  oddP <- rep(NA, 15)
  for (ii in 1:15) {
    QpredP <- predCLM(datHH, Qcoeff, Qzeta, Qaopt)
    Qcat <- simstudy:::matMultinom(as.matrix(QpredP))
    datHH1$QOrdY <- ordered(Qcat)
    Qclm2 <- clm(QOrdY ~ x11 + x11r + x12H + As1 + I(As1*x11) + Ar1 + I(Ar1*x11r) + I(As1 * Ar1) + I(As1 * Ar1 *x12H), data = datHH1)
    Qstage11[ii,] <- tail(Qclm2$coefficients,6)
    oddP[ii] <- table(datHH1$QOrdY)[3]/(table(datHH1$QOrdY)[1] + table(datHH1$QOrdY)[2])
    
    predP <- predCLM(datHH, coeffS2, zetaS2, aoptS2)
    cat <- simstudy:::matMultinom(as.matrix(predP))
    datHH1$OrdY <- ordered(cat)
    oddQ[ii] <- table(datHH1$OrdY)[3]/(table(datHH1$OrdY)[1] + table(datHH1$OrdY)[2])
    clmS1 <- clm(OrdY ~ x11 + x11r + x12H + As1 + I(As1*x11) + Ar1 + I(Ar1*x11r) + I(As1 * Ar1) + I(As1 * Ar1 *x12H) , weights = opwS1, data = datHH1)
    coeffS1 <- clmS1$coefficients[3:length(clmS1$coefficients)]
    zetaS1 <- clmS1$coefficients[1:2]
    w_dS1 <- kappa1(datHH1, coeffS1, zetaS1, opwS1)
    #######################
    clm2 <- clm(OrdY ~ x11 + x11r + x12H + As1 + I(As1*x11) + Ar1 + I(Ar1*x11r) + I(As1 * Ar1) + I(As1 * Ar1 *x12H), weights = w_dS1, data = datHH1)
    stage11[ii,] <- tail(clm2$coefficients,6)
  }
  stage1[i,] <- apply(stage11,2, mean)
  Qstage1[i,] <- apply(Qstage11,2, mean)
  Qodd.rat[i] <- mean(oddQ)/ob.add; Podd.rat[i] <- mean(oddP)/ob.add
}

apply(stage1,2, mean)
apply(stage2,2, mean)



library(vioplot)
par(mfrow=c(3,2))
vioplot(Qstage1[,1], stage1[,1], 
        names=c( "Q-learning","dWPOM" ),xlab= "H = 1000, Stage 1", ylab = expression(paste( xi[0], " estimates")))
abline(h = -0.25, col = "red")

vioplot(Qstage1[,2], stage1[,2], 
        names=c( "Q-learning","dWPOM" ),xlab= "H = 1000, Stage 1", ylab = expression(paste( xi[1], " estimates")))
abline(h = 0.5, col = "red")
vioplot(Qstage1[,3], stage1[,3], 
        names=c( "Q-learning","dWPOM" ),xlab= "H = 1000, Stage 1",  ylab = expression(paste( psi[0], " estimates")))
abline(h = -0.25, col = "red") 

vioplot(Qstage1[,4], stage1[,4], 
        names=c( "Q-learning","dWPOM" ), xlab= "H = 1000, Stage 1", ylab = expression(paste( psi[1], " estimates")))
abline(h = 0.5, col = "red") 

vioplot(Qstage1[,5], stage1[,5], 
        names=c( "Q-learning","dWPOM" ),xlab= "H = 1000, Stage 1",  ylab = expression(paste( phi[0], " estimates")))
abline(h = -0.5, col = "red") 

vioplot(Qstage1[,6], stage1[,6], 
        names=c( "Q-learning","dWPOM" ), xlab= "H = 1000, Stage 1", ylab = expression(paste( phi[1], " estimates")))
abline(h = 0.25, col = "red") 




library(vioplot)
par(mfrow=c(3,2))
vioplot(Qstage2[,1], stage2[,1], 
        names=c( "Q-learning","dWPOM" ),xlab= "H = 1000, Stage 2", ylab = expression(paste( xi[0], " estimates")))
abline(h = -0.25, col = "red")

vioplot(Qstage2[,2], stage2[,2], 
        names=c( "Q-learning","dWPOM" ),xlab= "H = 1000, Stage 2", ylab = expression(paste( xi[1], " estimates")))
abline(h = 0.5, col = "red")
vioplot(Qstage2[,3], stage2[,3], 
        names=c( "Q-learning","dWPOM" ),xlab= "H = 1000, Stage 2",  ylab = expression(paste( psi[0], " estimates")))
abline(h = -0.25, col = "red") 

vioplot(Qstage2[,4], stage2[,4], 
        names=c( "Q-learning","dWPOM" ), xlab= "H = 1000, Stage 2", ylab = expression(paste( psi[1], " estimates")))
abline(h = 0.5, col = "red") 

vioplot(Qstage2[,5], stage2[,5], 
        names=c( "Q-learning","dWPOM" ),xlab= "H = 1000, Stage 2",  ylab = expression(paste( phi[0], " estimates")))
abline(h = -0.5, col = "red") 

vioplot(Qstage2[,6], stage2[,6], 
        names=c( "Q-learning","dWPOM" ), xlab= "H = 1000, Stage 2", ylab = expression(paste( phi[1], " estimates")))
abline(h = 0.25, col = "red") 

mean(Qodd.rat)
mean(Podd.rat)