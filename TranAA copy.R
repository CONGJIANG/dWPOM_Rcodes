tranAA <- function(AS, AR){
  tran <- rep(3, length(AS))
  tran[AS == 1] <- 2
  tran[(AS + AR) == 2] <- 1
  tran[(AS + AR) == 0] <- 4
  return(tran)
}
