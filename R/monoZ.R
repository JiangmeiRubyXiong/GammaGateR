#' Make cfGMM probabilities monotonic
#'
#' Takes a cfGMM model and forces the right-most component
#' probabilities to be monotonic increasing with respect
#' to the input values.
#'
#' @param fit A cfGMM model object
monoZ <- function(fit){
  z <- fit$posterior
  z1 <- z[,1]; z2 <- z[,2]
  x <- fit$x
  n <- length(x)
  order.x <- order(x)
  rank.x <- rank(x)

  z1 <- z1[order.x]
  z1 <- rev(z1)
  diff1 <- z1[2:n]-z1[1:(n-1)]
  diff1 <- c(min(z1), diff1)
  if(any(diff1<0)){
    diff1[diff1<0] <- 0
    z1 <- cumsum(diff1)
  }
  z1 <- rev(z1)
  z1 <- z1[rank.x]


  z2 <- z2[order.x]
  diff2 <- z2[2:n]-z2[1:(n-1)]
  diff2 <- c(min(z2), diff2)
  if(any(diff2<0)){
    diff2[diff2<0] <- 0
    z2 <- cumsum(diff2)
  }
  z2 <- z2[rank.x]


  z <- cbind(z1, z2)
  fit$posterior <- z
  return(fit)
}
