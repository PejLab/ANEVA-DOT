#' This function calculates Binomial coefficients for choosing 0:n out of n,
#' in log scale. It is useful for precalculating the coefficients and
#' providing them to "Binom_test_fast.R" or "pdf_binom_fast.R" when they
#' are redone on the same N, over and over again. This function was
#' designed by Pejman Mohammadi, 2018, Scripps Research, San Diego, CA.
#'
#' @param n The desired n for which coefficients will be calculated
#' @return log_z Vector of log scale Binomial coefficients

log_BinCoeffs<-function(n){
  log_z<-numeric(n+1) #note first and last element will remain 0

  lgamm<-lgamma(1:(n+1)) #precalculate gamma function

  m1<-floor(n/2) #this is the middle. binom coeffs symmetric so we calculate only half
  m2<-ceiling(n/2) #this is the "other" middle
  x<- 1:m1
  xc<- n-x
  log_z[2:(m1+1)]<- (lgamm[n+1]-lgamm[x+1]-lgamm[xc+1])

  log_z[n:(m2+1)]<-log_z[2:(m1+1)] #fill in the other side
  return(log_z)
}

