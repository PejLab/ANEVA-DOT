#' This function generates the binomial pdf, with the option to provide
#' Binomial coefficients to it. It is useful when it is recalculated on
#' the same N, over and over again.
#' log_BinCoeffs can be precalculated using the function log_BinCoeffs(n)
#' NEED TO INSERT WARNING IF p OUT OF DOMAIN [0,1]
#'
#' @param n The desired n for which the pdf will be calculated
#' @param p The Binomial parameter probability of success
#' @param log_BinCoeffs Optional vector of log-scale binomial coefficients
#' @return Bin_p The pdf for 0:n, out of n

pdf_Binom_fast<-function(n,p,log_BinCoeffs){
  if (p==0){
    Bin_p<-numeric(n+1)
    Bin_p[1]<-1
    return(Bin_p)
  }
  else if (p==1){
    Bin_p<-numeric(n+1)
    Bin_p[n+1]<-1
    return(Bin_p)
  }
  else{
    log_p<-log(p)
    log_pc<-log(1-p)
    x<-0:n
    xc<-n-x
    Bin_p<-exp(log_p*x+log_pc*xc+log_BinCoeffs)
    Bin_p<-(Bin_p/sum(Bin_p))
    return(Bin_p)
  }
}
