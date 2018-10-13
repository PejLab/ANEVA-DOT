#' ANEVA Dosage Outlier Test
#'
#' This test is designed to detect if ASE data reveals sufficient allelic imbalance to induce
#' an outlier in total gene expression.
#'
#' @param ASEdat Dataframe containing reference and alternative count data.
#' @param output_columns A vector containing additional ASEdat column name strings which the user
#' wishes to duplicate in the output.
#' @param eh1 A string with the column name of the reference count data
#' @param eh2 A string with the column name of the alternative count data
#' @param Eg_std a vector containing standard deviation in the log2 transformed total gene
#' expression in a healthy population. Eg_std vector must be in one-to-one correspondence with
#' ASE count data, and must be ordered correctly. Can be a vector in ASEdat if present. P-values
#' will not be generated for records with missing or infinite variances.
#' @param r0 The ratio of the eh1 allele (i.e., eh1/(eh1+eh2)) in the absence of any regulatory
#' difference (reference bias due to alignment). The simplest way to get such an estimate would
#' be to get the median ratio between eh1 and eh2 across the entire library. Can be a single
#' value across entire library, or a SNP-wise vector.
#' @param p0 An average noise rate p(R->A) or p(A->R), i.e., the probability of seeing an allele
#' due to noise when it is essentially not there. (For v7: LAMP = 0.0003) Can be a single
#' value across entire library, or a SNP-wise vector.
#' @param FDR A numeric value between 0 and 1 indicating the desired false discovery rate.
#' Default FDR is 0.05.
#' @param coverage A numeric value such that if total allelic count is less than that value,
#' p-values will not be generated for that record. Default value is 8.
#' @param plot A logical T/F indicating whether plots should be generated for the test data
#' @return A table containing the user-specified output_columns, as well as unadjusted and
#' adjusted p-values for detection of potential dosage outlier. P-values are adjusted using
#' Benjamini-Hoschberg method. P-values are not generated for records with missing or infinite
#' variances.
#@examples result<-ANEVAdot("data/testdata.txt", output_columns = c("eh1","eh2"), eh1 = "eh1", eh2 = "eh2", Eg_std=Sg)
#' @export
ANEVA_DOT<-function(ASEdat, output_columns = c("refCount","altCount"), eh1 = "refCount",
                   eh2 = "altCount", Eg_std, r0 = 0.5, p0 = 0.0003, FDR = 0.05,
                   coverage = 8, plot = TRUE){
  output<-ASEdat[,output_columns]
  #insert user warning about r0 and p0 defaults
  if (length(r0)==1){
    r0<-rep(r0,nrow(ASEdat))
  }
  if (length(p0)==1){
    p0<-rep(p0,nrow(ASEdat))
  }
  for (i in 1:nrow(ASEdat)){
    if (!is.finite(Eg_std[i])){
      output$p.val[i]<-NA
      next
    }
    else if ((ASEdat[i,eh1]+ASEdat[i,eh2])>coverage){
      output$p.val[i]<-Test_ASE_Outliers(Eg_std[i],ASEdat[i,eh1],ASEdat[i,eh2],r0[i],p0[i])
    }
    else{
      output$p.val[i]<-NA #due to inadequate coverage
    }
  }
  #Carry out Benjamini-Hochberg procedure to get adjusted p-values
  output$adj.pval<-p.adjust(output$p.val,method = "BH")

  #Carry out plotting if called for
  if (plot==TRUE){
    result<-output
    result[which(result[,eh1] == 0),eh1]<-.1
    result[which(result[,eh2] == 0),eh2]<-.1
    plot(result[,eh1], result[,eh2], log="xy", #yaxt = "n", xaxt = "n",
         main = "Reference Count vs. Alternative Count",
         xlab = "Reference Count", ylab = "Alternative Count",
         col = ifelse(result$adj.pval<.05,'red','black'), pch = 19)
    abline(a = 0, b = 1, col = "blue")
  }
  return(output)
}


#' @param Eg_std Standard deviation in the log2 transformed total gene expression in a healthy
#' population.
#' @param eh1 Integer value for of expression of haplotype 1.
#' @param eh2 Integer value of expression of haplotype 2.
#' @param r0 The ratio of the eh1 allele (i.e., eh1/(eh1+eh2)) in the absence of any regulatory
#' difference (reference bias due to alignment). The simplest way to get such an estimate would
#' be to get the median ratio between eh1 and eh2 across the entire library.
#' @param p0 An average noise rate p(R->A) or p(A->R), i.e., the probability of seeing an allele
#' due to noise when it is essentially not there. (For v7: LAMP = 0.0003).
#' @return Undajusted 2-sided p-value, testing whether: H0:There does not exist a significantly
#' abnormal allelic imbalance for this SNP, vs. H1: Allelic ratio significantly deviates from
#' normal population.
Test_ASE_Outliers<-function(Eg_std, eh1, eh2, r0, p0){
  if (eh1==eh2){
    p.val<-1
  }
  Eg_std<-max(Eg_std,2^(-52)) #this variance is already zero when compared to binomial variance.
  rad<-Eg_std*4 #integration radius
  log_BinCoeff<-log_BinCoeffs(eh1+eh2) #pre-calculate binomial coefficients to avoid redoing them within the integral
  p.val<-integrate(integrand,-rad,rad,eh1,eh2, Eg_std, r0, p0, log_BinCoeff, abs.tol = 0, rel.tol = 1e-4)$value
}


#' Function generates the required integrand for the test.
#'
#' @param dE Variable representing the true change in gene dosage
#' (variable of integration)
#' @param eh1 Integer value for count of expression of haplotype1
#' @param eh2 Integer value for count of expression of haplotype2
#' @param Eg_std Numeric value for the estimated dosage standard deviation
#' in normal population (can be derived from GTEx)
#' @param r0 The ratio of the eh1 allele (i.e., eh1/(eh1+eh2)) in the absence of any regulatory
#' difference (reference bias due to alignment). The simplest way to get such an estimate would
#' be to get the median ratio between eh1 and eh2 across the entire library.
#' @param p0 An average noise rate p(R->A) or p(A->R), i.e., the probability of seeing an allele
#' due to noise when it is essentially not there. (For v7: LAMP = 0.0003).
#' @param log_BinCoeffs a vector of log transformed binomial coefficients 0:N=eh1+eh2.
#' @return The integrand over which we wish to integrate to get
#' desired p-values

integrand<-function(dE, eh1, eh2, Eg_std, r0, p0, log_BinCoeff){
  N<-eh1+eh2
  Prob.dE<-dnorm(dE,0,Eg_std) #probability of observed dE in population
  dE[dE<log(.5)]<-(log(.5)) #-1 is the minimum possible loss of expression in log2
  kr<-2*exp(dE)-1 #regulatory effect size of the mutation
  rr<-kr/(kr+1) #allelic ratio of the mutated haplotype
  rn<-rr+p0*(1-2*rr) #allelic ratio of the mutated haplotype after accounting for noise
  kn<-rn/(1-rn) #effective effect size after all is said and done
  k0<-r0/(1-r0) #artificial effect size due to bias

  #Mutation on reference haplotype
  k<-kn*k0 #total expected aFC of R hap. to A hap. after accounting for bias
  r_mR<-k/(k+1) #reference allelic expression ratio (ASE Ref hap)

  k<-k0/kn #total expected aFC of R hap. to A hap. after accounting for bias
  r_mA<-k/(k+1)
  ifelse(k==Inf,r_mA<-1,r_mA<-(k/(k+1)))
  ifelse(k0==0 && kn==0,r_mA<-0.5,r_mA<-(k/(k+1)))

  return(Binom_test_ctm_dbl(eh1,N,r_mR,r_mA,log_BinCoeff,r0)*Prob.dE)
}

#For the following test we assume X and N are scalars while p is a vector
Binom_test_ctm_dbl<-function(X,N,p1,p2,log_BinCoeff,r0){
  m<-round(r0*N)
  p.val<-numeric(length(p1))
  for (i in 1:(length(p1))){
    if (X==m){
      p.val[i]<-1
    }
    else{
      Bnp1<-pdf_Binom_fast(N,p1[i],log_BinCoeff)
      Bnp2<-pdf_Binom_fast(N,p2[i],log_BinCoeff)
      Bnp<-(Bnp1+Bnp2)/2

      tpl<-sum(Bnp[1:X])
      tpr<-sum(Bnp[(X+2):(length(Bnp))])
      p.val[i]<-2*min(tpl,tpr,Bnp[X+1])
    }
  }
  return(p.val)
}

#' #################################################################################################
#' Generate plot of ASE data showing significant outliers. (This Function not currently implemented)
#'
#' @param filepath A string with quotation marks around it indicating ASE data filepath
#' @param eh1 A string with the column name of the reference count data
#' @param eh2 A string with the column name of the alternative count data
#' @param Eg_std a vector with variance estimates for dosage distributions associated with count data.
#' Eg_std must be in one-to-one correspondence with ASE data, and ordered correctly.
#' @param FDR A numeric value between 0 and 1 indicating the desired false discovery rate
#' @examples (tbd)
#' @export
plot.ANEVAdot<-function(filepath, eh1 = "refCount", eh2 = "altCount", Eg_std, FDR = 0.05){
  result<-ANEVAdot(filepath = filepath, output_columns = c(eh1,eh2), Eg_std = Eg_std, FDR = FDR, plot = FALSE)
  result[which(result[,eh1] == 0),eh1]<-.1
  result[which(result[,eh2] == 0),eh2]<-.1
  plot(result[,eh1], result[,eh2], log="xy", main = "Reference Count vs. Alternative Count",
       xlab = "Reference Count", ylab = "Alternative Count",
       col = ifelse(result$adj.pval<.10,'red','black'), pch = 19)
  abline(a = 0, b = 1, col = "blue")
}
