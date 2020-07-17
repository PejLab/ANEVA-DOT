#' ANEVA Dosage Outlier Test
#'
#' This test is designed to detect if ASE data reveals sufficient allelic imbalance to induce
#' an outlier in total gene expression.
#'
#' @param ASEdat Dataframe containing columns with reference and alternative count data.
#' @param output_columns Vector containing any ASEdat column name strings which the user
#' wishes to duplicate in the output.
#' @param eh1 String containing the column name of the reference count data.
#' @param eh2 String with the column name of the alternative count data.
#' @param Eg_std Vector containing genetic standard deviation, the square root of genetic
#' variation in regulation in general population in natural log scale (numeric value). Variance
#' values from GTEx are available on Github.com/PejLab. For additional utilities and starter code
#' to match Ensembl IDs with population variance estimates from GTEx v7, please see vgdat package.
#' Square root transformation must be applied to GTEx variance estimates before running this test.
#' P-values will not be generated for records with missing or infinite standard deviations. Eg_std
#' vector must be in one-to-one correspondence with ASE count data, and must
#' be ordered correctly.
#' @param r0 The ratio of the eh1 allele (i.e., eh1/(eh1+eh2)) in the absence of any regulatory
#' difference (reference bias due to alignment). The simplest way to get such an estimate would
#' be to get the median ratio between eh1 and eh2 across the entire library (i.e., eh1/(eh1+eh2)).
#' This is done automatically if no user-specified r0 value is detected. Can be a single value
#' across entire library, or a SNP-wise vector.
#' @param p0 An average noise rate p(R->A) or p(A->R), i.e., the probability of seeing an allele
#' due to noise when it is essentially not there (for GTEx v7: LAMP = 0.0003). Can be a single
#' value across entire library, or a SNP-wise vector.
#' @param FDR A numeric value between 0 and 1 indicating the desired false discovery rate.
#' Default FDR is 0.05.
#' @param coverage A numeric value such that if total allelic count is less than this value,
#' p-values will not be generated for that record. Default value is 10.
#' @param plot Logical T/F indicating whether plots should be generated for the test data.
#' @param jobs Number of workers to use for parallelization. Default value is 1.
#' @return Data frame containing the user-specified output_columns, as well as unadjusted and
#' adjusted p-values for detection of potential dosage outlier, testing H0:There does not exist
#' a significantly abnormal allelic imbalance for this SNP, vs. H1: Allelic ratio significantly
#' deviates from normal population. P-values are adjusted using Benjamini-Hoschberg method.
#' P-values are not generated for records with missing or infinite standard deviations.
#' @examples
#' # Define the output columns
#' output_columns <- c("GENE_ID", "TISSUE_ID",  "REF_COUNT", "ALT_COUNT", "TOTAL_COUNT", "NULL_RATIO")
#'
#' # re organize the tables by:
#' # 1: Selecting only genes that have Vg scores available
#' # 2: Reordering ASE data and Vg scores so they align
#' tiss <- "MSCLSK" # The data comes from a skeletal muscle sample
#' covered_genes <- intersect(Vg_GTEx_v7$IDs, sample_ASE$GENE_ID)
#' covered_gene_Vgs <- Vg_GTEx_v7[match(covered_genes, Vg_GTEx_v7$IDs), tiss]
#' covered_gene_ASE_data <- sample_ASE[match(covered_genes, sample_ASE$GENE_ID),]
#'
#' # Take the square root of the Vg scores to the get the Standard Deviation (SDg)
#' covered_gene_SDgs <- sqrt(covered_gene_Vgs)
#'
#' # Run ANEVA-DOT
#' ANEVADOT_scores <- ANEVADOT_test(covered_gene_ASE_data, output_columns = output_columns,
#'                                  eh1 = "REF_COUNT", eh2 = "ALT_COUNT", coverage = 10,
#'                                  r0 = covered_gene_ASE_data$NULL_RATIO,
#'                                  Eg_std = covered_gene_SDgs, plot = TRUE)
#' @import foreach parallel doSNOW
#' @export

ANEVADOT_test<-function(ASEdat, output_columns = c("refCount","altCount"), eh1 = "refCount",
                   eh2 = "altCount", Eg_std, r0 = NULL, p0 = NULL, FDR = 0.05,
                   coverage = 10, plot = TRUE, jobs = 1){
  # FIXME: we need this explicit import to make %dopar% work
  library(foreach)
  output<-ASEdat[,output_columns]
  if (is.null(r0)){
    r0<-get_r0(ASEdat, eh1, eh2)
    warning(paste("There was no r0 provided, I estimate it to be",r0))
  }
  if (is.null(p0)){
    p0<-0.000326
    warning(paste("There was no p0 provided, I estimate it to be",p0,"(obtained from GTEx v7 data)."))
  }
  if (length(r0)==1){
    r0<-rep(r0,nrow(ASEdat))
  }
  if (length(p0)==1){
    p0<-rep(p0,nrow(ASEdat))
  }
  #create progress bar
  pb <- txtProgressBar(min = 0, max = nrow(ASEdat), style = 3)
  #get p-values
  cl <- parallel::makeCluster(jobs)
  doSNOW::registerDoSNOW(cl)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  output$p.val <- foreach(i=1:nrow(ASEdat), .combine='c', .options.snow = opts) %dopar% {
    p.val<-NA
    if (!is.finite(Eg_std[i])){ #check for finite standard deviation
      p.val <- NA;
    }
    else if ((ASEdat[i,eh1]+ASEdat[i,eh2])>=coverage){ #check that minimum coverage has been met
      p.val<-Test_ASE_Outliers(Eg_std[i],ASEdat[i,eh1],ASEdat[i,eh2],r0[i],p0[i])
    }
    else{
      p.val<-NA #NA due to inadequate coverage
    }
    return(p.val)
  }
  parallel::stopCluster(cl)
  close(pb)

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


#' Test ASE Outliers
#'
#' This is a black box function which performs the statistical test on a single SNP.
#'
#' @param Eg_std Genetic standard deviation, the square root of genetic
#' variation in regulation in general population in natural log scale (numeric value).
#' @param eh1 Integer value for of expression of haplotype 1.
#' @param eh2 Integer value of expression of haplotype 2.
#' @param r0 The ratio of the eh1 allele (i.e., eh1/(eh1+eh2)) in the absence of any regulatory
#' difference (reference bias due to alignment). The simplest way to get such an estimate would
#' be to get the median ratio between eh1 and eh2 across the entire library.
#' @param p0 An average noise rate p(R->A) or p(A->R), i.e., the probability of seeing an allele
#' due to noise when it is essentially not there (for GTEx v7: LAMP = 0.0003).
#' @return Undajusted 2-sided p-value, testing whether: H0:There does not exist a significantly
#' abnormal allelic imbalance for this SNP, vs. H1: Allelic ratio significantly deviates from
#' normal population.
Test_ASE_Outliers<-function(Eg_std, eh1, eh2, r0, p0){
  Eg_std<-max(Eg_std,(.Machine$double.eps)) #this variance is already zero when compared to binomial variance.
  rad<-Eg_std*4 #integration radius
  log_BinCoeff<-log_BinCoeffs(eh1+eh2) #pre-calculate binomial coefficients to avoid redoing them within the integral
  if (eh1==eh2){
    p.val<-1
  }
  else{
  p.val<-integrate(integrand,-rad,rad,eh1,eh2, Eg_std, r0, p0, log_BinCoeff, abs.tol = 0, rel.tol = 1e-4, stop.on.error = FALSE)$value
  }
  return(p.val)
}


#' Function generates the required integrand for the test.
#'
#' @param dE Variable representing the true change in gene dosage
#' (variable of integration)
#' @param eh1 Integer value for count of expression of haplotype1
#' @param eh2 Integer value for count of expression of haplotype2
#' @param Eg_std Vector containing genetic standard deviation, the square root of genetic
#' variation in regulation in general population in natural log scale (numeric value).
#' @param r0 The ratio of the eh1 allele (i.e., eh1/(eh1+eh2)) in the absence of any regulatory
#' difference (reference bias due to alignment). The simplest way to get such an estimate would
#' be to get the median ratio between eh1 and eh2 across the entire library.
#' @param p0 An average noise rate p(R->A) or p(A->R), i.e., the probability of seeing an allele
#' due to noise when it is essentially not there. (For v7: LAMP = 0.0003).
#' @param log_BinCoeff a vector of natural log (ln) transformed binomial coefficients 0:N=eh1+eh2.
#' @return The integrand over which we wish to integrate to get
#' desired p-values

integrand<-function(dE, eh1, eh2, Eg_std, r0, p0, log_BinCoeff){
  N<-eh1+eh2
  Prob.dE<-dnorm(dE,0,Eg_std) #probability of observed dE in population
  dE[dE<log(.5)]<-(log(.5)) #-1 is the minimum possible loss of expression in log2 (log base 2)
  kr<-2*exp(dE)-1 #regulatory effect size of the mutation
  rr<-kr/(kr+1) #allelic ratio of the mutated haplotype
  rn<-rr+p0*(1-2*rr) #allelic ratio of the mutated haplotype after accounting for noise
  kn<-rn/(1-rn) #effective effect size after all is said and done
  k0<-r0/(1-r0) #artificial effect size due to bias

  #Mutation on reference haplotype
  k<-kn*k0 #total expected aFC of R hap. to A hap. after accounting for bias
  r_mR<-k/(k+1) #reference allelic expression ratio (ASE Ref hap)

  k<-k0/kn #total expected aFC of R hap. to A hap. after accounting for bias
  r_mA<-ifelse(k==Inf,r_mA<-1,r_mA<-(k/(k+1)))
  #ifelse((k0==0 && kn==0),r_mA<-0.5,r_mA<-(k/(k+1))) #this line not returning properly

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
      Bnp  = Bnp/sum(Bnp) #Just to get rid of potential numerical issues
      tpl<-sum(Bnp[0:X])
      tpr<-ifelse(X==N,tpr<-0,tpr<-sum(Bnp[(X+2):(N+1)])) #to properly handly right tail indices
      p.val[i]<-2*min(tpl,tpr)
      p.val[i]<-p.val[i]+Bnp[X+1]
    }
  }
  return(p.val)
}
