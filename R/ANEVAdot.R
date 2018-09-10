
#' We implement a helper function to generate the desired integrand
#' for the test.
#'
#' @param dE Variable representing the true change in gene dosage
#' (variable of integration)
#' @param eh1 An integer value for count of expression of haplotype1
#' @param eh2 An integer value for count of expression of haplotype2
#' @param Eg_std A numeric value for the estimated dosage standard deviation
#' in normal population (can be derived from GTEx)
#' @return The integrand over which we wish to integrate to get
#' desired p-values

integrand<-function(dE, eh1, eh2, Eg_std){
  Hh<-max(eh1,eh2)
  Lh<-min(eh1,eh2)
  Prob.dE<-dnorm(dE,0,Eg_std) #probability of observed dE in population
  dE[dE<(-1)]<-(-1)#-1 is the minimum possible loss of expression in log2
  k<-2^(dE+1)-1 #regulatory effect size of the mutation
  r<-k/(k+1) #true allelic expression ratio of the mutated haplotype
  N<-Lh+Hh
  #The first term below is the probability of observing equal or less of
  #the lower expressed allele (Lh).
  #The second term below is the probability of observing equal or
  #more of the higher expressed allele (Hh)
  return((pbinom(Lh,N,r)+pbinom(Hh-1,N,r,lower.tail = FALSE))*Prob.dE)
}

#'
#' ANEVA Dosage Outlier Test
#'
#' This test is designed to detect if ASE data reveals sufficient imbalance to induce
#' an outlier in total gene expression.
#'
#' @param filepath A string with quotation marks around it indicating ASE count data filepath.
#' Data should be a raw, tab-delimited text file.
#' @param output_columnns A vector containing strings of the exact column names to be duplicated
#' from input file in the output
#' @param eh1 A string with the column name of the reference count data
#' @param eh2 A string with the column name of the alternative count data
#' @param Eg_std a vector containing standard deviation in the log2 transformed total gene
#' expression in a healthy population. Eg_std vector must be in one-to-one correspondence with
#' ASE count data, and must be ordered correctly.
#' @param FDR A numeric value between 0 and 1 indicating the desired false discovery rate
#' @param plot A logical T/F indicating whether plots should be generated for the test data
#' @return A table containing the user-specified output columns, as well as unadjusted and
#' Benjamini-Hoschberg adjusted p-values for detection of potential dosage outlier. P-values
#' will not be generated for records with missing or infinite variances, or for those for
#' which the total allelic count is less than 8.
#' @examples result<-ANEVAdot("data/testdata.txt", output_columns = c("eh1","eh2"), eh1 = "eh1", eh2 = "eh2", Eg_std=Sg)
#' @export
ANEVAdot<-function(filepath, output_columns = c("refCount","altCount"), eh1 = "refCount",
                     eh2 = "altCount", Eg_std, FDR = 0.05, plot = TRUE){
  dat<-read.table(filepath, header = TRUE, sep = "\t")
  output<-dat[,output_columns]
  for (i in 1:nrow(dat)){
    if (!is.finite(Eg_std[i])){
      output$p.val[i]<-NA
      next
    }
    if (dat[i,eh1]==dat[i,eh2]){
      output$p.val[i]<-1
    }
    else if ((dat[i,eh1]+dat[i,eh2])>8){
      Hh<-max(dat[i,eh1],dat[i,eh2]) #higher expressed haplotype
      Lh<-min(dat[i,eh1],dat[i,eh2]) #lower expressed haplotype
      rad<-Eg_std[i]*4 #integration radius
      output$p.val[i]<-integrate(integrand,-rad,rad,eh1=Lh,eh2=Hh, Eg_std = Eg_std[i], abs.tol = 0, rel.tol = 1e-4)$value
    }
    else{
      output$p.val[i]<-NA
    }
  }
  #Carry out Benjamini-Hochberg procedure to get adjusted p-values
  output$adj.pval<-p.adjust(output$p.val,method = "BH")
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

#' Generate plot of ASE data showing significant outliers.
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
