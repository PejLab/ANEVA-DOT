
#' We implement a helper function to generate the desired integrand
#' for the h1M test, that is, assuming that haplotype1 harbors the mutation.
#'
#' @param r Variable representing the true allelic ratio of haplotype 1
#' @param eh1 An integer value for count of expression of haplotype1
#' @param eh2 An integer value for count of expression of haplotype2
#' @param Vg A numeric value for the estimated dosage variance in normal population
#' (can be gotten from GTEx)
#' @return The integrand for which we wish to integrate over to get
#' desired p-values
#' @examples (tbd)

integrand.h1M<-function(r, eh1, eh2, Vg){
  delta<-r/(1-r) #comment about delta, is effect size of mutation.... define r
  delta.E<-log2(delta+1)-1
  exp((pnorm(abs(delta.E),mean = 0, sd=sqrt(Vg), lower.tail=FALSE, log.p = TRUE))+dbeta(r,eh1,eh2, log = TRUE))
}

#' We implement a helper function to generate the desired integrand
#' for the h2M test, that is, assuming haplotype 2 harbors the mutation.
#'
#' @param r Variable representing the true allelic ratio of haplotype 1
#' @param eh1 An integer value for count of expression of haplotype1
#' @param eh2 An integer value for count of expression of haplotype2
#' @param Vg A numeric value for the estimated dosage variance in normal population
#' (can be gotten from GTEx)
#' @return The integrand for which we wish to integrate over to get
#' desired p-values
#' @examples (tbd)

integrand.h2M<-function(r, eh1, eh2, Vg){
  delta<-(1-r)/r
  delta.E<-log2(delta+1)-1
  exp((pnorm(abs(delta.E),mean = 0, sd=sqrt(Vg), lower.tail=FALSE, log.p = TRUE))+dbeta(r,eh1,eh2, log = TRUE))
}


#'
#' Bayesian Change of Expression test on ASE count data.
#'
#' This test is designed to detect if ASE shows sufficient imbalance to induce
#' an outlier in total gene expression. We first test under the assumption that
#' the mutation lies on haplotype 1, then test again assuming that the
#' mutation lies on haplotype 2. These results are averaged, resulting in an
#' effective two-sided test.
#'
#' @param filepath_ASE A string with quotation marks around it indicating ASE count data filepath
#' @param filepath_Vg A string with quotation marks around it indicating filepath for Variance estimates
#' @param output_columnns A vector containing strings of the exact column names to be duplicated
#' from input file in the output
#' @param eh1 A string with the column name of the reference count data
#' @param eh2 A string with the column name of the alternative count data
#' @param Vg  A string with the column name of the variances
#' @param FDR A numeric value between 0 and 1 indicating the desired false discovery rate
#' @param plot A logical T/F indicating whether plots should be generated for the test data
#' @return A table containing the user-specified output columns, as well as unadjusted and
#' Benjamini-Hoschberg adjusted p-values for h1M, h2M and total change of expression.
#' @examples (tbd)
#' @export
ASEoutlier<-function(filepath, output_columns = c("refCount","altCount"), eh1 = "refCount",
                      eh2 = "altCount", Vg = .00837743, FDR = 0.05, plot = TRUE){
  dat<-read.table(filepath, header = TRUE)
  output<-dat[,output_columns]
  for (i in 1:nrow(dat)){
    if ((dat[i,eh1]+dat[i,eh2])>8){
      output$p.val.h1M[i]<-2*integrate(integrand.h1M,0,1,eh1=dat[i,eh1],eh2=dat[i,eh2],Vg =Vg)$value
      output$p.val.h2M[i]<-2*integrate(integrand.h2M,0,1,eh1=dat[i,eh1],eh2=dat[i,eh2],Vg =Vg)$value
      output$p.val.2sided[i]<-(0.5*output$p.val.h1M[i])+(0.5*output$p.val.h2M[i])
    }
    else{
      output$p.val.h1M[i]<-NA
      output$p.val.h2M[i]<-NA
      output$p.val.2sided[i]<-NA
      }
  }
  #Carry out Benjamini-Hochberg procedure to get adjusted p-values
  output$adj.pval.h1M<-p.adjust(output$p.val.h1M,method = "BH")
  output$adj.pval.h2M<-p.adjust(output$p.val.h2M,method = "BH")
  output$adj.pval.2sided<-p.adjust(output$p.val.2sided,method = "BH")
  #output$sign<-(output$adj.pval<FDR)
  if (plot==TRUE){
    result<-output
    result[which(result[,eh1] == 0),eh1]<-.1
    result[which(result[,eh2] == 0),eh2]<-.1
    plot(result[,eh1], result[,eh2], log="xy", #yaxt = "n", xaxt = "n",
         main = "Reference Count vs. Alternative Count",
         xlab = "Reference Count", ylab = "Alternative Count",
         col = ifelse(result$adj.pval.2sided<.05,'red','black'), pch = 19)
    abline(a = 0, b = 1, col = "blue")
  }
  return(output)
}

#' Generate plot of ASE data showing significant outliers.
#'
#' @param filepath A string with quotation marks around it indicating ASE data filepath
#' @param r The ?true reference allele ratio?
#' @param eh1 A string with the column name of the reference count data
#' @param eh2 A string with the column name of the alternative count data
#' @param Vg A numeric value for the estimated dosage variance in normal population
#' @param FDR A numeric value between 0 and 1 indicating the desired false discovery rate
#' @examples (tbd)
#' @export
plot.ASEoutlier<-function(filepath, r, eh1 = "refCount", eh2 = "altCount", Vg, FDR = 0.05){
  result<-ASEoutlier(filepath = filepath, output_columns = c(eh1,eh2), r=r, Vg = Vg, FDR = FDR, plot = FALSE)
  result[which(result[,eh1] == 0),eh1]<-.1
  result[which(result[,eh2] == 0),eh2]<-.1
  plot(result[,eh1], result[,eh2], log="xy", main = "Reference Count vs. Alternative Count",
       xlab = "Reference Count", ylab = "Alternative Count",
       col = ifelse(result$adj.pval<.10,'red','black'), pch = 19)
  abline(a = 0, b = 1, col = "blue")
}
