# This test is designed to detect if ASE shows sufficient imbalance to induce
# an outlier in total gene expression. We test for both loss of expression
# and gain of expression, as well as a two-sided test.


#' We first implement a helper function to generate the desired integrand
#' for the Loss of Expression test, assuming that the haplotype with the
#' higher count is the wild type.
#'
#' @param r The ?true reference allele ratio?
#' @param eh1 An integer value for one of the haplotype counts
#' @param eh2 An integer value for the second haplotype count
#' @param Vg A numeric value for the estimated dosage variance in normal population
#' (can be gotten from GTEx)
#' @return The integrand for which we wish to integrate over to get
#' desired p-values
#' @examples (tbd)

integrand.UndExp<-function(x, eh1, eh2, Vg){
  Hh<-max(eh1,eh2)
  Lh<-min(eh1,eh2)
  k<-x/(1-x)
  dt<-log2(k+1)-1
  pnorm(dt,mean = 0, sd=sqrt(Vg))*dbeta(x,Lh+1,Hh+1)
}

#' We next implement a helper function to generate the desired integrand
#' for the Gain of Expression test, assuming that the haplotype with the
#' lower count is the wild type.
#'
#' @param r The ?true reference allele ratio?
#' @param eh1 An integer value for one of the haplotype counts
#' @param eh2 An integer value for the second haplotype count
#' @param Vg A numeric value for the estimated dosage variance in normal population
#' (can be gotten from GTEx)
#' @return The integrand for which we wish to integrate over to get
#' desired p-values
#' @examples (tbd)

integrand.OvrExp<-function(x, eh1, eh2, Vg){
  Hh<-max(eh1,eh2)
  Lh<-min(eh1,eh2)
  k<-x/(1-x)
  dt<-1-log2(k+1)
  pnorm(dt,mean = 0, sd=sqrt(Vg))*dbeta(x,Lh+1,Hh+1)
}


#'
#' Bayesian Change of Expression test on ASE count data.
#'
#' @param filepath A string with quotation marks around it indicating ASE data filepath
#' @param output_columnns A vector containing strings of the exact column names to be duplicated
#' from input file in the output
#' @param r The ?true reference allele ratio?
#' @param eh1 A string with the column name of the reference count data
#' @param eh2 A string with the column name of the alternative count data
#' @param Vg A numeric value for the estimated dosage variance in normal population
#' @param FDR A numeric value between 0 and 1 indicating the desired false discovery rate
#' @param plot A logical T/F indicating whether plots should be generated for the test data
#' @return A table containing the user-specified output columns, as well as unadjusted and
#' Benjamini-Hoschberg adjusted p-values for Loss, Gain, and change of expression
#' @examples (tbd)
#' @export
ASEoutlier<-function(filepath, output_columns = c("refCount","altCount"), r, eh1 = "refCount",
                      eh2 = "altCount", Vg, FDR = 0.05, plot = TRUE){
  dat<-read.table(filepath, header = TRUE)
  output<-dat[,output_columns]
  for (i in 1:nrow(dat)){
    if ((dat[i,eh1]+dat[i,eh2])>8){
      #output$p.val.LoE[i]<-2*integrate(dp.UndExp(r,dat[i,eh1],dat[i,eh2],Vg),0,1)
      #output$p.val.GoE[i]<-2*integrate(dp.OvrExp(r,dat[i,eh1],dat[i,eh2],Vg),0,1)
      #output$p.val.2sided[i]<-(0.5*output$p.val.LoE)+(0.5*output$p.val.GoE)
      #integrand<-function(x){
      #  pnorm(log2((x/(1-x))+1)-1,0,sqrt(Vg))*dbeta(x,min(dat[i,eh1],dat[i,eh2]),max(dat[i,eh1],dat[i,eh2]))
      #  }
      integrand<-function(x){pnorm(log2((x/(1-x))+1)-1,0,sqrt(Vg))*dbeta(x,min(dat[i,eh1],dat[i,eh2]),max(dat[i,eh1],dat[i,eh2]))}
      #integral.result<-integrate(integrand,0,1)
      output$p.val.LoE[i]<-2*integrate(integrand,0,1)$value
    }
    else{
      output$p.val.LoE[i]<-NA
      #output$p.val.GoE[i]<-NA
      #output$p.val.2sided[i]<-NA
      }
  }
  #Carry out Benjamini-Hochberg procedure to get adjusted p-values
  output$adj.pval<-p.adjust(output$p.val,method = "BH")
  #output$sign<-(output$adj.pval<FDR)
  if (plot==TRUE){
    result<-output
    result[which(result[,eh1] == 0),eh1]<-.1
    result[which(result[,eh2] == 0),eh2]<-.1
    plot(result[,eh1], result[,eh2], log="xy", #yaxt = "n", xaxt = "n",
         main = "Reference Count vs. Alternate Count",
         xlab = "Reference Count", ylab = "Alternative Count",
         col = ifelse(result$adj.pval<.10,'red','black'), pch = 19)
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
