
#The following function performs binomial tests on each observation of reference and alternative
#count data (eh1 and eh2) taken from a raw .txt file of ASE data. Filepath for the raw file must 
#be entered in quotation marks (""). The default binomial parameter to test is prob=0.5, while the 
#default FDR=0.10. The function will return the user-specified output_columns, in addition to the
#accompanying p-values and adjusted p-values using Benjamini-Hoschberg method. Observations with
#a total count (eh1+eh2) less than 8 are assigned a p-value of NA. 
ASE.binom.test<-function(filepath, output_columns, eh1 = "refCount", eh2 = "altCount", prob = 0.5, 
                         FDR = 0.05, plot = TRUE){
  dat<-read.table(filepath, header = TRUE)
  output<-dat[,output_columns]
  for (i in 1:nrow(dat)){
    if ((dat[i,eh1]+dat[i,eh2])>8){
      test<-binom.test(c(dat[i,eh1],dat[i,eh2]), p = prob)
      output$p.val[i]<-test$p.value
    }
    else{output$p.val[i]<-NA}
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

#function returns a plot of the ASE reference and alternative counts. Counts have been log transformed
#Observations with 0 values for eh1 or eh2 are replaced with .1 to facilitate the log transformation.
plot.ASE.binom<-function(filepath, eh1 = "refCount", eh2 = "altCount", prob = 0.5, FDR = 0.05){
  result<-ASE.binom.test(filepath = filepath, output_columns = c(eh1,eh2), prob = prob, FDR = FDR)
  result[which(result[,eh1] == 0),eh1]<-.1
  result[which(result[,eh2] == 0),eh2]<-.1
  plot(result[,eh1], result[,eh2], log="xy", main = "Reference Count vs. Alternate Count",
       xlab = "Reference Count", ylab = "Alternative Count",
       col = ifelse(result$adj.pval<.10,'red','black'), pch = 19)
  abline(a = 0, b = 1, col = "blue")
}