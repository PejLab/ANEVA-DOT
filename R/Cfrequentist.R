#assume mutation is on eh1

ASEoutlierCfrequentist<-function(filepath, output_columns = c("refCount","altCount"), eh1 = "refCount",
                      eh2 = "altCount", Vg = 0.25, FDR = 0.05, plot = TRUE){
  dat<-read.table(filepath, header = TRUE)
  output<-dat[,output_columns]
  for (i in 1:nrow(dat)){
    if ((dat[i,eh1]+dat[i,eh2])>8){
      #assume M is on eh1
      r.obs<-(dat[i,eh1]/(dat[i,eh1]+dat[i,eh2]))
      f.r.obs<-(-1-log2(1-r.obs))/sqrt(Vg)
      if (f.r.obs<=0){
        output$p.val.eh1[i]<-2*pnorm(f.r.obs)
      }
      else{
        output$pval.eh1[i]<-2*pnorm(-f.r.obs)
      }
      g.r.obs<-(-1-log2(r.obs))/sqrt(Vg)
      if (g.r.obs>=0){
        output$p.val.eh2[i]<-2*pnorm(g.r.obs, lower.tail = FALSE)
      }
      else{
        output$pval.eh2[i]<-2*pnorm(-g.r.obs, lower.tail = FALSE)
      }
      output$pval.2sided[i]<-mean(c(output$pval.eh1,output$pval.eh2))
    }
    else{
      output$p.val.eh1[i]<-NA
      output$p.val.eh2[i]<-NA
      output$p.val.2sided[i]<-NA
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
         main = "Reference Count vs. Alternative Count",
         xlab = "Reference Count", ylab = "Alternative Count",
         col = ifelse(result$adj.pval<.05,'red','black'), pch = 19)
    abline(a = 0, b = 1, col = "blue")
  }
  return(output)
}
