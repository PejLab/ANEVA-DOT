#goal: write a function that retrieves correct ENS_id from messy data column

getENS_id<-function(name){
  name<-as.character(name)
  result<-unlist(strsplit(name,"\\."))[1]
}

#goal: write a function that retrieves Vg

getVg<-function(Vg, colname="IDs",ENS_id){
  index<-which(Vg[,colname]==ENS_id)
  return(Vg[index,"TESTIS"])
}
