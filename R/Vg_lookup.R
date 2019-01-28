#goal: write a function that retrieves correct ENS_id from messy data column

#input a vector of ensemble id's with junk data
#output the ensemble id's without the junk

getENS_id<-function(name){
  result<-numeric(length(name))
  for (i in 1:length(name)){
    tmp<-as.character(name[i])
    result[i]<-unlist(strsplit(tmp,"\\."))[1]
  }
  return(result)
}

#goal: write a function that retrieves Vg
#use TESTIS variances for now because the records are most complete.

#input an r object containing variance records, the column name in that
#file containing the ensemble id's to join on, and the vector of ensemble id's
#for which variances are desired
#output a vector of variances tied to those ensemble id's

getVg<-function(Vg, colname="IDs",ENS_id){
  result<-numeric(length(ENS_id))
  indices<-match(ENS_id,Vg[,colname])
  Vg[indices,"MSCLSK"]
  for (i in 1:length(ENS_id)){
    result[i]<-Vg[indices[i],"TESTIS"]
  }
  return(result)
}
