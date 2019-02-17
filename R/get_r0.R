#' Reference bias due to alignment
#'
#' This function is designed to estimate the reference ratio in the absence of any regulatory
#' difference; that is, the reference bias due to alignment. This is done by obtaining the
#' median ratio between eh1 and eh2 across the top 20% expressed loci in the library.
#'
#' @param ASEdat Dataframe containing reference and alternative count data.
#' @param eh1 A string with the column name of the reference count data
#' @param eh2 A string with the column name of the alternative count data
#' @return A numeric estimate of the reference bias r0 for the given library.
#' @examples
#' @export

get_r0<-function(ASEdat, eh1 = "refCount", eh2 = "altCount"){
  totalCount<-ASEdat[,eh1]+ASEdat[,eh2]
  indices<-which(totalCount>=quantile(totalCount,.8,na.rm=TRUE))
  r0<-median(ASEdat[indices,eh1]/totalCount[indices])
}
