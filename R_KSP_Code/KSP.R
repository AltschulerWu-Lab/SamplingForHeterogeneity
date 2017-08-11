# Define KS-Prime function ----

# This function quantifies the similarity  (in terms of KS-prime and KS measures) of 
# subsample distribution to that of the whole. 
# Input:
# 1) wholeData - the whole-sample data as a vector (each element representing a single cell measurement)
# 2) subsetIndexList - A list of subsampling indices: each element of the list represents a collection of cells sub-sample extracted 
# from the wholeData. These subsamples are supplied as vectors containing the indices (in wholeData) of the cells belonging to the subsamples.
# Output: 
# Res -  is a list containing the KS-primes and the KS for the distributions of the different subsamples.
# Each of these is a vector of the same length as subsetIndexList.
ksPrime<-function(wholeData,subsetIndexList) {
  wholeDataRanks<-ceiling(rank(wholeData, na.last = "keep",ties.method = "average"))+1
  wholeDataRanks[which(is.na(wholeDataRanks))]<-1;
  wholeHist<-hist(wholeDataRanks,(1:(length(wholeDataRanks)+2)-0.5),plot=F)
  #wholeHist<-tail(wholeHist,length(wholeHist)-1);
  wholePDF<-wholeHist$counts
  wholePDF<-tail(wholePDF,length(wholePDF)-1);
  wholeCDF<-cumsum(wholePDF)
  wholeCDF<-wholeCDF/tail(wholeCDF,1)
  scaleFactor<-sqrt(wholeCDF*(1-wholeCDF))
  scaleFactor[which(scaleFactor==0)]<-Inf;
  
  nList<-length(subsetIndexList);
  ks<- as.numeric(rep(NA, nList)); 
  ksp<- as.numeric(rep(NA, nList)); 
  for(listNumber in 1:nList){ 
  subsetHist<-hist(wholeDataRanks[subsetIndexList[[listNumber]]],(1:(length(wholeDataRanks)+2)-0.5),plot=F);
  #subsetHist<-tail(subsetHist,length(subsetHist)-1);
  subsetPDF<-subsetHist$counts;
  subsetPDF<-tail(subsetPDF,length(subsetPDF)-1);
  subsetPDF<-subsetPDF/sum(subsetPDF);
  subsetCDF<-cumsum(subsetPDF);
  ks[listNumber]=max(abs(subsetCDF-wholeCDF));
  ksp[listNumber]=0.5*max(abs(subsetCDF-wholeCDF)/scaleFactor);
  }
  res=list("ks"=ks,"ksp"=ksp)
  return(res)
}
