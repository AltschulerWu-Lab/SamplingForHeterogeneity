library("R.matlab")
source('KSP.R')
# Load Well Data from MATLAB ----
pathname <- file.path("Cell_Culture", "testData.mat")
data <- readMat(pathname)
wellNumbers<-data$plateInfo[[1]]
uW<-unique(wellNumbers);
numberOfWells=length(uW);

# Generate well subsets ----
# Each subset consists of one cells from a single well in this demonstration
subsetIdx <- vector(mode="list", length=numberOfWells);
for(i in 1:length(subsetIdx)){
  subsetIdx[[i]]=which(wellNumbers==uW[i])
}

# calculate KS & KS-prime

numberOfFeatures<-ncol(data$plateData)
ksMat<-matrix(0,nrow=numberOfWells,ncol=numberOfFeatures);
kspMat<-matrix(0,nrow=numberOfWells,ncol=numberOfFeatures);
for(featNum in 1:numberOfFeatures){
  res<-ksPrime(data$plateData[,featNum],subsetIdx);
  ksMat[,featNum]<-res$ks;
  kspMat[,featNum]<-res$ksp;
}


# Save results ----
write.table(kspMat,file="R_ksp.csv",row.names = F,col.names = F);
write.table(ksMat,file="R_ks.csv",row.names = F,col.names = F);