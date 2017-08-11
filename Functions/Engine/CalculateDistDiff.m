function [kspMat,ksMat,deltaMeanMat,deltaMediansMat]=CalculateDistDiff(wholeDataMat,sampleCellSets)

[numberOfCells,numberOfFeatures]=size(wholeDataMat);
numberOfSampleSets=length(sampleCellSets);

kspMat=zeros(numberOfSampleSets,numberOfFeatures);
ksMat=zeros(numberOfSampleSets,numberOfFeatures);
deltaMeanMat=zeros(numberOfSampleSets,numberOfFeatures);
deltaMediansMat=zeros(numberOfSampleSets,numberOfFeatures);
wholeDataMeans=nanmean(wholeDataMat,1);
wholeDataMedians=nanmedian(wholeDataMat,1);

for sampleSetCounter=1:numberOfSampleSets
       deltaMeanMat(sampleSetCounter,:)=bsxfun(@minus,...
           mean(wholeDataMat(sampleCellSets{sampleSetCounter},:),1),...
           wholeDataMeans); 
 end



for featCounter=1:numberOfFeatures
    data=wholeDataMat(:,featCounter);
    %badPoints=isnan(data);
    %data=data(~badPoints);
    tRank=ceil(tiedrank(data))+1;
    tRank(isnan(tRank))=1;
    wholePDF=accumarray(tRank,1,[length(tRank)+1,1]);
    wholeCDF=cumsum(wholePDF(2:end));
    wholeCDF=wholeCDF/wholeCDF(end);
    scaleFactor=sqrt(wholeCDF.*(1-wholeCDF));
    scaleFactor(scaleFactor==0)=Inf;
    
    setMedians=zeros(numberOfSampleSets,1);
    for sampleSetCounter=1:numberOfSampleSets
        setIdx=sampleCellSets{sampleSetCounter};
        setMedians(sampleSetCounter)=nanmedian(data(setIdx));
        
        nonNanData=data(setIdx);
        nonNanData=nonNanData(~isnan(nonNanData));
        nLess=nnz(nonNanData<wholeDataMedians(featCounter))/length(nonNanData);
        nLEq=nnz(nonNanData<=wholeDataMedians(featCounter))/length(nonNanData);
        
        
        if(nLess>0.5)
            deltaMediansMat(sampleSetCounter,featCounter)=nLess-0.5;
        elseif(nLEq<0.5)
            deltaMediansMat(sampleSetCounter,featCounter)=0.5-nLEq;
        else
             deltaMediansMat(sampleSetCounter,featCounter)=0;
        end

  
        subsetPDF=accumarray(tRank(setIdx),1,[length(tRank)+1,1]);
        subsetPDF=subsetPDF(2:end);
        %nonEmpty=subsetPDF>0;
        %subsetPDF(nonEmpty)=subsetPDF(nonEmpty)./sum(subsetPDF(nonEmpty));
        subsetPDF=subsetPDF/sum(subsetPDF);
        
        subsetCDF=cumsum(subsetPDF);

        %subsetCDF=subsetCDF/subsetCDF(end);
        kspMat(sampleSetCounter,featCounter)=0.5*max(abs(subsetCDF-wholeCDF)./scaleFactor);
        ksMat(sampleSetCounter,featCounter)=max(abs(subsetCDF-wholeCDF));
%      [~,~,ksstat]=kstest2(data,data(setIdx));
%         if(abs(ksMat(sampleSetCounter,featCounter)-ksstat)>1E-10)
%            disp('oops'); 
%            1+1;
%         end
    end
%     [f,x]=ecdf(wholeDataMat(:,featCounter));
%     crazyPoints=isinf(x)|isnan(x);
%     x=x(~crazyPoints);
%     f=f(~crazyPoints);
%     %f1=interp1(x(2:end),f(2:end),sampleMedians);
%     %deltaMediansMat(:,featCounter)=(tRank(knnsearch(data,setMedians))-1)/max(tRank-1);
%     if(length(x)<=2)
%         deltaMediansMat(:,featCounter)=0;
%     else
%           deltaMediansMat(:,featCounter)=abs(interp1(x(2:end),f(2:end),setMedians)-0.5);
%     end
  
   
end
