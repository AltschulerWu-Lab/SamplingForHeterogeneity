classdef KSP_Calculator
   properties
       wholeDataRanks;
       wholeCDF;
       scaleFactor;
   end
   
   methods
       function obj=KSP_Calculator(wholeData)
           obj.wholeDataRanks=ceil(tiedrank(wholeData))+1;
           obj.wholeDataRanks(isnan(obj.wholeDataRanks))=1;
           wholePDF=accumarray(obj.wholeDataRanks,1,[length(obj.wholeDataRanks)+1,1]);
           obj.wholeCDF=cumsum(wholePDF(2:end));
           obj.wholeCDF=obj.wholeCDF/obj.wholeCDF(end);
           obj.scaleFactor=sqrt(obj.wholeCDF.*(1-obj.wholeCDF));
           obj.scaleFactor(obj.scaleFactor==0)=Inf;
       end
       
       function [ksp,ks]=Calculate(obj,subset_idx)
           
           subsetPDF=accumarray(obj.wholeDataRanks(subset_idx),1,...
               [length(obj.wholeDataRanks)+1,1]);
           subsetPDF=subsetPDF(2:end);
           subsetPDF=subsetPDF/sum(subsetPDF);

           subsetCDF=cumsum(subsetPDF);
          
           ksp=0.5*max(abs(subsetCDF-obj.wholeCDF)./obj.scaleFactor);
           ks=max(abs(subsetCDF-obj.wholeCDF));
       end
       
   end
    
    
    
end