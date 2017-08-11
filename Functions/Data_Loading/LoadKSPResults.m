function kspResults=LoadKSPResults(isOnCluster,params)

    numberOfMarkerSets=2;
    if(nargin<2)
		if(nargin==1)
	        params=GetParams('LiverCancer',isOnCluster);
		else
	        params=GetParams('LiverCancer');
		end
    end
    
    kspResults=cell(numberOfMarkerSets,1);
    sampleInfo=load(params.sampleAnnotationFile);
    sampleInfo=sampleInfo.anno;
    
    for markerSetCounter=1:numberOfMarkerSets
        numberOfSamples=length(sampleInfo{markerSetCounter});
        kspResults{markerSetCounter}=struct;
        dataCounter=1;
        for sampleCounter=1:numberOfSamples
           sampleID=sampleInfo{markerSetCounter}(sampleCounter).sampleID;
           kspResultsFileName=fullfile(params.kspResDir,...
               ['MS' num2str(markerSetCounter)],[sampleID '.mat']);
           if(exist(kspResultsFileName,'file'))
               results=load(kspResultsFileName);
               annoFields=fieldnames(results.anno);
               for fieldCounter=1:length(annoFields)
                  field=annoFields{fieldCounter};
                  kspResults{markerSetCounter}(dataCounter).(field)=...
                      results.anno.(field);
                   
               end
               kspResults{markerSetCounter}(dataCounter).kspStats=results.kspStats;
               dataCounter=dataCounter+1;
           end
        end
        
    end
    


end
