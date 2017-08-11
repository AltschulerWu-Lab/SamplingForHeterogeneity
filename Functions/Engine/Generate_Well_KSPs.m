% Calculations for Fig2d: Cell Culture
% This script calculates the KS-Prime values that are then analyzed and
% plotted in Fig 2d.

% Make sure that root diectory and Functions/Data_Loading are in path
% addpath('../../');
% addpath('../Data_Loading');
params=GetParams('CellCulture');
plateList=params.plateList;

%%
nrWellsRange=1:20;
numberOfPlates=length(plateList);
plateData=cell(numberOfPlates,1);
plateInfo=cell(numberOfPlates,1);

kspMat=cell(numberOfPlates,length(nrWellsRange));
ksMat=cell(numberOfPlates,length(nrWellsRange));
deltaMeans=cell(numberOfPlates,length(nrWellsRange));
deltaMedians=cell(numberOfPlates,length(nrWellsRange));
plateCellSubsets=cell(numberOfPlates,length(nrWellsRange));

for plateCounter=1:numberOfPlates
    
    [plateData{plateCounter}, plateInfo{plateCounter},featureNames]=...
        LoadWellData(plateList(plateCounter),'IsOnCluster',false);

    
    
    wellsInPlate=unique(plateInfo{plateCounter}.wellNumbers);
    numberOfWellsInPlate=length(wellsInPlate);
    
    for numberOfWellsSubsampled=nrWellsRange
        tic;
        if(numberOfWellsSubsampled<3)
            
            wellsToUse=nchoosek(wellsInPlate,numberOfWellsSubsampled);
            
        else
            wellsToUse=zeros(1000,numberOfWellsSubsampled);
            for repNumber=1:1000
                wellsToUse(repNumber,:)=randsample(wellsInPlate,numberOfWellsSubsampled);
            end
        end
        
        numberOfReps=size(wellsToUse,1);
 
        cellIdx=cell(numberOfReps,1);

        for repCounter=1:numberOfReps
            cellIdx{repCounter}=find(ismember(plateInfo{plateCounter}.wellNumbers,...
                wellsToUse(repCounter,:)));

        end
         [kspMat{plateCounter,numberOfWellsSubsampled},...
            ksMat{plateCounter,numberOfWellsSubsampled},...
            deltaMeans{plateCounter,numberOfWellsSubsampled},...
            deltaMedians{plateCounter,numberOfWellsSubsampled}]=...
            CalculateDistDiff(plateData{plateCounter},...
           cellIdx);
        disp(['Plate ' num2str(plateCounter) ' NrWells:' ...
            num2str(numberOfWellsSubsampled)  'Done']);
        toc;
    end
    
    disp(['Plate ' num2str(plateCounter) ' Loaded!']);
    
end

% This workspace is what is loaded up in Fig2d.m 
