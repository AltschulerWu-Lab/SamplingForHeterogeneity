% This script generates Figures 2a, 2b and 2c in the paper based on the
% liver cancer data set.

% This script assumes pre-calculated KS prime scores for the various
% sample marker combinations, that are loaded by the LoadKSPResults
% function below.
addpath(genpath('Functions/')); % This contains all the functions used
%% Load Data and Annotations

params=GetParams('LiverCancer');

% Load Pre-calculated KS' values for all marker sample combinations, for
% 1000 random samplings each for 1 to 10 pooled cores
kspResults=LoadKSPResults();

% Find which samples are shared across the two marker sets.
sampleIDs=cellfun(@(x) {x.sampleID},kspResults,'Unif',false);
[commonSamples,commonIdx1,commonIdx2]=intersect(sampleIDs{1},sampleIDs{2});
%% Set Parameters
kspCutoff=0.2;
minConfidence=0.8;


markerSets={{'DAPI','YAP','LKB1'},{'DAPI','YAP','Beta-catenin'}};
numberOfMarkerSets=length(markerSets);

% The largest number of TMA cores that was used in sampling calculations
maxNumberOfCores=10;

%Whether KS' was defined with factor of 0.5 or not when calculation was performed
isHalfFactorInKSPrime=false; 

% IDs of adjacent normal samples (i.e. non-tumor) samples. 
nonTumorSamples={'773'}; %

% Plotting
lim=[0,maxNumberOfCores+0.5]; % Range of number of cores displayed
markerColors=0.7*[0 0 1; 0 1 0; 1 0 0; 0 0 1; 0 1 0; 1 0 1]; % Colors for markers in markerSets
%% Calculations


% Convenience function to calculate fraction of elements in row nC of 
% matrix x, that are less than a cutoff scoreThreshold. This will be used
% in the calculation of confidence of being within a KS-Prime threshold
% with a given number of cores (which will be stored in confCurves)
confCal =@(x,nC,scoreThreshold) nnz(x(nC,:)<scoreThreshold)/numel(x(nC,:));
confCurves=cell(numberOfMarkerSets,3);

% Will store the min number of cores needed to achieve KS-Prime < kspCutoff
% with confidence greater than minConfidence.
minCores=cell(numberOfMarkerSets,1); 

for markerSetCounter=1:numberOfMarkerSets % Loop over marker sets
    
    % Initialize Variables
    numberOfSamplesInMarkerSet=length(kspResults{markerSetCounter});
    numberOfMarkersInMarkerSet=length(markerSets{markerSetCounter});
    minCores{markerSetCounter}=zeros(numberOfSamplesInMarkerSet,numberOfMarkersInMarkerSet);
    confCurves(markerSetCounter,:)={zeros(numberOfSamplesInMarkerSet,maxNumberOfCores)};
    
    
    for sampleCounter=1: numberOfSamplesInMarkerSet %Loop over samples in marker set
        
        % Load results for specific sample
        sampleResults=kspResults{markerSetCounter}(sampleCounter);
        
        
        for markerCounter=1:numberOfMarkersInMarkerSet %Loop over markers in marker set
            
            % Pull out results for marker in order specified in markerSets
            markerName= markerSets{markerSetCounter}{markerCounter};
            markerIdx=find(strcmpi(sampleResults.markerNames,markerName));
            stats=sampleResults.kspStats{markerIdx};
            if(isHalfFactorInKSPrime)
                ksPrimes=stats.coreKsPrimes; 
            else
                %Multiply by 0.5 if KS prime calculation was performed
                %without this factor.
                ksPrimes=0.5*stats.coreKsPrimes; 
            end
            
            % Calculate minimum number of cores needed for sample-marker
            % combo (in marker set) to achieve KS-Prime < kspCutoff
            % with confidence greater than minConfidence.
            minCores{markerSetCounter}(sampleCounter,markerCounter)=...
                MinSamplesNeeded(ksPrimes,kspCutoff,minConfidence);
            
            
            % Calculate for sample-marker combo (in marker set) 
            % how likely we are to achieve KS-Prime < kspCutoff
            % with different number of cores
            confCurves{markerSetCounter,markerCounter}(sampleCounter,:)=...
                arrayfun(@(n) confCal(ksPrimes,n,kspCutoff) ,1:maxNumberOfCores);
            
        end
        
        
    end
end

%% Fig 2a: Comparing Number of Cores for YAP across marker sets

markerNum=2; % For YAP


% Load number of cores results
x=minCores{1}(commonIdx1,markerNum);
y=minCores{2}(commonIdx2,markerNum);

% Drop non-tumor samples
samplesToUse=~ismember(commonSamples,nonTumorSamples);
x=x(samplesToUse);
y=y(samplesToUse);

% Count how many times each pair of number of cores appears (this will
% become the point size)
[a,b,c]=unique([x,y],'rows');
if(length(a)<length(x))
    temp=accumarray(c,50*ones(length(c),1));
    s=temp(c);
else
    s=50;
end

% Make Plot
figure;
scatter(x,y,s,'k','filled');
hold on;
plot(lim,lim,'--k','LineWidth',2);
hold off;
set(gca,'XLim',lim,'YLim',lim);
set(gca,'XTick',[1,5,10],'YTick',[1,5,10],'XTickLabel',[],'YTickLabel',[]);
axis square;

comment=cell(size(x));comment(:)={''};
comment(isnan(x)|isnan(y))={'Empty value means 10 cores, the max tested, was not enough to gurantee KS-Prime was within bounds at desired level of confidence. This point was not included in the plot'};

Fig2aData=table([1:length(x)]',double(x),double(y),s/50,comment,'VariableNames',...
    {'SampleNumber','Cores_YAP_MS1','Cores_YAP_MS2','Point_Size_Equals_Num_Samples_With_NumCores_Pair','Comment'});
writetable(Fig2aData,'Fig2.xls','Sheet','Fig2a');



%% Comparing number of cores needed by DAPI and YAP in marker set 1

% Load number of cores results
markerSetNumber=1;
x=minCores{markerSetNumber}(:,strcmp(markerSets{markerSetNumber},'YAP'));
y=minCores{markerSetNumber}(:,strcmp(markerSets{markerSetNumber},'DAPI'));


% Drop non-tumor samples
samplesToUse=~ismember(sampleIDs{markerSetNumber},nonTumorSamples);
x=x(samplesToUse);
y=y(samplesToUse);


% Count how many times each pair of number of cores appears (this will
% become the point size)
[a,b,c]=unique([x,y],'rows');
if(length(a)<length(x))
    temp=accumarray(c,50*ones(length(c),1));
    s=temp(c);
else
    s=50;
end

% Make plot
figure;
scatter(x,y,s,'k','filled');

hold on;
plot(lim,lim,'--k','LineWidth',2); % Diagonal line
hold off;

set(gca,'XLim',lim,'YLim',lim);

set(gca,'XTick',[1,5,10],'YTick',[1,5,10],'XTickLabel',[],'YTickLabel',[]);
axis square;


comment=cell(size(x));comment(:)={''};
comment(isnan(x)|isnan(y))={'Empty value means 10 cores, the max tested, was not enough to gurantee KS-Prime was within bounds at desired level of confidence. This point was not included in the plot'};

Fig2bData=table([1:length(x)]',double(x),double(y),s/50,comment,'VariableNames',...
    {'Sample_Number','Cores_YAP_MS1','Cores_DAPI_MS1','Point_Size_Equals_Num_Samples_With_NumCores_Pair','Comment'});
writetable(Fig2bData,'Fig2.xls','Sheet','Fig2b');


%% Fig 2c: Confidence curves for all markers

figure;
counter=1; % This counter is used to map markers onto their colors
confidenceCurve=cell(6,1); %6 markers used;
for markerSetCounter=1:numberOfMarkerSets % Loop over marker sets
    % Only preserve tumor samples
    samplesToUse=~ismember(sampleIDs{markerSetCounter},nonTumorSamples);
    
    for markerCounter=1:length(markerSets{markerSetCounter}) % Loop over markers in marker set
        
        % Average confidence across all samples (that we want to use)
        % Note, the 
        confidenceCurve{counter}=100*mean(confCurves{markerSetCounter,markerCounter}(samplesToUse,:))';
        
        % Plot Confidence Curve 
        plot(confidenceCurve{counter},'Color',markerColors(counter,:),'LineWidth',3);
        
        hold on;
        counter=counter+1;
    end
end

hold off;
set(gca,'XLim',[1,10],'YLim',[30,100]);
set(gca,'XTickLabel',[],'YTickLabel',[]);

columnNames={'PctPatients_DAPI_MS1','PctPatients_YAP_MS1','PctPatients_LKB1_MS1',...
    'PctPatients_DAPI_MS2','PctPatients_YAP_MS2','PctPatients_Bcat_MS2'};
Fig2cData=table([1:maxNumberOfCores]',confidenceCurve{:},'VariableNames',...
   ['Number_Of_Cores',columnNames]);
writetable(Fig2cData,'Fig2.xls','Sheet','Fig2c');
