% This script generates Figures 2d, i.e. it analyzes the A549 cell line data
% to explore how many wells are needed to capture the heterogeneity of
% diffrent cellular features,

% This script assumes the following already have been extracted from
%  the raw image data:
%  1) Single cell feature values for cells extracted from multiple wells
%   and replicated across multiple plates. This data is loaded by the LoadWellData
%  2) Pre-calculated KS prime scores for the various features, for different
%    numbers of sampled wells and replicate plates as loaded from
%    params.kspResultFile
%
addpath(genpath('Functions/'));
%% Set Parameters
kspThreshold=0.05;
cutofConfidence=0.95;


%% Load Data

%Whether KS' was defined with factor of 0.5 or not when calculation was performed
isHalfFactorInKSPrime=false;

params=GetParams('CellCulture');
kspMat=load(params.kspResultFile,'kspMat');
kspMat=kspMat.kspMat;
if(~isHalfFactorInKSPrime)
    kspMat=cellfun(@(x) 0.5*x,kspMat,'Unif',false) ;
end

[plateData, plateInfo,featureNames, featInfo]=...
    LoadWellData(params.plateList);
featInfo=featInfo.features;
%% Perform Calculations

[numberOfPlates,maxWellsSampled]=size(kspMat);
numberOfFeatures=size(kspMat{1},2);

% Calculate fraction of random samplings that yield
% value less than the KS-Prime threshold (i.e. the confidence of good
% sampling), across all features, numbers of well-sampled, randomizations
% for each, and replicated across plates.
confLevels= cellfun(@(x) sum(x<kspThreshold,1)/size(x,1) ...
    ,kspMat,'Unif',false);
confMat=zeros(numberOfPlates,maxWellsSampled,numberOfFeatures);
for featCounter=1:numberOfFeatures
    confMat(:,:,featCounter)=cellfun(@(x) x(featCounter),confLevels);
    
end

% Find minimum number of cores that ensures KS-Prime < kspThreshold with
% desired confidence. Loop over all features/plate, find highest number of
% cores where the confidence requirement is not satisfied, and pick one
% core more than that.
nrCores=zeros(numberOfPlates,numberOfFeatures);
for featCounter=1:numberOfFeatures
    for plateCounter=1:numberOfPlates
        n=find(squeeze(confMat(plateCounter,:,featCounter))...
            <cutofConfidence,1,'last');
        if(isempty(n))
            nrCores(plateCounter,featCounter)=1;
        else
            nrCores(plateCounter,featCounter)=n+1;
        end
    end
end

% We test upto maxWellsSampled (=20) cores. So, if our confidence condition
% is not satisfied for any of these, we don't truly know how many will be
% needed. Thus mean&std are unreliable. We set to Inf to be conservative
nrCores(nrCores==(maxWellsSampled+1))=Inf;


% Calculate Mean and Standard deviation of number of cores across replicate
% plates
meanNrCores=mean(nrCores);
stdNrCores=std(nrCores);
%% Select Features

% The feature generator classifies features into 6 classes (featInfo.category_id)
% 0 = ID - e.g. cellID, frameID these are not biological and not used
% 1 = intensity: measures various aspects of marker pixel intensities
% 2 = object_morphology: objects are high intensity regions in an image.
%     This measured their size etc. Because of their non-intuitiveness we
%     don't use them in our analysis
% 3 = object_moment: similarly this measures various moments of high intensity pixels
%     Because of their non-intuitiveness we don't use them in our analysis
% 4 = zernike moment- zernike texture features
% 5 = haralick texture - Haralick textures
% 6 = morphology - measuring cell morphology

% Drop ID and object features
relevantFeatures=find(~ismember(featInfo.category_id,[0,2,3]));

meanNrCores=meanNrCores(relevantFeatures);
stdNrCores=stdNrCores(relevantFeatures);
featNames=featInfo.name(relevantFeatures);
featClass=featInfo.category_id(relevantFeatures);
nrCores=nrCores(:,relevantFeatures);
featureClass=featInfo.category_id(relevantFeatures);


%% Identify Low Contrast features by based on feature names

% Set up
nrRelF=length(relevantFeatures);
channelUsed=zeros(nrRelF,3);
cellCompLabels={'D','Y','C'}; % which cellular compartment is being measure: D-Nuclear/Y-Cytoplasmic/C-WholeCell
intQuantLabels={'MultiChannel','R','I5','I20','I50','Iav','I80','Itot'}; % Different quantifiers
intQuant=zeros(nrRelF,length(intQuantLabels));
%offTargetLabels={'H2B_Y','H2B_C','XRCC5_Y','XRCC5_C','Cyto_D','Cyto_C'};
offTargetLabels={'H2B_Y','H2B_C','XRCC5_Y','XRCC5_C'}; % These represent features measuring intensity of a marker in a region it is not-expected to be. e.g. H2B in nucleus
offTarget=zeros(nrRelF,length(offTargetLabels));




for featCounter=1:nrRelF
    featIdx=relevantFeatures(featCounter);
    channelUsed(featCounter,featInfo.channel_id{featIdx})=1;
    %featureClass(featCounter,featInfo.category_id(featIdx))=1;
    
    
    
    
    name=featInfo.name{featIdx};
    temp=strsplit(name,'_');
    cInfo=temp{1};
    fInfo=temp{2};
    if(featInfo.category_id(featIdx)==1)%Break down the intensity features
        if(length(featInfo.channel_id{featIdx})==1) % Single channel features
            if(fInfo(1)=='R') %ratio feature
                featCompartment={fInfo(2),fInfo(3)};
                featQuant={'R',fInfo(4:end)};
            else
                featCompartment={fInfo(1)};
                featQuant=fInfo(2:end);
            end
            
            
        else % Across channel ratio features
            featCompartment=fInfo(2);
            featQuant={'MultiChannel','R',fInfo(3:end)};
        end
        
        intQuant(featCounter,ismember(intQuantLabels,featQuant))=1;
        
        
        % Identify off target features
        offTargetStr={};
        if(channelUsed(featCounter,1) && ismember('Y',featCompartment))
            offTargetStr= [offTargetStr {'H2B_Y'}];
        end
        if(channelUsed(featCounter,1) && ismember('C',featCompartment))
            offTargetStr= [offTargetStr {'H2B_C'}];
        end
        if(channelUsed(featCounter,2) && ismember('Y',featCompartment))
            offTargetStr= [offTargetStr {'XRCC5_Y'}];
        end
        if(channelUsed(featCounter,2) && ismember('C',featCompartment))
            offTargetStr= [offTargetStr {'XRCC5_C'}];
        end
        if(channelUsed(featCounter,2) && ismember('D',featCompartment))
            offTargetStr= [offTargetStr {'Cyto_D'}];
        end
        if(channelUsed(featCounter,2) && ismember('C',featCompartment))
            offTargetStr= [offTargetStr {'Cyto_C'}];
        end
        
        offTarget(featCounter,ismember(offTargetLabels,offTargetStr))=1;
        
        
        
    end
end

% Map all offtarget features into a new class, and renumber classes so that
% Intensity - 1
% Low Contrast Intensity - 2
% Zernike & Haralick textures both map to 3, a texture feature class
% Morphology - 4
% Object features (which are already dropped) would map to 0 if present.
temp=featureClass;
temp(any(offTarget,2))=7;
classMapping=[1,0,0,3,3,4,2];
classNamesAlt={'Intensity','Intensity(low contrast)','Texture','Morphology'};
featureClassAlt=classMapping(temp);
%% Sort features by the number of cores needed
[nrCoresSorted,featOrder]=sort(meanNrCores);
stdCoresSorted=stdNrCores(featOrder);
featClassSorted=featureClassAlt(featOrder);
featNamesSorted=featNames(featOrder);

%% Make the bar plot

classColors=[0.7 0 0; 0.7 0.4 0.4 ;0.4 0.7 0.4; 0.4 0.4 0.7]; % Colors used for the different classes

% For plotting purpose we replace the features where at least one plate
% requires >20 cores by 21 instead of Inf. This is out of the picture and
% denoted by asterisk and explained in the figure caption.
tempMean=nrCoresSorted;
tempMean(isinf(tempMean))=maxWellsSampled+1;
tempStd=stdCoresSorted;
tempStd(isinf(nrCoresSorted))=0;

figure;
for i=1:length(featClassSorted)
    fClass=featClassSorted(i);

    barColor=classColors(fClass,:);
    bar(i,tempMean(i),1,...
        'FaceColor',barColor,'EdgeColor','k');
    hold on;
end
errorbar(tempMean,tempStd,'Color','k');
xl=get(gca,'XLim');
plot(xl,[2,2],'--','Color',[0 0 0],'LineWidth',3);
hold off;
set(gca,'XLim',[0,length(featClassSorted)+1],'YLim',[0,20]);
set(gca,'XTickLabel',[],'YTickLabel',[]);
%% Export data to xls file
colNames=[{'Feature_Name','Feature_Class','Mean_Number_Cores','Stdev_Number_Cores'},...
    arrayfun(@(x) strcat('Plate',num2str(x),'_NrCores') ,1:7,'Unif',false)];
% Excel will show Inf as 65535. So we conver to NaN which show up as empty
meanCoresNan=nrCoresSorted';
meanCoresNan(isinf(meanCoresNan))=NaN;
nrCoresNan=nrCores;
nrCoresNan(isinf(nrCores))=NaN;

Fig2DData=table(featNamesSorted',classNamesAlt(featClassSorted)',meanCoresNan,stdCoresSorted',...
    nrCoresNan(1,featOrder)',nrCoresNan(2,featOrder)',nrCoresNan(3,featOrder)',...
    nrCoresNan(4,featOrder)',nrCoresNan(5,featOrder)', nrCoresNan(6,featOrder)',...
    nrCoresNan(7,featOrder)','VariableNames',colNames);
writetable(Fig2DData,'Fig2.xls','Sheet','Fig2d');