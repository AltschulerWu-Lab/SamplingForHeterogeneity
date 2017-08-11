function [wellData,cellInfo,featureNames,featRes]=LoadWellData(plateList,varargin)

	p=inputParser;

	defaultRowList=[2:15,2:15]';
	addParameter(p,'RowList',defaultRowList,...
		@(x) validateattributes(x,{'numeric'},...
		{'positive','integer','vector','<=',16})); 

	
	defaultColList=repmat([2,23],14,1);
	defaultColList=defaultColList(:);
	addParameter(p,'ColList',defaultColList,...
		@(x) validateattributes(x,{'numeric'},...
		{'positive','integer','vector','<=',24})); 

	addParameter(p,'IsOnCluster',false,...
		@(x) validateattributes(x,{'logical'},{'scalar'}));
	
	parse(p,varargin{:});

	rowList=p.Results.RowList;
	colList=p.Results.ColList;
	if(length(rowList)~=length(colList))
		error('Row and Column list must have the same lengths');
	end
	
	params=GetParams('CellCulture');
	if(p.Results.IsOnCluster)
		projectDir=params.plateDataDir;
	else
		projectDir=params.altPlateDataDir;
	end
	
	vertcatcell=@(cellArray) vertcat(cellArray{:});

	numberOfWellsPerPlate=length(rowList);
	numberOfPlates=length(plateList);
	numberOfWells=numberOfPlates*numberOfWellsPerPlate;

	featData=cell(numberOfWells,1);
	cellPos=cell(numberOfWells,1);
	for plateCounter=1:numberOfPlates
		for wellCounter=1:numberOfWellsPerPlate
			
			% Define Variables
			wn=wellCounter+(plateCounter-1)*numberOfWellsPerPlate;
			rowNumber=rowList(wellCounter);
			colNumber=colList(wellCounter);
			
			% Open Data Files
			plateDir=fullfile(projectDir,...
				plateList{plateCounter});
			featuresFile=fullfile(plateDir,'/features/cbfeatures/',...
				['cbfeatures-' num2str(rowNumber) '-'...
    			num2str(colNumber) '.mat']);
			segmentationFile=fullfile(plateDir,'/segmentation/',...
				['segmentation-' num2str(rowNumber) '-'...
    			num2str(colNumber) '.mat']);
        	featRes=load(featuresFile);
        	segRes=load(segmentationFile);
        
			% Extract Feature and Position Data
        	cellPos{wn}=vertcatcell(cellfun(@mean,...
            	segRes.segmentation.dnabdry,'Unif',false));
        	featData{wn}=featRes.features.data;


		end

	end
	featureNames=featRes.features.name;
	wellData=vertcat(featData{:});
	cellPos=vertcat(cellPos{:});

	AnnotateWells=@(dataArray,annotArray) vertcatcell(arrayfun(@(x) ...
    ones(size(dataArray{x},1),1)*annotArray(x), 1:size(dataArray,1),...
    'Unif',false));

	cellInfo=struct;
	cellInfo.wellNumbers=AnnotateWells(featData,1:length(featData));
	cellInfo.rowNum=AnnotateWells(featData,repmat(rowList,numberOfPlates,1));
	cellInfo.colNum=AnnotateWells(featData,repmat(colList,numberOfPlates,1));
	cellInfo.yPos=cellPos(:,1);
	cellInfo.xPos=cellPos(:,2);
	plateNumbers=repmat(1:numberOfPlates,numberOfWellsPerPlate,1);
	cellInfo.plate=AnnotateWells(featData,plateNumbers(:));


end
