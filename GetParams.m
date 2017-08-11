function params=GetParams(paramSet,isOnCluster,varargin)

	if(nargin<2)
		isOnCluster=false;
	end

	if(isOnCluster)
		zDrive='/awlab/users/srajaram/';

	else
		zDrive='/home/z/';
	end

	params=struct;
	switch(paramSet)
		case 'Papadeno'
			params.tmaRadius=650; % comes from 0.6mm diameter core with ~ 0.4626 microns per pixel imaging.
			params.dataDir=fullfile(zDrive,'Data/Papadeno_Aperio');
			params.afiFile=fullfile(params.dataDir,'Papadeno_L34_MS1_Tumor3_2Channel.afi');
            %params.afiFile=fullfile(params.dataDir,'Papadeno_L34_MS1_Tumor3.afi');
			%params.channelNames={'DAPI','Ki67','TTF1','Bad'};
            params.channelNames={'DAPI','TTF1'};
			params.tumorOutlineFilename=fullfile(params.dataDir,'Papadeno_L34_MS1_Tumor3_TumorArea.xml');
			params.masksAndCoresFile=fullfile(params.dataDir,'cluster_nuclear1_Papadeno_L34_MS1_Tumor3_1.mat');
			params.pp_cutoff=10000;
			params.displayed_core_centers=8*[1421 705;1777 2235];
		case 'LiverCancer'
			% TMA Core Properties
			params.tmaRadius=300; %TMA radius in microns (300 corresponds to 0.6mm diameter cores)
			params.minCoreFgFraction=0.7; % Fraction of Core that must be cellular for core to be used
			
%             p=inputParser;
%             acceptableMaskTypes={'Auto','JG'};
%             addparameter(p,'MaskType','JG',...
%                 @(x) validatestring(x,acceptableMaskTypes));
%             
%             parse(p,varargin{:});
%             
%             params.maskType=p.Results.MaskType;

            params.maskType='Auto';'JG'; % Auto means only badly stained regions (identified automatically) will be dropped, whereas JG uses masks of
			% tissue regions drawn by John Gordan

			%Image Data Paths
			params.yapLkb1ImgDir=fullfile(zDrive,...
				'/Data/JG_images/AxioScan_June2016/Tiffs/');
			params.yapBcatVersaImgDir=fullfile(zDrive,...
				'/Data/JG_images/Versa_April2016/');
			params.yapBcatAxioImgDir=fullfile(zDrive,...
			'/Data/JG_images/Axio_Bcat/');
			switch(params.maskType)
				case 'JG'
					params.masksDir=fullfile(zDrive,...
					'/Data/JG_images/Masks/JG/');
					params.coreCoordsDir=fullfile(zDrive,...
					'/Data/JG_images/CoreCoords/JG/');
				case 'Auto'
					params.masksDir=fullfile(zDrive,...
					'/Data/JG_images/Masks/Auto/');
					params.coreCoordsDir=fullfile(zDrive,...
					'/Data/JG_images/CoreCoords/Auto/');

			end

			% Annotation
			params.clinicalInfoFileName=fullfile(zDrive,...
				'/Data/JG_images/160623 HCC specimen clinical info for Satwik.xlsx');
			params.pathInfoFileName=fullfile(zDrive,'/Data/JG_images/HCC specimen path comments 11-10-15.xlsx');
			params.mappingFile=fullfile(zDrive,'/Data/JG_images/Img2sample_alt.mat');	
			if(isOnCluster)
				params.sampleAnnotationFile=fullfile(zDrive,...
					'/Data/JG_images/SampleAnnotationCluster.mat');	
			else
				params.sampleAnnotationFile=fullfile(zDrive,...
					'/Data/JG_images/SampleAnnotation.mat');	
			end

			%Results Paths
			params.segResDir=fullfile(zDrive,...
				'Cluster/Tissue_Length_Scale/Nuclear_Segmentation/Results/');
			switch(params.maskType)
				case 'JG'
					params.kspResDir=fullfile(zDrive,...
						'Cluster/Tissue_Length_Scale/KSP_Calculation/Results/JG/');
				case 'Auto'
					params.kspResDir=fullfile(zDrive,...
						'Cluster/Tissue_Length_Scale/KSP_Calculation/Results/Auto/');

			end

		case 'CellCulture'
			params.kspResultFile='/home/srajaram/Work/Tissue/Code/tissue_length_scale/Cell_Culture/Well_Variation_Results.mat';
			params.plateDataDir='/home/project/2015_09_HTS_LE/2016/';
			params.plateDataDir='/home/project/2015_09_HTS_LE/plates/2016/';
			params.altPlateDataDir='/home/project/2015_09_HTS_LE/plates/2016/';
			
			params.plateList={'LE_20160202_InCell_plate_2016018019_10x_t24',...
				'LE_20160202_InCell_plate_2016018020_10x_t24',...
				'LE_20160202_InCell_plate_2016018021_10x_t24',...
				'LE_20160202_InCell_plate_2016018022_10x_t24',...
				'LE_20160202_InCell_plate_2016018023_10x_t24',...
				'LE_20160202_InCell_plate_2016018024_10x_t24',...
				'LE_20160202_InCell_plate_2016018025_10x_t24'};
		otherwise
			error('Unrecognized Parameter Set');

	end
end
