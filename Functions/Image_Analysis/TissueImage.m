classdef TissueImage
   properties
       numberOfChannels;
       numberOfMagLevels;
       dimensions;
       imgFileMat;
       imgLayerMat;
       micronsPerPixel;
       channelNames;
       bitDepth;
   end
   
   methods
       
       function obj=TissueImage(imgFileMat,imgLayerMat,varargin)
           if(isequal(size(imgFileMat),size(imgLayerMat)))
            obj.imgFileMat=imgFileMat;
            obj.imgLayerMat=imgLayerMat;
           else
               error('File and Layer Matrices must have matching sizes');
           end
           
           [obj.numberOfMagLevels,obj.numberOfChannels]=size(imgFileMat);
           
            p=inputParser;
            
            
            addParameter(p,'MicronsPerPixel',[],@(x) validateattributes(x,{'numeric'},...
                {'scalar','positive'}));
            defaultChannelNames=strcat('Channel',...
                arrayfun(@num2str,1:obj.numberOfChannels,'Unif',false));
            addParameter(p,'ChannelNames',defaultChannelNames,...
                @(x) iscell(x) & length(x)==obj.numberOfChannels& all(cellfun(@ischar,x)));
                  
            
           parse(p,varargin{:});
           obj.micronsPerPixel=p.Results.MicronsPerPixel;
           obj.channelNames=p.Results.ChannelNames;
           
           % Extract dimensions of the different mag levels
           dims=zeros(obj.numberOfMagLevels,2*obj.numberOfChannels);
           for magLevel=1:obj.numberOfMagLevels
              for channelCounter=1:obj.numberOfChannels
                 info=imfinfo(imgFileMat{magLevel,channelCounter});
                 info=info(imgLayerMat(magLevel,channelCounter));
                 dims(magLevel,2*(channelCounter-1)+1)=info.Height;
                 dims(magLevel,2*channelCounter)=info.Width;
              end
           end
           obj.bitDepth=info.BitDepth;
           
           if(all(std(dims(:,1:2:2*obj.numberOfChannels),[],2)==0)&&...
                   all(std(dims(:,2:2:2*obj.numberOfChannels),[],2)==0)) %Check if the different channels yield images of the same size
               obj.dimensions=dims(:,1:2);
               
           else
                error('The different channels dont have same image sizes');
           end
           if(~(isequal(obj.dimensions,sort(obj.dimensions,'descend'))))
               error('Mag levels need to be in decreasing order of size');
           end
           
           
       end
       
       function img=LoadImage(obj,varargin)
            p=inputParser;
            
            
            addParameter(p,'MagLevel',1,@(x) validateattributes(x,{'numeric'},...
                {'scalar','integer','positive','<=',obj.numberOfMagLevels}));
            addParameter(p,'Channels',1:obj.numberOfChannels,@(x) ...
                validateattributes(x,{'numeric'},{'vector','integer',...
                'positive','<=',obj.numberOfChannels}));
            addParameter(p,'PixelRegion',{}, @(x) validateattributes(x,...
                {'cell'},{'numel',2}));
            
            
            parse(p,varargin{:});
            
            %add checker function for PixelRegion
            if(isempty(p.Results.PixelRegion))
                img=zeros([obj.dimensions(p.Results.MagLevel,:),length(p.Results.Channels)]);
            else
                img=zeros([cellfun(@diff,p.Results.PixelRegion)+1,length(p.Results.Channels)]);
            end
            for channelCounter=1:length(p.Results.Channels)
                
                channelNumber=p.Results.Channels(channelCounter);
                img(:,:,channelCounter)=imread(obj.imgFileMat{p.Results.MagLevel,channelNumber},...
                    obj.imgLayerMat(p.Results.MagLevel,channelNumber),'PixelRegion',p.Results.PixelRegion);
            end
            
            
       end
       
   end
   
   
end
