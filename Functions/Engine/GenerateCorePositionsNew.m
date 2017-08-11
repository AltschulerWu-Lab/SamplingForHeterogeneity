function corePositions =GenerateCorePositionsNew(tImgList,fgMaskList,...
    maxNumberSamples, numberOfReplicates,params)

numberOfImages=length(tImgList);


if(length(fgMaskList) ~= numberOfImages)
    error('TissueImage and Mask List length inconsistent');
end


% Make sure fgmasks are at same maglevel/mpp
mppList=zeros(numberOfImages,1);
for imgCounter=1:numberOfImages
    dims=tImgList{imgCounter}.dimensions;
    fgLayer=find(dims(:,1)==size(fgMaskList{imgCounter},1) ...
        & dims(:,2)==size(fgMaskList{imgCounter},2));
    if(isempty(fgLayer))
        error('Cannot Match FG Maks and Image Layers');
    end
    mppList(imgCounter)=tImgList{imgCounter}.micronsPerPixel*dims(1,1)/dims(fgLayer,1);
    
end

micronsPerPixel=max(mppList); %Use the highest mpp (i.e. most zoomed out)
wrongMags=find(mppList/micronsPerPixel<0.95); % Rescale masks that are not at correct mpp
if(~isempty(wrongMags))
    warning('Masks are at different magnifications, attempting to correct');
end
for layerCounter=1:length(wrongMags)
    layerNum=wrongMags(layerCounter);
    fgMaskList{layerNum}=imresize(fgMaskList{layerNum},mppList(layerNum)/micronsPerPixel);
end


% Generate valid sampling positions and sampling blocks, organize as
% layers for different images
tmaRadius=params.tmaRadius/micronsPerPixel; %Convert to pixel units
blockWidth=round(tmaRadius/sqrt(2)); % This way cores not in neighboring blocks are at least 2 tmaRadii apart

fgPixelIdxByBlock=cell(numberOfImages,1);
validBlocks=cell(numberOfImages,1);

for imgCounter=1:numberOfImages
    
    % Find pixels that have more than the min frac of fg pixels within
    % TMA radius. These are potential TMA core positions. The two
    % reszing operations are to speed up calculations
    
    rescaleFactor=10/micronsPerPixel;
    fgMask=imresize(fgMaskList{imgCounter},1/rescaleFactor);
    isValidCorePos=imfilter(double(fgMask),...
        fspecial('disk',tmaRadius/rescaleFactor))>params.minCoreFgFraction;
    isValidCorePos=imresize(isValidCorePos,size(fgMaskList{imgCounter}));
    
    % Organize valid positions by the image block in which they are
    % contained
    fgPixelIdxByBlock{imgCounter}=GetTruePixelsByBlock(isValidCorePos,blockWidth);
    
    % Valid blocks are those that contain at least one potential core
    % position. For each valid block, image from which it came is also
    % stored.
    %[vBy,vBx]=find(~cellfun(@isempty,fgPixelIdxByBlock{imgCounter}));

    %validBlocks{imgCounter}=[vBy,vBx,imgCounter*ones(length(vBy),1)];
    %validBlocks{imgCounter}=sparse(vBy,vBx,ones(length(vBy),1),...
    %    size(fgPixelIdxByBlock{imgCounter},1),size(fgPixelIdxByBlock{imgCounter},2));
    validBlocks{imgCounter}=sparse(cellfun(@(x) size(x,1),fgPixelIdxByBlock{imgCounter}));
    
end
% validBlocks=vertcat(validBlocks{:}); % Combine blocks from all images

corePositions=cell(maxNumberSamples,numberOfReplicates);
for numberOfSamples=1:maxNumberSamples
    for repCounter=1:numberOfReplicates
        blocksAvailable=validBlocks;
        
       
        corePos=zeros(numberOfSamples,3);
 
        sampleCounter=1;    
        numberFailures=0;
        while sampleCounter<=numberOfSamples
        
            [blockList,pointsPerBox]=MakeBlockList(blocksAvailable);
            
            isValid=false;
            while(size(blockList,1)>0 & ~isValid)
                
              
                %chosenBox=randi(size(blockList,1));
                chosenBox=find(mnrnd(1,pointsPerBox/sum(pointsPerBox)));
                x=blockList(chosenBox,2);
                y=blockList(chosenBox,1);
                imgNum=blockList(chosenBox,3);
                
                possiblePositions=fgPixelIdxByBlock{imgNum}{y,x};
                if(sampleCounter==1)
                 
                    chosenPosition=possiblePositions(randi(...
                        size(possiblePositions,1)),:);
                
                    isValid=true;
                
                else
                    %Check if there are any valid positions far enough from
                    %previous points
                    % If yes, pick one of them and set valid
                    
                    % If not, remove this block from available list and keep invalid
                    prevPosInImg= corePos(1:(sampleCounter-1),:);
                    prevPosInImg=prevPosInImg(prevPosInImg(:,3)==imgNum,1:2);
                    
                    
                    
                    if(~isempty(prevPosInImg))
                        posIdx=find((min(pdist2(possiblePositions,prevPosInImg),...
                            [],2)>2*tmaRadius));
                    else
                        posIdx=1:size(possiblePositions);
                    end
                    
                    
                    
                    if(~isempty(posIdx))
                        chosenPosIdx=posIdx(randi(length(posIdx)));
                        chosenPosition=possiblePositions(chosenPosIdx,:);
                        isValid=true;
                    else
                        blocksAvailable{imgNum}(y,x)=false;
                        [blockList,pointsPerBox]=MakeBlockList(blocksAvailable);
                        
                    end
                    %                     if(min(pdi possiblePositions=fgPixelIdxByBlock{imgNum}{y,x};st(pos,prevPosInImg))<tmaRadius)
                end
        
                
           
                
                
                
            end
            
            if(~isempty(blockList))
                corePos(sampleCounter,:)=[chosenPosition,imgNum];
                
                neighborsX=x+(-1:1);
                neighborsX(neighborsX<1)=[];
                neighborsX(neighborsX>size(blocksAvailable{imgNum},2))=[];
                neighborsY=y+(-1:1);
                neighborsY(neighborsY<1)=[];
                neighborsY(neighborsY>size(blocksAvailable{imgNum},1))=[];
                
                blocksAvailable{imgNum}(neighborsY,neighborsX)=false;
                sampleCounter=sampleCounter+1;
            else
                corePos=zeros(numberOfSamples,3);
                sampleCounter=1;
                numberFailures=numberFailures+1;
                if(numberFailures>2000)
                   error('MyComponent:TooManyFailures',...
                       'I give up, cannot fit %d samples in',numberOfSamples);
                end
            end
        end
        % TODO: rescale core pos to max mag and add random stagger 
        magFactor=round(micronsPerPixel/tImgList{1}.micronsPerPixel); % all images supposed to be at same mpp, so can use first
        randStagger=randi(magFactor,size(corePos,1),2)-round(magFactor/2);
        corePos(:,1:2)=corePos(:,1:2)*magFactor+randStagger;
        corePositions{numberOfSamples,repCounter}=corePos;
     
    end
end


%sampling loop over number of samples replicates etc

end

function [blockList,pointsPerBlock]=MakeBlockList(blocksAvailable)
    
    blockList=cell(length(blocksAvailable),1);
    pointsPerBlock=cell(length(blocksAvailable),1);
    for i=1:length(blocksAvailable)
        idx=find(blocksAvailable{i});
%         [by,bx]=find(blocksAvailable{i});
        [by,bx]=ind2sub(size(blocksAvailable{i}),idx);
        blockList{i}=[by,bx,i*ones(length(by),1)];
        pointsPerBlock{i}=full(blocksAvailable{i}(idx));
    end
    blockList=vertcat(blockList{:});
    pointsPerBlock=vertcat(pointsPerBlock{:});
end

function pixIdx=GetTruePixelsByBlock(boolMat,blockWidth)


ny=ceil(size(boolMat,1)/blockWidth);
nx=ceil(size(boolMat,2)/blockWidth);

pixIdx=cell(ny,nx);
for i=1:ny
    yRange=(1:blockWidth)+(i-1)*blockWidth;
    yRange=yRange(yRange<=size(boolMat,1));
    for j=1:nx
        xRange=(1:blockWidth)+(j-1)*blockWidth;
        xRange=xRange(xRange<=size(boolMat,2));
        [yIdx,xIdx]=find(boolMat(yRange,xRange));
        pixIdx{i,j}=[yIdx+yRange(1)-1,xIdx+xRange(1)-1];
        
    end
end

end