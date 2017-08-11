% This script generates the Fig 1 in the paper

% These results are generated from the analysis of a single Papilarray
% adenocarcinoma image stained for DAPI and TTF1. 
% Data needed to run this script (locations of data files specified in GetParams.m
% under case 'Papadeno', and will need to be changed to reflect paths in your system) 
% 1) Image file(s) for tissue: Specified in GetParams.m as params.afiFile.
% Note: you will need to ensure tif files pointed to by the afi, are in the
% same directory as the afi file.
% 2) Tumor Outline: Hand drawn annotation of tumor boundary -green line in
% figure 1b. This is stored as an xml file pointed to by
% params.tumorOutlineFilename
% 3) Core positions/Mask: Mask identifying forgeround & locations of randomly
% placed cores were previously generated. For this tissue sample they are 
%stored in a mat file, pointed to by params.masksAndCoresFile

addpath(genpath('Functions/')); % This contains all the functions used


%% Parameter Setup

params=GetParams('Papadeno');

% Segmentation (identification of cells in the image)
% Because the tissue image is so large, it is broken into smaller images,
% to perform analysis. Here we spcify size of these images.
block_xres=2000; 
block_yres=2000;
min_nuclear_area=100; % Minimum size of a nucleus in pixels. 

% Main Analysis
core_radius=params.tmaRadius; % The radius of a TMA core. 
max_number_of_cores=10; % The max number of pooled cores to analyze in generating samples
number_of_reps=1000; % Number of "replicate" samplings for a given number of cores
number_of_ksp_bins=100; % Number of bins used in calculating the CDFs in Fig 1d.


% Plotting
hist_bins=2.^linspace(10,16,64); % intensity binning. Used for plotting not calculations
bin_widths=diff([hist_bins,hist_bins(end)])+1;
nr_cores_to_show=[1,3,6,10]; % values for number of cores used in Fig. 1c/d
kspRange=[0,0.4];
wholeColor=[0.4416,0.7490,0.4322]; % fyi, same as linspecer 3
kspColors=TwoColorSpanningCmap([0.9 0.9 0.1],[0 0 1],128);
nrCores_colors=(0:(max_number_of_cores-1))'*[(1/(max_number_of_cores-1)) 0 0];
randomSamplingColor=[0.7 0.7 0.7];
% For Fig 1.b 
chosen_core_centers=params.displayed_core_centers;
display_scaling_range=[0,15000];% Min/max intensities for display of tissue image
numberOfCores=size(chosen_core_centers,1); 
border_colors=linspecer(numberOfCores+1);
% Fig 1c
nrRands=50; % Number of curves to plot
%% Load Data

% Images
tImg=Afi2TImg(params.afiFile);
dapi_img=tImg.LoadImage('Channels',find(strcmp(params.channelNames,'DAPI')));
img=tImg.LoadImage('Channels',find(strcmp(params.channelNames,'TTF1')));
tumor_outline=Load_Tumor_Outline(params.tumorOutlineFilename);
img_lowres=tImg.LoadImage('MagLevel',tImg.numberOfMagLevels,...
    'Channels',find(strcmp(params.channelNames,'TTF1')));

%mask
mask_and_cores=load(params.masksAndCoresFile,'pixel_mask_to_use','raw_tma_positions');
pixel_mask=mask_and_cores.pixel_mask_to_use;

% Core coordinates

%% Segmentation Calculation
%Perform nuclear segmentation and extract mean intensity of each nucleus
[yres,xres]=size(img);
block_nx=ceil(xres/block_xres);
block_ny=ceil(yres/block_yres);
marker_vals=cell(block_ny,block_nx);
nuclear_centroids=cell(block_ny,block_nx);

% Note the image is broken into blocks for performing segmentation since 
% analysis on full tissue image is too memory intensive (results
% ultimately combined across blocks). Loop is over image blocks.
for nx_counter=1:block_nx
    for ny_counter=1:block_ny
        tic;
        
        % Load up data for block of image being analyzed
        block_yrange=(1:(block_yres)) + (ny_counter-1)*block_yres;
        block_xrange=(1:(block_xres)) + (nx_counter-1)*block_xres;
        block_yrange(block_yrange>yres)=[];
        block_xrange(block_xrange>xres)=[];
        dna_img=dapi_img(block_yrange,block_xrange);
        marker_img=img(block_yrange,block_xrange);
        mask_fg=pixel_mask(block_yrange,block_xrange);
        
        % Calculate Laplacian of Gaussian response on dna image. Higher
        % response expected at center of nuclei.
        h_log=fspecial('log',30,10);
        img_log=imfilter(double(dna_img),h_log);
        img_log_fg=img_log;
        img_log_fg(~mask_fg)=max(img_log_fg(:));
        
        % Find good minima to act as seeds for watershed
        img_min=imextendedmin(img_log_fg,0.1);
        img_log_fg=imimposemin(img_log_fg,img_min);
        
        % Perform watershed        
        L=watershed(img_log_fg);
        
        % Clean up
        border_mask=imclearborder(L>0);
        clean_labels=L;
        clean_labels(~border_mask)=0;
        clean_labels(~mask_fg)=0;
        
        props=regionprops(clean_labels,{'Centroid','Area','PixelIdxList'});
        
        if(~isempty(props))
         
            is_too_small=[props(:).Area]<min_nuclear_area;
            clean_labels(vertcat(props(is_too_small).PixelIdxList))=0;
            props(is_too_small)=[];
            if(~isempty(props))
           
                % Value of TTF1 intensity is calculated on each nucleus
                marker_vals{ny_counter,nx_counter}=cellfun(@(x) mean(marker_img(x)),{props(:).PixelIdxList});
                
                % Position of cells are stored,
                nuclear_centroids{ny_counter,nx_counter}=bsxfun(@plus,vertcat(props(:).Centroid),...
                    [block_xrange(1),block_yrange(1)]-1);
                
            end
            
        end
        toc;
        disp([num2str(nx_counter) ' x ' num2str(ny_counter) ' done!']);
    end
end

% Combine results across all blocks
centroids_mat=vertcat(nuclear_centroids{:});
marker_vals_mat=horzcat(marker_vals{:})';
%% KS-Prime Calculations

core_nuclear_vals=cell(max_number_of_cores,number_of_reps);
core_nuclear_hists=cell(max_number_of_cores,number_of_reps);

whole_tissue_hist=histc(marker_vals_mat,hist_bins);
whole_tissue_pdf=whole_tissue_hist/sum(whole_tissue_hist)./bin_widths';
% whole_tissue_cdf=cumsum(whole_tissue_pdf).*bin_widths';
% whole_tissue_cdf(whole_tissue_cdf>1)=1;


kspEngine=KSP_Calculator(marker_vals_mat); %initialize the KSP calculator

core_ks_primes=zeros(max_number_of_cores,number_of_reps);
core_ks=zeros(max_number_of_cores,number_of_reps);
for numberOfCores=1:max_number_of_cores % Loop over number of cores
    tic
    for rep_number=1:number_of_reps % loop over replicate samplings
        
        % Find cells within any of the sampled cores
        core_coords=mask_and_cores.raw_tma_positions{numberOfCores,rep_number};
        core_coords=fliplr(core_coords);
        dists_to_cores=pdist2(core_coords,centroids_mat);
        [tma_idx,cell_idx]=find(dists_to_cores<core_radius);
        
        % Get core distributions (for plots)
        core_nuclear_vals{numberOfCores,rep_number}=marker_vals_mat(cell_idx);
        core_nuclear_hists{numberOfCores,rep_number}=histc(core_nuclear_vals{numberOfCores,rep_number},hist_bins);
        
        % Calculate KS prime
        [core_ks_primes(numberOfCores,rep_number),...
            core_ks(numberOfCores,rep_number)]=kspEngine.Calculate(cell_idx);
    end
    toc;
end


%whole_tissue_mean=mean(marker_vals_mat);
whole_tissue_median=median(marker_vals_mat);
%core_means=cellfun(@mean,core_nuclear_vals);
core_delta_medians=abs(cellfun(@(x) nnz(x>whole_tissue_median)/numel(x), core_nuclear_vals)-0.5);

% Cells selected in a spatially random manner (rather than using randomly placed cores)
random_ksp=zeros(number_of_reps,1);
random_ks=zeros(number_of_reps,1);
random_delta_median=zeros(number_of_reps,1);
%random_delta_mean=zeros(number_of_reps,1);
for rCounter=1:number_of_reps
    
    nPoints=length(core_nuclear_vals{1,randi(1000)});
    cell_idx=randsample(length(marker_vals_mat),nPoints);
    
    [random_ksp(rCounter),random_ks(rCounter)]=kspEngine.Calculate(cell_idx);
    
    
    %random_delta_mean(rCounter)=...
    %    100*abs(mean(marker_vals_mat(cell_idx))-whole_tissue_mean)/whole_tissue_mean;
    random_delta_median(rCounter)=...
        abs((nnz(marker_vals_mat(cell_idx)>whole_tissue_median)/nPoints)-0.5);
end



%% KSP Distribution Calculations (for Fig 1.d)
min_ksp=nanmin(core_ks_primes(:));
max_ksp=nanmax(core_ks_primes(:));

ksp_bins=linspace(0,max_ksp,number_of_ksp_bins);
ksp_bin_width=ksp_bins(2)-ksp_bins(1);



ksp_pdf=cell(max_number_of_cores,1);
ksp_cdf=cell(max_number_of_cores,1);
for numberOfCores=1:max_number_of_cores
    
    ksp_pdf{numberOfCores}=histc(core_ks_primes(numberOfCores,:),ksp_bins);
    ksp_pdf{numberOfCores}=ksp_pdf{numberOfCores}/sum(ksp_pdf{numberOfCores});
    ksp_cdf{numberOfCores}=cumsum(ksp_pdf{numberOfCores});
    ksp_pdf{numberOfCores}=ksp_pdf{numberOfCores}/ksp_bin_width;
end


ks_pdf=cell(max_number_of_cores,1);
ks_cdf=cell(max_number_of_cores,1);
for numberOfCores=1:max_number_of_cores
    
    ks_pdf{numberOfCores}=histc(core_ks(numberOfCores,:),ksp_bins);
    ks_pdf{numberOfCores}=ks_pdf{numberOfCores}/sum(ks_pdf{numberOfCores});
    ks_cdf{numberOfCores}=cumsum(ks_pdf{numberOfCores});
    ks_pdf{numberOfCores}=ks_pdf{numberOfCores}/ksp_bin_width;
end


scatter_data=cell(max_number_of_cores,1);
scatter_data_ks=cell(max_number_of_cores,1);
for numberOfCores=1:max_number_of_cores
    scatter_data{numberOfCores}=[ core_ks_primes(numberOfCores,:)',...
        core_delta_medians(numberOfCores,:)'];
   scatter_data_ks{numberOfCores}=[ core_ks_primes(numberOfCores,:)',...
        core_ks(numberOfCores,:)'];
    %scatter_data{number_of_cores}=[ core_ks_primes(number_of_cores,:)',...
    %    100*abs(core_means(number_of_cores,:)-whole_tissue_mean)'/whole_tissue_mean];
    
end


%random core data

for numberOfCores=1%:max_number_of_cores
    
    random_ksp_pdf=histc(random_ksp,ksp_bins);
    random_ksp_pdf= random_ksp_pdf/sum( random_ksp_pdf);
    random_ksp_cdf=cumsum(random_ksp_pdf);
    random_ksp_pdf=random_ksp_pdf/ksp_bin_width;
    
    
    random_ks_pdf=histc(random_ks,ksp_bins);
    random_ks_pdf= random_ks_pdf/sum( random_ks_pdf);
    random_ks_cdf=cumsum(random_ks_pdf);
    random_ks_pdf=random_ks_pdf/ksp_bin_width;
    
end






%% Fig 1b top: Display Images Whole Tissue Image with tumor and core outline + individual cores


% Calculate colors used to plot core borders based on their KS-Prime values
numberOfCoresDisplayed=size(chosen_core_centers,1);
coreOutlineColors=zeros(numberOfCoresDisplayed,3);
for coreCounter=1:numberOfCoresDisplayed
    core_coords=chosen_core_centers(coreCounter,:);
    dists_to_cores=pdist2(core_coords,centroids_mat);
    [~,cell_idx]=find(dists_to_cores<core_radius);
    
    kVal=kspEngine.Calculate(cell_idx);
    cNum=round((size(kspColors,1)-1)*((kVal-kspRange(1))/(kspRange(2)-kspRange(1))).^1)+1;
    if(cNum<1)
        cNum=1;
    end
    if(cNum>size(kspColors,1))
        cNum=size(kspColors,1);
    end
    coreOutlineColors(coreCounter,:)=kspColors(cNum,:);
end


figure;
imagesc(img_lowres,display_scaling_range);colormap('gray');
hold on;


for coreCounter=1:numberOfCoresDisplayed
    
    rectangle('Position',[chosen_core_centers(coreCounter,1)/8 - core_radius/8, ...
        chosen_core_centers(coreCounter,2)/8 - core_radius/8, core_radius/4, core_radius/4],...
        'Curvature',[1,1],...
        'EdgeColor',coreOutlineColors(coreCounter,:),'LineWidth',4,'LineStyle','-');
end
plot(tumor_outline(:,1)/8,tumor_outline(:,2)/8,'LineWidth',3,'Color',border_colors(coreCounter+1,:));

% Add scale bar.
scale_bar_pos=[850,2250];
plot(scale_bar_pos(:,1)+[0,0],scale_bar_pos(2)+[1,core_radius/4],...
    'w','LineWidth',5); 
%note the core_radius by 4 is because core_radius is calculated on the high
%res image. the low res image shown is 8 fold smaller.

hold off;
axis off equal;


% Display Core Images separately

[core_imgs,core_stats]=Extract_Cores_From_Image(img,pixel_mask,chosen_core_centers,...
    core_radius,params.pp_cutoff,0:64:2^16); % Note the core_stats extracted here are pixel rather than cell based


for coreCounter=1:numberOfCoresDisplayed
    cmap1=gray(128);

    min_val=display_scaling_range(1);
    max_val=display_scaling_range(2);
    
    
    figure;
    temp_img= padarray(double(core_imgs{coreCounter}),[50,50],Inf);
    temp_img(temp_img==2^16-1)=Inf;
    
    temp_img(isfinite(temp_img)&temp_img>max_val)=max_val;
    temp_img=1+floor(double(temp_img)*128/(max_val+1));
    
    
    imagesc(temp_img,[1,128]);
    
    
    
    
    
    hold on;
    rectangle('Position',[size(temp_img,1)/2+core_radius*[-1 -1], 2*core_radius*[1,1] ]  ,...
        'Curvature',[1,1],...
        'EdgeColor',coreOutlineColors(coreCounter,:),'LineWidth',4,'LineStyle','-');
    hold off;
    colormap(cmap1);
    axis off equal;
end



%% Fig, 1b bottom: Histograms of selected cores
markerSize=150;
%hist_bins=2.^linspace(10,16,64);


wholeH=(histc(marker_vals_mat,hist_bins)/length(marker_vals_mat))./bin_widths';


plot(hist_bins,wholeH,'Color',wholeColor,'LineWidth',3);
hold on;
%plot(mean(marker_vals_mat),-0.0015,'Marker','^','Clipping',q'off','MarkerSize',18,...
%     'MarkerEdgeColor','k','MarkerFaceColor',border_colors(core_counter+1,:),'LineWidth',3);
scatter(mean(marker_vals_mat),-0.000015,markerSize,'filled','Marker','^','Clipping','off',...
    'MarkerEdgeColor','k','MarkerFaceColor',wholeColor,'MarkerFaceAlpha',0.5);
core_hist=cell(size(params.displayed_core_centers,1),1);
for core_counter=1:size(params.displayed_core_centers,1)
    
    core_coords=chosen_core_centers(core_counter,:);
    %core_coords=fliplr(core_coords);
    dists_to_cores=pdist2(core_coords,centroids_mat);
    [tma_idx,cell_idx]=find(dists_to_cores<core_radius);
    core_intensity_vals=marker_vals_mat(cell_idx);
    core_hist{core_counter}=histc(core_intensity_vals,hist_bins);
    core_hist{core_counter}=(core_hist{core_counter}/sum(core_hist{core_counter}))./bin_widths';
    %plot(mean(vals),-0.0015,'Marker','^','Clipping','off','MarkerSize',18,...
    %    'MarkerEdgeColor','k','MarkerFaceColor',border_colors(core_counter,:),'LineWidth',3);
    scatter(mean(core_intensity_vals),-0.000015,markerSize,'filled','Marker','^','Clipping','off',...
        'MarkerEdgeColor','k','MarkerFaceColor',coreOutlineColors(core_counter,:),'MarkerFaceAlpha',0.5);
    
    plot(hist_bins,core_hist{core_counter},'Color',coreOutlineColors(core_counter,:),...
        'LineWidth',3);
    
end
hold off;
ylim=get(gca,'YLim');
ylim(1)=0;
set(gca,'Xscale','log','YLim',ylim,'Clipping','off','XTickLabel',[],'YTickLabel',[]);
Fig1bData=table(hist_bins',wholeH,core_hist{1},core_hist{2},'VariableNames',...
    {'TTF1_Intensity','Whole_PDF','Core1_PDF','Core2_PDF'});
writetable(Fig1bData,'Fig1.xls','Sheet','Fig1b');
%% Fig 1c: Plots of Intensity histograms colored by KSP values for varying #Cores



for nrCoreCounter=1:length(nr_cores_to_show)

    figure;
    
    core_pdf_collection=cell(nrRands,1);
    for rCounter=1:nrRands
        coreIdx= randi(1000);
        core_pdf=core_nuclear_hists{nr_cores_to_show(nrCoreCounter),coreIdx};
        core_pdf=(core_pdf/sum(core_pdf))./bin_widths';
        kVal=core_ks_primes(nr_cores_to_show(nrCoreCounter),coreIdx);
        
        cNum=round((size(kspColors,1)-1)*((kVal-kspRange(1))/(kspRange(2)-kspRange(1))).^1)+1;
        if(cNum<1)
            cNum=1;
        end
        if(cNum>size(kspColors,1))
            cNum=size(kspColors,1);
        end
        
        
        plot(hist_bins,core_pdf,'Color',kspColors(cNum,:),'LineWidth',2);
        core_pdf_collection{rCounter}=core_pdf;
       
        hold on;
    end
    plot(hist_bins,whole_tissue_pdf,'Color',wholeColor(1,:),'LineWidth',5);
    hold off;
    set(gca,'YLim',[0,1E-3],'XTick',[],'YTick',[],'XScale','log');
    %title(['Core Based Sampling:' num2str(nrCoreVals(nrCoreCounter)) ' cores']);
    colNames=[{'TTF1_Intensity','Whole_PDF'},...
        strcat('Sampling',arrayfun(@num2str,1:nrRands,'Unif',false),'_PDF')];
    sheetName=['Fig1c_' num2str(nr_cores_to_show(nrCoreCounter)) 'Cores'];
    temp_table=table(hist_bins',whole_tissue_pdf,core_pdf_collection{:},...
        'VariableNames',colNames);
    writetable(temp_table,'Fig1.xls','Sheet',sheetName);
end


for nrCoreCounter=1%:length(nr_cores_to_show)

    figure;
    
    core_pdf_collection=cell(nrRands,1);
    for rCounter=1:nrRands
        
        nPoints=length(core_nuclear_vals{nr_cores_to_show(nrCoreCounter),randi(1000)});
        idx=randsample(length(marker_vals_mat),nPoints);
        core_pdf=histc(marker_vals_mat(idx),hist_bins);
        core_pdf=(core_pdf/sum(core_pdf))./bin_widths';
        
        
        [kVal,~]=kspEngine.Calculate(idx);
        
        cNum=round((size(kspColors,1)-1)*((kVal-kspRange(1))/(kspRange(2)-kspRange(1))).^1)+1;
        if(cNum<1)
            cNum=1;
        end
        if(cNum>size(kspColors,1))
            cNum=size(kspColors,1);
        end
        
        %plot(hist_bins,core_pdf,'Color',[0.5 0.5 0.5],'LineWidth',2);
        plot(hist_bins,core_pdf,'Color',kspColors(cNum,:),'LineWidth',2);
        core_pdf_collection{rCounter}=core_pdf;
        hold on;
    end
    plot(hist_bins, whole_tissue_pdf,'Color',wholeColor(1,:),'LineWidth',5);
    hold off;
    set(gca,'YLim',[0,1E-3],'XTick',[],'YTick',[],'XScale','log');
   
    colNames=[{'TTF1_Intensity','Whole_PDF'},...
        strcat('Sampling',arrayfun(@num2str,1:nrRands,'Unif',false),'_PDF')];
    sheetName='Fig1c_RandomSampling';
    temp_table=table(hist_bins',whole_tissue_pdf,core_pdf_collection{:},...
        'VariableNames',colNames);
    writetable(temp_table,'Fig1.xls','Sheet',sheetName);
    %title(['Random Sampling:' num2str(nrCoreVals(nrCoreCounter)) ' cores worth of cells']);
end
figure;
imagesc([1:length(kspColors)]');
colormap(kspColors);
axis off;

%% Fig.1d top:Scatter Plot of Delta-Median vs KSP colored by number of cores
showRandom=true;

xmax=kspRange(2);

for i=1:length(nr_cores_to_show)
    numberOfCores=nr_cores_to_show(i);
    
    alpha=1.0;
    scatter(scatter_data{numberOfCores}(:,1),scatter_data{numberOfCores}(:,2),...
        30,'filled','MarkerFaceColor',nrCores_colors(numberOfCores,:),...
        'MarkerEdgeColor','none', 'MarkerFaceAlpha',alpha);
    
    hold on;
end
if(showRandom)
    scatter(random_ksp,random_delta_median,...
        30,'filled','MarkerFaceColor',randomSamplingColor,...
        'MarkerEdgeColor','none', 'MarkerFaceAlpha',1.0);
 
end
set(gca,'XLim',[0,xmax],'YLim',[0,xmax]);
hold on;
set(gca,'XTick',[0,xmax/2,xmax],'YTick',[0,xmax/2,xmax],'XTickLabel',[],...
    'YTickLabel',[]);
hold off;
hold on;
axis square;

dataPlotted=scatter_data(nr_cores_to_show);
dataPlotted=vertcat(dataPlotted{:});
dataKsp=[dataPlotted(:,1);random_ksp];
dataMedianError=[dataPlotted(:,2);random_delta_median];
dataNrCoresUsed=repmat(nr_cores_to_show,number_of_reps,1);
dataNrCoresUsed=[dataNrCoresUsed(:); NaN*ones(number_of_reps,1)];
dataComment=cell(size(dataKsp));
dataComment(1:(length(nr_cores_to_show)*number_of_reps))={'Core Based Sampling'};
dataComment((length(nr_cores_to_show)*number_of_reps)+(1:number_of_reps))=...
    {'Random Sampling, 1 core equivalent'};
Fig1dTopTable=table(dataKsp,dataMedianError,dataNrCoresUsed,dataComment,...
    'VariableNames',{'KS_Prime','Median_Error','Number_Of_Cores','Sampling_Method'});
writetable(Fig1dTopTable,'Fig1.xls','Sheet','Fig1d_Top');
%temp=Fig1dTopTable.Number_Of_Cores;
%temp(isnan(temp))=20;
%scatter(Fig1dTopTable.KS_Prime,Fig1dTopTable.Median_Error,50,...
%    temp,'filled');
%% Fig.1d bottom: KS-Prime CDFs
figure
for i=1:length(nr_cores_to_show)
    numberOfCores=nr_cores_to_show(i);
    plot(ksp_bins,100*ksp_cdf{numberOfCores},'Color',nrCores_colors(numberOfCores,:),...
        'LineWidth',3) ;
    hold on;
end
if(showRandom)
    plot(ksp_bins,100*random_ksp_cdf,'Color',randomSamplingColor,'LineWidth',3) ;
end
set(gca,'XLim',[0,xmax],'YLim',[0,100],'XTickLabel',[],'YTickLabel',[]);

temp=cellfun(@(x) 100*x',ksp_cdf(nr_cores_to_show),'Unif',false);
   
colNames=[{'KS_Prime'},...
        strcat('CDF_',arrayfun(@num2str,nr_cores_to_show,'Unif',false),'Core'),{'CDF_Random'}];
Fig1dBottomTable=table(ksp_bins',temp{:},100*random_ksp_cdf,...
    'VariableNames',colNames);
writetable(Fig1dBottomTable,'Fig1.xls','Sheet','Fig1d_Bottom');
%plot(Fig1dBottomTable.KS_Prime,table2array(Fig1dBottomTable(:,2:end)));
%set(gca,'XLim',kspRange);