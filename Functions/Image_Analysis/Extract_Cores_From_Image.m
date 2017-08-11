function [core_imgs,core_stats]=Extract_Cores_From_Image(img,foreground_mask,...
        core_centers,core_radius,pp_cutoff,hist_bins,varargin)
    p=inputParser;
    
    default_outside_core_val=Inf;
    addParameter(p,'outside_core_val',default_outside_core_val);
    default_force_double_output=false;
    addParameter(p,'force_double_output',default_force_double_output);
    
    parse(p,varargin{:});
    
    outside_core_val=p.Results.outside_core_val;
    force_double_output=p.Results.force_double_output;
    
    number_of_cores=size(core_centers,1);
    
    [x_offsets,y_offsets]=meshgrid(-core_radius:core_radius,-core_radius:core_radius);
    x_offsets=x_offsets(:);
    y_offsets=y_offsets(:);
    dists=sqrt(x_offsets.^2+y_offsets.^2);%chosen_core_centers=8*[1421 705;1777 2220];round(8*[1421 705;1934 2300]);%1342,1692; % Coordinates of tma cores displayed in full mag image
    
    x_offsets=x_offsets(dists<=core_radius);
    y_offsets=y_offsets(dists<=core_radius);
    
    core_stats=struct;
    core_imgs=cell(number_of_cores,1);
    number_of_channels=size(img,3);
    
    for core_number=1:number_of_cores
        if(force_double_output)
            core_imgs{core_number}=zeros(2*core_radius+1,2*core_radius+1,number_of_channels);
           
        else
            core_imgs{core_number}=zeros(2*core_radius+1,2*core_radius+1,...
                number_of_channels,'like',img);

        end
        core_stats(core_number).mean=zeros(number_of_channels,1);
        core_stats(core_number).ppp=zeros(number_of_channels,1);
        core_stats(core_number).pdf=zeros(number_of_channels,length(hist_bins));
        core_stats(core_number).cdf=zeros(number_of_channels,length(hist_bins));
        for channel_number=1:number_of_channels
            channel_img=img(:,:,channel_number);
            
            idx_of_pixels_in_core=sub2ind(size(channel_img),...
                x_offsets+core_centers(core_number,2),...
                y_offsets+core_centers(core_number,1));
            
            
            pixel_intensities=channel_img(idx_of_pixels_in_core);
            fg_pixels_in_core=foreground_mask(idx_of_pixels_in_core);
            fg_pixel_intensities=pixel_intensities(fg_pixels_in_core);
            
            core_stats(core_number).mean(channel_number)=mean(fg_pixel_intensities);
            core_stats(core_number).ppp(channel_number)=nnz(fg_pixel_intensities>pp_cutoff)/...
                numel(fg_pixel_intensities);
            if(~force_double_output)
                temp=channel_img(core_centers(core_number,2)+...
                    (-core_radius:core_radius),core_centers(core_number,1)+(-core_radius:core_radius));
            else
                temp=double(channel_img(core_centers(core_number,2)+...
                    (-core_radius:core_radius),core_centers(core_number,1)+(-core_radius:core_radius)));
            end
            temp(dists>core_radius)=outside_core_val;
            core_imgs{core_number}(:,:,channel_number)=temp;
            core_stats(core_number).pdf(channel_number,:)=...
                histc(fg_pixel_intensities,hist_bins);
            core_stats(core_number).pdf(channel_number,:)=...
                core_stats(core_number).pdf(channel_number,:)...
                /sum(core_stats(core_number).pdf(channel_number,:));
            core_stats(core_number).cdf(channel_number,:)=...
                cumsum(core_stats(core_number).pdf(channel_number,:));
        end
    end
    
end
