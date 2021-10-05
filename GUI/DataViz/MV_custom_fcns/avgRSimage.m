function avgRSimage(stats_data)
% This function aggregates imaging data in "dat_filename" across subjects.
% Inputs:
% stats_data (struct) = structure from StatsManager class containing all
% information of each recording and observations.

disp('Go!')

%
acqIdx_list = [stats_data.acqIdx];
obs_list = arrayfun(@(x) {x.observations.ID}, stats_data,'UniformOutput',0);
obs_list = unique(cellfun(@(x) num2str(x),[obs_list{:}],'UniformOutput',0));
idx = cellfun(@(x) isempty(x), obs_list);
obs_list(idx) = [];
groups = unique({stats_data.groupName});


for i = 1:numel(obs_list)
    f = figure('Name', obs_list{i}, 'Visible', 'off');
    a=1;
    for j = 1:max(acqIdx_list)
        for k = 1:numel(groups)
            idx = [stats_data.acqIdx] == j & strcmp({stats_data.groupName},groups{k}) &...
            arrayfun(@(x) any(strcmp(obs_list{i},{x.observations.ID})),stats_data);
            [img,N] = align_and_average(stats_data(idx), obs_list{i},1);
            s = subplot(max(acqIdx_list),numel(groups),a);
            imagesc(s,img);
%             s.CLim = [-1 1];
            s.DataAspectRatio = [1 1 1]; colormap(s,'jet');
            title(s,['Group ' groups{k} ' Acq ' num2str(j) ' N = ' num2str(N)],...
                'Interpreter','none');
            a = a+1;
        end
    end
    
    % Normalize CLims:
    all_ax = findall(f, 'Type', 'Axes');
    data = arrayfun(@(x) x.Children.CData(:),all_ax, 'UniformOutput',0);
    data = vertcat(data{:});
    idx = isinf(data) | isnan(data);
    clims = [min(data(~idx)) max(data(~idx))];
    arrayfun(@(x) set(x, 'CLim', clims), all_ax);
    % Set figure visible:
    f.Visible = 'on'; drawnow
    % End!
    disp(['Done with ' obs_list{i}])
end
disp('Finished!')
end



% Local Function

function [out,N] = align_and_average(data, obsID, b_applyGauss)
% This function aligns all frames and averages the images.
disp(['Calculating average for ' obsID '...'])
% Retrieve ROI's centroids:
out = {};
refPt = {};
% Set first item of data as reference"
roifile = load(data(1).MatFile.ROIfile);
obs_indx = find(strcmp(obsID,{roifile.ROI_info.Name}));

[fY,fX]= find(bwmorph(roifile.ROI_info(obs_indx).Stats.ROI_binary_mask, 'shrink', Inf));
dist1 =  sqrt((abs((roifile.img_info.refPt(1)) - fX))^2 + (abs((roifile.img_info.refPt(2)) - fY))^2);
ref2d = imref2d(size(roifile.ROI_info(obs_indx).Stats.ROI_binary_mask));
imgs = zeros([ref2d.ImageSize, numel(data)], 'single');
for i = 1:length(data)
    roifile = load(data(i).MatFile.ROIfile);
    origFile = data(i).MatFile.datFile; origFile = origFile{1};
    spcm_file = mapDatFile(origFile);
    spcm_data = spcm_file.Data.CM;
    obs_indx = find(strcmp(obsID,{roifile.ROI_info.Name}));
    seed_fr = arrayfun(@(x) find(bwmorph(x.Stats.ROI_binary_mask, 'shrink', Inf)),...
        roifile.ROI_info(obs_indx));
    in = spcm_data(:,:,seed_fr);
    [mvY,mvX] = ind2sub([size(in,1),size(in,2)], seed_fr);
    % Try to infer scaling by calculating difference of distance between
    % Bregma and centroid:
%     distM =  sqrt((abs((roifile.img_info.refPt(1)) - mvX))^2 + (abs((roifile.img_info.refPt(2)) - mvY))^2);
%     scale_factor = distM/dist1
    scale_factor = 1;
    refPt{i} = roifile.img_info.refPt;    
    % Temporarily replace NaNs and Infs with zeros:
    msk = isnan(in);
    msk_inf = isinf(in);
    in(msk|msk_inf) = 0;
    % Translate/Scale to match bregma position:
    % 1- Scale:
    tf = affine2d([scale_factor 0 0;  0 scale_factor 0; 0 0 1]);
    tmp = imwarp(in,tf,'nearest', 'OutputView', ref2d);
    msk = imwarp(msk,tf,'nearest', 'OutputView', ref2d, 'FillValues',1);
    msk_inf = imwarp(msk_inf,tf,'nearest', 'OutputView', ref2d, 'FillValues',0);
    % 2- Find "scaled" bregma:
    [x,y] = transformPointsForward(tf,refPt{i}(1), refPt{i}(2));
    % Translate to Reference Bregma:
    txy = refPt{1} - [x,y];
    tmp = imtranslate(tmp,ref2d,txy,'nearest');
    msk = imtranslate(msk,ref2d,txy,'nearest', 'FillValues',1);
    msk_inf= imtranslate(msk_inf,ref2d,txy,'nearest', 'FillValues',0);
    % Replace NaNs and Infs:
    tmp(msk) = nan;
    tmp(msk_inf)= inf;
    % apply spatial filter:
    if b_applyGauss
        % Remove Infs:
        tmp(msk_inf) = max(tmp(~msk_inf));
        tmp = imgaussfilt(tmp,2);
    end
    imgs(:,:,i) = tmp;
end
% Average images:
N = size(imgs,3);
% out = any(imgs == Inf,3); % Checking alignment. 
out = mean(imgs,3,'omitnan');
end