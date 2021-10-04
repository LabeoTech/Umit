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
            [img,N] = align_and_average(stats_data(idx), obs_list{i});
            s = subplot(max(acqIdx_list),numel(groups),a);
            imagesc(s,img);
            s.CLim = [-1 1]; s.DataAspectRatio = [1 1 1]; colormap(s,'jet');
            title(s,['Acq ' num2str(j) ' N = ' num2str(N)]);
            a = a+1;
        end
    end
    f.Visible = 'on'; drawnow
    disp(['Done with ' obs_list{i}])
end
disp('Finished!')
end



% Local Function

function [out,N] = align_and_average(data, obsID)
% This function aligns all frames and averages the images.
disp(['Calculating average for ' obsID '...'])
% Retrieve ROI's centroids:
out = {};
refPt = {};
% Set first item of data as reference"
roifile = load(data(1).MatFile.ROIfile);
obs_indx = find(strcmp(obsID,{roifile.ROI_info.Name}));
[fY,fX]= find(bwmorph(roifile.ROI_info(obs_indx).Stats.ROI_binary_mask, 'shrink', Inf));
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
%     refPt{i} = data(i).MatFile.refPt;
    refPt{i} = roifile.img_info.refPt;
    tf = fitgeotrans([refPt{i};mvY, mvX], [refPt{1};fY,fX],'nonreflectivesimilarity');
    % Remove NaNs from image:
    msk = isnan(in);
    in(msk) = 0;
    tmp = imwarp(in,tf,'nearest','OutputView', ref2d);
    tmp(msk) = nan;
    imgs(:,:,i) = tmp;
end
% Average images:
N = size(imgs,3);
% out = any(imgs == 1,3);
out = mean(imgs,3,'omitnan');
end