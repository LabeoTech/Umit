function outDataStat = genCorrelationMatrix(dataStat)
% GENCORRELATIONMATRIX creates a correlation matrix from SPCM data extracted using
% getDataFromROI function.

% Inputs:
%   dataStat: structure containing data from a statistics function.
% Output:
%   outDataStat: structure containing stats-ready data extracted from ROIs.

% Defaults:
default_Output = 'CorrMatrix.mat'; %#ok
%%% Arguments parsing and validation %%%
p = inputParser;
% The input of the function must be a File , RawFolder or SaveFolder
% Imaging Object:
addRequired(p, 'dataStat', @(x) isstruct(x))
% Parse inputs:
parse(p,dataStat);
% Initialize Variables:
data = p.Results.dataStat;
clear p
%%%%

% Select observations in data.ROIfile:
idx_obs = ismember({data.ROIfile.ROI_info.Name}, data.obsID);
data.ROIfile.ROI_info(~idx_obs)= [];
% Find seed frames in data.ROIfile masks:
centroid_px = arrayfun(@(x) find(bwmorph(x.Stats.ROI_binary_mask,'shrink', Inf),...
    1, 'last'), data.ROIfile.ROI_info, 'UniformOutput', true);
if size(data.obsID,1)~= size(centroid_px,1)
    centroid_px = centroid_px';
end
% Rearrange centroid_px to match data.obsID list
[~,lb] = ismember({data.ROIfile.ROI_info.Name}, data.obsID);
centroid_px = centroid_px(lb);
% Fill the correlation matrix with values:
corrMat = cellfun(@(x) single(x(centroid_px)),data.data, 'UniformOutput',0);
% Create new dimension names as {'O', 'O'}:
dim_names = {'O','O'};
outDataStat = save2Mat([], corrMat,data.obsID, dim_names, 'label',data.obsID , 'appendMetaData', data,'genFile', true);
end
