function outFile = getDataFromROI(File, SaveFolder, varargin)
% GETDATAFROMROI extracts and aggregates data from regions of interest
% (ROIs) in imaging data using an "ROI_xxxxxx.mat" file located in
% subject's folder.
%
% Inputs:
%   File: Imaging data file with 2+ dimensions (.DAT file).
%   SaveFolder: Folder where data are stored.
%   Output: name of .DAT file saved in SAVEFOLDER.
% Outputs:
%   outFile: .MAT file containing ROI names and aggregate data.
%
% Defaults:
default_Output = 'ROI_data.mat'; 

%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p, 'File', @(x) isfile(x) && endsWith(x, '.dat'));
addRequired(p, 'SaveFolder', @isfolder);
% Optional Parameters:
% opts structure:
default_opts = struct('SpatialAggFcn', 'mean');
addOptional(p, 'opts', default_opts,@(x) isstruct(x) &&...
    ismember(x.SpatialAggFcn, {'none','mean', 'max', 'min', 'median', 'mode', 'sum', 'std'}));
% Parse inputs:
parse(p,File, SaveFolder, varargin{:});
% Initialize Variables:
File = p.Results.File;
SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;

% map File:
[mData, metaData] = mapDatFile(File);
datName = metaData.datName;
if iscell(datName)
    datName = datName{1};
end
data = mData.Data.(datName);
data_sz = size(data);
% Parse File path to find subject folder:
str = split(File, filesep);
subjFolder = strjoin(str(1:end-3), filesep);
ROIfile = dir([subjFolder filesep 'ROI_*.mat']);
if isempty(ROIfile)
    errID = 'Umitoolbox:getDataFromROI:FileNotFound';
    errMsg = ['ROI file not found in ' subjFolder];
    error(errID, errMsg);
end
% Load ROI file:
roi_data = load(fullfile(ROIfile.folder, ROIfile.name));
% locate "X" and "Y" dimensions in metaData:
dim_names = metaData.dim_names;
xLoc = find(strcmp(dim_names, 'X'));
yLoc = find(strcmp(dim_names, 'Y'));
% Check if frame size is the same as the one in ROI file:
errID = 'Umitoolbox:getDataFromROI:IncompatibleSizes';
errMsg = 'Data file frame size is different from the one in ROI file';
assert(isequal([data_sz(yLoc) data_sz(xLoc)], size(roi_data.img_info.data)), errID, errMsg)
% permute matrix:
orig_dim = 1:ndims(mData.Data.data);
new_dim = [xLoc yLoc setdiff(orig_dim, [xLoc yLoc])];
data = permute(data, new_dim);
dim_names = dim_names(new_dim);
% reshape matrix:
data_sz = data_sz(new_dim);
data = reshape(data, prod(data_sz([1 2])),[]);
% extract ROI pixel values from data:
roi_names = {roi_data.ROI_info.Name}';
roi_pixVals = cell(size(roi_names));
for i = 1:length(roi_pixVals)
    roi_msk = roi_data.ROI_info(i).Stats.ROI_binary_mask;
    pixVals = data(roi_msk(:),:);
    % Perform aggregate operation:
    pixVals = applyAggFcn(pixVals, opts.SpatialAggFcn);
    % reshape data back to retrieve dimensions:
    if length(data_sz) > 2
        pixVals = reshape(pixVals, [size(pixVals,1) data_sz(3:end)]);
    end
    roi_pixVals{i} = pixVals;
end
% Save data to .mat:
filename = fullfile(SaveFolder, default_Output);

new_dim_names ={'O', dim_names{3:end}};
save2Mat(filename, roi_pixVals, roi_names, new_dim_names, 'appendMetaData', metaData)
end

% Local function:
function out = applyAggFcn(vals, fcn_name)
% APPLYAGGFCN performs the aggregate function of name "fcn_name" on the 1st
% dimension of the data "vals". All aggregate functions EXCLUDE NaNs!

switch fcn_name
    case 'mean'
        out = nanmean(vals, 1);
    case 'median'
        out = median(vals, 1, 'omitnan');
    case 'mode'
        out = mode(vals, 1);
    case 'std'
        out = std(vals, 0, 1, 'omitnan');
    case 'max' 
        out = max(vals, [], 1, 'omitnan');
    case 'min'
        out = min(vals, [], 1, 'omitnan');
    case 'sum'
        out = sum(vals, 1, 'omitnan');
    otherwise 
        out = vals;
end
end










