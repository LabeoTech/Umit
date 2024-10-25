function outData = getDataFromROI(data, varargin)
% GETDATAFROMROI extracts and aggregates data from regions of interest
% (ROIs) in imaging data using an "ROI_xxxxxx.mat" file located in
% subject's folder.

% Inputs:
%   data: numerical matrix OR structure (see genDataMetaStructure.m) containing imaging data.
% Output:
%   outData: structure containing stats-ready data extracted from ROIs.

% Defaults:
default_Output = 'ROI_data.mat'; %#ok This line is here just for Pipeline management.
default_opts = struct('ROImasks_filename', 'myROIs.roimsk', 'SpatialAggFcn', 'mean');
opts_values = struct('ROImasks_filename', {{'myROIs.roimsk'}}, 'SpatialAggFcn',{{'none','mean', 'max', 'min', 'median', 'mode', 'sum', 'std'}});% This is here only as a reference for PIPELINEMANAGER.
%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'data',@(x) (isnumeric(x) && ndims(x)== 3) || isstruct(x)); % Validate if the input is a 3-D numerical matrix or a structure
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ...
    ismember(x.SpatialAggFcn, opts_values.SpatialAggFcn));
% Parse inputs:
parse(p,data, varargin{:});
% Initialize Variables:
data = p.Results.data;
opts = p.Results.opts;
clear p
%%%%%%%%%%%%%%%%
if isstruct(data)
    % Validate input:
    % Here, if the data is a structure, accept only if:
    %   1) There is a single observation
    %   2) The data from the observation is an image, image time series or both
    %       split by events.
    fn = fieldnames(data.data);
    assert(length(data.obsID) == 1 && ndims(data.data.(fn{1}))>= 2, 'Invalid input data! The data should be an Image or Image time series.');
    % Calculate on data saved in structure:
    outData = extractROIdata(data.data.(fn{1}),opts,data.b_hasEvents);
    fn = setdiff(fieldnames(data),fieldnames(outData));
    for ii = 1:length(fn)
        outData.(fn{ii}) = data.(fn{ii});
    end
else
    % Calculate on 3D array:
    outData = extractROIdata(data,opts,false);
end
disp('Done');
end
% Local functions:
function outStruct = extractROIdata(data,opts, b_hasEvents)
% Check if ROImasks file exist:
% opts.ROImasks_filename = findMyROIfile(opts.ROImasks_filename);
% Load ROI file:
roi_data = load(opts.ROImasks_filename,'-mat');
% locate "X" and "Y" dimensions:
if b_hasEvents
    yxLoc = [2,3];
else
    yxLoc = [1 2];
end

% Check if frame size is the same as the one in ROI file:
data_sz = size(data);
errID = 'Umitoolbox:getDataFromROI:IncompatibleSizes';
errMsg = 'Data file frame size is different from the one in ROI file.';
assert(isequal(data_sz(yxLoc), size(roi_data.img_info.imageData)), errID, errMsg)
% permute matrix:
orig_dim = 1:ndims(data);
new_dim = [yxLoc setdiff(orig_dim, yxLoc)];
data = permute(data, new_dim);
% reshape matrix:
data_sz = data_sz(new_dim);
data = reshape(data, prod(data_sz([1 2])),[]);
% extract ROI pixel values from data:
roi_names = {roi_data.ROI_info.Name}';
roi_pixVals = cell(1,length(roi_names));
disp(['Extracting ' opts.SpatialAggFcn ' values from ROIs...'])
for i = 1:length(roi_names)
    roi_msk = roi_data.ROI_info(i).Stats.ROI_binary_mask;
    pixVals = data(roi_msk(:),:);
    % Perform aggregate operation:
    pixVals = applyAggFcn(pixVals, opts.SpatialAggFcn);
    % reshape data back to retrieve dimensions:
    if length(data_sz) > 2
        pixVals = reshape(pixVals, [size(pixVals,1) data_sz(3:end)]);
    end
    roi_pixVals{i}= pixVals;
end
% Add ROImask to output data:
extraInfo.ROImasks_filename = opts.ROImasks_filename;
outStruct = genDataMetaStructure(roi_pixVals, roi_names,'hasEvents',b_hasEvents,'extraInfo',extraInfo);
end

function out = applyAggFcn(vals, fcn_name)
% APPLYAGGFCN performs the aggregate function of name "fcn_name" on the 1st
% dimension of the data "vals". All aggregate functions EXCLUDE NaNs!

switch fcn_name
    case 'mean'
        out = mean(vals, 1, 'omitnan');
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
