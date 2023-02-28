function outData = getDataFromROI(data, metaData, varargin)
% GETDATAFROMROI extracts and aggregates data from regions of interest
% (ROIs) in imaging data using an "ROI_xxxxxx.mat" file located in
% subject's folder.

% Inputs:
%   data: numerical matrix containing imaging data.
%   metaData: .mat file with meta data associated with "data".
% Output:
%   outData: structure containing stats-ready data extracted from ROIs.

% Defaults:
default_Output = 'ROI_data.mat'; %#ok This line is here just for Pipeline management.
default_opts = struct('ROImasks_filename', 'ROImasks_data.mat', 'SpatialAggFcn', 'mean');
opts_values = struct('ROImasks_filename', {{'ROImasks_data.mat'}}, 'SpatialAggFcn',{{'none','mean', 'max', 'min', 'median', 'mode', 'sum', 'std'}});% This is here only as a reference for PIPELINEMANAGER.m.
default_object = ''; % This line is here just for Pipeline management to be able to detect this input.
%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'data',@(x) isnumeric(x)); % Validate if the input is a 3-D numerical matrix:
addRequired(p,'metaData', @(x) isa(x,'matlab.io.MatFile') | isstruct(x)); % MetaData associated to "data".
% Optional Parameters:
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ...
    ismember(x.SpatialAggFcn, opts_values.SpatialAggFcn));
addOptional(p, 'object', default_object, @(x) isempty(x) || isa(x,'Acquisition') || isa(x,'Modality'));

% Parse inputs:
parse(p,data, metaData, varargin{:});
% Initialize Variables:
data = p.Results.data;
metaData = p.Results.metaData;
opts = p.Results.opts;
object = p.Results.object;
clear p
%%%%%%%%%%%%%%%%

% Parse File path to find subject folder:
opts.ROImasks_filename = findMyROIfile(opts.ROImasks_filename,object);

% Load ROI file:
roi_data = load(opts.ROImasks_filename);
% locate "X" and "Y" dimensions in metaData and in ROI info:
dim_names = metaData.dim_names;
[~,yxLoc] = ismember({'Y','X'}, dim_names);

% Check if frame size is the same as the one in ROI file:
data_sz = size(data);
errID = 'Umitoolbox:getDataFromROI:IncompatibleSizes';
errMsg = 'Data file frame size is different from the one in ROI file.';
assert(isequal(data_sz(yxLoc), size(roi_data.img_info.imageData)), errID, errMsg)
% permute matrix:

orig_dim = 1:ndims(data);
new_dim = [yxLoc setdiff(orig_dim, yxLoc)];
data = permute(data, new_dim);
dim_names = dim_names(new_dim);
% reshape matrix:
data_sz = data_sz(new_dim);
data = reshape(data, prod(data_sz([1 2])),[]);
% extract ROI pixel values from data:
roi_names = {roi_data.ROI_info.Name}';
roi_pixVals = cell(size(roi_names));
disp(['Extracting ' opts.SpatialAggFcn ' values from ROIs...'])
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
disp('Done');
% Create new dimension names:
if strcmp(opts.SpatialAggFcn, 'none')
    % This is an special case where the User wants to keep all individual
    % pixels from the ROI. In this case, the dimension "P" (for pixel) is
    % created after the "O" (observation).
    new_dim_names = {'O', 'P', dim_names{3:end}};
else
    new_dim_names ={'O', dim_names{3:end}};
end
if isa(metaData, 'matlab.io.MatFile')
    metaData.Properties.Writable = true;
end
metaData.ROIfile = opts.ROImasks_filename;
if isempty(object)
    outData = save2Mat('', roi_pixVals, roi_names,...
        new_dim_names, 'appendMetaData', metaData, 'genFile', false);
else
    outData = save2Mat('', roi_pixVals, roi_names,...
        new_dim_names, 'appendMetaData', metaData, 'genFile', false,...
        'appendObjectInfo',object);
end

end

% Local function:
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
