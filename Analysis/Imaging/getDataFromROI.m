function outFile = getDataFromROI(data, metaData, SaveFolder, varargin)
% GETDATAFROMROI extracts and aggregates data from regions of interest
% (ROIs) in imaging data using an "ROI_xxxxxx.mat" file located in
% subject's folder.

% Inputs:
%   data: numerical matrix containing imaging data.
%   metaData: .mat file with meta data associated with "data".
%   SaveFolder: Folder where data are stored.
% Outputs:
%   outFile: .MAT file containing ROI names and aggregate data.

% Defaults:
default_Output = 'ROI_data.mat'; %#ok
default_opts = struct('ROI_filename', 'ROI_data.mat', 'SpatialAggFcn', 'mean');
%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'data',@(x) isnumeric(x)); % Validate if the input is a 3-D numerical matrix:
addRequired(p,'metaData', @(x) isa(x,'matlab.io.MatFile') | isstruct(x)); % MetaData associated to "data".
addRequired(p, 'SaveFolder', @isfolder);
% Optional Parameters:
addOptional(p, 'opts', default_opts,@(x) isstruct(x) &&...
    ismember(x.SpatialAggFcn, {'none','mean', 'max', 'min', 'median', 'mode', 'sum', 'std'}));
addOptional(p, 'object', [], @(x) isempty(x) || isa(x,'Acquisition') || isa(x,'Modality')); 

% Parse inputs:
parse(p,data, metaData, SaveFolder, varargin{:});
% Initialize Variables:
data = p.Results.data;
metaData = p.Results.metaData;
SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;
object = p.Results.object;
clear p
%%%%%%%%%%%%%%%%

% For backward compatibility with previous versions of umIT when .DAT files contained 
% more than one matrix:
datName = metaData.datName;
dim_names = metaData.dim_names;
if iscell(datName)
    dim_names = dim_names{1};
end
%%%%%%%%%%%%%%%

% Parse File path to find subject folder:
if ~isempty(object)
    % If a umIT's valid object is provided, it means that the function is being run by 
    % PipelineManager. In this case, find the subject's folder an try to
    % find the ROI file there: 
    idx = false;
    while ~idx
        ParentObj = object.MyParent;
        idx = isa(ParentObj, 'Subject');
    end  
    if ~isfile(fullfile(ParentObj.SaveFolder, opts.ROI_filename))
        errID = 'Umitoolbox:getDataFromROI:FileNotFound';
        subjFolder = strrep(ParentObj.SaveFolder, '\', '\\');
        errMsg = ['ROI file not found in ' subjFolder];
        error(errID, errMsg);
    else
        opts.ROI_filename = fullfile(ParentObj.SaveFolder, opts.ROI_filename);
    end
end

% Load ROI file:
roi_data = load(opts.ROI_filename);
% locate "X" and "Y" dimensions in metaData and in ROI info:
xLoc = find(strcmp(dim_names, 'X'));
yLoc = find(strcmp(dim_names, 'Y'));
% Check if frame size is the same as the one in ROI file:
data_sz = size(data);
errID = 'Umitoolbox:getDataFromROI:IncompatibleSizes';
errMsg = 'Data file frame size is different from the one in ROI file';
assert(isequal([data_sz(yLoc) data_sz(xLoc)], size(roi_data.img_info.data)), errID, errMsg)
% permute matrix:
orig_dim = 1:ndims(data);
new_dim = [yLoc xLoc setdiff(orig_dim, [xLoc yLoc])];
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
[~, datFile,~] = fileparts(metaData.datFile);
[~,roi_filename,~] = fileparts(opts.ROI_filename);
outFile = [roi_filename '_' datFile '.mat'];
mFile = fullfile(SaveFolder, outFile);
new_dim_names ={'O', dim_names{3:end}};
if isa(metaData, 'matlab.io.MatFile')
    metaData.Properties.Writable = true;
else
    metaData.Writable = true;
end
metaData.ROIfile = opts.ROI_filename;
save2Mat(mFile, roi_pixVals, roi_names, new_dim_names, 'appendMetaData', metaData)
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










