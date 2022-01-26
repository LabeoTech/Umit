function [outData, metaData] = calculate_response_amplitude(data, metaData, varargin)
% CALCULATE_RESPONSE_AMPLITUDE works with data split into events ("E" dimension in 
% "dim_names" variable). It calculates the signal amplitude in time ("T").
% The calculation consists on a value "postEvent_value" calculated in a time window
% between the event trigger and a given time in seconds ("timeWindow") minus a value
% "preEvent_value" calculated between the beginning of the trial an the Event frame.
% Default values for post and preEventTimes are 'max' and
% 'median'.
 
% Inputs:
%   data: numerical matrix containing image time series split by event (dimensions "E", "Y", "X", "T").
%   metaData: .mat file with meta data associated with "data".
%   opts (optional) : structure containing extra parameters.

% Outputs: 
%   outData: numerical matrix with dimensions {E,Y,X}.   
%   metaData: .mat file with meta data associated with "outData".

% Defaults:
default_Output = 'amplitude_Map.dat'; %#ok This line is here just for Pipeline management.
default_opts = struct('preEvent_value', 'median', 'postEvent_value', 'max', 'timeWindow', -1);

%%% Arguments parsing and validation %%%
% Parse inputs:
p = inputParser;
addRequired(p,'data',@(x) isnumeric(x)); % Validate if the input is a 3-D numerical matrix:
addRequired(p,'metaData', @(x) isa(x,'matlab.io.MatFile') | isstruct(x)); % MetaData associated to "data".
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
% Parse inputs:
parse(p,data, metaData, varargin{:});
%Initialize Variables:
data = p.Results.data; 
metaData = p.Results.metaData;
opts = p.Results.opts;
clear p
%%%%
% Further validation of optional arguments:
errID = 'Umitoolbox:calculate_respose_amplitude:InvalidDataType';
errMsg = 'Invalid Input. Input must be member of {"mean", "median", "min","max"} or a numeric scalar';
valid_Opts1 = @(x) (ismember(char(x), {'mean', 'median', 'min','max'}) || (isscalar(x) && isnumeric(x)));
valid_Opts2 = @(x) (isscalar(x) && (x == -1 || x>0));
assert(valid_Opts1(opts.preEvent_value), errID, errMsg);
assert(valid_Opts1(opts.postEvent_value), errID, errMsg);
assert(valid_Opts2(opts.timeWindow), errID, 'Input must be a scalar positive numeric value or "-1".');
% Get dimension names:
dims = metaData.dim_names;
% Validate if data has the following dimension names "E" and "T":
errID = 'Umitoolbox:calculate_respose_amplitude:WrongInput';
errMsg = 'Input Data is invalid. It must have an "Event" and a "Time" dimensions.';
assert(all(ismember({'E', 'T'}, dims)), errID, errMsg);

% Identify "T" dimension and permute data so Time is the first dimension:
idxT = find(strcmp('T', dims));
orig_dim_indx = 1:numel(dims);
new_dim_indx = [idxT setdiff(orig_dim_indx, idxT)];
data = permute(data, new_dim_indx);

% Store data size:
data_sz = size(data);

% Reshape data:
data = reshape(data,data_sz(1), []);

% Perform amplitude calculation:
trigFrame = round(metaData.preEventTime_sec * metaData.Freq);
switch opts.timeWindow
    case -1
        endFrame  = size(data,1);
    otherwise
        endFrame  = trigFrame + round(metaData.Freq * opts.timeWindow);
        if endFrame > size(data,1)
            endFrame = size(data,1);
        end
end

% Get baseline (pre_trigger) and postTrigger data:
bsln = data(1:trigFrame-1,:);
postTrig = data(trigFrame:endFrame,:);

% Use aggregate function OR value defined by User:
bsln = applyAggFcn(bsln, opts.preEvent_value);
postTrig = applyAggFcn(postTrig, opts.postEvent_value);

% calculate amplitude:
outData = postTrig - bsln;

% Reshape outData to match data size:
new_sz = data_sz;
new_sz(1) = 1;
outData = reshape(outData, new_sz);

% Find singleton dimensions:
singletonDims = size(outData) == 1;
% Permute outData to original data size and remove singleton dimensions:
outData = squeeze(permute(outData,[2:numel(dims) 1])); 
% Do the same in dimension names:
singletonDims = singletonDims([2:numel(dims) 1]);
new_dim_names = dims(~singletonDims);

% Create new metaData:
metaData = genMetaData(outData, new_dim_names, metaData);

end

% Local function:
function out = applyAggFcn(vals, aggfcn)
% APPLYAGGFCN performs the aggregate function of name "fcn_name" on the 1st
% dimension of the data "vals". All aggregate functions EXCLUDE NaNs!

switch aggfcn
    case 'mean'
        out = nanmean(vals, 1);
    case 'median'
        out = median(vals, 1, 'omitnan');
%     case 'mode'
%         out = mode(vals, 1);
%     case 'std'
%         out = std(vals, 0, 1, 'omitnan');
    case 'max'
        out = max(vals, [], 1, 'omitnan');
    case 'min'
        out = min(vals, [], 1, 'omitnan');
%     case 'sum'
%         out = sum(vals, 1, 'omitnan');
    otherwise
        out = repmat(aggfcn,size(vals));
end
end
