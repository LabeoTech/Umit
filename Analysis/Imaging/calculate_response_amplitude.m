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

% Defaults: !Instantiate one default per line!
default_Output = 'amplitude_Map.dat'; %#ok This line is here just for Pipeline management.
default_opts = struct('preEvent_value', 'median', 'postEvent_value', 'max', 'TimeWindow_sec', 'all');
opts_values = struct('preEvent_value', {{'mean', 'median', 'min','max'}}, 'postEvent_value',{{'mean', 'median', 'min','max'}},'TimeWindow_sec',{{'all',Inf}});% This is here only as a reference for PIPELINEMANAGER.m. 

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
valid_Opts1 = @(x) (ismember(char(x), opts_values.preEvent_value) || (isscalar(x) && isnumeric(x)));
assert(valid_Opts1(opts.preEvent_value), errID, errMsg);
assert(valid_Opts1(opts.postEvent_value), errID, errMsg);
if numel(opts.TimeWindow_sec)~= 2 || diff(opts.TimeWindow_sec < 0)
    error('Invalid time window. The input must be a pair of increasing positive numbers.');
end
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

disp('Calculating response amplitude...')

% Store data size:
data_sz = size(data);

% Reshape data:
data = reshape(data,data_sz(1), []);

% Perform amplitude calculation:
trigFrame = round(metaData.preEventTime_sec * metaData.Freq);
if strcmpi(opts.TimeWindow_sec,'all')
    % Use all post-event time when the user selected "all":
    frOn = trigFrame + 1;
    frOff  = size(data,1);
else
    % Otherwise, select the post-event time window of choice in "opts":
    % Get frames values:
    frOn = round(opts.TimeWindow_sec(1)*metaData.Freq) + trigFrame;
    frOff = round(opts.TimeWindow_sec(2)*metaData.Freq) + trigFrame;
    
    % Reset to default if input values are out of range
    if frOn > size(data,1) || frOff > size(data,1) || frOn > frOff
        warning('TimeWindow onset is out of range! Reset to default ("all")')
        frOn = trigFrame + 1;
        frOff = size(data,1);
    end
end

% Get baseline (pre_trigger) and postTrigger data:
bsln = data(1:trigFrame,:);
postTrig = data(frOn:frOff,:);

% Use aggregate function OR value defined by User:
bsln = applyAggFcn(bsln, opts.preEvent_value);
postTrig = applyAggFcn(postTrig, opts.postEvent_value);
% Reshape data to match original data size:
new_sz = data_sz;
new_sz(1) = 1;
bsln = reshape(bsln, new_sz);
postTrig = reshape(postTrig, new_sz);
% calculate amplitude:
outData = postTrig - bsln;
% Find singleton dimensions:
singletonDims = size(outData) == 1;
% Permute outData to original data size and remove singleton dimensions:
outData = squeeze(permute(outData,[2:numel(dims) 1])); 
% Do the same in dimension names:
singletonDims = singletonDims([2:numel(dims) 1]);
new_dim_names = dims(~singletonDims);

% Create new metaData:
metaData = genMetaData(outData, new_dim_names, metaData);
disp('Done!')
end

% Local function:
function out = applyAggFcn(vals, aggfcn)
% APPLYAGGFCN performs the aggregate function of name "fcn_name" on the 1st
% dimension of the data "vals". All aggregate functions EXCLUDE NaNs!

switch aggfcn
    case 'mean'
        out = mean(vals, 1, 'omitnan');
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
