function outData = genAmplitudeMaps(data,SaveFolder,varargin)
% GENAMPLITUDEMAPS computes response amplitude maps by subtracting the baseline data 
% from the response time period data. An "events.mat" file is required.
%
% Inputs:
%   data (numeric array): Numerical matrix containing image time series with dimensions [Y, X, T].
%   SaveFolder (char): Folder containing "data" and the associated "AcqInfos.mat" file.
%   opts (struct, optional): Structure containing extra parameters. Defaults:
%       - baselineMeasure (char): Measure for baseline ('mean', 'median', 'min', 'max').
%       - responseMeasure (char): Measure for response ('mean', 'median', 'min', 'max').
%       - TimeWindow_sec (numeric array or char): Time window for response ('all' or [start, end]).
%
% Outputs:
%   outData (struct): Structure containing the amplitude maps and event-related information.
%
% Note:
%   - The "events.mat" file must be present in the SaveFolder.


% Defaults: 
default_Output = 'amplitude_Map.dat'; %#ok This line is here just for Pipeline management.
default_opts = struct('baselineMeasure', 'median', 'responseMeasure', 'max', 'TimeWindow_sec', 'all');
opts_values = struct('baselineMeasure', {{'mean', 'median', 'min','max'}}, 'responseMeasure',{{'mean', 'median', 'min','max'}},'TimeWindow_sec',{{'all',[0 Inf]}});% This is here only as a reference for PIPELINEMANAGER.m. 

%%% Arguments parsing and validation %%%

% Parse inputs:
p = inputParser;
addRequired(p,'data',@(x) isnumeric(x) && ndims(x) == 3); % Validate if the input is numerical or a structure.
addRequired(p,'SaveFolder',@(x) isfolder(x)); % Validate if the SaveFolder exists.
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
% Parse inputs:
parse(p,data,SaveFolder,varargin{:});
opts = p.Results.opts;
clear p
% Further validation of optional arguments:
errID = 'Umitoolbox:calculate_respose_amplitude:InvalidDataType';
errMsg = 'Invalid Input. Input must be member of {"mean", "median", "min","max"} or a numeric scalar';
valid_Opts1 = @(x) (ismember(char(x), opts_values.preEvent_value) || (isscalar(x) && isnumeric(x)));
valid_Opts2 = @(x) ( (~isempty(x) && isscalar(x) && x > 0) )  || strcmp(x, 'all');
assert(valid_Opts1(opts.preEvent_value), errID, errMsg);
assert(valid_Opts1(opts.postEvent_value), errID, errMsg);
assert(valid_Opts2(opts.postEvnt_timeWindow_sec), errID, 'Input must be a scalar positive numeric value or "all".');
% Validate if the "events.mat" file exists in the save folder:
assert(isfile(fullfile(SaveFolder,"events.mat")), 'The "events.mat" file wasn''t found in the Save Folder.');

% Split data by Events:
ev = EventsManager(SaveFolder);
data = ev.splitDataByEvents(data);
evInfo = ev.exportEventInfo;
% Get frames to analyse:
baselineFrames = 1:round(evInfo.baselinePeriod*evInfo.FrameRateHz); % Baseline
if strcmpi(opts.TimeWindow_sec,'all')
    frOn = baselineFrames(end)+1;
    frOff = size(data,4);
else       
    % Get frames values:
    frOn = round(opts.TimeWindow_sec(1)*evInfo.FrameRateHz) + baselineFrames(end);
    frOff = round(opts.TimeWindow_sec(2)*evInfo.FrameRateHz) + baselineFrames(end);
    % Reset to default if input values are out of range
    if any([frOn, frOff] > size(data,4)) || frOn > frOff
        warning('TimeWindow onset is out of range! Reset to default ("all")')
        frOn = baselineFrames(end) +1;
        frOff = size(data,4);
    end
end
responseFrames = frOn:frOff;
%
disp('Calculating response amplitude...')

% Store data size:
data_sz = size(data);
% Reshape data:
data = reshape(data,data_sz(1), []);
% Use aggregate function OR value defined by User:
baseline = applyAggFcn(data(baselineFrames,:), opts.baselineMeasure);
response = applyAggFcn(data(responseFrames,:), opts.responseMeasure);
% Calculate amplitude map:
amplitudeMap = response - baseline;
% Reshape amplitude maps:
amplitudeMap = reshape(amplitudeMap,data_sz([1 2 3]));
% Store data in structure with meta data:
evInfo.TimeWindow_sec = opts.TimeWindow_sec;
outData = genDataStructure(amplitudeMap,'hasEvents',true,'extraInfo',evInfo);
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
