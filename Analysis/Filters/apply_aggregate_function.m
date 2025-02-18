function outData = apply_aggregate_function(data, SaveFolder, varargin)
% APPLY_AGGREGATE_FUNCTION applies an aggregate function to one
% dimensions of a .DAT file. Works with image-time-series only!

% Inputs:
%   data: numerical matrix containing image-time-series data.
%   SaveFolder (char): folder containing "data" and the associated AcqInfos.mat file.
%   opts (optional) : structure containing the function's parameters:
%       aggregateFcn (default = "mean") : name of the aggregate function.
%       dimension (default = "Time") : name of the dimension to perform
%       the calculation:
% Output:
%   outData: structure containing stats-ready aggregated data.

% Defaults:
default_Output = ' dataAgg.dat'; %#ok This is here for PIPELINEMANAGER.M.
default_opts = struct('aggregateFcn', 'mean', 'dimension', 'Time');
opts_values = struct('aggregateFcn', {{'mean', 'max', 'min', 'median', 'mode', 'sum', 'std'}}, 'dimension',{{'Space','Time','Events'}});% This is here only as a reference for PIPELINEMANAGER.m.
%%% Arguments parsing and validation %%%
% Validator for data from structure:
validateDataStructure = @(x) isstruct(x) && isDatStat(x) && isDataImageTimeSeries(x);
% Parse inputs:
p = inputParser;
addRequired(p,'data',@(x) (isnumeric(x) && ndims(x) == 3) || validateDataStructure(x)); % Validate if the input is numerical or a structure
addRequired(p,'SaveFolder',@(x) isfolder(x)); % Validate if the SaveFolder exists.
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x) && ...
    ismember(x.aggregateFcn, opts_values.aggregateFcn) && ...
    ismember(x.dimension, opts_values.dimension));
% Parse inputs:
parse(p,data, SaveFolder, varargin{:});
%Initialize Variables:
opts = p.Results.opts;
clear p
%
if isstruct(data)
    % For data in structure, build new structure and store aggregated
    % values:
    dataAgg = struct();
    fn = fieldnames(data.data);
    evList = [];
    if data.b_hasEvents
        evList = data.eventID;
    end    
    for ii =1:length(fn)
        [dataAgg.(fn{ii}), newEvList] = applyAggFcn(data.data.(fn{ii}), opts.aggregateFcn,opts.dimension,data.b_hasEvents,evList);
    end
    if ~isempty(newEvList)
        data.eventID = newEvList;
    end
    % Create new stats-ready structure:
    outData = genDataStructure(dataAgg,data,obsID,'hasEvents',data.b_hasEvents,'extraInfo',data);
else
    % For 3D arrays:
    dataAgg = applyAggFcn(data, opts.aggregateFcn, opts.dimension, false,[]);
    % Create new stats-ready structure:
    % Get Frame Rate from Folder:
    info = load(fullfile(SaveFolder,'AcqInfos.mat'));
    extraInfo.FrameRateHz = info.AcqInfoStream.FrameRateHz;
    outData = genDataStructure(dataAgg, data,'extraInfo',extraInfo);
end

end

% Local functions:

function [out,new_eventID] = applyAggFcn(vals,aggFcn,dimName,b_hasEvents,eventID)
% applyAggFcn performs aggregation on the specified dimension of the input data.
%
% This function aggregates the input array `vals` based on the specified
% dimension `dimName` ('SPACE', 'TIME', or 'EVENTS') using the provided
% aggregation function `aggFcn`. If the data is split by events, it permutes
% the array before aggregation. The function also returns a new event ID
% array `new_eventID` when aggregating across events.
%
% Input:
%   vals        - Input array to be aggregated.
%   aggFcn      - Aggregation function handle (e.g., @mean, @sum).
%   dimName     - Name of the dimension to aggregate ('SPACE', 'TIME', 'EVENTS').
%   b_hasEvents - Boolean flag indicating if the data is split by events.
%   eventID     - Event ID array for aggregating across events.

% Permute array, if the data is split by events:
if b_hasEvents
    vals = permute(vals,[2 3 4 1]);
end
new_eventID = [];
% Perform aggregation of selected dimension:
switch upper(dimName)
    case 'SPACE'
        % Aggregates in X and Y dimensions:
        out = calcAgg(vals,aggFcn,[1 2]);
    case 'TIME'
        % Aggregates in time dimension:
        out = calcAgg(vals,aggFcn,3);
    case 'EVENTS'
        % Aggregates across events:
        new_eventID = unique(eventID,'stable');
        out = zeros([size(vals,1) size(vals,2), size(vals,3), length(new_eventID)],'single');
        for ii = 1:length(new_eventID)
            idxID = eventID == new_eventID(ii);
            out(:,:,:,ii) = calcAgg(vals(:,:,:,idxID),aggFcn,4);
        end
    otherwise
        error('Unknown dimension identifier!')
end
%Permute array back to original:
if b_hasEvents
    out = ipermute(out,[2 3 4 1]);
end
% Remove singleton dimensions:
out = squeeze(out);
end

function out = calcAgg(vals,aggfcn,idxDim)
% Aggregate the values in the selected dimension.

switch aggfcn
    case 'mean'
        out = mean(vals, idxDim,'omitnan');
    case 'median'
        out = median(vals, idxDim, 'omitnan');
    case 'mode'
        out = mode(vals, idxDim);
    case 'std'
        out = std(vals, 0, idxDim, 'omitnan');
    case 'max'
        out = max(vals, [], idxDim, 'omitnan');
    case 'min'
        out = min(vals, [], idxDim, 'omitnan');
    case 'sum'
        out = sum(vals, idxDim, 'omitnan');
    otherwise
        out = vals;
end
end
