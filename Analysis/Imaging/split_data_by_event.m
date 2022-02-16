function [outData, metaData] = split_data_by_event(data, metaData, SaveFolder, varargin)
% SPLIT_DATA_BY_EVENT reshapes an image time series dataset in a 4D matrix of 
% dimensions: {E,Y,X,T}%
% Inputs:
%   data: 3D numerical matrix with dimensions: {Y,X,T}
%   metaData: .mat file with meta data associated with "data".
%   opts (optional) : structure containing extra parameters.

% Outputs: 
%   outData: 4D numerical matrix with dimensions {E,Y,X,T}.   
%   metaData: .mat file with meta data associated with "outData".

% Defaults:
default_Output = 'data_splitByEvent.dat';  %#ok. This line is here just for Pipeline management.
default_opts = struct('preEventTime_sec',2, 'postEventTime_sec',4, 'PadWith', 'mean');
%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'data',@(x) isnumeric(x) & ndims(x) == 3); % Validate if the input is a 3-D numerical matrix:
addRequired(p,'metaData', @(x) isa(x,'matlab.io.MatFile') | isstruct(x)); % MetaData associated to "data".
addRequired(p, 'SaveFolder', @isfolder);
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x) && ...
    ismember(x.PadWith, {'mean', 'NaN'}) || isnumeric(x)); % Padding options for cases where movie snippets dont have the same length.
% Parse inputs:
parse(p,data, metaData, SaveFolder, varargin{:});
%Initialize Variables:
data = p.Results.data; 
metaData = p.Results.metaData;
folder = p.Results.SaveFolder;
opts = p.Results.opts;
clear p
%%%%

% Load "events.mat" file:
evFile = fullfile(folder, 'events.mat');
if ~isfile(evFile)
    errID = 'umIToolbox:split_data_by_event:FileNotFound';
    errMsg = ['Event file ("events.mat") not found in ' folder]; 
    errMsg = strrep(errMsg, filesep, [filesep filesep]);
    error(errID, errMsg);
else
    evDat = load(fullfile(folder, 'events.mat'));
end
%%%%%%%%%%%
szdat = size(data);
sr = metaData.Freq;
preFr = round(sr*opts.preEventTime_sec);
postFr = round(sr*opts.postEventTime_sec);
centralFr = preFr + 1;
len_trial = preFr + postFr;
n_trial = sum(evDat.state == 1);
timestamps = evDat.timestamps(evDat.state == 1);
new_dims = {'E', 'Y', 'X','T'};
[~, locB]= ismember(new_dims([2,3]), metaData.dim_names);
% Create empty matrix:
outData = nan([n_trial, szdat(locB), len_trial], 'single');
% Fill empty matrix with data segments
fix_snippet = false;
for i = 1:n_trial
    trialFr = round(sr*timestamps(i));
    start = trialFr - preFr;
    stop = trialFr + postFr - 1;
    if start < 1
        start = 1;
        fix_snippet = true;
    elseif stop > szdat(3)
        stop = szdat(3);
        fix_snippet = true;
    end
    snippet = data(:,:,start:stop);
    startFr = centralFr - (trialFr - start);
    stopFr = centralFr + (stop - trialFr);
    if fix_snippet
        warning(['Snippet size is out of bounds from Data.'...
            ' Missing data points will be replaced with ' opts.PadWith]);
        switch opts.PadWith
            case 'mean'
                avg = mean(snippet, 3, 'omitnan');
                avg = repmat(avg,1,1,size(outData,4));
                outData(i,:,:,:) = avg;
            case 'NaN'
                % empty
            otherwise
                outData(i,:,:,:) = opts.PadWith;
        end
    end
    outData(i,:,:,startFr:stopFr) = snippet;
end
% Add variables to metaData.
extraParams = metaData;
extraParams.preEventTime_sec = opts.preEventTime_sec;
extraParams.preEventTime_sec = opts.preEventTime_sec;
extraParams.postEventTime_sec = opts.postEventTime_sec;
extraParams.eventID = evDat.eventID(evDat.state == 1);
extraParams.eventNameList = evDat.eventNameList;
metaData = genMetaData(outData, new_dims, extraParams);
end