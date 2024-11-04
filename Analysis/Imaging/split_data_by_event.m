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
default_opts = struct('preEventTime_sec','auto', 'postEventTime_sec','auto', 'PadWith', 'mean');
opts_values = struct('preEventTime_sec', {{'auto',Inf}}, 'postEventTime_sec', {{'auto',Inf}}, 'PadWith',{{'mean','NaN',Inf}});%#ok. This is here only as a reference for PIPELINEMANAGER.m. 

%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'data',@(x) isnumeric(x) & ndims(x) == 3); % Validate if the input is a 3-D numerical matrix:
addRequired(p,'metaData', @(x) isa(x,'matlab.io.MatFile') | isstruct(x)); % MetaData associated to "data".
addRequired(p, 'SaveFolder', @isfolder);
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x) && ...
    ( any(strcmpi(x.PadWith, {'mean', 'NaN'})) || isnumeric(x.PadWith) )); % Padding options for cases where movie snippets dont have the same length.
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
% Get pre/post event times from metaData or try to find the best timing
% based on the timestamps of events:
if strcmp(opts.preEventTime_sec, 'auto') || strcmp(opts.postEventTime_sec, 'auto')
    if isfield(metaData, 'preEventTime_sec')
        % Grab info from file's metaData:        
        fprintf('Using info from data''s metadata:\n\tPre event time: %d seconds.\n\tPost event time: %d seconds.\n',...
            [metaData.preEventTime_sec, metaData.postEventTime_sec]);
        opts.preEventTime_sec = metaData.preEventTime_sec;
        opts.postEventTime_sec = metaData.postEventTime_sec;
    else        
       % Calculate pre/post times from "events" file:
       disp('Calculating from events...')       
       tmTrial= round(mean(diff(evDat.timestamps(evDat.state == 1)), 'omitnan'));
       % Use 20% of the time as pre and 80 % as post:
       opts.preEventTime_sec = round(.2*tmTrial,2);
       opts.postEventTime_sec = tmTrial - opts.preEventTime_sec;
       fprintf('Pre and post-event times calculated from "events" file:\n\tPre event time: %0.2f seconds.\n\tPost event time: %0.2f seconds.\n',...
            [opts.preEventTime_sec, opts.postEventTime_sec]);
    end
% else
%     % Given that the variables of pre and post event times accept the
%     % string "auto" as input, numbers is also interpreted as strings and 
%     % need to be transformed back to numerical data types:
%     opts.preEventTime_sec = str2double(opts.preEventTime_sec);
%     opts.postEventTime_sec = str2double(opts.postEventTime_sec);
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
disp('Splitting data by events...')
% Fill empty matrix with data segments
fix_snippet = false;
for i = 1:n_trial
    trialFr = round(sr*timestamps(i))+1;
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
disp('Done!');
% Add variables to metaData.
extraParams = metaData;
extraParams.preEventTime_sec = opts.preEventTime_sec;
extraParams.preEventTime_sec = opts.preEventTime_sec;
extraParams.postEventTime_sec = opts.postEventTime_sec;
extraParams.eventID = evDat.eventID(evDat.state == 1);
extraParams.eventNameList = evDat.eventNameList;
metaData = genMetaData(outData, new_dims, extraParams);
end