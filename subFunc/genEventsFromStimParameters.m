function genEventsFromStimParameters(saveFolder, varargin)
% GENEVENTSFROMSTIMPARAMETERS is a helper function for run_ImagesClassification
% for data obtained with LabeoTech Optical Imaging Systems.
% This function creates an "events.mat" file based on the STIM info stored in the
% "StimParameters.mat" file.

% Inputs:
%   saveFolder (char): Fullpath of the folder where to save the events.mat
%   file.
%   stimStruct (struct, optional): Structure containing the "Stim"
%   variables with names starting with "Stim_". If provided, this structure
%   will be used to generate the "events.mat" file instead of the
%   "StimParameters.mat" file.

if nargin == 2
    stimFile = varargin{:};
else
    % Create events.mat from StimParameters.mat file:
    stimFile = load(fullfile(saveFolder,'StimParameters.mat'));
end
if isempty(stimFile)
    warning('Aborted event file creation. StimParameters file not found in save folder!')
    return
elseif ~any(startsWith(fieldnames(stimFile), 'stim', 'IgnoreCase',true))
    warning('Aborted event file creation. Triggers not found in StimParameters file!')
    return
end

fn = fieldnames(stimFile);
stimVar = fn(startsWith(fn, 'stim_', 'IgnoreCase', true));
disp('Creating events file...');
% Find stim source:
eventID = [];
state = [];
timestamps = [];
% eventNameList = {};
for i = 1:length(stimVar)
    [ID,chanState,tmstmp] = getEventFromChannel(stimFile.(stimVar{i}), stimFile.Freq);
    eventID = [eventID;ID.*i];
    state = [state;chanState];
    timestamps = [timestamps;tmstmp];    
end
% Reorganize arrays chronologically:
[timestamps, order] = sort(timestamps);
eventID = eventID(order);
state = state(order);
% Save event File:
saveEventsFile(saveFolder, eventID, timestamps, state, num2cell(unique(eventID)))
end

% Local function
function [ID,state,timestamps] = getEventFromChannel(data, FrameRateHz)
ID = [];state = [];timestamps = [];
id_list = unique(data(:)); id_list(id_list == 0) = [];
for i = 1:length(id_list)
    on_indx = find(data(1:end-1)<.5 & data(2:end)>.5 & data(:,2:end) == id_list(i));
    off_indx = find(data(1:end-1)>.5 & data(2:end)<.5 & data(:,1:end-1) == id_list(i));
    timestamps =[timestamps; (sort([on_indx,off_indx]))'./FrameRateHz];
    state =[state; repmat([true;false], numel(on_indx),1)];
    ID = [ID; repmat(id_list(i),2,1)];
end
% Rearrange arrays by chronological order:
[timestamps,idxTime] = sort(timestamps);
state = state(idxTime);
ID = ID(idxTime);
end
