function genEventsFromStimParams(Folder)
% GENEVENTSFROMSTIMPARAMS is a helper function for run_ImagesClassification
% and other data import functions for data obtained with LabeoTech Optical
% Imaging Systems.
% This function creates an "events.mat" file based on the experiment info and
% Stim info. The former is stored in the "AcqInfos.mat" file and the latter
% in the "StimParameters.mat" file. Both files are created during the data
% import of .bin files (using ImagesClassification, for instance).
% Inputs:
%   Folder (char): Fullpath of the folder where to save the
%     "events.mat" file and where the files "AcqInfo.mat" and
%     "StimParameters.mat" are located.


% Create events.mat from StimParameters.mat file:
disp('Creating events file...');
acqInfoFile = fullfile(Folder, 'AcqInfos.mat');
stimFile = fullfile(Folder, 'StimParameters.mat');
% Check if files exist in Folder:
assert(isfile(acqInfoFile), 'umIToolbox:genEventsFromStimParams:MissingFile',...
    'Acquisition Info file is missing!')
assert(isfile(stimFile), 'umIToolbox:genEventsFromStimParams:MissingFile',...
    'Stimulation Parameters file is missing!')
% Load information:
exp_info = load(acqInfoFile);
stim_info = load(stimFile);
% Get On and off timestamps of stims:
fn = fieldnames(stim_info);
ID =[];
timestamps = [];
state = [];
for i = 1:length(fn)
    on_indx = find(stim_info.(fn{i})(1:end-1)<.5 & stim_info.(fn{i})(2:end)>.5);
    off_indx = find(stim_info.(fn{i})(1:end-1)>.5 & stim_info.(fn{i})(2:end)<.5);
    timestamps =[timestamps; (sort([on_indx;off_indx]))./exp_info.AcqInfoStream.FrameRateHz];
    state =[state; repmat([true;false], numel(on_indx),1)];
    ID = [ID; repmat(i,size(timestamps))];
end
% Rearrange arrays by chronological order:
[timestamps,idxTime] = sort(timestamps);
state = state(idxTime);
ID = ID(idxTime);
% Look for events:
if any(startsWith('event', fieldnames(exp_info.AcqInfoStream)))
    disp('Digital stimulation data found!')
    eventID = repelem(exp_info.AcqInfoStream.Events_Order,1,2);
    uniqID = unique(eventID);
    eventNameList = cell(1,numel(uniqID));
    for i = uniqID
        eventNameList{i} = exp_info.AcqInfoStream.(['Stim' num2str(i)]).name;
    end
else
    eventID = ID;
    eventNameList = fn;
end
saveEventsFile(Folder, eventID, timestamps, state, eventNameList)
end


