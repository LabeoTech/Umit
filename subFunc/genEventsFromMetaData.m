function genEventsFromMetaData(metaData, saveFolder)
% GENEVENTSFROMMETADATA is a helper function for run_ImagesClassification
% and other data import functions for data obtained with LabeoTech Optical
% Imaging Systems.
% This function creates an "events.mat" file based on the STIM info stored in the
% data's meta data .MAT file under any variable starting with "Stim". If,
% no Stimulation is found in those variables, the function will try to use
% the trial timestamps stored in "trialID" to split the data into trials,
% with a generic trigger timestamp. The "trialID" is created by the
% function "mergeRecordings.m" to merge separate trials into one recording.
% Inputs:
%   metaData (struct): structure with the meta data of a given file 
%     from which an "events.mat" file will be created.
%   saveFolder (char): Fullpath of the folder where to save the events.mat
%   file.


% Create events.mat from StimParameters.mat file:
disp('Creating events file...');
ID = [];
timestamps = [];
state = [];
fn = fieldnames(metaData);
stim_fn = startsWith(fn, 'Stim_');
% Try to find stim source:
if ~any(stim_fn)
    stim_fn = strcmpi(fn,'stim');
end
if isfield(metaData, 'TrialID') && ( ~any(stim_fn) || isempty(metaData.Stim) || sum(metaData.Stim) == 0 )
    stimVar = 'TrialID';
elseif any(stim_fn)
    stimVar = fn(stim_fn);
else
    disp('No Stimulation info found in meta data!')
    return
end

% Decide split method
if strcmpi(stimVar, 'trialid')
    splitMethod = 'trial';
elseif numel(stimVar) == 1
    stimVar = stimVar{:};
    splitMethod = 'singleStim';
else
    splitMethod = 'multiStim';
end


switch splitMethod
    case 'trial'
        trialID_list = unique(metaData.(stimVar));
        for i = 1:length(trialID_list)
            timestamps = [timestamps; find(metaData.(stimVar) == trialID_list(i), 1, 'first')/metaData.Freq;...
                find(metaData.(stimVar) == trialID_list(i), 1, 'last')/metaData.Freq];
            state = [state;true;false];
            ID = [ID;1;1]; % Add generic ID = 1;
        end        
    case 'singleStim'        
        on_indx = find(metaData.(stimVar)(1:end-1)<.5 & metaData.(stimVar)(2:end)>.5);
        off_indx = find(metaData.(stimVar)(1:end-1)>.5 & metaData.(stimVar)(2:end)<.5);
        timestamps =[timestamps; (sort([on_indx;off_indx]))./metaData.Freq];
        state =[state; repmat([true;false], numel(on_indx),1)];
        ID = [ID; ones(size(timestamps))];
        % Rearrange arrays by chronological order:
        [timestamps,idxTime] = sort(timestamps);
        state = state(idxTime);
        ID = ID(idxTime);
    otherwise
        disp('Code to be created...')
end
eventNameList = {'1'}; % for now...
% Save event File:
saveEventsFile(saveFolder, ID, timestamps, state, eventNameList)
end


