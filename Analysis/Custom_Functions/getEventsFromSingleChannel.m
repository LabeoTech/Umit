function getEventsFromSingleChannel(object, SaveFolder, varargin)
% GET_EVENTSFROMSINGLECHANNEL creates an event .MAT file in SAVEFOLDER
% from Analog signals recorded from ONE CHANNEL of LabeoTech Imaging systems.
% Inputs:
%   RawFolder: directory containing ai_xxxx.bin files.
%   SaveFolder: directory to save .MAT eventsfile.
%       optional parameters:
%       opts.Channel = Analog channel index with the TTL signal. If not
%       provided, the function will try to find the one with the largest Standard
%       Deviation value.
%       opts.threshold = Threshold  to be used in the signal detection.
%       Default = 2.5 (v).
% Output:
%   TTL_events.mat file containing channel ID, state and timestamps.
%   For details, see function CREATE_TTL_EVENTSFILE.m.
%
%%% Arguments parsing and validation %%%
p = inputParser;
% The input of the function must be a File , RawFolder or SaveFolder
% Imaging Object:
addRequired(p, 'object', @(x) isa(x,'Modality'));
% Save folder:
addRequired(p, 'SaveFolder', @isfolder);
% Optional Parameters:
% opts structure
default_opts = struct('channel', -1, 'threshold', 2.5);
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(fieldnames(x)));
% Parse inputs:
parse(p,object, SaveFolder, varargin{:});
% Initialize Variables:
object = p.Results.object;
SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;
%%%%
% Check if MetaDataFile exists:
if ~isfile(object.MetaDataFile)
    errID = 'MATLAB:UMIToolbox:getEventsFromSingleChannel:FileNotFound';
    errMsg = ['MetaDataFile not found in ' object.RawFolder];
    errMsg = strrep(errMsg, '\' , '\\');
    error(errID, errMsg);
end
disp('Reading Analog inputs...')
txt = fileread(fullfile(object.RawFolder, 'info.txt'));
sr = regexp(txt, '(?<=AISampleRate:\s*)\d+', 'match', 'once'); sr = str2double(sr);
tAIChan = regexp(txt, '(?<=AINChannels:\s*)\d+', 'match', 'once'); tAIChan = str2double(tAIChan);
aiFilesList = dir(fullfile(object.RawFolder,'ai_*.bin'));
AnalogIN = [];
for ind = 1:size(aiFilesList,1)
    data = memmapfile(fullfile(aiFilesList(ind).folder, aiFilesList(ind).name), 'Offset', 5*4, 'Format', 'double', 'repeat', inf);
    tmp = data.Data;
    tmp = reshape(tmp, sr, tAIChan, []);
    tmp = permute(tmp,[1 3 2]);
    tmp = reshape(tmp,[],tAIChan);
    AnalogIN = [AnalogIN; tmp];
end
% Crop to first and last camera triggers:
camT = diff(AnalogIN(:,1) > 2.5); camT = [camT;NaN];
camTOn = find(camT == 1,1,'first');
camTOff = find(camT == -1,1,'last');
AnalogIN = AnalogIN(camTOn:camTOff,:);

if opts.channel == -1
    STDev = std(AnalogIN(:,2:end), 0, 1);% exclude Cam triggers from search.
    sigChan = find(STDev == max(STDev)) + 1;
    sigChan = sigChan(1); % if find more than one, pick the first (arbitrary choice here...)
else
    sigChan = opts.channel;
end
signal = downsample(AnalogIN(:,sigChan),100); % I did this to try to eliminate fast artifacts due to the photodiode voltage fluctuations(BrunoO 23/03/2021).
sr = sr/100;
disp('Finding events...')
% For data created in PsychToolbox:
if endsWith(object.MetaDataFile, '.mat')
    a = load(object.MetaDataFile);
elseif endsWith(object.MetaDataFile, '.csv')
    % For data created in PsychoPy:
    a = readtable(object.MetaDataFile);
    a = a(~isnan(a.SF),:);
    condList = cell(1,height(a));
    colNames = a.Properties.VariableNames;
    indx = find(strcmp(colNames,'trials_thisRepN'));
    colNames = colNames(1:indx-1);
    for i = 1:height(a)
        str = {};
        for j = 1:numel(colNames)
            str = [str,{[colNames{j} num2str(a.(colNames{j})(i))]}];
        end
        condList{i} = strjoin(str,'-');
    end
end
%
if isfield(a, 'BarSize')
    Direction = a.DriftDirection;
    nSweeps = a.nTrials;
    if Direction == -1
        eventID = repmat([0:90:270]',nSweeps*2,1); % accounts for ON/OFF state.
        eventID = sort(eventID);
    else
        eventID = repmat(Direction, nSweeps*2,1);
    end
    condList = num2cell(unique(eventID),2);
    [~, state, timestamps] = getEventsFromTTL(signal, sr, opts.threshold);
elseif isfield(a,'TrialList')
    [condList, ~, eventID] = unique(a.TrialList, 'rows');
    condList = num2cell(condList,2);
    % duplicate event ID to account for ON/OFF states:
    eventID = repelem(eventID,2);
    [~, state, timestamps] = getEventsFromTTL(signal, sr, opts.threshold);
elseif exist('condList', 'var')
    [condList, ~, eventID] = unique(condList, 'rows');    
    % duplicate event ID to account for ON/OFF states:
    eventID = repelem(eventID,2);
    [~, state, timestamps] = getEventsFromTTL(signal, sr, opts.threshold);
else
    [eventID, state, timestamps] = getEventsFromTTL(signal, sr, opts.threshold);
    condList = unique(eventID);
    condList = num2cell(condList,2);
end
% This sections is to try to remove different artifacts from
% PsychToolbox...
deltaT = diff(timestamps);
idx = deltaT < (mean(deltaT) - 1.96*std(deltaT,1)); idx = [idx;false]; % Remove every pulse outside 1,96 STD below the mean (basic statistical outlier detection).
if sum(idx) > 0
    disp('Artifact Detected and removed from Photodiode signal')
end
stateOn = find(state(~idx) == 1, 1, 'first');
stateOff = find(state(~idx) == 0, 1, 'last');
state = state(stateOn:stateOff);
timestamps = timestamps(stateOn:stateOff);
%%%
% Save to EVENTS.MAT file:
saveEventsFile(SaveFolder, eventID, timestamps, state, condList)
disp('Done.')
end