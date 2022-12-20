function getEventsFromSingleChannel(object, SaveFolder, varargin)
% GET_EVENTSFROMSINGLECHANNEL creates an event .MAT file in SAVEFOLDER
% from Analog signals recorded from ONE CHANNEL of LabeoTech Imaging
% systems. 
% 
% Important!
% This function won't work when event signals exist in more than one channel!
%
% Inputs:
%   RawFolder: directory containing ai_xxxx.bin files.
%   SaveFolder: directory to save .MAT eventsfile.
%       optional parameters:
%       opts.Channel = Analog channel index with the TTL signal. If not
%       provided, the function will try to find the one with the largest Standard
%       Deviation value.
%       opts.threshold = Threshold to be used in the signal detection.
%       Default = 2.5 (v).
% Output:
%   TTL_events.mat file containing channel ID, state and timestamps.
%   For details, see function CREATE_TTL_EVENTSFILE.m.
%

% opts structure
default_opts = struct('channel', 'auto', 'threshold', 'auto', 'TriggerType','EdgeSet');
opts_values = struct('channel', {{'auto',Inf}}, 'threshold',{{'auto',Inf}}, 'TriggerType', {{'EdgeSet', 'EdgeToggle'}});%#ok  % This is here only as a reference for PIPELINEMANAGER.m.

%%% Arguments parsing and validation %%%
p = inputParser;
% Imaging Object:
addRequired(p, 'object', @(x) isa(x,'Modality'));
% Save folder:
addRequired(p, 'SaveFolder', @isfolder);
% Optional Parameters:

% Notes on TriggerType:
% Two possible modes:
% 1 - "EdgeSet": The signal stays on for the duration of the trial.
% 2 - "EdgeToggle": The signal turns on and off at the beginning and at the
% end of the trial.
% h = 
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
a = fileread(fullfile(object.RawFolder, 'info.txt'));
sr = regexp(a, '(?<=AISampleRate:\s*)\d+', 'match', 'once'); sr = str2double(sr);
tAIChan = regexp(a, '(?<=AINChannels:\s*)\d+', 'match', 'once'); tAIChan = str2double(tAIChan);
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

if strcmp(opts.channel, 'auto')
    STDev = std(AnalogIN(:,2:end), 0, 1);% exclude Cam triggers from search.
    sigChan = find(STDev == max(STDev)) + 1;
    sigChan = sigChan(1); % if find more than one, pick the first (arbitrary choice here...)
else
    sigChan = opts.channel;
end
f = fdesign.lowpass('N,F3dB', 4, 30, sr); % Apply low-pass filter @30Hz to remove high-frequency noise.
lpass = design(f,'butter');
signal = filtfilt(lpass.sosMatrix, lpass.ScaleValues, AnalogIN(:,sigChan)')';
% 
if strcmp(opts.threshold, 'auto')
    opts.threshold =  min(signal) + ((max(signal) - min(signal))/2);
  fprintf('Automatic threshold selected!\n\tMin signal: %0.2f\n\tMax signal: %0.2f\n\tThreshold: %0.2f\n',...
      [min(signal), max(signal), opts.threshold])
else
    opts.threshold = str2double(opts.threshold);
end
disp('Finding events...')
% Loading meta data file:
if endsWith(object.MetaDataFile, '.mat')
    % For data created in PsychToolbox:
    a = load(object.MetaDataFile);
elseif endsWith(object.MetaDataFile, '.csv')
    
    % For data created in PsychoPy:
    a = readtable(object.MetaDataFile);
    condList = cell(1,height(a));
    if isfield(a,'SF')
        % For Veronique:
        a = a(~isnan(a.SF),:);        
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
    elseif isfield(a, 'Contrast')
        % For Catherine:
        % Here the columns of interest are : "AngDir" and "Contrast". DO
        % NOT CHANGE THEM!
        a = a(~isnan(a.Contrast),:);        
        condList = arrayfun(@(x,y) ['Dir-', num2str(x), 'Cont-',num2str(y)], ...
            a.AngDir, a.Contrast, 'UniformOutput',false);        
    end
    
elseif endsWith(object.MetaDataFile, '.txt')
    % For data created in DMD system
    disp('CREATING SECTION ...');
    a = fileread(object.MetaDataFile);
    % Get Events Description table:
    idx_start = regexpi(a, 'ID\t\LP\t\Duration\tImage');
    a(idx_start:end);
    % Create event Name list:
    out = textscan(a(idx_start:end), '%s %s %s %s', 'Delimiter', '\t', ...
        'HeaderLines', 1); 
    out = horzcat(out{:});
    evntNameList = cell(size(out,1),1);    
    for i = 1:size(out,1)
        evntNameList{i} = ['LP_', out{i,2}, '-Dur_' out{i,3}, '-Img_' out{i,4}];
    end
    % Get list of event order:
    evntOrd = regexpi(a, '(?<=Events Order:).*', 'match','once','dotexceptnewline');
    evntOrd = strip(evntOrd);
    evntOrd = str2num(evntOrd);    
end
% Here, we identify the type of recording based on the variables in the
% meta data file (from PsychToolbox scripts!)
if isfield(a, 'BarSize')
    Direction = a.DriftDirection;
    nSweeps = a.nTrials;
    if Direction == -1
        eventID = repmat([1:4]',nSweeps*2,1); % accounts for ON/OFF state.
        eventID = sort(eventID);
    else
        eventID = repmat(Direction, nSweeps*2,1);
    end
    condList = {'0','90','180','270'};
    [~, state, timestamps] = getEventsFromTTL(signal, sr, opts.threshold, opts.TriggerType);
elseif isfield(a,'TrialList')
    [condList, ~, eventID] = unique(a.TrialList, 'rows');
    condList = num2cell(condList,2);
    % duplicate event ID to account for ON/OFF states:
    eventID = repelem(eventID,2);
    [~, state, timestamps] = getEventsFromTTL(signal, sr, opts.threshold, opts.TriggerType);
elseif exist('condList', 'var')
    [condList, ~, eventID] = unique(condList, 'rows');    
    % duplicate event ID to account for ON/OFF states:
    eventID = repelem(eventID,2);
    [~, state, timestamps] = getEventsFromTTL(signal, sr, opts.threshold, opts.TriggerType);
elseif exist('evntOrd', 'var')
    % For DMD system data
    [~, state, timestamps] = getEventsFromTTL(signal, sr, opts.threshold, opts.TriggerType);
    eventID = evntOrd;
    condList = evntNameList;
else
    [eventID, state, timestamps] = getEventsFromTTL(signal, sr, opts.threshold, opts.TriggerType);
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