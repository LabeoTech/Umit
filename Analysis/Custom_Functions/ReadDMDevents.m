function ReadDMDevents(object,SaveFolder, varargin)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% THIS FUNCTION IS NOW OBSOLETE AND WILL BE EVENTUALLY REPLACED BY 
% THE "EVENTSMANAGER" CLASS!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% READDMDEVENTS creates an event .MAT file in SAVEFOLDER
% from Analog signals using the DMD system from LabeoTech.
% 
% Inputs:
%   RawFolder: directory containing ai_xxxx.bin files.
%   SaveFolder: directory to save .MAT eventsfile.
%   opts (struct): structure containg optional parameters.

%  Defaults:
default_opts = struct('channel', 3, 'threshold', 'auto', 'TriggerType','EdgeSet');
opts_values = struct('channel', [2:11], 'threshold',{{'auto', Inf}}, 'TriggerType', {{'EdgeSet', 'EdgeToggle'}});%#ok  % This is here only as a reference for PIPELINEMANAGER.m.

% Notes on TriggerType:
% Two possible modes:
% 1 - "EdgeSet": The signal stays on for the duration of the trial.
% 2 - "EdgeToggle": The signal turns on and off at the beginning and at the
% end of the trial.

%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p, 'object', @(x) isa(x,'Modality'));
% addRequired(p, 'RawFolder', @isfolder);
addRequired(p, 'SaveFolder', @isfolder);
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(fieldnames(x)));
% Parse inputs:
parse(p,object,SaveFolder, varargin{:});
% Initialize Variables:
object = p.Results.object;
% RawFolder = p.Results.RawFolder;
SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;
%%%%
% Look for the DMD meta data file:

if ~isfile(object.MetaDataFile) || ~endsWith(object.MetaDataFile, 'OptoGen2DParams.txt')
    error('umIToolbox:ReadDMDevents:FileNotFound', ...
        'Failed to read DMD events! The file "OptoGen2DParams.txt" was not found in raw folder.')        
end
% Read DMD file:
txt = readcell(object.MetaDataFile, 'Delimiter',':', 'HeaderLines',1);
eventID = repelem(str2num(txt{strcmp('Events Order', txt(:,1)),2}),1,2);
idx_trigType = strcmpi('external trigger', txt(:,1));
% read event description:
indx_header = find(startsWith(txt(:,1), 'ID'));
evntDesc = regexp(txt(indx_header+1:end,1), '\t', 'split');
condList = cell(size(evntDesc));
for i = 1:numel(evntDesc)
    condList{i} = ['LP_',evntDesc{i}{2},'-Dur_', evntDesc{i}{3}, '-Img_', evntDesc{i}{4}];
end

if txt{idx_trigType,2} == 0 && opts.channel ~=3    
    warning('External channel was ignored because trigger was internal.')
    opts.channel = 3;
end

% Read analog IN:
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
% Try to guess channel if it was external while the user entry was "3"
% (internal):
if txt{idx_trigType,2} == 1 && opts.channel ==3    
    warning('Internal channel was ignored because trigger was external. Trying to find which one...')
    disp('Trying to guess channel with external triggers...');
    [~,idx] = max(std(AnalogIN,0,1));
    disp(['Found triggers on channel #' num2str(idx)]);
    opts.channel = idx;
end

signal = AnalogIN(:,opts.channel);
%
if strcmp(opts.threshold, 'auto')
    opts.threshold =  min(signal) + ((max(signal) - min(signal))/2);
  fprintf('Automatic threshold selected!\n\tMin signal: %0.2f\n\tMax signal: %0.2f\n\tThreshold: %0.2f\n',...
      [min(signal), max(signal), opts.threshold])
else
    opts.threshold = str2double(opts.threshold);
end
disp('Finding events...')
% Extracting time stamps and state from analog signal:
[~, state, timestamps] = getEventsFromTTL(signal, sr, opts.threshold, opts.TriggerType);

% TEMPORARY FIX FOR TRIGGER ERRORS FROM PSYCHOPY:
% Here we get only the first triggers recorded by LabeoTech system and
% described in the .txt file. Triggers sent from psychopy after these will
% be ignored.
if length(eventID) < length(state)
    warning(['DMD System recorded ' num2str(length(eventID)/2) ' events but '...
        num2str(sum(state)) ' triggers were detected. Extra triggers were ignored!'])
    state = state(1:length(eventID));
    timestamps = timestamps(1:length(eventID));
end
% Save to EVENTS.MAT file:
saveEventsFile(SaveFolder, eventID, timestamps, state, condList)
disp('Done.')
end