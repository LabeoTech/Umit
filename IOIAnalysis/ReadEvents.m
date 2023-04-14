function ReadEvents(DataFolder,SaveFolder,TriggerChannelName,EventInfoFile,FileParsingMethod)
% READEVENTS reads trigger timestamps from the analog inputs (ai_xxxx.bin file(s)) of
% LabeoTech's imaging systems or from a .csv/.txt file. The event IDs are
% read from the text file using one of the available parsing methods (ex.
% VPIXX file). In cases where there are no triggers in the analog inputs,
% the timestamps are read directly from the text file.
%
% Inputs:
% 1- DataFolder (path):
%   Path contaning dataset from Labeo's system and/or the event info file.
% 2- SaveFolder (path):
%   Path where to save the events.mat file.
% 3- TriggerChannelName (char|cell):
%   Analog Input channel name containing the triggers. Leave empty if the
%   trigger timestamps are already stored in the "EventInfoFile"
%   The available channel names are the following:
%       'Internal-main' (default): Main internal channel.
%       'Internal-Aux' : Auxiliary internal channel.
%       'AI1' : External Analog channel #1.
%       'AI2' : External Analog channel #2.
%       'AI3' : External Analog channel #3.
%       'AI4' : External Analog channel #4.
%       'AI5' : External Analog channel #5.
%       'AI6' : External Analog channel #6.
%       'AI7' : External Analog channel #7.
%       'AI8' : External Analog channel #8.
%  If the triggers are located in more than one channel, give a list of
%  channels as a cell array of characters.
% 4- EventInfoFile (char):
%   Name of the file containing the event indices and names. This file must
%   be located in the "DataFolder". Leave empty if no event info is
%   available. In this case, the triggers will be read from the analog
%   inputs and the data will be treated as a single condition.
% 5- FileParsingMethod (char):
%   Method for parsing the "EventInfoFile". The default method consists of
%   The other available methods are the following:
%   a) 'default: reads the values of the columns with headers TIMESTAMP (optional),ID,NAME
%       stored in a .csv file.
%   b) 'VPIXX': reads content of the .txt or .vpixx files. Extracts the list
%       of event IDs and Names from the columns "Condition" and "Stimulus"
%       respectively. This information is extracted from the section "RAW DATA".
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Read Analog inputs:

% Get timestamps from AnalogIN channels:

if ~isempty(TriggerChannelName)
    chanNameList = {'Internal-main', 'Internal-Aux','AI1', 'AI2','AI3','AI4','AI5','AI6','AI7','AI8'}; % List of existing Analog channel names.
    [~,stimChan] = ismember(upper(TriggerChannelName), upper(chanNameList));
    if stimChan == 0
        warning('Invalid channel name! The "Internal-main" channel will be read instead.');
        stimChan = 2;
    else
        stimChan = stimChan + 1 ; % Shift channel index to skip camera channel.
    end
    AcqInfoStream = ReadInfoFile(DataFolder);
    Fields = fieldnames(AcqInfoStream); %Recovers Stimulation information from info.txt file
    idx = contains(Fields, 'stimulation','IgnoreCase',true);
    Fields = Fields(idx);
    cnt = 0;
    for indS = 1:length(Fields)
        if(isnumeric(AcqInfoStream.(Fields{indS})) && ...
                AcqInfoStream.(Fields{indS}) > 0)
            if( AcqInfoStream.Stimulation == 0 )
                AcqInfoStream.Stimulation = 1;
            end
            triggerInfo = ReadAnalogsIn(DataFolder, SaveFolder, AcqInfoStream, stimChan); %Read analog inputs.
            cnt = cnt + 1;
            break;
        end
    end
    if( cnt == 0 ) %No stimulation found in info.txt file.
        fprintf('Stimulation not detected. \n');        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(EventInfoFile) && ~exist('triggerInfo','var')
    error('Triggers not found in analog inputs and no event info file was provided! Operation Aborted!');
end
if exist('triggerInfo','var')
    if ~isempty(EventInfoFile) && ~exist(fullfile(DataFolder,EventInfoFile),'file')
        warning('Event info file does not exist. Events.mat data will be saved as a single condition recording')
    end
end
%%% Parse Event file %%%%%%%%%%%%
if ~isempty(EventInfoFile)
    switch lower(FileParsingMethod)
        case 'default'
            tab = readtable(fullfile(DataFolder,EventInfoFile));
            tab.Properties.VariableNames = upper(tab.Properties.VariableNames);
        case 'vpixx'
            tab = readVpixxFile(fullfile(DataFolder,EventInfoFile));
        otherwise
            error('Unknown file parsing method!')
    end
end
% Get triggers from the EventInfoFile, if applicable:
if ~exist('triggerInfo','var')
    if ~any(strcmpi(tab.Properties.VariableNames,'timestamp'))
        error('Missing timestamps in EventInfoFile! Operation aborted.')
    end
    timestamp = tab.TIMESTAMP;
else
    % Create timestamps (in seconds) from triggerInfo:   
    disp('creating timestamps from AnalogIN data...')
    fn = setdiff(fieldnames(triggerInfo),'Stim');
    fn = fn(startsWith(fn,'Stim','IgnoreCase',true));    
    for ii = 1:length(fn)        
        
        % TDB
        
    end
end
% Get condition IDs and names from the EventInfoFile:
if exist('tab','var')
    eventID = tab.ID;
    disp('Getting list of event names...')
    IDlist = unique(eventID);
    eventNameList = cell(size(IDlist));
    for ii = 1:length(IDlist)
        row = find(eventID == IDlist(ii),1,'first');
        eventNameList(ii) = tab.NAME(row);
    end
else
    % Check for event info in triggerInfo (Stim digital!)
    disp('Getting event info from LabeoSystem...');
    % TBD
end

% Save the info to the "events.mat" file:
saveEventsFile(SaveFolder,eventID,timestamp,[],eventNameList)
end
%%%%% Local functions for parsing event info from text files:
function out = readVpixxFile(file)
% READVPIXXFILE extracts the condition ID, name and timestamps (if applicable)
% from the RAW DATA section of a .txt or .vpixx file. 

% Extract only the part of the text containing the real sequence of
% conditions (RAW DATA section):
filetext = fileread(file);
[~,idx_start] = regexp(filetext, 'RAW DATA');
idx_stop = regexp(filetext, 'SORTED');
filetext = strip(filetext(idx_start+1:idx_stop-1));
tab = strsplit(filetext, '\n')';
% Get stimulus ID and order. Ignore Event and Time columns for now!
out = {};
myCols = [2 3]; % Keep just Condition and Stimulus.
for i = 2:length(tab)
    str = strsplit(tab{i},'\t');
    if isempty(str{1})
        % Skip Event rows.
        continue
    end
    out = [out;str(myCols)];
end

out(:,1) = cellfun(@str2double, out(:,1), 'UniformOutput',false);
out = cell2table(out);
out.Properties.VariableNames = {'ID','NAME'};
end