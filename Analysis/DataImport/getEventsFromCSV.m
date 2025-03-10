function getEventsFromCSV(RawFolder,SaveFolder,varargin)
% GETEVENTSFROMCSV generates an "events.mat" file from events information stored in a .CSV file.
% This function reads a .CSV file from the specified RawFolder, which contains events information
% with columns 'timestamps', 'state', and 'eventID'. It processes the data and saves it into an
% "events.mat" file in the specified SaveFolder. The function assumes there is a single .CSV file
% in the RawFolder unless a specific csvFile name is provided.
%
% Inputs:
%    RawFolder - Path to the folder containing the .CSV file with events information.
%    SaveFolder - Path to the folder where the "events.mat" file will be saved.
%    csvFile (optional) - Name of the .CSV file. If not provided,
%                           the function will automatically search for a .CSV file in the RawFolder. 
%                           Important: If there are multiple .CSV files in the folder,
%                                      the function will pick the first one found.
%

% Parse inputs
p = inputParser;
addRequired(p, 'RawFolder', @isfolder);
addRequired(p, 'SaveFolder', @isfolder);
addOptional(p, 'csvFile','auto',@ischar);
parse(p, RawFolder, SaveFolder,varargin{:});

if strcmp(p.Results.csvFile, 'auto')
    csvTag = '*.csv';
else
    [~,fname,~] = fileparts(p.Results.csvFile);
    csvTag = [fname '.csv'];
end
% Find .CSV file in folder:
csvFile = dir(fullfile(RawFolder, csvTag));
assert(~isempty(csvFile),'CSV file not found in folder "%s"!',RawFolder)    
csvFile = csvFile(1);
% Read the file:
evInfo  = readtable(fullfile(RawFolder, csvFile.name));
% Check if all columns exist:
evInfo.Properties.VariableNames = lower(evInfo.Properties.VariableNames);
assert(all(ismember(evInfo.Properties.VariableNames,{'timestamps','state','eventid'})),...
    'Invalid column names! The following columns should exist: "timestamps", "state" and "eventID')
% Create event indices from "eventID" column:
eventNameList = unique(evInfo.eventid);
eventID = zeros(height(evInfo),1,'uint16');
for ii = 1:height(evInfo)
    eventID(ii) = find(strcmp(eventNameList,evInfo.eventid(ii)));
end
state = uint8(evInfo.state);
timestamps = single(evInfo.timestamps);
% Save data to events.mat file:
saveEventsFile(SaveFolder,eventID,timestamps,state,eventNameList);
end
