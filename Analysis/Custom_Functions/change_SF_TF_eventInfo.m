function change_SF_TF_eventInfo(File, SaveFolder)
% CHANGE_SF_TF_EVENTINFO modifies the metadata .MAT file from data_splitByEvent.dat
% to consider different directions of the same grating as repetitions (not
% distinct conditions). 
% This is a custom function from Matthieu Vanni's
% lab.
% Inputs:
%   RawFolder: directory containing ai_xxxx.bin files.
%   SaveFolder: directory to save .MAT eventsfile.
% Output:
%   events.mat file containing channel ID, state and timestamps.
%   For details, see function CREATE_TTL_EVENTSFILE.m.
%
%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'File',@isfile)
addRequired(p, 'SaveFolder', @isfolder);
% Parse inputs:
parse(p,File, SaveFolder);
%Initialize Variables:
File = p.Results.File; 
SaveFolder = p.Results.SaveFolder;
%%%%
% Check if events.mat exists:
evntFile = fullfile(SaveFolder, 'events.mat');
if ~isfile(evntFile)
    errID = 'MATLAB:UMIToolbox:change_SF_TF_eventFile:FileNotFound';
    errMsg = ['Events file not found in ' SaveFolder];
    errMsg = strrep(errMsg, '\' , '\\');
    error(errID, errMsg);
end
mD = matfile(evntFile);
state = mD.state;
evntID = mD.eventID;
evntID = evntID(state == 1);

evntNames = mD.eventNameList;

evntNames = cell2mat(evntNames);
uniq_gratings = unique(evntNames(:,1:3), 'rows');
newIDs = zeros(size(evntNames,1),1);
for i = 1:size(uniq_gratings,1)
    idx = ismember(evntNames(:,1:3), uniq_gratings(i,:),'rows');
    newIDs(idx) = i;
end

newID_list = zeros(size(evntID));
for i = 1:length(evntID)
    idx = evntID == evntID(i);
    newID = newIDs(evntID(i));
    newID_list(idx) = newID;
end

new_eventID = newID_list;
new_eventNameList = num2cell(uniq_gratings,2);
% Change info in file metadata:
[~, metaDat] = mapDatFile(File);
metaDat.Properties.Writable = true;
metaDat.eventID = new_eventID;
metaDat.eventNameList = new_eventNameList;
end