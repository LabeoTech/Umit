function mergeDirs_in_eventFile(SaveFolder)
% MERGEDIRS_IN_EVENTFILE modifies the variables in events.mat file 
% by merging the grating directions as "repetitions" of the same stimulus.
% This is a custom function used in Matthieu Vanni's lab.

% Inputs:
%   SaveFolder: directory to save .MAT eventsfile.

%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p, 'SaveFolder', @isfolder);
% Parse inputs:
parse(p, SaveFolder);
%Initialize Variables:
SaveFolder = p.Results.SaveFolder;
clear p
%%%%
disp('Merging Directions...')
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
% Get Directions from strings:
evntNums = cellfun(@(x) str2num(x), mD.eventNameList, 'UniformOutput',0);
if isempty(evntNums{1})
    % FOR CATHERINE
    list = mD.eventNameList;
    angDirs = regexpi(list, 'AngDir-\d+\-(?=Cont)','match');angDirs = [angDirs{:}]';    
    [new_eventNameList, ~,newIDs] = unique(erase(list,angDirs),'stable');    
else
    % FOR EVERYBODY ELSE
    evntNums = cell2mat(evntNums);
    uniq_gratings = unique(evntNums(:,1:3), 'rows');
    newIDs = zeros(size(evntNums,1),1);
    for i = 1:size(uniq_gratings,1)
        idx = ismember(evntNums(:,1:3), uniq_gratings(i,:),'rows');
        newIDs(idx) = i;
    end
    % Recreate eventNameList
    new_eventNameList = num2cell(uniq_gratings,2);
    % Transform items back to strings:
    new_eventNameList = cellfun(@(x) num2str(x), new_eventNameList, 'UniformOutput',0);
end
newID_list = zeros(size(evntID));
for i = 1:length(evntID)
    idx = evntID == evntID(i);    
    newID_list(idx) = newIDs(evntID(i));
end
% Change "events.mat" file:
mD.Properties.Writable = true;
mD.eventID = repelem(newID_list,2); % Replace IDs for state == 0
mD.eventNameList = new_eventNameList;
disp('Events file updated!')
end