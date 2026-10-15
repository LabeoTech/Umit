function [labels, rigUUIDs, defaultIndex] = buildLibraryRigChoices(rigs, preselectUUID)
%BUILDLIBRARYRIGCHOICES Choices for dialogs that pick one library rig.
%
%   [labels, rigUUIDs, defaultIndex] = buildLibraryRigChoices(rigs)
%   [labels, rigUUIDs, defaultIndex] = buildLibraryRigChoices(rigs, preselectUUID)
%
%   RIGS is the table returned by UMITRigStore.listRigs. Only readable
%   Active and Available rigs are offered (archived rigs cannot receive new
%   work). DEFAULTINDEX selects PRESELECTUUID when it is offered, otherwise
%   the Active rig (the rig library's only role in dialogs is this
%   pre-selection), otherwise the first choice; it is 0 when nothing can be
%   offered. Labels read "Display name (rigID)", with " - Active" appended
%   for the Active rig.

if nargin < 2
    preselectUUID = '';
end
labels = {};
rigUUIDs = {};
defaultIndex = 0;

requiredVariables = {'RigUUID', 'RigID', 'DisplayName', 'Status', 'IsReadable'};
if ~istable(rigs) || ~all(ismember(requiredVariables, rigs.Properties.VariableNames))
    return
end
rigs = rigs(rigs.IsReadable & ...
    ismember(lower(string(rigs.Status)), ["active", "available"]), :);

activeIndex = 0;
for iRig = 1:height(rigs)
    rigID = char(rigs.RigID(iRig));
    displayName = char(rigs.DisplayName(iRig));
    if isempty(displayName) || strcmp(displayName, rigID)
        label = rigID;
    else
        label = sprintf('%s (%s)', displayName, rigID);
    end
    if strcmpi(char(rigs.Status(iRig)), 'active')
        label = [label ' - Active']; %#ok<AGROW>
        activeIndex = iRig;
    end
    labels{end+1} = label; %#ok<AGROW>
    rigUUIDs{end+1} = char(rigs.RigUUID(iRig)); %#ok<AGROW>
end

if isempty(rigUUIDs)
    return
end
preselected = find(strcmpi(rigUUIDs, char(string(preselectUUID))), 1);
if ~isempty(preselected)
    defaultIndex = preselected;
elseif activeIndex > 0
    defaultIndex = activeIndex;
else
    defaultIndex = 1;
end
end
