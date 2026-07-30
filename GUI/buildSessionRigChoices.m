function [labels, ids, defaultRigID] = buildSessionRigChoices(rigs)
%BUILDSESSIONRIGCHOICES Build Add Session rig dropdown choices.
%
% Keeps readable, active rigs; prepends an explicit no-rig choice; and
% selects the first available default rig.

labels = {'(No rig)'};
ids = {''};
defaultRigID = '';

requiredVariables = { ...
    'RigID', 'DisplayName', 'IsDefault', 'Status', 'IsReadable'};
if ~istable(rigs) || ...
        ~all(ismember(requiredVariables, ...
        rigs.Properties.VariableNames))
    return
end

isAvailable = rigs.IsReadable & ...
    strcmpi(string(rigs.Status), "active");
rigs = rigs(isAvailable, :);

for iRig = 1:height(rigs)
    rigID = char(rigs.RigID(iRig));
    displayName = char(rigs.DisplayName(iRig));
    if isempty(displayName) || strcmp(displayName, rigID)
        label = rigID;
    else
        label = sprintf('%s (%s)', displayName, rigID);
    end

    labels{end+1} = label; %#ok<AGROW>
    ids{end+1} = rigID; %#ok<AGROW>
    if isempty(defaultRigID) && rigs.IsDefault(iRig)
        defaultRigID = rigID;
    end
end
end
