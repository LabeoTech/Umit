function [labels, ids, defaultRigID] = buildSessionRigChoices(rigs, currentRigID)
%BUILDSESSIONRIGCHOICES Build assignable Rig dropdown choices.
%
%   New assignments include only Active and Available Rigs. An Archived
%   currentRigID may be included solely to preserve a historical selection.

if nargin < 2
    currentRigID = '';
end
currentRigID = char(string(currentRigID));
labels = {};
ids = {};
defaultRigID = '';

requiredVariables = {'RigID', 'DisplayName', 'Status', 'IsReadable'};
if ~istable(rigs) || ...
        ~all(ismember(requiredVariables, ...
        rigs.Properties.VariableNames))
    return
end

isAvailable = rigs.IsReadable & ...
    ismember(lower(string(rigs.Status)), ["active", "available"]);
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
    if isempty(defaultRigID) && strcmpi(char(rigs.Status(iRig)), 'active')
        defaultRigID = rigID;
    end
end

if ~isempty(currentRigID) && ~any(strcmpi(ids, currentRigID))
    try
        historical = UMITRigStore.openByRigID(currentRigID).getRigInfo();
        if strcmpi(historical.status, 'archived')
            displayName = char(string(historical.displayName));
            if isempty(displayName) || strcmp(displayName, historical.rigID)
                label = historical.rigID;
            else
                label = sprintf('%s (%s)', displayName, historical.rigID);
            end
            ids{end+1} = historical.rigID;
            labels{end+1} = sprintf('%s [archived - historical]', label);
        end
    catch
        % Keep unresolved historical IDs out of normal selectors.
    end
end
end
