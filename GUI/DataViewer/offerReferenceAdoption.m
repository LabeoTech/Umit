function result = offerReferenceAdoption(parentFigure, saveFolder, opts)
%OFFERREFERENCEADOPTION Offer to add an aligned folder's reference to its subject library.
%
%   result = offerReferenceAdoption(parentFigure, saveFolder)
%   result = offerReferenceAdoption(..., 'Quiet', true)
%   result = offerReferenceAdoption(..., 'Interactive', false, 'Choice', choice)
%
%   Project redesign Phase 7 (section 4 Q7). A folder aligned without a
%   project carries its committed reference in DataParams.registration.
%   Once the folder is bound to a project session (planReferenceAdoption):
%     - the reference is already in the subject library (active or
%       archived): the match is reported only, nothing is recorded;
%     - it is not: the user may add it (adoptFolderReference), optionally
%       as the subject's active reference. The folder is never changed.
%
%   Name-Value options:
%       Quiet       - true: say nothing when there is nothing to offer (not
%                     bound, not aligned). Use after a binding. Default
%                     false (menu: every outcome is shown).
%       Interactive - false: no dialogs (tests). Default true.
%       Choice      - with Interactive false: 'add', 'addActive' or 'skip'
%                     (default 'skip').
%
%   RESULT fields: changed (a reference was added), state (the plan state,
%   or 'error'), resourceUUID, message (for the status line) and severity
%   ('info', 'success', 'warning' or 'error').
%
%   See also planReferenceAdoption, adoptFolderReference, offerStaleLockClearing.

arguments
    parentFigure
    saveFolder (1,1) string
    opts.Quiet (1,1) logical = false
    opts.Interactive (1,1) logical = true
    opts.Choice (1,1) string {mustBeMember(opts.Choice, ["add", "addActive", "skip"])} = "skip"
end

saveFolder = char(saveFolder);
result = struct('changed', false, 'state', 'error', 'resourceUUID', '', ...
    'message', '', 'severity', 'info');
title = 'Folder Reference';

try
    plan = planReferenceAdoption(saveFolder);
catch ME
    result.message = sprintf('Could not compare the folder reference with the library: %s', ...
        ME.message);
    result.severity = 'error';
    iAlert(parentFigure, opts, result.message, title, 'error');
    return
end
result.state = plan.state;
result.message = plan.message;

switch plan.state
    case {'notBound', 'notAligned'}
        if ~opts.Quiet
            iAlert(parentFigure, opts, plan.message, title, 'info');
        end
        return

    case 'inLibrary'
        % Shown only (section 4 Q7): nothing is recorded.
        if ~opts.Quiet
            iAlert(parentFigure, opts, plan.message, title, 'info');
        end
        return
end

% canAdd
choice = char(opts.Choice);
if opts.Interactive
    question = sprintf(['%s\n\nAdding copies the committed reference image into the ' ...
        'library of subject %s (project %s), so the subject''s other recordings can ' ...
        'be aligned to it. This folder is not changed. "Add and Make Active" also ' ...
        'makes it the subject''s active reference.'], plan.message, plan.subjectID, ...
        plan.projectName);
    answer = uiconfirm(parentFigure, question, 'Add Folder Reference to Library', ...
        'Options', {'Add to Library', 'Add and Make Active', 'Not Now'}, ...
        'DefaultOption', 1, 'CancelOption', 3, 'Icon', 'question');
    choice = iChoiceFromAnswer(answer);
end
if strcmp(choice, 'skip')
    result.message = ['The folder reference was not added to the library. Add it later: ' ...
        'Management > Project > Add Folder Reference to Library.'];
    return
end

setActive = strcmp(choice, 'addActive');
try
    result.resourceUUID = adoptFolderReference(saveFolder, 'SetActive', setActive);
catch ME
    result.severity = 'error';
    result.message = sprintf('Could not add the folder reference to the library: %s', ...
        ME.message);
    [~, lockMessage, isLockError] = offerStaleLockClearing(parentFigure, ME, ...
        'Folder', saveFolder, 'Interactive', opts.Interactive);
    if isLockError
        result.message = lockMessage;
        result.severity = 'warning';
    else
        iAlert(parentFigure, opts, result.message, title, 'error');
    end
    return
end
result.changed = true;
result.severity = 'success';
if setActive
    result.message = sprintf(['Reference "%s" added to the library of subject %s and ' ...
        'made active.'], plan.referenceName, plan.subjectID);
else
    result.message = sprintf('Reference "%s" added to the library of subject %s.', ...
        plan.referenceName, plan.subjectID);
end
end

function choice = iChoiceFromAnswer(answer)
switch answer
    case 'Add to Library'
        choice = 'add';
    case 'Add and Make Active'
        choice = 'addActive';
    otherwise
        choice = 'skip';
end
end

function iAlert(parentFigure, opts, message, title, icon)
if opts.Interactive
    uialert(parentFigure, message, title, 'Icon', icon);
end
end
