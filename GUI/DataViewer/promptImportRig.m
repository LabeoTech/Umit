function result = promptImportRig(parentFigure, opts)
%PROMPTIMPORTRIG DataViewer Import Rig workflow.
%
%   result = promptImportRig(parentFigure)
%   result = promptImportRig(..., 'Interactive', false, 'FilePath', filePath, ...
%       'IfExists', 'cancel' | 'copy')
%
%   Lets the user choose a .umitrig file (uigetfile) written by Export Rig
%   (promptExportRig) and installs it in the rig library with
%   UMITRigStore.importRig. The rig keeps its identity, its spectra and its
%   coregistration resources. It becomes Active only when the library has no
%   Active rig, otherwise Available; no folder is bound to it.
%
%   When the rig is already in the library (same UUID), or another rig uses
%   its rig ID, the user is asked whether to import it anyway: a separate
%   copy with a new identity in the first case, or the same rig under a new
%   rig ID in the second. The new ID is "<rigID>_copy", "<rigID>_copy2", ...
%
%   Filter sets carried by the file are added to the user library when
%   missing (and named in the message); one already present with the same
%   content is left alone, one with different content under the same ID
%   makes the import fail and nothing is installed.
%
%   RESULT has fields changed (true when a rig was installed), message,
%   severity ('info' or 'error'), rigUUID and rigID of the installed rig ('' when
%   nothing was installed), and filterSetsAdded (IDs written to the library). Errors (an unreadable file, a spectrum ID that
%   exists with different content, ...) are shown in a uialert and reported
%   in RESULT, never rethrown; nothing is installed on failure.
%
%   Name-Value options:
%       Interactive - false: no dialogs (tests); FilePath is then required.
%                     Default true.
%       FilePath    - the .umitrig file, with Interactive false.
%       IfExists    - with Interactive false, what to do when the rig or its
%                     ID already exists: 'cancel' (default) or 'copy'.

arguments
    parentFigure
    opts.Interactive (1,1) logical = true
    opts.FilePath (1,:) char = ''
    opts.IfExists (1,1) string {mustBeMember(opts.IfExists, ["cancel", "copy"])} = "cancel"
end

result = struct('changed', false, 'message', '', 'severity', 'info', ...
    'rigUUID', '', 'rigID', '', 'filterSetsAdded', {{}});
try
    filePath = opts.FilePath;
    if opts.Interactive
        [fileName, folder] = uigetfile( ...
            {'*.umitrig', 'umIT rig export (*.umitrig)'}, 'Import Rig');
        if isequal(fileName, 0)
            result.message = 'Rig import cancelled.';
            return
        end
        filePath = fullfile(folder, fileName);
    elseif isempty(filePath)
        error('Umitoolbox:promptImportRig:missingInput', ...
            'FilePath is required when Interactive is false.');
    end

    try
        [store, summary] = UMITRigStore.importRig(filePath);
    catch ME
        isSameRig = strcmp(ME.identifier, 'Umitoolbox:UMITRigStore:importRigExists');
        isIDTaken = strcmp(ME.identifier, 'Umitoolbox:UMITRigStore:importRigIDExists');
        if ~isSameRig && ~isIDTaken
            rethrow(ME)
        end
        if ~iAcceptCopy(parentFigure, opts, isSameRig, ME.message)
            result.message = 'Rig import cancelled; the rig library was left unchanged.';
            return
        end
        [store, summary] = iImportUnderFreeID(filePath, isSameRig);
    end

    rigInfo = store.getRigInfo();
    result.changed = true;
    result.rigUUID = rigInfo.uuid;
    result.rigID = rigInfo.rigID;
    result.filterSetsAdded = summary.filterSetsAdded;
    result.message = sprintf(['Imported rig "%s" (%s); it is %s in the rig library. ' ...
        'Bind it to a folder with Management > Rig > Folder Rig.'], ...
        rigInfo.displayName, rigInfo.rigID, rigInfo.status);
    if ~isempty(summary.filterSetsAdded)
        result.message = sprintf('%s Added filter set(s): %s.', ...
            result.message, strjoin(summary.filterSetsAdded, ', '));
    end
catch ME
    result.message = sprintf('Could not import the rig: %s', ME.message);
    result.severity = 'error';
    if opts.Interactive
        uialert(parentFigure, ME.message, 'Import Rig failed', 'Icon', 'error');
    end
end
end

function accepted = iAcceptCopy(parentFigure, opts, isSameRig, storeMessage)
%IACCEPTCOPY Ask (or take the test option) whether to import under a new ID.
if ~opts.Interactive
    accepted = opts.IfExists == "copy";
    return
end
if isSameRig
    question = ['This rig is already in the rig library. Import a separate ' ...
        'copy with a new identity (rig ID "<rigID>_copy")?'];
    button = 'Import as copy';
else
    question = sprintf(['%s\n\nImport it under a new rig ID ' ...
        '("<rigID>_copy")?'], storeMessage);
    button = 'Import under new ID';
end
choice = uiconfirm(parentFigure, question, 'Import Rig', ...
    'Options', {button, 'Cancel'}, 'DefaultOption', 2, 'CancelOption', 2, ...
    'Icon', 'question');
accepted = strcmp(choice, button);
end

function [store, summary] = iImportUnderFreeID(filePath, newIdentity)
%IIMPORTUNDERFREEID Import under the first free "<rigID>_copy[N]" rig ID.
baseID = iPeekRigID(filePath);
maxLength = getUMITRigSchema().namingRules.maxIDLength;
for suffixIndex = 1:99
    if suffixIndex == 1
        suffix = '_copy';
    else
        suffix = sprintf('_copy%d', suffixIndex);
    end
    candidate = [baseID(1:min(numel(baseID), maxLength - numel(suffix))) suffix];
    try
        [store, summary] = UMITRigStore.importRig(filePath, 'RigID', candidate, ...
            'NewIdentity', newIdentity);
        return
    catch ME
        if ~strcmp(ME.identifier, 'Umitoolbox:UMITRigStore:importRigIDExists')
            rethrow(ME)
        end
    end
end
error('Umitoolbox:promptImportRig:noFreeRigID', ...
    'Could not find a free rig ID for "%s" after 99 attempts.', baseID);
end

function rigID = iPeekRigID(filePath)
%IPEEKRIGID Rig ID recorded in a .umitrig file ('Rig' when it cannot be read).
rigID = 'Rig';
try
    loaded = load(filePath, 'UMITRigExport', '-mat');
    candidate = char(string(loaded.UMITRigExport.rig.rigID));
    if ~isempty(candidate)
        rigID = candidate;
    end
catch
    % importRig reports an unreadable file itself.
end
end
