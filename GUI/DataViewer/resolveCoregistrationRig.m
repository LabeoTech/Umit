function target = resolveCoregistrationRig(dataFolder, varargin)
%RESOLVECOREGISTRATIONRIG Library rig that a dual-camera calibration is saved into.
%
%   target = resolveCoregistrationRig(dataFolder)
%   target = resolveCoregistrationRig(dataFolder, 'RigUUID', rigUUID)
%
%   Used by DataViewer_Coreg2Cams. Read-only (never writes a ledger or the
%   rig library). Never falls back to the machine's default rig.
%   Resolution order:
%       1) the rig bound in DATAFOLDER's rig ledger (UMITRigLedger), when
%          the folder is bound;
%       2) otherwise RigUUID, when the caller supplies one;
%       3) otherwise target.needsSelection is true: the caller lets the user
%          pick a library rig (promptSelectLibraryRig, Active rig
%          pre-selected).
%
%   TARGET has fields needsSelection, rigUUID, rigID, displayName, origin
%   ('folderLedger', 'caller' or '') and message (why a selection is
%   needed, for the dialog).
%
%   Errors (Umitoolbox:Coreg2Cams:*):
%       boundRigNotInLibrary - the folder's bound rig is not in this rig
%                              library (import it from its exported
%                              .umitrig file with UMITRigStore.importRig)
%       boundRigArchived     - the folder's bound rig is archived (restore it
%                              in the Rig Manager to save calibrations)
%   An invalid folder ledger throws Umitoolbox:RigLedger:invalidLedger.

p = inputParser;
p.FunctionName = 'resolveCoregistrationRig';
addParameter(p, 'RigUUID', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
parse(p, varargin{:});

target = struct('needsSelection', false, 'rigUUID', '', 'rigID', '', ...
    'displayName', '', 'origin', '', 'message', '');
dataFolder = char(string(dataFolder));

if ~isempty(dataFolder) && isfolder(dataFolder)
    ledger = UMITRigLedger.read(dataFolder);
    if ~isempty(ledger.binding)
        sourceRig = ledger.binding.sourceRig;
        try
            rigStore = UMITRigStore.open(sourceRig.uuid);
        catch
            error('Umitoolbox:Coreg2Cams:boundRigNotInLibrary', ...
                ['The folder "%s" is bound to rig "%s" (%s), which is not in this ' ...
                'rig library. Import the rig from its exported .umitrig file in the ' ...
                'MATLAB Command Window, UMITRigStore.importRig(<file>) (the Rig ' ...
                'Manager has no rig import), to save its calibration here.'], ...
                dataFolder, sourceRig.displayName, sourceRig.rigID);
        end
        RigInfo = rigStore.getRigInfo();
        if strcmpi(RigInfo.status, 'archived')
            error('Umitoolbox:Coreg2Cams:boundRigArchived', ...
                ['The folder "%s" is bound to rig "%s" (%s), which is archived. ' ...
                'Restore it in the Rig Manager to save a new calibration for it.'], ...
                dataFolder, RigInfo.displayName, RigInfo.rigID);
        end
        target = iTargetFromRigInfo(target, RigInfo, 'folderLedger');
        return
    end
end

callerUUID = char(string(p.Results.RigUUID));
if ~isempty(callerUUID)
    target = iTargetFromRigInfo(target, UMITRigStore.open(callerUUID).getRigInfo(), 'caller');
    return
end

target.needsSelection = true;
target.message = ['The current data folder is not bound to a rig. Choose the library ' ...
    'rig whose dual-camera calibration you want to create or edit.'];
end

function target = iTargetFromRigInfo(target, RigInfo, origin)
target.rigUUID = RigInfo.uuid;
target.rigID = RigInfo.rigID;
target.displayName = RigInfo.displayName;
target.origin = origin;
end
