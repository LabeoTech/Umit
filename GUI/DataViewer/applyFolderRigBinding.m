function result = applyFolderRigBinding(saveFolder, rigUUID)
%APPLYFOLDERRIGBINDING Bind or rebind a folder to a library rig (DataViewer).
%
%   result = applyFolderRigBinding(saveFolder, rigUUID)
%
%   Re-plans with planFolderRigBinding (so the write never uses a stale
%   preview), then writes UMITRigLedger.bind or rebind with Origin
%   'DataViewer'. Nothing is written when the plan action is 'none'.
%   RESULT has fields action, rigLabel, conflict (the rebind result) and
%   message (a one-line status for DataViewer).

plan = planFolderRigBinding(saveFolder, rigUUID);
result = struct('action', plan.action, 'rigLabel', plan.rigLabel, ...
    'conflict', plan.conflict, 'message', '');
switch plan.action
    case 'bind'
        UMITRigLedger.bind(saveFolder, plan.binding, 'Origin', 'DataViewer');
        result.message = sprintf('Bound rig %s to the folder.', plan.rigLabel);
    case 'rebind'
        result.conflict = UMITRigLedger.rebind(saveFolder, plan.binding, 'Origin', 'DataViewer');
        if result.conflict.hasConflict
            result.message = sprintf(['Rebound the folder to rig %s. The folder is now ' ...
                'CONFLICTED (%s).'], plan.rigLabel, strjoin(result.conflict.parts, ', '));
        else
            result.message = sprintf('Rebound the folder to rig %s.', plan.rigLabel);
        end
    otherwise
        result.message = plan.confirmText;
end
end
