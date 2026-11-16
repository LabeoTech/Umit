function result = promptExportRig(parentFigure, opts)
%PROMPTEXPORTRIG DataViewer Export Rig workflow.
%
%   result = promptExportRig(parentFigure)
%   result = promptExportRig(..., 'Interactive', false, 'RigUUID', uuid, ...
%       'FilePath', filePath)
%
%   Lets the user choose a library rig (promptSelectLibraryRig) and a
%   destination file (uiputfile, default "<rigID>.umitrig"), then writes the
%   rig, its spectra and its coregistration resources to a self-contained
%   .umitrig file with UMITRigStore.exportRig. The rig library is not
%   changed. Install the file on another machine with Management > Rig >
%   Import Rig (promptImportRig).
%
%   Filter sets are not part of a rig. When the library has user-defined
%   filter sets, the user is asked which of them (all pre-selected, see
%   promptSelectFilterSets) to carry in the file together with the filter
%   spectra they use; built-in filter sets exist in every installation and
%   are never offered.
%
%   RESULT has fields changed (always false: nothing in the app changes),
%   message, severity ('info' or 'error'), filePath (the written file,
%   '' when nothing was written) and filterSetIDs (the filter sets written).
%   Errors are shown in a uialert and reported in RESULT, never rethrown.
%
%   Name-Value options:
%       Interactive  - false: no dialogs (tests); RigUUID and FilePath are
%                      then required. Default true.
%       RigUUID      - rig to export. With Interactive true it skips the rig
%                      chooser (the Rig Manager passes its selected rig);
%                      with Interactive false it is required.
%       FilePath     - destination, with Interactive false.
%       FilterSetIDs - with Interactive false, the filter sets to include
%                      (string array or cellstr). Default none.

arguments
    parentFigure
    opts.Interactive (1,1) logical = true
    opts.RigUUID (1,:) char = ''
    opts.FilePath (1,:) char = ''
    opts.FilterSetIDs = {}
end

result = struct('changed', false, 'message', '', 'severity', 'info', ...
    'filePath', '', 'filterSetIDs', {{}});
try
    rigUUID = opts.RigUUID;
    filePath = opts.FilePath;
    filterSetIDs = cellstr(string(opts.FilterSetIDs));

    if opts.Interactive
        if isempty(rigUUID)
            rigUUID = promptSelectLibraryRig(parentFigure, ...
                'Title', 'Export Rig', ...
                'Message', 'Choose the library rig to export to a .umitrig file.');
        end
        if isempty(rigUUID)
            result.message = 'Rig export cancelled.';
            return
        end
        userFilterSets = UMITRigStore.listFilterSets();
        if ~isempty(userFilterSets)
            userFilterSets = userFilterSets(userFilterSets.Origin == "user", :);
        end
        [filterSetIDs, accepted] = promptSelectFilterSets(parentFigure, userFilterSets);
        if ~accepted
            result.message = 'Rig export cancelled.';
            return
        end
        rigInfo = UMITRigStore.open(rigUUID).getRigInfo();
        [fileName, folder] = uiputfile( ...
            {'*.umitrig', 'umIT rig export (*.umitrig)'}, 'Export Rig', ...
            [rigInfo.rigID '.umitrig']);
        if isequal(fileName, 0)
            result.message = 'Rig export cancelled.';
            return
        end
        filePath = fullfile(folder, fileName);
    elseif isempty(rigUUID) || isempty(filePath)
        error('Umitoolbox:promptExportRig:missingInput', ...
            'RigUUID and FilePath are required when Interactive is false.');
    end

    store = UMITRigStore.open(rigUUID);
    rigInfo = store.getRigInfo();
    writtenPath = store.exportRig(filePath, 'FilterSets', filterSetIDs);

    result.filePath = writtenPath;
    result.filterSetIDs = filterSetIDs;
    result.message = sprintf('Exported rig "%s" (%s) to %s.', ...
        rigInfo.displayName, rigInfo.rigID, writtenPath);
    if ~isempty(filterSetIDs)
        result.message = sprintf('%s Included filter set(s): %s.', ...
            result.message, strjoin(filterSetIDs, ', '));
    end
catch ME
    result.message = sprintf('Could not export the rig: %s', ME.message);
    result.severity = 'error';
    if opts.Interactive
        uialert(parentFigure, ME.message, 'Export Rig failed', 'Icon', 'error');
    end
end
end
