function plan = classifyAlignmentDatFiles(targetFolder, doTransform)
%CLASSIFYALIGNMENTDATFILES Decide which .dat files of a folder an alignment will transform.
%
%   plan = classifyAlignmentDatFiles(targetFolder)
%   plan = classifyAlignmentDatFiles(targetFolder, doTransform)
%
%   Shared by applyImageAlignmentToFolder (the apply backend) and the Alignment Tool's
%   confirmation dialog, so the file list the user confirms is the file list that is
%   rewritten.
%
%   Rule:
%       - Every .dat is aligned, temporary pipeline outputs (PMTMP_*) included: they can
%         be saved as permanent files later, and a left-behind old-size copy would be a
%         wrong-size file in an aligned folder.
%       - A permanent .dat that cannot be aligned (not single precision, first two axes
%         not Y and X, a frame size different from the rest) is an error, as before.
%       - A temporary PMTMP_* file that cannot be aligned (unreadable header, not single
%         precision, first two axes not Y and X, or a frame size different from the
%         folder's image size, for example the output of a crop) is skipped with a reason
%         instead of blocking the whole folder.
%
%   The folder's image size is the one shared by the permanent .dat files; only when the
%   folder has no readable permanent .dat does it come from the temporary ones.
%
%   Output (struct):
%       files        - dir entries of the .dat files that will be transformed.
%       infos        - loadMetaData of each (cell, aligned with files); empty cells when
%                      doTransform is false.
%       isTemporary  - logical column, aligned with files.
%       skipped      - struct array (name, reason) of temporary files left out.
%       sizeYX       - [Y X] shared by the files to transform; [] when there is none.
%
%   Errors (identifiers as in applyImageAlignmentToFolder):
%       ImageAlignmentTool:MissingDatFile, ImageAlignmentTool:UnsupportedLayout,
%       ImageAlignmentTool:InconsistentDatSizes, and the single-precision error of
%       assertSingleDat, raised for permanent files only.
%
%   See also applyImageAlignmentToFolder.

if nargin < 2 || isempty(doTransform)
    doTransform = true;
end

targetFolder = char(string(targetFolder));
datFiles = dir(fullfile(targetFolder, '*.dat'));
nFiles = numel(datFiles);

isTemporary = false(nFiles, 1);
for iFile = 1:nFiles
    isTemporary(iFile) = startsWith(datFiles(iFile).name, 'PMTMP_', 'IgnoreCase', true);
end

plan = struct( ...
    'files', datFiles, ...
    'infos', {cell(1, nFiles)}, ...
    'isTemporary', isTemporary, ...
    'skipped', struct('name', {}, 'reason', {}), ...
    'sizeYX', []);

if ~doTransform || nFiles == 0
    return
end

% Describe every file. A permanent file that cannot be described is an error; a
% temporary one is skipped.
keep = true(nFiles, 1);
for iFile = 1:nFiles
    datPath = fullfile(datFiles(iFile).folder, datFiles(iFile).name);
    if ~isfile(datPath)
        error('ImageAlignmentTool:MissingDatFile', ...
            'DAT file disappeared during preflight: %s', datPath);
    end

    try
        info = loadMetaData(datPath);
        % umIT supports only single-precision .dat files: refuse before anything is written.
        assertSingleDat(info, 'Context', sprintf('"%s"', datFiles(iFile).name));
        dimNames = cellstr(string(info.dimNames));
        dimNames = dimNames(:).';
        if numel(dimNames) < 2 || ~isequal(dimNames(1:2), {'Y', 'X'})
            error('ImageAlignmentTool:UnsupportedLayout', ...
                ['DAT file "%s" has dimensions {%s}. Alignment needs image data ' ...
                 'whose first two axes are Y and X.'], ...
                datFiles(iFile).name, strjoin(dimNames, ','));
        end
        plan.infos{iFile} = info;
    catch ME
        if ~isTemporary(iFile)
            rethrow(ME);
        end
        plan.skipped(end + 1) = struct('name', datFiles(iFile).name, ...
            'reason', strtrim(strrep(ME.message, newline, ' ')));
        keep(iFile) = false;
    end
end

% The folder's image size comes from the permanent files; the temporary files only when
% there are none.
described = keep & ~cellfun(@isempty, plan.infos(:));
permanent = described & ~isTemporary;
if any(permanent)
    reference = permanent;
elseif any(described)
    reference = described;
else
    reference = false(nFiles, 1);
end

sizeYX = [];
if any(reference)
    sizes = cellfun(@(info) [datAxisSize(info, 'Y'), datAxisSize(info, 'X')], ...
        plan.infos(reference), 'UniformOutput', false);
    sizes = vertcat(sizes{:});
    if size(unique(sizes, 'rows'), 1) > 1
        names = {datFiles(reference).name};
        error('ImageAlignmentTool:InconsistentDatSizes', ...
            ['The .dat files in the folder have different frame sizes (%s). ' ...
             'No file was changed.'], ...
            strjoin(cellfun(@(n, k) sprintf('%s: %d x %d', n, sizes(k, 1), sizes(k, 2)), ...
            names, num2cell(1:numel(names)), 'UniformOutput', false), '; '));
    end
    sizeYX = sizes(1, :);
end

% A temporary file with a different frame size is skipped, not an error.
for iFile = 1:nFiles
    if ~described(iFile) || ~isTemporary(iFile) || isempty(sizeYX)
        continue
    end
    info = plan.infos{iFile};
    fileSizeYX = [datAxisSize(info, 'Y'), datAxisSize(info, 'X')];
    if ~isequal(fileSizeYX, sizeYX)
        plan.skipped(end + 1) = struct('name', datFiles(iFile).name, 'reason', sprintf( ...
            'frames of %d x %d, but the folder''s image size is %d x %d', ...
            fileSizeYX(1), fileSizeYX(2), sizeYX(1), sizeYX(2)));
        keep(iFile) = false;
    end
end

plan.files = datFiles(keep);
plan.infos = plan.infos(keep);
plan.isTemporary = isTemporary(keep);
plan.sizeYX = sizeYX;
end
