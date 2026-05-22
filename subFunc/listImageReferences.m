function referenceTable = listImageReferences(opts)
%LISTIMAGEREFERENCES List ImageReference files under the UMIT reference tree.
%
%   referenceTable = listImageReferences()
%   referenceTable = listImageReferences('projectName', projectName)
%   referenceTable = listImageReferences(..., 'includeInvalid', true)
%
%   Scans:
%       getUmitFolder('referenceImages', 'create', false)
%
%   recursively for ImageReference_*.mat files and returns a table with
%   metadata suitable for an Image Reference manager GUI.

arguments
    opts.projectName (1,1) string = ""
    opts.groupName (1,1) string = ""
    opts.mouseID (1,1) string = ""
    opts.sessionID (1,1) string = ""
    opts.includeInvalid (1,1) logical = false
end

referenceRoot = getUmitFolder('referenceImages', 'create', false);

if isempty(referenceRoot)
    referenceTable = iEmptyReferenceTable();
    return
end

fileList = dir(fullfile(referenceRoot, '**', 'ImageReference_*.mat'));

if isempty(fileList)
    referenceTable = iEmptyReferenceTable();
    return
end

rows = cell(numel(fileList), 1);
nRows = 0;

umitRoot = getUmitFolder();

for iFile = 1:numel(fileList)

    fileAbs = fullfile(fileList(iFile).folder, fileList(iFile).name);
    fileRel = iMakeRelativeToRoot(fileAbs, umitRoot);

    isValid = true;
    errorMessage = "";

    try
        ImageReference = loadImageReference(fileAbs);

        projectName = string(ImageReference.projectName);
        groupName = string(ImageReference.groupName);
        mouseID = string(ImageReference.mouseID);
        sessionID = string(ImageReference.sessionID);
        referenceName = string(ImageReference.name);
        description = string(ImageReference.description);
        sourceFolder = string(ImageReference.sourceFolder);
        createdOn = ImageReference.createdOn;
        imageSizeYX = string(sprintf('%d x %d', ...
            ImageReference.imageSizeYX(1), ImageReference.imageSizeYX(2)));

    catch ME
        isValid = false;
        errorMessage = string(ME.message);

        if ~opts.includeInvalid
            continue
        end

        projectName = "";
        groupName = "";
        mouseID = "";
        sessionID = "";
        referenceName = "";
        description = "";
        sourceFolder = "";
        createdOn = NaT;
        imageSizeYX = "";
    end

    if ~iMatchesFilter(projectName, opts.projectName)
        continue
    end
    if ~iMatchesFilter(groupName, opts.groupName)
        continue
    end
    if ~iMatchesFilter(mouseID, opts.mouseID)
        continue
    end
    if ~iMatchesFilter(sessionID, opts.sessionID)
        continue
    end

    nRows = nRows + 1;
    rows{nRows} = { ...
        createdOn, ...
        projectName, ...
        groupName, ...
        mouseID, ...
        sessionID, ...
        referenceName, ...
        imageSizeYX, ...
        description, ...
        sourceFolder, ...
        string(fileList(iFile).name), ...
        string(fileRel), ...
        string(fileAbs), ...
        isValid, ...
        errorMessage};
end

rows = rows(1:nRows);

if isempty(rows)
    referenceTable = iEmptyReferenceTable();
    return
end

referenceTable = cell2table(vertcat(rows{:}), ...
    'VariableNames', { ...
        'CreatedOn', ...
        'ProjectName', ...
        'GroupName', ...
        'MouseID', ...
        'SessionID', ...
        'ReferenceName', ...
        'ImageSizeYX', ...
        'Description', ...
        'SourceFolder', ...
        'FileName', ...
        'RelativePath', ...
        'AbsolutePath', ...
        'IsValid', ...
        'ErrorMessage'});

end

function T = iEmptyReferenceTable()
%IEMPTYREFERENCETABLE Return empty table with expected columns.

T = table( ...
    NaT(0,1), ...
    strings(0,1), ...
    strings(0,1), ...
    strings(0,1), ...
    strings(0,1), ...
    strings(0,1), ...
    strings(0,1), ...
    strings(0,1), ...
    strings(0,1), ...
    strings(0,1), ...
    strings(0,1), ...
    strings(0,1), ...
    false(0,1), ...
    strings(0,1), ...
    'VariableNames', { ...
        'CreatedOn', ...
        'ProjectName', ...
        'GroupName', ...
        'MouseID', ...
        'SessionID', ...
        'ReferenceName', ...
        'ImageSizeYX', ...
        'Description', ...
        'SourceFolder', ...
        'FileName', ...
        'RelativePath', ...
        'AbsolutePath', ...
        'IsValid', ...
        'ErrorMessage'});
end

function tf = iMatchesFilter(value, filterValue)
%IMATCHESFILTER Return true if text value matches optional filter.

filterValue = strtrim(string(filterValue));

if filterValue == ""
    tf = true;
    return
end

tf = strcmpi(strtrim(string(value)), filterValue);
end

function relPath = iMakeRelativeToRoot(absPath, rootPath)
%IMAKERELATIVETOROOT Return path relative to rootPath when possible.

absPath = char(string(absPath));
rootPath = char(string(rootPath));

if startsWith(absPath, rootPath)
    relPath = erase(absPath, [rootPath filesep]);
else
    relPath = absPath;
end
end
