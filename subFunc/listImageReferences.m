function referenceTable = listImageReferences(opts)
%LISTIMAGEREFERENCES List active ImageReference files.
%
%   referenceTable = listImageReferences()
%   referenceTable = listImageReferences('projectName', projectName)
%   referenceTable = listImageReferences(..., 'includeInvalid', true)
%   referenceTable = listImageReferences(..., 'includeArchived', true)
%
%   Scans the UMIT ImageReference tree:
%
%       getUmitFolder('referenceImages', 'create', false)
%
%   and returns metadata from ImageReference_*.mat files.
%
%   By default, files inside:
%
%       referenceImages/_archive/
%
%   are excluded. This means that after the ImageReferenceManager archives a
%   file and refreshes the table, the archived row disappears from the active
%   table.
%
%   Name-Value options:
%       projectName      - Optional project-name filter.
%       groupName        - Optional group-name filter.
%       mouseID          - Optional mouse-ID filter.
%       sessionID        - Optional session-ID filter.
%       includeInvalid   - If true, include invalid files with error text.
%                          Default: false.
%       includeArchived  - If true, include archived ImageReference files.
%                          Default: false.
%
%   Output:
%       referenceTable - Table with ImageReference metadata.

arguments
    opts.projectName (1,1) string = ""
    opts.groupName (1,1) string = ""
    opts.mouseID (1,1) string = ""
    opts.sessionID (1,1) string = ""
    opts.includeInvalid (1,1) logical = false
    opts.includeArchived (1,1) logical = false
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

% Archived files live under referenceImages/_archive. They should not appear
% in the active manager table unless explicitly requested.
if ~opts.includeArchived
    keepFile = true(numel(fileList), 1);

    for iFile = 1:numel(fileList)
        keepFile(iFile) = ~iIsArchivedReferenceFolder( ...
            fileList(iFile).folder, referenceRoot);
    end

    fileList = fileList(keepFile);
end

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
    isArchived = iIsArchivedReferenceFolder(fileList(iFile).folder, referenceRoot);

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
        isArchived, ...
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
        'IsArchived', ...
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
        'IsArchived', ...
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

absPathNorm = strrep(absPath, '/', filesep);
absPathNorm = strrep(absPathNorm, '\', filesep);

rootPathNorm = strrep(rootPath, '/', filesep);
rootPathNorm = strrep(rootPathNorm, '\', filesep);

if endsWith(rootPathNorm, filesep)
    rootPrefix = rootPathNorm;
else
    rootPrefix = [rootPathNorm filesep];
end

if startsWith(absPathNorm, rootPrefix)
    relPath = extractAfter(string(absPathNorm), strlength(rootPrefix));
    relPath = char(relPath);
else
    relPath = absPath;
end
end

function tf = iIsArchivedReferenceFolder(folderPath, referenceRoot)
%IISARCHIVEDREFERENCEFOLDER True if folder is inside referenceImages/_archive.

folderPath = char(string(folderPath));
referenceRoot = char(string(referenceRoot));

folderPath = strrep(folderPath, '/', filesep);
folderPath = strrep(folderPath, '\', filesep);

referenceRoot = strrep(referenceRoot, '/', filesep);
referenceRoot = strrep(referenceRoot, '\', filesep);

if startsWith(folderPath, referenceRoot)
    relFolder = erase(folderPath, referenceRoot);
else
    relFolder = folderPath;
end

relFolder = strip(string(relFolder), filesep);
pathParts = split(relFolder, filesep);
pathParts = pathParts(strlength(pathParts) > 0);

tf = any(strcmpi(pathParts, "_archive"));
end