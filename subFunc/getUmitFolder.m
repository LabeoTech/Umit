function folderPath = getUmitFolder(subFolder, opts)
%GETUMITFOLDER Return or create a folder under the UMIT config root.
%
%   folderPath = getUmitFolder()
%   folderPath = getUmitFolder(subFolder)
%   folderPath = getUmitFolder(subFolder, 'create', true)
%   folderPath = getUmitFolder(subFolder, 'create', false)
%
%   Returns the UMIT configuration root folder:
%
%       <USER>/Documents/LabeoTech/umIT
%
%   or a requested relative subfolder below it. The root holds the project
%   store ('projects'), the rig library ('rigs') and the preferences.
%
%   Root override: when the environment variable UMIT_ROOT is set to a
%   non-empty absolute folder path, that folder is the root instead. Tests
%   set it through UMITTemporaryRootFixture so that they never read or
%   write the real stores. An empty value means "not set".
%
%   Inputs:
%       subFolder - Optional relative folder path below the UMIT root.
%                   Folder levels can be separated with "/" or "\".
%
%   Name-Value options:
%       create    - If true, missing folders are created.
%                   If false and the folder does not exist, returns ''.
%                   Default: true.
%
%   Examples:
%       rootFolder = getUmitFolder();
%
%       refFolder = getUmitFolder('referenceImages');
%
%       sessionRefFolder = getUmitFolder( ...
%           'referenceImages/ProjectA/Group1/Mouse01/Session03');
%
%       maybeFolder = getUmitFolder('referenceImages/ProjectA', ...
%           'create', false);

arguments
    subFolder (1,1) string = ""
    opts.create (1,1) logical = true
end

folderPath = strtrim(getenv('UMIT_ROOT'));
if isempty(folderPath)
    if ispc
        root = getenv('USERPROFILE');
    else
        root = getenv('HOME');
    end
    folderPath = fullfile(root, 'Documents', 'LabeoTech','umIT');
elseif ~iIsAbsolutePath(folderPath)
    error('Umitoolbox:UmitFolder:invalidRootOverride', ...
        'UMIT_ROOT must be an absolute folder path, got "%s".', folderPath);
end

subFolder = strtrim(subFolder);

if subFolder ~= ""
    % Normalize Windows and Unix-style separators.
    subFolder = replace(subFolder, "\", filesep);
    subFolder = replace(subFolder, "/", filesep);

    if startsWith(subFolder, filesep) || contains(subFolder, ":")
        error('getUmitFolder:InvalidPath', ...
            'subFolder must be a relative path below the UMIT root.');
    end

    pathParts = split(subFolder, filesep);
    pathParts = pathParts(strlength(pathParts) > 0);

    if any(pathParts == "." | pathParts == "..")
        error('getUmitFolder:InvalidPathPart', ...
            'subFolder cannot contain "." or "..".');
    end

    folderPath = fullfile(folderPath, pathParts{:});
end

if ~isfolder(folderPath)
    if opts.create
        mkdir(folderPath);
    else
        folderPath = '';
    end
end

end

function tf = iIsAbsolutePath(folderPath)
%IISABSOLUTEPATH True for a drive, UNC or root-anchored path.
tf = startsWith(folderPath, {'/', '\'}) || ...
    ~isempty(regexp(folderPath, '^[A-Za-z]:[\\/]', 'once'));
end
