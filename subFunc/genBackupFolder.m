function backupPath = genBackupFolder(SaveFolder, backupChoice, varargin)
%GENBACKUPFOLDER Manage backup and optional deletion of files in a save folder.
%
%   backupPath = genBackupFolder(SaveFolder)
%   backupPath = genBackupFolder(SaveFolder, backupChoice)
%   backupPath = genBackupFolder(..., 'eraseFolder', true)
%
%   This function ignores the raw files created by LabeoTech imaging
%   systems. Therefore, only the files created after one of the data-import
%   functions are managed here.
%
%   Inputs:
%       SaveFolder    - Folder where files are saved and managed.
%       backupChoice  - Optional backup action:
%           'ERASE'      : Delete managed files from SaveFolder.
%           'GENBACKUP'  : Create a timestamped .zip backup in SaveFolder.
%           <char>       : Create a .zip backup using the provided base name.
%
%   Name-Value options:
%       eraseFolder   - Logical scalar. If true, delete managed files after
%                       the selected action is executed.
%                       Default: false.
%
%   Output:
%       backupPath    - Path to the created backup .zip file. Empty if no
%                       backup was created.
%
%   Notes:
%       - If backupChoice is omitted or empty, a dialog asks the user what
%         to do.
%       - Managed files are backed up as a .zip archive stored inside
%         SaveFolder.
%       - Managed files are deleted only when eraseFolder is true or when
%         backupChoice is explicitly set to 'ERASE'.
%
%   Examples:
%       genBackupFolder('C:\Data', 'ERASE')
%       genBackupFolder('C:\Data', 'GENBACKUP')
%       genBackupFolder('C:\Data', 'GENBACKUP', 'eraseFolder', true)
%       backupPath = genBackupFolder('C:\Data', 'MyBackupName')

backupPath = '';

% Validate input
if ~(ischar(SaveFolder) || (isstring(SaveFolder) && isscalar(SaveFolder)))
    error('Umitoolbox:genBackupFolder:InvalidSaveFolder', ...
        'SaveFolder must be a character vector or string scalar.');
end
SaveFolder = char(string(SaveFolder));

if ~isfolder(SaveFolder)
    error('Umitoolbox:genBackupFolder:MissingFolder', ...
        'SaveFolder does not exist: "%s".', SaveFolder);
end

if nargin < 2
    backupChoice = [];
end

p = inputParser;
p.FunctionName = 'genBackupFolder';
addParameter(p, 'eraseFolder', false, @(x) islogical(x) && isscalar(x));
parse(p, varargin{:});

eraseFolder = p.Results.eraseFolder;

% Check for fixed files that should never be managed here
fixedFiles = [ ...
    dir(fullfile(SaveFolder, '*.bin')); ...          % raw binary files
    dir(fullfile(SaveFolder, '*.tif')); ...          % TIFF files from TIFF import
    dir(fullfile(SaveFolder, '*.tiff')); ...         % TIFF files from TIFF import
    dir(fullfile(SaveFolder, '*info.txt')); ...      % info.txt / ExperimentInfo.txt
    dir(fullfile(SaveFolder, 'Comments.txt')); ...   % comments
    dir(fullfile(SaveFolder, 'Snapshot*.png'))];     % snapshots

% Get list of files to manage
movFiles = [...
    dir(fullfile(SaveFolder, '*.dat')); ...          % .dat files
    dir(fullfile(SaveFolder, '*.mat')); ...          % .mat files
    dir(fullfile(SaveFolder, '*.roi')); ...          % .roi files
    dir(fullfile(SaveFolder, '*.umt')); ...          % .umt files    
    ];

if ~isempty(fixedFiles)
    movFiles(ismember({movFiles.name}, {fixedFiles.name})) = [];
end

if isempty(movFiles)
    return
end

% Resolve backup choice
if isempty(backupChoice)
    choice = questdlg( ...
        'The save folder already contains files. Please choose an option:', ...
        'Folder Contains Files', ...
        'Erase all', 'Create backup', 'Cancel', ...
        'Create backup');

    switch choice
        case 'Erase all'
            backupChoice = 'ERASE';

        case 'Create backup'
            answer = inputdlg( ...
                'Type backup zip base name:', ...
                'backupChoice', ...
                [1 60], ...
                {['bkp_' datestr(now(), 'yyyymmddHHMMSS')]});

            if isempty(answer)
                disp('Operation cancelled by User')
                return
            elseif ~isempty(answer{1})
                backupChoice = answer{1};
            else
                backupChoice = ['bkp_' datestr(now(), 'yyyymmddHHMMSS')];
            end

        otherwise
            disp('Operation cancelled by User')
            return
    end
else
    backupChoice = char(string(backupChoice));
end

% Explicit ERASE always deletes managed files and does not create a backup.
if strcmpi(backupChoice, 'ERASE')
    arrayfun(@(x) delete(fullfile(x.folder, x.name)), movFiles);
    return
end

% Normalize GENBACKUP behavior to a timestamped zip base name
if strcmpi(backupChoice, 'GENBACKUP')
    backupChoice = ['bkp_' datestr(now(), 'yyyymmddHHMMSS')];
end

% Create zip backup
zipBaseName = backupChoice;
if endsWith(zipBaseName, '.zip', 'IgnoreCase', true)
    zipFilePath = fullfile(SaveFolder, zipBaseName);
else
    zipFilePath = fullfile(SaveFolder, [zipBaseName '.zip']);
end

zipFilePath = iMakeUniqueZipPath(zipFilePath);
fileList = fullfile({movFiles.folder}, {movFiles.name});

try
    zip(zipFilePath, fileList);
    backupPath = zipFilePath;
catch ME
    error('Umitoolbox:genBackupFolder:ZipFailed', ...
        'Failed to create backup zip "%s".\nOriginal error: %s', ...
        zipFilePath, ME.message);
end

% Delete managed files only when explicitly requested.
if eraseFolder
    arrayfun(@(x) delete(fullfile(x.folder, x.name)), movFiles);
end

end

function uniquePath = iMakeUniqueZipPath(zipFilePath)
%IMAKEUNIQUEZIPPATH Return a non-existing zip path.

uniquePath = zipFilePath;

if ~isfile(uniquePath)
    return
end

[folderName, baseName, ext] = fileparts(zipFilePath);

for iFile = 1:1000
    candidatePath = fullfile(folderName, sprintf('%s_%03d%s', baseName, iFile, ext));

    if ~isfile(candidatePath)
        uniquePath = candidatePath;
        return
    end
end

error('Umitoolbox:genBackupFolder:CouldNotCreateUniqueZip', ...
    'Could not create a unique backup zip file name for "%s".', zipFilePath);
end