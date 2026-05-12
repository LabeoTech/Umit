function genBackupFolder(SaveFolder, backupChoice)
%GENBACKUPFOLDER Manage backup and deletion of files in a save folder.
%
%   genBackupFolder(SaveFolder)
%   genBackupFolder(SaveFolder, backupChoice)
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
%   Notes:
%       - If backupChoice is omitted or empty, a dialog asks the user what
%         to do.
%       - Managed files are backed up as a .zip archive stored inside
%         SaveFolder.
%       - After backup creation, managed files are deleted from SaveFolder.
%
%   Examples:
%       genBackupFolder('C:\Data', 'ERASE')
%       genBackupFolder('C:\Data', 'GENBACKUP')
%       genBackupFolder('C:\Data', 'MyBackupName')

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

% Check for fixed files that should never be managed here
fixedFiles = [ ...
    dir(fullfile(SaveFolder, '*.bin')); ...          % raw binary files
    dir(fullfile(SaveFolder, '*.tif')); ...          % TIFF files from TIFF import
    dir(fullfile(SaveFolder, '*info.txt')); ...      % info.txt / ExperimentInfo.txt
    dir(fullfile(SaveFolder, 'Comments.txt')); ...   % comments
    dir(fullfile(SaveFolder, 'Snapshot*.png'))];     % snapshots

% Get list of files to manage
movFiles = dir(SaveFolder);
movFiles([movFiles.isdir] == 1) = [];

if ~isempty(fixedFiles)
    movFiles(ismember({movFiles.name}, {fixedFiles.name})) = [];
end

if isempty(movFiles)
    return
end

% Resolve backup choice
if nargin < 2 || isempty(backupChoice)
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

% Normalize GENBACKUP behavior to a timestamped zip base name
if strcmpi(backupChoice, 'GENBACKUP')
    backupChoice = ['bkp_' datestr(now(), 'yyyymmddHHMMSS')];
end

% Create zip backup when requested
if ~strcmpi(backupChoice, 'ERASE')
    zipBaseName = backupChoice;
    if endsWith(zipBaseName, '.zip', 'IgnoreCase', true)
        zipFilePath = fullfile(SaveFolder, zipBaseName);
    else
        zipFilePath = fullfile(SaveFolder, [zipBaseName '.zip']);
    end

    fileList = fullfile({movFiles.folder}, {movFiles.name});

    try
        zip(zipFilePath, fileList);
    catch ME
        error('Umitoolbox:genBackupFolder:ZipFailed', ...
            'Failed to create backup zip "%s".\nOriginal error: %s', ...
            zipFilePath, ME.message);
    end
end

% Delete managed files from the current folder
arrayfun(@(x) delete(fullfile(x.folder, x.name)), movFiles);

end