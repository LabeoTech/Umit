function genBackupFolder(SaveFolder, backupChoice)
% genBackupFolder - Manages backup and deletion of files in a specified folder.
% This function ignores the raw files created by LabeoTech's imaging
% systems. Therefore, only the files created after one of the data import
% functions will be managed by this function.
%
% Inputs:
%    SaveFolder - The folder where files are saved and managed.
%    backupChoice - Optional parameter to specify the backup action:
%       'ERASE' - Erases all files in the SaveFolder.
%       'GENBACKUP' - Automatically generates a backup folder with a timestamp.
%       Custom folder name - Moves files to a user-specified backup folder.
% Important notes:
%   - If backupChoice is not provided, a prompt will appear to allow the user to set it interactively.
%   - All backup folder options are created as a subfolder inside SaveFolder.
%
% Examples:
%    genBackupFolder('C:\Data', 'ERASE')
%    genBackupFolder('C:\Data', 'GENBACKUP')
%    genBackupFolder('C:\Data', 'MyBackupFolder')

% Check for existing .dat files in the SaveFolder
fixedFiles = [dir(fullfile(SaveFolder, '*.bin'))... % binary files
    ; dir(fullfile(SaveFolder,'*.tif'))... % TIF files (in case of data import from TIF)
    ; dir(fullfile(SaveFolder, '*info.txt'))... % info.txt and/or ExperimentInfo.txt
    ; dir(fullfile(SaveFolder, 'Comments.txt'))... % Comments.txt
    ; dir(fullfile(SaveFolder, 'Snapshot*.png'))]; % Snapshots

% Get list of files to move
movFiles = dir(SaveFolder);
movFiles([movFiles.isdir] == 1) = [];
if ~isempty(fixedFiles)
    movFiles(ismember({movFiles.name}, {fixedFiles.name})) = [];
end
if isempty(movFiles);return;end

% Create prompt if the backupChoice parameter was not set
if nargin < 2 || isempty(backupChoice)
    choice = questdlg('The save folder already contains files. Please choose an option:', ...
        'Folder Contains Files', 'Erase all', 'Create backup', 'Cancel', 'Create backup');
    
    % Process the user's choice
    switch choice
        case 'Erase all'
            backupChoice = 'ERASE';
        case 'Create backup'
            answer = inputdlg('Type backup folder name:', 'backupChoice', ...
                [1 60], {['bkp_' datestr(now(), 'yyyymmddHHMMSS')]});
            if isempty(answer)
                disp('Operation cancelled by User')
                return
            elseif ~isempty(answer{:})
                backupChoice = answer{:};
            else
                backupChoice = ['bkp_' datestr(now(), 'yyyymmddHHMMSS')];
            end
        otherwise
            disp('Operation cancelled by User')
            return
    end
elseif strcmpi(backupChoice, 'GENBACKUP')
    backupChoice = ['bkp_' datestr(now(), 'yyyymmddHHMMSS')];
end

if ~strcmpi(backupChoice, 'ERASE')
    % Create subfolder
    if ~isfolder(fullfile(SaveFolder, backupChoice))
        mkdir(fullfile(SaveFolder, backupChoice));
    end
    % Copy data to subfolder
    arrayfun(@(x) copyfile(fullfile(x.folder, x.name), fullfile(x.folder, backupChoice, x.name)), movFiles)
end
% Delete all files from the current folder (except .bin and .txt)
arrayfun(@(x) delete(fullfile(x.folder, x.name)), movFiles);

end
