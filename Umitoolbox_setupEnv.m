function Umitoolbox_setupEnv
% PROTOTYPE FOR AN EVENTUAL INSTALLATION FUNCTION FOR THE TOOLBOX.
% For now, it creates an environment variable with the folder where the
% toolbox is located. This is used by the GUI scripts to find folders and
% files in the computer.

% Get UMIT's env variable:
myenv = getenv('Umitoolbox');
if isfolder(myenv)
    answer = questdlg('Environment Variable Umitoolbox already exists. Choose an option:',...
        'Env. variable exists!', 'Finish Setup', 'Redo Setup', 'Cancel', 'Finish Setup');
    switch answer
        case 'Finish Setup'
            b_finish = checkEnvVar;
        case 'Redo Setup'
            setEnvVar;
            return
        otherwise
            disp('Operation cancelled by User')
            return
    end
else
    setEnvVar;
    return
end

if ~b_finish
    error('Failed to complete Umit Setup!')
end

% Edit info.xml file to be able to access the documentation through Matlab "doc":
info_file = fileread(fullfile(myenv,'docs','info.xml'));
path_str = regexp(info_file,'(?<=<help_location>).*(?=</help_location>)', 'match');
% Replace default path string with user's path:
info_file = strrep(info_file,path_str{:}, fullfile(myenv,'docs'));
% Save info_file:
fid = fopen(fullfile(myenv,'docs','info.xml'),'w');
fprintf(fid,'%s',info_file);
fclose(fid);
% Build serch database:
folders = dir(fullfile(myenv, 'docs'));
help_dir = folders([folders.isdir] & startsWith({folders.name}, 'helpsearch-v'));
status = rmdir(fullfile(help_dir.folder, help_dir.name), 's'); %#ok
builddocsearchdb(fullfile(myenv,'docs'))
% Add toolbox to Path for the current Matlab session:
addpath(genpath(getenv('Umitoolbox')));
disp('Everything is set! You can start using the toolbox now!');
end
% Local functions:
function setEnvVar
% This function automatically sets an environment variable named
% "Umitoolbox" for a windows user account.
saveDir = uigetdir(pwd, 'Select Toolbox Folder');
if saveDir == 0
    disp('Operation cancelled by User')
    return
end
switch computer
    case 'PCWIN64'
        [~, cmdOut] = system(['SETX Umitoolbox ' saveDir]);
        if contains(cmdOut, 'success','IgnoreCase',true)
            disp('done!');
            disp('Restart MATLAB to apply the changes and rerun this function!');
        else
            error('Failed to create environment variable! Try to do it manually.')
        end
    otherwise
        disp('For this computer, you have to set manually the environment variable "Umitoolbox"');
end
end

function b_isOk = checkEnvVar
% This function checks if the current Matlab sessions "sees" UMIT's env.
% variable.
b_isOk = false;
if isempty(getenv('Umitoolbox'))
    return
end
[~,cmdOut] = system(['echo %Umitoolbox%']);
if contains(cmdOut, getenv('Umitoolbox'))
    b_isOk = true;
end
end
