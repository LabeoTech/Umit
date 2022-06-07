function Umitoolbox_setupEnv
% PROTOTYPE FOR AN EVENTUAL INSTALLATION FUNCTION FOR THE TOOLBOX.
% For now, it creates an environment variable with the folder where the
% toolbox is located. This is used by the GUI scripts to find folders and
% files in the computer.

myenv = getenv('Umitoolbox');
if isfolder(myenv)
    disp('Environment Variable Umitoolbox already exists!')
else
    saveDir = uigetdir(pwd, 'Select Toolbox Folder');
    if saveDir == 0
        disp('Operation cancelled by User')
        return
    end
    sys = computer;
    switch sys
        case 'PCWIN64'
            system(['SETX Umitoolbox ' saveDir]);
            disp('done!');
            disp('Restart MATLAB to apply the changes and rerun this function!');
        otherwise
            disp('For this computer, you have to set manually the environment variable "Umitoolbox"');
    end
end
% Edit info.xml file to be able to access the documentation through Matlab "doc":
info_file = fileread(fullfile(myenv,'html','info.xml'));
path_str = regexp(info_file,'(?<=<help_location>)(\S+)(?=</help_location>)', 'match');
% Replace default path string with user's path:
info_file = strrep(info_file,path_str{:}, fullfile(myenv,'html'));
% Save info_file:
fid = fopen(fullfile(myenv,'html','info.xml'),'w');
fprintf(fid,'%s',info_file);
fclose(fid);
% Add toolbox to Path:
addpath(genpath(saveDir)); 
end