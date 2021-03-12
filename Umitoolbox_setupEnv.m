function Umitoolbox_setupEnv
% PROTOTYPE FOR AN EVENTUAL INSTALLATION FUNCTION FOR THE TOOLBOX.
% For now, it creates an environment variable with the folder where the
% toolbox is located. This is used by the GUI scripts to find folders and
% files in the computer.
saveDir = uigetdir(pwd, 'Select Toolbox Folder');
sys = computer;
myenv = getenv('Umitoolbox');
if isfolder(myenv)
    disp('Environment Variable Umitoolbox already exists!')
    return
end
if strcmp(sys, 'PCWIN64')
    system(['SETX Umitoolbox ' saveDir]);
else
    disp('Cant set an environment variable in this machine. This is temporary. Contact Labeo for details');
end
disp('done!')
end