function Umitoolbox_setupEnv
% PROTOTYPE FOR AN EVENTUAL INSTALLATION FUNCTION FOR THE TOOLBOX.
% For now, it creates an environment variable with the folder where the
% toolbox is located. This is used by the GUI scripts to find folders and
% files in the computer.
saveDir = uigetdir(pwd, 'Select Toolbox Folder');
setenv('Umitoolbox', saveDir);
disp('done!')
end