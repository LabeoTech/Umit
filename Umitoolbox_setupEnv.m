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
switch sys
    case 'PCWIN64'
        system(['SETX Umitoolbox ' saveDir]);
        disp('done!');
        disp('Restart MATLAB to apply the changes!');
    otherwise
        disp('For this computer, you have to set manually the environment variable "Umitoolbox"');
end

end