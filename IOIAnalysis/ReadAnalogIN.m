function AnalogIN = ReadAnalogIN(FolderPath,Infos)
% READANALOGIN reads the Analog inputs from the ai_xxxx.bin files created by the
% the OiS200 systems.
%
% Inputs:
%   FolderPath (char): Folder containing the binary files.
%   Infos (struct): Structure containing the acquisition information
%   (created by the "ReadInfoFile.m" function.
% Output:
%   AnalogIN (matrix) : Time x Channel array containing the analog inputs
%   data (in Volts).

aiFilesList = dir(fullfile(FolderPath,'ai_*.bin'));
if isempty(aiFilesList)
    warning(['Analog Input files (ai_xxxx.bin) not found in "' FolderPath '"'])
    return
end
disp('Reading analog inputs...')
% Opening of the files:
AnalogIN = [];
for ind = 1:size(aiFilesList,1)
    data = memmapfile(fullfile(FolderPath,aiFilesList(ind).name),...
        'Offset', 5*4, 'Format', 'double', 'repeat', inf);
    tmp = data.Data;
    tmp = reshape(tmp, 1e4, Infos.AINChannels, []);
    tmp = permute(tmp,[1 3 2]);
    tmp = reshape(tmp,[],Infos.AINChannels);
    AnalogIN = [AnalogIN; tmp];
end
% Crop to first and last camera triggers:
camT = diff(AnalogIN(:,1) > 2.5); camT = [camT;NaN];
camTOn = find(camT == 1,1,'first');
camTOff = find(camT == -1,1,'last');
AnalogIN = AnalogIN(camTOn:camTOff,:);
disp('Analog Inputs read!')
end