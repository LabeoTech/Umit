function info = loadAcqInfo(folder)
% Loads the "AcqInfo.mat" file in the folder
% Raises an error if the file doesn't exist

try
    info = load(fullfile(folder, 'AcqInfos.mat'));
catch 
    ME = MException('MATLAB:loadError','"AcqInfos.mat" file not found in "%s"',folder);
    throwAsCaller(ME)
end


