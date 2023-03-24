function updateVersionTag(figHandle)
    % This function updates the version tag in the mlapps.    
    path = fileparts(mfilename('fullpath'));
    version = fileread(fullfile(path, 'umit_version.txt'));
    figHandle.Name = [figHandle.Name ' (' strip(version) ')'] ; % append version name    
end