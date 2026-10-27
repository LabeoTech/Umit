function updateVersionTag(figHandle)
% UPDATEVERSIONTAG Append the toolbox version to a figure name.
%
%   updateVersionTag(figHandle) reads the version string stored in
%   subFunc/umit_version.txt and appends it, in parentheses, to the Name
%   property of figHandle (e.g. 'DataViewer (Bipolar-v2.0.0)').

path = fileparts(mfilename('fullpath'));
version = fileread(fullfile(path, 'umit_version.txt'));
figHandle.Name = [figHandle.Name ' (' strip(version) ')'];
end
