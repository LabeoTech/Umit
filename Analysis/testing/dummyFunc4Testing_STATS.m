function out = dummyFunc4Testing_STATS(object,SaveFolder, varargin)
% FUNCTEMPLATE is a dummy function that serves as a template for creating
% functions compatible with the toolbox.
%
% The arguments have to be stated in the following order:
% 1 - Input : keywords are {File, RawFolder or SaveFolder}
% 2 - SaveIn: save directory (fullpath).
% 3 - Output: Output file name or cell array of file names.
% 4 - opts: structure containing optional parameters for the function.
%
% Defaults:
default_Output = 'dummyStatsFile.mat'; 
%%% Arguments parsing and validation %%%
p = inputParser;
% The input of the function must be a File , RawFolder or SaveFolder
% Imaging Object:
addRequired(p, 'object')
addRequired(p, 'SaveFolder', @isfolder);
addOptional(p, 'Output', default_Output)
% Parse inputs:
parse(p,object, SaveFolder, varargin{:});
% Initialize Variables:
SaveFolder = p.Results.SaveFolder;
outFile = p.Results.Output;
%%%%
% Run your code here:
disp('This is a dummy function for Stats Module testing!');
n = randi(200);
data = randn(5,n);
% Create Meta Data variables:
labels{1} = {'var1', 'var2', 'var3', 'var4', 'var5'};
dim_names = {'O','X'};
% Save to .MAT file:
save2Mat(fullfile(SaveFolder, outFile), single(data), labels, dim_names);
out = outFile;
end
