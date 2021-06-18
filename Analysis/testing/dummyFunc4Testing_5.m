function out = dummyFunc4Testing_5(SaveFolder, varargin)
% FUNCTEMPLATE is a dummy function that serves as a template for creating
% functions compatible with the toolbox.
%
% The arguments have to be stated in the following order:
% 1 - Input : keywords are {File, RawFolder or SaveFolder}
% 2 - SaveIn: save directory (fullpath).
% 3 - Output: Output file name or cell array of file names.
% 4 - opts: structure containing optional parameters for the function.
%
default_Output = 'dummyFile5.dat';
%%% Arguments parsing and validation %%%
p = inputParser;
% The input of the function must be a File , RawFolder or SaveFolder
addRequired(p, 'RawFolder', @isfolder)% For Raw Folder as input
% Save folder:
addRequired(p, 'SaveFolder', @isfolder);
% Output:
addOptional(p, 'Output', default_Output);
% Parse inputs:
parse(p,RawFolder, SaveFolder);
%Initialize Variables:
RawFolder = p.Results.RawFolder;
SaveFolder = p.Results.SaveFolder;
Output = p.Results.Output;

%%%%
% Run your code here:
disp('This is a dummy function for Pipeline testing!');
a = zeros(3,3, 'single');
save2Dat(fullfile(SaveFolder, Output), a);
out = Output;
end
