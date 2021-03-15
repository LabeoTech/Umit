function out = dummyFunc4Testing_3(File, SaveFolder, varargin)
% FUNCTEMPLATE is a dummy function that serves as a template for creating
% functions compatible with the toolbox.
%
% The arguments have to be stated in the following order:
% 1 - Input : keywords are {File, RawFolder or SaveFolder}
% 2 - SaveIn: save directory (fullpath).
% 3 - Output: Output file name or cell array of file names.
% 4 - opts: structure containing optional parameters for the function.
%

default_Output = 'dummyFile3.dat';
%%% Arguments parsing and validation %%%
p = inputParser;
% The input of the function must be a File , RawFolder or SaveFolder
addRequired(p, 'File',@isfile)% For a File as input
% Save folder:
addRequired(p, 'SaveFolder', @isfolder);
% Output:
addOptional(p, 'Output', default_Output);
% Parse inputs:
parse(p,File, SaveFolder);
%Initialize Variables:
File = p.Results.File;
SaveFolder = p.Results.SaveFolder;
Output = p.Results.Output;

%%%%
% Run your code here:
disp('This is a dummy function for Pipeline testing!');
a = zeros(3,3, 'single');
save2Dat(fullfile(SaveFolder, Output), a);
out = Output;
end



