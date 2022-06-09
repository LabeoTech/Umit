function dummyFunc4Testing_4(object, SaveFolder, varargin)
% FUNCTEMPLATE is a dummy function that serves as a template for creating
% functions compatible with the toolbox.
%
% The arguments have to be stated in the following order:
% 1 - Input : keywords are {File, RawFolder or SaveFolder}
% 2 - SaveIn: save directory (fullpath).
% 3 - Output: Output file name or cell array of file names.
% 4 - opts: structure containing optional parameters for the function.
%
%Defaults:
default_opts = struct('field1','val1','field2',10);
opts_values = struct('field1',{{'val1', 'val2','val3'}}, 'field2',[1:10]);

%%% Arguments parsing and validation %%%
p = inputParser;
% The input of the function must be a File , RawFolder or SaveFolder
addRequired(p, 'object')
% Save folder:
addRequired(p, 'SaveFolder', @isfolder);
addOptional(p, 'opts',default_opts, @isstruct)
% Parse inputs:
parse(p,object, SaveFolder, varargin{:});
%Initialize Variables:
object = p.Results.object;
SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;
%%%%
% Run your code here:
disp('This is a dummy function for Pipeline testing!');
end


