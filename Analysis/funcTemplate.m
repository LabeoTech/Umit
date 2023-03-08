%% Templates for creating functions in umIT

% Below are N dummy functions that can be used as templates to build your 
% own function in umIT:
% 

function outFile = funcTemplate(File, SaveFolder, opts)
% FUNCTEMPLATE is a dummy function that serves as a template for creating
% functions compatible with the toolbox.
%
% The arguments have to be stated in the following order:
% 1 - Input : keywords are {File, RawFolder or SaveFolder}
% 2 - SaveIn: save directory (fullpath).
% 3 - opts: structure containing optional parameters for the function.
%
%%% Arguments parsing and validation %%%
p = inputParser;
% The input of the function must be a File , RawFolder or SaveFolder
addRequired(p,'File',@isfile)% For a file as input.
addRequired(p, 'RawFolder', @isfolder)% For Raw Folder as input
% Save folder:
addRequired(p, 'SaveFolder', @isfolder);
% Optional Parameters:
% opts structure
default_opts = struct('Output', default_Output, 'Option1', 'value_Option1', 'OptionN', 'value_OptionN');
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
% Parse inputs:
parse(p,File/RawFolder, SaveIn, varargin{:});
%Initialize Variables:
File = p.Results.File; 
% OR
RawFolder = p.Results.RawFolder;
SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;
Output = opts.Output;
%%%%

% Run your code here:
%

%Save data using save2dat.m function
save2Dat(datFile, data);
% Output file names
outFile = Output;
end



