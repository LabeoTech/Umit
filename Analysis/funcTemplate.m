function outFile = funcTemplate({File, RawFolder}, SaveFplder, opts, Output)
% FUNCTEMPLATE is a dummy function that serves as a template for creating
% functions compatible with the toolbox.
%
% The arguments have to be stated in the following order:
% 1 - Input : keywords are {File, RawFolder or SaveFolder}
% 2 - SaveIn: save directory (fullpath).
% 3 - Output: Output file name or cell array of file names.
% 4 - opts: structure containing optional parameters for the function.
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
default_opts = struct('Option1', 'value_Option1', 'OptionN', 'value_OptionN');
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
% Output File: 
default_Output = {'OutputFileName'};
addOptional(p, 'Output', default_Output, @(x) ischar(x) || isstring(x) || iscell(x));
% Parse inputs:
parse(p,File/RawFolder, SaveIn, varargin{:});
%Initialize Variables:
File = p.Results.File; 
% OR
RawFolder = p.Results.RawFolder;
SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;
Output = p.Results.Output;
%%%%

% Run your code here:
%

%Save data using save2dat.m function
save2Dat(datFile, data);
% Output file names
outFile = Output;
end



