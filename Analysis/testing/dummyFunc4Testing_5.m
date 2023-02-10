function outData = dummyFunc4Testing_5(data, varargin)
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
addRequired(p, 'data');
% Parse inputs:
parse(p,data, varargin{:});
%Initialize Variables:
data= p.Results.data;

a =  getEventsFromTTL(data,10); % INDUCE ERROR HERE!!
%%%%
% Run your code here:
disp('This is a dummy function for Pipeline testing!');
a = zeros(3,3, 'single');
save2Dat(fullfile(SaveFolder, Output), a);
out = Output;
end
