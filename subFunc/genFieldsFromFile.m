function varNames = genFieldsFromFile(filename, varargin)
% GENFIELDSFROMFILE returns the variable names from the specified .mat file.
% This function is used by some of the analysis function (e.g.
% run_HemoCompute.m) to dynamically update the default options without the
% need to modify the funcion code.
% In order to work, this function should replace the field value of the
% default_opts and opts_values variables inside the analysis function.
% See run_HemoCompute.m for an example.
%
% Inputs:
%   filename (string): The name of the .mat file to read.
%   numVars (scalar, optional): The number of variable names to retrieve from the file.
%
% Outputs:
%   varNames (string or cell array of strings): The variable names from the .mat file.

% Parse input arguments
p = inputParser;
addRequired(p, 'filename', @(x) ischar(x) && endsWith(x, '.mat'));
addOptional(p, 'numVars', [], @(x) isscalar(x) && x > 0);
parse(p, filename, varargin{:});
numVars = p.Results.numVars;
clear p
% Load the .mat file
file = load(filename);%#ok

% Get the variable names
vars = fieldnames(file);
if isempty(numVars)
    numVars = length(vars);
end
% Return the variable names
if numVars == 1
    varNames = vars{1};
else
    varNames = vars(1:numVars);
end

end
