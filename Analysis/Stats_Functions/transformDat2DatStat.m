function outData = transformDat2DatStat(data)
% TRANSFORMDAT2DATSTAT packages the image time series data in a structure 
% to be compatibles with with "Analysis" tab of the
% "umIToolbox" app and saved as a .DATSTAT file.

default_Output = 'datFile_transf.mat'; %#ok This line is here just for Pipeline management.
%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'data',@(x) isnumeric(x) && ndims(x) == 3); % Validate if the input is a 3-D numerical matrix:
% Initialize Variables:
parse(p, data);
%
outData = genDataStructure(data);
end