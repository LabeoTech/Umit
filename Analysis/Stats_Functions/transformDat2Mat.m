function outDataStat = transformDat2Mat(data, metaData)
% TRANSFORMDAT2MAT creates a .MAT file compatible with the "Analysis" tab of the
% "umIToolbox" app from a .DAT file containing any Imaging data. 
% The data is stored under a "dummy observation".

default_Output = 'datFile_transf.mat'; %#ok This line is here just for Pipeline management.
%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'data',@(x) isnumeric(x)); % Validate if the input is a 3-D numerical matrix:
addRequired(p,'metaData', @(x) isa(x,'matlab.io.MatFile') | isstruct(x)); % MetaData associated to "data".
% Initialize Variables:
parse(p, data, metaData);
data = p.Results.data;
metaData = p.Results.metaData;
%%%%
dim_names = [metaData.dim_names, {'O'}];
outDataStat = save2Mat([], {data} ,{'genericObs'}, dim_names, 'appendMetaData', metaData,'genFile', false);
end