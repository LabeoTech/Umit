function outData = calculateDF_F0(data, metaData)
% CALCULATEDF_F0 generates DeltaF/F0 movie from FILE.
% Inputs:
%   data: numerical matrix containing imaging data with a "T"ime dimension.
%   metaData: .mat file with meta data associated with "data".
% Output:
%   outData: numerical matrix containing aggregated imaging data.   
%   metaData: .mat file with meta data associated with "outData".
% Defaults:
default_Output = 'deltaF_F0.dat'; %#ok This is here only as a reference for PIPELINEMANAGER.m. 
%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'data',@(x) isnumeric(x)); % Validate if the input is a 3-D numerical matrix:
addRequired(p,'metaData', @(x) isa(x,'matlab.io.MatFile') | isstruct(x)); % MetaData associated to "data".
% Parse inputs:
parse(p,data, metaData);
%Initialize Variables:
data = p.Results.data; 
metaData = p.Results.metaData;
clear p
%%%%

% Find "T"ime dimension:
indxT = find(strcmp('T', metaData.dim_names));
% Calculate baseline over time
bsln = mean(data,indxT);
% Calculate DeltaF/F0:
outData = (data - bsln)./ bsln;
end