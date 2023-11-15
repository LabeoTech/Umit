function [outData, metaData] = normalizeCOND(data, metaData, varargin)
% NORMALIZECOND normalizes event-triggered and image time series data by the
% average of a baseline condition.
% The data must be a 4-D matrix containing image time series separated by trials.

% Inputs:
%   data (4D numerical matrix): Image time series or image time series separated by events.
%   metaData (struct): structure containing meta data associated with "data".
%   opts (optional) : structure containing extra parameters:
%       baselineConditionIndex (positive number): baseline condition index
%       stored in the "events.mat"file.

%       b_centerAtOne (bool): Set to TRUE to center the normalized data at
%           one. Otherwise, the data will be centered at zero.
% Output:
%   outData(4D numerical matrix): "data" with values transformed to express
%   DeltaR/R0.
%   metaData: .mat file with meta data associated with "outData".
% Defaults:
default_Output = 'normByCond.dat'; %#ok This is here only as a reference for PIPELINEMANAGER.m.
default_opts = struct('baselineConditionIndex',1, 'b_centerAtOne', false);
opts_values = struct('baselineConditionIndex', {[1:100]}, 'b_centerAtOne', [true,false]);%#ok. This is here only as a reference for PIPELINEMANAGER.m.
%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'data',@(x) isnumeric(x)); % Validate if the input is a numerical matrix:
addRequired(p,'metaData', @(x) isa(x,'matlab.io.MatFile') | isstruct(x)); % MetaData associated to "data".
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
% Parse inputs:
parse(p,data, metaData, varargin{:});
% Initialize Variables:
outData = p.Results.data;
metaData = p.Results.metaData;
opts = p.Results.opts;
clear p
%%%%
% Validate if the data is an Image time series by Events:
errMsg = 'Input Data must be an Image time series with dimensions "E","Y","X","T".';
errID = 'umIToolbox:normalizeBSLN:WrongInput';
assert(all(ismember({'E','Y','X','T'}, metaData.dim_names)), errID, errMsg)

disp('Normalizing data by condition...');
% Calculate baseline:
idx = metaData.eventID == opts.baselineConditionIndex;
bsln = mean(outData(idx,:,:,:), 1,'omitnan');
% Normalize data from baseline:
outData = (outData - bsln)./bsln;
% Center data at ONE:
if opts.b_centerAtOne
    outData = outData + 1;
end
end