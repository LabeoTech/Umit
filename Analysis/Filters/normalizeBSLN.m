function [outData, metaData] = normalizeBSLN(data, metaData, varargin)
% NORMALIZEBSLN normalizes event-triggered data by a baseline. Here, a
% baseline corresponds to the median of the pixel values before the event
% (pre-event time).
% The data must be a 4-D matrix containing image time series separated by
% trials.
% Inputs:
%   data (4D numerical matrix): Image time series separated by events with dimensions {'E','Y','X','T}.
%   metaData (struct): structure containing smeta data associated with "data".
%   opts (optional) : structure containing extra parameters.
% Output:
%   outData(4D numerical matrix): "data" with values transformed to express
%   DeltaR/R0.
%   metaData: .mat file with meta data associated with "outData".
% Defaults:
default_Output = 'normBSLN.dat'; %#ok This is here only as a reference for PIPELINEMANAGER.m. 
default_opts = struct('baseline_sec','auto', 'b_centerAtOne', false);
opts_values = struct('baseline_sec', {{'auto',Inf}}, 'b_centerAtOne', [true,false]);%#ok. This is here only as a reference for PIPELINEMANAGER.m. 
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
errMsg = 'Input Data must be an Image time series separated by events with dimensions {"E","Y","X",T"}.';
errID = 'umIToolbox:normalizeBSLN:WrongInput';
assert(isequal(metaData.dim_names, {'E','Y','X','T'}), errID, errMsg)
% 
% Update event timings in metaData:
if ~strcmpi(opts.baseline_sec,'auto')   
    % Update event timings in metaData:
    metaData.postEventTime_sec = metaData.postEventTime_sec - opts.baseline_sec;
    metaData.preEventTime_sec = single(opts.baseline_sec);    
end
% Calculate baseline:
bsln = median(outData(:,:,:,1:round(metaData.preEventTime_sec*metaData.Freq)), ...
    4,'omitnan');
% Normalize data from baseline:
disp('Normalizing data by baseline...');
outData = (outData - bsln)./bsln;
if opts.b_centerAtOne
    outData = outData + 1;
end
end