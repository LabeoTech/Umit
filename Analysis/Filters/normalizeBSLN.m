function [outData, metaData] = normalizeBSLN(data, metaData, varargin)
% NORMALIZEBSLN normalizes event-triggered and image time series data by a baseline.
% Here, a baseline corresponds to the median of the pixel values before the event
% (pre-event time) or a value in seconds stored in the variable "baseline_sec".
% If the data is an image time series, the baseline value "auto" corresponds 
% to the 20% first frames of the data.
% The data must be a 3-D matrix containg image time series or a 4-D matrix 
% containing image time series separated by trials.
% Inputs:
%   data (3D or 4D numerical matrix): Image time series or image time series separated by events.
%   metaData (struct): structure containing meta data associated with "data".
%   opts (optional) : structure containing extra parameters:
%       baseline_sec("auto" or a number): baseline period(in seconds). If
%           set to "auto" the value stored in the metaData as
%           "preEventTime_sec" will be applied. If the input data is a 3D
%           matrix, the "auto" value corresponds to the 20% first frames of the
%           data.
%       b_centerAtOne (bool): Set to TRUE to center the normalized data at
%           one. Otherwise, the data will be centered at zero.
% Output:
%   outData(3D or 4D numerical matrix): "data" with values transformed to express
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
errMsg = 'Input Data must be an Image time series with dimensions "Y","X",T" or "E","Y","X","T".';
errID = 'umIToolbox:normalizeBSLN:WrongInput';
assert(all(ismember({'Y','X','T'}, metaData.dim_names)), errID, errMsg)
if strcmpi('E',metaData.dim_names)
    % Update event timings in metaData:
    if ~strcmpi(opts.baseline_sec,'auto')
        % Update event timings in metaData:
        metaData.postEventTime_sec = metaData.postEventTime_sec - opts.baseline_sec;
        metaData.preEventTime_sec = single(opts.baseline_sec);
        bsln_sec = metaData.preEventTime_sec;
    else
        bsln_sec = metaData.pre_EventTime_sec;
    end
else
   if strcmpi(opts.baseline_sec, 'auto')
       bsln_sec = .2*size(outData,3)/metaData.Freq; % Set baseline as first 20% of the frames for image time series.
       warning('Setting arbitrary baseline time of %0.2f seconds for image time series input.', bsln_sec);
   else
       bsln_sec = single(opts.baseline_sec);
   end              
end

disp('Normalizing data by baseline...');
% Calculate baseline:
if strcmpi('E',metaData.dim_names)
    % for 4D data:
    bsln = median(outData(:,:,:,1:round(bsln_sec*metaData.Freq)), ...
        4,'omitnan');
else
    % for 3D data:
    bsln = median(outData(:,:,1:round(bsln_sec*metaData.Freq)), ...
        3,'omitnan');
end
% Normalize data from baseline:
outData = (outData - bsln)./bsln;
% Center data at ONE:
if opts.b_centerAtOne
    outData = outData + 1;
end
end