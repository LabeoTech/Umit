function outData = calculateDF_F_byEvent(data, metaData, varargin)
% CALCULATEDF_F_BYEVENT generates DeltaF/F0 movie from DATA.
% The data must be a 4-D matrix containing image time series separated by
% trials.
% Inputs:
%   data (4D numerical matrix): Image time series separated by events with dimensions {'E','Y','X','T}.
%   metaData: .mat file with meta data associated with "data".
%   opts (struct): opts structure containing :
%       b_applyDetrend (bool): If true, applies "detrend" function over
%       time.
% Output:
%   outData(4D numerical matrix): "data" with values transformed to express
%   DeltaF/F0.
% Defaults:
default_Output = 'deltaF_F_ByEvent.dat'; %#ok This is here only as a reference for PIPELINEMANAGER.m. 
default_opts = struct('b_applyDetrend',false);
%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'data',@(x) isnumeric(x)); % Validate if the input is a 3-D numerical matrix:
addRequired(p,'metaData', @(x) isa(x,'matlab.io.MatFile') | isstruct(x)); % MetaData associated to "data".
addOptional(p, 'opts', default_opts, @isstruct);
% Parse inputs:
parse(p,data, metaData, varargin{:});
%Initialize Variables:
outData = p.Results.data; 
metaData = p.Results.metaData;
opts = p.Results.opts;
clear p
%%%%
% Validate if the data is an Image time series by Events:
errMsg = 'Input Data must be an Image time series separated by events with dimensions {"E","Y","X",T"}.';
errID = 'umIToolbox:calculateDF_F0_byEvent:WrongInput';
assert(isequal(metaData.dim_names, {'E','Y','X','T'}), errID, errMsg)

% Perform linear detrending on data:
if opts.b_applyDetrend    
    medFr = 7;
    disp('Detrending...');    
    delta_y = mean(outData(:,:,:,end-medFr:end),4, 'omitnan') - mean(outData(:,:,:,1:medFr),4,'omitnan');
    delta_x =(size(outData,4)- medFr);
    M = delta_y./delta_x; clear delta_*    
    trend = bsxfun(@times,M,permute(linspace(-2,size(outData,4)-medFr+4,...
        size(outData,4)),[4 3 1 2])) + outData(:,:,:,3);       
    disp('Calculating DeltaF/F ...');
    outData = bsxfun(@rdivide,bsxfun(@minus,outData,trend),trend);
else
    % Calculate baseline:
    bsln = median(outData(:,:,:,1:round(metaData.preEventTime_sec*metaData.Freq)), ...
        4,'omitnan');
    % Normalize data to get DeltaF/F values
    disp('Calculating DeltaF/F ...');
    outData = bsxfun(@rdivide,bsxfun(@minus,outData,bsln),bsln);
end
disp('Done!')
end