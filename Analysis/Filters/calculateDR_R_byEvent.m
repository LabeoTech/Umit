function outData = calculateDR_R_byEvent(data, metaData)
% CALCULATEDR_R_BYEVENT generates DeltaR/R0 movie from DATA.
% The data must be a 4-D matrix containing image time series separated by
% trials.
% Inputs:
%   data (4D numerical matrix): Image time series separated by events with dimensions {'E','Y','X','T}.
%   metaData (struct): structure containing smeta data associated with "data".
% Output:
%   outData(4D numerical matrix): "data" with values transformed to express
%   DeltaR/R0.
% Defaults:
default_Output = 'deltaR_R_ByEvent.dat'; %#ok This is here only as a reference for PIPELINEMANAGER.m. 
%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'data',@(x) isnumeric(x)); % Validate if the input is a numerical matrix:
addRequired(p,'metaData', @(x) isa(x,'matlab.io.MatFile') | isstruct(x)); % MetaData associated to "data".
% Parse inputs:
parse(p,data, metaData);
%Initialize Variables:
outData = p.Results.data; 
metaData = p.Results.metaData;
clear p
%%%%
% Validate if the data is an Image time series by Events:
errMsg = 'Input Data must be an Image time series separated by events with dimensions {"E","Y","X",T"}.';
errID = 'umIToolbox:calculateDF_F0_byEvent:WrongInput';
assert(isequal(metaData.dim_names, {'E','Y','X','T'}), errID, errMsg)
% Calculate baseline:
bsln = median(outData(:,:,:,1:round(metaData.preEventTime_sec*metaData.Freq)), ...
    4,'omitnan');
% Normalize data to get DeltaR/R values
disp('Calculating DeltaR/R ...');
outData = (outData - bsln)./bsln;
disp('Finished with DeltaR/R!')
end