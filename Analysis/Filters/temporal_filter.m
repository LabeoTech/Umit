function outFile = temporal_filter(File, SaveFolder, varargin)
% This function applies a band-pass 4th-order Butterworth filter to the Time
% dimension of an image time series (Y,X,T) dataset. 
% The filtering algorithm consists in creating two low-passed versions of the 
% signal with a given cut-off frequency ("LowCutOff" and "HighCutOff") and 
% subsequently subtracting the two filtered signals. 
% In this case, the signal(R) is expressed as DeltaR.
% Optionally, the subtracted signals can be normalized to the low cut-off signal 
% to express the signal as DeltaR/R.

% Limitations:
% The data must be an Image time series with dimensions
% {Y,X,T}.

% This function is a wrapper of the IOI library function "NormalisationFiltering.m".
% For more information on the algorithm, refer to the function's documentation.





%%% Arguments parsing and validation %%%
p = inputParser;
% The input of the function must be a File , RawFolder or SaveFolder
addRequired(p,'File',@isfile)% For a file as input.
% Save folder:
addRequired(p, 'SaveFolder', @isfolder);
% Optional Parameters:
% opts structure
default_opts = struct('LowCutOffHz', .0033, 'HighCutOffHz', 0, 'Normalize', true); 
% Some notes on the CutOff values:
    % 1) The HighCutOffHz value of 0 will be translated as the Nyquist of the sample rate
    % 2) For the LowCutOff, values equal or less than zero will give a low-passed signal at "HighCutOff". 
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
% Output File:
default_Output = 'BPtemporalFilter.dat';
addOptional(p, 'Output', default_Output, @(x) ischar(x) || isstring(x) || iscell(x));
% Parse inputs:
parse(p,File, SaveFolder, varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialize Variables:
File = p.Results.File;
SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;
Output = p.Results.Output;
%%%%
% memory map DataFile and MetaDataFile:
[mDat, mDt] = mapDatFile(File);
% Validata if Data is an Image Time Series:
if ~all(ismember(mDt.dim_names, {'Y', 'X', 'T'}))
    error('umIToolbox:temporal_filter:InvalidInput', ...
        'Wrong Input Data type. Data must be an Image time series with dimensions "X", "Y" and "T".');
end
if opts.HighCutOffHz == -1
    opts.HighCutOffHz = mDt.Freq/2; % Translate "-1" to Nyquist.
end
% Load Data:
data = mDat.Data.data;
% Find NaNs and replace them with zeros:
idx_nan = isnan(data);
data(idx_nan) = 0;
% Run Temporal filter function
data = NormalisationFiltering(fileparts(File), data, opts.LowCutOffHz, opts.HighCutOffHz, ...
    opts.Normalize, mDt.Freq);
% Put NaNs back to data:
data(idx_nan) = NaN;
% SAVING DATA :
% Generate .DAT and .MAT file Paths:
datFile = fullfile(SaveFolder, Output);
% Save to .DAT file and create .MAT file with metaData:
save2Dat(datFile, data, mDt.dim_names);
outFile = Output;
end
