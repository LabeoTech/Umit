function outFile = temporal_filter(File, SaveFolder, varargin)
% This function applies a 4th-order Butterworth filter to the Time dimension of
% imaging data. This is achieved by the subtraction of two low-pass
% filters.



% The data must be an Image time series with dimensions
% {Y,X,T}.
% This function is a wrapper of the IOI library function
% "NormalisationFiltering.m".

%%% Arguments parsing and validation %%%
p = inputParser;
% The input of the function must be a File , RawFolder or SaveFolder
addRequired(p,'File',@isfile)% For a file as input.
% Save folder:
addRequired(p, 'SaveFolder', @isfolder);
% Optional Parameters:
% opts structure
default_opts = struct('LowCutOffHz', .0033, 'HighCutOffHz', -1, 'Normalize', true); 
% Some notes on the CutOff values:
    % 1) The HighCutOffHz value of -1 will be translated as the Nyquist of the sample rate
    % 2) Other than the case above, all cut-off values equal or less than
    % zero will disable the Low-pass filter. 
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
if opts.HighCutOffHz == -1
    opts.HighCutOffHz = mDt.Freq/2; % Translate "-1" to Nyquist.
end
% Load Data:
data = mDat.Data.data;
% Find NaNs and replace them with zeros:
idx_nan = isnan(data);
data(idx_nan) = 0;
% Run Temporal filter function
data = NormalisationFiltering(data, mDt.Freq, opts.LowCutOffHz, opts.HighCutOffHz, ...
    opts.Normalize);
% Put NaNs back to data:
data(idx_nan) = NaN;
% SAVING DATA :
% Generate .DAT and .MAT file Paths:
datFile = fullfile(SaveFolder, Output);
% Save to .DAT file and create .MAT file with metaData:
save2Dat(datFile, data, mDt.dim_names);
outFile = Output;
end
