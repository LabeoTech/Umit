function outData = temporal_filter(data, metaData, varargin)
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

% Output:
% outData (3D numerical array) : filtered data from "File".

% Defaults:
default_Output = 'BPtemporalFilter.dat';  %#ok. This line is here just for Pipeline management.
default_opts = struct('LowCutOffHz', '1/120', 'HighCutOffHz', '0', 'Normalize', true);
% Some notes on the CutOff values:
% 1) The HighCutOffHz value of 0 will be translated as the Nyquist of the sample rate
% 2) For the LowCutOff, values equal or less than zero will give a low-passed signal at "HighCutOff".


%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'data',@(x) isnumeric(x) & ndims(x) == 3); % Validate if the input is a 3-D numerical matrix:
addRequired(p,'metaData', @(x) isa(x,'matlab.io.MatFile') | isstruct(x)); % MetaData associated to "data".
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
% Parse inputs:
parse(p,data, metaData, varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize Variables and remove inputParser object:
outData = p.Results.data;
metaData = p.Results.metaData;
opts = p.Results.opts;
clear p
%%%%
% Transform strings in opts to numbers:
eval(['opts.LowCutOffHz = ' opts.LowCutOffHz, ';'])
eval(['opts.HighCutOffHz= ' opts.HighCutOffHz, ';'])
% Validate if "data" is an Image Time Series:
errID = 'umIToolbox:temporal_filter:InvalidInput';
errMsg = 'Wrong Input Data type. Data must be an Image time series with dimensions "X", "Y" and "T".';
assert(all(ismember(metaData.dim_names,{'Y', 'X', 'T'})), errID, errMsg);

% Find NaNs and replace them with zeros:
idx_nan = isnan(outData);
outData(idx_nan) = 0;
% Run Temporal filter function
disp('Filtering data...')
outData = NormalisationFiltering(pwd, outData, opts.LowCutOffHz, opts.HighCutOffHz, ...
    opts.Normalize,0, metaData.Freq);
disp('Finished with temporal filter.')
% Put NaNs back to data:
outData(idx_nan) = NaN;
end