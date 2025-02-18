function outData = normalizeLPF(data, SaveFolder, varargin)
% NORMALIZELPF performs a data normalization by low pass filtering.
% The filtering algorithm consists in creating two low-passed versions of the
% signal with a given cut-off frequency ("LowCutOff" and "HighCutOff") and
% subsequently subtracting the two filtered signals.
% In this case, the signal(R) is expressed as DeltaR.
% Optionally, the subtracted signals can be normalized to the low cut-off signal
% to express the signal as DeltaR/R.
% This function is a wrapper of the IOI library function "NormalisationFiltering.m".
% For more information on the algorithm, refer to the function's documentation.
%
% Limitations:
% The data must be an Image time series with dimensions {Y,X,T} or split by events ('E','Y','X','T').
%
% Inputs:
%   data: numerical matrix containing image time series (with dimensions "Y", "X", "T" or "E","Y", "X", "T").
%   SaveFolder (char): folder containing "data" and the associated AcqInfos.mat file.
%   opts (optional): structure containing extra parameters. See "default_opts" variable below for details!
% Outputs:
%   outData: numerical matrix with dimensions {Y,X,T}.

% Defaults:
default_Output = 'normLPF.dat';  %#ok. This line is here just for Pipeline management.
default_opts = struct('LowCutOffHz', 0.0083, 'HighCutOffHz', 1, 'Normalize', true, 'bApplyExpFit', false);
opts_values = struct('LowCutOffHz', [0,Inf], 'HighCutOffHz',[eps,Inf],'Normalize',[false, true],'bApplyExpFit', [true,false]);%#ok  % This is here only as a reference for PIPELINEMANAGER.m.
% Some notes on the CutOff values:
% 1) The HighCutOffHz value of 0 will be translated as the Nyquist of the sample rate
% 2) For the LowCutOff, values equal to zero will give a low-passed signal at "HighCutOff".


%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'data',@(x) (isnumeric(x) && ndims(x) == 3) || isstruct(x)); % Validate if the input is numerical or a structure.
addRequired(p,'SaveFolder',@(x) isfolder(x)); % Validate if the SaveFolder exists.
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
% Parse inputs:
parse(p,data, SaveFolder, varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize Variables and remove inputParser object:
outData = data;
opts = p.Results.opts;
clear p data
%%%%
if isstruct(outData)
    Freq = outData.FrameRateHz;
else
    info = load(fullfile(SaveFolder,'AcqInfos.mat'));
    Freq = info.AcqInfoStream.FrameRateHz;
end    
% Check if cut-off frequencies are in the acceptable range:
errID = 'umIToolbox:normalizeLPF:InvalidInput';
% Check Low cut-off frequency:
if opts.LowCutOffHz < 0 || opts.LowCutOffHz > Freq/2
    error(errID,['Invalid cut off value! LowCutOffHz must be between 0 and ' num2str(Freq/2) '!'])
end
% Check High cut-off frequency:
if opts.HighCutOffHz < opts.LowCutOffHz
    error(errID,'Invalid cut off value! HighCutOffHz must be higher than LowCutOffHz!');
end

if opts.HighCutOffHz == 0 ||  opts.HighCutOffHz > Freq/2
    error(errID,['Invalid cut off value! HighCutOffHz must be a positive number less than or equal to ' num2str(Freq/2) '!']);
end

% Process data:
if isstruct(outData)
    % For data from .DATASTAT files
    % Check if the data is an Image time series
    errID = 'umIToolbox:normalizeLPF:InvalidInput';
    errMsg = 'Wrong Input Data type. Data must be an Image time series.';
    assert(strcmpi(outData.dataCategory.data,'Image-time-series'),errID,errMsg);    
    if outData.b_hasEvents
        % The data is split by events       
        for ii = 1:size(outData.data.data,1)
            disp(['Filtering Event #' num2str(ii) '/' size(outData.data.data,1) '...']);
            outData.data.data(ii,:,:,:) = NormalisationFiltering(pwd, squeeze(outData.data.data(ii,:,:,:)), opts.LowCutOffHz, opts.HighCutOffHz, ...
                opts.Normalize,opts.bApplyExpFit, Freq);
        end
    else
        % The data is an image-time-series.
        outData.data.data = NormalisationFiltering(SaveFolder, outData.data.data, opts.LowCutOffHz, opts.HighCutOffHz, ...
        opts.Normalize,opts.bApplyExpFit, Freq);
    end            
else
    % For Image time series from .DAT files:        
    outData = NormalisationFiltering(SaveFolder, outData, opts.LowCutOffHz, opts.HighCutOffHz, ...
        opts.Normalize,opts.bApplyExpFit, Freq);
    disp('Finished with normalization filtering.')
end

disp('Finished with temporal filter.')

end