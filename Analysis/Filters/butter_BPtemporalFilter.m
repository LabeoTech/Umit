function outFile = butter_BPtemporalFilter(File, SaveFolder, varargin)
% This function applies a 4th-order Butterworth filter to the Time dimension of
% imaging data.


%%% Arguments parsing and validation %%%
p = inputParser;
% The input of the function must be a File , RawFolder or SaveFolder
addRequired(p,'File',@isfile)% For a file as input.
% Save folder:
addRequired(p, 'SaveFolder', @isfolder);
% Optional Parameters:
% opts structure
default_opts = struct('LowCutOffHz', .3, 'HighCutOffHz', 3);
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
% memmap DataFile and MetaDataFile:
[mDat, mDt] = mapDatFile(File);
Freq = mDt.Freq;
% Locate "T"ime dimension:
indxT = strcmp('T', mDt.dim_names);
% Initialize Butterworth filter Parameters:
if( opts.LowCutOffHz > 0 )
    UseLPFilt = 1;
    f = fdesign.lowpass('N,F3dB', 4, opts.LowCutOffHz, Freq);    %Fluo lower Freq
    lpass = design(f,'butter');
else
    UseLPFilt = 0;
end
if( opts.HighCutOffHz > 0 )
    UseHPFilt = 1;
    f = fdesign.lowpass('N,F3dB', 4, opts.HighCutOffHz, Freq);   %Fluo Higher Freq
    hpass = design(f,'butter');
else
    UseHPFilt = 0;
end
% Load data:
data = mDat.Data.data;
% Find NaNs and replace them with zeros:
idx_nan = isnan(data);
data(idx_nan) = 0;
% Filter:
szData = size(data);
data = reshape(data,[], szData(indxT));
for i = 1:size(data,1)
    Sig = double(data(i,:));
    if( UseLPFilt )
        LP = filtfilt(lpass.sosMatrix, lpass.ScaleValues, Sig')';
    else
        LP = ones(size(Signal));
    end
    if( UseHPFilt )
        HP = filtfilt(hpass.sosMatrix, hpass.ScaleValues, Sig')';
    else
        HP = Signal;
    end
    data(i,:) = HP - LP + mean(Sig, 'omitnan');
end
data = reshape(data,szData);
% Put NaNs back to data:
data(idx_nan) = NaN;
% SAVING DATA :
% Generate .DAT and .MAT file Paths:
datFile = fullfile(SaveFolder, Output);
% Save to .DAT file and create .MAT file with metaData:
save2Dat(datFile, data, mDt.dim_names);
outFile = Output;
end
