function outFile = cheby1_BPtemporalFilter(File, SaveFolder, varargin)

%%% Arguments parsing and validation %%%
p = inputParser;
% The input of the function must be a File , RawFolder or SaveFolder
addRequired(p,'File',@isfile)% For a file as input.
% Save folder:
addRequired(p, 'SaveFolder', @isfolder);
% Optional Parameters:
% opts structure
default_opts = struct('LowPassHz', .3, 'HighPassHz', 3);
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
% Output File:
default_Output = 'cheby1Filt.dat';
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
% Initialize Chebyshev Filter Parameters:
bp = [opts.LowPassHz opts.HighPassHz]/(Freq/2); % Divide LowPass and HighPass values by Nyquist.
[b,a] = cheby1(2,0.5,bp);
% Load data:
data = mDat.Data.data;
% Find NaNs and replace them with zeros:
idx_nan = isnan(data);
data(idx_nan) = 0;
% Filter:
szData = size(data);
data = reshape(data,[], szData(indxT));
for i = 1:size(data,1)
    data(i,:) = single(filtfilt(b,a,double(data(i,:))));
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


