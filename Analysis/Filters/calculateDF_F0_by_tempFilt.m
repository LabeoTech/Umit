function outFile = calculateDF_F0_by_tempFilt(File, SaveFolder, varargin)
% CALCULATEDF_F0_BY_TEMPFILT normalizes imaging data using temporal filtering method.


%%% Arguments parsing and validation %%%
p = inputParser;
% The input of the function must be a File , RawFolder or SaveFolder
addRequired(p,'File',@isfile)% For a file as input.
% Save folder:
addRequired(p, 'SaveFolder', @isfolder);
% Optional Parameters:
% opts structure
default_opts = struct('HighPassHz', -1); % This default value will be translated as the Nyquist of the sample rate
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
% Output File:
default_Output = 'deltaF_F0_by_tempFilt.dat';
addOptional(p, 'Output', default_Output, @(x) ischar(x) || isstring(x) || iscell(x));
% Parse inputs:
parse(p,File, SaveFolder, varargin{:});
%Initialize Variables:
File = p.Results.File;
SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;
Output = p.Results.Output;
%%%%
% Open memMapfile and metaData:
mmData = mapDatFile(File);
mDt = matfile(strrep(File, '.dat','_info.mat'));
% Load Data and Sample Rate:
data = mmData.Data.data;
if (mean(reshape(data,[],size(data,3)),1) < 0.5)
    data = data + 1;
end
Freq = mDt.Freq;
if opts.HighPassHz == -1
    opts.HighPassHz = Freq/2; % Translate "-1" to Nyquist.
end
% Set temporal Filter params:
szData = size(data);
f = fdesign.lowpass('N,F3dB', 4, 1/120, Freq);
lpass = design(f,'butter'); lp_sosM = lpass.sosMatrix; lp_SV = lpass.ScaleValues;
f = fdesign.lowpass('N,F3dB', 4, opts.HighPassHz, Freq);
hpass = design(f,'butter'); hp_sosM = hpass.sosMatrix; hp_SV = hpass.ScaleValues;
for ind = 1:szData(1)
    Sig = double(squeeze(data(ind,:,:)))';
    % Apply temporal filter:
    dl = single(filtfilt(lp_sosM, lp_SV, Sig));
    dh = single(filtfilt(hp_sosM, hp_SV, Sig));
    data(ind,:,:) = (dh./dl)';
end
% Replace NaNs with zeros:
idx = isnan(data);
data(idx) = 0;
% SAVING DATA :
% Generate .DAT and .MAT file Paths:
datFile = fullfile(SaveFolder, Output);
% Save to .DAT file and create .MAT file with metaData:
save2Dat(datFile, data);
outFile = Output;
end


