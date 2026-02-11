function writeTrial_YXTE(fidOut, trialIdx, trialData, szYXTE, datatype)
% WRITETRIAL_YXTE Writes a single trial with E as the last dimension
%
% File layout on disk: [Y, X, T, E]
%
% Inputs:
%   fidOut     : file ID opened with 'w' or 'r+'
%   trialIdx   : event index (1-based)
%   trialData  : [Y,X,T] numeric array
%   szYXTE     : full dataset size [Y,X,T,E]
%   datatype   : e.g. 'single', 'double'

% Sizes
nY = szYXTE(1);
nX = szYXTE(2);
nT = szYXTE(3);

bytesPerElem = getByteSize(datatype);
nElemPerTrial = nY * nX * nT;

% Offset to beginning of this trial block
offset = (trialIdx-1) * nElemPerTrial * bytesPerElem;

% Seek and write in one shot (FAST)
fseek(fidOut, offset, 'bof');
fwrite(fidOut, trialData, datatype);
end
