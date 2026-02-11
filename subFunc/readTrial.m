function trialData = readTrial(fidIn, trialIdx, sz, datatype)
% READTRIAL Reads a single trial (event) from an interleaved E,Y,X,T .dat file
% Faster version using 'fread' with skip
%
% Inputs:
%   fidIn      : file ID of open .dat file
%   trialIdx   : index of the event to read (1-based)
%   sz         : [E,Y,X,T] size of the dataset
%   datatype   : data type string, e.g. 'single', 'double'
%
% Output:
%   trialData  : [Y,X,T] matrix for this event

nE = sz(1);
nY = sz(2);
nX = sz(3);
nT = sz(4);

nElemPerEvent = nY*nX*nT;
bytesPerElem  = getByteSize(datatype);

% Compute offset of first element for this trial
startOffset = (trialIdx-1) * bytesPerElem;

% Compute skip between consecutive elements of this event
skipBytes = (nE-1) * bytesPerElem;

% Move to first element of trial
fseek(fidIn, startOffset, 'bof');

% Read all elements of this trial with automatic skip
trialData = fread(fidIn, nElemPerEvent, ['*' datatype], skipBytes);

% Reshape to [Y,X,T]
trialData = reshape(trialData, [nY, nX, nT]);
end
