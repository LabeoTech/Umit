function [outData, metaData] = trim_movie(data, metaData, varargin)
% TRIM_MOVIE crops an image time series by N seconds or frames
% at the beginning and/or end of the movie.
%
% Supports automatic low-RAM mode when DATA is a filename (char).

% Defaults
default_Output = 'cropped_mov.dat'; %#ok
default_opts = struct('crop_start', 0, 'crop_end', 0, 'TimeUnit', 'Frames');
opts_values = struct('crop_start', [0 Inf], 'crop_end', [0 Inf],'TimeUnit', {{'Seconds','Frames'}}); %#ok

%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'data',@(x) isnumeric(x) || ischar(x));
addRequired(p,'metaData', @(x) isa(x,'matlab.io.MatFile') || isstruct(x));
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
parse(p,data, metaData, varargin{:});

outData  = p.Results.data;
metaData = p.Results.metaData;
opts     = p.Results.opts;
clear p

dim_names = metaData.dim_names;

%%% Validation %%%
errID  = 'umIToolbox:trim_movie:InvalidDataType';
errMsg = 'Input must be a 3D matrix with dimension "T"';
assert(all(ismember(dim_names, {'Y','X','T'})), errID, errMsg);

assert(opts.crop_start >= 0, errID, 'crop_start must be >= 0');
assert(opts.crop_end   >= 0, errID, 'crop_end must be >= 0');

%%% Identify T dimension %%%
idxT = find(strcmp('T', dim_names));

%%% Determine total frames %%%
if isnumeric(outData)
    nT = size(outData, idxT);
else
    sz = [metaData.datSize, metaData.datLength];
    nT = sz(idxT);
end

%%% Convert crop units to frames %%%
if startsWith(lower(opts.TimeUnit),'sec')
    fr_start = round(metaData.Freq * opts.crop_start);
    fr_stop  = round(nT - metaData.Freq * opts.crop_end);
else
    fr_start = opts.crop_start;
    fr_stop  = nT - opts.crop_end;
end

fr_start = max(fr_start, 1);
fr_stop  = min(fr_stop, nT);

assert(fr_start <= fr_stop, errID, ...
    sprintf('Invalid crop range [%d %d]', fr_start, fr_stop));

%% ==========================================================
% STANDARD MODE (numeric input) — UNCHANGED
% ==========================================================
if isnumeric(outData)

    orig_dim_indx = 1:numel(dim_names);
    new_dim_indx  = [idxT setdiff(orig_dim_indx, idxT)];

    outData = permute(outData, new_dim_indx);

    data_sz = size(outData);
    outData = reshape(outData, data_sz(1), []);

    disp('Cropping movie...')
    new_sz    = data_sz;
    new_sz(1)= fr_stop - fr_start + 1;

    outData = outData(fr_start:fr_stop, :);
    outData = reshape(outData, new_sz);
    outData = permute(outData, [2:numel(dim_names) 1]);

    metaData = genMetaData(outData, dim_names, metaData);
    disp('Done')
    return
end

%% ==========================================================
% RAM-SAFE MODE (disk-backed input)
% ==========================================================
disp('Cropping movie (RAM-Safe mode)...')

% Open input file
fidIn = fopen(outData, 'r');
cIn = onCleanup(@() safeFclose(fidIn));
assert(fidIn > 0, errID, 'Failed to open input file');

% Output file
outFile = fullfile(fileparts(data),'TRIMMED_DATA.dat');
fidOut = fopen(outFile, 'w');
cOut = onCleanup(@() safeFclose(fidOut));
assert(fidOut > 0, errID, 'Failed to create output file');

% Data geometry

nY = sz(1);
nX = sz(2);

frameStride  = nY * nX * 4;

% Seek to first kept frame
fseek(fidIn, (fr_start-1) * frameStride, 'bof');

% Stream frames
for t = fr_start:fr_stop
    frame = fread(fidIn, nY*nX, ['*' metaData.Datatype]);
    fwrite(fidOut, frame, metaData.Datatype);
end

fclose(fidIn);
fclose(fidOut);

% Update metadata
metaData.datLength = fr_stop - fr_start + 1;
outData = outFile;
save(strrep(outFile,'.dat','.mat'), '-struct','metaData')
disp('Done')
end
