function outData = run_HemoCorrection(SaveFolder, data, metaData, varargin)
% RUN_HEMOCORRECTION Apply hemodynamic correction to a fluorescence channel.
%
% This function supports two correction approaches:
%
%   1) LinearRegression
%      Pixel-wise linear regression of fluorescence onto one or more
%      reflectance channels.
%
%   2) Ratiometric
%      Single-reference correction based on normalized fluorescence and
%      normalized reference reflectance.
%
% The function supports both:
%   - standard mode (numeric fluorescence array in RAM)
%   - low-RAM mode (fluorescence provided as a .dat filename)
%
% Inputs
% ------
% SaveFolder : char
%   Folder containing fluorescence and reflectance files.
%
% data : numeric array or char
%   Either:
%       - numeric fluorescence data [Y, X, T]
%       - fluorescence .dat filename
%
% metaData : struct or matlab.io.MatFile
%   Metadata associated with the fluorescence data.
%
% opts (optional struct)
%   Structure with fields:
%       - Algorithm : 'LinearRegression' or 'Ratiometric'
%       - Red       : logical
%       - Green     : logical
%       - Amber     : logical
%       - Other     : custom channel filename (char), used only in
%                     Ratiometric mode
%
% Output
% ------
% outData :
%   - Standard mode: corrected fluorescence array
%   - Low-RAM mode : corrected fluorescence filename
%
% Notes
% -----
% In Ratiometric mode, exactly one reference channel must be selected.

% Defaults:
default_Output = 'hemoCorr_fluo.dat'; %#ok This line is here just for Pipeline management.
default_opts = struct( ...
    'Algorithm', 'LinearRegression', ...
    'Red', true, ...
    'Green', true, ...
    'Amber', true, ...
    'Other', '');

opts_values = struct( ... %#ok This is here only as a reference for PIPELINEMANAGER.m.
    'Algorithm', {{'LinearRegression','Ratiometric'}}, ...
    'Red', [false, true], ...
    'Green', [false, true], ...
    'Amber', [false, true], ...
    'Other', {{''}});

% -------------------------------------------------------------------------
% Input parsing
% -------------------------------------------------------------------------
p = inputParser;
addRequired(p, 'SaveFolder', @isfolder);
addRequired(p, 'data', @(x) (isnumeric(x) && ndims(x) == 3) || ischar(x));
addRequired(p, 'metaData', @(x) isa(x,'matlab.io.MatFile') || isstruct(x));
addOptional(p, 'opts', default_opts, @(x) isstruct(x) && ~isempty(x));
parse(p, SaveFolder, data, metaData, varargin{:});

SaveFolder = p.Results.SaveFolder;
data = p.Results.data;
metaData = p.Results.metaData;
opts = p.Results.opts;
clear p

% -------------------------------------------------------------------------
% Build selected channel list
% -------------------------------------------------------------------------
list = {};
if opts.Red
    list{end+1} = 'Red'; %#ok<AGROW>
end
if opts.Green
    list{end+1} = 'Green'; %#ok<AGROW>
end
if opts.Amber
    list{end+1} = 'Amber'; %#ok<AGROW>
end
if ~isempty(opts.Other)
    list{end+1} = opts.Other; %#ok<AGROW>
end

fprintf('Performing hemodynamic correction in fluo channel using %s algorithm...\n', ...
    opts.Algorithm);

% -------------------------------------------------------------------------
% Dispatch by algorithm
% -------------------------------------------------------------------------
switch lower(opts.Algorithm)
    case 'linearregression'
        outData = HemoCorrection(SaveFolder, data, metaData, list);

    case 'ratiometric'
        assert(numel(list) == 1, ...
            'umIToolbox:run_HemoCorrection:InvalidInput', ...
            'Ratiometric correction requires exactly one reference channel.');

        refFile = localResolveReferenceFile(SaveFolder, list{1});
        fprintf('Using channel "%s" in hemodynamic correction...\n', refFile);

        if ischar(data)
            outData = localRatiometricLowRAM(SaveFolder, data, metaData, refFile);
        else
            outData = localRatiometricStandard(SaveFolder, data, metaData, refFile);
        end

    otherwise
        error('umIToolbox:run_HemoCorrection:InvalidInput', ...
            'Unknown correction algorithm "%s".', opts.Algorithm);
end

disp('Finished hemodynamic correction.')
end


% =========================================================================
% Local functions
% =========================================================================
function refFile = localResolveReferenceFile(SaveFolder, channelTag)
% Resolve channel tag or filename to an existing .dat file in SaveFolder.

tag = lower(channelTag);

switch tag
    case 'red'
        if isfile(fullfile(SaveFolder, 'rChan.dat'))
            refFile = 'rChan.dat';
        else
            refFile = 'red.dat';
        end

    case {'amber', 'yellow'}
        if isfile(fullfile(SaveFolder, 'yChan.dat'))
            refFile = 'yChan.dat';
        else
            refFile = 'yellow.dat';
        end

    case 'green'
        if isfile(fullfile(SaveFolder, 'gChan.dat'))
            refFile = 'gChan.dat';
        else
            refFile = 'green.dat';
        end

    otherwise
        [~, name, ext] = fileparts(channelTag);
        if isempty(ext)
            ext = '.dat';
        end
        refFile = [name, ext];
end

assert(isfile(fullfile(SaveFolder, refFile)), ...
    'umIToolbox:run_HemoCorrection:FileNotFound', ...
    'Reference channel file "%s" not found.', refFile);
end


function outData = localRatiometricStandard(SaveFolder, data, metaData, refFile)
% Ratiometric correction in standard in-memory mode.

refData = loadDatFile(fullfile(SaveFolder, refFile));

assert(all(size(data) == size(refData)), ...
    'umIToolbox:run_HemoCorrection:InvalidInput', ...
    'Reference and fluorescence channels must have the same size.');

datSize = size(data);

% Fluorescence normalization
data2D = reshape(data, [], datSize(end));
mData = mean(data2D, 2, 'omitnan');
data2D = (data2D - mData) ./ mData;

% Reference normalization
ref2D = reshape(refData, [], datSize(end));
mRef = mean(ref2D, 2, 'omitnan');
ref2D = (ref2D - mRef) ./ mRef;

% Ratiometric correction
data2D = data2D - ref2D;

% Restore fluorescence mean
data2D = (data2D .* mData) + mData;

outData = reshape(data2D, datSize);
end


function outFile = localRatiometricLowRAM(SaveFolder, fluoFile, metaData, refFile)
% Ratiometric correction in low-RAM mode using spatial slabs.

outFile = fullfile(SaveFolder, 'hemoCorr_fluo.dat');

Ny = metaData.datSize(1);
Nx = metaData.datSize(2);
Nt = metaData.datLength;

fidFluo = fopen(fullfile(SaveFolder, fluoFile), 'r');
assert(fidFluo ~= -1, ...
    'umIToolbox:run_HemoCorrection:FileOpenError', ...
    'Could not open fluorescence file "%s".', fluoFile);
cFluo = onCleanup(@() safeFclose(fidFluo)); %#ok<NASGU>

fidRef = fopen(fullfile(SaveFolder, refFile), 'r');
assert(fidRef ~= -1, ...
    'umIToolbox:run_HemoCorrection:FileOpenError', ...
    'Could not open reference file "%s".', refFile);
cRef = onCleanup(@() safeFclose(fidRef)); %#ok<NASGU>

preallocateDatFile(outFile, metaData);

fidOut = fopen(outFile, 'r+');
assert(fidOut ~= -1, ...
    'umIToolbox:run_HemoCorrection:FileOpenError', ...
    'Could not create output file "%s".', outFile);
cOut = onCleanup(@() safeFclose(fidOut)); %#ok<NASGU>

% Chunk along X to control RAM
dataBytes = prod([metaData.datSize, metaData.datLength, getByteSize(metaData.Datatype)]);
nChunks = calculateMaxChunkSize(dataBytes, 3, .1);
chunkX = ceil(Nx / nChunks);

for ii = 1:nChunks
    xStart = (ii - 1) * chunkX + 1;
    xEnd   = min(xStart + chunkX - 1, Nx);
    xIdx   = xStart:xEnd;

    % Read slabs
    fSlab = spatialSlabIO('read', fidFluo, Ny, Nx, Nt, xIdx, metaData.Datatype);
    rSlab = spatialSlabIO('read', fidRef,  Ny, Nx, Nt, xIdx, metaData.Datatype);

    assert(all(size(fSlab) == size(rSlab)), ...
        'umIToolbox:run_HemoCorrection:InvalidInput', ...
        'Reference and fluorescence slabs must have the same size.');

    slabSz = size(fSlab);

    % Reshape to [pixels x time]
    fSlab = reshape(fSlab, [], Nt);
    rSlab = reshape(rSlab, [], Nt);

    % Normalize each pixel trace
    mFluo = mean(fSlab, 2, 'omitnan');
    fSlab = (fSlab - mFluo) ./ mFluo;

    mRef = mean(rSlab, 2, 'omitnan');
    rSlab = (rSlab - mRef) ./ mRef;

    % Ratiometric correction
    fSlab = fSlab - rSlab;

    % Restore fluorescence mean
    fSlab = (fSlab .* mFluo) + mFluo;

    % Write corrected slab
    fSlab = reshape(fSlab, slabSz);
    spatialSlabIO('write', fidOut, Ny, Nx, Nt, xIdx, metaData.Datatype, fSlab);
end

[~, outName, outExt] = fileparts(outFile);
outFile = [outName outExt];
end