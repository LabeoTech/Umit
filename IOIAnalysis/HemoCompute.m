function [HbO, HbR] = HemoCompute(DataFolder, SaveFolder, FilterSet, Illumination, b_normalize, varargin)
%HEMOCOMPUTE Approximate HbO and HbR concentration changes from intrinsic imaging.
%
%   [HbO, HbR] = HemoCompute(DataFolder, SaveFolder, FilterSet, Illumination, b_normalize)
%   [HbO, HbR] = HemoCompute(..., 'HbT_uM', 100, 'O2_sat', 60, 'RAMSafeMode', false)
%
%   This function approximates concentration changes of oxygenated (HbO)
%   and deoxygenated (HbR) hemoglobin from two or three intrinsic imaging
%   wavelengths.
%
%   Inputs:
%       DataFolder    - Folder containing red.dat, green.dat, and/or yellow.dat.
%       SaveFolder    - Folder where HbO/HbR outputs are saved. If empty,
%                       outputs are returned only in RAM in standard mode.
%       FilterSet     - 'gCaMP', 'jrGECO', or 'none'.
%       Illumination  - Cell array containing at least two wavelengths from
%                       {'red','green','yellow','amber'}.
%       b_normalize   - Logical scalar. If true, normalize input channels
%                       when needed.
%
%   Name-Value parameters:
%       'HbT_uM'      - Total hemoglobin concentration in micromolar.
%                       Default = 100.
%       'O2_sat'      - Oxygen saturation percentage. Default = 60.
%       'RAMSafeMode' - Logical scalar. If true, write HbO/HbR directly to
%                       disk instead of holding them in RAM. Default = false.
%
%   Outputs:
%       HbO           - In standard mode, numeric HbO array. In RAM-safe
%                       mode, output filename.
%       HbR           - In standard mode, numeric HbR array. In RAM-safe
%                       mode, output filename.

% Default output names for pipeline management.
default_Output_HbO = 'HbO.dat'; %#ok<NASGU>
default_Output_HbR = 'HbR.dat'; %#ok<NASGU>

if nargin == 1 && (ischar(DataFolder) || (isstring(DataFolder) && isscalar(DataFolder))) && ...
        strcmpi(strtrim(char(string(DataFolder))), 'pipelineInfo')
    HbO = localPipelineInfo();
    HbR = [];
    return
end

p = inputParser;
p.FunctionName = 'HemoCompute';
addRequired(p, 'DataFolder', @(x) (ischar(x) || (isstring(x) && isscalar(x))) && isfolder(x));
addRequired(p, 'SaveFolder', @(x) isempty(x) || isfolder(x) || ischar(x) || (isstring(x) && isscalar(x)));
addRequired(p, 'FilterSet', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addRequired(p, 'Illumination', @(x) iscell(x) && ~isempty(x));
addRequired(p, 'b_normalize', @(x) islogical(x) && isscalar(x));
addParameter(p, 'HbT_uM', 100, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'O2_sat', 60, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && x <= 100);
addParameter(p, 'RAMSafeMode', false, @(x) islogical(x) && isscalar(x));
parse(p, DataFolder, SaveFolder, FilterSet, Illumination, b_normalize, varargin{:});

DataFolder = char(string(p.Results.DataFolder));
SaveFolder = char(string(p.Results.SaveFolder));
FilterSet = char(string(p.Results.FilterSet));
Illumination = p.Results.Illumination;
b_normalize = p.Results.b_normalize;
HbT_uM = double(p.Results.HbT_uM);
O2_sat = double(p.Results.O2_sat);
b_RAMsafeMode = p.Results.RAMSafeMode;

bSave = ~isempty(SaveFolder);
if bSave
    assert(isfolder(SaveFolder), ...
        'Umitoolbox:HemoCompute:InvalidSaveFolder', ...
        'SaveFolder must be an existing folder or empty.');
end

assert(ismember(lower(FilterSet), {'gcamp','jrgeco','none'}), ...
    'Umitoolbox:HemoCompute:InvalidFilterSet', ...
    'FilterSet must be one of: gCaMP, jrGECO, none.');

% Keep this exact normalization of amber to yellow.
idx = contains(lower(Illumination), 'amber');
if any(idx)
    Illumination{idx} = 'yellow';
end

illuminationLower = lower(string(Illumination));
validIllum = ismember(illuminationLower, {'red','green','yellow'});
assert(all(validIllum), ...
    'Umitoolbox:HemoCompute:InvalidIllumination', ...
    'Illumination entries must be red, green, yellow, or amber.');
assert(numel(unique(illuminationLower)) >= 2, ...
    'Umitoolbox:HemoCompute:InvalidIllumination', ...
    'At least two different illumination wavelengths are needed for Hb computation.');

acqFile = fullfile(DataFolder, 'AcqInfos.mat');
assert(isfile(acqFile), ...
    'Umitoolbox:HemoCompute:MissingAcqInfos', ...
    'AcqInfos.mat was not found in "%s".', DataFolder);

infos = load(acqFile, 'AcqInfoStream');
assert(isfield(infos, 'AcqInfoStream') && isstruct(infos.AcqInfoStream), ...
    'Umitoolbox:HemoCompute:InvalidAcqInfos', ...
    'AcqInfos.mat does not contain a valid AcqInfoStream structure.');
acq = infos.AcqInfoStream;

assert(isfield(acq, 'Height') && isfield(acq, 'Width') && isfield(acq, 'Length') && ...
    isfield(acq, 'FrameRateHz') && isfield(acq, 'Camera_Model'), ...
    'Umitoolbox:HemoCompute:InvalidAcqInfos', ...
    'AcqInfoStream must contain Height, Width, Length, FrameRateHz, and Camera_Model.');

Ny = double(acq.Height);
Nx = double(acq.Width);
Nt = double(acq.Length);
Freq = double(acq.FrameRateHz);
Datatype = 'single'; %#ok<NASGU>

% File opening and availability checks.
fidR = 0;
fidG = 0;
fidY = 0;

if any(strcmp(illuminationLower, 'red'))
    redFile = fullfile(DataFolder, 'red.dat');
    assert(isfile(redFile), ...
        'Umitoolbox:HemoCompute:MissingChannelFile', ...
        'Required channel file was not found: %s', redFile);
    fidR = fopen(redFile, 'r');
    assert(fidR ~= -1, ...
        'Umitoolbox:HemoCompute:FileOpenError', ...
        'Could not open channel file "%s".', redFile);
    c_r = onCleanup(@() safeFclose(fidR)); %#ok<NASGU>
end

if any(strcmp(illuminationLower, 'green'))
    greenFile = fullfile(DataFolder, 'green.dat');
    assert(isfile(greenFile), ...
        'Umitoolbox:HemoCompute:MissingChannelFile', ...
        'Required channel file was not found: %s', greenFile);
    fidG = fopen(greenFile, 'r');
    assert(fidG ~= -1, ...
        'Umitoolbox:HemoCompute:FileOpenError', ...
        'Could not open channel file "%s".', greenFile);
    c_g = onCleanup(@() safeFclose(fidG)); %#ok<NASGU>
end

if any(strcmp(illuminationLower, 'yellow'))
    yellowFile = fullfile(DataFolder, 'yellow.dat');
    assert(isfile(yellowFile), ...
        'Umitoolbox:HemoCompute:MissingChannelFile', ...
        'Required channel file was not found: %s', yellowFile);
    fidY = fopen(yellowFile, 'r');
    assert(fidY ~= -1, ...
        'Umitoolbox:HemoCompute:FileOpenError', ...
        'Could not open channel file "%s".', yellowFile);
    c_y = onCleanup(@() safeFclose(fidY)); %#ok<NASGU>
end

% Data size is strictly continuous YXT.
datsz = [Ny, Nx, Nt];
NbPix = [Ny, Nx];

% Check whether the data are already normalized.
indxNorm = [-2 -2 -2];
fTags = {'fidR', 'fidG', 'fidY'};

fprintf('Checking channel data...\n');
fidList = [fidR, fidG, fidY];
bytesPerFrame = prod(datsz(1:2)) * 4;

fprintf('Checking channel data...\n');
for i = 1:3
    thisFid = fidList(i);
    if thisFid == 0
        continue
    end

    % Sample 10 frames to estimate normalization state.
    frIdx = unique(floor(linspace(1, datsz(3), 10)));
    frIdx(frIdx < 1) = [];
    tmp = zeros(datsz(1), datsz(2), 'single');

    for jj = 1:numel(frIdx)
        fseek(thisFid, (frIdx(jj) - 1) * bytesPerFrame, 'bof');
        frame = fread(thisFid, datsz([1 2]), '*single');
        tmp = tmp + frame';
    end

    tmp = tmp ./ numel(frIdx);
    Mdat = mean(tmp, 'all', 'omitnan');

    if Mdat > .75 && Mdat < 1.25
        indxNorm(i) = 1;
    elseif Mdat > -.25 && Mdat < .25
        indxNorm(i) = 0;
    else
        indxNorm(i) = -1;
    end
end

if ~all(indxNorm)
    error('Umitoolbox:HemoCompute:HeterogeneousInput', ...
        'The input data is heterogeneous. All channels must be preprocessed in the same way.');
end

indxNorm = indxNorm(find(indxNorm ~= -2, 1, 'first'));

if indxNorm == -1 && ~b_normalize
    error('Umitoolbox:HemoCompute:NeedsNormalization', ...
        'The channels must be normalized or b_normalize must be true.');
end

if indxNorm == 0
    warning('Umitoolbox:HemoCompute:CenteredAtZero', ...
        'The channels are centered at zero. They will be shifted to be centered at one.');
end

if b_normalize && (indxNorm == 1 || indxNorm == 0)
    b_normalize = false;
    warning('Umitoolbox:HemoCompute:AlreadyNormalized', ...
        'The input data is already normalized. Normalization will be skipped.');
end

fprintf('Data checked!\n');

% Compute extinction/pathlength model.
A = ioi_epsilon_pathlength( ...
    'Hillman', ...
    HbT_uM, ...
    O2_sat, ...
    100 - O2_sat, ...
    FilterSet, ...
    acq.Camera_Model);

% Temporal filters used for input normalization.
f = fdesign.lowpass('N,F3dB', 4, 1, Freq);
lpass_high = design(f, 'butter');

f = fdesign.lowpass('N,F3dB', 4, 1/120, Freq);
lpass_low = design(f, 'butter');

% Chunking strategy for continuous YXT data only.
nChunks = calculateMaxChunkSize(prod(datsz) * 4, 12, .1);
chunkX  = ceil(NbPix(2) / nChunks);

% Output allocation / preallocation.
if ~b_RAMsafeMode
    HbO = zeros(datsz, 'single');
    HbR = zeros(datsz, 'single');
else
    assert(bSave, ...
        'Umitoolbox:HemoCompute:MissingSaveFolder', ...
        'SaveFolder must be provided when RAMSafeMode is true.');

    hboPath = fullfile(SaveFolder, default_Output_HbO);
    if isfile(hboPath)
        [folderPath, baseName, ext] = fileparts(hboPath);
        hboPath = fullfile(folderPath, [baseName '_preallocData' ext]);
    end

    hbrPath = fullfile(SaveFolder, default_Output_HbR);
    if isfile(hbrPath)
        [folderPath, baseName, ext] = fileparts(hbrPath);
        hbrPath = fullfile(folderPath, [baseName '_preallocData' ext]);
    end

    preallocateDatFile(hboPath, [Ny, Nx, Nt], 'single');
    fid_hbo = fopen(hboPath, 'r+');
    c_hbo = onCleanup(@() safeFclose(fid_hbo)); %#ok<NASGU>

    preallocateDatFile(hbrPath, [Ny, Nx, Nt], 'single');
    fid_hbr = fopen(hbrPath, 'r+');
    c_hbr = onCleanup(@() safeFclose(fid_hbr)); %#ok<NASGU>
end

% Computation loop.
h = waitbar(0, 'Computing');
for indP = 1:nChunks
    xStart = (indP - 1) * chunkX + 1;
    xEnd   = min(xStart + chunkX - 1, NbPix(2));
    xIdx   = xStart:xEnd;

    if b_RAMsafeMode
        h.Name = ['HemoCompute (chunk ' num2str(indP) '/' num2str(nChunks) ')'];
        drawnow()
    end

    if fidR
        waitbar(indP/nChunks, h, 'Red channel [Reading file...]')
        Red = spatialSlabIO('read', fidR, NbPix(1), NbPix(2), Nt, xIdx, 'single');
        Red = reshape(Red, [], Nt);

        if b_normalize
            waitbar(indP/nChunks, h, 'Red channel [Normalizing data...]')
            Red = single(filtfilt(lpass_high.sosMatrix, lpass_high.ScaleValues, double(Red)'))';
            tmp = single(filtfilt(lpass_low.sosMatrix, lpass_low.ScaleValues, double(Red)'))';
            tmp(tmp < min(Red(:))) = min(Red(:));
            Red = Red ./ tmp;
        end

        if indxNorm == 0
            Red = Red + 1;
        end
        Red = -log(Red);
    end

    if fidG
        waitbar(indP/nChunks, h, 'Green channel [Reading file...]')
        Green = spatialSlabIO('read', fidG, NbPix(1), NbPix(2), Nt, xIdx, 'single');
        Green = reshape(Green, [], Nt);

        if b_normalize
            waitbar(indP/nChunks, h, 'Green channel [Normalizing data...]')
            Green = single(filtfilt(lpass_high.sosMatrix, lpass_high.ScaleValues, double(Green)'))';
            tmp = single(filtfilt(lpass_low.sosMatrix, lpass_low.ScaleValues, double(Green)'))';
            tmp(tmp < min(Green(:))) = min(Green(:));
            Green = Green ./ tmp;
        end

        if indxNorm == 0
            Green = Green + 1;
        end
        Green = -log(Green);
    end

    if fidY
        waitbar(indP/nChunks, h, 'Yellow channel [Reading file...]')
        Yel = spatialSlabIO('read', fidY, NbPix(1), NbPix(2), Nt, xIdx, 'single');
        Yel = reshape(Yel, [], Nt);

        if b_normalize
            waitbar(indP/nChunks, h, 'Yellow channel [Normalizing data...]')
            Yel = single(filtfilt(lpass_high.sosMatrix, lpass_high.ScaleValues, double(Yel)'))';
            tmp = single(filtfilt(lpass_low.sosMatrix, lpass_low.ScaleValues, double(Yel)'))';
            tmp(tmp < min(Yel(:))) = min(Yel(:));
            Yel = Yel ./ tmp;
        end

        if indxNorm == 0
            Yel = Yel + 1;
        end
        Yel = -log(Yel);
    end
    clear tmp

    waitbar(indP/nChunks, h, 'Computing [HbO] and [HbR]...')

    if fidR * fidG * fidY > 0
        Ainv = pinv(A);
        Hbs = Ainv * ([Red(:), Green(:), Yel(:)]') .* 1e6;
        clear Red Green Yel
    elseif fidR * fidG > 0
        Ainv = pinv(A(1:2,:));
        Hbs = Ainv * ([Red(:), Green(:)]') .* 1e6;
        clear Red Green
    elseif fidG * fidY > 0
        Ainv = pinv(A(2:3,:));
        Hbs = Ainv * ([Green(:), Yel(:)]') .* 1e6;
        clear Green Yel
    else
        Ainv = pinv(A([1 3],:));
        Hbs = Ainv * ([Red(:), Yel(:)]') .* 1e6;
        clear Red Yel
    end

    Hbs = reshape(Hbs, 2, NbPix(1), numel(xIdx), []);
    Hbs = real(Hbs);

    if ~b_RAMsafeMode
        HbO(:,xIdx,:) = squeeze(Hbs(1,:,:,:));
        HbR(:,xIdx,:) = squeeze(Hbs(2,:,:,:));
    else
        waitbar(indP/nChunks, h, 'Writing to files...')
        spatialSlabIO('write', fid_hbo, NbPix(1), NbPix(2), Nt, xIdx, 'single', squeeze(Hbs(1,:,:,:)));
        spatialSlabIO('write', fid_hbr, NbPix(1), NbPix(2), Nt, xIdx, 'single', squeeze(Hbs(2,:,:,:)));
    end
end

close(h);

if b_RAMsafeMode
    fclose(fid_hbr);
    fclose(fid_hbo);

    [~, hboName, hboExt] = fileparts(hboPath);
    [~, hbrName, hbrExt] = fileparts(hbrPath);
    HbO = [hboName hboExt];
    HbR = [hbrName hbrExt];
    return
end

% Save file management (standard mode).
if bSave
    hboPath = fullfile(SaveFolder, default_Output_HbO);
    hbrPath = fullfile(SaveFolder, default_Output_HbR);

    if isfile(hboPath)
        delete(hboPath);
    end
    if isfile(hbrPath)
        delete(hbrPath);
    end

    fidHbO = fopen(hboPath, 'W');
    fwrite(fidHbO, HbO, '*single');
    fclose(fidHbO);

    fidHbR = fopen(hbrPath, 'W');
    fwrite(fidHbR, HbR, '*single');
    fclose(fidHbR);
end

    function info = localPipelineInfo()
        info = PipelineManager.createPipelineInfo( ...
            'HemoCompute', ...
            'Approximate HbO and HbR concentration changes from intrinsic imaging.');

        info = PipelineManager.addInput(info, ...
            'DataFolder', ...
            {'parameter'}, ...
            'Folder containing reflectance channels and AcqInfos.mat.', ...
            'kind', 'parameter', ...
            'position', 1, ...
            'callType', 'positional', ...
            'default', '', ...
            'dataType', 'char');

        info = PipelineManager.addInput(info, ...
            'SaveFolder', ...
            {'parameter'}, ...
            'Folder where HbO/HbR files are saved.', ...
            'kind', 'parameter', ...
            'position', 2, ...
            'callType', 'positional', ...
            'default', '', ...
            'dataType', 'char');

        info = PipelineManager.addInput(info, ...
            'FilterSet', ...
            'parameter', ...
            'Intrinsic imaging filter set.', ...
            'kind', 'parameter', ...
            'position', 3, ...
            'callType', 'positional', ...
            'default', 'none', ...
            'allowed', {'gCaMP','jrGECO','none'}, ...
            'dataType', 'char');

        info = PipelineManager.addInput(info, ...
            'Illumination', ...
            'parameter', ...
            'Illumination wavelength list.', ...
            'kind', 'parameter', ...
            'position', 4, ...
            'callType', 'positional', ...
            'default', {{'red','green'}}, ...
            'dataType', 'cell');

        info = PipelineManager.addInput(info, ...
            'b_normalize', ...
            'parameter', ...
            'Whether to normalize channels before computation.', ...
            'kind', 'parameter', ...
            'position', 5, ...
            'callType', 'positional', ...
            'default', false, ...
            'allowed', {false,true}, ...
            'dataType', 'logical');

        info = PipelineManager.addInput(info, ...
            'HbT_uM', ...
            'parameter', ...
            'Total hemoglobin concentration in micromolar.', ...
            'kind', 'parameter', ...
            'position', 6, ...
            'callType', 'namevalue', ...
            'default', 100, ...
            'dataType', 'numeric');

        info = PipelineManager.addInput(info, ...
            'O2_sat', ...
            'parameter', ...
            'Oxygen saturation percentage.', ...
            'kind', 'parameter', ...
            'position', 7, ...
            'callType', 'namevalue', ...
            'default', 60, ...
            'dataType', 'numeric');

        info = PipelineManager.addInput(info, ...
            'RAMSafeMode', ...
            'parameter', ...
            'Whether to compute Hb outputs in file-backed mode.', ...
            'kind', 'parameter', ...
            'position', 8, ...
            'callType', 'namevalue', ...
            'default', false, ...
            'allowed', {false,true}, ...
            'dataType', 'logical');

        info = PipelineManager.addOutput(info, ...
            'HbO', ...
            {'ImageTimeSeries'}, ...
            'data', ...
            'Oxygenated hemoglobin signal.', ...
            default_Output_HbO, ...
            1, ...
            'isData', true);

        info = PipelineManager.addOutput(info, ...
            'HbR', ...
            {'ImageTimeSeries'}, ...
            'data', ...
            'Deoxygenated hemoglobin signal.', ...
            default_Output_HbR, ...
            2, ...
            'isData', true);
    end
end
