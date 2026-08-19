function [HbO, HbR] = HemoCompute(DataFolder, SaveFolder, FilterSet, Illumination, b_normalize, varargin)
%HEMOCOMPUTE Approximate HbO and HbR concentration changes from intrinsic imaging.
%
%   [HbO, HbR] = HemoCompute(DataFolder, SaveFolder, FilterSet, Illumination, b_normalize)
%   [HbO, HbR] = HemoCompute(..., 'HbT_uM', 100, 'O2_sat', 60, 'RAMSafeMode', false)
%
%   This function approximates concentration changes of oxygenated (HbO)
%   and deoxygenated (HbR) hemoglobin from two or three intrinsic imaging
%   wavelengths.
%   When UMITRigStore is available but AcqInfos does not identify a Rig, the
%   existing active default Rig is used with a warning. AcqInfos is not
%   modified by this fallback.
%
%   If selected channels have different temporal sampling rates or lengths,
%   higher-frequency channels are anti-aliased and resampled to match the
%   lowest-frequency selected channel. Same-frequency channels with small
%   length differences are resampled without anti-aliasing. Lower-frequency
%   channels are not upsampled.
%
%   Inputs:
%       DataFolder    - Folder containing red.dat, green.dat, and/or yellow.dat.
%       SaveFolder    - Folder where HbO/HbR outputs are saved. If empty,
%                       outputs are returned only in RAM in standard mode.
%       FilterSet     - UMITRigStore filter-set ID when the Rig service is
%                       available, or a legacy FilterSets.mat name otherwise.
%       Illumination  - Cell array containing at least two wavelengths from
%                       {'red','green','yellow','amber'}.
%       b_normalize   - Logical scalar. If true, normalize input channels
%                       when needed.
%
%   Name-Value parameters:
%       'HbT_uM'      - Total hemoglobin concentration in micromolar.
%                       Default = 100.
%       'O2_sat'      - Oxygen saturation percentage. Default = 60.
%       'OpticalInfo' - Optional resolved optical-information structure.
%                       It must contain the spectral fields documented by
%                       ioi_epsilon_pathlength plus a channels struct array
%                       with name, datFile, and camIdx for each requested
%                       illumination. When supplied, no optical library
%                       lookup is made.
%                       Otherwise UMITRigStore is used when it is on the
%                       MATLAB path; if it is absent, the legacy MAT optical
%                       definitions under IOIAnalysis are used.
%       'RAMSafeMode' - Logical scalar. If true, write HbO/HbR directly to
%                       disk instead of holding them in RAM. Default = false.
%
%   Outputs:
%       HbO           - In standard mode, numeric HbO array. In RAM-safe
%                       mode, output filename.
%       HbR           - In standard mode, numeric HbR array. In RAM-safe
%                       mode, output filename.

% Default output names for pipeline management.
default_Output_HbO = 'HbO.dat';
default_Output_HbR = 'HbR.dat';

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
addParameter(p, 'OpticalInfo', [], @(x) isempty(x) || (isstruct(x) && isscalar(x)));
addParameter(p, 'RAMSafeMode', false, @(x) islogical(x) && isscalar(x));
parse(p, DataFolder, SaveFolder, FilterSet, Illumination, b_normalize, varargin{:});

DataFolder = char(string(p.Results.DataFolder));
SaveFolder = char(string(p.Results.SaveFolder));
FilterSet = char(string(p.Results.FilterSet));
Illumination = p.Results.Illumination;
b_normalize = p.Results.b_normalize;
HbT_uM = double(p.Results.HbT_uM);
O2_sat = double(p.Results.O2_sat);
opticalInfo = p.Results.OpticalInfo;
b_RAMsafeMode = p.Results.RAMSafeMode;

if ~strcmp(DataFolder(end), filesep)
    DataFolder = [DataFolder filesep];
end

bSave = ~isempty(SaveFolder);
if bSave
    assert(isfolder(SaveFolder), ...
        'Umitoolbox:HemoCompute:InvalidSaveFolder', ...
        'SaveFolder must be an existing folder or empty.');
    if ~strcmp(SaveFolder(end), filesep)
        SaveFolder = [SaveFolder filesep];
    end
end

illuminationLower = string(cellfun(@localNormalizeIlluminationName, ...
    Illumination, 'UniformOutput', false));
assert(numel(unique(illuminationLower)) >= 2, ...
    'Umitoolbox:HemoCompute:InvalidIllumination', ...
    'At least two different illumination wavelengths are needed for Hb computation.');

useRigOptics = ~isempty(opticalInfo);
if ~useRigOptics
    hasRigStore = exist('UMITRigStore', 'class') == 8 || ...
        exist('UMITRigStore', 'file') == 2;
    if hasRigStore
        opticalInfo = localResolveRigOpticalInfo( ...
            DataFolder, cellstr(illuminationLower), FilterSet);
        useRigOptics = true;
    end
end

% Resolve selected channel files and file-specific metadata.
channelNames = {'red', 'green', 'yellow'};
if useRigOptics
    channelFiles = cell(1, 3);
    for iResolved = 1:numel(opticalInfo.channels)
        row = find(strcmp(channelNames, opticalInfo.channels(iResolved).name), 1);
        channelFiles{row} = fullfile(DataFolder, opticalInfo.channels(iResolved).datFile);
    end
else
    channelFiles = { ...
        fullfile(DataFolder, 'red.dat'), ...
        fullfile(DataFolder, 'green.dat'), ...
        fullfile(DataFolder, 'yellow.dat')};
    acq = localLoadAcqInfo(DataFolder);
    assert(isfield(acq, 'Camera_Model') && ~isempty(acq.Camera_Model), ...
        'Umitoolbox:HemoCompute:InvalidAcqInfos', ...
        'AcqInfoStream must contain Camera_Model for legacy optical lookup.');
    assert(any(strcmpi(FilterSet, {'gcamp','jrgeco','none'})), ...
        'Umitoolbox:HemoCompute:InvalidFilterSet', ...
        'FilterSet must be one of: gCaMP, jrGECO, none.');
end
selectedChannels = ismember(channelNames, cellstr(illuminationLower));

channelInfo = cell(1, 3);
fidList = zeros(1, 3);
cleanupList = cell(1, 3);

for ii = find(selectedChannels)
    assert(isfile(channelFiles{ii}), ...
        'Umitoolbox:HemoCompute:MissingChannelFile', ...
        'Required channel file was not found: %s', channelFiles{ii});

    channelInfo{ii} = loadMetaData(channelFiles{ii});

    assert(isfield(channelInfo{ii}, 'Height') && isfield(channelInfo{ii}, 'Width') && ...
            isfield(channelInfo{ii}, 'Length') && isfield(channelInfo{ii}, 'FrameRateHz'), ...
        'Umitoolbox:HemoCompute:InvalidChannelMetadata', ...
        'Metadata for "%s" must contain Height, Width, Length, and FrameRateHz.', ...
        channelFiles{ii});

    fidList(ii) = fopen(channelFiles{ii}, 'r');
    assert(fidList(ii) ~= -1, ...
        'Umitoolbox:HemoCompute:FileOpenError', ...
        'Could not open channel file "%s".', channelFiles{ii});
    thisFid = fidList(ii);
    cleanupList{ii} = onCleanup(@() safeFclose(thisFid)); 
end

idxSelected = find(selectedChannels);
nativeFreq = nan(1,3);
nativeLength = nan(1,3);
for ii = idxSelected
    nativeFreq(ii) = double(channelInfo{ii}.FrameRateHz);
    nativeLength(ii) = double(channelInfo{ii}.Length);
end

% Assert that all selected channels span the same recording duration before
% any temporal resampling. Interpolation is only valid for channels covering
% the same acquisition interval; it should not silently repair cropped or
% truncated files.
durationTolSec = 1e-3;
nativeDurationSec = nativeLength ./ nativeFreq;
targetDurationForCheck = nativeDurationSec(idxSelected(1));
for ii = idxSelected
    assert(isfinite(nativeDurationSec(ii)) && nativeDurationSec(ii) > 0, ...
        'Umitoolbox:HemoCompute:InvalidChannelDuration', ...
        'Channel "%s" has invalid duration metadata: Length=%g, FrameRateHz=%g.', ...
        channelNames{ii}, nativeLength(ii), nativeFreq(ii));

    assert(abs(nativeDurationSec(ii) - targetDurationForCheck) <= durationTolSec, ...
        'Umitoolbox:HemoCompute:DurationMismatch', ...
        ['Selected channels do not span the same recording duration. ' ...
         'Channel "%s": Length=%g, FrameRateHz=%g, Duration=%0.6f s. ' ...
         'Reference duration=%0.6f s.'], ...
        channelNames{ii}, nativeLength(ii), nativeFreq(ii), ...
        nativeDurationSec(ii), targetDurationForCheck);
end

% Use the lowest-frequency selected channel as the target timeline.
targetFreq = min(nativeFreq(idxSelected));
freqTol = max(1e-6, abs(targetFreq) * 1e-6);
idxTargetCandidates = idxSelected(abs(nativeFreq(idxSelected) - targetFreq) <= freqTol);
[~, idxMinLength] = min(nativeLength(idxTargetCandidates));
idxTarget = idxTargetCandidates(idxMinLength);

Ny = double(channelInfo{idxTarget}.Height);
Nx = double(channelInfo{idxTarget}.Width);
Nt = double(channelInfo{idxTarget}.Length);
Freq = double(channelInfo{idxTarget}.FrameRateHz);
Datatype = char(string(channelInfo{idxTarget}.Datatype)); %#ok<NASGU>

for ii = idxSelected
    assert(isequal(double(channelInfo{ii}.Height), Ny) && ...
           isequal(double(channelInfo{ii}.Width), Nx), ...
        'Umitoolbox:HemoCompute:SpatialMismatch', ...
        'All selected channels must have matching spatial dimensions.');
end

% Data size is strictly continuous YXT.
datsz = [Ny, Nx, Nt];
NbPix = [Ny, Nx];

% Check whether the selected data are already normalized.
indxNorm = [-2 -2 -2];

fprintf('Checking channel data...\n');
for ii = idxSelected
    thisFid = fidList(ii);
    thisInfo = channelInfo{ii};
    thisNt = double(thisInfo.Length);
    bytesPerFrame = double(thisInfo.Height) * double(thisInfo.Width) * getByteSize(thisInfo.Datatype);

    % Sample 10 frames from the native channel timeline.
    frIdx = unique(floor(linspace(1, thisNt, 10)));
    frIdx(frIdx < 1) = [];
    tmp = zeros(double(thisInfo.Height), double(thisInfo.Width), 'single');

    for jj = 1:numel(frIdx)
        fseek(thisFid, (frIdx(jj) - 1) * bytesPerFrame, 'bof');
        frame = fread(thisFid, [double(thisInfo.Height), double(thisInfo.Width)], '*single');
        tmp = tmp + frame;
    end

    tmp = tmp ./ numel(frIdx);
    Mdat = mean(tmp, 'all', 'omitnan');

    if Mdat > .75 && Mdat < 1.25
        indxNorm(ii) = 1;
    elseif Mdat > -.25 && Mdat < .25
        indxNorm(ii) = 0;
    else
        indxNorm(ii) = -1;
    end
end

selectedNorm = indxNorm(idxSelected);
if isempty(selectedNorm) || any(selectedNorm ~= selectedNorm(1))
    error('Umitoolbox:HemoCompute:HeterogeneousInput', ...
        'The input data is heterogeneous. All selected channels must be preprocessed in the same way.');
end

indxNorm = selectedNorm(1);

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
if useRigOptics
    A = ioi_epsilon_pathlength( ...
        'Hillman', HbT_uM, O2_sat, 100 - O2_sat, opticalInfo);
else
    A = ioi_epsilon_pathlength( ...
        'Hillman', HbT_uM, O2_sat, 100 - O2_sat, FilterSet, acq.Camera_Model);
end

% Temporal filters used for input normalization after all channels are on
% the target timeline.
f = fdesign.lowpass('N,F3dB', 4, 1, Freq);
lpass_high = design(f, 'butter');

f = fdesign.lowpass('N,F3dB', 4, 1/120, Freq);
lpass_low = design(f, 'butter');

% Chunking strategy for continuous YXT data. Use the largest native length
% so memory is not underestimated when repeated illuminations are present.
maxNativeLength = max(nativeLength(idxSelected));
nChunks = calculateMaxChunkSize(Ny * Nx * maxNativeLength * 4, 12, .1);
chunkX  = ceil(NbPix(2) / nChunks);

% Output allocation / preallocation.
if ~b_RAMsafeMode
    HbO = zeros(datsz, 'single');
    HbR = zeros(datsz, 'single');
else
    assert(bSave, ...
        'Umitoolbox:HemoCompute:MissingSaveFolder', ...
        'SaveFolder must be provided when RAMSafeMode is true.');

    % Write through scratch files, then move them onto the declared
    % outputs. Renaming the output when it already exists would make every
    % pipeline re-run write to a different file and leave the stale
    % original in place.
    hboOutPath = fullfile(SaveFolder, default_Output_HbO);
    hbrOutPath = fullfile(SaveFolder, default_Output_HbR);
    [~, hboBaseName, hboExt] = fileparts(default_Output_HbO);
    [~, hbrBaseName, hbrExt] = fileparts(default_Output_HbR);
    hboPath = fullfile(SaveFolder, [hboBaseName '_writing' hboExt]);
    hbrPath = fullfile(SaveFolder, [hbrBaseName '_writing' hbrExt]);

    preallocateDatFile(hboPath, [Ny, Nx, Nt], 'single');
    fid_hbo = fopen(hboPath, 'r+');
    c_hbo = onCleanup(@() safeFclose(fid_hbo)); 

    preallocateDatFile(hbrPath, [Ny, Nx, Nt], 'single');
    fid_hbr = fopen(hbrPath, 'r+');
    c_hbr = onCleanup(@() safeFclose(fid_hbr)); 
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

    if fidList(1)
        waitbar(indP/nChunks, h, 'Red channel [Reading file...]')
        Red = spatialSlabIO('read', fidList(1), NbPix(1), NbPix(2), ...
            channelInfo{1}.Length, xIdx, channelInfo{1}.Datatype);
        Red = iResampleChannelToLowestFrequency(Red, channelInfo{1}, Nt, Freq);
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

    if fidList(2)
        waitbar(indP/nChunks, h, 'Green channel [Reading file...]')
        Green = spatialSlabIO('read', fidList(2), NbPix(1), NbPix(2), ...
            channelInfo{2}.Length, xIdx, channelInfo{2}.Datatype);
        Green = iResampleChannelToLowestFrequency(Green, channelInfo{2}, Nt, Freq);
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

    if fidList(3)
        waitbar(indP/nChunks, h, 'Yellow channel [Reading file...]')
        Yel = spatialSlabIO('read', fidList(3), NbPix(1), NbPix(2), ...
            channelInfo{3}.Length, xIdx, channelInfo{3}.Datatype);
        Yel = iResampleChannelToLowestFrequency(Yel, channelInfo{3}, Nt, Freq);
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

    if fidList(1) * fidList(2) * fidList(3) > 0
        Ainv = pinv(A);
        Hbs = Ainv * ([Red(:), Green(:), Yel(:)]') .* 1e6;
        clear Red Green Yel
    elseif fidList(1) * fidList(2) > 0
        Ainv = pinv(A(1:2,:));
        Hbs = Ainv * ([Red(:), Green(:)]') .* 1e6;
        clear Red Green
    elseif fidList(2) * fidList(3) > 0
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

    [moveOk, moveMsg] = movefile(hboPath, hboOutPath, 'f');
    assert(moveOk, 'Umitoolbox:HemoCompute:OutputMoveFailed', ...
        'Failed to move "%s" onto "%s": %s', hboPath, hboOutPath, moveMsg);
    [moveOk, moveMsg] = movefile(hbrPath, hbrOutPath, 'f');
    assert(moveOk, 'Umitoolbox:HemoCompute:OutputMoveFailed', ...
        'Failed to move "%s" onto "%s": %s', hbrPath, hbrOutPath, moveMsg);

    HbO = default_Output_HbO;
    HbR = default_Output_HbR;
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

    saveData(hboPath, single(HbO));
    saveData(hbrPath, single(HbR));
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

function canonicalName = localNormalizeIlluminationName(name)
%LOCALNORMALIZEILLUMINATIONNAME Normalize legacy illumination names locally.

if ~(ischar(name) || (isstring(name) && isscalar(name)))
    error('Umitoolbox:HemoCompute:InvalidIllumination', ...
        'Illumination entries must be scalar text.');
end
canonicalName = lower(strtrim(char(string(name))));
if strcmp(canonicalName, 'amber')
    canonicalName = 'yellow';
end
if ~ismember(canonicalName, {'red', 'green', 'yellow'})
    error('Umitoolbox:HemoCompute:InvalidIllumination', ...
        'Illumination entries must be red, green, yellow, or amber.');
end
end

function acqInfo = localLoadAcqInfo(dataFolder)
%LOCALLOADACQINFO Load and validate the folder-level acquisition structure.

acqFile = fullfile(char(string(dataFolder)), 'AcqInfos.mat');
if ~isfile(acqFile)
    error('Umitoolbox:HemoCompute:MissingAcqInfos', ...
        'AcqInfos.mat was not found in "%s".', char(string(dataFolder)));
end
loaded = load(acqFile, 'AcqInfoStream');
if ~isfield(loaded, 'AcqInfoStream') || ...
        ~isstruct(loaded.AcqInfoStream) || ~isscalar(loaded.AcqInfoStream)
    error('Umitoolbox:HemoCompute:InvalidAcqInfos', ...
        'AcqInfos.mat does not contain a scalar AcqInfoStream structure.');
end
acqInfo = loaded.AcqInfoStream;
end

function opticalInfo = localResolveRigOpticalInfo(dataFolder, illuminations, filterSet)
%LOCALRESOLVERIGOPTICALINFO Resolve optics through an available Rig store.

acqInfo = localLoadAcqInfo(dataFolder);
[acqInfo, ~] = resolveImportedChannelFallback(acqInfo, dataFolder, ...
    'RequiredChannels', illuminations);
if isfield(acqInfo, 'rigUUID') && ~isempty(acqInfo.rigUUID)
    rigStore = UMITRigStore.open(acqInfo.rigUUID);
elseif isfield(acqInfo, 'rigID') && ~isempty(acqInfo.rigID)
    rigStore = UMITRigStore.openByRigID(acqInfo.rigID);
else
    rigStore = localResolveActiveRig();
    rigInfo = rigStore.getRigInfo();
    warning('Umitoolbox:HemoCompute:MissingRigUsingActiveRig', ...
        ['AcqInfoStream does not identify a Rig. Using active Rig "%s" ' ...
        'for HemoCompute without modifying AcqInfos.mat.'], rigInfo.rigID);
end
opticalInfo = rigStore.resolveOpticalConfiguration( ...
    acqInfo, illuminations, filterSet, dataFolder);
end

function rigStore = localResolveActiveRig()
%LOCALRESOLVEACTIVERIG Resolve an existing active Rig without creating one.

rigStore = UMITRigStore.getDefaultRig();
if ~isempty(rigStore)
    return
end

rigs = UMITRigStore.listRigs();
candidates = rigs(rigs.Status == "active" & rigs.IsReadable, :);
if height(candidates) == 1
    rigStore = UMITRigStore.open(char(candidates.RigUUID(1)));
    return
end

if isempty(candidates)
    detail = 'No active readable Rig exists.';
else
    detail = 'Multiple active Rigs exist and none is selected as default.';
end
error('Umitoolbox:HemoCompute:MissingRig', ...
    ['AcqInfoStream does not identify a Rig. %s Select an active default ' ...
    'Rig before running HemoCompute.'], detail);
end

function data = iResampleChannelToLowestFrequency(data, sourceMeta, targetNt, targetFreq)
%IRESAMPLECHANNELTOLOWESTFREQUENCY Match one channel to the target timeline.
%
%   data = iResampleChannelToLowestFrequency(data, sourceMeta, targetNt,
%   targetFreq) resamples data along the last dimension so it has targetNt
%   samples. If the source frequency is higher than the target frequency, an
%   anti-aliasing low-pass filter is applied before interpolation. If the
%   source frequency matches the target but the lengths differ, interpolation
%   is applied without anti-aliasing. Lower-frequency sources are not
%   upsampled.

sourceNt = double(sourceMeta.Length);
sourceFreq = double(sourceMeta.FrameRateHz);

if sourceNt == targetNt
    data = single(data);
    return
end

freqTol = max(1e-6, abs(targetFreq) * 1e-6);
inputSize = size(data);
if inputSize(end) ~= sourceNt
    error(['Temporal resampling failed. The input data length does not match ' ...
        'the source metadata Length.']);
end

data = reshape(data, [], sourceNt);

if sourceFreq > targetFreq + freqTol
    aaCutoff = 0.45 * targetFreq;
    if aaCutoff > 0 && aaCutoff < sourceFreq/2
        f = fdesign.lowpass('N,F3dB', 4, aaCutoff, sourceFreq);
        lpass = design(f,'butter');
        data = single(filtfilt(lpass.sosMatrix, lpass.ScaleValues, double(data)'))';
    end
elseif sourceFreq < targetFreq - freqTol
    error(['Cannot upsample a lower-frequency channel for Hb computation. ' ...
        'Select channels with compatible frequencies or downsample all ' ...
        'channels to the lowest-frequency timeline.']);
end

xSource = linspace(0, 1, sourceNt);
xTarget = linspace(0, 1, targetNt);
data = interp1(xSource, single(data)', xTarget, 'linear', 'extrap')';

inputSize(end) = targetNt;
data = reshape(single(data), inputSize);

end
