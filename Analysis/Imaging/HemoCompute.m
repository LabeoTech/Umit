function varargout = HemoCompute(SaveFolder, data, varargin)
%HEMOCOMPUTE Approximate HbO and HbR concentration changes from intrinsic imaging.
%
%   outFiles = HemoCompute(SaveFolder, data)
%   [HbOFile, HbRFile] = HemoCompute(SaveFolder, data)
%   outFiles = HemoCompute(SaveFolder, data, 'FilterSet', 'none', ...)
%   HemoCompute('preflight', SaveFolder, parameters)
%
%   This function approximates concentration changes of oxygenated (HbO)
%   and deoxygenated (HbR) hemoglobin from two or three intrinsic imaging
%   wavelengths, and writes them as HbO.dat and HbR.dat in SaveFolder. The
%   channels are streamed in X slabs (low-RAM). It normalizes stable
%   canonical illumination names and resolves the optics ONLY from the
%   folder rig ledger (UMITRigLedger.json, see UMITRigLedger): the camera
%   and illumination spectra come from the ledger's optics snapshot, never
%   from the rig library or the Active rig. The filter set comes from the
%   library, or, with the warning Umitoolbox:HemoCompute:FilterSetFromLedger,
%   from the folder's committed snapshot when the library no longer has it.
%
%   Rig guards (errors, never warnings; preflight runs the same guards):
%       - no ledger, a "no rig" ledger, a Conflicted or an invalid ledger
%         is an error whose message explains how to bind a rig, resolve the
%         conflict or repair the folder;
%       - a filter-set ID committed in this folder with different content
%         (modified in the library) is an error.
%   After a successful computation the optics and filter-set snapshots are
%   committed in the ledger (only when not already committed).
%
%   If the selected channels have different temporal sampling rates or
%   lengths, higher-frequency channels are anti-aliased and resampled to
%   match the lowest-frequency selected channel. Same-frequency channels
%   with small length differences are resampled without anti-aliasing.
%   Lower-frequency channels are not upsampled.
%
%   Inputs:
%       SaveFolder - Folder containing AcqInfos.mat and intrinsic channels.
%                    Receives HbO.dat and HbR.dat.
%       data       - Name of one of the Y-X-T illumination channel .dat files
%                    in SaveFolder (for example 'red.dat'). Arrays, other
%                    extensions, other layouts, and files that are not a
%                    resolved illumination channel are not supported, and
%                    every resolved channel file must be Y-X-T.
%
%                    WHY "data" EXISTS: it is a pipeline dependency, not
%                    the computation's input. HemoCompute does not read
%                    this file in particular: it reads ALL the illumination
%                    channel files resolved from AcqInfos.mat and the Rig
%                    optical configuration in SaveFolder, whichever one is
%                    passed. The input is kept because PipelineManager
%                    schedules a step from its DATA inputs. Without one, the
%                    step would be classified as a source (no DATA input)
%                    and run in the first phase, before setup steps such as
%                    applyRegistrationTformOnFolder and without being
%                    ordered after the step that produces the channels.
%                    Wiring one channel file into "data" makes the step wait
%                    for that producer and run in the main schedule. For
%                    that reason "data" is validated (it must exist, be
%                    Y-X-T, and be one of the resolved channels) but its
%                    values are never used.
%
%   Name-Value parameters:
%       'FilterSet'              - Filter set name.
%       'b_normalize'            - Whether to normalize input channels.
%       'Illuminations'          - Two or three canonical names from
%                                  {'red','green','yellow'}. Historical
%                                  alias 'amber' is normalized to 'yellow'.
%       'HbT_concentration_uM'   - Total hemoglobin concentration.
%       'StO2perc'               - Oxygen saturation percentage.
%
%   Outputs:
%       outFiles          - With one requested output, returns a
%                           {HbOFile,HbRFile} manifest of the Y-X-T .dat
%                           files written into SaveFolder.
%       HbOFile, HbRFile  - With two requested outputs, returns the two
%                           backing filenames. Load those files with
%                           loadData when numeric arrays are required.
%
%   Notes:
%       - Normalization (b_normalize) smooths each channel with a 1 Hz
%         low-pass and divides it by its 1/120 Hz baseline. When 1 Hz is not
%         below the Nyquist frequency of the target timeline (frame rate
%         <= 2 Hz), the 1 Hz smoothing is skipped with a warning.
%       - Samples whose (normalized) intensity is not positive and finite
%         have no optical density: they give NaN in HbO and HbR, with a
%         warning. NaN samples stay NaN.

% Default outputs for pipeline management.
default_Output = {'HbO.dat', 'HbR.dat'};

if nargin == 1 && (ischar(SaveFolder) || (isstring(SaveFolder) && isscalar(SaveFolder))) && ...
        strcmpi(strtrim(char(string(SaveFolder))), 'pipelineInfo')
    varargout{1} = localPipelineInfo();
    return
end

if nargin >= 1 && (ischar(SaveFolder) || (isstring(SaveFolder) && isscalar(SaveFolder))) && ...
        strcmpi(strtrim(char(string(SaveFolder))), 'preflight')
    if nargin ~= 3 || numel(varargin) ~= 1
        error('Umitoolbox:HemoCompute:InvalidPreflightCall', ...
            'Use HemoCompute(''preflight'', SaveFolder, parameters).');
    end
    localPreflight(data, varargin{1});
    if nargout > 0
        varargout = cell(1, nargout);
    end
    return
end

p = inputParser;
p.FunctionName = 'HemoCompute';
addRequired(p, 'SaveFolder', @(x) (ischar(x) || (isstring(x) && isscalar(x))) && isfolder(x));
addRequired(p, 'data');
addParameter(p, 'FilterSet', 'none', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'b_normalize', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Illuminations', {'red','green','yellow'}, ...
    @(x) (iscell(x) && (isempty(x) || isvector(x))) || ...
    (isstring(x) && (isempty(x) || isvector(x))));
addParameter(p, 'HbT_concentration_uM', 100, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'StO2perc', 60, @iValidateStO2Percentage);
parse(p, SaveFolder, data, varargin{:});

SaveFolder = char(string(p.Results.SaveFolder));
filterSet = char(string(p.Results.FilterSet));
b_normalize = p.Results.b_normalize;
illuminationInput = cellstr(string(p.Results.Illuminations));
HbT_uM = double(p.Results.HbT_concentration_uM);
O2_sat = double(p.Results.StO2perc);

illumination = cellfun(@UMITRigStore.normalizeIlluminationName, ...
    illuminationInput, 'UniformOutput', false);

assert(numel(illumination) >= 2, ...
    'Umitoolbox:HemoCompute:InvalidInput', ...
    'At least two different illumination wavelengths must be selected.');

% The wired-in file: a Y-X-T .dat in SaveFolder.
dataName = localResolveDataFile(data, SaveFolder);
[opticalInfo, commitPlan] = localResolveOpticalInfo(SaveFolder, illumination, filterSet);

% The illumination channels are read from SaveFolder, so the wired file
% must be one of them and every channel file must be Y-X-T.
localAssertChannelFiles(dataName, opticalInfo, SaveFolder);

iComputeHb(SaveFolder, illumination, b_normalize, HbT_uM, O2_sat, ...
    opticalInfo, default_Output);

writtenOutputs = localConfirmWrittenOutputs(SaveFolder, default_Output);
localCommitConsumedSnapshots(SaveFolder, commitPlan);
if nargout <= 1
    varargout{1} = writtenOutputs;
elseif nargout == 2
    varargout = writtenOutputs;
else
    error('Umitoolbox:HemoCompute:TooManyOutputs', ...
        'HemoCompute supports at most two outputs.');
end

fprintf('Finished HemoCompute.\n');

    function info = localPipelineInfo()
        info = PipelineManager.createPipelineInfo( ...
            'HemoCompute', ...
            'Approximate HbO and HbR concentration changes from intrinsic imaging.');
        info.version = '2.0.0';

        info = PipelineManager.addInput(info, ...
            'SaveFolder', ...
            'SaveFolder', ...
            'Folder containing AcqInfos.mat and intrinsic channels.', ...
            'kind', 'input', ...
            'position', 1, ...
            'callType', 'positional', ...
            'isData', false);

        info = PipelineManager.addInput(info, ...
            'data', ...
            'ImageTimeSeries', ...
            ['One of the Y-X-T illumination channel .dat files (red, green or ' ...
             'yellow) in SaveFolder. Dependency input only: HemoCompute reads ' ...
             'all resolved channels, and this input orders the step after the ' ...
             'one that produces them (a step with no DATA input would be ' ...
             'scheduled as a source, before setup steps).'], ...
            'position', 2, ...
            'callType', 'positional', ...
            'isData', true, ...
            'supportsFile', true, ...
            'dataMode', 'file');

        info = PipelineManager.addInput(info, ...
            'FilterSet', ...
            'parameter', ...
            'Intrinsic imaging filter set.', ...
            'kind', 'parameter', ...
            'position', 3, ...
            'callType', 'namevalue', ...
            'default', 'none', ...
            'allowed', {'gCaMP','jrGECO','none'}, ...
            'dataType', 'char');

        info = PipelineManager.addInput(info, ...
            'b_normalize', ...
            'parameter', ...
            'Whether to normalize channels before computation.', ...
            'kind', 'parameter', ...
            'position', 4, ...
            'callType', 'namevalue', ...
            'default', true, ...
            'allowed', [true false], ...
            'dataType', 'logical');

        info = PipelineManager.addInput(info, ...
            'Illuminations', ...
            'parameter', ...
            ['Canonical semantic illumination names. Stored values stay stable ' ...
            'when Rig display names or physical spectrum IDs change.'], ...
            'kind', 'parameter', ...
            'position', 5, ...
            'callType', 'namevalue', ...
            'default', {'red','green','yellow'}, ...
            'allowed', {'red';'green';'yellow'}, ...
            'dataType', 'char');

        info = PipelineManager.addInput(info, ...
            'HbT_concentration_uM', ...
            'parameter', ...
            'Total hemoglobin concentration in micromolar.', ...
            'kind', 'parameter', ...
            'position', 6, ...
            'callType', 'namevalue', ...
            'default', 100, ...
            'allowed', [eps, Inf], ...
            'dataType', 'numeric');

        info = PipelineManager.addInput(info, ...
            'StO2perc', ...
            'parameter', ...
            'Oxygen saturation percentage in the inclusive range [0, 100].', ...
            'kind', 'parameter', ...
            'position', 7, ...
            'callType', 'namevalue', ...
            'default', 60, ...
            'allowed', [0 100], ...
            'dataType', 'numeric');

        info = PipelineManager.addOutput(info, ...
            'HbO', ...
            'ImageTimeSeries', ...
            'data', ...
            'Oxygenated hemoglobin backing filename.', ...
            default_Output{1}, ...
            1, ...
            'isData', true);

        info = PipelineManager.addOutput(info, ...
            'HbR', ...
            'ImageTimeSeries', ...
            'data', ...
            'Deoxygenated hemoglobin backing filename.', ...
            default_Output{2}, ...
            2, ...
            'isData', true);

        info.supportsPreflight = true;
        % Results depend on the folder's committed optics: step skipping must
        % not reuse outputs computed before an acknowledged rig conflict.
        info.rigDependent = true;
    end
end

function dataName = localResolveDataFile(data, SaveFolder)
%LOCALRESOLVEDATAFILE Validate the wired-in file: a Y-X-T .dat in SaveFolder.

assert(ischar(data) || (isstring(data) && isscalar(data)), ...
    'Umitoolbox:HemoCompute:UnsupportedInputType', ...
    'Input "data" must be the name of a Y-X-T channel .dat file (arrays are not supported).');

[~, stem, ext] = fileparts(char(string(data)));
assert(isempty(ext) || strcmpi(ext, '.dat'), ...
    'Umitoolbox:HemoCompute:UnsupportedInputFile', ...
    'Unsupported input file extension "%s". Only .dat files are supported.', ext);
dataName = [stem '.dat'];

dataPath = fullfile(SaveFolder, dataName);
assert(isfile(dataPath), ...
    'Umitoolbox:HemoCompute:FileNotFound', ...
    'Input .dat file not found in SaveFolder: "%s".', dataName);
assertDatLayout(loadMetaData(dataPath), {{'Y','X','T'}}, 'HemoCompute');
end

function localAssertChannelFiles(dataName, opticalInfo, SaveFolder)
%LOCALASSERTCHANNELFILES The wired file is a resolved channel; all are Y-X-T.

channelFiles = {};
if isstruct(opticalInfo) && isfield(opticalInfo, 'channels') && ...
        isfield(opticalInfo.channels, 'datFile')
    channelFiles = cellstr(string({opticalInfo.channels.datFile}));
end

assert(any(strcmpi(channelFiles, dataName)), ...
    'Umitoolbox:HemoCompute:DataInputNotAChannel', ...
    ['Input file "%s" is not one of the resolved illumination channels (%s). ' ...
     'HemoCompute processes those channels, so the wired file must be one of them.'], ...
    dataName, strjoin(channelFiles, ', '));

for iFile = 1:numel(channelFiles)
    channelPath = fullfile(SaveFolder, channelFiles{iFile});
    if isfile(channelPath)
        assertDatLayout(loadMetaData(channelPath), {{'Y','X','T'}}, 'HemoCompute');
    end
end
end

function outFiles = localConfirmWrittenOutputs(SaveFolder, expectedFiles)
%LOCALCONFIRMWRITTENOUTPUTS Verify the declared outputs reached disk.

outFiles = cellstr(string(expectedFiles));
missing = outFiles(~cellfun(@(f) isfile(fullfile(SaveFolder, f)), outFiles));

if ~isempty(missing)
    error('Umitoolbox:HemoCompute:MissingOutputFile', ...
        'HemoCompute did not write the expected output file(s): %s.', ...
        strjoin(missing, ', '));
end

end

function tf = iValidateStO2Percentage(value)
%IVALIDATESTO2PERCENTAGE Validate physiological oxygen saturation input.

tf = isnumeric(value) && isscalar(value) && isfinite(value) && ...
    value >= 0 && value <= 100;

if ~tf
    error('Umitoolbox:HemoCompute:InvalidStO2perc', ...
        'StO2perc must be a finite scalar oxygen saturation percentage in the range [0, 100].');
end

end

function localPreflight(saveFolder, parameters)
%LOCALPREFLIGHT Validate Rig optics before pipeline execution.

if ~isstruct(parameters)
    error('Umitoolbox:HemoCompute:InvalidPreflightParameters', ...
        'Pipeline preflight parameters must be a struct array.');
end
names = {parameters.name};
illuminationIdx = find(strcmpi(names, 'Illuminations'), 1);
filterSetIdx = find(strcmpi(names, 'FilterSet'), 1);
if isempty(illuminationIdx) || isempty(filterSetIdx)
    error('Umitoolbox:HemoCompute:InvalidPreflightParameters', ...
        'HemoCompute preflight requires Illuminations and FilterSet parameters.');
end
localResolveOpticalInfo( ...
    saveFolder, parameters(illuminationIdx).value, parameters(filterSetIdx).value);
end

function [opticalInfo, commitPlan] = localResolveOpticalInfo( ...
        saveFolder, illuminations, filterSet)
%LOCALRESOLVEOPTICALINFO Resolve optics from the folder rig ledger.
%   Runs the ledger guards first (no ledger, NoRig, Conflicted, invalid and
%   modified filter set are errors), then maps the acquired channels to the
%   ledger's optics snapshot. Never reads the rig library's rigs or the
%   Active rig, and never writes. COMMITPLAN lists what a successful run
%   must commit.

saveFolder = char(string(saveFolder));
operation = 'HemoCompute';
opticsUse = UMITRigLedger.assertCanUseOptics(saveFolder, 'Operation', operation);
filterSetSnapshot = localFilterSetSnapshot(saveFolder, filterSet);
filterSetNeedsCommit = UMITRigLedger.assertCanUseFilterSet( ...
    saveFolder, filterSetSnapshot, 'Operation', operation);

acqFile = fullfile(saveFolder, 'AcqInfos.mat');
if ~isfile(acqFile)
    error('Umitoolbox:HemoCompute:MissingAcqInfos', ...
        'AcqInfos.mat was not found in "%s".', saveFolder);
end
loaded = load(acqFile, 'AcqInfoStream');
if ~isfield(loaded, 'AcqInfoStream') || ...
        ~isstruct(loaded.AcqInfoStream) || ~isscalar(loaded.AcqInfoStream)
    error('Umitoolbox:HemoCompute:InvalidAcqInfos', ...
        'AcqInfos.mat does not contain a scalar AcqInfoStream structure.');
end
acqInfo = loaded.AcqInfoStream;
[acqInfo, ~] = resolveImportedChannelFallback(acqInfo, saveFolder, ...
    'RequiredChannels', illuminations);

opticalInfo = UMITRigStore.opticalConfigurationFromSnapshots( ...
    opticsUse.snapshot, filterSetSnapshot, acqInfo, illuminations, saveFolder);
opticalInfo.rigUUID = opticsUse.sourceRig.uuid;
opticalInfo.rigID = opticsUse.sourceRig.rigID;

commitPlan = struct( ...
    'opticsHash', opticsUse.hash, ...
    'opticsNeedsCommit', opticsUse.needsCommit, ...
    'filterSet', filterSetSnapshot, ...
    'filterSetNeedsCommit', filterSetNeedsCommit);
end

function snapshot = localFilterSetSnapshot(saveFolder, filterSetID)
%LOCALFILTERSETSNAPSHOT Library filter set, or the folder's committed copy.
%   The library definition is used when it exists, so a modified set is
%   detected by the ledger guard. When the library no longer has the set
%   (another machine, removed user set), the snapshot this folder committed
%   to is used instead.

try
    snapshot = UMITRigLedger.buildFilterSetSnapshot(filterSetID);
catch ME
    if ~ismember(ME.identifier, {'Umitoolbox:UMITRigStore:filterSetNotFound', ...
            'Umitoolbox:UMITRigStore:spectrumNotFound'})
        rethrow(ME)
    end
    snapshot = UMITRigLedger.committedFilterSet(saveFolder, filterSetID);
    if isempty(snapshot)
        rethrow(ME)
    end
    warning('Umitoolbox:HemoCompute:FilterSetFromLedger', ...
        ['Filter set "%s" could not be loaded from this machine''s rig library (%s). ' ...
        'HemoCompute uses the copy committed in the rig ledger of "%s". ' ...
        'To use the library again, import the filter set (with its spectra) into ' ...
        'the rig library under the same ID and content.'], ...
        snapshot.id, ME.message, saveFolder);
end
end

function localCommitConsumedSnapshots(saveFolder, commitPlan)
%LOCALCOMMITCONSUMEDSNAPSHOTS Commit the optics and filter set just used.

if commitPlan.opticsNeedsCommit
    UMITRigLedger.commit(saveFolder, 'optics', commitPlan.opticsHash, ...
        'Origin', 'HemoCompute', 'FirstUsedBy', 'HemoCompute');
end
if commitPlan.filterSetNeedsCommit
    UMITRigLedger.commit(saveFolder, 'filterSet', commitPlan.filterSet, ...
        'Origin', 'HemoCompute', 'FirstUsedBy', 'HemoCompute');
end
end

% =========================================================================
% Hb computation (low-RAM, file-backed)
% =========================================================================
function iComputeHb(SaveFolder, Illumination, b_normalize, HbT_uM, O2_sat, ...
        opticalInfo, outNames)
%ICOMPUTEHB Stream the resolved channels into HbO/HbR files in SaveFolder.
%
% ILLUMINATION holds canonical names; OPTICALINFO is the resolved optical
% configuration (channel files and spectra); OUTNAMES = {HbO, HbR} file
% names. Outputs are written through *_writing.dat scratch files and then
% moved onto the declared names.

if ~strcmp(SaveFolder(end), filesep)
    SaveFolder = [SaveFolder filesep];
end

illuminationLower = string(Illumination);
assert(numel(unique(illuminationLower)) >= 2, ...
    'Umitoolbox:HemoCompute:InvalidIllumination', ...
    'At least two different illumination wavelengths are needed for Hb computation.');

% Resolve selected channel files and file-specific metadata.
channelNames = {'red', 'green', 'yellow'};
channelFiles = cell(1, 3);
for iResolved = 1:numel(opticalInfo.channels)
    row = find(strcmp(channelNames, opticalInfo.channels(iResolved).name), 1);
    channelFiles{row} = fullfile(SaveFolder, opticalInfo.channels(iResolved).datFile);
end
selectedChannels = ismember(channelNames, cellstr(illuminationLower));
assert(isequal(find(opticalInfo.activeRows(:))', find(selectedChannels)), ...
    'Umitoolbox:HemoCompute:OpticalInfoChannelMismatch', ...
    ['OpticalInfo marks illumination rows [%s] as active, but the ' ...
     'selected channels are [%s]. The extinction matrix would contain ' ...
     'NaN rows for the selected channels.'], ...
    num2str(find(opticalInfo.activeRows(:))'), ...
    num2str(find(selectedChannels)));

channelInfo = cell(1, 3);
slabList = cell(1, 3);
hasChannel = false(1, 3);
cleanupList = cell(1, 3); % holds the slab onCleanup objects

for ii = find(selectedChannels)
    assert(isfile(channelFiles{ii}), ...
        'Umitoolbox:HemoCompute:MissingChannelFile', ...
        'Required channel file was not found: %s', channelFiles{ii});

    channelInfo{ii} = loadMetaData(channelFiles{ii});

    assert(isfield(channelInfo{ii}, 'dimNames') && isfield(channelInfo{ii}, 'dimSizes') && ...
            datAxisSize(channelInfo{ii}, 'T') > 0 && isfield(channelInfo{ii}, 'frameRateHz'), ...
        'Umitoolbox:HemoCompute:InvalidChannelMetadata', ...
        'Metadata for "%s" must describe the Y, X, and T sizes and the frame rate.', ...
        channelFiles{ii});

    slabList{ii} = spatialSlabIO('open', channelFiles{ii}, 'Info', channelInfo{ii});
    hasChannel(ii) = true;
    thisSlab = slabList{ii};
    cleanupList{ii} = onCleanup(@() spatialSlabIO('close', thisSlab));
end

idxSelected = find(selectedChannels);
nativeFreq = nan(1,3);
nativeLength = nan(1,3);
for ii = idxSelected
    nativeFreq(ii) = double(channelInfo{ii}.frameRateHz);
    nativeLength(ii) = double(datAxisSize(channelInfo{ii}, 'T'));
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

Ny = double(datAxisSize(channelInfo{idxTarget}, 'Y'));
Nx = double(datAxisSize(channelInfo{idxTarget}, 'X'));
Nt = double(datAxisSize(channelInfo{idxTarget}, 'T'));
Freq = double(channelInfo{idxTarget}.frameRateHz);

% Header exposure of HbO/HbR: the selected channels' common exposure, or
% NaN when they differ (the outputs combine the channels).
selectedExposures = cellfun(@(c) double(c.exposureMsec), channelInfo(idxSelected));
if all(arrayfun(@(e) isequaln(e, selectedExposures(1)), selectedExposures))
    outExposureMsec = selectedExposures(1);
else
    outExposureMsec = NaN;
end

for ii = idxSelected
    assert(isequal(double(datAxisSize(channelInfo{ii}, 'Y')), Ny) && ...
           isequal(double(datAxisSize(channelInfo{ii}, 'X')), Nx), ...
        'Umitoolbox:HemoCompute:SpatialMismatch', ...
        'All selected channels must have matching spatial dimensions.');
end

NbPix = [Ny, Nx];

% Check whether the selected data are already normalized.
indxNorm = [-2 -2 -2];

fprintf('Checking channel data...\n');
for ii = idxSelected
    thisInfo = channelInfo{ii};
    thisNt = double(datAxisSize(thisInfo, 'T'));

    % Sample 10 frames from the native channel timeline.
    frIdx = unique(floor(linspace(1, thisNt, 10)));
    frIdx(frIdx < 1) = [];
    tmp = zeros(double(datAxisSize(thisInfo, 'Y')), double(datAxisSize(thisInfo, 'X')), 'single');

    for jj = 1:numel(frIdx)
        frame = single(spatialSlabIO('read', slabList{ii}, 1:datAxisSize(thisInfo, 'X'), frIdx(jj)));
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
A = hemooptics.ioi_epsilon_pathlength( ...
    'Hillman', HbT_uM, O2_sat, 100 - O2_sat, opticalInfo);

% Temporal filters used for input normalization after all channels are on
% the target timeline: a 1 Hz smoothing (lpass_high) and a 1/120 Hz
% baseline (lpass_low). They are only needed when normalizing. The 1 Hz
% smoothing cannot be designed at or above the Nyquist frequency; it is
% then skipped and only the baseline division is applied.
lpass_high = [];
lpass_low = [];
if b_normalize
    if 1 < Freq / 2
        f = fdesign.lowpass('N,F3dB', 4, 1, Freq);
        lpass_high = design(f, 'butter');
    else
        warning('Umitoolbox:HemoCompute:NormalizationSmoothingSkipped', ...
            ['The 1 Hz smoothing of the normalization is skipped: the ' ...
             'channels are at %g Hz (Nyquist %g Hz), so 1 Hz is not below ' ...
             'the Nyquist frequency. Each channel is only divided by its ' ...
             '1/120 Hz baseline.'], Freq, Freq / 2);
    end

    f = fdesign.lowpass('N,F3dB', 4, 1/120, Freq);
    lpass_low = design(f, 'butter');
end

% Chunking strategy for continuous YXT data. Use the largest native length
% so memory is not underestimated when repeated illuminations are present.
maxNativeLength = max(nativeLength(idxSelected));
nChunks = calculateMaxChunkSize(Ny * Nx * maxNativeLength * 4, 12, .1);
chunkX  = ceil(NbPix(2) / nChunks);
nChunks = ceil(NbPix(2) / chunkX);

% Write through scratch files, then move them onto the declared outputs.
% Renaming the output when it already exists would make every pipeline
% re-run write to a different file and leave the stale original in place.
hboOutPath = fullfile(SaveFolder, outNames{1});
hbrOutPath = fullfile(SaveFolder, outNames{2});
[~, hboBaseName, hboExt] = fileparts(outNames{1});
[~, hbrBaseName, hbrExt] = fileparts(outNames{2});
hboPath = fullfile(SaveFolder, [hboBaseName '_writing' hboExt]);
hbrPath = fullfile(SaveFolder, [hbrBaseName '_writing' hbrExt]);

% Outputs are on the target timeline (lowest-frequency channel).
outHeaderArgs = {'dataClass', 'single', 'dimNames', {'Y', 'X', 'T'}, ...
    'dimSizes', [Ny, Nx, Nt], 'frameRateHz', Freq, 'exposureMsec', outExposureMsec};
slab_hbo = spatialSlabIO('create', hboPath, ...
    datHeaderFromInfo(channelInfo{idxTarget}, hboBaseName, outHeaderArgs{:}));
c_hbo = onCleanup(@() spatialSlabIO('close', slab_hbo));

slab_hbr = spatialSlabIO('create', hbrPath, ...
    datHeaderFromInfo(channelInfo{idxTarget}, hbrBaseName, outHeaderArgs{:}));
c_hbr = onCleanup(@() spatialSlabIO('close', slab_hbr));

% Clamp-floor pre-pass. min(Red/Green/Yel(:)) must be a whole-image
% reduction, not a chunk-local one, or the clamp floor (and therefore the
% normalized output) depends on chunkX, which itself varies with available
% RAM via calculateMaxChunkSize. Fold the minimum across chunks one slab at
% a time so this stays RAM-safe.
globalMin = Inf(1, 3);
if b_normalize
    for indP = 1:nChunks
        xStart = (indP - 1) * chunkX + 1;
        xEnd   = min(xStart + chunkX - 1, NbPix(2));
        xIdx   = xStart:xEnd;

        for ii = idxSelected
            chan = iReadChannel(slabList{ii}, channelInfo{ii}, xIdx, Nt, Freq);
            chan = iSmooth(chan, lpass_high);
            globalMin(ii) = min(globalMin(ii), min(chan(:)));
        end
    end
    clear chan
end

% Computation loop.
nInvalid = 0;
nSamples = 0;
channelLabels = {'Red', 'Green', 'Yellow'};
h = waitbar(0, 'Computing');
for indP = 1:nChunks
    xStart = (indP - 1) * chunkX + 1;
    xEnd   = min(xStart + chunkX - 1, NbPix(2));
    xIdx   = xStart:xEnd;

    h.Name = ['HemoCompute (chunk ' num2str(indP) '/' num2str(nChunks) ')'];
    drawnow()

    OD = cell(1, 3);
    for ii = idxSelected
        waitbar(indP/nChunks, h, [channelLabels{ii} ' channel [Reading file...]'])
        chan = iReadChannel(slabList{ii}, channelInfo{ii}, xIdx, Nt, Freq);

        if b_normalize
            waitbar(indP/nChunks, h, [channelLabels{ii} ' channel [Normalizing data...]'])
            chan = iSmooth(chan, lpass_high);
            tmp = single(filtfilt(lpass_low.sosMatrix, lpass_low.ScaleValues, double(chan)'))';
            tmp(tmp < globalMin(ii)) = globalMin(ii);
            chan = chan ./ tmp;
        end

        if indxNorm == 0
            chan = chan + 1;
        end

        % A sample that is not positive and finite has no optical density:
        % it becomes NaN in HbO and HbR. NaN samples stay NaN.
        isInvalid = chan <= 0 | isinf(chan);
        nInvalid = nInvalid + nnz(isInvalid);
        nSamples = nSamples + numel(chan);
        chan(isInvalid) = NaN;

        OD{ii} = -log(chan);
    end
    clear tmp chan

    waitbar(indP/nChunks, h, 'Computing [HbO] and [HbR]...')

    if hasChannel(1) * hasChannel(2) * hasChannel(3) > 0
        Ainv = pinv(A);
        Hbs = Ainv * ([OD{1}(:), OD{2}(:), OD{3}(:)]') .* 1e6;
    elseif hasChannel(1) * hasChannel(2) > 0
        Ainv = pinv(A(1:2,:));
        Hbs = Ainv * ([OD{1}(:), OD{2}(:)]') .* 1e6;
    elseif hasChannel(2) * hasChannel(3) > 0
        Ainv = pinv(A(2:3,:));
        Hbs = Ainv * ([OD{2}(:), OD{3}(:)]') .* 1e6;
    else
        Ainv = pinv(A([1 3],:));
        Hbs = Ainv * ([OD{1}(:), OD{3}(:)]') .* 1e6;
    end
    clear OD

    Hbs = reshape(Hbs, 2, NbPix(1), numel(xIdx), []);

    waitbar(indP/nChunks, h, 'Writing to files...')
    spatialSlabIO('write', slab_hbo, xIdx, reshape(Hbs(1,:,:,:), NbPix(1), numel(xIdx), []));
    spatialSlabIO('write', slab_hbr, xIdx, reshape(Hbs(2,:,:,:), NbPix(1), numel(xIdx), []));
end

close(h);

if nInvalid > 0
    warning('Umitoolbox:HemoCompute:NonPositiveIntensity', ...
        ['%d of %d channel samples (%.3g%%) were not positive and finite ' ...
         'after normalization, so they have no optical density: HbO and ' ...
         'HbR are NaN at those samples. Masked pixels stored as 0 are a ' ...
         'common cause.'], nInvalid, nSamples, 100 * nInvalid / nSamples);
end

spatialSlabIO('finalize', slab_hbr);
spatialSlabIO('finalize', slab_hbo);
clear c_hbo c_hbr

[moveOk, moveMsg] = movefile(hboPath, hboOutPath, 'f');
assert(moveOk, 'Umitoolbox:HemoCompute:OutputMoveFailed', ...
    'Failed to move "%s" onto "%s": %s', hboPath, hboOutPath, moveMsg);
[moveOk, moveMsg] = movefile(hbrPath, hbrOutPath, 'f');
assert(moveOk, 'Umitoolbox:HemoCompute:OutputMoveFailed', ...
    'Failed to move "%s" onto "%s": %s', hbrPath, hbrOutPath, moveMsg);
end

function chan = iReadChannel(slab, info, xIdx, Nt, Freq)
%IREADCHANNEL Read X columns XIDX of a channel on the target timeline (pixels x Nt).
chan = spatialSlabIO('read', slab, xIdx);
chan = iResampleChannelToLowestFrequency(chan, info, Nt, Freq);
chan = reshape(chan, [], Nt);
end

function chan = iSmooth(chan, lpass_high)
%ISMOOTH Apply the 1 Hz normalization smoothing (pixels x T); none when empty.
if ~isempty(lpass_high)
    chan = single(filtfilt(lpass_high.sosMatrix, lpass_high.ScaleValues, double(chan)'))';
end
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

sourceNt = double(datAxisSize(sourceMeta, 'T'));
sourceFreq = double(sourceMeta.frameRateHz);

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
