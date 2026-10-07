function [outData, metaData] = split_data_by_event(data, SaveFolder, varargin)
%SPLIT_DATA_BY_EVENT Split continuous image time series into event trials.
%
%   [outData, metaData] = split_data_by_event(data, SaveFolder)
%   [outData, metaData] = split_data_by_event(data, SaveFolder, 'FrameRateHz', rate)
%
%   This function uses the event definitions stored in "events.mat" to
%   split a continuous image time series into event instances. The output
%   has dimensions Y X T E, where E is the event-instance axis. The output
%   format follows the input: .dat for array and .dat inputs, .umt for UMT
%   inputs (see Output).
%
%   Inputs:
%       data       - One of (the data must have a T axis):
%                    * Numeric image time series with dimensions Y X T
%                    * Raw .dat filename containing continuous Y X T data
%                      (LOW-RAM MODE, see Notes)
%                    * UMT struct (any kind, image or non-image) with one
%                      entry that has a T axis as its last axis, e.g.
%                      Y-X-T or ROI-T
%                    * .umt filename containing such a UMT struct
%       SaveFolder - Folder containing events.mat.
%
%   Name-Value parameters:
%       FrameRateHz - Frame rate of DATA (Hz), used to convert event times to
%                     frames. PipelineManager injects it from the data; a .dat
%                     input's header provides it otherwise, and a UMT entry's
%                     meta.FrameRateHz. In-RAM arrays need it explicitly;
%                     AcqInfos.mat is not used.
%
%   Output:
%       outData    - For array and .dat inputs: numeric Y x X x T x E array,
%                    one slice per event instance of events.mat (ignored ones
%                    included), saved as .dat by PipelineManager (.dat header
%                    Phase 8c); in LOW-RAM MODE (.dat filename input) the full
%                    path of the Y-X-T-E .dat file written in SaveFolder
%                    instead. The .dat stores no labels:
%                    resolveDatEventMapping matches its E axis to events.mat,
%                    whose selection flags apply at display and processing
%                    time. For UMT inputs: a UMT struct of the same kind whose
%                    single entry keeps its name and gets an E axis appended
%                    (Y-X-T -> Y-X-T-E, ROI-T -> ROI-T-E), with the shared
%                    eventInfo; saved as .umt. Labels of the other axes are
%                    kept; T labels are dropped (the trial length differs).
%       metaData   - For array and .dat inputs: struct with dimNames
%                    {'Y','X','T','E'}, used by PipelineManager to save the
%                    .dat. Empty for UMT inputs.
%
%   Notes:
%       - LOW-RAM MODE: for a .dat filename the continuous recording is
%         never loaded. The trial frames are computed once (same rules as
%         the in-RAM split, including the crop to the shortest trial), then
%         each trial is read from the input and written to the output .dat,
%         so only one trial is in memory. The output is written to
%         "dataByEv.dat" in SaveFolder and its path is returned.
%       - In-RAM mode: the split algorithm is delegated to
%         EventsManager.splitDataByEvents(...).
%       - UMT inputs (struct or .umt file) are processed in RAM and keep a
%         UMT output. Exactly one entry must have a T axis (as its last
%         axis); only that entry is split and written. Already event-split
%         data (an entry with T and E) are rejected.
%       - Every event instance of events.mat is saved, including the ones the
%         user removed (ignored), so the E axis maps one to one onto
%         events.mat (.dat header Phase 8b). Removal affects display and
%         processing, not saving: functions that reduce this output over
%         events exclude the ignored instances (Phase 8c).

% Default output for pipeline management:
default_Output = 'dataByEv.dat';
metaData = struct();

if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) && ...
        strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outData = localPipelineInfo();
    return
end

assert(nargin >= 2, ...
    'Umitoolbox:split_data_by_event:notEnoughInputs', ...
    'split_data_by_event requires both data and SaveFolder inputs.');

assert((ischar(SaveFolder) || (isstring(SaveFolder) && isscalar(SaveFolder))) && isfolder(SaveFolder), ...
    'Umitoolbox:split_data_by_event:invalidSaveFolder', ...
    'SaveFolder must be an existing folder.');

assert(isfile(fullfile(SaveFolder, 'events.mat')), ...
    'Umitoolbox:split_data_by_event:missingEventsFile', ...
    'The "events.mat" file is missing in "%s".', SaveFolder);

opts = inputParser;
opts.FunctionName = mfilename;
addParameter(opts, 'FrameRateHz', []);
parse(opts, varargin{:});

% LOW-RAM MODE: a .dat filename is streamed trial by trial; the continuous
% recording is never loaded and the result is written straight to a .dat file.
datFile = iResolveDatFile(data, SaveFolder);
if ~isempty(datFile)
    frameRateHz = resolveDataInfoValue('frameRateHz', opts.Results.FrameRateHz, datFile, ...
        mfilename);
    outData = iSplitDatFile(datFile, SaveFolder, EventsManager(SaveFolder), ...
        frameRateHz, default_Output);
    metaData = struct('dimNames', {{'Y','X','T','E'}});
    return
end

src = iResolveInput(data, SaveFolder);
frameRateHz = resolveDataInfoValue('frameRateHz', opts.Results.FrameRateHz, [], ...
    mfilename, 'OwnValue', src.ownRate, 'OwnSource', 'the UMT entry meta.FrameRateHz');
ev = EventsManager(SaveFolder);

dataByEv = ev.splitDataByEvents(src.dataYXT, 'FrameRateHz', frameRateHz, ...
    'IncludeIgnored', true);
assert(ndims(dataByEv) == 4, ...
    'Umitoolbox:split_data_by_event:invalidSplitOutput', ...
    'EventsManager.splitDataByEvents returned unexpected data dimensions.');

if src.isUMT
    % UMT input keeps a UMT output (same kind and entry name, with an E axis
    % appended), carrying the event metadata and the labels of the other axes.
    dataByEv = reshape(dataByEv, [src.leadSize, size(dataByEv, 3), size(dataByEv, 4)]);
    outData = genUMTStruct(single(dataByEv), ...
        'kind', src.kind, ...
        'entryName', src.entryName, ...
        'dimNames', [src.dims, {'E'}], ...
        'labels', src.labels, ...
        'meta', struct('FrameRateHz', frameRateHz));
    outData = appendUMTEventInfo(outData, ...
        'eventInfo', ev.exportEventInfo('FrameRateHz', frameRateHz, 'IncludeIgnored', true), ...
        'overwrite', true);
    validateUMTStruct(outData);
    metaData = struct();
    return
end

% Array input: event-split image data are saved as .dat (Phase 8c).
outData = single(dataByEv);
metaData = struct('dimNames', {{'Y','X','T','E'}});

    function info = localPipelineInfo()
        info = PipelineManager.createPipelineInfo( ...
            'split_data_by_event', ...
            'Split continuous image data into event instances using events.mat.');

        info = PipelineManager.addInput(info, ...
            'data', ...
            {'ImageTimeSeries', 'ProcessedData'}, ...
            'Continuous image time series input.', ...
            'kind', 'input', ...
            'position', 1, ...
            'callType', 'positional', ...
            'isData', true, ...
            'supportsFile', true, ...
            'dataMode', 'either');

        info = PipelineManager.addInput(info, ...
            'SaveFolder', ...
            'SaveFolder', ...
            'Folder containing events.mat.', ...
            'kind', 'input', ...
            'position', 2, ...
            'callType', 'positional', ...
            'isData', false);

        info = PipelineManager.addInput(info, ...
            'FrameRateHz', ...
            'sourceInfo', ...
            'Frame rate of the input data (Hz), injected from the data.', ...
            'kind', 'sourceInfo', ...
            'sourceField', 'frameRateHz', ...
            'required', false);

        info = PipelineManager.addOutput(info, ...
            'outData', ...
            'ProcessedData', ...
            'data', ...
            ['Event-split image data (Y-X-T-E), one slice per event instance: ' ...
             '.dat for array/.dat input, .umt for UMT input.'], ...
            default_Output, ...
            1, ...
            'isData', true);

        info = PipelineManager.addOutput(info, ...
            'metaData', ...
            'metaData', ...
            'data', ...
            'Axes of the .dat output (empty for UMT input).', ...
            '', ...
            2, ...
            'isData', false);
    end
end

function datFile = iResolveDatFile(data, SaveFolder)
%IRESOLVEDATFILE Full path of a .dat filename input, or '' for any other input.

datFile = '';
if ~(ischar(data) || (isstring(data) && isscalar(data)))
    return
end

inPath = char(string(data));
if ~isfile(inPath) && isfile(fullfile(SaveFolder, inPath))
    inPath = fullfile(SaveFolder, inPath);
end
[~, ~, ext] = fileparts(inPath);
if strcmpi(ext, '.dat') && isfile(inPath)
    datFile = inPath;
end
end

function outFile = iSplitDatFile(inFile, SaveFolder, ev, frameRateHz, defaultOutput)
%ISPLITDATFILE Split a continuous Y-X-T .dat file into a Y-X-T-E .dat, one trial at a time.
%
% The trial frames come from EventsManager.getFrameMatrix and are cropped to
% the shortest trial exactly as EventsManager.splitDataByEvents does, so the
% result equals the in-RAM split. Only one trial is held in memory.

info = loadMetaData(inFile);
assertDatLayout(info, {{'Y','X','T'}}, 'split_data_by_event');

frMat = ev.getFrameMatrix(datAxisSize(info, 'T'), '', [], ...
    'FrameRateHz', frameRateHz, 'IncludeIgnored', true);

firstNaNCol = find(any(isnan(frMat), 1), 1, 'first');
if ~isempty(firstNaNCol)
    frMat = frMat(:, 1:firstNaNCol - 1);
end
assert(~isempty(frMat), ...
    'Umitoolbox:split_data_by_event:noEventFrames', ...
    'No event instance of events.mat has frames inside the recording.');

nTrials = size(frMat, 1);
nFramesPerTrial = size(frMat, 2);

slabIn = spatialSlabIO('open', inFile, 'Info', info);
cIn = onCleanup(@() spatialSlabIO('close', slabIn));
Ny = slabIn.Ny;
Nx = slabIn.Nx;

% Write through a scratch file so the declared output only appears once the
% run has completed, and so the input may be the file it would overwrite.
outFile = fullfile(SaveFolder, defaultOutput);
[~, outStem, outExt] = fileparts(defaultOutput);
tmpFile = fullfile(SaveFolder, [outStem, '_writing', outExt]);
hdr = datHeaderFromInfo(info, outStem, 'dataClass', 'single', ...
    'dimNames', {'Y','X','T','E'}, ...
    'dimSizes', [Ny, Nx, nFramesPerTrial, nTrials]);
slabOut = spatialSlabIO('create', tmpFile, hdr);
cOut = onCleanup(@() spatialSlabIO('close', slabOut));

for iTrial = 1:nTrials
    slab = single(spatialSlabIO('read', slabIn, 1:Nx, frMat(iTrial, :)));
    spatialSlabIO('write', slabOut, 1:Nx, slab, ...
        (iTrial - 1) * nFramesPerTrial + (1:nFramesPerTrial));
end

spatialSlabIO('finalize', slabOut);
clear cIn cOut; % close both files before the move below

[moveOk, moveMsg] = movefile(tmpFile, outFile, 'f');
assert(moveOk, 'Umitoolbox:split_data_by_event:outputMoveFailed', ...
    'Failed to move "%s" onto "%s": %s', tmpFile, outFile, moveMsg);
end

function src = iResolveInput(data, SaveFolder)
%IRESOLVEINPUT Resolve in-RAM input forms into data with a trailing time axis.
%
% A .dat filename never gets here (LOW-RAM MODE handles it). Accepted: a
% numeric Y X T array, a UMT struct, or a .umt filename. A UMT may be of any
% kind (image, roi); the entry to split is the one entry that has a T axis,
% which must be its last axis. The entry is reshaped to Y' x X' x T, the
% layout EventsManager.splitDataByEvents takes (leadSize restores it).

src = struct();
src.dataYXT = [];
src.ownRate = [];    % UMT entry meta.FrameRateHz
src.isUMT = false;
src.kind = '';
src.entryName = '';
src.dims = {};
src.leadSize = [];
src.labels = struct();

if isnumeric(data)
    assert(ndims(data) == 3, ...
        'Umitoolbox:split_data_by_event:invalidNumericInput', ...
        'Numeric input data must be a 3D array with dimensions Y, X, T.');
    src.dataYXT = single(data);
    return
end

if ischar(data) || (isstring(data) && isscalar(data))
    inPath = char(string(data));
    if ~isfile(inPath)
        altPath = fullfile(SaveFolder, inPath);
        if isfile(altPath)
            inPath = altPath;
        end
    end
    assert(isfile(inPath), ...
        'Umitoolbox:split_data_by_event:missingInputFile', ...
        'Input file was not found: %s', inPath);

    [~,~,ext] = fileparts(inPath);
    assert(strcmpi(ext, '.umt'), ...
        'Umitoolbox:split_data_by_event:unsupportedInputFile', ...
        'Unsupported input file extension: %s', ext);
    data = loadData(inPath);
end

assert(isstruct(data) && isscalar(data), ...
    'Umitoolbox:split_data_by_event:invalidInputType', ...
    'Input data must be numeric YXT, raw .dat, UMT struct, or .umt file.');

validateUMTStruct(data, 'requireEventInfo', false);

[entry, entryName] = iSelectUMTTimeEntry(data);
dims = cellstr(string(entry.dimNames));
value = entry.value;

% Every UMT pattern with a T axis ends in T (before E); the split output
% then appends E, matching the schema patterns {..., 'T', 'E'}.
assert(strcmp(dims{end}, 'T'), ...
    'Umitoolbox:split_data_by_event:unsupportedTimeAxis', ...
    'The T axis of UMT entry "%s" must be its last axis (dimNames: %s).', ...
    entryName, strjoin(dims, ', '));

nDimsLead = numel(dims) - 1;
sz = size(value, 1:numel(dims));
src.isUMT = true;
src.kind = char(string(data.kind));
src.entryName = entryName;
src.dims = dims;
src.leadSize = sz(1:nDimsLead);
src.dataYXT = single(reshape(value, sz(1), prod(sz(2:nDimsLead)), sz(end)));

if isfield(entry, 'meta') && isstruct(entry.meta) && ...
        isfield(entry.meta, 'FrameRateHz') && ~isempty(entry.meta.FrameRateHz)
    src.ownRate = double(entry.meta.FrameRateHz);
end

% Labels describe the other axes and stay valid; T labels do not (the trial
% length differs from the recording length) and E has none yet.
if isfield(data, 'labels') && isstruct(data.labels)
    src.labels = rmfield(data.labels, intersect(fieldnames(data.labels), {'T', 'E'}));
end
end

function [entry, entryName] = iSelectUMTTimeEntry(umt)
%ISELECTUMTTIMEENTRY Select the one entry of a UMT that has a T axis.

entryNames = fieldnames(umt.data);
validNames = {};
sawEventSplit = false;

for iEntry = 1:numel(entryNames)
    thisEntry = umt.data.(entryNames{iEntry});
    if ~isfield(thisEntry, 'dimNames')
        continue
    end
    dims = cellstr(string(thisEntry.dimNames));
    if ~any(strcmp(dims, 'T'))
        continue
    end
    if any(strcmp(dims, 'E'))
        sawEventSplit = true;
        continue
    end
    validNames{end+1} = entryNames{iEntry}; %#ok<AGROW>
end

assert(~isempty(validNames) || ~sawEventSplit, ...
    'Umitoolbox:split_data_by_event:alreadyEventSplit', ...
    'The UMT input is already event-split (its entry with a T axis also has an E axis).');

assert(~isempty(validNames), ...
    'Umitoolbox:split_data_by_event:noTimeDimension', ...
    'The UMT input has no entry with a T (time) axis.');

assert(isscalar(validNames), ...
    'Umitoolbox:split_data_by_event:multipleCompatibleUMTEntries', ...
    ['Multiple entries with a T axis were found in the UMT input (%s). ' ...
     'The current version can split only one entry.'], strjoin(validNames, ', '));

entryName = validNames{1};
entry = umt.data.(entryName);
end
