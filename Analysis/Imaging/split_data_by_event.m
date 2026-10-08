function outData = split_data_by_event(data, SaveFolder, varargin)
%SPLIT_DATA_BY_EVENT Split continuous data with a time axis into event trials.
%
%   outData = split_data_by_event(data, SaveFolder)
%   outData = split_data_by_event(data, SaveFolder, 'FrameRateHz', rate)
%
%   This function uses the event definitions stored in "events.mat" to
%   split a continuous recording into event instances: an E (event
%   instance) axis is appended after the T axis. The output format follows
%   the input: .dat for a .dat input, a UMT for a UMT input (see Output).
%
%   Inputs:
%       data       - One of (the data must have a T axis):
%                    * .dat filename (or path) of a continuous Y-X-T
%                      recording (LOW-RAM MODE, see Notes)
%                    * UMT struct (any kind, image or non-image) with a
%                      single entry that has a T axis as its last axis, e.g.
%                      Y-X-T or ROI-T
%                    * .umt filename containing such a UMT struct
%                    Arrays are not supported.
%       SaveFolder - Folder containing events.mat; also the folder of the
%                    .dat output.
%
%   Name-Value parameters:
%       FrameRateHz - Frame rate of DATA (Hz), used to convert event times to
%                     frames. PipelineManager injects it from the data; a .dat
%                     input's header provides it otherwise, and a UMT entry's
%                     meta.FrameRateHz. AcqInfos.mat is not used.
%
%   Output:
%       outData    - For a .dat input: the full path of the Y-X-T-E .dat file
%                    written in SaveFolder ("dataByEv.dat"), one slice per
%                    event instance of events.mat (ignored ones included).
%                    The .dat stores no labels: resolveDatEventMapping
%                    matches its E axis to events.mat, whose selection flags
%                    apply at display and processing time. For UMT inputs: a
%                    UMT struct of the same kind whose single entry keeps its
%                    name and gets an E axis appended (Y-X-T -> Y-X-T-E,
%                    ROI-T -> ROI-T-E), with the shared eventInfo; saved as
%                    .umt. Labels of the other axes are kept; T labels are
%                    dropped (the trial length differs).
%
%   Notes:
%       - LOW-RAM MODE: for a .dat filename the continuous recording is
%         never loaded. The trial frames are computed once (same rules as
%         EventsManager.splitDataByEvents, including the crop to the
%         shortest trial), then each trial is read from the input and
%         written to the output .dat in X slabs sized from the available
%         RAM. The layout is checked first: only Y-X-T is accepted
%         (event-split data are already split).
%       - UMT inputs (struct or .umt file) are processed in RAM, by
%         EventsManager.splitDataByEvents, and keep a UMT output. The UMT
%         must hold a single entry, with a T axis as its last axis. Already
%         event-split data (an entry with T and E) are rejected.
%       - Every event instance of events.mat is saved, including the ones the
%         user removed (ignored), so the E axis maps one to one onto
%         events.mat (.dat header Phase 8b). Removal affects display and
%         processing, not saving: functions that reduce this output over
%         events exclude the ignored instances (Phase 8c).

% Default output for pipeline management:
default_Output = 'dataByEv.dat';

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
    % The layout is checked before anything else is resolved: a layout
    % without a T axis, or an already event-split one, has no frame rate or
    % events to look up.
    datInfo = loadMetaData(datFile);
    assertDatLayout(datInfo, {{'Y','X','T'}}, 'split_data_by_event');
    frameRateHz = resolveDataInfoValue('frameRateHz', opts.Results.FrameRateHz, datFile, ...
        mfilename);
    outData = iSplitDatFile(datFile, datInfo, SaveFolder, EventsManager(SaveFolder), ...
        frameRateHz, default_Output);
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

    function info = localPipelineInfo()
        info = PipelineManager.createPipelineInfo( ...
            'split_data_by_event', ...
            'Split continuous data with a time axis into event instances using events.mat.');
        info.version = '2.0.0';

        info = PipelineManager.addInput(info, ...
            'data', ...
            {'ImageTimeSeries', 'ProcessedData'}, ...
            'Continuous Y-X-T .dat file, or a .umt file / UMT struct with one entry that has a T axis.', ...
            'kind', 'input', ...
            'position', 1, ...
            'callType', 'positional', ...
            'isData', true, ...
            'supportsFile', true, ...
            'dataMode', 'file');

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
            ['Event-split data (E axis appended after T), one slice per event ' ...
             'instance: .dat (Y-X-T-E) for .dat input, .umt for UMT input.'], ...
            default_Output, ...
            1, ...
            'isData', true);
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

function outFile = iSplitDatFile(inFile, info, SaveFolder, ev, frameRateHz, defaultOutput)
%ISPLITDATFILE Split a continuous Y-X-T .dat file into a Y-X-T-E .dat, one trial at a time.
%
% The trial frames come from EventsManager.getFrameMatrix and are cropped to
% the shortest trial exactly as EventsManager.splitDataByEvents does, so the
% result equals the in-RAM split. Only one trial, in X slabs sized from the
% available RAM, is held in memory.

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

% One trial's frames are read and written in X slabs; a slab and its single
% copy are resident at once.
nChunks = calculateMaxChunkSize(double(Ny) * double(Nx) * nFramesPerTrial * ...
    double(getByteSize('single')), 2, 0.2);
chunkX = max(1, ceil(Nx / nChunks));
nChunks = ceil(Nx / chunkX);

for iTrial = 1:nTrials
    outFrames = (iTrial - 1) * nFramesPerTrial + (1:nFramesPerTrial);
    for c = 1:nChunks
        xIdx = ((c-1) * chunkX + 1):min(c * chunkX, Nx);
        slab = single(spatialSlabIO('read', slabIn, xIdx, frMat(iTrial, :)));
        spatialSlabIO('write', slabOut, xIdx, slab, outFrames);
    end
    fprintf('Trial %i/%i written\n', iTrial, nTrials)
end

spatialSlabIO('finalize', slabOut);
clear cIn cOut; % close both files before the move below

[moveOk, moveMsg] = movefile(tmpFile, outFile, 'f');
assert(moveOk, 'Umitoolbox:split_data_by_event:outputMoveFailed', ...
    'Failed to move "%s" onto "%s": %s', tmpFile, outFile, moveMsg);
end

function src = iResolveInput(data, SaveFolder)
%IRESOLVEINPUT Resolve a UMT input (struct or .umt file) into data with a trailing time axis.
%
% A .dat filename never gets here (LOW-RAM MODE handles it). Accepted: a UMT
% struct or a .umt filename. A UMT may be of any kind (image, roi) but holds
% a single entry, which has a T axis as its last axis. The entry is reshaped
% to Y' x X' x T, the layout EventsManager.splitDataByEvents takes (leadSize
% restores it).

src = struct();
src.dataYXT = [];
src.ownRate = [];    % UMT entry meta.FrameRateHz
src.kind = '';
src.entryName = '';
src.dims = {};
src.leadSize = [];
src.labels = struct();

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
    'Input data must be a .dat file, a UMT struct, or a .umt file (arrays are not supported).');

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
%ISELECTUMTTIMEENTRY Return the single entry of a UMT, which must have a T axis.

entryNames = fieldnames(umt.data);

assert(isscalar(entryNames), ...
    'Umitoolbox:split_data_by_event:multipleCompatibleUMTEntries', ...
    ['The UMT input has %d entries (%s). The current version can split ' ...
     'only a UMT with a single entry.'], numel(entryNames), strjoin(entryNames, ', '));

entryName = entryNames{1};
entry = umt.data.(entryName);
dims = cellstr(string(entry.dimNames));

assert(any(strcmp(dims, 'T')), ...
    'Umitoolbox:split_data_by_event:noTimeDimension', ...
    'The UMT input has no T (time) axis.');

assert(~any(strcmp(dims, 'E')), ...
    'Umitoolbox:split_data_by_event:alreadyEventSplit', ...
    'The UMT input is already event-split (its entry has both T and E axes).');
end
