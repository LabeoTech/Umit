function outData = split_data_by_event(data, SaveFolder)
%SPLIT_DATA_BY_EVENT Split continuous image time series into event trials.
%
%   outData = split_data_by_event(data, SaveFolder)
%
%   This function uses the event definitions stored in "events.mat" to
%   split a continuous image time series into event instances. The output is
%   returned as a UMT structure containing event-split image data with
%   dimensions Y X T E, where E is the event-instance axis.
%
%   Inputs:
%       data       - One of:
%                    * Numeric image time series with dimensions Y X T
%                    * Raw .dat filename containing continuous Y X T data
%                    * UMT struct containing one continuous image entry
%                    * .umt filename containing one continuous image entry
%       SaveFolder - Folder containing events.mat and any metadata needed to
%                    resolve file-backed inputs.
%
%   Output:
%       outData    - UMT structure containing one image entry with dimNames
%                    {'Y','X','T','E'} and shared top-level eventInfo.
%
%   Notes:
%       - This function does not implement low-RAM mode.
%       - The underlying split algorithm is delegated to
%         EventsManager.splitDataByEvents(...).
%       - Continuous UMT inputs are supported. Already event-split UMT data
%         are rejected.

% Default output for pipeline management:
default_Output = 'dataByEv.umt'; %#ok<NASGU>

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

src = iResolveInput(data, SaveFolder);
ev = EventsManager(SaveFolder);

dataByEv = ev.splitDataByEvents(src.dataYXT);
assert(ndims(dataByEv) == 4, ...
    'Umitoolbox:split_data_by_event:invalidSplitOutput', ...
    'EventsManager.splitDataByEvents returned unexpected data dimensions.');

entryMeta = struct();
if isfield(src, 'entryMeta') && isstruct(src.entryMeta) && ...
        isfield(src.entryMeta, 'FrameRateHz') && ~isempty(src.entryMeta.FrameRateHz)
    entryMeta.FrameRateHz = src.entryMeta.FrameRateHz;
elseif isfield(src, 'frameRateHz') && ~isempty(src.frameRateHz)
    entryMeta.FrameRateHz = src.frameRateHz;
end

outData = genUMTStruct(single(dataByEv), ...
    'kind', 'image', ...
    'entryName', 'main', ...
    'dimNames', {'Y','X','T','E'}, ...
    'meta', entryMeta, ...
    'SaveFolder', SaveFolder);

outData = appendUMTEventInfo(outData, ...
    'eventInfo', ev.exportEventInfo(), ...
    'overwrite', true);

validateUMTStruct(outData);

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
            'Folder containing events.mat and metadata.', ...
            'kind', 'input', ...
            'position', 2, ...
            'callType', 'positional', ...
            'isData', false);

        info = PipelineManager.addOutput(info, ...
            'outData', ...
            'ProcessedData', ...
            'data', ...
            'Event-split image data in UMT format.', ...
            default_Output, ...
            1, ...
            'isData', true);
    end
end

function src = iResolveInput(data, SaveFolder)
%IRESOLVEINPUT Resolve supported input forms into continuous YXT data.

src = struct();
src.dataYXT = [];
src.entryMeta = struct();
src.frameRateHz = [];

if isnumeric(data)
    assert(ndims(data) == 3, ...
        'Umitoolbox:split_data_by_event:invalidNumericInput', ...
        'Numeric input data must be a 3D array with dimensions Y, X, T.');
    src.dataYXT = single(data);
    return
end

if ischar(data) || (isstring(data) && isscalar(data))
    inPath = char(string(data));
    assert(isfile(inPath), ...
        'Umitoolbox:split_data_by_event:missingInputFile', ...
        'Input file was not found: %s', inPath);

    [~,~,ext] = fileparts(inPath);
    ext = lower(ext);

    switch ext
        case '.dat'
            loaded = loadData(inPath);
            assert(isnumeric(loaded) && ndims(loaded) == 3, ...
                'Umitoolbox:split_data_by_event:invalidDatInput', ...
                'Raw .dat input must resolve to continuous YXT data.');
            src.dataYXT = single(loaded);

            md = load(fullfile(SaveFolder,'AcqInfos.mat'),'AcqInfoStream');
            md = md.AcqInfoStream;
            if isstruct(md)
                if isfield(md, 'FrameRateHz') && ~isempty(md.FrameRateHz)
                    src.frameRateHz = double(md.FrameRateHz);
                elseif isfield(md, 'Freq') && ~isempty(md.Freq)
                    src.frameRateHz = double(md.Freq);
                end
            end
            return

        case '.umt'
            data = loadData(inPath);

        otherwise
            error('Umitoolbox:split_data_by_event:unsupportedInputFile', ...
                'Unsupported input file extension: %s', ext);
    end
end

assert(isstruct(data) && isscalar(data), ...
    'Umitoolbox:split_data_by_event:invalidInputType', ...
    'Input data must be numeric YXT, raw .dat, UMT struct, or .umt file.');

validateUMTStruct(data, 'requireEventInfo', false);
assert(strcmpi(data.kind, 'image'), ...
    'Umitoolbox:split_data_by_event:invalidUMTKind', ...
    'UMT input must be of kind "image".');

entry = iSelectUMTImageEntry(data);
dimNames = entry.dimNames;

assert(isequal(dimNames, {'Y','X','T'}), ...
    'Umitoolbox:split_data_by_event:invalidUMTDims', ...
    'UMT input must contain one continuous image entry with dimNames {''Y'',''X'',''T''}.');

src.dataYXT = single(entry.value);

if isfield(entry, 'meta') && isstruct(entry.meta)
    src.entryMeta = entry.meta;
end

if isfield(src.entryMeta, 'FrameRateHz') && ~isempty(src.entryMeta.FrameRateHz)
    src.frameRateHz = double(src.entryMeta.FrameRateHz);
elseif isfolder(SaveFolder)
    md = load(fullfile(SaveFolder,'AcqInfos.mat'),'AcqInfoStream');
    md = md.AcqInfoStream;
    if isstruct(md)
        if isfield(md, 'FrameRateHz') && ~isempty(md.FrameRateHz)
            src.frameRateHz = double(md.FrameRateHz);
        elseif isfield(md, 'Freq') && ~isempty(md.Freq)
            src.frameRateHz = double(md.Freq);
        end
    end
end
end

function entry = iSelectUMTImageEntry(umt)
%ISELECTUMTIMAGEENTRY Select one compatible continuous image entry from UMT.

entryNames = fieldnames(umt.data);
validNames = {};

for iEntry = 1:numel(entryNames)
    thisEntry = umt.data.(entryNames{iEntry});
    if isfield(thisEntry, 'dimNames') && isequal(thisEntry.dimNames, {'Y','X','T'})
        validNames{end+1} = entryNames{iEntry}; %#ok<AGROW>
    end
end

assert(~isempty(validNames), ...
    'Umitoolbox:split_data_by_event:noCompatibleUMTEntry', ...
    'No compatible continuous image entry with dimNames {''Y'',''X'',''T''} was found in the UMT input.');

assert(numel(validNames) == 1, ...
    'Umitoolbox:split_data_by_event:multipleCompatibleUMTEntries', ...
    ['Multiple compatible continuous image entries were found in the UMT input. ' ...
     'The current version can process only one image entry.']);

entry = umt.data.(validNames{1});
end
