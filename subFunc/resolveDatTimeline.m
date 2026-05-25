function timelineInfo = resolveDatTimeline(fileOrLength, AcqInfoStream, varargin)
%RESOLVEDATTIMELINE Resolve a .dat temporal length against AcqInfos metadata.
%
%   timelineInfo = resolveDatTimeline(datFile, AcqInfoStream)
%   timelineInfo = resolveDatTimeline(T, AcqInfoStream)
%   timelineInfo = resolveDatTimeline(T, AcqInfoStream, 'DatFile', datFile)
%
%   A .dat file is valid only when its temporal length matches
%   one known imported/base timeline described by AcqInfos.mat. This helper
%   validates that rule and returns the resolved frame rate.
%
%   Resolution order:
%       1) Exact DatFile match in AcqInfoStream.ImportedChannels
%       2) Unique/safe Length match in ImportedChannels or top-level Length
%       3) Error when no match or ambiguous conflicting frame rates
%
%   Input fileOrLength:
%       - .dat file path: T is inferred from file size using top-level
%         Height/Width and Datatype, defaulting to single precision.
%       - Numeric scalar: interpreted directly as T.

p = inputParser;
p.FunctionName = 'resolveDatTimeline';

addRequired(p, 'fileOrLength', @(x) (isnumeric(x) && isscalar(x)) || ischar(x) || isstring(x));
addRequired(p, 'AcqInfoStream', @(x) isstruct(x) && isscalar(x));
addParameter(p, 'DatFile', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'ThrowError', true, @(x) islogical(x) && isscalar(x));

parse(p, fileOrLength, AcqInfoStream, varargin{:});

throwError = p.Results.ThrowError;
datFile = char(string(p.Results.DatFile));

try
    if isnumeric(fileOrLength)
        actualLength = double(fileOrLength);
    else
        filePath = char(string(fileOrLength));
        if isempty(datFile)
            datFile = filePath;
        end
        actualLength = iInferDatLengthFromFile(filePath, AcqInfoStream);
    end

    validateattributes(actualLength, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive', 'integer'}, ...
        'resolveDatTimeline', 'Length');

    if ~isempty(datFile)
        [~, datBase, datExt] = fileparts(datFile);
        if isempty(datExt)
            datExt = '.dat';
        end
        datFile = [datBase, datExt];
    end

    importedChannels = iGetImportedChannels(AcqInfoStream);

    % Exact file match is authoritative for imported channels.
    if ~isempty(datFile) && ~isempty(importedChannels)
        idxFile = find(strcmpi({importedChannels.DatFile}, datFile));
        if ~isempty(idxFile)
            if numel(idxFile) > 1
                error('Umitoolbox:resolveDatTimeline:duplicateImportedChannel', ...
                    'Multiple ImportedChannels entries were found for "%s".', datFile);
            end

            entry = importedChannels(idxFile);
            if double(entry.Length) ~= actualLength
                error('Umitoolbox:resolveDatTimeline:lengthMismatch', ...
                    ['File "%s" is registered as an imported channel with Length=%d, ' ...
                     'but the actual/inferred length is %d.'], ...
                    datFile, double(entry.Length), actualLength);
            end

            timelineInfo = iBuildTimelineInfo(actualLength, double(entry.FrameRateHz), ...
                'ImportedChannels.DatFile', idxFile, datFile);
            return
        end
    end

    % Length-based match against base and imported timelines.
    candidates = iBuildTimelineCandidates(AcqInfoStream);
    if isempty(candidates)
        error('Umitoolbox:resolveDatTimeline:noTimelines', ...
            'AcqInfoStream does not define any base/imported timelines.');
    end

    candidateLengths = [candidates.Length];
    idxMatch = find(candidateLengths == actualLength);

    if isempty(idxMatch)
        error('Umitoolbox:resolveDatTimeline:noMatch', ...
            ['Length %d does not match any imported/base timeline in ' ...
             'AcqInfos.mat. Save this output as .umt instead of .dat, ' ...
             'or provide explicit metadata.'], ...
            actualLength);
    end

    freqList = [candidates(idxMatch).FrameRateHz];
    uniqueFreq = unique(freqList);
    if numel(uniqueFreq) ~= 1
        error('Umitoolbox:resolveDatTimeline:ambiguousTimeline', ...
            ['Length %d matches multiple known timelines with different ' ...
             'FrameRateHz values. Provide explicit source metadata or use .umt.'], ...
            actualLength);
    end

    timelineInfo = iBuildTimelineInfo(actualLength, uniqueFreq, ...
        'LengthMatch', idxMatch, datFile);

catch ME
    if throwError
        rethrow(ME)
    end

    timelineInfo = struct();
    timelineInfo.IsValid = false;
    timelineInfo.ErrorIdentifier = ME.identifier;
    timelineInfo.ErrorMessage = ME.message;
end

end

% =========================================================================
% Local helpers
% =========================================================================
function actualLength = iInferDatLengthFromFile(filePath, AcqInfoStream)
%IINFERDATLENGTHFROMFILE Infer T from file size and AcqInfoStream Y/X.

if ~isfile(filePath)
    error('Umitoolbox:resolveDatTimeline:fileNotFound', ...
        'File not found: "%s".', filePath);
end

if ~isfield(AcqInfoStream, 'Height') || ~isfield(AcqInfoStream, 'Width')
    error('Umitoolbox:resolveDatTimeline:missingFrameSize', ...
        'AcqInfoStream must contain Height and Width to infer .dat length.');
end

height = double(AcqInfoStream.Height);
width = double(AcqInfoStream.Width);
validateattributes(height, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive', 'integer'}, ...
    'resolveDatTimeline', 'Height');
validateattributes(width, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive', 'integer'}, ...
    'resolveDatTimeline', 'Width');

if isfield(AcqInfoStream, 'Datatype') && ~isempty(AcqInfoStream.Datatype)
    datatype = char(string(AcqInfoStream.Datatype));
else
    datatype = 'single';
end

bytesPerElement = iGetByteSize(datatype);
fileInfo = dir(filePath);
frameBytes = height * width * bytesPerElement;

if mod(fileInfo.bytes, frameBytes) ~= 0
    error('Umitoolbox:resolveDatTimeline:invalidFileLength', ...
        'File size is incompatible with AcqInfoStream Height/Width/Datatype.');
end

actualLength = fileInfo.bytes / frameBytes;

end

function importedChannels = iGetImportedChannels(AcqInfoStream)
%IGETIMPORTEDCHANNELS Return normalized ImportedChannels entries.

importedChannels = struct('DatFile', {}, 'Length', {}, 'FrameRateHz', {});

if ~isfield(AcqInfoStream, 'ImportedChannels') || isempty(AcqInfoStream.ImportedChannels)
    return
end

raw = AcqInfoStream.ImportedChannels(:).';
for iEntry = 1:numel(raw)
    if ~isfield(raw(iEntry), 'DatFile') || ...
            ~isfield(raw(iEntry), 'Length') || ...
            ~isfield(raw(iEntry), 'FrameRateHz')
        error('Umitoolbox:resolveDatTimeline:invalidImportedChannels', ...
            ['Each ImportedChannels entry must contain DatFile, Length, ' ...
             'and FrameRateHz.']);
    end

    datFile = char(string(raw(iEntry).DatFile));
    [~, datBase, datExt] = fileparts(datFile);
    if isempty(datExt)
        datExt = '.dat';
    end

    importedChannels(end+1).DatFile = [datBase, datExt]; %#ok<AGROW>
    importedChannels(end).Length = double(raw(iEntry).Length);
    importedChannels(end).FrameRateHz = double(raw(iEntry).FrameRateHz);
end

end

function candidates = iBuildTimelineCandidates(AcqInfoStream)
%IBUILDTIMELINECANDIDATES Build base/imported timeline candidates.

candidates = struct('Length', {}, 'FrameRateHz', {}, 'SourceType', {}, 'SourceIndex', {});

if isfield(AcqInfoStream, 'Length') && ~isempty(AcqInfoStream.Length) && ...
        isfield(AcqInfoStream, 'FrameRateHz') && ~isempty(AcqInfoStream.FrameRateHz)
    candidates(end+1) = struct( ... %#ok<AGROW>
        'Length', double(AcqInfoStream.Length), ...
        'FrameRateHz', double(AcqInfoStream.FrameRateHz), ...
        'SourceType', 'AcqInfoStream.Length', ...
        'SourceIndex', 0);
end

importedChannels = iGetImportedChannels(AcqInfoStream);
for iEntry = 1:numel(importedChannels)
    candidates(end+1) = struct( ... %#ok<AGROW>
        'Length', double(importedChannels(iEntry).Length), ...
        'FrameRateHz', double(importedChannels(iEntry).FrameRateHz), ...
        'SourceType', 'ImportedChannels.Length', ...
        'SourceIndex', iEntry);
end

end

function timelineInfo = iBuildTimelineInfo(len, freq, sourceType, sourceIndex, datFile)
%IBUILDTIMELINEINFO Create timeline-resolution output structure.

timelineInfo = struct();
timelineInfo.IsValid = true;
timelineInfo.Length = double(len);
timelineInfo.datLength = double(len);
timelineInfo.FrameRateHz = double(freq);
timelineInfo.Freq = double(freq);
timelineInfo.SourceType = sourceType;
timelineInfo.SourceIndex = sourceIndex;
timelineInfo.DatFile = datFile;

end

function nBytes = iGetByteSize(datatype)
%IGETBYTESIZE Return byte size for supported numeric datatypes.

switch lower(char(string(datatype)))
    case {'single', 'float32'}
        nBytes = 4;
    case {'double', 'float64'}
        nBytes = 8;
    case {'uint8', 'int8', 'logical'}
        nBytes = 1;
    case {'uint16', 'int16'}
        nBytes = 2;
    case {'uint32', 'int32'}
        nBytes = 4;
    case {'uint64', 'int64'}
        nBytes = 8;
    otherwise
        error('Umitoolbox:resolveDatTimeline:unsupportedDatatype', ...
            'Unsupported datatype: "%s".', char(string(datatype)));
end

end
