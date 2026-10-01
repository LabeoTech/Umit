function timelineInfo = resolveDatTimeline(fileOrLength, AcqInfoStream, varargin)
%RESOLVEDATTIMELINE Resolve a temporal length against AcqInfos metadata.
%
%   timelineInfo = resolveDatTimeline(T, AcqInfoStream)
%   timelineInfo = resolveDatTimeline(T, AcqInfoStream, 'DatFile', datFile)
%
%   Matches a temporal length T against the imported/base timelines
%   described by AcqInfos.mat and returns the resolved frame rate. Used for
%   data without a file source (in-memory arrays); .dat files describe
%   their own frame rate (see loadMetaData).
%
%   Resolution order:
%       1) Exact DatFile match in AcqInfoStream.ImportedChannels
%       2) Unique/safe Length match in ImportedChannels or top-level Length
%       3) Error when no match or ambiguous conflicting frame rates
%
%   Input fileOrLength:
%       - Numeric scalar T. (Passing a .dat file path, which inferred T
%         from AcqInfos.mat Height/Width, is no longer supported.)

p = inputParser;
p.FunctionName = 'resolveDatTimeline';

addRequired(p, 'fileOrLength', @(x) isnumeric(x) && isscalar(x));
addRequired(p, 'AcqInfoStream', @(x) isstruct(x) && isscalar(x));
addParameter(p, 'DatFile', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'ThrowError', true, @(x) islogical(x) && isscalar(x));

parse(p, fileOrLength, AcqInfoStream, varargin{:});

throwError = p.Results.ThrowError;
datFile = char(string(p.Results.DatFile));

try
    actualLength = double(fileOrLength);

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
function importedChannels = iGetImportedChannels(AcqInfoStream)
%IGETIMPORTEDCHANNELS Return normalized ImportedChannels entries.

importedChannels = struct('DatFile', {}, 'Length', {}, 'FrameRateHz', {});

if ~isfield(AcqInfoStream, 'ImportedChannels') || isempty(AcqInfoStream.ImportedChannels)
    return
end

raw = AcqInfoStream.ImportedChannels(:).';
importedChannels = repmat(importedChannels, 1, numel(raw));
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

    importedChannels(iEntry).DatFile = [datBase, datExt];
    importedChannels(iEntry).Length = double(raw(iEntry).Length);
    importedChannels(iEntry).FrameRateHz = double(raw(iEntry).FrameRateHz);
end

end

function candidates = iBuildTimelineCandidates(AcqInfoStream)
%IBUILDTIMELINECANDIDATES Build base/imported timeline candidates.

hasBaseTimeline = isfield(AcqInfoStream, 'Length') && ...
    ~isempty(AcqInfoStream.Length) && isfield(AcqInfoStream, 'FrameRateHz') && ...
    ~isempty(AcqInfoStream.FrameRateHz);
importedChannels = iGetImportedChannels(AcqInfoStream);
candidates = repmat(struct('Length', 0, 'FrameRateHz', 0, ...
    'SourceType', '', 'SourceIndex', 0), 1, ...
    double(hasBaseTimeline) + numel(importedChannels));
candidateIdx = 0;

if hasBaseTimeline
    candidateIdx = candidateIdx + 1;
    candidates(candidateIdx) = struct( ...
        'Length', double(AcqInfoStream.Length), ...
        'FrameRateHz', double(AcqInfoStream.FrameRateHz), ...
        'SourceType', 'AcqInfoStream.Length', ...
        'SourceIndex', 0);
end

for iEntry = 1:numel(importedChannels)
    candidateIdx = candidateIdx + 1;
    candidates(candidateIdx) = struct( ...
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
