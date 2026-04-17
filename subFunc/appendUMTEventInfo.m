function out = appendUMTEventInfo(umt, varargin)
%APPENDUMTEVENTINFO Append or replace shared event metadata in a .umt structure.
%
%   Supported call patterns:
%       out = appendUMTEventInfo(umt, ...
%           'eventID', eventID, ...
%           'repetitionIndex', repetitionIndex, ...
%           'eventName', eventName, ...
%           'eventAxisMode', eventAxisMode)
%
%       out = appendUMTEventInfo(umt, 'eventInfo', eventInfoStruct)
%
%   Notes:
%       - baselinePeriod may be stored in top-level eventInfo.
%       - When an eventInfo struct is provided, missing fields commonly
%         omitted by EventsManager.exportEventInfo are derived when possible.

errID = 'Umitoolbox:appendUMTEventInfo:invalidInput';

validateUMTStruct(umt, 'requireEventInfo', false);
schema = getUMTSchema(umt.version);

p = inputParser;
p.FunctionName = 'appendUMTEventInfo';
addParameter(p, 'eventID', [], @(x) isnumeric(x) && isvector(x));
addParameter(p, 'repetitionIndex', [], @(x) isnumeric(x) && isvector(x));
addParameter(p, 'eventName', [], @(x) isstring(x) || iscell(x));
addParameter(p, 'eventAxisMode', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'eventInfo', struct(), @(x) isstruct(x) && isscalar(x));
addParameter(p, 'overwrite', false, @(x) islogical(x) && isscalar(x));
parse(p, varargin{:});

entryNames = fieldnames(umt.data);
eLengths = [];
for iEntry = 1:numel(entryNames)
    entryName = entryNames{iEntry};
    entryData = umt.data.(entryName);
    if ~isstruct(entryData) || ~isscalar(entryData)
        error(errID, ...
            'Operation aborted. Entry "%s" must be a scalar struct.', entryName);
    end
    if ~isfield(entryData, 'value') || ~isfield(entryData, 'dimNames')
        error(errID, ...
            ['Operation aborted. Entry "%s" must contain fields "value" ' ...
             'and "dimNames".'], ...
            entryName);
    end
    dimNames = iNormalizeDimNames(entryData.dimNames, schema, errID, entryName);
    dimSizes = iGetDeclaredDimensionSizes(entryData.value, dimNames, errID, entryName);
    if isempty(dimSizes)
        continue
    end
    idxE = find(strcmp(dimNames, 'E'), 1, 'first');
    if ~isempty(idxE)
        eLengths(end+1) = dimSizes(idxE); %#ok<AGROW>
    end
end

if isempty(eLengths)
    error(errID, ...
        ['Operation aborted. appendUMTEventInfo can only be used when at ' ...
         'least one entry uses the "E" dimension.']);
end
if numel(unique(eLengths)) ~= 1
    error(errID, ...
        ['Operation aborted. All entries that use the "E" dimension must ' ...
         'share the same E length.']);
end

if iHasNameValue(varargin, 'eventInfo') && ~isempty(fieldnames(p.Results.eventInfo))
    eventInfo = iNormalizeEventInfoFromStruct(p.Results.eventInfo, eLengths(1), schema, errID);
else
    requiredNV = {'eventID','repetitionIndex','eventName','eventAxisMode'};
    for iReq = 1:numel(requiredNV)
        if ~iHasNameValue(varargin, requiredNV{iReq})
            error(errID, ...
                'Operation aborted. "%s" must be provided.', requiredNV{iReq});
        end
    end
    eventInfo = iNormalizeEventInfoFromNameValue( ...
        p.Results.eventID, p.Results.repetitionIndex, p.Results.eventName, ...
        p.Results.eventAxisMode, eLengths(1), schema, errID);
end

out = umt;
if isfield(out, 'eventInfo') && ~p.Results.overwrite
    if ~isequal(out.eventInfo, eventInfo)
        error(errID, ...
            ['Operation aborted. "eventInfo" already exists with a ' ...
             'different value. Use overwrite=true to replace it.']);
    end
end
out.eventInfo = eventInfo;
validateUMTStruct(out, 'requireEventInfo', true);
end

function tf = iHasNameValue(args, name)
tf = false;
if mod(numel(args), 2) ~= 0
    return
end
for iArg = 1:2:numel(args)
    key = args{iArg};
    if ischar(key) || (isstring(key) && isscalar(key))
        if strcmpi(char(string(key)), name)
            tf = true;
            return
        end
    end
end
end

function dimNames = iNormalizeDimNames(dimNamesIn, schema, errID, entryName)
allowedDims = schema.allowedDims;
if isempty(dimNamesIn)
    dimNames = {};
    return
end
if isstring(dimNamesIn)
    if ~isvector(dimNamesIn)
        error(errID, ...
            'Operation aborted. Entry "%s.dimNames" must be a 1D text array.', ...
            entryName);
    end
    rawDims = cellstr(dimNamesIn(:).');
elseif iscell(dimNamesIn) && isvector(dimNamesIn) && ...
        all(cellfun(@(c) ischar(c) || (isstring(c) && isscalar(c)), dimNamesIn))
    rawDims = cellstr(string(dimNamesIn(:).'));
else
    error(errID, ...
        ['Operation aborted. Entry "%s.dimNames" must be a cell array of ' ...
         'character vectors or a string vector.'], ...
        entryName);
end

dimNames = cell(1, numel(rawDims));
for iDim = 1:numel(rawDims)
    idx = find(strcmpi(rawDims{iDim}, allowedDims), 1, 'first');
    if isempty(idx)
        error(errID, ...
            ['Operation aborted. Entry "%s" contains invalid dimension name ' ...
             '"%s".'], ...
            entryName, rawDims{iDim});
    end
    dimNames{iDim} = allowedDims{idx};
end
end

function dimSizes = iGetDeclaredDimensionSizes(value, dimNames, errID, entryName)
nDimsExpected = numel(dimNames);
if nDimsExpected == 0
    if ~isscalar(value)
        error(errID, ...
            'Operation aborted. Scalar entry "%s" must use dimNames = {}.', ...
            entryName);
    end
    dimSizes = [];
    return
end
if isscalar(value)
    error(errID, ...
        'Operation aborted. Scalar entry "%s" must use dimNames = {}.', ...
        entryName);
end
sz = size(value);
nonSingletonDims = find(sz ~= 1);
if numel(nonSingletonDims) == 1 && nonSingletonDims ~= 1
    error(errID, ...
        ['Operation aborted. Entry "%s" is one-dimensional but not stored ' ...
         'as a column vector.'], ...
        entryName);
end
if numel(sz) < nDimsExpected
    sz(end+1:nDimsExpected) = 1;
elseif numel(sz) > nDimsExpected
    if any(sz(nDimsExpected+1:end) ~= 1)
        error(errID, ...
            ['Operation aborted. Entry "%s.value" has more dimensions than ' ...
             'declared in dimNames, and the extra dimensions are not singleton.'], ...
            entryName);
    end
    sz = sz(1:nDimsExpected);
end
dimSizes = sz;
end

function eventInfo = iNormalizeEventInfoFromStruct(eventInfoIn, eLen, schema, errID)
filtered = struct();
rawFields = fieldnames(eventInfoIn);
keepFields = intersect(rawFields, ...
    [schema.requiredEventInfoFields, schema.optionalEventInfoFields, {'eventNameList'}], ...
    'stable');
for iField = 1:numel(keepFields)
    filtered.(keepFields{iField}) = eventInfoIn.(keepFields{iField});
end

if ~isfield(filtered, 'eventID')
    error(errID, ...
        'Operation aborted. Provided eventInfo struct must contain eventID.');
end

% Default axis mode for exported instance-wise event info.
if ~isfield(filtered, 'eventAxisMode') || isempty(filtered.eventAxisMode)
    filtered.eventAxisMode = 'instances';
end

% Derive repetition indices when missing.
if ~isfield(filtered, 'repetitionIndex') || isempty(filtered.repetitionIndex)
    ids = double(filtered.eventID(:));
    repIdx = zeros(size(ids));
    uniqueIDs = unique(ids, 'stable');
    for iID = 1:numel(uniqueIDs)
        idx = find(ids == uniqueIDs(iID));
        repIdx(idx) = 1:numel(idx);
    end
    filtered.repetitionIndex = repIdx;
end

% Derive names when missing.
if ~isfield(filtered, 'eventName') || isempty(filtered.eventName)
    if isfield(filtered, 'eventNameList') && ~isempty(filtered.eventNameList)
        ids = double(filtered.eventID(:));
        nameList = filtered.eventNameList;
        if isstring(nameList)
            nameList = cellstr(nameList(:));
        elseif iscell(nameList)
            nameList = cellstr(string(nameList(:)));
        else
            nameList = {};
        end

        if ~isempty(nameList) && all(ids >= 1) && all(ids <= numel(nameList))
            filtered.eventName = reshape(string(nameList(ids)), [], 1);
        end
    end

    if ~isfield(filtered, 'eventName') || isempty(filtered.eventName)
        ids = double(filtered.eventID(:));
        filtered.eventName = reshape(string(arrayfun(@(x) sprintf('Event%d', x), ids, 'UniformOutput', false)), [], 1);
    end
end

if isfield(filtered, 'eventNameList')
    filtered = rmfield(filtered, 'eventNameList');
end

eventInfo = iNormalizeEventInfoFromNameValue( ...
    filtered.eventID, filtered.repetitionIndex, filtered.eventName, filtered.eventAxisMode, ...
    eLen, schema, errID);

if isfield(filtered, 'baselinePeriod')
    bp = filtered.baselinePeriod;
    if ~isnumeric(bp) || ~isscalar(bp) || ~isfinite(bp) || bp <= 0
        error(errID, ...
            ['Operation aborted. eventInfo.baselinePeriod must be a ' ...
             'positive numeric scalar when provided.']);
    end
    eventInfo.baselinePeriod = double(bp);
end
end

function eventInfo = iNormalizeEventInfoFromNameValue(eventID, repIdx, evName, axisMode, eLen, schema, errID)
if ~isnumeric(eventID) || ~isvector(eventID) || isempty(eventID) || ...
        any(~isfinite(eventID(:))) || any(eventID(:) < 1) || ...
        any(mod(eventID(:),1) ~= 0)
    error(errID, ...
        ['Operation aborted. "eventID" must be a non-empty vector of ' ...
         'positive integers.']);
end
eventID = eventID(:);
if ~isnumeric(repIdx) || ~isvector(repIdx) || isempty(repIdx) || ...
        any(~isfinite(repIdx(:))) || any(repIdx(:) < 0) || ...
        any(mod(repIdx(:),1) ~= 0)
    error(errID, ...
        ['Operation aborted. "repetitionIndex" must be a non-empty vector ' ...
         'of non-negative integers.']);
end
repIdx = repIdx(:);
if isstring(evName) && isvector(evName)
    evName = cellstr(evName(:));
elseif iscell(evName) && isvector(evName) && ...
        all(cellfun(@(c) ischar(c) || (isstring(c) && isscalar(c)), evName))
    evName = cellstr(string(evName(:)));
else
    error(errID, ...
        ['Operation aborted. "eventName" must be a string vector or a cell ' ...
         'array of character vectors.']);
end
if ~(ischar(axisMode) || (isstring(axisMode) && isscalar(axisMode)))
    error(errID, ...
        'Operation aborted. "eventAxisMode" must be a text scalar.');
end
axisMode = lower(char(string(axisMode)));
if ~ismember(axisMode, schema.allowedEventAxisModes)
    error(errID, ...
        'Operation aborted. Invalid eventAxisMode "%s".', axisMode);
end
nE = numel(eventID);
if nE ~= eLen
    error(errID, ...
        ['Operation aborted. Event info length (%d) does not match the ' ...
         'shared E length (%d).'], ...
        nE, eLen);
end
if numel(repIdx) ~= nE || numel(evName) ~= nE
    error(errID, ...
        ['Operation aborted. eventID, repetitionIndex, and eventName must ' ...
         'all have the same number of elements.']);
end
switch axisMode
    case 'instances'
        if any(repIdx < 1)
            error(errID, ...
                ['Operation aborted. repetitionIndex must contain positive ' ...
                 'integers when eventAxisMode="instances".']);
        end
    case 'aggregated_repetitions'
        if any(repIdx ~= 0)
            error(errID, ...
                ['Operation aborted. repetitionIndex must be all zeros ' ...
                 'when eventAxisMode="aggregated_repetitions".']);
        end
end
eventInfo = struct();
eventInfo.eventID = eventID;
eventInfo.repetitionIndex = repIdx;
eventInfo.eventName = evName;
eventInfo.eventAxisMode = axisMode;
end
