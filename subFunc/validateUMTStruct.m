function validateUMTStruct(umt, varargin)
%VALIDATEUMTSTRUCT Validate UMIT processed-data structure for .umt files.
%
%   validateUMTStruct(umt)
%   validateUMTStruct(umt, 'requireEventInfo', false)
%
%   Required top-level fields:
%       version
%       kind
%       data
%
%   Optional top-level fields:
%       labels
%       eventInfo
%
%   Required fields for each entry in "umt.data":
%       value
%       dimNames
%
%   Optional entry fields:
%       meta
%
%   Notes:
%       - labels are shared top-level display/reference metadata only.
%       - eventInfo is shared top-level event metadata.
%       - entry.meta stores entry-specific processing metadata, such as
%         FrameRateHz.

errID = 'Umitoolbox:validateUMTStruct:invalidInput';

p = inputParser;
p.FunctionName = 'validateUMTStruct';
addParameter(p, 'requireEventInfo', true, @(x) islogical(x) && isscalar(x));
parse(p, varargin{:});

requireEventInfo = p.Results.requireEventInfo;

if ~isstruct(umt) || ~isscalar(umt)
    error(errID, ...
        'Operation aborted. UMT data must be a scalar struct.');
end

if ~isfield(umt, 'version')
    error(errID, ...
        'Operation aborted. UMT data is missing required field "version".');
end

if ~isnumeric(umt.version) || ~isscalar(umt.version)
    error(errID, ...
        'Operation aborted. "version" must be a numeric scalar.');
end

schema = getUMTSchema(umt.version);

topFields = fieldnames(umt);
allowedTopFields = [schema.requiredTopFields, schema.optionalTopFields];

if ~all(ismember(schema.requiredTopFields, topFields))
    error(errID, ...
        'Operation aborted. UMT data is missing one or more required top-level fields.');
end

unknownTopFields = setdiff(topFields, allowedTopFields);
if ~isempty(unknownTopFields)
    error(errID, ...
        'Operation aborted. Unsupported top-level field(s): %s.', ...
        strjoin(unknownTopFields, ', '));
end

kind = iNormalizeKind(umt.kind, schema, errID);

if ~isstruct(umt.data) || ~isscalar(umt.data)
    error(errID, ...
        'Operation aborted. "data" must be a scalar struct.');
end

entryNames = fieldnames(umt.data);
if isempty(entryNames)
    error(errID, ...
        'Operation aborted. "data" must contain at least one entry.');
end

dimUsage = struct();
eLengths = [];
bUsesEvents = false;

for iEntry = 1:numel(entryNames)
    entryName = entryNames{iEntry};
    entryData = umt.data.(entryName);

    if ~isstruct(entryData) || ~isscalar(entryData)
        error(errID, ...
            'Operation aborted. Entry "%s" must be a scalar struct.', entryName);
    end

    entryFields = fieldnames(entryData);
    allowedEntryFields = [schema.requiredEntryFields, schema.optionalEntryFields];

    if ~all(ismember(schema.requiredEntryFields, entryFields))
        error(errID, ...
            'Operation aborted. Entry "%s" is missing one or more required fields.', ...
            entryName);
    end

    unknownEntryFields = setdiff(entryFields, allowedEntryFields);
    if ~isempty(unknownEntryFields)
        error(errID, ...
            'Operation aborted. Entry "%s" contains unsupported field(s): %s.', ...
            entryName, strjoin(unknownEntryFields, ', '));
    end

    value = entryData.value;
    if ~(isnumeric(value) || islogical(value))
        error(errID, ...
            'Operation aborted. Entry "%s.value" must be numeric or logical.', ...
            entryName);
    end

    dimNames = iNormalizeDimNames(entryData.dimNames, schema, errID, entryName);
    dimSizes = iGetDeclaredDimensionSizes(value, dimNames, errID, entryName);

    if ~iIsAllowedPattern(kind, dimNames, schema)
        error(errID, ...
            ['Operation aborted. Entry "%s" has invalid dimNames for kind "%s".'], ...
            entryName, kind);
    end

    if isfield(entryData, 'meta')
        iValidateEntryMeta(entryData.meta, schema, errID, entryName);
    end

    for iDim = 1:numel(dimNames)
        thisDim = dimNames{iDim};
        thisLen = dimSizes(iDim);

        if isfield(dimUsage, thisDim)
            dimUsage.(thisDim)(end+1) = thisLen; %#ok<AGROW>
        else
            dimUsage.(thisDim) = thisLen;
        end
    end

    idxE = find(strcmp(dimNames, 'E'), 1, 'first');
    if ~isempty(idxE)
        bUsesEvents = true;
        eLengths(end+1) = dimSizes(idxE); %#ok<AGROW>
    end
end

if isfield(umt, 'labels')
    labels = iNormalizeLabelsStruct(umt.labels, schema, errID);
    labelDims = fieldnames(labels);
    usedDims  = fieldnames(dimUsage);

    if ~all(ismember(labelDims, usedDims))
        invalidDims = labelDims(~ismember(labelDims, usedDims));
        error(errID, ...
            ['Operation aborted. "labels" contains dimension(s) that are not ' ...
             'used by any data entry: %s.'], ...
            strjoin(invalidDims, ', '));
    end

    for iDim = 1:numel(labelDims)
        thisDim = labelDims{iDim};
        thisLabels = labels.(thisDim);
        thisLengths = dimUsage.(thisDim);

        if any(thisLengths ~= numel(thisLabels))
            error(errID, ...
                ['Operation aborted. labels.%s does not match the size of all ' ...
                 'entries that use dimension "%s".'], ...
                thisDim, thisDim);
        end
    end
end

if bUsesEvents
    if numel(unique(eLengths)) ~= 1
        error(errID, ...
            ['Operation aborted. All entries that use the "E" dimension ' ...
             'must have the same E length when shared top-level eventInfo ' ...
             'is used.']);
    end

    if isfield(umt, 'eventInfo')
        iValidateEventInfoStruct(umt.eventInfo, eLengths(1), schema, errID);
    elseif requireEventInfo
        error(errID, ...
            ['Operation aborted. Top-level "eventInfo" is required when at ' ...
             'least one entry uses the "E" dimension.']);
    end
else
    if isfield(umt, 'eventInfo')
        error(errID, ...
            ['Operation aborted. Top-level "eventInfo" is only allowed when ' ...
             'at least one entry uses the "E" dimension.']);
    end
end

end

function kind = iNormalizeKind(kindIn, schema, errID)
if ~(ischar(kindIn) || (isstring(kindIn) && isscalar(kindIn)))
    error(errID, ...
        'Operation aborted. "kind" must be a character vector or string scalar.');
end
kind = lower(char(string(kindIn)));
if ~ismember(kind, schema.allowedKinds)
    error(errID, ...
        'Operation aborted. Invalid kind "%s".', kind);
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
if nDimsExpected == 1
    % Only a genuinely 1-D declared entry is ambiguous between a row and a
    % column orientation. For nDimsExpected >= 2, the array's own shape
    % positionally assigns a size to every declared dimension (e.g. a
    % single ROI paired with N timepoints is legitimately [1, N]), so this
    % check must not run there.
    nonSingletonDims = find(sz ~= 1);
    if numel(nonSingletonDims) == 1 && nonSingletonDims ~= 1
        error(errID, ...
            ['Operation aborted. Entry "%s" is one-dimensional but not stored ' ...
             'as a column vector.'], ...
            entryName);
    end
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

function tf = iIsAllowedPattern(kind, dimNames, schema)
if isempty(dimNames)
    tf = true;
    return
end
allowedPatterns = schema.allowedPatterns.(kind);
tf = any(cellfun(@(c) isequal(dimNames, c), allowedPatterns));
end

function labels = iNormalizeLabelsStruct(labelsIn, schema, errID)
allowedDims = schema.allowedDims;
if ~isstruct(labelsIn) || ~isscalar(labelsIn)
    error(errID, ...
        'Operation aborted. "labels" must be a scalar struct.');
end
labels = struct();
rawFields = fieldnames(labelsIn);
if isempty(rawFields)
    return
end
canonicalFields = cell(size(rawFields));
for iField = 1:numel(rawFields)
    idx = find(strcmpi(rawFields{iField}, allowedDims), 1, 'first');
    if isempty(idx)
        error(errID, ...
            'Operation aborted. Invalid label dimension "%s".', rawFields{iField});
    end
    canonicalFields{iField} = allowedDims{idx};
end
if numel(unique(canonicalFields)) ~= numel(canonicalFields)
    error(errID, ...
        ['Operation aborted. "labels" contains duplicate dimensions after ' ...
         'case normalization.']);
end
for iField = 1:numel(rawFields)
    labels.(canonicalFields{iField}) = ...
        iNormalizeLabelVector(labelsIn.(rawFields{iField}), errID, canonicalFields{iField});
end
end

function labels = iNormalizeLabelVector(labelsIn, errID, dimName)
if isstring(labelsIn) && isvector(labelsIn)
    labels = cellstr(labelsIn(:).');
elseif iscell(labelsIn) && isvector(labelsIn) && ...
        all(cellfun(@(c) ischar(c) || (isstring(c) && isscalar(c)), labelsIn))
    labels = cellstr(string(labelsIn(:).'));
else
    error(errID, ...
        ['Operation aborted. labels.%s must be a cell array of character ' ...
         'vectors or a string vector.'], ...
        dimName);
end
end

function iValidateEntryMeta(metaIn, schema, errID, entryName)
if ~isstruct(metaIn) || ~isscalar(metaIn)
    error(errID, ...
        'Operation aborted. Entry "%s.meta" must be a scalar struct.', entryName);
end

metaFields = fieldnames(metaIn);
reservedFields = [{'value','dimNames','meta'}, schema.requiredEntryFields, schema.optionalEntryFields];
colliding = intersect(metaFields, reservedFields);
if ~isempty(colliding)
    error(errID, ...
        ['Operation aborted. Entry "%s.meta" contains reserved field name(s): %s.'], ...
        entryName, strjoin(colliding, ', '));
end

if isfield(metaIn, 'FrameRateHz')
    val = metaIn.FrameRateHz;
    if ~isnumeric(val) || ~isscalar(val) || ~isfinite(val) || val <= 0
        error(errID, ...
            ['Operation aborted. Entry "%s.meta.FrameRateHz" must be a ' ...
             'positive numeric scalar.'], ...
            entryName);
    end
end
end

function iValidateEventInfoStruct(eventInfoIn, eLen, schema, errID)
if ~isstruct(eventInfoIn) || ~isscalar(eventInfoIn)
    error(errID, ...
        'Operation aborted. "eventInfo" must be a scalar struct.');
end

rawFields = fieldnames(eventInfoIn);
reqFields = schema.requiredEventInfoFields;
allowedFields = [schema.requiredEventInfoFields, schema.optionalEventInfoFields];

if ~all(ismember(reqFields, rawFields))
    error(errID, ...
        ['Operation aborted. "eventInfo" is missing one or more required ' ...
         'fields: %s.'], ...
        strjoin(reqFields, ', '));
end

unknownFields = setdiff(rawFields, allowedFields);
if ~isempty(unknownFields)
    error(errID, ...
        'Operation aborted. Unsupported field(s) in "eventInfo": %s.', ...
        strjoin(unknownFields, ', '));
end

eventID = eventInfoIn.eventID;
repIdx = eventInfoIn.repetitionIndex;
evName = eventInfoIn.eventName;
axisMode = eventInfoIn.eventAxisMode;

if ~isnumeric(eventID) || ~isvector(eventID) || isempty(eventID) || ...
        any(~isfinite(eventID(:))) || any(eventID(:) < 1) || ...
        any(mod(eventID(:),1) ~= 0)
    error(errID, ...
        ['Operation aborted. eventInfo.eventID must be a non-empty vector ' ...
         'of positive integers.']);
end
eventID = eventID(:);

if ~isnumeric(repIdx) || ~isvector(repIdx) || isempty(repIdx) || ...
        any(~isfinite(repIdx(:))) || any(repIdx(:) < 0) || ...
        any(mod(repIdx(:),1) ~= 0)
    error(errID, ...
        ['Operation aborted. eventInfo.repetitionIndex must be a non-empty ' ...
         'vector of non-negative integers.']);
end
repIdx = repIdx(:);

if isstring(evName) && isvector(evName)
    evName = cellstr(evName(:));
elseif iscell(evName) && isvector(evName) && ...
        all(cellfun(@(c) ischar(c) || (isstring(c) && isscalar(c)), evName))
    evName = cellstr(string(evName(:)));
else
    error(errID, ...
        ['Operation aborted. eventInfo.eventName must be a string vector ' ...
         'or a cell array of character vectors.']);
end

if ~(ischar(axisMode) || (isstring(axisMode) && isscalar(axisMode)))
    error(errID, ...
        'Operation aborted. eventInfo.eventAxisMode must be a text scalar.');
end
axisMode = lower(char(string(axisMode)));
if ~ismember(axisMode, schema.allowedEventAxisModes)
    error(errID, ...
        'Operation aborted. Invalid eventAxisMode "%s".', axisMode);
end

nE = numel(eventID);
if numel(repIdx) ~= nE || numel(evName) ~= nE
    error(errID, ...
        ['Operation aborted. eventInfo fields eventID, repetitionIndex, ' ...
         'and eventName must all have the same number of elements.']);
end
if nE ~= eLen
    error(errID, ...
        ['Operation aborted. Shared "eventInfo" length (%d) does not match ' ...
         'the size of the "E" dimension (%d).'], ...
        nE, eLen);
end

switch axisMode
    case 'instances'
        if any(repIdx(:) < 1)
            error(errID, ...
                ['Operation aborted. eventInfo.repetitionIndex must contain ' ...
                 'positive integers when eventAxisMode="instances".']);
        end
    case 'aggregated_repetitions'
        if any(repIdx(:) ~= 0)
            error(errID, ...
                ['Operation aborted. eventInfo.repetitionIndex must be all ' ...
                 'zeros when eventAxisMode="aggregated_repetitions".']);
        end
end

if isfield(eventInfoIn, 'baselinePeriod')
    bp = eventInfoIn.baselinePeriod;
    if ~isnumeric(bp) || ~isscalar(bp) || ~isfinite(bp) || bp <= 0
        error(errID, ...
            ['Operation aborted. eventInfo.baselinePeriod must be a ' ...
             'positive numeric scalar when provided.']);
    end
end
end
