function out = appendUMTEventInfo(umt, varargin)
%APPENDUMTEVENTINFO Append or replace shared event metadata in a .umt structure.
%
%   out = appendUMTEventInfo(umt, ...
%       'eventID', eventID, ...
%       'repetitionIndex', repetitionIndex, ...
%       'eventName', eventName, ...
%       'eventAxisMode', eventAxisMode)
%
%   out = appendUMTEventInfo(..., 'overwrite', true)
%
%   Inputs:
%       umt - UMT structure to update.
%
%   Name-Value options:
%       eventID         - Numeric vector of positive integer event IDs.
%       repetitionIndex - Numeric vector of non-negative integers.
%       eventName       - Text vector/cellstr with one name per E element.
%       eventAxisMode   - One of:
%                           'instances'
%                           'aggregated_repetitions'
%       overwrite       - Logical scalar. If true, existing eventInfo is
%                         replaced. Default: false
%
%   Output:
%       out             - Updated and fully validated UMT structure.
%
%   Notes:
%       - eventInfo is shared globally across all entries that use the
%         "E" dimension.
%       - All entries that use "E" must have the same E length.
%       - If no entry uses "E", this function errors.
%       - Trailing singleton dimensions declared in dimNames are allowed
%         even if MATLAB suppresses them in ndims(value).

errID = 'Umitoolbox:appendUMTEventInfo:invalidInput';

validateUMTStruct(umt, 'requireEventInfo', false);

schema = getUMTSchema(umt.version);

p = inputParser;
p.FunctionName = 'appendUMTEventInfo';

addParameter(p, 'eventID', [], @(x) isnumeric(x) && isvector(x));
addParameter(p, 'repetitionIndex', [], @(x) isnumeric(x) && isvector(x));
addParameter(p, 'eventName', [], @(x) isstring(x) || iscell(x));
addParameter(p, 'eventAxisMode', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'overwrite', false, @(x) islogical(x) && isscalar(x));

parse(p, varargin{:});

requiredNV = {'eventID','repetitionIndex','eventName','eventAxisMode'};
for iReq = 1:numel(requiredNV)
    if ~iHasNameValue(varargin, requiredNV{iReq})
        error(errID, ...
            'Operation aborted. "%s" must be provided.', requiredNV{iReq});
    end
end

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

eventInfo = iNormalizeEventInfoFromNameValue( ...
    p.Results.eventID, ...
    p.Results.repetitionIndex, ...
    p.Results.eventName, ...
    p.Results.eventAxisMode, ...
    eLengths(1), ...
    schema, ...
    errID);

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

% =========================================================================
% Local helpers
% =========================================================================

function tf = iHasNameValue(args, name)
%IHASNAMEVALUE Return true when a name-value list contains the given name.

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
%INORMALIZEDIMNAMES Normalize and validate dimNames.

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
%IGETDECLAREDDIMENSIONSIZES Return dimension sizes compatible with dimNames.
%
% This helper does not rely on ndims(value). It allows trailing singleton
% dimensions declared in dimNames even when MATLAB suppresses them in ndims.

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

function eventInfo = iNormalizeEventInfoFromNameValue(eventID, repIdx, evName, axisMode, eLen, schema, errID)
%INORMALIZEEVENTINFOFROMNAMEVALUE Normalize name-value event info inputs.

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