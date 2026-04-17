function out = genUMTStruct(data, varargin)
%GENUMTSTRUCT Create or append one measurement entry to a .umt structure.
%
%   Create a new .umt structure:
%       out = genUMTStruct(value, 'kind', kind, 'dimNames', dimNames)
%       out = genUMTStruct(value, 'kind', kind, 'entryName', name, ...
%                          'dimNames', dimNames, 'labels', labels, ...
%                          'meta', meta, 'eventInfo', eventInfo, ...
%                          'SaveFolder', SaveFolder)
%
%   Append a new entry to an existing .umt structure:
%       out = genUMTStruct(umt, 'value', value, 'entryName', name, ...
%                          'dimNames', dimNames)
%       out = genUMTStruct(umt, 'value', value, 'entryName', name, ...
%                          'dimNames', dimNames, 'labels', labels, ...
%                          'meta', meta, 'eventInfo', eventInfo, ...
%                          'SaveFolder', SaveFolder, 'overwrite', true)
%
%   Notes:
%       - labels remain shared top-level display/reference metadata.
%       - entry-specific processing metadata is stored in entry.meta.
%       - shared event metadata is stored only at the top level in
%         out.eventInfo.
%       - When entry.meta.FrameRateHz is not provided and SaveFolder is
%         given, genUMTStruct attempts to auto-populate FrameRateHz from
%         AcqInfos.mat via loadMetaData, with a direct AcqInfos.mat
%         fallback if needed.

errID = 'Umitoolbox:genUMTStruct:invalidInput';

if iLooksLikeUMT(data)
    out = iAppendToExistingUMT(data, varargin{:});
else
    out = iCreateNewUMT(data, varargin{:});
end

validateUMTStruct(out, 'requireEventInfo', false);

end

function out = iCreateNewUMT(value, varargin)
errID = 'Umitoolbox:genUMTStruct:invalidInput';
schema = getUMTSchema(1);

if ~(isnumeric(value) || islogical(value))
    error(errID, ...
        'Operation aborted. New entry value must be numeric or logical.');
end

p = inputParser;
p.FunctionName = 'genUMTStruct';
addParameter(p, 'kind', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'entryName', 'main', @iValidateEntryNameInput);
addParameter(p, 'dimNames', {}, @iValidateDimNamesInput);
addParameter(p, 'labels', struct(), @(x) isstruct(x) && isscalar(x));
addParameter(p, 'meta', struct(), @(x) isstruct(x) && isscalar(x));
addParameter(p, 'eventInfo', struct(), @(x) isstruct(x) && isscalar(x));
addParameter(p, 'SaveFolder', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
parse(p, varargin{:});

if ~iHasNameValue(varargin, 'kind')
    error(errID, ...
        'Operation aborted. "kind" must be provided when creating a new UMT structure.');
end
if ~iHasNameValue(varargin, 'dimNames')
    error(errID, ...
        'Operation aborted. "dimNames" must be provided when creating a new UMT structure.');
end

kind = iNormalizeKind(p.Results.kind, schema, errID);
entryName = iNormalizeEntryName(p.Results.entryName, errID);
dimNames = iNormalizeDimNames(p.Results.dimNames, schema, errID, entryName);
labels = iNormalizeLabelsStruct(p.Results.labels, schema, errID);
meta = iNormalizeEntryMetaStruct(p.Results.meta, schema, errID, entryName);

% Auto-populate critical per-entry processing metadata from AcqInfos.mat
% when the caller supplied SaveFolder but did not explicitly provide the
% metadata in entry.meta.
meta = iAutofillEntryMeta(meta, p.Results.SaveFolder);
meta = iNormalizeEntryMetaStruct(meta, schema, errID, entryName);

value = iForceDeclaredShape(value, dimNames, errID, entryName);

out = struct();
out.version = schema.version;
out.kind = kind;
out.data = struct();
out.data.(entryName) = struct('value', value, 'dimNames', {dimNames});
if ~isempty(fieldnames(meta))
    out.data.(entryName).meta = meta;
end

if ~isempty(fieldnames(labels))
    out.labels = labels;
end
if iHasNameValue(varargin, 'eventInfo') && ~isempty(fieldnames(p.Results.eventInfo))
    out.eventInfo = iNormalizeEventInfoStruct(p.Results.eventInfo, schema, errID);
end
end

function out = iAppendToExistingUMT(umt, varargin)
errID = 'Umitoolbox:genUMTStruct:invalidInput';
validateUMTStruct(umt, 'requireEventInfo', false);
out = umt;
schema = getUMTSchema(out.version);

p = inputParser;
p.FunctionName = 'genUMTStruct';
addParameter(p, 'value', [], @(x) isnumeric(x) || islogical(x));
addParameter(p, 'entryName', '', @iValidateEntryNameInput);
addParameter(p, 'dimNames', {}, @iValidateDimNamesInput);
addParameter(p, 'labels', struct(), @(x) isstruct(x) && isscalar(x));
addParameter(p, 'meta', struct(), @(x) isstruct(x) && isscalar(x));
addParameter(p, 'eventInfo', struct(), @(x) isstruct(x) && isscalar(x));
addParameter(p, 'SaveFolder', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'overwrite', false, @(x) islogical(x) && isscalar(x));
parse(p, varargin{:});

if ~iHasNameValue(varargin, 'value')
    error(errID, ...
        'Operation aborted. "value" must be provided when appending to an existing UMT structure.');
end
if ~iHasNameValue(varargin, 'entryName')
    error(errID, ...
        'Operation aborted. "entryName" must be provided when appending to an existing UMT structure.');
end
if ~iHasNameValue(varargin, 'dimNames')
    error(errID, ...
        'Operation aborted. "dimNames" must be provided when appending to an existing UMT structure.');
end

value = p.Results.value;
entryName = iNormalizeEntryName(p.Results.entryName, errID);
dimNames = iNormalizeDimNames(p.Results.dimNames, schema, errID, entryName);
labels = iNormalizeLabelsStruct(p.Results.labels, schema, errID);
meta = iNormalizeEntryMetaStruct(p.Results.meta, schema, errID, entryName);
overwrite = p.Results.overwrite;

% Auto-populate critical per-entry processing metadata from AcqInfos.mat
% when the caller supplied SaveFolder but did not explicitly provide the
% metadata in entry.meta.
meta = iAutofillEntryMeta(meta, p.Results.SaveFolder);
meta = iNormalizeEntryMetaStruct(meta, schema, errID, entryName);

value = iForceDeclaredShape(value, dimNames, errID, entryName);

if isfield(out.data, entryName) && ~overwrite
    error(errID, ...
        ['Operation aborted. Entry "%s" already exists. Use overwrite=true ' ...
         'to replace it.'], ...
        entryName);
end

out.data.(entryName) = struct('value', value, 'dimNames', {dimNames});
if ~isempty(fieldnames(meta))
    out.data.(entryName).meta = meta;
end

out = iMergeTopLevelLabels(out, labels, overwrite, errID);

if iHasNameValue(varargin, 'eventInfo') && ~isempty(fieldnames(p.Results.eventInfo))
    newEventInfo = iNormalizeEventInfoStruct(p.Results.eventInfo, schema, errID);
    out = iMergeTopLevelEventInfo(out, newEventInfo, overwrite, errID);
end
end

function valueOut = iForceDeclaredShape(valueIn, dimNames, errID, entryName)
nDimsExpected = numel(dimNames);
if nDimsExpected == 0
    if ~isscalar(valueIn)
        error(errID, ...
            'Operation aborted. Scalar entry "%s" must use dimNames = {}.', ...
            entryName);
    end
    valueOut = valueIn;
    return
end
if isscalar(valueIn)
    error(errID, ...
        'Operation aborted. Scalar entry "%s" must use dimNames = {}.', ...
        entryName);
end
sz = size(valueIn);
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
valueOut = reshape(valueIn, sz);
end

function out = iMergeTopLevelLabels(out, newLabels, overwrite, errID)
newFields = fieldnames(newLabels);
if isempty(newFields)
    return
end
if ~isfield(out, 'labels') || isempty(fieldnames(out.labels))
    out.labels = newLabels;
    return
end
for iField = 1:numel(newFields)
    thisField = newFields{iField};
    if isfield(out.labels, thisField)
        if ~isequal(out.labels.(thisField), newLabels.(thisField))
            if overwrite
                out.labels.(thisField) = newLabels.(thisField);
            else
                error(errID, ...
                    ['Operation aborted. labels.%s already exists with a ' ...
                     'different value. Use overwrite=true to replace it.'], ...
                    thisField);
            end
        end
    else
        out.labels.(thisField) = newLabels.(thisField);
    end
end
end

function out = iMergeTopLevelEventInfo(out, newEventInfo, overwrite, errID)
if isempty(fieldnames(newEventInfo))
    return
end
if ~isfield(out, 'eventInfo') || isempty(fieldnames(out.eventInfo))
    out.eventInfo = newEventInfo;
    return
end
if isequal(out.eventInfo, newEventInfo)
    return
end
if overwrite
    out.eventInfo = newEventInfo;
else
    error(errID, ...
        ['Operation aborted. "eventInfo" already exists with a different ' ...
         'value. Use overwrite=true to replace it.']);
end
end

function tf = iLooksLikeUMT(x)
tf = isstruct(x) && isscalar(x) && ...
    all(ismember({'version', 'kind', 'data'}, fieldnames(x)));
end

function tf = iValidateEntryNameInput(x)
tf = ischar(x) || (isstring(x) && isscalar(x));
end

function entryName = iNormalizeEntryName(entryNameIn, errID)
entryName = char(string(entryNameIn));
if isempty(entryName) || ~isvarname(entryName)
    error(errID, ...
        ['Operation aborted. "entryName" must be a valid MATLAB field ' ...
         'name.']);
end
end

function tf = iValidateDimNamesInput(x)
tf = isempty(x) || isstring(x) || iscell(x);
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
            'Operation aborted. "%s.dimNames" must be a 1D text array.', ...
            entryName);
    end
    rawDims = cellstr(dimNamesIn(:).');
elseif iscell(dimNamesIn) && isvector(dimNamesIn) && ...
        all(cellfun(@(c) ischar(c) || (isstring(c) && isscalar(c)), dimNamesIn))
    rawDims = cellstr(string(dimNamesIn(:).'));
else
    error(errID, ...
        ['Operation aborted. "%s.dimNames" must be a cell array of ' ...
         'character vectors or a string vector.'], ...
        entryName);
end

dimNames = cell(1, numel(rawDims));
for iDim = 1:numel(rawDims)
    idx = find(strcmpi(rawDims{iDim}, allowedDims), 1, 'first');
    if isempty(idx)
        error(errID, ...
            'Operation aborted. Invalid dimension name "%s".', rawDims{iDim});
    end
    dimNames{iDim} = allowedDims{idx};
end
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
elif iscell(labelsIn) && isvector(labelsIn) && ...
        all(cellfun(@(c) ischar(c) || (isstring(c) && isscalar(c)), labelsIn))
    labels = cellstr(string(labelsIn(:).'));
else
    error(errID, ...
        ['Operation aborted. labels.%s must be a cell array of character ' ...
         'vectors or a string vector.'], ...
        dimName);
end
end

function meta = iNormalizeEntryMetaStruct(metaIn, schema, errID, entryName)
if ~isstruct(metaIn) || ~isscalar(metaIn)
    error(errID, ...
        'Operation aborted. "%s.meta" must be a scalar struct.', entryName);
end
meta = metaIn;
if isempty(fieldnames(meta))
    return
end
reservedFields = [{'value','dimNames','meta'}, schema.requiredEntryFields, schema.optionalEntryFields];
colliding = intersect(fieldnames(meta), reservedFields);
if ~isempty(colliding)
    error(errID, ...
        ['Operation aborted. "%s.meta" contains reserved field name(s): %s.'], ...
        entryName, strjoin(colliding, ', '));
end
if isfield(meta, 'FrameRateHz')
    val = meta.FrameRateHz;
    if ~isnumeric(val) || ~isscalar(val) || ~isfinite(val) || val <= 0
        error(errID, ...
            ['Operation aborted. "%s.meta.FrameRateHz" must be a positive ' ...
             'numeric scalar.'], ...
            entryName);
    end
end
end

function metaOut = iAutofillEntryMeta(metaIn, saveFolderIn)
%IAUTOFILLENTRYMETA Populate critical per-entry metadata from SaveFolder.
%
% FrameRateHz is first resolved through loadMetaData for consistency with
% existing toolbox loaders. If that does not return a usable value, the
% function falls back to reading AcqInfos.mat directly.

metaOut = metaIn;
saveFolder = char(string(saveFolderIn));

if isempty(saveFolder) || ~isfolder(saveFolder)
    return
end

if isfield(metaOut, 'FrameRateHz') && ~isempty(metaOut.FrameRateHz)
    return
end

frameRateHz = [];

% First try the toolbox compatibility loader.
try
    md = loadMetaData(saveFolder);
    if isstruct(md)
        if isfield(md, 'FrameRateHz') && ~isempty(md.FrameRateHz)
            frameRateHz = md.FrameRateHz;
        elseif isfield(md, 'Freq') && ~isempty(md.Freq)
            frameRateHz = md.Freq;
        end
    end
catch
end

% Fall back to reading AcqInfos.mat directly when needed.
if isempty(frameRateHz)
    acqFile = fullfile(saveFolder, 'AcqInfos.mat');
    if isfile(acqFile)
        try
            s = load(acqFile);
            fn = fieldnames(s);
            for iField = 1:numel(fn)
                thisVal = s.(fn{iField});
                if ~isstruct(thisVal)
                    continue
                end
                if isfield(thisVal, 'FrameRateHz') && ~isempty(thisVal.FrameRateHz)
                    frameRateHz = thisVal.FrameRateHz;
                    break
                elseif isfield(thisVal, 'Freq') && ~isempty(thisVal.Freq)
                    frameRateHz = thisVal.Freq;
                    break
                end
            end
        catch
        end
    end
end

if ~isempty(frameRateHz)
    frameRateHz = double(frameRateHz);
    if isscalar(frameRateHz) && isfinite(frameRateHz) && frameRateHz > 0
        metaOut.FrameRateHz = frameRateHz;
    end
end
end

function eventInfo = iNormalizeEventInfoStruct(eventInfoIn, schema, errID)
if ~isstruct(eventInfoIn) || ~isscalar(eventInfoIn)
    error(errID, ...
        'Operation aborted. "eventInfo" must be a scalar struct.');
end

eventInfo = struct();
rawFields = fieldnames(eventInfoIn);
if isempty(rawFields)
    return
end

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
repIdx  = eventInfoIn.repetitionIndex;
evName  = eventInfoIn.eventName;
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
elif iscell(evName) && isvector(evName) && ...
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

eventInfo.eventID = eventID;
eventInfo.repetitionIndex = repIdx;
eventInfo.eventName = evName;
eventInfo.eventAxisMode = axisMode;

if isfield(eventInfoIn, 'baselinePeriod')
    bp = eventInfoIn.baselinePeriod;
    if ~isnumeric(bp) || ~isscalar(bp) || ~isfinite(bp) || bp <= 0
        error(errID, ...
            ['Operation aborted. eventInfo.baselinePeriod must be a ' ...
             'positive numeric scalar when provided.']);
    end
    eventInfo.baselinePeriod = double(bp);
end
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
