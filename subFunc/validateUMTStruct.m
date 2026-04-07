function validateUMTStruct(umt)
%VALIDATEUMTSTRUCT Validate UMIT processed-data structure for .umt files.
%
%   validateUMTStruct(umt)
%
%   Validates whether the input structure follows the current .umt schema.
%   The schema definition is obtained from getUMTSchema using umt.version.
%
%   Required top-level fields:
%       version
%       kind
%       data
%
%   Optional top-level fields:
%       labels
%
%   Required fields for each entry in "umt.data":
%       value
%       dimNames
%
%   Rules:
%       - umt.version must be numeric scalar and supported by getUMTSchema.
%       - umt.kind must be one of the allowed kinds from the schema.
%       - umt.data must be a non-empty scalar struct.
%       - Each entry value must be numeric or logical.
%       - Scalars must use dimNames = {}.
%       - One-dimensional data must be stored as column vectors.
%       - Labels, if provided, must be shared top-level labels and must
%         match the corresponding dimension sizes in all entries that use
%         that dimension.
%
%   Error:
%       Throws an error if validation fails.

errID = 'Umitoolbox:validateUMTStruct:invalidInput';

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

% Track dimension lengths by name for shared top-level label validation.
dimUsage = struct();

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
    dimSizes = iGetEffectiveDimensionSizes(value, errID, entryName);

    if isempty(dimSizes)
        if ~isempty(dimNames)
            error(errID, ...
                'Operation aborted. Scalar entry "%s" must use dimNames = {}.', ...
                entryName);
        end
    else
        if numel(dimNames) ~= numel(dimSizes)
            error(errID, ...
                ['Operation aborted. The number of dimNames in entry "%s" ' ...
                 'does not match the effective dimensionality of its value.'], ...
                entryName);
        end
    end

    if ~iIsAllowedPattern(kind, dimNames, schema)
        error(errID, ...
            ['Operation aborted. Entry "%s" has invalid dimNames for kind "%s".'], ...
            entryName, kind);
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

end

% =========================================================================
% Local helpers
% =========================================================================

function kind = iNormalizeKind(kindIn, schema, errID)
%INORMALIZEKIND Normalize and validate top-level kind.

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

function dimSizes = iGetEffectiveDimensionSizes(value, errID, entryName)
%IGETEFFECTIVEDIMENSIONSIZES Return effective dimension sizes for validation.
%
% Scalars return [].
% One-dimensional data must be stored along the first dimension.
% Higher-dimensional data return size(value).

if isscalar(value)
    dimSizes = [];
    return
end

sz = size(value);
nonSingletonDims = find(sz ~= 1);

if isempty(nonSingletonDims)
    dimSizes = [];
    return
end

if numel(nonSingletonDims) == 1
    if nonSingletonDims ~= 1
        error(errID, ...
            ['Operation aborted. Entry "%s" is one-dimensional but not stored ' ...
             'as a column vector.'], ...
            entryName);
    end
    dimSizes = sz(1);
    return
end

dimSizes = sz;

end

function tf = iIsAllowedPattern(kind, dimNames, schema)
%IISALLOWEDPATTERN Return true when dimNames is allowed for the given kind.

if isempty(dimNames)
    tf = true;
    return
end

allowedPatterns = schema.allowedPatterns.(kind);
tf = any(cellfun(@(c) isequal(dimNames, c), allowedPatterns));

end

function labels = iNormalizeLabelsStruct(labelsIn, schema, errID)
%INORMALIZELABELSSTRUCT Normalize and validate top-level labels struct.

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
%INORMALIZELABELVECTOR Normalize one label vector.

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