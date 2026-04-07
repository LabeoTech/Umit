function out = genUMTStruct(data, varargin)
%GENUMTSTRUCT Create or append one measurement entry to a .umt structure.
%
%   Create a new .umt structure:
%       out = genUMTStruct(value, 'kind', kind, 'dimNames', dimNames)
%       out = genUMTStruct(value, 'kind', kind, 'entryName', name, ...
%                          'dimNames', dimNames, 'labels', labels)
%
%   Append a new entry to an existing .umt structure:
%       out = genUMTStruct(umt, 'value', value, 'entryName', name, ...
%                          'dimNames', dimNames)
%       out = genUMTStruct(umt, 'value', value, 'entryName', name, ...
%                          'dimNames', dimNames, 'labels', labels, ...
%                          'overwrite', true)
%
%   Inputs:
%       data    - Either:
%                 1) Numeric/logical payload used to create a new .umt
%                    structure, or
%                 2) Existing .umt structure to which one entry will be
%                    appended (or overwritten if requested).
%
%   Name-Value options for create mode:
%       kind      - Required. Allowed values are defined by getUMTSchema.
%       entryName - Name of the created entry. Default: 'main'
%       dimNames  - Required. Allowed values are defined by getUMTSchema.
%       labels    - Optional shared top-level labels struct. Default: struct()
%
%   Name-Value options for append mode:
%       value     - Required. Numeric/logical payload to append.
%       entryName - Required. Name of the entry to append.
%       dimNames  - Required. Allowed values are defined by getUMTSchema.
%       labels    - Optional shared top-level labels struct to merge into
%                   the output. Default: struct()
%       overwrite - Logical scalar. If true, existing entry names and
%                   conflicting labels can be replaced. Default: false
%
%   Output:
%       out       - Validated .umt structure using the current schema.
%
%   Notes:
%       - This function appends one entry per call. Multiple measurements
%         can be added by calling genUMTStruct repeatedly on the same
%         output structure.
%       - Shared labels are stored only at the top level as out.labels.
%       - Final output is validated with validateUMTStruct before return.

errID = 'Umitoolbox:genUMTStruct:invalidInput';

if iLooksLikeUMT(data)
    out = iAppendToExistingUMT(data, varargin{:});
else
    out = iCreateNewUMT(data, varargin{:});
end

validateUMTStruct(out);

end

% =========================================================================
% Local helpers
% =========================================================================

function out = iCreateNewUMT(value, varargin)
%ICREATENEWUMT Create a new .umt structure with one entry.

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

out = struct();
out.version = schema.version;
out.kind = kind;
out.data = struct();
out.data.(entryName) = struct('value', value, 'dimNames', {dimNames});

if ~isempty(fieldnames(labels))
    out.labels = labels;
end

end

function out = iAppendToExistingUMT(umt, varargin)
%IAPPENDTOEXISTINGUMT Append one entry to an existing .umt structure.

errID = 'Umitoolbox:genUMTStruct:invalidInput';

validateUMTStruct(umt);
out = umt;
schema = getUMTSchema(out.version);

p = inputParser;
p.FunctionName = 'genUMTStruct';

addParameter(p, 'value', [], @(x) isnumeric(x) || islogical(x));
addParameter(p, 'entryName', '', @iValidateEntryNameInput);
addParameter(p, 'dimNames', {}, @iValidateDimNamesInput);
addParameter(p, 'labels', struct(), @(x) isstruct(x) && isscalar(x));
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
overwrite = p.Results.overwrite;

if isfield(out.data, entryName) && ~overwrite
    error(errID, ...
        ['Operation aborted. Entry "%s" already exists. Use overwrite=true ' ...
         'to replace it.'], ...
        entryName);
end

out.data.(entryName) = struct('value', value, 'dimNames', {dimNames});
out = iMergeTopLevelLabels(out, labels, overwrite, errID);

end

function out = iMergeTopLevelLabels(out, newLabels, overwrite, errID)
%IMERGETOPLEVELLABELS Merge shared labels into an output structure.

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

function tf = iLooksLikeUMT(x)
%ILOOKSLIKEUMT Return true when input resembles a .umt structure.

tf = isstruct(x) && isscalar(x) && ...
    all(ismember({'version', 'kind', 'data'}, fieldnames(x)));

end

function tf = iValidateEntryNameInput(x)
%IVALIDENTRYNAMEINPUT Lightweight validator for parser entryName input.

tf = ischar(x) || (isstring(x) && isscalar(x));

end

function entryName = iNormalizeEntryName(entryNameIn, errID)
%INORMALIZEENTRYNAME Normalize and validate one entry name.

entryName = char(string(entryNameIn));

if isempty(entryName) || ~isvarname(entryName)
    error(errID, ...
        ['Operation aborted. "entryName" must be a valid MATLAB field ' ...
         'name.']);
end

end

function tf = iValidateDimNamesInput(x)
%IVALIDATEDIMNAMESINPUT Lightweight validator for parser dimNames input.

tf = isempty(x) || isstring(x) || iscell(x);

end

function kind = iNormalizeKind(kindIn, schema, errID)
%INORMALIZEKIND Normalize and validate kind.

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