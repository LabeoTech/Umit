function out = genUMTStruct(data, varargin)
%GENUMTSTRUCT Validate and package derived data into a single structure.
%
%   out = genUMTStruct(data)
%   out = genUMTStruct(data, obsID)
%   out = genUMTStruct(data, obsID, 'hasEvents', true)
%   out = genUMTStruct(data, obsID, 'extraInfo', S)
%
%   Inputs:
%       data    - Numeric array, struct array, or cell array with one
%                 element per observation.
%
%                 Accepted forms:
%                 1) numeric array
%                    A single observation containing one measure named "data".
%                 2) struct array
%                    One element per observation, with one or more measure
%                    fields per observation.
%                 3) cell array
%                    One element per observation. Each cell becomes a struct
%                    element with a single measure field named "data".
%
%   Optional positional input:
%       obsID   - 1D cell array of character vectors or strings describing
%                 each observation. Default: {'genericObs'}
%
%   Name-Value options:
%       hasEvents  - Logical scalar. True if data are split by events.
%                    Default: false
%       extraInfo  - Struct containing additional metadata to merge into
%                    the output. Default: struct()
%
%   Output:
%       out     - Structure containing:
%                 out.obsID
%                 out.data
%                 out.b_hasEvents
%                 out.b_hasMultipleMeasures
%                 out.dataCategory
%                 plus any non-reserved fields from extraInfo
%
%   Notes:
%       Data categories are assigned per measure as:
%           'scalar'
%           'time-vector'
%           'map'
%           'image-time-series'
%           'unknown'
%
%       Reserved field names in extraInfo are ignored with a warning:
%           'obsID'
%           'data'
%           'b_hasEvents'
%           'b_hasMultipleMeasures'
%           'dataCategory'

% -------------------------------------------------------------------------
% Parse inputs
% -------------------------------------------------------------------------
p = inputParser;
p.FunctionName = 'genUMTStruct';

addRequired(p, 'data', @(x) ...
    (isnumeric(x) && ~isempty(x)) || ...
    (isstruct(x) && ~isempty(x)) || ...
    (iscell(x) && ~isempty(x)));
addOptional(p, 'obsID', {'genericObs'}, @iValidateObsID);
addParameter(p, 'hasEvents', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'extraInfo', struct(), @(x) isstruct(x) && isscalar(x));

parse(p, data, varargin{:});

data = p.Results.data;
obsID = iNormalizeObsID(p.Results.obsID);
b_hasEvents = p.Results.hasEvents;
extraInfo = p.Results.extraInfo;

% -------------------------------------------------------------------------
% Remove reserved fields from extraInfo
% -------------------------------------------------------------------------
reservedFields = {'obsID', 'data', 'b_hasEvents', 'b_hasMultipleMeasures', 'dataCategory'};
extraFields = fieldnames(extraInfo);
conflictFields = extraFields(ismember(extraFields, reservedFields));

if ~isempty(conflictFields)
    warning('Umitoolbox:genUMTStruct:reservedExtraInfoFields', ...
        ['Ignoring reserved field(s) in extraInfo: %s.'], strjoin(conflictFields, ', '));
    extraInfo = rmfield(extraInfo, conflictFields);
end

% -------------------------------------------------------------------------
% Normalize input data into a struct array
% -------------------------------------------------------------------------
if isnumeric(data)
    data = struct('data', data);

elseif iscell(data)
    tmp = repmat(struct('data', []), numel(data), 1);
    for ii = 1:numel(data)
        tmp(ii).data = data{ii};
    end
    data = tmp;
end

% At this point, data must be a struct array.
if ~isstruct(data) || isempty(data)
    error('Umitoolbox:genUMTStruct:invalidInput', ...
        '"data" must resolve to a non-empty struct array.');
end

% -------------------------------------------------------------------------
% Validate observation count
% -------------------------------------------------------------------------
errID = 'Umitoolbox:genUMTStruct:IncompatibleSize';
errMsg = 'The length of data is different from the number of observations.';
assert(numel(data) == numel(obsID), errID, errMsg);

% -------------------------------------------------------------------------
% Validate struct-array field consistency
% -------------------------------------------------------------------------
fn = fieldnames(data);
if isempty(fn)
    error('Umitoolbox:genUMTStruct:invalidInput', ...
        'Input struct data must contain at least one field.');
end

for ii = 2:numel(data)
    if ~isequal(fieldnames(data(ii)), fn)
        error('Umitoolbox:genUMTStruct:invalidInput', ...
            'All elements of the input struct array must have the same fields.');
    end
end

% -------------------------------------------------------------------------
% Categorize each measure and validate consistency across observations
% -------------------------------------------------------------------------
dataCategory = struct();

for ii = 1:numel(fn)
    refCat = iClassifyMeasure(data(1).(fn{ii}), b_hasEvents);

    for jj = 2:numel(data)
        thisCat = iClassifyMeasure(data(jj).(fn{ii}), b_hasEvents);
        if ~strcmp(thisCat, refCat)
            error('Umitoolbox:genUMTStruct:inconsistentMeasureType', ...
                'Measure "%s" does not have consistent dimensionality across observations.', fn{ii});
        end
    end

    dataCategory.(fn{ii}) = refCat;
end

% -------------------------------------------------------------------------
% Assemble output
% -------------------------------------------------------------------------
out = extraInfo;
out.obsID = obsID;
out.data = data;
out.b_hasEvents = b_hasEvents;
out.b_hasMultipleMeasures = numel(fn) > 1;
out.dataCategory = dataCategory;

% Final format validation for .umt compatibility
validateUMTStruct(out);
end

% =========================================================================
% Local helpers
% =========================================================================

function tf = iValidateObsID(x)
%IVALIDATEOBSID Validate obsID cell array.
tf = iscell(x) && isvector(x) && ...
    all(cellfun(@(c) ischar(c) || (isstring(c) && isscalar(c)), x));
end

function obsID = iNormalizeObsID(obsID)
%INORMALIZEOBSID Convert obsID entries to char row cell array.
obsID = cellfun(@char, cellstr(string(obsID)), 'UniformOutput', false);
obsID = reshape(obsID, 1, []);
end

function dataCat = iClassifyMeasure(x, b_hasEvents)
%ICLASSIFYMEASURE Classify one measure by dimensionality.

if b_hasEvents
    x = iRemoveEventDimensionForClassification(x);
end

if isscalar(x)
    dataCat = 'scalar';
elseif isvector(x)
    dataCat = 'time-vector';
elseif ndims(x) == 2
    dataCat = 'map';
elseif ndims(x) == 3
    dataCat = 'image-time-series';
else
    dataCat = 'unknown';
end
end

function x = iRemoveEventDimensionForClassification(x)
%IREMOVEEVENTDIMENSIONFORCLASSIFICATION Remove leading event dimension.

if isempty(x)
    return
end

nDims = ndims(x);
idx = repmat({':'}, 1, nDims);
idx{1} = 1;
x = squeeze(x(idx{:}));
end