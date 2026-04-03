function DataParams = updateDataParam(folderPath, fieldPath, fieldValue, varargin)
%UPDATEDATAPARAM Update one field in folder-global DataParams.mat.
%
%   DataParams = updateDataParam(folderPath, fieldPath, fieldValue)
%   DataParams = updateDataParam(folderPath, 'view.origin_xy_px', [120 85])
%   DataParams = updateDataParam(folderPath, 'mask.logical', BW)
%
%   Loads DataParams.mat from the specified folder, updates the requested
%   nested field, validates the result, and saves the updated structure
%   back to disk.
%
%   Inputs:
%       folderPath - Path to the folder containing DataParams.mat
%       fieldPath  - Dotted field path, e.g. 'view.imageSizeYX'
%       fieldValue - Value to assign to the specified field
%
%   Name-Value options:
%       validateAfterSet - Logical scalar. If true, validate after update.
%                          Default: true
%
%   Output:
%       DataParams - Updated DataParams structure.
%
%   Notes:
%       - Intermediate fields must already exist.
%       - This function enforces folder-centric updates for the unique
%         DataParams.mat file in the folder.

p = inputParser;
p.FunctionName = 'updateDataParam';

addRequired(p, 'folderPath', @(x) ischar(x) || isstring(x));
addRequired(p, 'fieldPath', @(x) ischar(x) || isstring(x));
addRequired(p, 'fieldValue');
addParameter(p, 'validateAfterSet', true, @(x) islogical(x) && isscalar(x));

parse(p, folderPath, fieldPath, fieldValue, varargin{:});
R = p.Results;

fieldPath = char(string(fieldPath));
parts = strsplit(fieldPath, '.');

if isempty(fieldPath) || any(cellfun(@isempty, parts))
    error('updateDataParam:InvalidFieldPath', ...
        'fieldPath must be a non-empty dotted field path.');
end

DataParams = loadDataParams(folderPath);
DataParams = iSetNestedField(DataParams, parts, fieldValue);

if strcmp(fieldPath, 'mask.logical')
    if ~isempty(fieldValue) && isempty(DataParams.mask.createdOn)
        DataParams.mask.createdOn = datetime('now');
    end
end

if strcmp(fieldPath, 'registration.tform')
    if ~isempty(fieldValue)
        DataParams.registration.isRegistered = true;
        if isempty(DataParams.registration.createdOn)
            DataParams.registration.createdOn = datetime('now');
        end
    else
        DataParams.registration.isRegistered = false;
    end
end

DataParams.lastModified = datetime('now');

if R.validateAfterSet
    validateDataParams(DataParams);
end

filePath = fullfile(char(string(folderPath)), 'DataParams.mat');
save(filePath, 'DataParams', '-mat');
end

function S = iSetNestedField(S, parts, value)
%ISETNESTEDFIELD Recursively assign value to an existing nested struct field.

thisField = parts{1};

if ~isfield(S, thisField)
    error('updateDataParam:UnknownField', ...
        'Unknown field: %s', strjoin(parts, '.'));
end

if numel(parts) == 1
    S.(thisField) = value;
    return
end

if ~isstruct(S.(thisField)) || ~isscalar(S.(thisField))
    error('updateDataParam:NonStructIntermediate', ...
        'Intermediate field "%s" is not a scalar struct.', thisField);
end

S.(thisField) = iSetNestedField(S.(thisField), parts(2:end), value);
end