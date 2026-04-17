function schema = getUMTSchema(version)
%GETUMTSCHEMA Return centralized UMT schema definition.
%
%   schema = getUMTSchema()
%   schema = getUMTSchema(version)
%
%   Returns the toolbox schema definition used by genUMTStruct,
%   validateUMTStruct, and appendUMTEventInfo.
%
%   Inputs:
%       version - Numeric scalar schema version. Default: 1
%
%   Output:
%       schema  - Structure containing:
%                 schema.version
%                 schema.allowedKinds
%                 schema.allowedDims
%                 schema.requiredTopFields
%                 schema.optionalTopFields
%                 schema.requiredEntryFields
%                 schema.optionalEntryFields
%                 schema.requiredEventInfoFields
%                 schema.allowedEventAxisModes
%                 schema.allowedPatterns
%
%   Notes:
%       - This function centralizes all allowed kinds, dimension names, and
%         dimension layouts for .umt files.
%       - Store only "version" inside the .umt file. The actual schema
%         definition is maintained here in code.
%       - Shared event metadata is stored at the top level in "eventInfo"
%         and is required only for fully validated UMT structures when at
%         least one entry uses the "E" dimension.

if nargin < 1 || isempty(version)
    version = 1;
end

if ~isnumeric(version) || ~isscalar(version)
    error('Umitoolbox:getUMTSchema:invalidVersion', ...
        'Schema version must be a numeric scalar.');
end

switch version
    case 1
        schema = struct();

        schema.version = 1;

        schema.allowedKinds = {'image', 'roi'};
        schema.allowedDims  = {'Y', 'X', 'T', 'E', 'ROI', 'Measure'};

        schema.requiredTopFields = {'version', 'kind', 'data'};
        schema.optionalTopFields = {'labels', 'eventInfo'};

        schema.requiredEntryFields = {'value', 'dimNames'};
        schema.optionalEntryFields = {};

        schema.requiredEventInfoFields = { ...
            'eventID', ...
            'repetitionIndex', ...
            'eventName', ...
            'eventAxisMode'};

        schema.allowedEventAxisModes = {'instances', 'aggregated_repetitions'};

        schema.allowedPatterns = struct();
        schema.allowedPatterns.image = { ...
            {'Y','X'}, ...
            {'Y','X','T'}, ...
            {'Y','X','E'}, ...
            {'Y','X','T','E'}};

        schema.allowedPatterns.roi = { ...
            {'ROI'}, ...
            {'ROI','T'}, ...
            {'ROI','E'}, ...
            {'ROI','Measure'}, ...
            {'ROI','T','E'}, ...
            {'ROI','ROI'}};

    otherwise
        error('Umitoolbox:getUMTSchema:unsupportedVersion', ...
            'Unsupported UMT schema version: %g.', version);
end
end