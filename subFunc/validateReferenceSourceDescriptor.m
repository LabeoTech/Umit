function validateReferenceSourceDescriptor(descriptor, label)
%VALIDATEREFERENCESOURCEDESCRIPTOR Check the shape of a reference source descriptor.
%
%   validateReferenceSourceDescriptor(descriptor)
%   validateReferenceSourceDescriptor(descriptor, label)
%
%   Checks the field set and types only. The values are human aids and are
%   never resolved, so their content is not checked.
%
%   Error: validateDataParams:InvalidReferenceSource.
%
%   See also referenceSourceDescriptor.

if nargin < 2
    label = 'reference source';
end
errID = 'validateDataParams:InvalidReferenceSource';

if ~isstruct(descriptor) || ~isscalar(descriptor)
    error(errID, 'The %s must be a scalar struct.', label);
end
expected = fieldnames(referenceSourceDescriptor('legacy', ''));
missing = setdiff(expected, fieldnames(descriptor));
extra = setdiff(fieldnames(descriptor), expected);
if ~isempty(missing) || ~isempty(extra)
    error(errID, 'The %s has missing field(s) {%s} or unexpected field(s) {%s}.', ...
        label, strjoin(missing, ', '), strjoin(extra, ', '));
end
if ~ismember(char(string(descriptor.kind)), {'project', 'umitDat', 'file', 'legacy'})
    error(errID, 'The %s kind "%s" is unknown.', label, char(string(descriptor.kind)));
end
for iField = 1:numel(expected)
    name = expected{iField};
    value = descriptor.(name);
    switch name
        case 'extractionRuleVersion'
            ok = isempty(value) || (isnumeric(value) && isscalar(value));
        case 'sourceAligned'
            ok = islogical(value) && isscalar(value);
        otherwise
            ok = ischar(value) || (isstring(value) && isscalar(value)) || isempty(value);
    end
    if ~ok
        error(errID, 'The %s field "%s" has an invalid type.', label, name);
    end
end
end
