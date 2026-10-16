function validateRegistrationRecord(registration, view)
%VALIDATEREGISTRATIONRECORD Validate DataParams.registration (version 2).
%
%   validateRegistrationRecord(registration, view)
%
%   Shape, types and cross-field rules of the version-2 record (see
%   newRegistrationRecord). VIEW is DataParams.view. The reference image is
%   not re-hashed here (loading stays cheap); estimate and apply re-verify
%   the stored image against contentHash (verifyRegistrationReference).
%
%   Cross-field rules:
%       - reference.imageSizeYX equals the size of reference.image.
%       - Hash fields are '' or 'sha256:' + 64 lowercase hex digits.
%       - isReviewed requires an estimate.
%       - Not aligned: appliedOn, appliedBy, confirmationMode and
%         appliedReferenceHash are empty.
%       - Aligned: an estimate exists; when a reference exists,
%         appliedReferenceHash equals reference.contentHash and the
%         reference size equals DataParams.view.imageSizeYX.
%
%   Errors use validateDataParams:<condition> identifiers.
%
%   See also validateDataParams, newRegistrationRecord.

[empty, referenceTemplate, estimateTemplate] = newRegistrationRecord();

iRequireExactFields(registration, fieldnames(empty), 'DataParams.registration', ...
    'validateDataParams:InvalidRegistration');
if ~isequal(registration.version, 2)
    error('validateDataParams:InvalidRegistration', ...
        'DataParams.registration.version must be 2.');
end

% Reference
reference = registration.reference;
if ~isempty(reference)
    iRequireExactFields(reference, fieldnames(referenceTemplate), ...
        'DataParams.registration.reference', 'validateDataParams:InvalidRegistrationReference');
    image = reference.image;
    if ~(isnumeric(image) || islogical(image)) || isempty(image) || ~ismatrix(image)
        error('validateDataParams:InvalidReferenceImage', ...
            'DataParams.registration.reference.image must be a non-empty 2-D numeric array.');
    end
    if ~isequal(double(reference.imageSizeYX(:).'), double(size(image)))
        error('validateDataParams:ReferenceImageSizeMismatch', ...
            'DataParams.registration.reference.imageSizeYX must equal the image size.');
    end
    iRequireHash(reference.contentHash, 'reference.contentHash', false);
    iRequireText(reference, {'name', 'description', 'axisConvention', 'selectedOn'}, ...
        'DataParams.registration.reference');
    validateReferenceSourceDescriptor(reference.source, 'DataParams.registration.reference.source');
end

% Estimate
estimate = registration.estimate;
if ~isempty(estimate)
    iRequireExactFields(estimate, fieldnames(estimateTemplate), ...
        'DataParams.registration.estimate', 'validateDataParams:InvalidRegistrationEstimate');
    iRequireTform(estimate.tform);
    iRequireText(estimate, {'transformType', 'method', 'qcStatus', 'qcWarning', ...
        'qcFigureFile', 'qcPreviewImageFile', 'createdOn', 'movingSource'}, ...
        'DataParams.registration.estimate');
    iRequireHash(estimate.referenceHash, 'estimate.referenceHash', true);
    iValidateQcMetrics(estimate.qcMetrics, fieldnames(estimateTemplate.qcMetrics));
end

% Apply record
apply = registration.apply;
iRequireExactFields(apply, fieldnames(empty.apply), ...
    'DataParams.registration.apply', 'validateDataParams:InvalidRegistrationApply');
for name = {'isRegistered', 'isReviewed'}
    if ~islogical(apply.(name{1})) || ~isscalar(apply.(name{1}))
        error('validateDataParams:InvalidRegistrationApply', ...
            'DataParams.registration.apply.%s must be a scalar logical.', name{1});
    end
end
iRequireText(apply, {'appliedOn', 'appliedBy', 'confirmationMode', ...
    'appliedReferenceHash', 'notes'}, 'DataParams.registration.apply');
iRequireHash(apply.appliedReferenceHash, 'apply.appliedReferenceHash', true);
if ~isstruct(apply.manualAdjustment) || ~isscalar(apply.manualAdjustment)
    error('validateDataParams:InvalidRegistrationApply', ...
        'DataParams.registration.apply.manualAdjustment must be a scalar struct.');
end
if apply.isReviewed && isempty(estimate)
    error('validateDataParams:ReviewedWithoutTform', ...
        'Registration cannot be marked as reviewed without an estimate.');
end
if ~apply.isRegistered
    for name = {'appliedOn', 'appliedBy', 'confirmationMode', 'appliedReferenceHash'}
        if ~isempty(apply.(name{1}))
            error('validateDataParams:AppliedFieldWithoutRegistration', ...
                'DataParams.registration.apply.%s must be empty when isRegistered is false.', ...
                name{1});
        end
    end
else
    if isempty(estimate)
        error('validateDataParams:MissingTform', ...
            'An aligned folder must keep the estimate (transform) it applied.');
    end
    if ~isempty(reference)
        if ~strcmp(apply.appliedReferenceHash, reference.contentHash)
            error('validateDataParams:AppliedReferenceMismatch', ...
                ['DataParams.registration.apply.appliedReferenceHash must equal the ' ...
                 'reference contentHash of an aligned folder.']);
        end
        if isstruct(view) && isfield(view, 'imageSizeYX') && ~isempty(view.imageSizeYX) && ...
                ~isequal(double(view.imageSizeYX(:).'), double(reference.imageSizeYX(:).'))
            error('validateDataParams:ReferenceImageSizeMismatch', ...
                ['DataParams.registration.reference size must match ' ...
                 'DataParams.view.imageSizeYX once the folder is aligned.']);
        end
    end
end
end

function iRequireExactFields(s, expected, label, errID)
if ~isstruct(s) || ~isscalar(s)
    error(errID, '%s must be a scalar struct.', label);
end
missing = setdiff(expected, fieldnames(s));
extra = setdiff(fieldnames(s), expected);
if ~isempty(missing)
    error(errID, 'Missing field(s) in %s: %s.', label, strjoin(missing, ', '));
end
if ~isempty(extra)
    error(errID, 'Unexpected field(s) in %s: %s.', label, strjoin(extra, ', '));
end
end

function iRequireText(s, names, label)
for iName = 1:numel(names)
    value = s.(names{iName});
    if ~(isempty(value) || ischar(value) || (isstring(value) && isscalar(value)))
        error('validateDataParams:InvalidRegistrationTextField', ...
            '%s.%s must be empty, a char vector, or a string scalar.', label, names{iName});
    end
end
end

function iRequireHash(value, label, allowEmpty)
text = char(string(value));
if isempty(text) && allowEmpty
    return
end
if isempty(regexp(text, '^sha256:[0-9a-f]{64}$', 'once'))
    error('validateDataParams:InvalidContentHash', ...
        'DataParams.registration.%s must be ''sha256:'' followed by 64 lowercase hex digits.', label);
end
end

function iRequireTform(tf)
isValid = isa(tf, 'affine2d') || isa(tf, 'projective2d') || isa(tf, 'rigidtform2d') || ...
    isa(tf, 'simtform2d') || isa(tf, 'affinetform2d') || ...
    isa(tf, 'images.geotrans.GeometricTransformation2D');
if ~isValid
    error('validateDataParams:InvalidTform', ...
        'DataParams.registration.estimate.tform is not a recognized 2D geometric transform object.');
end
end

function iValidateQcMetrics(qc, required)
if ~isstruct(qc) || ~isscalar(qc)
    error('validateDataParams:InvalidQCMetrics', ...
        'DataParams.registration.estimate.qcMetrics must be a scalar struct.');
end
missing = setdiff(required, fieldnames(qc));
if ~isempty(missing)
    error('validateDataParams:MissingQCMetricsField', ...
        'Missing field: DataParams.registration.estimate.qcMetrics.%s', missing{1});
end
for name = {'MIBefore', 'MIAfter', 'MIDelta', 'rotationDeg', 'determinant'}
    value = qc.(name{1});
    if ~(isempty(value) || (isnumeric(value) && isscalar(value) && isfinite(value)))
        error('validateDataParams:InvalidQCMetricScalar', ...
            'DataParams.registration.estimate.qcMetrics.%s must be empty or a finite numeric scalar.', ...
            name{1});
    end
end
for name = {'translationXY_px', 'scaleXY'}
    value = qc.(name{1});
    if ~(isempty(value) || (isnumeric(value) && numel(value) == 2 && all(isfinite(value(:)))))
        error('validateDataParams:InvalidQCMetricVector', ...
            ['DataParams.registration.estimate.qcMetrics.%s must be empty or a ' ...
             '2-element finite numeric vector.'], name{1});
    end
end
if ~isempty(qc.scaleXY) && any(qc.scaleXY(:) <= 0)
    error('validateDataParams:InvalidScaleXY', ...
        'DataParams.registration.estimate.qcMetrics.scaleXY must contain positive values.');
end
end
