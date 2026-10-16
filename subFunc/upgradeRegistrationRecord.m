function registration = upgradeRegistrationRecord(registration, view)
%UPGRADEREGISTRATIONRECORD Return DataParams.registration as version 2 (in memory).
%
%   registration = upgradeRegistrationRecord(registration, view)
%
%   A version-2 record is returned unchanged. A version-1 (flat) record is
%   converted IN MEMORY so existing DataParams.mat files stay loadable;
%   nothing is written here (loadDataParams never writes, and the file is
%   only rewritten as version 2 when a writer saves it):
%       - referenceImage   -> reference (contentHash computed; name from
%                             referenceDescription; source kind 'legacy'
%                             with the old referenceFile as its path and the
%                             old imageReferenceUUID as resourceUUID)
%       - tform and QC     -> estimate (referenceHash = the reference hash)
%       - isRegistered...  -> apply (appliedReferenceHash = the reference
%                             hash when the folder was aligned)
%   An empty or missing record becomes newRegistrationRecord().
%
%   See also newRegistrationRecord, validateRegistrationRecord.

[empty, referenceTemplate, estimateTemplate] = newRegistrationRecord();
if ~isstruct(registration) || ~isscalar(registration) || isempty(fieldnames(registration))
    registration = empty;
    return
end
if isfield(registration, 'version') && isequal(registration.version, 2)
    return
end

old = registration;
registration = empty;
axisConvention = 'imageXY_topLeft';
if isstruct(view) && isfield(view, 'axisConvention') && ~isempty(view.axisConvention)
    axisConvention = char(string(view.axisConvention));
end

if isfield(old, 'referenceImage') && ~isempty(old.referenceImage)
    reference = referenceTemplate;
    [reference.contentHash, reference.image] = computeImageContentHash(old.referenceImage);
    reference.description = iText(old, 'referenceDescription');
    reference.name = reference.description;
    if isempty(reference.name)
        reference.name = 'Legacy reference';
    end
    reference.axisConvention = axisConvention;
    reference.imageSizeYX = double(size(reference.image));
    reference.selectedOn = iText(old, 'createdOn');
    reference.source = referenceSourceDescriptor('legacy', iText(old, 'referenceFile'));
    reference.source.resourceUUID = iText(old, 'imageReferenceUUID');
    registration.reference = reference;
end
referenceHash = '';
if ~isempty(registration.reference)
    referenceHash = registration.reference.contentHash;
end

if isfield(old, 'tform') && ~isempty(old.tform)
    estimate = estimateTemplate;
    estimate.tform = old.tform;
    for name = {'transformType', 'method', 'qcStatus', 'qcWarning', 'qcFigureFile', ...
            'qcPreviewImageFile', 'createdOn'}
        estimate.(name{1}) = iText(old, name{1});
    end
    if isfield(old, 'qcMetrics') && isstruct(old.qcMetrics)
        for name = fieldnames(estimate.qcMetrics)'
            if isfield(old.qcMetrics, name{1})
                estimate.qcMetrics.(name{1}) = old.qcMetrics.(name{1});
            end
        end
    end
    estimate.movingSource = iText(old, 'source');
    estimate.referenceHash = referenceHash;
    registration.estimate = estimate;
end

apply = registration.apply;
if isfield(old, 'isRegistered') && islogical(old.isRegistered) && isscalar(old.isRegistered)
    apply.isRegistered = old.isRegistered;
end
if isfield(old, 'isReviewed') && islogical(old.isReviewed) && isscalar(old.isReviewed)
    apply.isReviewed = old.isReviewed && ~isempty(registration.estimate);
end
if apply.isRegistered
    apply.appliedOn = iText(old, 'appliedOn');
    apply.appliedBy = iText(old, 'appliedBy');
    apply.confirmationMode = iText(old, 'confirmationMode');
    apply.appliedReferenceHash = referenceHash;
end
if isfield(old, 'manualAdjustment') && isstruct(old.manualAdjustment)
    apply.manualAdjustment = old.manualAdjustment;
end
apply.notes = iText(old, 'notes');
registration.apply = apply;
end

function text = iText(s, name)
text = '';
if isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
    if isdatetime(value)
        value.Format = 'yyyy-MM-dd HH:mm:ss';
    end
    text = char(string(value));
end
end
