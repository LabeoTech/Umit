function reference = resolveImageReferenceSource(source, opts)
%RESOLVEIMAGEREFERENCESOURCE Resolve a reference source to image, hash and descriptor.
%
%   reference = resolveImageReferenceSource(source)
%   reference = resolveImageReferenceSource(source, 'Interactive', false)
%
%   The single reference resolver used by the GUI and by pipelines. Given a
%   source (one of three kinds), it returns the reference image, its content
%   hash (computeImageContentHash) and a source descriptor, the human-aid
%   record stored in DataParams.registration.reference.source. Reading never
%   writes: not to the source folder, not to the target, not to the store.
%
%   SOURCE is a scalar struct with field 'kind' and:
%       kind 'project' - a reference of a subject library:
%           store       UMITProjectStore of the project
%           subjectID   the subject (ID)
%           resourceUUID (optional) a library reference; default: the
%                        subject's active reference
%       kind 'umitDat' - a UMIT .dat channel, directly or by its SaveFolder:
%           datPath     path of the .dat file, or
%           saveFolder + datFile (channel file name in that folder)
%       kind 'file'    - a standalone reference file:
%           path        a reference exported by UMIT (MAT holding an
%                       ImageReference variable; exportImageReferenceFile)
%                       or a plain image file (.tif/.tiff, .mat, .umt) read
%                       by readImageReferenceSource (§4 Q6)
%           variable, entry (optional) MAT variable / UMT entry for a plain
%                       image file with several candidates
%
%   Rules:
%       - umitDat: the .dat file's folder must contain a present, readable
%         AcqInfos.mat (with AcqInfoStream); nothing else is required.
%         Otherwise Umitoolbox:ImageRegistration:invalidDatSource.
%       - umitDat: the image is extracted by the fixed rule
%         'evenlySpacedMean' version 1 (ExtractionRule/ExtractionRuleVersion):
%         the mean of up to 100 evenly spaced frames, as
%         readImageReferenceSource reads a .dat.
%       - umitDat: if the source folder is already aligned
%         (DataParams.registration.apply.isRegistered), the selection is
%         ALLOWED with the warning Umitoolbox:ImageRegistration:alignedSource
%         (the source data are already transformed, the target is
%         interpolated a second time, and the source's committed reference).
%         The warning text is also returned in reference.source.warning so
%         it can be recorded. It is the only warning-level path.
%       - Paths in the descriptor are human aids only; UMIT never resolves
%         them again. Identity and verification use the image and its hash.
%
%   Output REFERENCE (scalar struct):
%       image, contentHash, name, description, axisConvention, imageSizeYX,
%       source (descriptor, see referenceSourceDescriptor), and sourceFile:
%       the file the image was made from when known (the payload's
%       sourceFile for a library or exported reference, the .dat or image
%       file otherwise); used by createRegistrationTform UseFile='auto'.
%
%   Name-Value options:
%       Interactive - false never opens a selection dialog (pipelines).
%                     Default: true.
%
%   See also computeImageContentHash, referenceSourceDescriptor,
%   readImageReferenceSource, exportImageReferenceFile.

arguments
    source (1,1) struct
    opts.Interactive (1,1) logical = true
end

if ~isfield(source, 'kind')
    error('Umitoolbox:ImageRegistration:invalidSource', ...
        'A reference source needs a ''kind'' (project, umitDat or file).');
end

switch char(string(source.kind))
    case 'project'
        reference = iResolveProject(source);
    case 'umitDat'
        reference = iResolveUmitDat(source);
    case 'file'
        reference = iResolveFile(source, opts.Interactive);
    otherwise
        error('Umitoolbox:ImageRegistration:invalidSource', ...
            'Unknown reference source kind "%s" (expected project, umitDat or file).', ...
            char(string(source.kind)));
end
end

% =========================================================================
% project
% =========================================================================
function reference = iResolveProject(source)
iRequire(source, {'store', 'subjectID'}, 'project');
store = source.store;
if ~isa(store, 'UMITProjectStore')
    error('Umitoolbox:ImageRegistration:invalidSource', ...
        'A project reference source needs a UMITProjectStore in ''store''.');
end
subjectID = char(string(source.subjectID));
SubjectInfo = store.getSubjectInfo(subjectID);

try
    if isfield(source, 'resourceUUID') && ~isempty(source.resourceUUID)
        resource = store.getResource(char(string(source.resourceUUID)));
    else
        resource = store.getActiveImageReference(subjectID);
    end
catch ME
    wrapped = MException('Umitoolbox:ImageRegistration:referenceNotFound', ...
        'Could not resolve the image reference of subject "%s": %s', subjectID, ME.message);
    throw(addCause(wrapped, ME));
end
if isempty(resource)
    error('Umitoolbox:ImageRegistration:referenceNotFound', ...
        ['Subject "%s" has no active image reference. Select one in Image ' ...
         'Reference Manager, or pass a reference source explicitly.'], subjectID);
end
if ~strcmp(resource.type, 'imageReference') || ~strcmpi(resource.ownerUUID, SubjectInfo.uuid)
    error('Umitoolbox:ImageRegistration:invalidSource', ...
        'Resource %s is not an image reference of subject "%s".', resource.uuid, subjectID);
end
if ~resource.fileExists || ~isfile(resource.absolutePath)
    error('Umitoolbox:ImageRegistration:referenceNotFound', ...
        'The image reference file is missing from the project store: %s', resource.absolutePath);
end
if ~strcmp(computeFileChecksum(resource.absolutePath), resource.checksum)
    error('Umitoolbox:ImageRegistration:referenceChecksumMismatch', ...
        'The image reference file does not match its registered checksum: %s', ...
        resource.absolutePath);
end
ImageReference = iLoadPayload(resource.absolutePath);

ProjectInfo = store.getProjectInfo();
if (~isempty(ImageReference.projectUUID) && ...
        ~strcmpi(ImageReference.projectUUID, ProjectInfo.projectUUID)) || ...
        (~isempty(ImageReference.subjectUUID) && ...
        ~strcmpi(ImageReference.subjectUUID, SubjectInfo.uuid))
    error('Umitoolbox:ImageRegistration:referenceIdentityMismatch', ...
        ['The image reference payload names another project or subject than ' ...
         'the library that holds it: %s'], resource.absolutePath);
end

descriptor = referenceSourceDescriptor('project', resource.absolutePath);
descriptor.projectUUID = ProjectInfo.projectUUID;
descriptor.projectName = ProjectInfo.projectName;
descriptor.subjectUUID = SubjectInfo.uuid;
descriptor.subjectID = SubjectInfo.subjectID;
descriptor.resourceUUID = resource.uuid;

reference = iReferenceFromPayload(ImageReference, descriptor);
reference.description = iFirstNonEmpty({resource.description, ImageReference.description, ...
    resource.displayName, ImageReference.name});
end

% =========================================================================
% umitDat
% =========================================================================
function reference = iResolveUmitDat(source)
if isfield(source, 'datPath') && ~isempty(source.datPath)
    datPath = char(string(source.datPath));
else
    iRequire(source, {'saveFolder', 'datFile'}, 'umitDat');
    datPath = fullfile(char(string(source.saveFolder)), char(string(source.datFile)));
end
[saveFolder, name, ext] = fileparts(datPath);
if ~strcmpi(ext, '.dat') || ~isfile(datPath)
    error('Umitoolbox:ImageRegistration:invalidDatSource', ...
        'The reference source "%s" is not an existing .dat file.', datPath);
end

acqInfosPath = fullfile(saveFolder, 'AcqInfos.mat');
acqOK = false;
if isfile(acqInfosPath)
    try
        acqOK = ~isempty(whos('-file', acqInfosPath, 'AcqInfoStream'));
    catch
        acqOK = false;
    end
end
if ~acqOK
    error('Umitoolbox:ImageRegistration:invalidDatSource', ...
        ['The reference source "%s" is not in a UMIT save folder: "%s" has no ' ...
         'readable AcqInfos.mat. Select a .dat imported by UMIT, or use a file source.'], ...
        datPath, saveFolder);
end

[image, ~] = readImageReferenceSource(datPath, 'maxFramesToAverage', 100, ...
    'Interactive', false);

descriptor = referenceSourceDescriptor('umitDat', datPath);
descriptor.saveFolder = saveFolder;
descriptor.extractionRule = 'evenlySpacedMean';
descriptor.extractionRuleVersion = 1;
[descriptor.sourceAligned, descriptor.sourceCommittedHash, descriptor.sourceCommittedName] = ...
    iSourceAlignment(saveFolder);

[contentHash, image] = computeImageContentHash(image);
reference = struct( ...
    'image', image, ...
    'contentHash', contentHash, ...
    'name', [name ext], ...
    'description', sprintf('Mean of up to 100 evenly spaced frames of %s', datPath), ...
    'axisConvention', 'imageXY_topLeft', ...
    'imageSizeYX', double(size(image)), ...
    'source', descriptor, ...
    'sourceFile', datPath);

if descriptor.sourceAligned
    message = sprintf(['The reference source "%s" comes from a folder whose data are ' ...
        'already aligned (to reference "%s", %s). Its image is already transformed, so ' ...
        'the target data will be interpolated a second time.'], datPath, ...
        descriptor.sourceCommittedName, UMITLedgerCore.shortHash(descriptor.sourceCommittedHash));
    reference.source.warning = message;
    warning('Umitoolbox:ImageRegistration:alignedSource', '%s', message);
end
end

function [aligned, committedHash, committedName] = iSourceAlignment(saveFolder)
% Read-only: loadDataParams never writes.
aligned = false;
committedHash = '';
committedName = '';
if ~isfile(fullfile(saveFolder, 'DataParams.mat'))
    return
end
DataParams = loadDataParams(saveFolder);
registration = DataParams.registration;
if ~registration.apply.isRegistered
    return
end
aligned = true;
committedHash = registration.apply.appliedReferenceHash;
if ~isempty(registration.reference)
    committedName = registration.reference.name;
end
end

% =========================================================================
% file
% =========================================================================
function reference = iResolveFile(source, interactive)
iRequire(source, {'path'}, 'file');
filePath = char(string(source.path));
if ~isfile(filePath)
    error('Umitoolbox:ImageRegistration:invalidSource', ...
        'The reference file does not exist: %s', filePath);
end
[~, name, ext] = fileparts(filePath);
if strcmpi(ext, '.dat')
    error('Umitoolbox:ImageRegistration:invalidSource', ...
        'A .dat reference source uses the source kind ''umitDat'': %s', filePath);
end

if strcmpi(ext, '.mat') && ~isempty(whos('-file', filePath, 'ImageReference'))
    ImageReference = iLoadPayload(filePath);
    descriptor = referenceSourceDescriptor('file', filePath);
    descriptor.fileFormat = 'umitReference';
    for fieldName = {'projectUUID', 'projectName', 'subjectUUID', 'subjectID', ...
            'sessionUUID', 'sessionID'}
        descriptor.(fieldName{1}) = char(string(ImageReference.(fieldName{1})));
    end
    reference = iReferenceFromPayload(ImageReference, descriptor);
    return
end

options = {'Interactive', interactive};
if isfield(source, 'variable') && ~isempty(source.variable)
    options = [options, {'Variable', char(string(source.variable))}];
end
if isfield(source, 'entry') && ~isempty(source.entry)
    options = [options, {'Entry', char(string(source.entry))}];
end
[image, sourceInfo] = readImageReferenceSource(filePath, options{:});
descriptor = referenceSourceDescriptor('file', filePath);
descriptor.fileFormat = 'image';
if isfield(sourceInfo, 'selectedEntry') && ~isempty(sourceInfo.selectedEntry)
    descriptor.entry = char(string(sourceInfo.selectedEntry));
elseif isfield(sourceInfo, 'selectedVariable') && ~isempty(sourceInfo.selectedVariable)
    descriptor.entry = char(string(sourceInfo.selectedVariable));
end
[contentHash, image] = computeImageContentHash(image);
reference = struct( ...
    'image', image, ...
    'contentHash', contentHash, ...
    'name', [name ext], ...
    'description', char(string(sourceInfo.processingNotes)), ...
    'axisConvention', 'imageXY_topLeft', ...
    'imageSizeYX', double(size(image)), ...
    'source', descriptor, ...
    'sourceFile', filePath);
end

% =========================================================================
% Helpers
% =========================================================================
function reference = iReferenceFromPayload(ImageReference, descriptor)
reference = struct( ...
    'image', ImageReference.image, ...
    'contentHash', ImageReference.contentHash, ...
    'name', char(string(ImageReference.name)), ...
    'description', char(string(ImageReference.description)), ...
    'axisConvention', char(string(ImageReference.axisConvention)), ...
    'imageSizeYX', double(size(ImageReference.image)), ...
    'source', descriptor, ...
    'sourceFile', char(string(ImageReference.sourceFile)));
end

function ImageReference = iLoadPayload(filePath)
try
    loaded = load(filePath, 'ImageReference', '-mat');
catch ME
    wrapped = MException('Umitoolbox:ImageRegistration:invalidSource', ...
        'Could not read the image reference file "%s": %s', filePath, ME.message);
    throw(addCause(wrapped, ME));
end
if ~isfield(loaded, 'ImageReference')
    error('Umitoolbox:ImageRegistration:invalidSource', ...
        'The file holds no ImageReference variable: %s', filePath);
end
try
    ImageReference = normalizeImageReferenceStruct(loaded.ImageReference);
catch ME
    wrapped = MException('Umitoolbox:ImageRegistration:invalidSource', ...
        'The image reference in "%s" is invalid: %s', filePath, ME.message);
    throw(addCause(wrapped, ME));
end
end

function iRequire(source, fields, kind)
missing = fields(~isfield(source, fields));
if ~isempty(missing)
    error('Umitoolbox:ImageRegistration:invalidSource', ...
        'A %s reference source needs field(s): %s.', kind, strjoin(missing, ', '));
end
end

function text = iFirstNonEmpty(candidates)
text = '';
for iCandidate = 1:numel(candidates)
    candidate = strtrim(char(string(candidates{iCandidate})));
    if ~isempty(candidate)
        text = candidate;
        return
    end
end
end
