function varargout = applyRegistrationTformOnFolder(SaveFolder, varargin)
%APPLYREGISTRATIONTFORMONFOLDER Apply stored registration to all .dat files in a folder.
%
%   applyRegistrationTformOnFolder(SaveFolder)
%   applyRegistrationTformOnFolder(SaveFolder, 'Name', Value, ...)
%   info = applyRegistrationTformOnFolder('pipelineInfo')
%
%   This function reads the currently stored registration transform from the
%   folder DataParams file and destructively applies it to all .dat files
%   in the folder.
%
%   IMPORTANT:
%       This operation rewrites the folder .dat files in place.
%       Due to interpolation and potential out-of-frame data loss, the
%       change should be considered irreversible in practice.
%
%   Inputs:
%       SaveFolder - Folder containing DataParams.mat and the target .dat
%                    files.
%
%   Name-Value parameters:
%       RequireUserConfirmation   - Logical scalar. If true, the saved QC
%                                   figure is opened when available and the
%                                   user is asked to confirm before the
%                                   destructive apply. Default: true
%       OpenQCFigure              - Logical scalar. If true and the QC
%                                   figure exists, it is opened before the
%                                   confirmation dialog. Default: true
%       AllowReapply              - Logical scalar. If false, the function
%                                   errors when the folder is already marked
%                                   as registered. Default: false
%
%   Notes:
%       - This function has no normal runtime output.
%       - After completion, it prints the modified .dat files and the
%         updated DataParams file to the command window.
%       - A pipelineInfo structure can still be queried with:
%             info = applyRegistrationTformOnFolder('pipelineInfo')

if nargin == 1 && (ischar(SaveFolder) || (isstring(SaveFolder) && isscalar(SaveFolder))) ...
        && strcmpi(strtrim(char(string(SaveFolder))), 'pipelineInfo')
    if nargout > 0
        varargout{1} = localPipelineInfo();
    else
        disp(localPipelineInfo());
    end
    return
end

p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'SaveFolder', @isfolder);
addParameter(p, 'RequireUserConfirmation', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'OpenQCFigure', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'AllowReapply', false, @(x) islogical(x) && isscalar(x));
parse(p, SaveFolder, varargin{:});

SaveFolder = char(string(p.Results.SaveFolder));
requireUserConfirmation = p.Results.RequireUserConfirmation;
openQCFigure = p.Results.OpenQCFigure;
allowReapply = p.Results.AllowReapply;

dataParamsPath = fullfile(SaveFolder, 'DataParams.mat');
assert(isfile(dataParamsPath), ...
    'Umitoolbox:applyRegistrationTformOnFolder:MissingDataParams', ...
    'DataParams file not found: "%s".', dataParamsPath);

S = load(dataParamsPath, 'DataParams');
assert(isfield(S, 'DataParams') && isstruct(S.DataParams) && isscalar(S.DataParams), ...
    'Umitoolbox:applyRegistrationTformOnFolder:InvalidDataParams', ...
    'DataParams file does not contain a valid scalar DataParams struct.');
DataParams = S.DataParams;
validateDataParams(DataParams);

assert(isfield(DataParams, 'registration') && isstruct(DataParams.registration), ...
    'Umitoolbox:applyRegistrationTformOnFolder:MissingRegistration', ...
    'DataParams.registration is missing or invalid.');
assert(~isempty(DataParams.registration.tform), ...
    'Umitoolbox:applyRegistrationTformOnFolder:MissingTform', ...
    'No registration tform is stored in DataParams.registration.tform.');

if DataParams.registration.isRegistered && ~allowReapply
    error('Umitoolbox:applyRegistrationTformOnFolder:AlreadyRegistered', ...
        ['This folder is already marked as registered. Reapplication is blocked ' ...
         'by default because the operation is destructive.']);
end

datList = dir(fullfile(SaveFolder, '*.dat'));
assert(~isempty(datList), ...
    'Umitoolbox:applyRegistrationTformOnFolder:NoDatFiles', ...
    'No .dat files were found in "%s".', SaveFolder);

% Validate every target before touching any of them, and before prompting.
% This operation is destructive, so it must be all-or-nothing: gating inside
% the rewrite loop would abort partway, leaving some files transformed,
% others not, and DataParams still marked unregistered.
datPlan = iPreflightDatFiles(datList, SaveFolder, DataParams);

if requireUserConfirmation
    if openQCFigure && isfield(DataParams.registration, 'qcFigureFile') && ...
            ~isempty(DataParams.registration.qcFigureFile)

        qcFigPath = char(string(DataParams.registration.qcFigureFile));
        if isempty(fileparts(qcFigPath))
            qcFigPath = fullfile(SaveFolder, qcFigPath);
        end

        if isfile(qcFigPath)
            openfig(qcFigPath, 'new', 'visible');
        end
    end

    msg = sprintf(['Registration will irreversibly rewrite all .dat files in:\n\n%s\n\n' ...
        'Continue?'], SaveFolder);
    choice = questdlg(msg, 'Apply registration?', 'Continue', 'Cancel', 'Cancel');
    if ~strcmp(choice, 'Continue')
        fprintf('Registration cancelled by user.\n');
        return
    end
    confirmationMode = 'manual_confirmed';
    reviewedFlag = true;
else
    confirmationMode = 'forced_no_confirmation';
    reviewedFlag = DataParams.registration.isReviewed;
end

tform = DataParams.registration.tform;
modifiedFiles = cell(0,1);

for iFile = 1:numel(datPlan)
    fileName = datPlan(iFile).name;
    datPath = datPlan(iFile).path;

    ny = datPlan(iFile).ny;
    nx = datPlan(iFile).nx;
    nt = datPlan(iFile).nt;
    datatype = datPlan(iFile).datatype;
    Rfixed = imref2d([ny nx]);
    bytesPerElem = getByteSize(datatype);

    fidIn = fopen(datPath, 'r');
    assert(fidIn ~= -1, ...
        'Umitoolbox:applyRegistrationTformOnFolder:FileOpenFailed', ...
        'Could not open "%s" for reading.', datPath);
    cIn = onCleanup(@() safeFclose(fidIn));

    tmpPath = fullfile(SaveFolder, [fileName '.tmp']);
    fidOut = fopen(tmpPath, 'w');
    assert(fidOut ~= -1, ...
        'Umitoolbox:applyRegistrationTformOnFolder:FileOpenFailed', ...
        'Could not open temporary output file "%s".', tmpPath);
    cOut = onCleanup(@() safeFclose(fidOut));

    fseek(fidIn, 0, 'bof');
    firstFrame = fread(fidIn, ny * nx, ['*' datatype]);
    assert(numel(firstFrame) == ny * nx, ...
        'Umitoolbox:applyRegistrationTformOnFolder:InvalidDatFile', ...
        'Could not read the first frame from "%s".', datPath);
    firstFrame = reshape(firstFrame, ny, nx);
    nanMask = isnan(firstFrame);
    nanMaskWarped = imwarp(nanMask, tform, 'nearest', 'OutputView', Rfixed);

    for t = 1:nt
        fseek(fidIn, (t-1) * ny * nx * bytesPerElem, 'bof');
        frame = fread(fidIn, ny * nx, ['*' datatype]);
        assert(numel(frame) == ny * nx, ...
            'Umitoolbox:applyRegistrationTformOnFolder:InvalidDatFile', ...
            'Could not read frame %d from "%s".', t, datPath);

        frame = reshape(frame, ny, nx);
        frame(nanMask) = 0;
        frame = imwarp(frame, tform, 'nearest', 'OutputView', Rfixed);
        frame(nanMaskWarped) = NaN;
        fwrite(fidOut, frame, datatype);
    end

    clear cIn cOut
    iReplaceFileSafely(tmpPath, datPath);
    modifiedFiles{end+1,1} = fileName; %#ok<AGROW>
end

DataParams.registration.isRegistered = true;
DataParams.registration.isReviewed = reviewedFlag;
DataParams.registration.appliedOn = char( ...
    datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
DataParams.registration.appliedBy = mfilename;
DataParams.registration.confirmationMode = confirmationMode;
saveDataParams(SaveFolder, DataParams);

fprintf('\nRegistration applied to folder:\n%s\n', SaveFolder);
fprintf('Modified .dat files:\n');
for iFile = 1:numel(modifiedFiles)
    fprintf('  - %s\n', modifiedFiles{iFile});
end
fprintf('Updated DataParams:\n');
fprintf('  - %s\n\n', dataParamsPath);

end

function datPlan = iPreflightDatFiles(datList, SaveFolder, DataParams)
%IPREFLIGHTDATFILES Validate every target .dat before any file is rewritten.
%
% Returns one struct per file carrying the metadata the rewrite loop needs,
% so loadMetaData is not called twice per file. Any unsupported file aborts
% the whole operation: this path rewrites data in place, so a partial run is
% worse than no run.

refSizeYX = iResolveReferenceSizeYX(DataParams);

datPlan = struct('name', {}, 'path', {}, 'ny', {}, 'nx', {}, 'nt', {}, ...
    'datatype', {});

for iFile = 1:numel(datList)
    fileName = datList(iFile).name;
    datPath = fullfile(SaveFolder, fileName);
    md = loadMetaData(datPath);

    if ~isfield(md, 'Height') || ~isfield(md, 'Width') || ~isfield(md, 'datLength')
        error('Umitoolbox:applyRegistrationTformOnFolder:InvalidMetadata', ...
            'Could not resolve Height/Width/datLength for "%s".', datPath);
    end

    % Datatype is not gated here: loadMetaData already refuses any .dat that
    % is not single precision. It is still read from the metadata rather than
    % hardcoded, so the read/write calls below stay correct by construction
    % if that invariant is ever relaxed.
    datatype = 'single';
    if isfield(md, 'Datatype') && ~isempty(md.Datatype)
        datatype = char(string(md.Datatype));
    end

    % Layout gate. Event-split .dat is unsupported legacy on this path; its
    % frame count is not md.datLength, so it would be rewritten incorrectly.
    dimNames = {'Y','X','T'};
    if isfield(md, 'dim_names') && ~isempty(md.dim_names)
        dimNames = cellstr(string(md.dim_names));
    end
    if ~isequal(dimNames(:).', {'Y','X','T'})
        error('Umitoolbox:applyRegistrationTformOnFolder:UnsupportedLayout', ...
            ['File "%s" has dimensions {%s}. Folder registration only ' ...
             'supports continuous Y-X-T .dat files; event-split .dat data ' ...
             'are not supported on this path.'], ...
            fileName, strjoin(dimNames(:).', ','));
    end

    ny = double(md.Height);
    nx = double(md.Width);

    % Size gate. createRegistrationTform resizes the moving image to the
    % reference before estimating, so the stored transform is expressed in
    % reference pixel units. Applying it at a different native size would
    % silently shift and scale the data by the wrong amount.
    if ~isequal([ny nx], refSizeYX)
        error('Umitoolbox:applyRegistrationTformOnFolder:ReferenceSizeMismatch', ...
            ['File "%s" is %dx%d, but the stored transform was estimated ' ...
             'against a %dx%d reference. Re-run createRegistrationTform for ' ...
             'this folder before applying.'], ...
            fileName, ny, nx, refSizeYX(1), refSizeYX(2));
    end

    datPlan(end+1) = struct( ...
        'name', fileName, ...
        'path', datPath, ...
        'ny', ny, ...
        'nx', nx, ...
        'nt', double(md.datLength), ...
        'datatype', datatype); %#ok<AGROW>
end

end

function refSizeYX = iResolveReferenceSizeYX(DataParams)
%IRESOLVEREFERENCESIZEYX Frame size the stored transform was estimated against.

refSizeYX = [];

if isfield(DataParams, 'registration') && ...
        isfield(DataParams.registration, 'referenceImage') && ...
        ~isempty(DataParams.registration.referenceImage)
    refSizeYX = double(size(DataParams.registration.referenceImage, [1 2]));
elseif isfield(DataParams, 'view') && ...
        isfield(DataParams.view, 'imageSizeYX') && ...
        ~isempty(DataParams.view.imageSizeYX)
    refSizeYX = double(DataParams.view.imageSizeYX(:).');
end

if numel(refSizeYX) ~= 2
    error('Umitoolbox:applyRegistrationTformOnFolder:UnknownReferenceSize', ...
        ['Neither DataParams.registration.referenceImage nor ' ...
         'DataParams.view.imageSizeYX records the frame size the stored ' ...
         'transform was estimated against, so the transform cannot be ' ...
         'verified. Re-run createRegistrationTform for this folder.']);
end

end

function iReplaceFileSafely(tmpPath, destPath)
%IREPLACEFILESAFELY Install tmpPath as destPath without a data-loss window.
%
% Mirrors the sequence used by saveMatAtomic: move the existing file aside,
% install the replacement, then drop the backup. The destination is never
% deleted before its replacement is known to be in place, so a failed move
% cannot leave the folder without the data.

backupPath = [destPath '.bak'];
backupCreated = false;

if isfile(destPath)
    [ok, message] = movefile(destPath, backupPath, 'f');
    if ~ok
        error('Umitoolbox:applyRegistrationTformOnFolder:BackupFailed', ...
            'Could not back up "%s" before replacing it: %s', destPath, message);
    end
    backupCreated = true;
end

[ok, message] = movefile(tmpPath, destPath, 'f');
if ~ok
    if backupCreated && isfile(backupPath)
        movefile(backupPath, destPath, 'f');
    end
    error('Umitoolbox:applyRegistrationTformOnFolder:ReplaceFailed', ...
        'Could not install the registered data as "%s": %s', destPath, message);
end

if backupCreated && isfile(backupPath)
    delete(backupPath);
end

end

function info = localPipelineInfo()
%LOCALPIPELINEINFO Return PipelineManager metadata for applyRegistrationTformOnFolder.

info = PipelineManager.createPipelineInfo( ...
    mfilename, ...
    'Destructively apply the stored registration transform to all .dat files in a folder.');

info = PipelineManager.addInput( ...
    info, ...
    'SaveFolder', ...
    'SaveFolder', ...
    'Folder containing DataParams.mat and the target .dat files.', ...
    'kind', 'input', ...
    'position', 1, ...
    'callType', 'positional', ...
    'isData', false);

info = PipelineManager.addInput( ...
    info, ...
    'RequireUserConfirmation', ...
    'parameter', ...
    'If true, asks the user to confirm before destructively applying the stored registration.', ...
    'kind', 'parameter', ...
    'default', true, ...
    'allowed', [true false], ...
    'callType', 'namevalue');

info = PipelineManager.addInput( ...
    info, ...
    'OpenQCFigure', ...
    'parameter', ...
    'If true, opens the saved QC figure before the confirmation dialog when available.', ...
    'kind', 'parameter', ...
    'default', true, ...
    'allowed', [true false], ...
    'callType', 'namevalue');

info = PipelineManager.addInput( ...
    info, ...
    'AllowReapply', ...
    'parameter', ...
    'If true, allows applying registration again even if the folder is already marked as registered.', ...
    'kind', 'parameter', ...
    'default', false, ...
    'allowed', [true false], ...
    'callType', 'namevalue');

info = PipelineManager.addOutput( ...
    info, ...
    'rewrittenDatFiles', ...
    {'ImageTimeSeries','ProcessedData'}, ...
    'file', ...
    'Every .dat file in SaveFolder, rewritten in place with the registration applied.', ...
    '', ...
    1, ...
    'isData', false);

info = PipelineManager.addOutput( ...
    info, ...
    'dataParams', ...
    'SaveFolder', ...
    'file', ...
    'DataParams.mat, updated to record that registration was applied.', ...
    'DataParams.mat', ...
    2, ...
    'isData', false);
end
