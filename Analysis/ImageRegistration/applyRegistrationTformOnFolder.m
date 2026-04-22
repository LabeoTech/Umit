function outFile = applyRegistrationTformOnFolder(SaveFolder, varargin)
%APPLYREGISTRATIONTFORMONFOLDER Apply stored registration to all .dat files in a folder.
%
%   outFile = applyRegistrationTformOnFolder(SaveFolder)
%   outFile = applyRegistrationTformOnFolder(SaveFolder, 'Name', Value, ...)
%   info    = applyRegistrationTformOnFolder('pipelineInfo')
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
%       DataParamsFile            - Name of the DataParams MAT-file stored
%                                   in the folder. Default: 'DataParams.mat'
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
%   Output:
%       outFile - File manifest containing the modified .dat file names and
%                 the updated DataParams MAT-file.

if nargin == 1 && (ischar(SaveFolder) || (isstring(SaveFolder) && isscalar(SaveFolder))) ...
        && strcmpi(strtrim(char(string(SaveFolder))), 'pipelineInfo')
    outFile = localPipelineInfo();
    return
end

p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'SaveFolder', @isfolder);
addParameter(p, 'DataParamsFile', 'DataParams.mat', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'RequireUserConfirmation', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'OpenQCFigure', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'AllowReapply', false, @(x) islogical(x) && isscalar(x));
parse(p, SaveFolder, varargin{:});

SaveFolder = char(string(p.Results.SaveFolder));
dataParamsFile = char(string(p.Results.DataParamsFile));
requireUserConfirmation = p.Results.RequireUserConfirmation;
openQCFigure = p.Results.OpenQCFigure;
allowReapply = p.Results.AllowReapply;

if ~endsWith(dataParamsFile, '.mat', 'IgnoreCase', true)
    dataParamsFile = [dataParamsFile '.mat'];
end

dataParamsPath = fullfile(SaveFolder, dataParamsFile);
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

if requireUserConfirmation
    figHandle = [];
    if openQCFigure && isfield(DataParams.registration, 'qcFigureFile') && ...
            ~isempty(DataParams.registration.qcFigureFile)
        qcFigPath = fullfile(SaveFolder, char(string(DataParams.registration.qcFigureFile)));
        if isfile(qcFigPath)
            figHandle = openfig(qcFigPath, 'new', 'visible'); %#ok<NASGU>
        end
    end

    msg = sprintf(['Registration will irreversibly rewrite all .dat files in:\n\n%s\n\n' ...
        'Continue?'], SaveFolder);
    choice = questdlg(msg, 'Apply registration?', 'Continue', 'Cancel', 'Cancel');
    if ~strcmp(choice, 'Continue')
        outFile = {};
        return
    end
    confirmationMode = 'manual_confirmed';
    reviewedFlag = true;
else
    confirmationMode = 'forced_no_confirmation';
    reviewedFlag = DataParams.registration.isReviewed;
end

tform = DataParams.registration.tform;
outFile = {};

for iFile = 1:numel(datList)
    fileName = datList(iFile).name;
    datPath = fullfile(SaveFolder, fileName);
    md = loadMetaData(datPath);

    assert(isfield(md, 'Height') && isfield(md, 'Width') && isfield(md, 'Length'), ...
        'Umitoolbox:applyRegistrationTformOnFolder:InvalidMetadata', ...
        'Could not resolve Height/Width/Length for "%s".', datPath);

    ny = double(md.Height);
    nx = double(md.Width);
    nt = double(md.Length);
    Rfixed = imref2d([ny nx]);
    bytesPerElem = 4; % single

    fidIn = fopen(datPath, 'r');
    assert(fidIn ~= -1, ...
        'Umitoolbox:applyRegistrationTformOnFolder:FileOpenFailed', ...
        'Could not open "%s" for reading.', datPath);
    cIn = onCleanup(@() fclose(fidIn)); %#ok<NASGU>

    tmpPath = fullfile(SaveFolder, [fileName '.tmp']);
    fidOut = fopen(tmpPath, 'w');
    assert(fidOut ~= -1, ...
        'Umitoolbox:applyRegistrationTformOnFolder:FileOpenFailed', ...
        'Could not open temporary output file "%s".', tmpPath);
    cOut = onCleanup(@() fclose(fidOut)); %#ok<NASGU>

    fseek(fidIn, 0, 'bof');
    firstFrame = fread(fidIn, ny * nx, '*single');
    assert(numel(firstFrame) == ny * nx, ...
        'Umitoolbox:applyRegistrationTformOnFolder:InvalidDatFile', ...
        'Could not read the first frame from "%s".', datPath);
    firstFrame = reshape(firstFrame, ny, nx);
    nanMask = isnan(firstFrame);
    nanMaskWarped = imwarp(nanMask, tform, 'nearest', 'OutputView', Rfixed);

    for t = 1:nt
        fseek(fidIn, (t-1) * ny * nx * bytesPerElem, 'bof');
        frame = fread(fidIn, ny * nx, '*single');
        assert(numel(frame) == ny * nx, ...
            'Umitoolbox:applyRegistrationTformOnFolder:InvalidDatFile', ...
            'Could not read frame %d from "%s".', t, datPath);

        frame = reshape(frame, ny, nx);
        frame(nanMask) = 0;
        frame = imwarp(frame, tform, 'nearest', 'OutputView', Rfixed);
        frame(nanMaskWarped) = NaN;
        fwrite(fidOut, frame, 'single');
    end

    clear cIn cOut
    delete(datPath);
    movefile(tmpPath, datPath, 'f');
    outFile{end+1} = fileName; %#ok<AGROW>
end

DataParams.registration.isRegistered = true;
DataParams.registration.isReviewed = reviewedFlag;
DataParams.registration.appliedOn = datestr(now, 'yyyy-mm-dd HH:MM:SS');
DataParams.registration.appliedBy = mfilename;
DataParams.registration.confirmationMode = confirmationMode;
validateDataParams(DataParams);
save(dataParamsPath, 'DataParams');
outFile{end+1} = dataParamsFile;

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
    'DataParamsFile', ...
    'parameter', ...
    'Name of the DataParams MAT-file stored in the folder.', ...
    'kind', 'parameter', ...
    'default', 'DataParams.mat', ...
    'callType', 'namevalue');

info = PipelineManager.addInput( ...
    info, ...
    'RequireUserConfirmation', ...
    'parameter', ...
    'If true, asks the user to confirm before destructively applying the stored registration.', ...
    'kind', 'parameter', ...
    'default', true, ...
    'callType', 'namevalue');

info = PipelineManager.addInput( ...
    info, ...
    'OpenQCFigure', ...
    'parameter', ...
    'If true, opens the saved QC figure before the confirmation dialog when available.', ...
    'kind', 'parameter', ...
    'default', true, ...
    'callType', 'namevalue');

info = PipelineManager.addInput( ...
    info, ...
    'AllowReapply', ...
    'parameter', ...
    'If true, allows applying registration again even if the folder is already marked as registered.', ...
    'kind', 'parameter', ...
    'default', false, ...
    'callType', 'namevalue');

info = PipelineManager.addOutput( ...
    info, ...
    'outFile', ...
    'Unknown', ...
    'file', ...
    'File manifest containing the modified .dat files and updated DataParams MAT-file.', ...
    {'*.dat','DataParams.mat'}, ...
    1, ...
    'isData', false, ...
    'saveFileName', '');
end

end
