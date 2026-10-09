function fixture = setupSyntheticRegistrationFolder(rootFolder, varargin)
%SETUPSYNTHETICREGISTRATIONFOLDER Create a synthetic folder fixture.
%
%   fixture = setupSyntheticRegistrationFolder(rootFolder)
%   fixture = setupSyntheticRegistrationFolder(rootFolder, 'DatFormat', fmt)
%
%   'DatFormat' is 'header' (default: headered green.dat/red.dat) or
%   'legacySidecar' (headerless .dat plus legacy sidecar .mat), as written
%   by writeTestDat (.dat header Phase 5a).

p = inputParser;
addParameter(p, 'DatFormat', 'header', @(x) any(strcmpi(x, {'header', 'legacySidecar'})));
parse(p, varargin{:});
datFormat = p.Results.DatFormat;

if ~isfolder(rootFolder)
    mkdir(rootFolder);
end

Ny = 64;
Nx = 64;
Nt = 5;
frameRateHz = 10;
exposureMsec = 20;

[xg, yg] = meshgrid(1:Nx, 1:Ny);
refImage = exp(-((xg-20).^2 + (yg-20).^2)/80) + 0.7*exp(-((xg-45).^2 + (yg-40).^2)/120);
refImage = single(refImage);

shiftXY = [4 -3]; % [x y]
movingFrame = imtranslate(refImage, shiftXY, 'FillValues', 0);

greenData = repmat(movingFrame, 1, 1, Nt);
redData = repmat(movingFrame * 0.8, 1, 1, Nt);

writeTestDat(fullfile(rootFolder, 'green.dat'), greenData, frameRateHz, exposureMsec, 'Format', datFormat);
writeTestDat(fullfile(rootFolder, 'red.dat'), redData, frameRateHz, exposureMsec, 'Format', datFormat);

AcqInfoStream = struct();
AcqInfoStream.Width = Nx;
AcqInfoStream.Height = Ny;
AcqInfoStream.Length = Nt;
AcqInfoStream.FrameRateHz = frameRateHz;
AcqInfoStream.ExposureMsec = exposureMsec;
% Current-schema ImportedChannels (.dat header Phase 5a): the folder's
% channels are listed explicitly instead of inferred from file sizes.
for channel = {'green.dat', 'red.dat'}
    AcqInfoStream = appendImportedChannelInfo(AcqInfoStream, struct( ...
        'DatFile', channel{1}, 'Length', Nt, 'FrameRateHz', frameRateHz, ...
        'ExposureMsec', exposureMsec));
end
save(fullfile(rootFolder, 'AcqInfos.mat'), 'AcqInfoStream');

DataParams = createDataParams(rootFolder, 'overwrite', true);
DataParams.view.pixelSize_px_per_mm = 1;
saveDataParams(rootFolder, DataParams);

% DFR-20260819-010: pin an isolated fixture Rig to this dataset before
% addSession runs, so UMITRigStore.ensureDatasetRigAssociation resolves the
% pinned rigUUID instead of falling through to getActiveRig() -- the
% fixture never depends on (or mutates) which Rig is ambiently Active.
rigInfo = struct('rigID', ['ImgRegTestRig_' ...
    strrep(char(java.util.UUID.randomUUID()), '-', '')]);
rigStore = UMITRigStore.create(rigInfo);
UMITRigStore.assignDatasetRig(rootFolder, rigStore.getRigInfo().uuid);

projectInfo = struct( ...
    'projectName', ['Image Registration Test ' ...
    char(java.util.UUID.randomUUID())], ...
    'description', 'Temporary managed-reference test project.');
store = UMITProjectStore.create(projectInfo);
ProjectInfo = store.getProjectInfo();

subjectID = 'Subject_01';
subjectUUID = store.addSubject(struct('subjectID', subjectID));
sessionID = 'Session_01';
sessionUUID = store.addSession(subjectID, struct( ...
    'sessionID', sessionID, ...
    'processedDataFolder', rootFolder));

ImageReference = genImageReferenceStruct( ...
    refImage, ...
    'Name', 'Synthetic managed reference', ...
    'ProjectUUID', ProjectInfo.projectUUID, ...
    'ProjectName', ProjectInfo.projectName, ...
    'SubjectUUID', subjectUUID, ...
    'SubjectID', subjectID, ...
    'SessionUUID', sessionUUID, ...
    'SessionID', sessionID, ...
    'Description', 'Managed Image Reference payload description', ...
    'SourceFolder', rootFolder, ...
    'SourceFile', fullfile(rootFolder, 'green.dat'), ...
    'SourceFrame', 1, ...
    'SourceType', 'synthetic-test', ...
    'CreatedBy', 'setupSyntheticRegistrationFolder');
sourceReferenceFile = fullfile(rootFolder, 'reference_source.mat');
save(sourceReferenceFile, 'ImageReference', '-mat');

resourceInfo = struct( ...
    'displayName', 'Synthetic active reference', ...
    'description', 'Active managed Image Reference resource description');
imageReferenceUUID = store.addImageReference( ...
    subjectID, sourceReferenceFile, resourceInfo);
managedReference = store.getResource(imageReferenceUUID);

fixture = struct();
fixture.Ny = Ny;
fixture.Nx = Nx;
fixture.Nt = Nt;
fixture.ReferenceImage = refImage;
fixture.ShiftXY = shiftXY;
fixture.GreenData = greenData;
fixture.RedData = redData;
fixture.SaveFolder = rootFolder;
fixture.RigRoot = rigStore.RigRoot;
fixture.Store = store;
fixture.ProjectRoot = store.ProjectRoot;
fixture.ProjectUUID = ProjectInfo.projectUUID;
fixture.SubjectID = subjectID;
fixture.SubjectUUID = subjectUUID;
fixture.SessionID = sessionID;
fixture.SessionUUID = sessionUUID;
fixture.ImageReference = ImageReference;
fixture.ImageReferenceUUID = imageReferenceUUID;
fixture.ManagedReference = managedReference;
fixture.ManagedDescription = resourceInfo.description;
end
