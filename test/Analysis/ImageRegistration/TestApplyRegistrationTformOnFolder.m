classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestApplyRegistrationTformOnFolder < matlab.unittest.TestCase
    %TESTAPPLYREGISTRATIONTFORMONFOLDER Unit tests for applyRegistrationTformOnFolder.
    %
    %   Since project redesign Phase 1 the pipeline apply is a wrapper around
    %   the canonical backend applyImageAlignmentToFolder (ImageAlignmentTool):
    %   same folder state, same guards, staged all-or-nothing commit.

    properties
        TempFolder
        Fixture
    end

    methods (TestMethodSetup)
        function setup(testCase)
            testCase.TempFolder = fullfile(tempdir, ['TestApplyRegistrationTformOnFolder_' char(java.util.UUID.randomUUID)]);
            testCase.Fixture = setupSyntheticRegistrationFolder(testCase.TempFolder);
            createRegistrationTform(testCase.TempFolder, 'ShowFigure', false);
        end
    end

    methods (TestMethodTeardown)
        function teardown(testCase)
            if isstruct(testCase.Fixture) && ...
                    isfield(testCase.Fixture, 'ProjectRoot') && ...
                    isfolder(testCase.Fixture.ProjectRoot)
                rmdir(testCase.Fixture.ProjectRoot, 's');
            end
            if isstruct(testCase.Fixture) && ...
                    isfield(testCase.Fixture, 'RigRoot') && ...
                    isfolder(testCase.Fixture.RigRoot)
                rmdir(testCase.Fixture.RigRoot, 's');
            end
            if isfolder(testCase.TempFolder)
                rmdir(testCase.TempFolder, 's');
            end
        end
    end

    methods (Test)
        function testPipelineManagerExecutesRegistrationApplyNode(testCase)
            beforeGreen = loadData(fullfile(testCase.TempFolder, 'green.dat'));
            parameters = struct('RequireUserConfirmation', false, ...
                'OpenQCFigure', false);
            pm = buildPMForScenario(testCase.TempFolder, ...
                'applyRegistrationTformOnFolder', 'auto', ...
                'Parameters', parameters);
            result = pm.executePipeline('PrintSummary', false);

            testCase.verifyEqual(result.status, "completed");
            testCase.verifyEqual(sort(result.createdFiles.FileName), ...
                sort(["green.dat"; "red.dat"; "DataParams.mat"]));
            testCase.verifyTrue(all(result.createdFiles.FileExists));
            afterGreen = loadData(fullfile(testCase.TempFolder, 'green.dat'));
            testCase.verifyFalse(isequal(beforeGreen, afterGreen));
            S = load(fullfile(testCase.TempFolder, 'DataParams.mat'), 'DataParams');
            testCase.verifyTrue(S.DataParams.registration.apply.isRegistered);
        end

        function testApplyUpdatesDatFilesAndDataParams(testCase)
            beforeGreen = loadData(fullfile(testCase.TempFolder, 'green.dat'));
            beforeRed = loadData(fullfile(testCase.TempFolder, 'red.dat'));

            applyRegistrationTformOnFolder(testCase.TempFolder, ...
                'RequireUserConfirmation', false);

            afterGreen = loadData(fullfile(testCase.TempFolder, 'green.dat'));
            afterRed = loadData(fullfile(testCase.TempFolder, 'red.dat'));

            testCase.verifyFalse(isequal(beforeGreen, afterGreen));
            testCase.verifyFalse(isequal(beforeRed, afterRed));

            S = load(fullfile(testCase.TempFolder, 'DataParams.mat'), 'DataParams');
            DP = S.DataParams;
            validateDataParams(DP);
            testCase.verifyTrue(DP.registration.apply.isRegistered);
            testCase.verifyEqual(char(string(DP.registration.apply.appliedBy)), 'applyRegistrationTformOnFolder');
            testCase.verifyEqual(char(string(DP.registration.apply.confirmationMode)), 'forced_no_confirmation');
            % The commitment is the content hash the estimate was made against.
            testCase.verifyEqual(DP.registration.apply.appliedReferenceHash, ...
                computeImageContentHash(testCase.Fixture.ReferenceImage));
            testCase.verifyEqual(DP.registration.apply.appliedReferenceHash, ...
                DP.registration.estimate.referenceHash);
            testCase.verifyEqual(DP.registration.reference.source.kind, 'project');
            testCase.verifyEqual(DP.registration.reference.source.resourceUUID, ...
                testCase.Fixture.ImageReferenceUUID);
        end

        function testAlignedFolderRefusesEveryReapply(testCase)
            % Brief section 4 Q1: an aligned folder is committed; re-applying
            % is refused (no AllowReapply any more), before any change.
            applyRegistrationTformOnFolder(testCase.TempFolder, 'RequireUserConfirmation', false);
            before = snapshotRegistrationFiles(testCase.TempFolder);

            try
                applyRegistrationTformOnFolder(testCase.TempFolder, 'RequireUserConfirmation', false);
                testCase.verifyFail('Expected the re-apply to be refused.');
            catch ME
                testCase.verifyEqual(ME.identifier, 'Umitoolbox:ImageRegistration:alreadyAligned');
                testCase.verifySubstring(ME.message, testCase.TempFolder);
                testCase.verifySubstring(ME.message, 'fresh import');
            end
            testCase.verifyError(@() applyRegistrationTformOnFolder(testCase.TempFolder, ...
                'RequireUserConfirmation', false, 'AllowReapply', true), ...
                'MATLAB:InputParser:UnmatchedParameter');
            verifyRegistrationSnapshot(testCase, testCase.TempFolder, before);
        end

        function testMissingTformErrors(testCase)
            S = load(fullfile(testCase.TempFolder, 'DataParams.mat'), 'DataParams');
            S.DataParams.registration.estimate = [];
            save(fullfile(testCase.TempFolder, 'DataParams.mat'), '-struct', 'S');

            testCase.verifyError(@() applyRegistrationTformOnFolder(testCase.TempFolder, ...
                'RequireUserConfirmation', false), ...
                'Umitoolbox:applyRegistrationTformOnFolder:MissingTform');
        end

        function testProcessesHeaderedDat(testCase)
            % .dat header Phase 4e-1: a headered file is registered like a
            % headerless one, and every output is headered.
            fx = testCase.Fixture;
            iRewriteWithHeader(fullfile(testCase.TempFolder, 'red.dat'), [fx.Ny fx.Nx fx.Nt], 10);
            before = loadData(fullfile(testCase.TempFolder, 'red.dat'));

            applyRegistrationTformOnFolder(testCase.TempFolder, ...
                'RequireUserConfirmation', false, 'OpenQCFigure', false);

            for name = {'red.dat', 'green.dat'}
                f = fullfile(testCase.TempFolder, name{1});
                testCase.verifyTrue(isDatWithHeader(f), [name{1} ' must be headered']);
                testCase.verifyTrue(readDatHeader(f).writeComplete);
            end
            testCase.verifyEqual(readDatHeader(fullfile(testCase.TempFolder, 'red.dat')).frameRateHz, 10);
            testCase.verifyFalse(isequal(loadData(fullfile(testCase.TempFolder, 'red.dat')), before));
        end

        function testCancellationRaisesBeforeMutation(testCase)
            before = snapshotRegistrationFiles(testCase.TempFolder);
            installCancelQuestdlgMock(testCase);

            testCase.verifyError(@() applyRegistrationTformOnFolder( ...
                testCase.TempFolder, 'OpenQCFigure', false), ...
                'Umitoolbox:applyRegistrationTformOnFolder:Cancelled');

            verifyRegistrationSnapshot(testCase, testCase.TempFolder, before);
        end

        function testPipelineManagerReportsConfirmationCancellationAsFailure(testCase)
            before = snapshotRegistrationFiles(testCase.TempFolder);
            installCancelQuestdlgMock(testCase);
            pm = buildPMForScenario(testCase.TempFolder, ...
                'applyRegistrationTformOnFolder', 'auto', ...
                'Parameters', struct('OpenQCFigure', false));
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture( ...
                'PipelineManager:StepFailed'));

            result = pm.executePipeline('PrintSummary', false);

            testCase.verifyEqual(result.status, "failed");
            verifyRegistrationSnapshot(testCase, testCase.TempFolder, before);
        end

        function testLeavesNoTemporaryResidue(testCase)
            applyRegistrationTformOnFolder(testCase.TempFolder, ...
                'RequireUserConfirmation', false);

            testCase.verifyEmpty(dir(fullfile(testCase.TempFolder, '*.tmp')));
            testCase.verifyEmpty(dir(fullfile(testCase.TempFolder, '*.bak')));
            testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, 'green.dat')));
            testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, 'red.dat')));
        end

        function testRejectsInconsistentViewSize(testCase)
            % DataParams.view must describe the current frame; a stale view
            % would give a wrong view/mask after alignment.
            DataParams = loadDataParams(testCase.TempFolder);
            smallerYX = [32 32];
            DataParams.view.imageSizeYX = smallerYX;
            DataParams.mask.logical = [];
            saveDataParams(testCase.TempFolder, DataParams);
            before = snapshotRegistrationFiles(testCase.TempFolder);

            testCase.verifyError(@() applyRegistrationTformOnFolder( ...
                testCase.TempFolder, 'RequireUserConfirmation', false), ...
                'ImageAlignmentTool:ViewSizeMismatch');

            % Pre-flight must abort before any file is touched.
            verifyRegistrationSnapshot(testCase, testCase.TempFolder, before);
            testCase.verifyEmpty(dir(fullfile(testCase.TempFolder, '*.tmp')));
            testCase.verifyEmpty(dir(fullfile(testCase.TempFolder, '*.bak')));
        end

        function testRejectsLayoutsThatDoNotStartWithYX(testCase)
            % Only data whose first two axes are Y, X can be registered with a
            % 2-D transform; anything else is refused before any file is touched.
            fx = testCase.Fixture;
            writeLegacyDat(testCase.TempFolder, 'rogueLayout', ...
                zeros(2, fx.Ny, fx.Nx, 'single'), 'single', {'T','Y','X'});
            before = snapshotRegistrationFiles(testCase.TempFolder);

            testCase.verifyError(@() applyRegistrationTformOnFolder( ...
                testCase.TempFolder, 'RequireUserConfirmation', false), ...
                'ImageAlignmentTool:UnsupportedLayout');

            verifyRegistrationSnapshot(testCase, testCase.TempFolder, before);
            testCase.verifyEmpty(dir(fullfile(testCase.TempFolder, '*.tmp')));
            testCase.verifyEmpty(dir(fullfile(testCase.TempFolder, '*.bak')));
        end

        function testRegistersEveryYXLayout(testCase)
            % Y-X, Y-X-E and Y-X-F files are registered plane by plane, with
            % the same result as the corresponding planes of the recording.
            fx = testCase.Fixture;
            green = single(loadData(fullfile(testCase.TempFolder, 'green.dat')));
            writeTestDat(fullfile(testCase.TempFolder, 'frame.dat'), green(:, :, 1), 10, ...
                'DimNames', {'Y','X'});
            writeTestDat(fullfile(testCase.TempFolder, 'perEvent.dat'), green(:, :, 1:3), 10, ...
                'DimNames', {'Y','X','E'});
            writeTestDat(fullfile(testCase.TempFolder, 'maps.dat'), green(:, :, 1:2), 10, ...
                'DimNames', {'Y','X','F'});

            applyRegistrationTformOnFolder(testCase.TempFolder, ...
                'RequireUserConfirmation', false, 'OpenQCFigure', false);

            registeredGreen = single(loadData(fullfile(testCase.TempFolder, 'green.dat')));
            testCase.verifyEqual(single(loadData(fullfile(testCase.TempFolder, 'frame.dat'))), ...
                registeredGreen(:, :, 1));
            testCase.verifyEqual(single(loadData(fullfile(testCase.TempFolder, 'perEvent.dat'))), ...
                registeredGreen(:, :, 1:3));
            testCase.verifyEqual(single(loadData(fullfile(testCase.TempFolder, 'maps.dat'))), ...
                registeredGreen(:, :, 1:2));
            testCase.verifyEqual(double(loadMetaData( ...
                fullfile(testCase.TempFolder, 'maps.dat')).dimSizes(:).'), [fx.Ny fx.Nx 2]);
        end

        function testRegistersImageUMTAndLeavesOtherUMTsAlone(testCase)
            fx = testCase.Fixture;
            green = single(loadData(fullfile(testCase.TempFolder, 'green.dat')));
            imageUMT = genUMTStruct(green(:, :, 1:2), 'kind', 'image', ...
                'entryName', 'main', 'dimNames', {'Y','X','T'}, ...
                'meta', struct('FrameRateHz', 10));
            imageUMT = genUMTStruct(imageUMT, 'value', green(:, :, 1), ...
                'entryName', 'frame', 'dimNames', {'Y','X'});
            saveData(fullfile(testCase.TempFolder, 'maps.umt'), imageUMT);
            roiUMT = genUMTStruct(rand(2, 6, 'single'), 'kind', 'roi', ...
                'entryName', 'traces', 'dimNames', {'ROI','T'});
            saveData(fullfile(testCase.TempFolder, 'traces.umt'), roiUMT);
            roiBefore = computeFileChecksum(fullfile(testCase.TempFolder, 'traces.umt'));

            applyRegistrationTformOnFolder(testCase.TempFolder, ...
                'RequireUserConfirmation', false, 'OpenQCFigure', false);

            registeredGreen = single(loadData(fullfile(testCase.TempFolder, 'green.dat')));
            after = loadData(fullfile(testCase.TempFolder, 'maps.umt'));
            validateUMTStruct(after, 'requireEventInfo', false);
            testCase.verifyEqual(single(after.data.main.value), registeredGreen(:, :, 1:2));
            testCase.verifyEqual(single(after.data.frame.value), registeredGreen(:, :, 1));
            testCase.verifyEqual(after.data.main.meta.FrameRateHz, 10);
            testCase.verifyEqual(computeFileChecksum(fullfile(testCase.TempFolder, 'traces.umt')), ...
                roiBefore, 'a roi UMT is not image data');
            testCase.verifyEmpty(dir(fullfile(testCase.TempFolder, '*.tmp')));
            testCase.verifyEqual(fx.Ny, size(after.data.main.value, 1));
        end

        function testRejectsImageUMTOfAnotherSize(testCase)
            bad = genUMTStruct(rand(16, 16, 'single'), 'kind', 'image', ...
                'entryName', 'small', 'dimNames', {'Y','X'});
            saveData(fullfile(testCase.TempFolder, 'small.umt'), bad);
            before = snapshotRegistrationFiles(testCase.TempFolder);

            testCase.verifyError(@() applyRegistrationTformOnFolder( ...
                testCase.TempFolder, 'RequireUserConfirmation', false), ...
                'ImageAlignmentTool:UmtSizeMismatch');
            verifyRegistrationSnapshot(testCase, testCase.TempFolder, before);
        end

        function testRegistersEventSplitDatFrameByFrame(testCase)
            % A Y-X-T-E file is registered like the continuous recording:
            % every frame of every trial gets the same transform, so trials
            % made of the green frames equal the registered green frames.
            fx = testCase.Fixture;
            green = single(loadData(fullfile(testCase.TempFolder, 'green.dat')));
            trials = cat(4, green(:, :, 1:2), green(:, :, 3:4));
            writeTestDat(fullfile(testCase.TempFolder, 'trials.dat'), trials, 10, ...
                'DimNames', {'Y','X','T','E'});

            applyRegistrationTformOnFolder(testCase.TempFolder, ...
                'RequireUserConfirmation', false, 'OpenQCFigure', false);

            info = loadMetaData(fullfile(testCase.TempFolder, 'trials.dat'));
            testCase.verifyEqual(cellstr(string(info.dimNames(:).')), {'Y','X','T','E'});
            testCase.verifyEqual(double(info.dimSizes(:).'), [fx.Ny fx.Nx 2 2]);
            after = single(loadData(fullfile(testCase.TempFolder, 'trials.dat')));
            registeredGreen = single(loadData(fullfile(testCase.TempFolder, 'green.dat')));
            testCase.verifyEqual(after, ...
                cat(4, registeredGreen(:, :, 1:2), registeredGreen(:, :, 3:4)));
        end

        function testNaNPixelsAreRestoredPerFrame(testCase)
            % A NaN region that exists in only one trial is zeroed for the
            % interpolation and restored (warped) in that frame only.
            green = single(loadData(fullfile(testCase.TempFolder, 'green.dat')));
            trials = cat(4, green(:, :, 1:2), green(:, :, 3:4));
            trials(1:3, :, 2, 2) = NaN;
            writeTestDat(fullfile(testCase.TempFolder, 'trials.dat'), trials, 10, ...
                'DimNames', {'Y','X','T','E'});

            applyRegistrationTformOnFolder(testCase.TempFolder, ...
                'RequireUserConfirmation', false, 'OpenQCFigure', false);

            after = single(loadData(fullfile(testCase.TempFolder, 'trials.dat')));
            testCase.verifyTrue(any(isnan(after(:, :, 2, 2)), 'all'));
            testCase.verifyFalse(any(isnan(after(:, :, 1, 2)), 'all'));
            testCase.verifyFalse(any(isnan(after(:, :, :, 1)), 'all'));
        end

        function testGuiAndPipelineApplyGiveIdenticalFolderState(testCase)
            % Brief section 1.7: one apply path. The same transform applied
            % through the pipeline wrapper and through the GUI backend call
            % gives identical data, ROIs, mask and view.
            fx = testCase.Fixture;
            roi = iMakeRectROI('square', [fx.Ny fx.Nx], 20:40, 15:35);
            saveROIFile(fullfile(testCase.TempFolder, 'myROI.roi'), ...
                createROIFile([fx.Ny fx.Nx], 'ROIs', roi));
            DataParams = loadDataParams(testCase.TempFolder);
            DataParams.mask.logical(1:10, :) = false;
            saveDataParams(testCase.TempFolder, DataParams);
            guiFolder = [testCase.TempFolder '_gui'];
            copyfile(testCase.TempFolder, guiFolder);
            testCase.addTeardown(@() rmdir(guiFolder, 's'));
            DataParams = loadDataParams(testCase.TempFolder);

            applyRegistrationTformOnFolder(testCase.TempFolder, ...
                'RequireUserConfirmation', false, 'OpenQCFigure', false);
            evalc(['applyImageAlignmentToFolder(guiFolder, DataParams.registration.estimate.tform, ' ...
                '''referenceImage'', DataParams.registration.reference.image)']);

            for name = {'green.dat', 'red.dat'}
                testCase.verifyEqual(loadData(fullfile(guiFolder, name{1})), ...
                    loadData(fullfile(testCase.TempFolder, name{1})), name{1});
            end
            pipelineROI = load(fullfile(testCase.TempFolder, 'myROI.roi'), 'ROIFile', '-mat');
            guiROI = load(fullfile(guiFolder, 'myROI.roi'), 'ROIFile', '-mat');
            testCase.verifyEqual(guiROI.ROIFile.ROIs.geometry.verticesXY_px, ...
                pipelineROI.ROIFile.ROIs.geometry.verticesXY_px);
            testCase.verifyEqual(guiROI.ROIFile.ROIs.mask, pipelineROI.ROIFile.ROIs.mask);
            testCase.verifyNotEqual(pipelineROI.ROIFile.ROIs.geometry.verticesXY_px, ...
                roi.geometry.verticesXY_px, 'the ROI must have been transformed');
            pipelineDP = loadDataParams(testCase.TempFolder);
            guiDP = loadDataParams(guiFolder);
            testCase.verifyEqual(guiDP.mask.logical, pipelineDP.mask.logical);
            testCase.verifyNotEqual(pipelineDP.mask.logical, DataParams.mask.logical);
            testCase.verifyEqual(guiDP.view, pipelineDP.view);
            testCase.verifyTrue(guiDP.registration.apply.isRegistered);
            testCase.verifyTrue(pipelineDP.registration.apply.isRegistered);
            testCase.verifyEqual(guiDP.registration.reference.image, ...
                pipelineDP.registration.reference.image);
            % Same reference image, same commitment, whatever the path.
            testCase.verifyEqual(guiDP.registration.apply.appliedReferenceHash, ...
                pipelineDP.registration.apply.appliedReferenceHash);
            % The GUI call carried no managed provenance: nothing stale is
            % carried over from the earlier estimate.
            testCase.verifyEqual(guiDP.registration.reference.source.kind, 'file');
            testCase.verifyEmpty(guiDP.registration.reference.source.resourceUUID);
            testCase.verifyEmpty(guiDP.registration.estimate.qcFigureFile);
            testCase.verifyEqual(pipelineDP.registration.reference.source.resourceUUID, ...
                fx.ImageReferenceUUID);
            % Both write an alignment report and no backup zip.
            for folder = {testCase.TempFolder, guiFolder}
                testCase.verifyNumElements(dir(fullfile(folder{1}, 'manualAlignmentReport_*.mat')), 1);
                testCase.verifyEmpty(dir(fullfile(folder{1}, '*.zip')));
            end
        end

        function testFailedCommitRestoresEveryOriginal(testCase)
            % A failure while installing the aligned files (here: red.dat
            % held open, so it cannot be moved on Windows) restores every
            % original and DataParams; no temporary or backup file remains.
            testCase.assumeTrue(ispc, 'Needs Windows file locking.');
            fx = testCase.Fixture;
            roi = iMakeRectROI('square', [fx.Ny fx.Nx], 20:40, 15:35);
            saveROIFile(fullfile(testCase.TempFolder, 'myROI.roi'), ...
                createROIFile([fx.Ny fx.Nx], 'ROIs', roi));
            before = snapshotRegistrationFiles(testCase.TempFolder);
            roiBefore = computeFileChecksum(fullfile(testCase.TempFolder, 'myROI.roi'));
            fid = fopen(fullfile(testCase.TempFolder, 'red.dat'), 'r');
            closer = onCleanup(@() fclose(fid));

            testCase.verifyError(@() applyRegistrationTformOnFolder(testCase.TempFolder, ...
                'RequireUserConfirmation', false, 'OpenQCFigure', false), ...
                'ImageAlignmentTool:CommitFailed');
            clear closer

            verifyRegistrationSnapshot(testCase, testCase.TempFolder, before);
            testCase.verifyEqual(computeFileChecksum(fullfile(testCase.TempFolder, 'myROI.roi')), ...
                roiBefore);
            testCase.verifyEmpty(dir(fullfile(testCase.TempFolder, '*.tmp')));
            testCase.verifyEmpty(dir(fullfile(testCase.TempFolder, '*.bak')));
            testCase.verifyEmpty(dir(fullfile(testCase.TempFolder, 'manualAlignmentReport_*')));
            testCase.verifyNumElements(dir(fullfile(testCase.TempFolder, ...
                'ImageAlignmentFailedReport_*.txt')), 1);
            testCase.verifyFalse(loadDataParams(testCase.TempFolder).registration.apply.isRegistered);

            % Nothing was committed, so the folder can still be aligned.
            applyRegistrationTformOnFolder(testCase.TempFolder, ...
                'RequireUserConfirmation', false, 'OpenQCFigure', false);
            testCase.verifyTrue(loadDataParams(testCase.TempFolder).registration.apply.isRegistered);
        end

        function testPipelineInfoUsesCanonicalDataParamsFile(testCase)
            info = applyRegistrationTformOnFolder('pipelineInfo');

            testCase.verifyFalse(any(strcmp( ...
                {info.parameters.name}, 'DataParamsFile')));
            testCase.verifySubstring(info.inputs(1).description, ...
                'DataParams.mat');
            testCase.verifyEqual({info.outputs.name}, ...
                {'rewrittenDatFiles', 'dataParams'});
            testCase.verifyEqual({info.outputs.defOutfilename}, ...
                {'*.dat', 'DataParams.mat'});
            testCase.verifyFalse(any([info.outputs.isData]));
            testCase.verifyFalse(any([info.outputs.returnsValue]));
        end
    end
end

function snapshot = snapshotRegistrationFiles(folder)
files = dir(fullfile(folder, '*.dat'));
names = string({files.name}).';
names(end+1,1) = "DataParams.mat";
names = sort(names);

checksums = strings(size(names));
for iFile = 1:numel(names)
    checksums(iFile) = string(computeFileChecksum( ...
        fullfile(folder, char(names(iFile)))));
end
snapshot = table(names, checksums, ...
    'VariableNames', {'FileName', 'Checksum'});
end

function verifyRegistrationSnapshot(testCase, folder, expected)
actual = snapshotRegistrationFiles(folder);
testCase.verifyEqual(actual, expected);
end

function installCancelQuestdlgMock(testCase)
mockFolder = fullfile(testCase.TempFolder, 'questdlg_mock');
mkdir(mockFolder);
mockFile = fullfile(mockFolder, 'questdlg.m');
fid = fopen(mockFile, 'w');
assert(fid ~= -1, 'Could not create questdlg test double.');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'function choice = questdlg(varargin)\n');
fprintf(fid, 'choice = ''Cancel'';\n');
fprintf(fid, 'end\n');
clear cleanup
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(mockFolder));
rehash
end

function writeLegacyDat(folder, baseName, data, datatype, dimNames)
%WRITELEGACYDAT Write a .dat plus a legacy metadata sidecar.
%
% loadMetaData recognizes a sidecar when it carries at least two of
% dim_names/datSize/datLength/Freq/Datatype, and legacy metadata takes
% precedence over AcqInfos.mat. That lets a test declare a datatype or a
% dimension layout that differs from the folder default.

sz = size(data);
fid = fopen(fullfile(folder, [baseName '.dat']), 'w');
fwrite(fid, cast(data, datatype), datatype);
fclose(fid);

meta = struct();
meta.dim_names = dimNames;
meta.datSize = sz(1:end-1);
meta.datLength = sz(end);
meta.Datatype = datatype;
meta.Freq = 10;
save(fullfile(folder, [baseName '.mat']), '-struct', 'meta');
end

function iRewriteWithHeader(filePath, dimSizes, rate)
% Put a v1 header in front of the file's existing single values
% (.dat header Phase 4c-1 guard test).
values = loadData(filePath);   % the input may already be headered (Phase 5a)
[~, base] = fileparts(filePath);
hdr = struct('dataClass', 'single', 'frameRateHz', rate, 'exposureMsec', NaN, ...
    'channelName', base, 'dimNames', {{'Y', 'X', 'T'}}, 'dimSizes', dimSizes, ...
    'writeComplete', true);
fid = fopen(filePath, 'w', 'ieee-le');
fwrite(fid, encodeDatHeader(hdr), 'uint8');
fwrite(fid, values, 'single');
fclose(fid);
end

function roi = iMakeRectROI(name, imageSizeYX, rowRange, colRange)
%IMAKERECTROI Build one schema-valid rectangle ROI (as TestGetDataFromROI).

mask = false(imageSizeYX);
mask(rowRange, colRange) = true;
x0 = min(colRange) - 0.5; x1 = max(colRange) + 0.5;
y0 = min(rowRange) - 0.5; y1 = max(rowRange) + 0.5;
vertices = [x0 y0; x1 y0; x1 y1; x0 y1];
stats = struct('computedOn', datetime('now'), 'NPixels', nnz(mask), ...
    'areaPx2', nnz(mask), 'areaMM2', [], 'centroidXY_px', [], 'centroidXY_mm', [], ...
    'distanceFromOrigin_px', [], 'distanceFromOrigin_mm', [], 'spatialMean', [], ...
    'spatialStd', [], 'spatialMedian', [], 'spatialMin', [], 'spatialMax', []);
roi = struct('name', name, 'type', 'polygon', 'DOC', datetime('now'), ...
    'modifiedOn', datetime('now'), 'color', [1 0 0], 'notes', '', ...
    'geometry', struct('polyshape', polyshape(vertices(:,1), vertices(:,2)), ...
        'verticesXY_px', vertices, 'ROIType', 'polygon', ...
        'ROIParameters', struct('ROIType', 'polygon', 'Position', vertices, ...
        'Vertices', vertices)), ...
    'mask', mask, 'stats', stats);
end
