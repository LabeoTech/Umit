classdef TestApplyRegistrationTformOnFolder < matlab.unittest.TestCase
    %TESTAPPLYREGISTRATIONTFORMONFOLDER Unit tests for applyRegistrationTformOnFolder.

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
            testCase.verifyTrue(S.DataParams.registration.isRegistered);
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
            testCase.verifyTrue(DP.registration.isRegistered);
            testCase.verifyEqual(char(string(DP.registration.appliedBy)), 'applyRegistrationTformOnFolder');
            testCase.verifyEqual(char(string(DP.registration.confirmationMode)), 'forced_no_confirmation');
            testCase.verifyEqual( ...
                DP.registration.imageReferenceUUID, ...
                testCase.Fixture.ImageReferenceUUID);
            testCase.verifyEqual( ...
                DP.registration.imageReferenceChecksum, ...
                testCase.Fixture.ManagedReference.checksum);
        end

        function testAlreadyRegisteredBlockedByDefault(testCase)
            applyRegistrationTformOnFolder(testCase.TempFolder, 'RequireUserConfirmation', false);

            testCase.verifyError(@() applyRegistrationTformOnFolder(testCase.TempFolder, ...
                'RequireUserConfirmation', false), ...
                'Umitoolbox:applyRegistrationTformOnFolder:AlreadyRegistered');
        end

        function testAllowReapplyOverride(testCase)
            applyRegistrationTformOnFolder(testCase.TempFolder, 'RequireUserConfirmation', false);
            applyRegistrationTformOnFolder(testCase.TempFolder, ...
                'RequireUserConfirmation', false, 'AllowReapply', true);

            DataParams = loadDataParams(testCase.TempFolder);
            testCase.verifyTrue(DataParams.registration.isRegistered);
        end

        function testMissingTformErrors(testCase)
            S = load(fullfile(testCase.TempFolder, 'DataParams.mat'), 'DataParams');
            S.DataParams.registration.tform = [];
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

        function testRejectsReferenceSizeMismatch(testCase)
            % The transform is estimated against the reference frame size, so
            % applying it to differently sized data would silently shift and
            % scale by the wrong amount.
            DataParams = loadDataParams(testCase.TempFolder);
            smallerYX = [32 32];
            DataParams.view.imageSizeYX = smallerYX;
            DataParams.registration.referenceImage = zeros(smallerYX, 'single');
            DataParams.mask.logical = [];
            saveDataParams(testCase.TempFolder, DataParams);
            before = snapshotRegistrationFiles(testCase.TempFolder);

            testCase.verifyError(@() applyRegistrationTformOnFolder( ...
                testCase.TempFolder, 'RequireUserConfirmation', false), ...
                'Umitoolbox:applyRegistrationTformOnFolder:ReferenceSizeMismatch');

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
                'Umitoolbox:applyRegistrationTformOnFolder:UnsupportedLayout');

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
            testCase.verifyEmpty(dir(fullfile(testCase.TempFolder, '*_registering.umt')));
            testCase.verifyEqual(fx.Ny, size(after.data.main.value, 1));
        end

        function testRejectsImageUMTOfAnotherSize(testCase)
            bad = genUMTStruct(rand(16, 16, 'single'), 'kind', 'image', ...
                'entryName', 'small', 'dimNames', {'Y','X'});
            saveData(fullfile(testCase.TempFolder, 'small.umt'), bad);
            before = snapshotRegistrationFiles(testCase.TempFolder);

            testCase.verifyError(@() applyRegistrationTformOnFolder( ...
                testCase.TempFolder, 'RequireUserConfirmation', false), ...
                'Umitoolbox:applyRegistrationTformOnFolder:ReferenceSizeMismatch');
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
