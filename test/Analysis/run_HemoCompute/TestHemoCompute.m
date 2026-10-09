classdef TestHemoCompute < matlab.unittest.TestCase
    %TESTHEMOCOMPUTE Unit tests for HemoCompute.
    %
    %   These tests use small synthetic AcqInfos-driven datasets so the
    %   expected timelines are explicit. They cover equal-length imported
    %   channels and repeated-illumination-like datasets where one imported
    %   channel has a higher frame rate and more frames than the others.

    properties
        ProjectRoot
        TempFolder
        AcqInfo
        RigRoot
    end

    methods (TestMethodSetup)
        function createTempFolder(testCase)
            thisFile = mfilename('fullpath');
            testFolder = fileparts(thisFile);
            projectRoot = extractBefore(testFolder, [filesep 'test']);
            if isempty(projectRoot)
                projectRoot = fileparts(fileparts(testFolder));
            end

            testCase.ProjectRoot = char(projectRoot);
            addpath(genpath(testCase.ProjectRoot));

            testCase.TempFolder = fullfile(tempdir, ...
                ['TestHemoCompute_' char(java.util.UUID.randomUUID)]);
            mkdir(testCase.TempFolder);

            iCreateSyntheticHemoDataset(testCase.TempFolder, false);
            rigStore = iCreateConfiguredRig();
            testCase.RigRoot = rigStore.RigRoot;
            iBindDatasetToRig(testCase.TempFolder, rigStore);
            s = load(fullfile(testCase.TempFolder, 'AcqInfos.mat'), 'AcqInfoStream');
            testCase.AcqInfo = s.AcqInfoStream;
        end
    end

    methods (TestMethodTeardown)
        function removeTempFolder(testCase)
            if ~isempty(testCase.TempFolder) && isfolder(testCase.TempFolder)
                try %#ok<TRYNC>
                    rmdir(testCase.TempFolder, 's');
                end
            end
            if ~isempty(testCase.RigRoot) && isfolder(testCase.RigRoot)
                rmdir(testCase.RigRoot, 's');
            end
        end
    end

    methods (Test)
        function testTwoChannelStandard(testCase)
            [HbO, HbR] = HemoCompute(testCase.TempFolder, testCase.TempFolder, ...
                'none', {'red','green'}, true);

            testCase.verifyTrue(isnumeric(HbO));
            testCase.verifyTrue(isnumeric(HbR));
            testCase.verifySize(HbO, [testCase.AcqInfo.Height, testCase.AcqInfo.Width, testCase.AcqInfo.Length]);
            testCase.verifySize(HbR, [testCase.AcqInfo.Height, testCase.AcqInfo.Width, testCase.AcqInfo.Length]);
            testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, 'HbO.dat')));
            testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, 'HbR.dat')));
        end

        function testExplicitResolvedOpticalInfo(testCase)
            rigStore = UMITRigStore.open(testCase.AcqInfo.rigUUID);
            opticalInfo = rigStore.resolveOpticalConfiguration( ...
                testCase.AcqInfo, {'red','green'}, 'none', testCase.TempFolder);

            [HbO, HbR] = HemoCompute(testCase.TempFolder, '', ...
                'none', {'red','green'}, true, 'OpticalInfo', opticalInfo);

            testCase.verifyTrue(isnumeric(HbO));
            testCase.verifyTrue(isnumeric(HbR));
            testCase.verifySize(HbO, ...
                [testCase.AcqInfo.Height, testCase.AcqInfo.Width, testCase.AcqInfo.Length]);
        end

        function testRejectsMismatchedOpticalInfoActiveRows(testCase)
            rigStore = UMITRigStore.open(testCase.AcqInfo.rigUUID);
            opticalInfo = rigStore.resolveOpticalConfiguration( ...
                testCase.AcqInfo, {'red','green'}, 'none', testCase.TempFolder);
            opticalInfo.activeRows = [true; false; true];

            testCase.verifyError(@() HemoCompute( ...
                testCase.TempFolder, testCase.TempFolder, ...
                'none', {'red','green'}, true, 'OpticalInfo', opticalInfo), ...
                'Umitoolbox:HemoCompute:OpticalInfoChannelMismatch');
            testCase.verifyFalse(isfile(fullfile(testCase.TempFolder, 'HbO.dat')));
            testCase.verifyFalse(isfile(fullfile(testCase.TempFolder, 'HbR.dat')));
        end

        function testLegacyOpticsWhenRigStoreIsOffPath(testCase)
            rigPath = fullfile(testCase.ProjectRoot, 'RigManagement');
            testCase.assertNotEmpty(which('UMITRigStore'));
            rmpath(rigPath);
            pathCleanup = onCleanup(@() addpath(rigPath));
            testCase.assertEqual(exist('UMITRigStore', 'class'), 0);

            [HbO, HbR] = HemoCompute(testCase.TempFolder, '', ...
                'none', {'red','green'}, true);

            testCase.verifyTrue(isnumeric(HbO));
            testCase.verifyTrue(isnumeric(HbR));
            testCase.verifySize(HbO, ...
                [testCase.AcqInfo.Height, testCase.AcqInfo.Width, testCase.AcqInfo.Length]);
        end

        function testThreeChannelStandard(testCase)
            [HbO, HbR] = HemoCompute(testCase.TempFolder, testCase.TempFolder, ...
                'none', {'red','green','yellow'}, true);

            testCase.verifyTrue(isnumeric(HbO));
            testCase.verifyTrue(isnumeric(HbR));
            testCase.verifySize(HbO, [testCase.AcqInfo.Height, testCase.AcqInfo.Width, testCase.AcqInfo.Length]);
            testCase.verifySize(HbR, [testCase.AcqInfo.Height, testCase.AcqInfo.Width, testCase.AcqInfo.Length]);
        end

        function testAmberMapsToYellow(testCase)
            [HbO, HbR] = HemoCompute(testCase.TempFolder, testCase.TempFolder, ...
                'none', {'red','amber'}, true);

            testCase.verifyTrue(isnumeric(HbO));
            testCase.verifyTrue(isnumeric(HbR));
        end

        function testRepeatedIlluminationDownsamplesToLowestFrequency(testCase)
            repeatedFolder = fullfile(testCase.TempFolder, 'repeatedTimeline');
            mkdir(repeatedFolder);
            iCreateSyntheticHemoDataset(repeatedFolder, true);
            iBindDatasetToRig(repeatedFolder, ...
                UMITRigStore.openByRigID(iCurrentRigID(testCase.TempFolder)));

            [HbO, HbR] = HemoCompute(repeatedFolder, repeatedFolder, ...
                'none', {'red','green'}, true);

            infoGreen = loadMetaData(fullfile(repeatedFolder, 'green.dat'));
            infoRed = loadMetaData(fullfile(repeatedFolder, 'red.dat'));
            infoHbO = loadMetaData(fullfile(repeatedFolder, 'HbO.dat'));

            testCase.verifyGreaterThan(infoRed.FrameRateHz, infoGreen.FrameRateHz);
            testCase.verifyGreaterThan(infoRed.Length, infoGreen.Length);
            testCase.verifySize(HbO, [infoGreen.Height, infoGreen.Width, infoGreen.Length]);
            testCase.verifySize(HbR, [infoGreen.Height, infoGreen.Width, infoGreen.Length]);
            testCase.verifyEqual(infoHbO.Length, infoGreen.Length);
            testCase.verifyEqual(infoHbO.FrameRateHz, infoGreen.FrameRateHz);
        end

        function testRepeatedIlluminationRAMSafeMode(testCase)
            repeatedFolder = fullfile(testCase.TempFolder, 'repeatedTimelineRAMSafe');
            mkdir(repeatedFolder);
            iCreateSyntheticHemoDataset(repeatedFolder, true);
            iBindDatasetToRig(repeatedFolder, ...
                UMITRigStore.openByRigID(iCurrentRigID(testCase.TempFolder)));

            [HbO, HbR] = HemoCompute(repeatedFolder, repeatedFolder, ...
                'none', {'red','green'}, true, ...
                'RAMSafeMode', true);

            testCase.verifyEqual(HbO, 'HbO.dat');
            testCase.verifyEqual(HbR, 'HbR.dat');

            out = loadData(fullfile(repeatedFolder, 'HbO.dat'));
            infoGreen = loadMetaData(fullfile(repeatedFolder, 'green.dat'));
            testCase.verifySize(out, [infoGreen.Height, infoGreen.Width, infoGreen.Length]);
        end

        function testInvalidFilterSet(testCase)
            testCase.verifyError(@() HemoCompute(testCase.TempFolder, testCase.TempFolder, ...
                'badFilter', {'red','green'}, true), ...
                'Umitoolbox:UMITRigStore:filterSetNotFound');
        end

        function testInvalidIlluminationSelection(testCase)
            testCase.verifyError(@() HemoCompute(testCase.TempFolder, testCase.TempFolder, ...
                'none', {'red'}, true), ...
                'Umitoolbox:HemoCompute:InvalidIllumination');
        end

        function testMissingChannelFile(testCase)
            delete(fullfile(testCase.TempFolder, 'green.dat'));

            testCase.verifyError(@() HemoCompute(testCase.TempFolder, testCase.TempFolder, ...
                'none', {'red','green'}, true), ...
                'Umitoolbox:UMITRigStore:missingChannelFile');
        end

        function testRAMSafeMode(testCase)
            [HbO, HbR] = HemoCompute(testCase.TempFolder, testCase.TempFolder, ...
                'none', {'red','green'}, true, ...
                'RAMSafeMode', true);

            testCase.verifyEqual(HbO, 'HbO.dat');
            testCase.verifyEqual(HbR, 'HbR.dat');
            testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, 'HbO.dat')));
            testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, 'HbR.dat')));
        end

        function testRAMExhaustionFallsBackToNonEmptyChunking(testCase)
            % TASK 5.1 regression guard: calculateMaxChunkSize used to
            % return [] when no usable RAM remained, which made every
            % chunked caller's `for c = 1:nChunks` loop execute zero
            % iterations and silently leave a zero-filled output. Simulate
            % RAM exhaustion and verify HemoCompute still produces real,
            % non-zero-filled output.
            testCase.assumeEqual(computer, 'PCWIN64', ...
                'RAM-exhaustion simulation relies on shadowing the PCWIN64 memory() built-in.');

            mocksFolder = fullfile(testCase.ProjectRoot, ...
                'test', 'subFunc', 'calculateMaxChunkSize', 'mocks');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(mocksFolder));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_MODE', 'fixed'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_TOTAL', num2str(1e9, '%d')));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_AVAILABLE', '0'));

            lastwarn('');
            [HbO, HbR] = HemoCompute(testCase.TempFolder, testCase.TempFolder, ...
                'none', {'red','green'}, true);
            [~, warnId] = lastwarn();

            testCase.verifyEqual(warnId, 'Umitoolbox:calculateMaxChunkSize:NoUsableRAM');
            testCase.verifyTrue(isnumeric(HbO));
            testCase.verifyTrue(isnumeric(HbR));
            testCase.verifyGreaterThan(nnz(HbO), 0);
            testCase.verifyGreaterThan(nnz(HbR), 0);
        end

        function testClampFloorIsChunkCountInvariant(testCase)
            % F-8 regression: min(Red/Green(:)) used to be a chunk-local
            % reduction, so the clamp floor (and therefore the normalized
            % output) depended on chunkX, which itself depends on
            % available RAM at call time via calculateMaxChunkSize. Force
            % two different chunk counts (1 vs. 8, for an 8-column test
            % image) via the memory() mock and confirm identical output.
            testCase.assumeEqual(computer, 'PCWIN64', ...
                'Chunk-count forcing relies on shadowing the PCWIN64 memory() built-in.');

            mocksFolder = fullfile(testCase.ProjectRoot, ...
                'test', 'subFunc', 'calculateMaxChunkSize', 'mocks');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(mocksFolder));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_MODE', 'fixed'));

            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_TOTAL', num2str(1e12, '%d')));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_AVAILABLE', num2str(1e12, '%d')));
            [HbOsingleChunk, HbRsingleChunk] = HemoCompute(testCase.TempFolder, '', ...
                'none', {'red','green'}, true);

            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_TOTAL', num2str(1e9, '%d')));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_AVAILABLE', num2str(100005000, '%d')));
            [HbOmanyChunks, HbRmanyChunks] = HemoCompute(testCase.TempFolder, '', ...
                'none', {'red','green'}, true);

            testCase.verifyEqual(HbOmanyChunks, HbOsingleChunk);
            testCase.verifyEqual(HbRmanyChunks, HbRsingleChunk);
        end

        function testRAMSafeRerunOverwritesDeclaredOutputs(testCase)
            % A pre-existing output must be replaced, not sidestepped under a
            % different name: the declared outputs are HbO.dat/HbR.dat, so a
            % re-run has to keep writing there and must not leave a stale
            % original or scratch file behind (DFR-20260819-008).
            hboFile = fullfile(testCase.TempFolder, 'HbO.dat');
            hbrFile = fullfile(testCase.TempFolder, 'HbR.dat');
            fid = fopen(hboFile, 'w'); fclose(fid);
            fid = fopen(hbrFile, 'w'); fclose(fid);
            hboStaleBytes = dir(hboFile).bytes;
            hbrStaleBytes = dir(hbrFile).bytes;

            [HbO, HbR] = HemoCompute(testCase.TempFolder, testCase.TempFolder, ...
                'none', {'red','green'}, true, ...
                'RAMSafeMode', true);

            testCase.verifyEqual(HbO, 'HbO.dat');
            testCase.verifyEqual(HbR, 'HbR.dat');
            testCase.verifyTrue(isfile(hboFile));
            testCase.verifyTrue(isfile(hbrFile));
            testCase.verifyGreaterThan(dir(hboFile).bytes, hboStaleBytes);
            testCase.verifyGreaterThan(dir(hbrFile).bytes, hbrStaleBytes);
            testCase.verifyFalse(isfile(fullfile(testCase.TempFolder, 'HbO_preallocData.dat')));
            testCase.verifyFalse(isfile(fullfile(testCase.TempFolder, 'HbR_preallocData.dat')));
            testCase.verifyFalse(isfile(fullfile(testCase.TempFolder, 'HbO_writing.dat')));
            testCase.verifyFalse(isfile(fullfile(testCase.TempFolder, 'HbR_writing.dat')));
        end

        function testSaveDisabledStandardMode(testCase)
            [HbO, HbR] = HemoCompute(testCase.TempFolder, '', ...
                'none', {'red','green'}, true);

            testCase.verifyTrue(isnumeric(HbO));
            testCase.verifyTrue(isnumeric(HbR));
            testCase.verifyFalse(isfile(fullfile(testCase.TempFolder, 'HbO.dat')));
            testCase.verifyFalse(isfile(fullfile(testCase.TempFolder, 'HbR.dat')));
        end

        function testPipelineInfo(testCase)
            info = HemoCompute('pipelineInfo');
            testCase.verifyEqual(info.name, 'HemoCompute');
            testCase.verifyFalse(info.legacyOpts);
            testCase.verifyEqual(info.outputs(1).type, {'ImageTimeSeries'});
            testCase.verifyEqual(info.outputs(2).type, {'ImageTimeSeries'});
            testCase.verifyEqual(info.outputs(1).defOutfilename, 'HbO.dat');
            testCase.verifyEqual(info.outputs(2).defOutfilename, 'HbR.dat');
        end
    end
end

function store = iCreateConfiguredRig()
suffix = strrep(char(java.util.UUID.randomUUID()), '-', '');
cameras = struct( ...
    'index', 1, 'displayName', 'Test Camera', 'manufacturer', '', ...
    'model', 'D1024', 'serialNumber', '', 'spectrumID', 'PF1024');
illuminations = struct( ...
    'name', {'red','green','yellow'}, ...
    'displayName', {'Red','Green','Yellow'}, ...
    'manufacturer', {'','',''}, 'model', {'','',''}, ...
    'spectrumID', {'LED_632nm','LED_521nm','LED_593nm'});
store = UMITRigStore.create(struct( ...
    'rigID', ['HemoCore_' suffix(1:10)], ...
    'cameras', cameras, 'illuminations', illuminations));
end

function iBindDatasetToRig(folderPath, store)
loaded = load(fullfile(folderPath, 'AcqInfos.mat'), 'AcqInfoStream');
AcqInfoStream = loaded.AcqInfoStream;
rigInfo = store.getRigInfo();
AcqInfoStream.rigUUID = rigInfo.uuid;
AcqInfoStream.rigID = rigInfo.rigID;
save(fullfile(folderPath, 'AcqInfos.mat'), 'AcqInfoStream', '-mat');
end

function rigID = iCurrentRigID(folderPath)
loaded = load(fullfile(folderPath, 'AcqInfos.mat'), 'AcqInfoStream');
rigID = loaded.AcqInfoStream.rigID;
end

function iCreateSyntheticHemoDataset(folderPath, repeatedRed)
%ICREATESYNTHETICHEMODATASET Create normalized intrinsic-channel .dat files.

Ny = 8;
Nx = 8;
baseLength = 12;
baseFreq = 10;

if repeatedRed
    redLength = baseLength * 2;
    redFreq = baseFreq * 2;
else
    redLength = baseLength;
    redFreq = baseFreq;
end

AcqInfoStream = struct();
AcqInfoStream.Height = Ny;
AcqInfoStream.Width = Nx;
AcqInfoStream.Length = baseLength;
AcqInfoStream.FrameRateHz = baseFreq;
AcqInfoStream.Camera_Model = 'D1024';
AcqInfoStream.MultiCam = false;
AcqInfoStream.Datatype = 'single';
AcqInfoStream.ExposureMsec = 10;
AcqInfoStream.ImportedChannels = struct( ...
    'DatFile', {}, ...
    'Length', {}, ...
    'FrameRateHz', {}, ...
    'ExposureMsec', {}, ...
    'CamIdx', {});

channels = { ...
    'red.dat', redLength, redFreq, 0.010; ...
    'green.dat', baseLength, baseFreq, 0.015; ...
    'yellow.dat', baseLength, baseFreq, 0.020};

for iChan = 1:size(channels, 1)
    datFile = channels{iChan, 1};
    nFrames = channels{iChan, 2};
    frameRate = channels{iChan, 3};
    amp = channels{iChan, 4};

    data = iMakeNormalizedMovie(Ny, Nx, nFrames, amp);
    iWriteSingleDat(fullfile(folderPath, datFile), data, frameRate, 10);

    AcqInfoStream.ImportedChannels(end+1).DatFile = datFile;
    AcqInfoStream.ImportedChannels(end).Length = nFrames;
    AcqInfoStream.ImportedChannels(end).FrameRateHz = frameRate;
    AcqInfoStream.ImportedChannels(end).ExposureMsec = 10;
    AcqInfoStream.ImportedChannels(end).CamIdx = 1;
end

save(fullfile(folderPath, 'AcqInfos.mat'), 'AcqInfoStream');
end

function data = iMakeNormalizedMovie(Ny, Nx, Nt, amp)
baseMap = single(reshape(linspace(0, 1, Ny * Nx), Ny, Nx));
t = reshape(single(1:Nt), 1, 1, []);
data = single(1 + amp * sin(2*pi*t/max(Nt, 1)) + 0.001 * baseMap);
end

function iWriteSingleDat(filePath, data, frameRateHz, exposureMsec)
% Headered input (.dat header Phase 5a).
writeTestDat(filePath, single(data), frameRateHz, exposureMsec);
end
