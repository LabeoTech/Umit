classdef TestAnaSpeckle < matlab.unittest.TestCase
    properties
        ProjectRoot
        SampleDataFolder
        TempFolder
    end

    methods (TestClassSetup)
        function addProjectPathAndLocateFixture(testCase)
            thisFile = mfilename('fullpath');
            testFolder = fileparts(thisFile);
            projectRoot = extractBefore(testFolder, [filesep 'test']);
            if isempty(projectRoot)
                projectRoot = fileparts(fileparts(testFolder));
            end

            testCase.ProjectRoot = char(projectRoot);
            addpath(genpath(testCase.ProjectRoot));

            cfg = [];
            if isappdata(0, 'AnaSpeckleTestConfig')
                cfg = getappdata(0, 'AnaSpeckleTestConfig');
            end

            if ~isempty(cfg) && isfield(cfg, 'sampleDataFolder') && isfolder(cfg.sampleDataFolder)
                testCase.SampleDataFolder = char(string(cfg.sampleDataFolder));
            else
                testCase.SampleDataFolder = fullfile( ...
                    testCase.ProjectRoot, ...
                    'test', ...
                    'Analysis', ...
                    'TestingData_speckle');
            end

            testCase.assertTrue(isfile(fullfile(testCase.SampleDataFolder, 'speckle.dat')));
            testCase.assertTrue(isfile(fullfile(testCase.SampleDataFolder, 'AcqInfos.mat')));
        end
    end

    methods (TestMethodSetup)
        function createTempFolder(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture
            fx = testCase.applyFixture(TemporaryFolderFixture);
            testCase.TempFolder = fx.Folder;

            copyfile(fullfile(testCase.SampleDataFolder, 'speckle.dat'), ...
                fullfile(testCase.TempFolder, 'speckle.dat'));
            copyfile(fullfile(testCase.SampleDataFolder, 'AcqInfos.mat'), ...
                fullfile(testCase.TempFolder, 'AcqInfos.mat'));

            deleteIfExists(fullfile(testCase.TempFolder, 'Flow.dat'));
            deleteIfExists(fullfile(testCase.TempFolder, 'Flow_raw.dat'));
            deleteIfExists(fullfile(testCase.TempFolder, 'Flow_compute.dat'));
        end
    end

    methods (TestMethodTeardown)
        function cleanup(testCase)
            deleteIfExists(fullfile(testCase.TempFolder, 'Flow.dat'));
            deleteIfExists(fullfile(testCase.TempFolder, 'Flow_raw.dat'));
            deleteIfExists(fullfile(testCase.TempFolder, 'Flow_compute.dat'));
        end
    end

    methods (Test)
        function testPipelineInfo(testCase)
            info = Ana_Speckle('pipelineInfo');
            testCase.verifyEqual(info.name, 'Ana_Speckle');
            testCase.verifyEqual(info.outputs(1).defOutfilename, 'Flow.dat');
            testCase.verifyEqual(info.outputs(1).saveFileName, 'Flow.dat');
        end

        function testLoadMetaDataReturnsSpeckleExposure(testCase)
            % The speckle fixture is headered since .dat header Phase 5a:
            % its header carries the speckle exposure (the AcqInfos-bound
            % ExposureSpeckleMsec alias no longer applies).
            md = loadMetaData(fullfile(testCase.TempFolder, 'speckle.dat'));

            testCase.verifyEqual(md.format, 'header');
            testCase.verifyEqual(double(md.exposureMsec), 24);
        end

        function testStandardModeReturnsYXT(testCase)
            [out, meta] = iStandard(testCase, false);
            testCase.verifyClass(out, 'single');
            md = loadMetaData(fullfile(testCase.TempFolder, 'speckle.dat'));
            testCase.verifyEqual(size(out), [datAxisSize(md, 'Y') datAxisSize(md, 'X') datAxisSize(md, 'T')]);
            % .dat header Phase 4c-2b: without a written file, metaData is
            % the .dat Info schema the Flow.dat output would have.
            testCase.verifyEqual(meta.filePath, fullfile(testCase.TempFolder, 'Flow.dat'));
            testCase.verifyEqual(meta.format, 'header');
            testCase.verifyEqual(meta.dataOffset, 512);
            testCase.verifyEqual(meta.dataClass, 'single');
            testCase.verifyEqual(meta.dimNames, {'Y', 'X', 'T'});
            testCase.verifyEqual(meta.dimSizes, size(out));
            testCase.verifyEqual(meta.frameRateHz, double(md.frameRateHz));
            testCase.verifyEqual(meta.exposureMsec, double(md.exposureMsec));
            testCase.verifyEqual(meta.channelName, 'Flow');
            testCase.verifyFalse(isfield(meta, 'datLength'), 'no deprecated names');
        end

        function testLowRAMModeReturnsFile(testCase)
            [outFile, meta] = Ana_Speckle('speckle.dat', testCase.TempFolder, false);
            testCase.verifyTrue(ischar(outFile) || (isstring(outFile) && isscalar(outFile)));
            outFile = char(string(outFile));
            testCase.verifyEqual(outFile, fullfile(testCase.TempFolder, 'Flow.dat'));
            testCase.verifyTrue(isfile(outFile));
            md = loadMetaData(fullfile(testCase.TempFolder, 'speckle.dat'));
            outMd = loadMetaData(outFile);
            testCase.verifyEqual(datAxisSize(outMd, 'Y'), datAxisSize(md, 'Y'));
            testCase.verifyEqual(datAxisSize(outMd, 'X'), datAxisSize(md, 'X'));
            testCase.verifyEqual(datAxisSize(outMd, 'T'), datAxisSize(md, 'T'));
            % .dat header Phase 4c-2b: metaData is the written file's Info.
            testCase.verifyEqual(meta, outMd);
            testCase.verifyEqual(outMd.format, 'header');
            testCase.verifyEqual(outMd.exposureMsec, double(single(md.exposureMsec)));
            testCase.verifyEqual(outMd.frameRateHz, double(single(md.frameRateHz)));
            testCase.verifyEqual(outMd.channelName, 'Flow');
            outFileInfo = dir(outFile);
            testCase.verifyEqual(outFileInfo.bytes, outMd.dataOffset + ...
                double(datAxisSize(md, 'Y')) * double(datAxisSize(md, 'X')) * double(datAxisSize(md, 'T')) * 4);
        end

        function testLowRAMOverwritesDeclaredOutputOnRerun(testCase)
            % A pre-existing output must be replaced, not sidestepped under a
            % different name: the declared output is Flow.dat, so a re-run
            % has to keep writing there and must not leave a stale original
            % or scratch file behind (DFR-20260819-008 sibling).
            staleFile = fullfile(testCase.TempFolder, 'Flow.dat');
            fid = fopen(staleFile, 'w');
            fclose(fid);
            staleBytes = dir(staleFile).bytes;

            outFile = Ana_Speckle('speckle.dat', testCase.TempFolder, false);
            testCase.verifyEqual(char(string(outFile)), staleFile);
            testCase.verifyTrue(isfile(staleFile));
            testCase.verifyGreaterThan(dir(staleFile).bytes, staleBytes);
            testCase.verifyFalse( ...
                isfile(fullfile(testCase.TempFolder, 'Flow_compute.dat')));
        end

        function testNoOutputWritesDefaultFile(testCase)
            iStandard(testCase, false);
            flowFile = fullfile(testCase.TempFolder, 'Flow.dat');
            testCase.verifyTrue(isfile(flowFile));
            md = loadMetaData(flowFile);
            inputMd = loadMetaData(fullfile(testCase.TempFolder, 'speckle.dat'));
            testCase.verifyEqual(datAxisSize(md, 'T'), datAxisSize(inputMd, 'T'));
        end

        function testNormalizedOutputHasUnitTemporalMean(testCase)
            % Acceptance: output-level normalization must yield a per-pixel
            % temporal mean of ~1 when bNormalize is true.
            out = iStandard(testCase, true);
            pixelMean = mean(out, 3);
            testCase.verifyEqual(pixelMean, ones(size(pixelMean), 'single'), 'RelTol', single(1e-3));
        end

        function testStandardNormalizeMatchesUnnormalizedOwnMean(testCase)
            % The static-structure (MeanMap) correction must be identical
            % regardless of bNormalize: normalized output must equal the
            % unnormalized output divided by its own temporal mean.
            outRaw  = iStandard(testCase, false);
            outNorm = iStandard(testCase, true);
            expectedNorm = outRaw ./ mean(outRaw, 3);
            testCase.verifyEqual(outNorm, expectedNorm, 'RelTol', single(1e-3));
        end

        function testRAMSafeNormalizeMatchesUnnormalizedOwnMean(testCase)
            % Same equivalence as above, but for the RAMSafe streaming path.
            rawFile = Ana_Speckle('speckle.dat', testCase.TempFolder, false);
            rawDat = readFlowDat(char(string(rawFile)));
            movefile(char(string(rawFile)), fullfile(testCase.TempFolder, 'Flow_raw.dat'));

            normFile = Ana_Speckle('speckle.dat', testCase.TempFolder, true);
            normDat = readFlowDat(char(string(normFile)));

            expectedNorm = rawDat ./ mean(rawDat, 3);
            testCase.verifyEqual(normDat, expectedNorm, 'RelTol', single(1e-3));

            deleteIfExists(fullfile(testCase.TempFolder, 'Flow_raw.dat'));
        end

        function testArrayNeedsFrameRateAndExposure(testCase)
            arr = loadData(fullfile(testCase.TempFolder, 'speckle.dat'));

            testCase.verifyError( ...
                @() Ana_Speckle(arr, testCase.TempFolder, false, 'ExposureMsec', 24), ...
                'Umitoolbox:Ana_Speckle:missingFrameRateHz');
            testCase.verifyError( ...
                @() Ana_Speckle(arr, testCase.TempFolder, false, 'FrameRateHz', 10), ...
                'Umitoolbox:Ana_Speckle:missingExposureMsec');
        end

        function testExplicitExposureOverridesTheHeader(testCase)
            % An explicit value wins over the file's header, with a warning.
            md = loadMetaData(fullfile(testCase.TempFolder, 'speckle.dat'));
            newExposure = double(md.exposureMsec) * 2;

            testCase.verifyWarning( ...
                @() Ana_Speckle('speckle.dat', testCase.TempFolder, false, ...
                'ExposureMsec', newExposure), ...
                'Umitoolbox:Ana_Speckle:sourceInfoConflict');
            outMd = loadMetaData(fullfile(testCase.TempFolder, 'Flow.dat'));
            testCase.verifyEqual(outMd.exposureMsec, newExposure);
        end

        function testRejectsUnsupportedInputs(testCase)
            testCase.verifyError( ...
                @() Ana_Speckle(zeros(4, 4, 'single'), testCase.TempFolder, false), ...
                'Ana_Speckle:UnsupportedLayout');
            testCase.verifyError( ...
                @() Ana_Speckle({1}, testCase.TempFolder, false), ...
                'Ana_Speckle:UnsupportedInputType');
            testCase.verifyError( ...
                @() Ana_Speckle('missing.dat', testCase.TempFolder, false), ...
                'Ana_Speckle:InputFileNotFound');
        end

        function testStandardModeMatchesRAMSafeMode(testCase)
            % Standard mode and RAMSafe mode must agree for the same input
            % and the same bNormalize value (both true and false).
            for bNorm = [false true]
                stdOut = iStandard(testCase, bNorm);
                ramFile = Ana_Speckle('speckle.dat', testCase.TempFolder, bNorm);
                ramOut = readFlowDat(char(string(ramFile)));
                testCase.verifyEqual(stdOut, ramOut, 'RelTol', single(1e-3));
            end
        end
    end
end

function deleteIfExists(filePath)
if isfile(filePath)
    delete(filePath);
end
end

function dat = readFlowDat(filePath)
% Flow.dat is headered since .dat header Phase 4c-2b: read with loadData.
dat = single(loadData(filePath));
end

function varargout = iStandard(testCase, bNormalize)
%ISTANDARD Ana_Speckle in Standard mode: the fixture recording as an array.
%   An array carries no metadata, so the header values are passed explicitly.
datFile = fullfile(testCase.TempFolder, 'speckle.dat');
md = loadMetaData(datFile);
[varargout{1:nargout}] = Ana_Speckle(loadData(datFile), testCase.TempFolder, bNormalize, ...
    'FrameRateHz', double(md.frameRateHz), 'ExposureMsec', double(md.exposureMsec));
end
