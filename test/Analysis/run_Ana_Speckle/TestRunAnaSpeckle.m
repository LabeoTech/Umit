classdef TestRunAnaSpeckle < matlab.unittest.TestCase
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
            if isappdata(0, 'RunAnaSpeckleTestConfig')
                cfg = getappdata(0, 'RunAnaSpeckleTestConfig');
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
            deleteIfExists(fullfile(testCase.TempFolder, 'Flow_compute.dat'));
        end
    end

    methods (TestMethodTeardown)
        function cleanup(testCase)
            deleteIfExists(fullfile(testCase.TempFolder, 'Flow.dat'));
            deleteIfExists(fullfile(testCase.TempFolder, 'Flow_compute.dat'));
        end
    end

    methods (Test)
        function testPipelineInfo(testCase)
            info = run_Ana_Speckle('pipelineInfo');
            testCase.verifyEqual(info.name, 'run_Ana_Speckle');
            testCase.verifyEqual(info.outputs(1).defOutfilename, 'Flow.dat');
            testCase.verifyEqual(info.outputs(1).saveFileName, 'Flow.dat');

            dataInput = info.inputs(strcmp({info.inputs.name}, 'data'));
            testCase.verifyTrue(dataInput.supportsFile);
            testCase.verifyEqual(dataInput.dataMode, 'either');
            testCase.verifyFalse(any(strcmp({info.parameters.name}, 'SpeckleFileName')));
            for name = {'FrameRateHz', 'ExposureMsec'}
                testCase.verifyTrue(any(strcmp({info.sourceInfo.name}, name{1})), ...
                    sprintf('%s must be a sourceInfo input.', name{1}));
            end
        end

        function testArrayInputReturnsTheFlowArray(testCase)
            % An array runs in Standard mode and returns the flow array;
            % the metadata cannot come from a header, so they are explicit.
            datFile = fullfile(testCase.TempFolder, 'speckle.dat');
            md = loadMetaData(datFile);
            arr = loadData(datFile);

            out = run_Ana_Speckle(testCase.TempFolder, arr, ...
                'FrameRateHz', double(md.frameRateHz), ...
                'ExposureMsec', double(md.exposureMsec));
            fromFile = loadData(run_Ana_Speckle(testCase.TempFolder, 'speckle.dat'));

            testCase.verifyClass(out, 'single');
            testCase.verifyEqual(size(out), size(arr));
            testCase.verifyEqual(out, fromFile, 'RelTol', single(1e-3));
        end

        function testArrayWithoutMetadataErrors(testCase)
            arr = loadData(fullfile(testCase.TempFolder, 'speckle.dat'));

            testCase.verifyError( ...
                @() run_Ana_Speckle(testCase.TempFolder, arr), ...
                'Umitoolbox:Ana_Speckle:missingFrameRateHz');
        end

        function testFileOutsideSaveFolderIsRead(testCase)
            elsewhere = fullfile(testCase.TempFolder, 'elsewhere');
            mkdir(elsewhere);
            copyfile(fullfile(testCase.TempFolder, 'speckle.dat'), ...
                fullfile(elsewhere, 'speckle.dat'));

            out = run_Ana_Speckle(testCase.TempFolder, fullfile(elsewhere, 'speckle.dat'));

            testCase.verifyEqual(char(string(out)), fullfile(testCase.TempFolder, 'Flow.dat'));
            testCase.verifyTrue(isfile(out));
        end

        function testDatInputWritesAHeaderedYXTFile(testCase)
            out = run_Ana_Speckle(testCase.TempFolder, 'speckle.dat', 'bNormalize', false);

            testCase.verifyEqual(char(string(out)), fullfile(testCase.TempFolder, 'Flow.dat'));
            testCase.verifyTrue(isfile(out));
            testCase.verifyFalse(isfile(fullfile(testCase.TempFolder, 'Flow_compute.dat')));
            outInfo = loadMetaData(out);
            inInfo = loadMetaData(fullfile(testCase.TempFolder, 'speckle.dat'));
            testCase.verifyEqual(outInfo.dimNames, {'Y','X','T'});
            testCase.verifyEqual(outInfo.dimSizes, inInfo.dimSizes);
            testCase.verifyEqual(outInfo.dataClass, 'single');
            testCase.verifyEqual(outInfo.frameRateHz, inInfo.frameRateHz);
        end

        function testNormalizeDividesByTemporalMean(testCase)
            raw = loadData(run_Ana_Speckle(testCase.TempFolder, 'speckle.dat'));
            out = loadData(run_Ana_Speckle(testCase.TempFolder, 'speckle.dat', ...
                'bNormalize', true));

            expected = single(double(raw) ./ mean(double(raw), 3, 'omitnan'));
            testCase.verifyEqual(out, expected, 'RelTol', single(1e-4));
        end

        function testBareNameWithoutExtensionIsAccepted(testCase)
            out = run_Ana_Speckle(testCase.TempFolder, 'speckle');

            testCase.verifyTrue(isfile(out));
        end

        function testRejectsUnsupportedInputs(testCase)
            sv = testCase.TempFolder;

            testCase.verifyError(@() run_Ana_Speckle(sv, struct('a', 1)), ...
                'Umitoolbox:run_Ana_Speckle:UnsupportedInputType');
            testCase.verifyError(@() run_Ana_Speckle(sv, zeros(4, 4, 'single')), ...
                'Umitoolbox:run_Ana_Speckle:unsupportedLayout');
            testCase.verifyError(@() run_Ana_Speckle(sv, 'speckle.tif'), ...
                'Umitoolbox:run_Ana_Speckle:UnsupportedInputFile');
            testCase.verifyError(@() run_Ana_Speckle(sv, 'missing.dat'), ...
                'Umitoolbox:run_Ana_Speckle:FileNotFound');

            % Layouts other than Y-X-T are refused before anything is written.
            layouts = {{'Y','X','T','E'}, [6 5 4 2]; {'Y','X','E'}, [6 5 3]; {'Y','X'}, [6 5]};
            for k = 1:size(layouts, 1)
                writeTestDat(fullfile(sv, 'layout.dat'), rand([layouts{k, 2}, 1], 'single'), ...
                    10, 5, 'DimNames', layouts{k, 1});
                testCase.verifyError(@() run_Ana_Speckle(sv, 'layout.dat'), ...
                    'Umitoolbox:run_Ana_Speckle:unsupportedLayout');
            end
            testCase.verifyFalse(isfile(fullfile(sv, 'Flow.dat')));
        end
    end
end

function deleteIfExists(filePath)
if isfile(filePath)
    delete(filePath);
end
end
