classdef TestRunSpeckleMapping < matlab.unittest.TestCase
    %TESTRUNSPECKLEMAPPING Unit tests for run_SpeckleMapping.
    %
    % Uses the real fixture folder:
    %   Analysis/TestingData_speckle

    properties
        FixtureFolder = ''
        TempFolder = ''
    end

    methods (TestMethodSetup)
        function setup(testCase)
            thisFile = which(class(testCase));
            testRoot = fileparts(thisFile);
            analysisRoot = fileparts(testRoot);
            testCase.FixtureFolder = fullfile(analysisRoot, 'TestingData_speckle');

            testCase.assumeTrue(isfolder(testCase.FixtureFolder), ...
                sprintf('Fixture folder not found: %s', testCase.FixtureFolder));

            testCase.TempFolder = fullfile(tempdir, ...
                ['TestRunSpeckleMapping_' char(java.util.UUID.randomUUID)]);
            mkdir(testCase.TempFolder);

            copyfile(fullfile(testCase.FixtureFolder, 'speckle.dat'), ...
                fullfile(testCase.TempFolder, 'speckle.dat'));
            copyfile(fullfile(testCase.FixtureFolder, 'AcqInfos.mat'), ...
                fullfile(testCase.TempFolder, 'AcqInfos.mat'));
        end
    end

    methods (TestMethodTeardown)
        function teardown(testCase)
            if ~isempty(testCase.TempFolder) && isfolder(testCase.TempFolder)
                rmdir(testCase.TempFolder, 's');
            end
        end
    end

    methods (Test)
        function testPipelineInfo(testCase)
            info = run_SpeckleMapping('pipelineInfo');

            testCase.verifyTrue(isstruct(info) && isscalar(info));
            testCase.verifyTrue(isfield(info, 'inputs'));
            testCase.verifyTrue(isfield(info, 'outputs'));

            dataInput = info.inputs(strcmp({info.inputs.name}, 'data'));
            testCase.verifyTrue(dataInput.supportsFile);
            testCase.verifyEqual(dataInput.dataMode, 'either');
            testCase.verifyFalse(any(strcmp({info.parameters.name}, 'channel')), ...
                'The channel parameter was removed: the input is "data".');
        end

        function testFileInputLowRAMMode(testCase)
            out = run_SpeckleMapping( ...
                testCase.TempFolder, ...
                'speckle.dat', ...
                'sType', 'Spatial', ...
                'bSaveMap', false, ...
                'bLogScale', false);

            validateUMTStruct(out, 'requireEventInfo', false);
            testCase.verifyEqual(char(string(out.kind)), 'image');
            testCase.verifyTrue(isfield(out.data, 'SpeckleMap'));
        end

        function testSaveMapOption(testCase)
            tifFile = fullfile(testCase.TempFolder, 'std_speckle.tiff');
            if isfile(tifFile)
                delete(tifFile);
            end

            out = run_SpeckleMapping( ...
                testCase.TempFolder, ...
                'speckle.dat', ...
                'sType', 'Temporal', ...
                'bSaveMap', true, ...
                'bLogScale', false);

            validateUMTStruct(out, 'requireEventInfo', false);
            testCase.verifyTrue(isfile(tifFile));
        end

        function testArrayInputMatchesFileInput(testCase)
            % An array runs in Standard mode, a filename in Low-RAM mode;
            % both give the same map.
            arr = loadData(fullfile(testCase.TempFolder, 'speckle.dat'));

            for sType = ["Spatial", "Temporal"]
                fromArray = run_SpeckleMapping(testCase.TempFolder, arr, ...
                    'sType', char(sType), 'bLogScale', true);
                fromFile = run_SpeckleMapping(testCase.TempFolder, 'speckle.dat', ...
                    'sType', char(sType), 'bLogScale', true);

                validateUMTStruct(fromArray, 'requireEventInfo', false);
                testCase.verifyEqual(fromArray.data.SpeckleMap.value, ...
                    fromFile.data.SpeckleMap.value, 'AbsTol', single(1e-5), ...
                    sprintf('Array and file maps differ for sType=%s.', sType));
            end
        end

        function testFileNameIsNotRestrictedToChannelNames(testCase)
            % Any .dat name works: the file is the input, not a channel key.
            copyfile(fullfile(testCase.TempFolder, 'speckle.dat'), ...
                fullfile(testCase.TempFolder, 'other.dat'));

            out = run_SpeckleMapping(testCase.TempFolder, 'other.dat', ...
                'sType', 'Spatial', 'bLogScale', false);
            ref = run_SpeckleMapping(testCase.TempFolder, 'speckle.dat', ...
                'sType', 'Spatial', 'bLogScale', false);

            testCase.verifyEqual(out.data.SpeckleMap.value, ref.data.SpeckleMap.value);
        end

        function testFileOutsideSaveFolderIsRead(testCase)
            sv = testCase.TempFolder;
            elsewhere = fullfile(sv, 'elsewhere');
            mkdir(elsewhere);
            copyfile(fullfile(sv, 'speckle.dat'), fullfile(elsewhere, 'speckle.dat'));

            out = run_SpeckleMapping(sv, fullfile(elsewhere, 'speckle.dat'), ...
                'sType', 'Temporal');
            ref = run_SpeckleMapping(sv, 'speckle.dat', 'sType', 'Temporal');

            testCase.verifyEqual(out.data.SpeckleMap.value, ref.data.SpeckleMap.value);
        end

        function testRejectsUnsupportedInputs(testCase)
            sv = testCase.TempFolder;

            % Inputs that are neither a filename nor a numeric array.
            testCase.verifyError(@() run_SpeckleMapping(sv, struct('a', 1)), ...
                'Umitoolbox:run_SpeckleMapping:UnsupportedInputType');
            testCase.verifyError(@() run_SpeckleMapping(sv, 'speckle.tif'), ...
                'Umitoolbox:run_SpeckleMapping:UnsupportedInputFile');

            % Arrays other than Y-X-T.
            testCase.verifyError(@() run_SpeckleMapping(sv, zeros(4, 4, 4, 2, 'single')), ...
                'Umitoolbox:run_SpeckleMapping:unsupportedLayout');
            testCase.verifyError(@() run_SpeckleMapping(sv, zeros(4, 4, 'single')), ...
                'Umitoolbox:run_SpeckleMapping:unsupportedLayout');

            % Layouts other than Y-X-T.
            layouts = {{'Y','X','T','E'}, [6 5 4 2]; {'Y','X','E'}, [6 5 3]; {'Y','X'}, [6 5]};
            for k = 1:size(layouts, 1)
                writeTestDat(fullfile(sv, 'red.dat'), rand([layouts{k, 2}, 1], 'single'), ...
                    10, 5, 'DimNames', layouts{k, 1});
                testCase.verifyError(@() run_SpeckleMapping(sv, 'red.dat'), ...
                    'Umitoolbox:run_SpeckleMapping:unsupportedLayout');
            end
        end

        function testMissingFileInputErrors(testCase)
            delete(fullfile(testCase.TempFolder, 'speckle.dat'));

            testCase.verifyError( ...
                @() run_SpeckleMapping(testCase.TempFolder, 'speckle.dat'), ...
                'Umitoolbox:run_SpeckleMapping:FileNotFound');
        end
    end
end