classdef TestSpeckleMapping < matlab.unittest.TestCase
    %TESTSPECKLEMAPPING Unit tests for SpeckleMapping.
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
                ['TestSpeckleMapping_' char(java.util.UUID.randomUUID)]);
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
        function testStandardTemporalReturnsUMT(testCase)
            out = SpeckleMapping(testCase.TempFolder, iArray(testCase), 'sType', 'Temporal', 'bSaveMap', false, 'bLogScale', false);

            testCase.verifyTrue(isstruct(out) && isscalar(out));
            validateUMTStruct(out, 'requireEventInfo', false);
            testCase.verifyEqual(char(string(out.kind)), 'image');
            testCase.verifyTrue(isfield(out.data, 'SpeckleMap'));

            entry = out.data.SpeckleMap;
            testCase.verifyEqual(cellstr(string(entry.dimNames)), {'Y','X'});
            testCase.verifyClass(entry.value, 'single');
        end

        function testStandardSpatialReturnsUMT(testCase)
            out = SpeckleMapping(testCase.TempFolder, iArray(testCase), 'sType', 'Spatial', 'bSaveMap', false, 'bLogScale', false);

            validateUMTStruct(out, 'requireEventInfo', false);
            entry = out.data.SpeckleMap;
            testCase.verifyEqual(cellstr(string(entry.dimNames)), {'Y','X'});
            testCase.verifyClass(entry.value, 'single');
        end

        function testLowRAMTemporalReturnsUMT(testCase)
            out = SpeckleMapping(testCase.TempFolder, 'speckle.dat', 'sType', 'Temporal', 'bSaveMap', false, 'bLogScale', false);

            validateUMTStruct(out, 'requireEventInfo', false);
            testCase.verifyEqual(char(string(out.kind)), 'image');
            testCase.verifyTrue(isfield(out.data, 'SpeckleMap'));
            testCase.verifyClass(out.data.SpeckleMap.value, 'single');
        end

        function testLogScaleOption(testCase)
            outNoLog = SpeckleMapping(testCase.TempFolder, iArray(testCase), 'sType', 'Temporal', 'bSaveMap', false, 'bLogScale', false);
            outLog = SpeckleMapping(testCase.TempFolder, iArray(testCase), 'sType', 'Temporal', 'bSaveMap', false, 'bLogScale', true);

            a = outNoLog.data.SpeckleMap.value;
            b = outLog.data.SpeckleMap.value;

            testCase.verifySize(a, size(b));
            testCase.verifyFalse(isequaln(a, b));
        end

        function testSaveMapWritesTiff(testCase)
            tifFile = fullfile(testCase.TempFolder, 'std_speckle.tiff');
            if isfile(tifFile)
                delete(tifFile);
            end

            SpeckleMapping(testCase.TempFolder, iArray(testCase), 'sType', 'Temporal', 'bSaveMap', true, 'bLogScale', false);

            testCase.verifyTrue(isfile(tifFile));
        end

        function testInvalidSTypeErrors(testCase)
            testCase.verifyError( ...
                @() SpeckleMapping(testCase.TempFolder, 'speckle.dat', 'sType', 'InvalidMode', 'bSaveMap', false, 'bLogScale', false), ...
                'Umitoolbox:SpeckleMapping:InvalidSType');
        end

        function testMissingFileErrors(testCase)
            delete(fullfile(testCase.TempFolder, 'speckle.dat'));

            testCase.verifyError( ...
                @() SpeckleMapping(testCase.TempFolder, 'speckle.dat', 'sType', 'Temporal', 'bSaveMap', false, 'bLogScale', false), ...
                'Umitoolbox:SpeckleMapping:FileNotFound');
        end

        function testFilenameSelectsLowRAMAndArraySelectsStandard(testCase)
            % The type of "data" is the only mode switch: a filename and the
            % same data as an array give the same map.
            for sType = {'Spatial', 'Temporal'}
                fromFile = SpeckleMapping(testCase.TempFolder, 'speckle.dat', ...
                    'sType', sType{1}, 'bSaveMap', false, 'bLogScale', true);
                fromArray = SpeckleMapping(testCase.TempFolder, iArray(testCase), ...
                    'sType', sType{1}, 'bSaveMap', false, 'bLogScale', true);

                testCase.verifyEqual(fromArray.data.SpeckleMap.value, ...
                    fromFile.data.SpeckleMap.value, 'AbsTol', single(1e-5), ...
                    sprintf('Array and file maps differ for sType=%s.', sType{1}));
            end
        end

        function testArrayThatFitsInRAMIsFilteredWhole(testCase)
            % With ample RAM an array takes the single-STDFILT path: the map
            % is the plain whole-array algorithm.
            arr = iArray(testCase);
            dat = arr ./ mean(arr, 3, 'omitnan');
            kernels = struct('Spatial', single(fspecial('disk', 2) > 0), ...
                'Temporal', ones(1, 1, 5, 'single'));

            for sType = {'Spatial', 'Temporal'}
                expected = single(mean(stdfilt(dat, kernels.(sType{1})), 3, 'omitnan'));
                out = SpeckleMapping(testCase.TempFolder, arr, 'sType', sType{1}, 'bSaveMap', false, 'bLogScale', false);

                testCase.verifyEqual(out.data.SpeckleMap.value, expected, ...
                    sprintf('Whole-array map differs for sType=%s.', sType{1}));
            end
        end

        function testDoubleArrayMatchesSingleArray(testCase)
            % Array input is converted to single chunk by chunk.
            arr = iArray(testCase);

            fromSingle = SpeckleMapping(testCase.TempFolder, arr, 'sType', 'Temporal', 'bSaveMap', false, 'bLogScale', false);
            fromDouble = SpeckleMapping(testCase.TempFolder, double(arr), 'sType', 'Temporal', 'bSaveMap', false, 'bLogScale', false);

            testCase.verifyEqual(fromDouble.data.SpeckleMap.value, ...
                fromSingle.data.SpeckleMap.value);
        end

        function testRejectsUnsupportedArrays(testCase)
            testCase.verifyError( ...
                @() SpeckleMapping(testCase.TempFolder, zeros(4, 4, 'single'), 'sType', 'Temporal', 'bSaveMap', false, 'bLogScale', false), ...
                'Umitoolbox:SpeckleMapping:unsupportedLayout');
            testCase.verifyError( ...
                @() SpeckleMapping(testCase.TempFolder, {1}, 'sType', 'Temporal', 'bSaveMap', false, 'bLogScale', false), ...
                'Umitoolbox:SpeckleMapping:UnsupportedInputType');
        end
    end
end

function arr = iArray(testCase)
%IARRAY The fixture recording as an in-RAM Y-X-T array.
arr = loadData(fullfile(testCase.TempFolder, 'speckle.dat'));
end