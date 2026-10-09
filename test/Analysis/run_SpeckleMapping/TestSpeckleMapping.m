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
            out = SpeckleMapping(iArray(testCase), testCase.TempFolder, 'Temporal', false, false);

            testCase.verifyTrue(isstruct(out) && isscalar(out));
            validateUMTStruct(out, 'requireEventInfo', false);
            testCase.verifyEqual(char(string(out.kind)), 'image');
            testCase.verifyTrue(isfield(out.data, 'SpeckleMap'));

            entry = out.data.SpeckleMap;
            testCase.verifyEqual(cellstr(string(entry.dimNames)), {'Y','X'});
            testCase.verifyClass(entry.value, 'single');
        end

        function testStandardSpatialReturnsUMT(testCase)
            out = SpeckleMapping(iArray(testCase), testCase.TempFolder, 'Spatial', false, false);

            validateUMTStruct(out, 'requireEventInfo', false);
            entry = out.data.SpeckleMap;
            testCase.verifyEqual(cellstr(string(entry.dimNames)), {'Y','X'});
            testCase.verifyClass(entry.value, 'single');
        end

        function testLowRAMTemporalReturnsUMT(testCase)
            out = SpeckleMapping('speckle.dat', testCase.TempFolder, 'Temporal', false, false);

            validateUMTStruct(out, 'requireEventInfo', false);
            testCase.verifyEqual(char(string(out.kind)), 'image');
            testCase.verifyTrue(isfield(out.data, 'SpeckleMap'));
            testCase.verifyClass(out.data.SpeckleMap.value, 'single');
        end

        function testLogScaleOption(testCase)
            outNoLog = SpeckleMapping(iArray(testCase), testCase.TempFolder, 'Temporal', false, false);
            outLog = SpeckleMapping(iArray(testCase), testCase.TempFolder, 'Temporal', false, true);

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

            SpeckleMapping(iArray(testCase), testCase.TempFolder, 'Temporal', true, false);

            testCase.verifyTrue(isfile(tifFile));
        end

        function testInvalidSTypeErrors(testCase)
            testCase.verifyError( ...
                @() SpeckleMapping('speckle.dat', testCase.TempFolder, 'InvalidMode', false, false), ...
                'Umitoolbox:SpeckleMapping:InvalidSType');
        end

        function testMissingFileErrors(testCase)
            delete(fullfile(testCase.TempFolder, 'speckle.dat'));

            testCase.verifyError( ...
                @() SpeckleMapping('speckle.dat', testCase.TempFolder, 'Temporal', false, false), ...
                'Umitoolbox:SpeckleMapping:FileNotFound');
        end

        function testFilenameSelectsLowRAMAndArraySelectsStandard(testCase)
            % The type of "data" is the only mode switch: a filename and the
            % same data as an array give the same map.
            for sType = {'Spatial', 'Temporal'}
                fromFile = SpeckleMapping('speckle.dat', testCase.TempFolder, ...
                    sType{1}, false, true);
                fromArray = SpeckleMapping(iArray(testCase), testCase.TempFolder, ...
                    sType{1}, false, true);

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
                out = SpeckleMapping(arr, testCase.TempFolder, sType{1}, false, false);

                testCase.verifyEqual(out.data.SpeckleMap.value, expected, ...
                    sprintf('Whole-array map differs for sType=%s.', sType{1}));
            end
        end

        function testDoubleArrayMatchesSingleArray(testCase)
            % Array input is converted to single chunk by chunk.
            arr = iArray(testCase);

            fromSingle = SpeckleMapping(arr, testCase.TempFolder, 'Temporal', false, false);
            fromDouble = SpeckleMapping(double(arr), testCase.TempFolder, 'Temporal', false, false);

            testCase.verifyEqual(fromDouble.data.SpeckleMap.value, ...
                fromSingle.data.SpeckleMap.value);
        end

        function testRejectsUnsupportedArrays(testCase)
            testCase.verifyError( ...
                @() SpeckleMapping(zeros(4, 4, 'single'), testCase.TempFolder, 'Temporal', false, false), ...
                'Umitoolbox:SpeckleMapping:UnsupportedLayout');
            testCase.verifyError( ...
                @() SpeckleMapping({1}, testCase.TempFolder, 'Temporal', false, false), ...
                'Umitoolbox:SpeckleMapping:UnsupportedInputType');
        end
    end
end

function arr = iArray(testCase)
%IARRAY The fixture recording as an in-RAM Y-X-T array.
arr = loadData(fullfile(testCase.TempFolder, 'speckle.dat'));
end