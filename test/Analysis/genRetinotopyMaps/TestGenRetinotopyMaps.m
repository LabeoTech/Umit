classdef TestGenRetinotopyMaps < matlab.unittest.TestCase
    properties
        ProjectRoot
        SampleDataFolder
        TempFolder
        RawData
        MetaData
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
            if isappdata(0, 'GenRetinotopyMapsTestConfig')
                cfg = getappdata(0, 'GenRetinotopyMapsTestConfig');
            end

            if ~isempty(cfg) && isfield(cfg, 'sampleDataFolder') && isfolder(cfg.sampleDataFolder)
                testCase.SampleDataFolder = char(string(cfg.sampleDataFolder));
            else
                testCase.SampleDataFolder = fullfile( ...
                    testCase.ProjectRoot, ...
                    'test', ...
                    'Analysis', ...
                    'TestingData_retinotopy');
            end

            testCase.assertTrue(isfile(fullfile(testCase.SampleDataFolder, 'fluo.dat')));
            testCase.assertTrue(isfile(fullfile(testCase.SampleDataFolder, 'AcqInfos.mat')));
            testCase.assertTrue(isfile(fullfile(testCase.SampleDataFolder, 'events.mat')));
        end
    end

    methods (TestMethodSetup)
        function createTempFolder(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture
            fx = testCase.applyFixture(TemporaryFolderFixture);
            testCase.TempFolder = fx.Folder;

            copyfile(fullfile(testCase.SampleDataFolder, 'fluo.dat'), fullfile(testCase.TempFolder, 'fluo.dat'));
            copyfile(fullfile(testCase.SampleDataFolder, 'AcqInfos.mat'), fullfile(testCase.TempFolder, 'AcqInfos.mat'));
            copyfile(fullfile(testCase.SampleDataFolder, 'events.mat'), fullfile(testCase.TempFolder, 'events.mat'));

            testCase.RawData = single(loadData(fullfile(testCase.TempFolder, 'fluo.dat')));
            testCase.MetaData = loadMetaData(fullfile(testCase.TempFolder, 'fluo.dat'));
        end
    end

    methods (Test)
        function testPipelineInfo(testCase)
            info = genRetinotopyMaps('pipelineInfo');
            testCase.verifyEqual(info.name, 'genRetinotopyMaps');
            testCase.verifyEqual(info.outputs(1).type, {'ProcessedData'});

            for parameterName = {'ViewingDist_cm', 'ScreenXsize_cm', 'ScreenYsize_cm'}
                parameter = info.parameters(strcmp( ...
                    {info.parameters.name}, parameterName{1}));
                testCase.verifyEqual(parameter.allowed, [0 Inf]);
            end
        end

        function testNumericInputAllDirections(testCase)
            out = genRetinotopyMaps(testCase.RawData, testCase.TempFolder, 'FrameRateHz', testCase.MetaData.frameRateHz, 'Direction', 'All');
            testCase.verifyTrue(isstruct(out));
            testCase.verifyTrue(isfield(out.data, 'AzimuthMap'));
            testCase.verifyTrue(isfield(out.data, 'ElevationMap'));
            testCase.verifyEqual(size(out.data.AzimuthMap.value, 3), 2);
            testCase.verifyEqual(size(out.data.ElevationMap.value, 3), 2);
        end

        function testIgnoredSweepIsExcluded(testCase)
            % .dat header Phase 8c: an ignored sweep keeps its place in the
            % sweep structure but no longer contributes to the maps; the
            % numeric and chunked .dat paths agree.
            rate = testCase.MetaData.frameRateHz;
            full = genRetinotopyMaps(testCase.RawData, testCase.TempFolder, 'FrameRateHz', rate, 'Direction', 'All');

            ev = EventsManager(testCase.TempFolder);
            inst = ev.getEventInstances();
            ev.removeRepetition(char(inst.eventName(1)), inst.repetitionIndex(1));
            ev.saveEvents(testCase.TempFolder);

            out = genRetinotopyMaps(testCase.RawData, testCase.TempFolder, 'FrameRateHz', rate, 'Direction', 'All');
            testCase.verifySize(out.data.AzimuthMap.value, size(full.data.AzimuthMap.value));
            testCase.verifyFalse(isequaln(out.data.AzimuthMap.value, full.data.AzimuthMap.value) && ...
                isequaln(out.data.ElevationMap.value, full.data.ElevationMap.value), ...
                'ignoring a sweep must change the maps');

            outDat = genRetinotopyMaps(fullfile(testCase.TempFolder, 'fluo.dat'), testCase.TempFolder, 'Direction', 'All');
            testCase.verifyNumericEquivalent(out.data.AzimuthMap.value, outDat.data.AzimuthMap.value);
            testCase.verifyNumericEquivalent(out.data.ElevationMap.value, outDat.data.ElevationMap.value);
        end

        function testDatInputMatchesNumeric(testCase)
            outNumeric = genRetinotopyMaps(testCase.RawData, testCase.TempFolder, 'FrameRateHz', testCase.MetaData.frameRateHz, 'Direction', 'All');
            outDat = genRetinotopyMaps(fullfile(testCase.TempFolder, 'fluo.dat'), testCase.TempFolder, 'Direction', 'All');
            testCase.verifyNumericEquivalent(outNumeric.data.AzimuthMap.value, outDat.data.AzimuthMap.value);
            testCase.verifyNumericEquivalent(outNumeric.data.ElevationMap.value, outDat.data.ElevationMap.value);
        end

        function testRejectsUMTStructUMTFileAndEventSplitInput(testCase)
            rate = testCase.MetaData.frameRateHz;
            umt = genUMTStruct(testCase.RawData, 'kind', 'image', ...
                'entryName', 'main', 'dimNames', {'Y','X','T'});
            umtFile = fullfile(testCase.TempFolder, 'input.umt');
            fclose(fopen(umtFile, 'w'));
            byEvent = fullfile(testCase.TempFolder, 'byEvent.dat');
            saveData(byEvent, cat(4, testCase.RawData(:,:,1:10), testCase.RawData(:,:,1:10)), ...
                'DimNames', {'Y','X','T','E'}, 'FrameRateHz', rate);

            testCase.verifyError(@() genRetinotopyMaps(umt, testCase.TempFolder, ...
                'FrameRateHz', rate), 'Umitoolbox:genRetinotopyMaps:WrongInput');
            testCase.verifyError(@() genRetinotopyMaps(umtFile, testCase.TempFolder, ...
                'FrameRateHz', rate), 'Umitoolbox:genRetinotopyMaps:WrongInput');
            testCase.verifyError(@() genRetinotopyMaps(byEvent, testCase.TempFolder), ...
                'Umitoolbox:genRetinotopyMaps:unsupportedLayout');
            testCase.verifyError(@() genRetinotopyMaps( ...
                cat(4, testCase.RawData(:,:,1:10), testCase.RawData(:,:,1:10)), ...
                testCase.TempFolder, 'FrameRateHz', rate), ...
                'Umitoolbox:genRetinotopyMaps:WrongInput');
        end

        function testForcedMultiSlabAverageMovieMatchesNumeric(testCase)
            testCase.assumeEqual(computer, 'PCWIN64', ...
                'Chunk forcing relies on shadowing the PCWIN64 memory() built-in.');
            mocksFolder = fullfile(testCase.ProjectRoot, ...
                'test', 'subFunc', 'calculateMaxChunkSize', 'mocks');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(mocksFolder));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_MODE', 'fixed'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_TOTAL', '10000'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_AVAILABLE', '5000'));

            S = load(fullfile(testCase.TempFolder, 'events.mat'));
            framestamps = round(S.timestamps * testCase.MetaData.frameRateHz);
            indxOn = find(S.eventID == 1 & S.state == 1);
            indxOff = find(S.eventID == 1 & S.state == 0);
            baselineLength = round(mean( ...
                framestamps(indxOn(2:end)) - framestamps(indxOff(1:end-1))));
            trialLength = round(mean(framestamps(indxOff) - framestamps(indxOn)));
            totalLength = trialLength + baselineLength;
            bytesPerX = double(datAxisSize(testCase.MetaData, 'Y')) * 4 * ...
                (double(totalLength) + double(baselineLength) + ...
                 2 * double(trialLength) + 4);
            testCase.verifyGreaterThan( ...
                calculateMaxChunkSize( ...
                    bytesPerX * double(datAxisSize(testCase.MetaData, 'X')), 1, .1), ...
                1, 'The regression fixture must force more than one X slab.');

            outNumeric = genRetinotopyMaps( ...
                testCase.RawData, testCase.TempFolder, ...
                'FrameRateHz', testCase.MetaData.frameRateHz, ...
                'Direction', 'All', 'b_useAverageMovie', true);
            outDat = genRetinotopyMaps( ...
                fullfile(testCase.TempFolder, 'fluo.dat'), testCase.TempFolder, ...
                'Direction', 'All', 'b_useAverageMovie', true);

            testCase.verifyNumericEquivalent( ...
                outNumeric.data.AzimuthMap.value, outDat.data.AzimuthMap.value);
            testCase.verifyNumericEquivalent( ...
                outNumeric.data.ElevationMap.value, outDat.data.ElevationMap.value);
        end

        function testForcedMultiSlabConcatenatedSweepsMatchesNumeric(testCase)
            % Concatenated-sweeps mode reads only each direction's stimulus
            % frames; with several X slabs per direction the maps must still
            % equal the in-RAM result.
            testCase.assumeEqual(computer, 'PCWIN64', ...
                'Chunk forcing relies on shadowing the PCWIN64 memory() built-in.');
            mocksFolder = fullfile(testCase.ProjectRoot, ...
                'test', 'subFunc', 'calculateMaxChunkSize', 'mocks');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(mocksFolder));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_MODE', 'fixed'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_TOTAL', '10000'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_AVAILABLE', '5000'));

            outNumeric = genRetinotopyMaps( ...
                testCase.RawData, testCase.TempFolder, ...
                'FrameRateHz', testCase.MetaData.frameRateHz, 'Direction', 'All');
            progress = evalc(['outDat = genRetinotopyMaps(' ...
                'fullfile(testCase.TempFolder, ''fluo.dat''), testCase.TempFolder, ' ...
                '''Direction'', ''All'');']);

            slabCounts = cellfun(@(c) str2double(c{1}), ...
                regexp(progress, 'slab \[\d+/(\d+)\]', 'tokens'));
            testCase.verifyGreaterThan(max(slabCounts), 1, ...
                'The fixture must force more than one X slab per direction.');

            testCase.verifyNumericEquivalent( ...
                outNumeric.data.AzimuthMap.value, outDat.data.AzimuthMap.value);
            testCase.verifyNumericEquivalent( ...
                outNumeric.data.ElevationMap.value, outDat.data.ElevationMap.value);
        end

        function testAzimuthOnly(testCase)
            testCase.rewriteEventsForMode('Azimuth_only');
            out = genRetinotopyMaps(testCase.RawData, testCase.TempFolder, 'FrameRateHz', testCase.MetaData.frameRateHz, 'Direction', 'Azimuth_only');
            testCase.verifyTrue(isfield(out.data, 'AzimuthMap'));
            testCase.verifyFalse(isfield(out.data, 'ElevationMap'));
        end

        function testElevationOnly(testCase)
            testCase.rewriteEventsForMode('Elevation_only');
            out = genRetinotopyMaps(testCase.RawData, testCase.TempFolder, 'FrameRateHz', testCase.MetaData.frameRateHz, 'Direction', 'Elevation_only');
            testCase.verifyTrue(isfield(out.data, 'ElevationMap'));
            testCase.verifyFalse(isfield(out.data, 'AzimuthMap'));
        end

        function testGenericLabelsWarnAndRun(testCase)
            testCase.rewriteEventsForMode('Generic');
            testCase.verifyWarning(@() genRetinotopyMaps(testCase.RawData, testCase.TempFolder, 'FrameRateHz', testCase.MetaData.frameRateHz, 'Direction', 'All'), ...
                'Umitoolbox:genRetinotopyMaps:GenericDirectionLabels');
        end

        function testGenericLabelsInvalidSweepsError(testCase)
            testCase.rewriteEventsForMode('Generic', 3);
            testCase.verifyError(@() genRetinotopyMaps(testCase.RawData, testCase.TempFolder, 'FrameRateHz', testCase.MetaData.frameRateHz, 'Direction', 'All'), ...
                'Umitoolbox:genRetinotopyMaps:InvalidSweeps');
        end

        function testMissingEventsMat(testCase)
            delete(fullfile(testCase.TempFolder, 'events.mat'));
            testCase.verifyError(@() genRetinotopyMaps(testCase.RawData, testCase.TempFolder, 'FrameRateHz', testCase.MetaData.frameRateHz), ...
                'Umitoolbox:genRetinotopyMaps:MissingInput');
        end

        function testRejectUnmatchedOnsetAndOffset(testCase)
            S = load(fullfile(testCase.TempFolder, 'events.mat'));
            offsetIdx = find(S.state == 0, 1, 'last');
            S.state(offsetIdx) = [];
            S.timestamps(offsetIdx) = [];
            S.eventID(offsetIdx) = [];
            testCase.saveEvents(S);

            testCase.verifyError(@() genRetinotopyMaps(testCase.RawData, testCase.TempFolder, 'FrameRateHz', testCase.MetaData.frameRateHz), ...
                'Umitoolbox:genRetinotopyMaps:InvalidEvents');
        end

        function testRejectsInitialBaselineOutsideRecordingInDatMode(testCase)
            S = load(fullfile(testCase.TempFolder, 'events.mat'));
            onsetIdx = find(S.state == 1, 1, 'first');
            S.timestamps(onsetIdx) = 1;
            testCase.saveEvents(S);

            testCase.verifyError(@() genRetinotopyMaps( ...
                fullfile(testCase.TempFolder, 'fluo.dat'), testCase.TempFolder, ...
                'b_useAverageMovie', true), ...
                'Umitoolbox:genRetinotopyMaps:InvalidBaseline');
        end

        function testVisualAngleRescaling(testCase)
            out = genRetinotopyMaps(testCase.RawData, testCase.TempFolder, 'FrameRateHz', testCase.MetaData.frameRateHz, ...
                'Direction', 'All', ...
                'ViewingDist_cm', 20, ...
                'ScreenXsize_cm', 30, ...
                'ScreenYsize_cm', 20);

            azPhase = out.data.AzimuthMap.value(:,:,2);
            elPhase = out.data.ElevationMap.value(:,:,2);

            testCase.verifyGreaterThanOrEqual(min(azPhase, [], 'all'), -90);
            testCase.verifyLessThanOrEqual(max(azPhase, [], 'all'), 90);
            testCase.verifyGreaterThanOrEqual(min(elPhase, [], 'all'), -90);
            testCase.verifyLessThanOrEqual(max(elPhase, [], 'all'), 90);
        end

        function testVisualAngleMapsFixedPhaseRangeLinearly(testCase)
            % The calibration is the linear map of the theoretical phase
            % range [0, 2*pi] onto [-VA, +VA] (pi -> 0 degrees), not a rescale
            % of the observed extrema.
            rate = testCase.MetaData.frameRateHz;
            radians = genRetinotopyMaps(testCase.RawData, testCase.TempFolder, ...
                'FrameRateHz', rate, 'Direction', 'All');
            degrees = genRetinotopyMaps(testCase.RawData, testCase.TempFolder, ...
                'FrameRateHz', rate, 'Direction', 'All', ...
                'ViewingDist_cm', 20, 'ScreenXsize_cm', 30, 'ScreenYsize_cm', 20);

            vaAz = atand(30 / (2 * 20));
            vaEl = atand(20 / (2 * 20));
            azRad = radians.data.AzimuthMap.value(:,:,2);
            elRad = radians.data.ElevationMap.value(:,:,2);
            testCase.assertGreaterThanOrEqual(min(azRad, [], 'all'), 0);
            testCase.assertLessThanOrEqual(max(azRad, [], 'all'), 2*pi);

            % (double() because verifyEqual ignores tolerances for singles)
            testCase.verifyEqual(double(degrees.data.AzimuthMap.value(:,:,2)), ...
                -vaAz + double(azRad) / (2*pi) * 2 * vaAz, 'AbsTol', 1e-3);
            testCase.verifyEqual(double(degrees.data.ElevationMap.value(:,:,2)), ...
                -vaEl + double(elRad) / (2*pi) * 2 * vaEl, 'AbsTol', 1e-3);
            testCase.verifyEqual(degrees.data.AzimuthMap.value(:,:,1), ...
                radians.data.AzimuthMap.value(:,:,1));
            testCase.verifyClass(degrees.data.AzimuthMap.value, 'single');
        end

        function testVisualAngleIndependentOfPhaseSubset(testCase)
            % Every pixel is transformed independently, so cropping the field
            % of view (which changes the observed phase extrema) must leave
            % each remaining pixel's visual angle unchanged.
            rate = testCase.MetaData.frameRateHz;
            geometry = {'ViewingDist_cm', 20, 'ScreenXsize_cm', 30, 'ScreenYsize_cm', 20};
            full = genRetinotopyMaps(testCase.RawData, testCase.TempFolder, ...
                'FrameRateHz', rate, 'Direction', 'All', geometry{:});

            rows = 1:6;
            cols = 1:5;
            cropped = genRetinotopyMaps(testCase.RawData(rows, cols, :), ...
                testCase.TempFolder, 'FrameRateHz', rate, 'Direction', 'All', ...
                geometry{:});
            radiansFull = genRetinotopyMaps(testCase.RawData, testCase.TempFolder, ...
                'FrameRateHz', rate, 'Direction', 'All');
            radiansCrop = genRetinotopyMaps(testCase.RawData(rows, cols, :), ...
                testCase.TempFolder, 'FrameRateHz', rate, 'Direction', 'All');

            % Precondition: the crop really changes the observed phase range,
            % which is what made the previous min/max calibration data-dependent.
            fullRange = [min(radiansFull.data.AzimuthMap.value(:,:,2), [], 'all'), ...
                max(radiansFull.data.AzimuthMap.value(:,:,2), [], 'all')];
            cropRange = [min(radiansCrop.data.AzimuthMap.value(:,:,2), [], 'all'), ...
                max(radiansCrop.data.AzimuthMap.value(:,:,2), [], 'all')];
            testCase.assertNotEqual(fullRange, cropRange);

            for entry = {'AzimuthMap', 'ElevationMap'}
                name = entry{1};
                testCase.verifyEqual(double(cropped.data.(name).value(:,:,2)), ...
                    double(full.data.(name).value(rows, cols, 2)), 'AbsTol', 1e-3);
            end
        end
    end

    methods (Access = private)
        function rewriteEventsForMode(testCase, modeName, maxOnsets)
            if nargin < 3
                maxOnsets = [];
            end

            S = load(fullfile(testCase.TempFolder, 'events.mat'));

            onsetIdx = find(S.state == 1);
            offsetIdx = find(S.state == 0);
            testCase.assertEqual(numel(onsetIdx), numel(offsetIdx));

            switch lower(modeName)
                case 'azimuth_only'
                    keepNames = {'0','180'};
                    keepMask = ismember(S.eventNameList(S.eventID(onsetIdx)), keepNames);
                    keepOn = onsetIdx(keepMask);
                    keepOff = offsetIdx(keepMask);

                    S.state = [ones(numel(keepOn),1); zeros(numel(keepOff),1)];
                    S.timestamps = [S.timestamps(keepOn); S.timestamps(keepOff)];

                    % Important: in Azimuth_only mode, IDs must be 1 and 2.
                    oldIDs = S.eventID(keepOn);
                    newIDsOn = zeros(size(oldIDs));
                    newIDsOn(oldIDs == find(strcmp(S.eventNameList,'0'),1,'first')) = 1;
                    newIDsOn(oldIDs == find(strcmp(S.eventNameList,'180'),1,'first')) = 2;
                    S.eventID = [newIDsOn; newIDsOn];
                    S.eventNameList = {'0','180'};

                case 'elevation_only'
                    keepNames = {'90','270'};
                    keepMask = ismember(S.eventNameList(S.eventID(onsetIdx)), keepNames);
                    keepOn = onsetIdx(keepMask);
                    keepOff = offsetIdx(keepMask);

                    S.state = [ones(numel(keepOn),1); zeros(numel(keepOff),1)];
                    S.timestamps = [S.timestamps(keepOn); S.timestamps(keepOff)];

                    % Important: in Elevation_only mode, IDs must be 1 and 2.
                    oldIDs = S.eventID(keepOn);
                    newIDsOn = zeros(size(oldIDs));
                    newIDsOn(oldIDs == find(strcmp(S.eventNameList,'90'),1,'first')) = 1;
                    newIDsOn(oldIDs == find(strcmp(S.eventNameList,'270'),1,'first')) = 2;
                    S.eventID = [newIDsOn; newIDsOn];
                    S.eventNameList = {'90','270'};

                case 'generic'
                    if isempty(maxOnsets)
                        maxOnsets = numel(onsetIdx);
                    end
                    keepOn = onsetIdx(1:maxOnsets);
                    keepOff = offsetIdx(1:maxOnsets);
                    S.state = [ones(numel(keepOn),1); zeros(numel(keepOff),1)];
                    S.timestamps = [S.timestamps(keepOn); S.timestamps(keepOff)];
                    S.eventID = ones(numel(S.state),1);
                    S.eventNameList = {'Main'};

                otherwise
                    error('Unknown rewrite mode');
            end

            eventID = S.eventID; %#ok<NASGU>
            timestamps = S.timestamps; %#ok<NASGU>
            state = S.state; %#ok<NASGU>
            eventNameList = S.eventNameList; %#ok<NASGU>
            save(fullfile(testCase.TempFolder, 'events.mat'), 'eventID', 'timestamps', 'state', 'eventNameList');
        end

        function verifyNumericEquivalent(testCase, a, b)
            testCase.verifyEqual(size(a), size(b));
            if isequaln(a,b)
                return
            end
            diffVals = double(a(:)) - double(b(:));
            diffVals = diffVals(isfinite(diffVals));
            if isempty(diffVals)
                testCase.verifyEqual(a,b);
                return
            end
            testCase.verifyLessThanOrEqual(std(diffVals, 0, 'omitnan'), 1e-4);
        end

        function saveEvents(testCase, S)
            eventID = S.eventID; %#ok<NASGU>
            timestamps = S.timestamps; %#ok<NASGU>
            state = S.state; %#ok<NASGU>
            eventNameList = S.eventNameList; %#ok<NASGU>
            save(fullfile(testCase.TempFolder, 'events.mat'), ...
                'eventID', 'timestamps', 'state', 'eventNameList');
        end
    end
end
