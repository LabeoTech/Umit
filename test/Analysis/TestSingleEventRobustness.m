classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestSingleEventRobustness < matlab.unittest.TestCase
    %TESTSINGLEEVENTROBUSTNESS Data with exactly one event instance (E = 1).
    %
    %   MATLAB drops a trailing size-1 dimension, so Y-X-T-E data with E = 1 is
    %   held as a Y-X-T array. The dimension metadata (dimNames of a .dat
    %   header or of a UMT entry) is the authority for the axis count; no
    %   consumer may infer it from ndims. Every function below runs on the
    %   single-event fixtures of makeSingleEventFixture (non-square, frame-
    %   coded values) and must keep the E axis.

    properties
        Folder char
    end

    properties (TestParameter)
        Variant = {'plain', 'shortBaseline', 'nearEnd'}
        Input = {'dat', 'umt'}
        Filter = {'normalizeBSLN', 'normalizeZScore', 'apply_detrend', 'GSR', ...
            'spatialGaussFilt', 'normalizeLPF'}
    end

    methods (TestMethodSetup)
        function createFolder(testCase)
            testFolder = fileparts(mfilename('fullpath'));
            projectRoot = extractBefore(testFolder, [filesep 'test']);
            if isempty(projectRoot)
                projectRoot = fileparts(fileparts(testFolder));
            end
            addpath(genpath(char(projectRoot)));
            testCase.Folder = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
        end
    end

    methods (Test)
        % ---- the fixture itself -------------------------------------------

        function testFixtureHoldsOneInstanceWithMatchingMetadata(testCase, Variant)
            fx = makeSingleEventFixture(testCase.Folder, 'Variant', Variant);

            inst = EventsManager(testCase.Folder).getEventInstances();
            testCase.verifyNumElements(inst.eventID, 1);
            testCase.verifyEqual(loadMetaData(fx.splitDatFile).dimNames, {'Y','X','T','E'});
            testCase.verifyEqual(double(loadMetaData(fx.splitDatFile).dimSizes(4)), 1);
            testCase.verifyEqual(size(fx.expectedSplit, 1:4), ...
                [6, 5, size(fx.expectedSplit, 3), 1]);
            testCase.verifyEqual(fx.imageUMTSplit.data.main.dimNames, {'Y','X','T','E'});
            testCase.verifyEqual(fx.roiUMTSplit.data.traces.dimNames, {'ROI','T','E'});
            testCase.verifyNumElements(fx.imageUMTSplit.eventInfo.eventID, 1);
        end

        % ---- split_data_by_event ------------------------------------------

        function testSplitKeepsTheEAxis(testCase, Variant, Input)
            fx = makeSingleEventFixture(testCase.Folder, 'Variant', Variant);
            warning('off', 'Umitoolbox:split_data_by_event:firstTrialBaselinePadded');
            testCase.addTeardown(@() warning('on', ...
                'Umitoolbox:split_data_by_event:firstTrialBaselinePadded'));

            if strcmp(Input, 'dat')
                out = split_data_by_event(fx.datFile, testCase.Folder);
                info = loadMetaData(out);
                testCase.verifyEqual(info.dimNames, {'Y','X','T','E'});
                testCase.verifyEqual(double(info.dimSizes(:)).', ...
                    [6, 5, size(fx.expectedSplit, 3), 1]);
                value = loadData(out);
            else
                out = split_data_by_event(fx.imageUMT, testCase.Folder);
                testCase.verifyEqual(out.data.main.dimNames, {'Y','X','T','E'});
                value = out.data.main.value;
                testCase.verifyNumElements(out.eventInfo.eventID, 1);
                validateUMTStruct(out, 'requireEventInfo', true);
            end
            testCase.verifyEqual(size(value, 1:4), [6, 5, size(fx.expectedSplit, 3), 1]);
            testCase.verifyEqual(single(value), fx.expectedSplit);
        end

        function testRoiUMTSplitKeepsTheEAxis(testCase)
            fx = makeSingleEventFixture(testCase.Folder);

            out = split_data_by_event(fx.roiUMT, testCase.Folder);

            testCase.verifyEqual(out.data.traces.dimNames, {'ROI','T','E'});
            testCase.verifyEqual(size(out.data.traces.value, 1:3), ...
                [3, size(fx.expectedSplit, 3), 1]);
            testCase.verifyEqual(out.data.traces.value, ...
                reshape(fx.expectedSplit(1:3, 1, :), 3, size(fx.expectedSplit, 3)));
        end

        % ---- EventsManager ------------------------------------------------

        function testConditionSelectionGivesASingleTrial(testCase)
            fx = makeSingleEventFixture(testCase.Folder, 'Variant', 'multiCondition');
            ev = EventsManager(testCase.Folder);

            [frMat, condList, repList] = ev.getFrameMatrix(fx.nT, fx.condition, ...
                'FrameRateHz', fx.rate);
            [byEv, condB] = ev.splitDataByEvents(fx.data, 'condition', fx.condition, ...
                'FrameRateHz', fx.rate);

            testCase.verifySize(frMat, [1, 30]);
            testCase.verifyEqual(double(condList), 2);
            testCase.verifyEqual(double(repList), 1);
            testCase.verifyEqual(size(byEv, 1:4), [6, 5, 30, 1]);
            testCase.verifyEqual(byEv, fx.expectedSplit);
            testCase.verifyEqual(double(condB), 2);
        end

        function testSplitDataByEventsHandlesOneInstance(testCase)
            fx = makeSingleEventFixture(testCase.Folder);
            ev = EventsManager(testCase.Folder);

            byEv = ev.splitDataByEvents(fx.data, 'FrameRateHz', fx.rate, 'IncludeIgnored', true);

            testCase.verifyEqual(size(byEv, 1:4), [6, 5, size(fx.expectedSplit, 3), 1]);
            testCase.verifyEqual(byEv, fx.expectedSplit);
        end

        function testExportEventInfoOfOneInstance(testCase)
            fx = makeSingleEventFixture(testCase.Folder);
            info = EventsManager(testCase.Folder).exportEventInfo('FrameRateHz', fx.rate);

            testCase.verifyNumElements(info.eventID, 1);
            testCase.verifyNumElements(info.selected, 1);
            testCase.verifyEqual(info.onsetSec, 3);
        end

        function testAggregationPlanOfOneInstance(testCase)
            fx = makeSingleEventFixture(testCase.Folder);
            ev = EventsManager(testCase.Folder);
            info = ev.exportEventInfo('FrameRateHz', fx.rate);

            plan = EventsManager.conditionAggregationPlan(info);
            reduced = EventsManager.reduceByCondition(fx.expectedSplit, plan, ...
                @(x) mean(x, 4, 'omitnan'), 4);

            testCase.verifyEqual(size(reduced, 1:4), [6, 5, size(fx.expectedSplit, 3), 1]);
            testCase.verifyEqual(reduced, fx.expectedSplit);
        end

        % ---- event mapping of a single-slice .dat -------------------------

        function testMappingOfASingleSliceDat(testCase)
            fx = makeSingleEventFixture(testCase.Folder);

            m = resolveDatEventMapping(loadMetaData(fx.splitDatFile), testCase.Folder);

            testCase.verifyEqual(m.status, 'matched');
            testCase.verifyEqual(m.nE, 1);
            testCase.verifyNumElements(m.eventInfo.eventID, 1);
        end

        % ---- consumers of a single-slice Y-X-T-E .dat ---------------------

        function testPerTrialFiltersKeepTheEAxis(testCase, Filter)
            fx = makeSingleEventFixture(testCase.Folder);
            iSuppressWarnings(testCase);

            out = feval(Filter, fx.splitDatFile, testCase.Folder);

            info = loadMetaData(iFullPath(testCase.Folder, out));
            testCase.verifyEqual(info.dimNames, {'Y','X','T','E'});
            testCase.verifyEqual(double(info.dimSizes(:)).', ...
                [6, 5, size(fx.expectedSplit, 3), 1]);
            value = loadData(iFullPath(testCase.Folder, out));
            testCase.verifyEqual(size(value, 1:4), [6, 5, size(fx.expectedSplit, 3), 1]);
            testCase.verifyFalse(all(isnan(value), 'all'));
        end

        function testAggregationOverTKeepsTheEAxis(testCase)
            fx = makeSingleEventFixture(testCase.Folder);

            out = apply_aggregate_function(fx.splitDatFile, testCase.Folder, ...
                'aggregateFcn', 'mean', 'dimensionName', 'T');

            % Aggregating over T removes the T axis: Y-X-E with E = 1.
            info = loadMetaData(iFullPath(testCase.Folder, out));
            testCase.verifyEqual(info.dimNames, {'Y','X','E'});
            testCase.verifyEqual(double(info.dimSizes(:)).', [6, 5, 1]);
            testCase.verifyEqual(loadData(iFullPath(testCase.Folder, out)), ...
                mean(fx.expectedSplit, 3), 'AbsTol', single(1e-3));
        end

        function testAggregationOverEOfASingleSliceIsTheSlice(testCase)
            fx = makeSingleEventFixture(testCase.Folder);

            out = apply_aggregate_function(fx.splitDatFile, testCase.Folder, ...
                'aggregateFcn', 'mean', 'dimensionName', 'E');

            info = loadMetaData(iFullPath(testCase.Folder, out));
            testCase.verifyEqual(info.dimNames, {'Y','X','T','E'});
            testCase.verifyEqual(double(info.dimSizes(:)).', ...
                [6, 5, size(fx.expectedSplit, 3), 1]);
            testCase.verifyEqual(loadData(iFullPath(testCase.Folder, out)), ...
                fx.expectedSplit, 'AbsTol', single(1e-3));
        end

        function testAggregationOfAnImageUMTWithOneEvent(testCase)
            fx = makeSingleEventFixture(testCase.Folder);

            outT = apply_aggregate_function(fx.imageUMTSplit, testCase.Folder, ...
                'dimensionName', 'T');
            outE = apply_aggregate_function(fx.imageUMTSplit, testCase.Folder, ...
                'dimensionName', 'E');

            testCase.verifyEqual(outT.data.main.dimNames, {'Y','X','E'});
            testCase.verifyEqual(size(outT.data.main.value, 1:3), [6, 5, 1]);
            testCase.verifyEqual(outE.data.main.dimNames, {'Y','X','T','E'});
            testCase.verifyEqual(size(outE.data.main.value, 1:4), ...
                [6, 5, size(fx.expectedSplit, 3), 1]);
        end

        function testAggregationOfContinuousDataSplitsOneInstance(testCase)
            fx = makeSingleEventFixture(testCase.Folder);

            out = apply_aggregate_function(fx.datFile, testCase.Folder, ...
                'aggregateFcn', 'mean', 'dimensionName', 'E');

            info = loadMetaData(iFullPath(testCase.Folder, out));
            testCase.verifyEqual(info.dimNames, {'Y','X','T','E'});
            testCase.verifyEqual(double(info.dimSizes(4)), 1);
            testCase.verifyEqual(loadData(iFullPath(testCase.Folder, out)), ...
                fx.expectedSplit, 'AbsTol', single(1e-3));
        end

        function testAmplitudeMapsOfASingleSlice(testCase)
            fx = makeSingleEventFixture(testCase.Folder);

            out = genAmplitudeMaps(fx.splitDatFile, testCase.Folder);

            info = loadMetaData(iFullPath(testCase.Folder, out));
            testCase.verifyEqual(info.dimNames(1:2), {'Y','X'});
            testCase.verifyEqual(double(info.dimSizes(1:2)), [6, 5]);
        end

        function testTiffExportOfASingleSlice(testCase)
            fx = makeSingleEventFixture(testCase.Folder);

            out = ConvertToTiff(fx.splitDatFile, testCase.Folder);

            testCase.verifyTrue(isfile(iFullPath(testCase.Folder, out)) || ...
                ~isempty(dir(fullfile(testCase.Folder, '*.tif*'))));
        end

        function testHemoCorrectionOfASingleSlice(testCase)
            fx = makeSingleEventFixture(testCase.Folder);
            writeTestDat(fullfile(testCase.Folder, 'ref.dat'), ...
                fx.data .* single(1.01) + single(0.2), fx.rate);

            out = HemoCorrection(fx.splitDatFile, testCase.Folder, ...
                'Algorithm', 'Ratiometric', 'Red', false, 'Green', false, ...
                'Amber', false, 'Other', 'ref.dat');

            info = loadMetaData(iFullPath(testCase.Folder, out));
            testCase.verifyEqual(info.dimNames, {'Y','X','T','E'});
            testCase.verifyEqual(double(info.dimSizes(4)), 1);
        end

        function testResponseFeaturesOfASingleInstance(testCase, Variant)
            fx = makeSingleEventFixture(testCase.Folder, 'Variant', Variant);
            iSuppressWarnings(testCase);
            file = fullfile(testCase.Folder, 'roiOne.umt');
            saveData(file, fx.roiUMTSplit);

            out = calculateResponseFeatures(file);

            name = fieldnames(out.data);
            entry = out.data.(name{1});
            testCase.verifyEqual(cellstr(string(entry.dimNames)), {'ROI','Measure','E'});
            testCase.verifyEqual(size(entry.value, 1:3), [3, 6, 1]);
            validateUMTStruct(out, 'requireEventInfo', true);
        end

        % ---- UMT helpers --------------------------------------------------

        function testUMTWithOneEventValidatesAndRoundTrips(testCase)
            fx = makeSingleEventFixture(testCase.Folder);

            validateUMTStruct(fx.imageUMTSplit, 'requireEventInfo', true);
            validateUMTStruct(fx.roiUMTSplit, 'requireEventInfo', true);
            file = fullfile(testCase.Folder, 'one.umt');
            saveData(file, fx.imageUMTSplit);
            back = loadData(file);

            testCase.verifyEqual(back.data.main.dimNames, {'Y','X','T','E'});
            testCase.verifyEqual(size(back.data.main.value, 1:4), ...
                size(fx.imageUMTSplit.data.main.value, 1:4));
            testCase.verifyNumElements(back.eventInfo.eventID, 1);
        end
    end
end

function path = iFullPath(folder, name)
%IFULLPATH A function output (file name or full path) as a full path.
if iscell(name)
    name = name{1};
end
name = char(string(name));
if isfile(name)
    path = name;
else
    path = fullfile(folder, name);
end
end

function iSuppressWarnings(testCase)
%ISUPPRESSWARNINGS Silence the known warnings of a filled first trial for one test.
testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture( ...
    {'Umitoolbox:split_data_by_event:firstTrialBaselinePadded', ...
     'Umitoolbox:calculateResponseFeatures:firstTrialBaselineFilled'}));
end
