classdef TestEventsManagerSelectionAndSplit < matlab.unittest.TestCase
    %TESTEVENTSMANAGERSELECTIONANDSPLIT Tests for event selection and splitting.

    properties
        TempFolder char
        Obj
    end

    methods (TestMethodSetup)
        function createTestObject(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture
            fx = testCase.applyFixture(TemporaryFolderFixture);
            testCase.TempFolder = fx.Folder;
            addpath(fullfile(fileparts(mfilename('fullpath')), 'EMTestHelpers'));
            writeMinimalAcqInfosMat(testCase.TempFolder, 'FrameRateHz', 10, 'AISampleRate', 100);
            testCase.Obj = EventsManager(testCase.TempFolder, '', 'csv');
            testCase.Obj.EventFileParseMethod = 'csv';
        end
    end

    methods (Test)
        function testRemoveRepetitionIgnoresOnlyTargetTransitions(testCase)
            prepareLabeledExternalObject(testCase);

            testCase.Obj.removeRepetition('A', 2);

            testCase.verifyEqual(testCase.Obj.selectedEvents, logical([1;1;1;1;0;0]));

            [evIdx, condList, repList] = testCase.Obj.getEventIndex('A');
            testCase.verifyEqual(evIdx, logical([1;0;0]));
            testCase.verifyEqual(condList, uint16(1));
            testCase.verifyEqual(repList, uint16(1));
        end

        function testRemoveConditionIgnoresAllConditionTransitions(testCase)
            prepareLabeledExternalObject(testCase);

            testCase.Obj.removeCondition('B');

            testCase.verifyEqual(testCase.Obj.selectedEvents, logical([1;1;0;0;1;1]));

            [tm, st] = testCase.Obj.getConditionTimestamps('A');
            testCase.verifyEqual(numel(tm), 4);
            testCase.verifyEqual(st, logical([1;0;1;0]));
        end

        function testClearIgnoredEventsResetsSelection(testCase)
            prepareLabeledExternalObject(testCase);

            testCase.Obj.removeCondition('B');
            testCase.Obj.removeRepetition('A', 2);
            testCase.Obj.clearIgnoredEvents();

            testCase.verifyTrue(all(testCase.Obj.selectedEvents));
            testCase.verifySize(testCase.Obj.selectedEvents, size(testCase.Obj.eventID));
        end

        function testGetFrameMatrixAllEvents(testCase)
            prepareLabeledExternalObject(testCase);

            [frMat, conditionList, repetitionList] = testCase.Obj.getFrameMatrix(40, 'FrameRateHz', testCase.Obj.AcqInfo.FrameRateHz);

            expFrMat = [4:13; 14:23; 24:33];
            testCase.verifyEqual(frMat, expFrMat);
            testCase.verifyEqual(conditionList, uint16([1;2;1]));
            testCase.verifyEqual(repetitionList, uint16([1;1;2]));
        end

        function testGetFrameMatrixConditionAndRepetitionFilter(testCase)
            prepareLabeledExternalObject(testCase);

            [frMat, conditionList, repetitionList] = testCase.Obj.getFrameMatrix(40, 'A', 2, 'FrameRateHz', testCase.Obj.AcqInfo.FrameRateHz);

            testCase.verifyEqual(frMat, 24:33);
            testCase.verifyEqual(conditionList, uint16(1));
            testCase.verifyEqual(repetitionList, uint16(2));
        end

        function testExplicitBaselineHasNoFramePeriodBound(testCase)
            % .dat header Phase 7a: the baseline is in seconds, with no
            % frame-period bound from AcqInfos.mat (10 Hz here, i.e. 0.1 s).
            prepareLabeledExternalObject(testCase);   % triggers 1 s apart
            testCase.Obj.setBaselinePeriod(0.05);      % shorter than one 10 Hz frame
            testCase.verifyEqual(double(testCase.Obj.baselinePeriod), 0.05, 'AbsTol', 1e-6);
            testCase.Obj.setBaselinePeriod(0.95);      % within one frame of the interval
            testCase.verifyEqual(double(testCase.Obj.baselinePeriod), 0.95, 'AbsTol', 1e-6);
            testCase.verifyError(@() testCase.Obj.setBaselinePeriod(1.0), ?MException);
        end

        function testBaselineIgnoresAcqInfosFrameRate(testCase)
            % An absurd AcqInfos.mat frame rate has no effect on the baseline.
            writeMinimalAcqInfosMat(testCase.TempFolder, 'FrameRateHz', 0.001, 'AISampleRate', 100);
            testCase.Obj = EventsManager(testCase.TempFolder, '', 'csv');
            testCase.Obj.EventFileParseMethod = 'csv';
            prepareLabeledExternalObject(testCase);
            testCase.Obj.setBaselinePeriod();
            testCase.verifyEqual(double(testCase.Obj.baselinePeriod), 0.2, 'AbsTol', 1e-6);
            testCase.Obj.setBaselinePeriod(0.3);
            testCase.verifyEqual(double(testCase.Obj.baselinePeriod), 0.3, 'AbsTol', 1e-6);
        end

        function testGetFrameMatrixUsesTheDataFrameRate(testCase)
            % .dat header Phase 6b-2: event times convert to frames with the
            % data's own rate. AcqInfos.mat says 10 Hz; the data is 20 Hz.
            prepareLabeledExternalObject(testCase);

            [frMat, conditionList, repetitionList] = testCase.Obj.getFrameMatrix(80, ...
                'FrameRateHz', 20);
            testCase.verifyEqual(frMat, [7:26; 27:46; 47:66]);
            testCase.verifyEqual(conditionList, uint16([1;2;1]));
            testCase.verifyEqual(repetitionList, uint16([1;1;2]));

            % Positional condition/repetition filters still work with it.
            frMat2 = testCase.Obj.getFrameMatrix(80, 'A', 2, 'FrameRateHz', 20);
            testCase.verifyEqual(frMat2, 47:66);
        end

        function testCallsWithoutRateError(testCase)
            % .dat header Phase 8a: the AcqInfos.mat rate fallback is removed;
            % every event-to-frame conversion needs the data's own rate.
            prepareLabeledExternalObject(testCase);
            data = reshape(single(1:2*2*40), 2, 2, 40);
            testCase.verifyError(@() testCase.Obj.getFrameMatrix(40), ...
                'Umitoolbox:EventsManager:missingFrameRate');
            testCase.verifyError(@() testCase.Obj.getFrameMatrix(40, 'A', 2), ...
                'Umitoolbox:EventsManager:missingFrameRate');
            testCase.verifyError(@() testCase.Obj.splitDataByEvents(data), ...
                'Umitoolbox:EventsManager:missingFrameRate');
            testCase.verifyError(@() testCase.Obj.exportEventInfo(), ...
                'Umitoolbox:EventsManager:missingFrameRate');
        end

        function testSplitDataByEventsAndExportUseTheDataFrameRate(testCase)
            prepareLabeledExternalObject(testCase);
            data = reshape(single(1:2*2*80), 2, 2, 80);

            [dataByEv, conditionList] = testCase.Obj.splitDataByEvents(data, 'FrameRateHz', 20);
            testCase.verifySize(dataByEv, [2 2 20 3]);
            testCase.verifyEqual(dataByEv(:, :, :, 2), data(:, :, 27:46));
            testCase.verifyEqual(conditionList, uint16([1;2;1]));

            evInfo = testCase.Obj.exportEventInfo('FrameRateHz', 20);
            testCase.verifyEqual(evInfo.FrameRateHz, 20);
        end

        function testSplitDataByEventsUsesFrameMatrixExactly(testCase)
            prepareLabeledExternalObject(testCase);

            T = 40;
            data = zeros(2, 3, T, 'single');
            for t = 1:T
                data(:,:,t) = single(t);
            end

            [frMat, conditionList1, repetitionList1] = testCase.Obj.getFrameMatrix(size(data,3), 'FrameRateHz', testCase.Obj.AcqInfo.FrameRateHz);
            [dataByEv, conditionList2, repetitionList2] = testCase.Obj.splitDataByEvents(data, 'FrameRateHz', testCase.Obj.AcqInfo.FrameRateHz);

            testCase.verifyEqual(conditionList2, conditionList1);
            testCase.verifyEqual(repetitionList2, repetitionList1);
            testCase.verifyEqual(size(dataByEv), [2 3 10 3]);

            for ii = 1:size(frMat,1)
                testCase.verifyEqual(dataByEv(:,:,:,ii), data(:,:,frMat(ii,:)));
            end
        end
        function testExportEventInfoReflectsIgnoredOnEvents(testCase)
            prepareLabeledExternalObject(testCase);

            testCase.Obj.removeCondition('B');
            evInfo = testCase.Obj.exportEventInfo('FrameRateHz', testCase.Obj.AcqInfo.FrameRateHz);

            testCase.verifyEqual(evInfo.eventNameList(:)', {'A','B'});
            testCase.verifyEqual(evInfo.eventID, uint16([1;2;1]));
            testCase.verifyEqual(evInfo.selectedEvents, logical([1;0;1]));
            testCase.verifyEqual(evInfo.selected, logical([1;0;1]));
            testCase.verifyEqual(evInfo.repetitionIndex, uint16([1;1;2]));
            testCase.verifyEqual(evInfo.durationSec, [0.1;0.1;0.1], 'AbsTol', 1e-6);
            testCase.verifyFalse(isfield(evInfo, 'PurgedEvents'));
            testCase.verifyEqual(evInfo.baselinePeriod, testCase.Obj.baselinePeriod, 'AbsTol', 1e-6);

            % 'IncludeIgnored', false keeps only the selected instances.
            evSel = testCase.Obj.exportEventInfo('FrameRateHz', 10, 'IncludeIgnored', false);
            testCase.verifyEqual(evSel.eventID, uint16([1;1]));
            testCase.verifyEqual(evSel.selected, true(2, 1));
            testCase.verifyEqual(evInfo.FrameRateHz, testCase.Obj.AcqInfo.FrameRateHz);
        end

        function testSplitDataByEventsReturnsYXTEAndMatchesFrameMatrix(testCase)
            prepareLabeledExternalObject(testCase);

            data = reshape(single(1:(4 * 3 * 40)), 4, 3, 40);

            [frMat, conditionIDlist, repetitionList] = testCase.Obj.getFrameMatrix(size(data,3), 'FrameRateHz', testCase.Obj.AcqInfo.FrameRateHz);
            [dataByEv, conditionIDlist2, repetitionList2] = testCase.Obj.splitDataByEvents(data, 'FrameRateHz', testCase.Obj.AcqInfo.FrameRateHz);

            testCase.verifyEqual(conditionIDlist2, conditionIDlist);
            testCase.verifyEqual(repetitionList2, repetitionList);

            frMatExpected = frMat;
            if any(isnan(frMatExpected(:)))
                firstNaNCol = find(any(isnan(frMatExpected),1), 1, 'first');
                if ~isempty(firstNaNCol)
                    frMatExpected(:, firstNaNCol:end) = [];
                end
            end

            testCase.verifySize(dataByEv, ...
                [size(data,1), size(data,2), size(frMatExpected,2), size(frMatExpected,1)]);

            for ii = 1:size(frMatExpected,1)
                testCase.verifyFalse(any(isnan(frMatExpected(ii,:))));
                testCase.verifyEqual(dataByEv(:,:,:,ii), data(:,:,frMatExpected(ii,:)));
            end
        end

        function testIncludeIgnoredKeepsAllInstancesWithFlags(testCase)
            % .dat header Phase 8b: saved split data keeps ignored instances;
            % the selected ones are identical to the default split.
            prepareLabeledExternalObject(testCase);
            testCase.Obj.removeRepetition('A', 2);
            rate = testCase.Obj.AcqInfo.FrameRateHz;
            data = reshape(single(1:2*2*40), 2, 2, 40);

            [frAll, condAll, repAll, selAll] = testCase.Obj.getFrameMatrix(40, ...
                'FrameRateHz', rate, 'IncludeIgnored', true);
            [frSel, condSel, ~, selSel] = testCase.Obj.getFrameMatrix(40, 'FrameRateHz', rate);
            testCase.verifyEqual(selAll, logical([1;1;0]));
            testCase.verifyEqual(condAll, uint16([1;2;1]));
            testCase.verifyEqual(repAll, uint16([1;1;2]));
            testCase.verifyEqual(selSel, true(2, 1));
            testCase.verifyEqual(condSel, uint16([1;2]));
            testCase.verifyEqual(frAll(selAll, :), frSel);

            [byEvAll, ~, ~, sel] = testCase.Obj.splitDataByEvents(data, ...
                'FrameRateHz', rate, 'IncludeIgnored', true);
            byEvSel = testCase.Obj.splitDataByEvents(data, 'FrameRateHz', rate);
            testCase.verifySize(byEvAll, [2 2 size(byEvAll, 3) 3]);
            testCase.verifyEqual(sel, logical([1;1;0]));
            nT = min(size(byEvAll, 3), size(byEvSel, 3));
            testCase.verifyEqual(byEvAll(:, :, 1:nT, sel), byEvSel(:, :, 1:nT, :));

            inst = testCase.Obj.getEventInstances();
            testCase.verifyEqual(inst.eventID, [1;2;1]);
            testCase.verifyEqual(inst.selected, logical([1;1;0]));
            testCase.verifyEqual(inst.eventName, ["A";"B";"A"]);
            testCase.verifyEqual(diff(inst.onsetSec), [1;1], 'AbsTol', 1e-6);
            testCase.verifyEqual(inst.durationSec, [0.1;0.1;0.1], 'AbsTol', 1e-6);
        end

        function testEventTimeVectorZeroIsTheOnsetColumn(testCase)
            % .dat header Phase 8b: time 0 of the event view is the frame
            % getFrameMatrix uses as the onset, including a fractional
            % baseline (the 3.1133 s x 5 Hz case had a one-frame shift).
            prepareLabeledExternalObject(testCase);
            testCase.Obj.setBaselinePeriod(0.31);
            for rate = [10 5 7]
                frMat = testCase.Obj.getFrameMatrix(80, 'FrameRateHz', rate);
                inst = testCase.Obj.getEventInstances();
                t = eventTimeVector(size(frMat, 2), rate, testCase.Obj.baselinePeriod);
                onsetFrame = ceil(inst.onsetSec(1) * rate);
                col = find(frMat(1, :) == onsetFrame, 1);
                testCase.assertNotEmpty(col);
                testCase.verifyEqual(t(col), 0, 'AbsTol', 1e-12);
            end
        end

        function testPurgeBlockedByEventSplitFiles(testCase)
            % .dat header Phase 8b: no purge while event-split data exists.
            prepareBaselinePseudoEventObject(testCase);
            saveData(fullfile(testCase.TempFolder, 'byEv.dat'), zeros(2, 2, 3, 2, 'single'), ...
                'DimNames', {'Y', 'X', 'T', 'E'}, 'FrameRateHz', 10);
            before = testCase.Obj.eventID;

            testCase.verifyError(@() testCase.Obj.removeCondition('baseline', 'Purge', true), ...
                'Umitoolbox:EventsManager:purgeBlockedByEventSplitFiles');
            testCase.verifyError(@() testCase.Obj.removeRepetition('baseline', 1, 'Purge', true), ...
                'Umitoolbox:EventsManager:purgeBlockedByEventSplitFiles');
            testCase.verifyEqual(testCase.Obj.eventID, before);
            testCase.verifyEqual(EventsManager.findEventSplitFiles(testCase.TempFolder), {'byEv.dat'});

            % Ignoring (no purge) is still allowed.
            testCase.Obj.removeCondition('baseline');
            testCase.verifyEqual(testCase.Obj.selectedEvents, logical([1;1;0;0;1;1]));

            delete(fullfile(testCase.TempFolder, 'byEv.dat'));
            testCase.Obj.clearIgnoredEvents();
            testCase.Obj.removeCondition('baseline', 'Purge', true);
            testCase.verifyEqual(testCase.Obj.eventID, uint16([1;1;1;1]));
        end

        function testRemoveConditionWithoutPurgeKeepsBoundaryClipped(testCase)
            prepareBaselinePseudoEventObject(testCase);

            testCase.Obj.removeCondition('baseline');
            [frMat, conditionList, repetitionList] = testCase.Obj.getFrameMatrix(60, 'A', 'FrameRateHz', testCase.Obj.AcqInfo.FrameRateHz);

            % Default (non-purge) ignore: "baseline"'s timestamp still bounds the
            % first "A" trial, clipping it to 2 frames.
            testCase.verifyEqual(conditionList, uint16([1;1]));
            testCase.verifyEqual(repetitionList, uint16([1;2]));
            testCase.verifyEqual(frMat(1,:), [6 7 nan(1,16)]);
        end

        function testPurgeConditionExtendsAdjacentTrialLength(testCase)
            prepareBaselinePseudoEventObject(testCase);

            testCase.Obj.removeCondition('baseline', 'Purge', true);

            % .dat header Phase 8b: a purge removes the events from the list.
            testCase.verifyEqual(testCase.Obj.eventID, uint16([1;1;1;1]));
            testCase.verifyEqual(testCase.Obj.selectedEvents, true(4, 1));
            testCase.verifyEqual(testCase.Obj.PurgedEvents, false(4, 1));

            [frMat, conditionList, repetitionList] = testCase.Obj.getFrameMatrix(60, 'A', 'FrameRateHz', testCase.Obj.AcqInfo.FrameRateHz);

            % With "baseline" purged, the first "A" trial now extends across the
            % gap that used to be bounded by "baseline", instead of stopping at 2
            % frames.
            testCase.verifyEqual(conditionList, uint16([1;1]));
            testCase.verifyEqual(repetitionList, uint16([1;2]));
            testCase.verifyEqual(frMat(1,:), 6:25);
            testCase.verifyEqual(frMat(2,:), 26:45);
        end

        function testPurgeRepetitionExtendsAdjacentTrialLength(testCase)
            prepareBaselinePseudoEventObject(testCase);

            testCase.Obj.removeRepetition('baseline', 1, 'Purge', true);

            % .dat header Phase 8b: a purge removes the events from the list.
            testCase.verifyEqual(testCase.Obj.eventID, uint16([1;1;1;1]));
            testCase.verifyEqual(testCase.Obj.PurgedEvents, false(4, 1));

            [frMat] = testCase.Obj.getFrameMatrix(60, 'A', 'FrameRateHz', testCase.Obj.AcqInfo.FrameRateHz);
            testCase.verifyEqual(frMat(1,:), 6:25);
        end

        function testClearIgnoredEventsResetsPurgedEvents(testCase)
            prepareBaselinePseudoEventObject(testCase);

            testCase.Obj.removeCondition('baseline', 'Purge', true);
            testCase.Obj.clearIgnoredEvents();

            testCase.verifyTrue(all(testCase.Obj.selectedEvents));
            testCase.verifyTrue(all(~testCase.Obj.PurgedEvents));
            testCase.verifySize(testCase.Obj.PurgedEvents, size(testCase.Obj.eventID));
        end

        function testSplitDataByEventsFilteredSelectionUsesYXTE(testCase)
            prepareLabeledExternalObject(testCase);

            data = reshape(single(1:(4 * 3 * 40)), 4, 3, 40);

            [frMat, conditionIDlist, repetitionList] = ...
                testCase.Obj.getFrameMatrix(size(data,3), 'A', 2, 'FrameRateHz', testCase.Obj.AcqInfo.FrameRateHz);

            [dataByEv, conditionIDlist2, repetitionList2] = ...
                testCase.Obj.splitDataByEvents(data, 'condition', 'A', 'repetition', 2, 'FrameRateHz', testCase.Obj.AcqInfo.FrameRateHz);

            frMatExpected = frMat;
            if any(isnan(frMatExpected(:)))
                firstNaNCol = find(any(isnan(frMatExpected),1), 1, 'first');
                if ~isempty(firstNaNCol)
                    frMatExpected(:, firstNaNCol:end) = [];
                end
            end

            testCase.verifyEqual(conditionIDlist2, conditionIDlist);
            testCase.verifyEqual(repetitionList2, repetitionList);
            testCase.verifySize(dataByEv, ...
                [size(data,1), size(data,2), size(frMatExpected,2), size(frMatExpected,1)]);

            for ii = 1:size(frMatExpected,1)
                testCase.verifyEqual(dataByEv(:,:,:,ii), data(:,:,frMatExpected(ii,:)));
            end
        end


    end


    methods (Access = private)
        function prepareLabeledExternalObject(testCase)
            sr = 100;
            signal = makePulseSignal(400, [50 150 250], 10, 'Amplitude', 5);
            testCase.Obj.getTriggersFromSignal(signal, sr, false);

            csvFile = writeCSVConditionFile(testCase.TempFolder, 'labels.csv', {'A','B','A'});
            status = testCase.Obj.readConditionFile(csvFile);
            testCase.assertTrue(status);
            testCase.assertEqual(testCase.Obj.eventNameList(:)', {'A','B'});
            testCase.assertEqual(testCase.Obj.eventID, uint16([1;1;2;2;1;1]));
            testCase.assertEqual(testCase.Obj.repetitionID, uint16([1;1;1;1;2;2]));
            testCase.assertEqual(testCase.Obj.baselinePeriod, single(0.2), 'AbsTol', 1e-6);
        end

        function prepareBaselinePseudoEventObject(testCase)
            %PREPAREBASELINEPSEUDOEVENTOBJECT Two "A" repetitions with an unevenly
            % spaced "baseline" pseudo-event sitting close after the first one, so
            % boundary-clipping vs. purge-based exclusion is distinguishable.
            sr = 100;
            signal = makePulseSignal(400, [50 70 250], 10, 'Amplitude', 5);
            % The default minInterStim (2 s) would otherwise merge the closely
            % spaced "A"/"baseline" pulses into a single burst trigger.
            testCase.Obj.minInterStim = 0.05;
            testCase.Obj.getTriggersFromSignal(signal, sr, false);

            csvFile = writeCSVConditionFile(testCase.TempFolder, 'baselineLabels.csv', {'A','baseline','A'});
            status = testCase.Obj.readConditionFile(csvFile);
            testCase.assertTrue(status);
            testCase.assertEqual(testCase.Obj.eventNameList(:)', {'A','baseline'});
            testCase.assertEqual(testCase.Obj.eventID, uint16([1;1;2;2;1;1]));
            testCase.assertEqual(testCase.Obj.repetitionID, uint16([1;1;1;1;2;2]));
        end
    end
end
