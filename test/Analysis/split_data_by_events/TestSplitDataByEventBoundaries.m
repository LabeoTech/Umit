classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestSplitDataByEventBoundaries < matlab.unittest.TestCase
    %TESTSPLITDATABYEVENTBOUNDARIES split_data_by_event when trials touch the recording edges.
    %
    %   The trial of an event starts baselinePeriod before its onset, so an
    %   event that comes sooner than the baseline after the recording start
    %   (or an onset after the recording end) leaves frames outside the
    %   recording. The FIRST event may do so: its missing baseline frames
    %   are filled with the mean of the baseline frames it has (with a
    %   warning). Any other trial outside the recording is an error. Every
    %   case runs for a .dat input (LOW-RAM MODE) and a UMT input, which must
    %   behave and report identically.
    %
    %   All fixtures are synthetic: a 10 Hz recording of 6 x 5 pixels (non
    %   square on purpose) whose value encodes its frame, with an events.mat
    %   written directly (onsets/baseline in seconds, exactly representable
    %   at 10 Hz so the frame arithmetic is exact). A trial holds
    %   round(baseline*rate) baseline columns, then the onset column; the
    %   first trial's column k is frame (start + k - 1).

    properties (TestParameter)
        Input = {'dat', 'umt'}
    end

    properties
        TempFolder char
    end

    properties (Constant)
        Rate = 10
        NY = 6
        NX = 5
        ErrId = 'Umitoolbox:split_data_by_event:noEventFrames'
        PadId = 'Umitoolbox:split_data_by_event:firstTrialBaselinePadded'
    end

    methods (TestMethodSetup)
        function createTempFolder(testCase)
            testFolder = fileparts(mfilename('fullpath'));
            projectRoot = extractBefore(testFolder, [filesep 'test']);
            if isempty(projectRoot)
                projectRoot = fileparts(fileparts(testFolder));
            end
            addpath(genpath(char(projectRoot)));

            fx = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.TempFolder = fx.Folder;
        end
    end

    methods (Test)
        % ---- first event sooner than the baseline: baseline is filled ----

        function testFirstEventSoonerThanBaselineIsFilledWithTheBaselineMean(testCase, Input)
            % Onset at 1 s (frame 10), baseline 2 s (20 columns): the first
            % trial starts at frame -9, so its first 10 columns are missing
            % and take the mean of the 10 baseline frames it has (frames 1-10).
            data = iData(testCase, 100);
            iWriteEvents(testCase.TempFolder, [1 4 7], 2);

            [out, warnId, warnMsg] = iSplitWarn(testCase, Input, data);

            testCase.verifyEqual(warnId, testCase.PadId);
            testCase.verifySubstring(warnMsg, 'starts 10 frames (1 s) before the recording');
            testCase.verifySubstring(warnMsg, 'mean of its 10 valid baseline');
            fill = mean(data(:, :, 1:10), 3);
            first = cat(3, repmat(fill, 1, 1, 10), data(:, :, 1:20));
            expected = cat(4, first, data(:, :, 21:50), data(:, :, 51:80));
            testCase.verifyEqual(iResult(Input, out), expected);
        end

        function testFirstEventOneFrameSoonerThanBaselineIsFilled(testCase, Input)
            % Baseline 2.1 s = 21 columns, onset frame 20: start frame 0, one
            % missing column, filled with the mean of frames 1-20.
            data = iData(testCase, 100);
            iWriteEvents(testCase.TempFolder, [2 5 8], 2.1);

            [out, warnId] = iSplitWarn(testCase, Input, data);

            testCase.verifyEqual(warnId, testCase.PadId);
            first = cat(3, mean(data(:, :, 1:20), 3), data(:, :, 1:29));
            expected = cat(4, first, data(:, :, 30:59), data(:, :, 60:89));
            testCase.verifyEqual(iResult(Input, out), expected);
        end

        function testFirstEventWithASingleBaselineFrameIsFilledFromIt(testCase, Input)
            % Onset at 0.0625 s (frame 1): only frame 1 exists before the
            % onset, so the 19 missing columns repeat it. The first trial is
            % one frame shorter than the others, so all are cropped to 29.
            data = iData(testCase, 100);
            iWriteEvents(testCase.TempFolder, [0.0625 3 6], 2);

            [out, warnId, warnMsg] = iSplitWarn(testCase, Input, data);

            testCase.verifyEqual(warnId, testCase.PadId);
            testCase.verifySubstring(warnMsg, 'mean of its 1 valid baseline');
            first = cat(3, repmat(data(:, :, 1), 1, 1, 19), data(:, :, 1:10));
            expected = cat(4, first, data(:, :, 11:39), data(:, :, 41:69));
            testCase.verifyEqual(iResult(Input, out), expected);
        end

        function testFirstEventAtExactlyTheBaselineHasNoPaddingAndNoWarning(testCase, Input)
            % Onset at 2 s, baseline 2 s: the first trial starts at frame 1.
            data = iData(testCase, 100);
            iWriteEvents(testCase.TempFolder, [2 5 8], 2);

            [out, warnId] = iSplitWarn(testCase, Input, data);

            testCase.verifyNotEqual(warnId, testCase.PadId);
            testCase.verifyEqual(iResult(Input, out), iExpected(data, {1:30, 31:60, 61:90}));
        end

        function testSingleEventSoonerThanBaselineIsFilledAndRunsToTheEnd(testCase, Input)
            % A lone event's trial extends from its start to the end of the
            % data; its first 10 baseline columns are missing.
            data = iData(testCase, 100);
            iWriteEvents(testCase.TempFolder, 1, 2);

            [out, warnId] = iSplitWarn(testCase, Input, data);

            testCase.verifyEqual(warnId, testCase.PadId);
            expected = cat(3, repmat(mean(data(:, :, 1:10), 3), 1, 1, 10), data);
            testCase.verifyEqual(iResult(Input, out), expected);
        end

        function testSingleEventAfterTheBaselineRunsToTheEnd(testCase, Input)
            data = iData(testCase, 100);
            iWriteEvents(testCase.TempFolder, 3, 2);

            out = iResult(Input, iSplit(testCase, Input, data));

            testCase.verifyEqual(out, iExpected(data, {11:100}));
        end

        function testIgnoredFirstEventIsStillSplitAndFilled(testCase, Input)
            % Ignored instances are saved too (IncludeIgnored): the ignored
            % first event keeps its E slice, filled like any first trial.
            data = iData(testCase, 100);
            iWriteEvents(testCase.TempFolder, [1 4 7], 2, 'Selected', [false true true]);

            [out, warnId] = iSplitWarn(testCase, Input, data);

            testCase.verifyEqual(warnId, testCase.PadId);
            testCase.verifySize(iResult(Input, out), [testCase.NY, testCase.NX, 30, 3]);
        end

        % ---- trials that cannot be filled ---------------------------------

        function testFirstEventAtTheRecordingStartErrors(testCase, Input)
            % Onset at 0 s: no recording frame belongs to the first baseline.
            data = iData(testCase, 100);
            iWriteEvents(testCase.TempFolder, [0 3 6], 2);

            ME = iCapture(@() iSplit(testCase, Input, data));

            testCase.assertNotEmpty(ME, 'The split must fail.');
            testCase.verifyEqual(ME.identifier, testCase.ErrId);
            testCase.verifySubstring(ME.message, 'starts outside the recording');
            testCase.verifySubstring(ME.message, 'Recording: 100 frames at 10 Hz');
            testCase.verifySubstring(ME.message, 'Baseline period: 2 s');
            testCase.verifySubstring(ME.message, 'event 1: starts at frame -19');
            testCase.verifySubstring(ME.message, '(1 of 3 event instances)');
            testCase.verifySubstring(ME.message, 'Only the first event may start before');
        end

        function testFillableFirstEventIsNotListedAsOffending(testCase, Input)
            % The first event (onset 1 s) could be filled, but the last one
            % is after the recording end: the split fails and only the last
            % one is reported.
            data = iData(testCase, 100);
            iWriteEvents(testCase.TempFolder, [1 4 12], 2);

            ME = iCapture(@() iSplit(testCase, Input, data));

            testCase.assertNotEmpty(ME);
            testCase.verifyEqual(ME.identifier, testCase.ErrId);
            testCase.verifySubstring(ME.message, '(1 of 3 event instances)');
            testCase.verifySubstring(ME.message, 'event 3: starts at frame');
            testCase.verifyFalse(contains(ME.message, 'event 1:'));
        end

        function testOffendingEventsAreListedWithACap(testCase, Input)
            % Seven events 1 s apart with an 8 s baseline: every trial starts
            % before the recording, and even the first has no frame left to
            % fill it from. Only five are listed.
            data = iData(testCase, 100);
            iWriteEvents(testCase.TempFolder, 1:7, 8);

            ME = iCapture(@() iSplit(testCase, Input, data));

            testCase.assertNotEmpty(ME);
            testCase.verifyEqual(ME.identifier, testCase.ErrId);
            testCase.verifySubstring(ME.message, '(7 of 7 event instances)');
            testCase.verifySubstring(ME.message, 'event 5: starts at frame');
            testCase.verifyFalse(contains(ME.message, 'event 6:'));
            testCase.verifySubstring(ME.message, '... and 2 more');
        end

        function testShortIntervalAndShortFirstDelayErrorsOnTheFirstEvent(testCase, Input)
            % Events 1 s apart with a 2 s baseline and the first onset at
            % 1 s: the first trial is only 10 frames long, all of them before
            % the recording, so nothing is left to fill it from.
            data = iData(testCase, 100);
            iWriteEvents(testCase.TempFolder, [1 2 3], 2);

            ME = iCapture(@() iSplit(testCase, Input, data));

            testCase.assertNotEmpty(ME);
            testCase.verifyEqual(ME.identifier, testCase.ErrId);
            testCase.verifySubstring(ME.message, 'event 1: starts at frame -9');
            testCase.verifySubstring(ME.message, '(1 of 3 event instances)');
        end

        function testEventAfterTheRecordingEndErrors(testCase, Input)
            % The last onset (12 s = frame 120) lies past the 100-frame
            % recording: its trial starts at frame 111.
            data = iData(testCase, 100);
            iWriteEvents(testCase.TempFolder, [2 5 12], 1);

            ME = iCapture(@() iSplit(testCase, Input, data));

            testCase.assertNotEmpty(ME);
            testCase.verifyEqual(ME.identifier, testCase.ErrId);
            testCase.verifySubstring(ME.message, 'event 3: starts at frame 111');
            testCase.verifySubstring(ME.message, 'Starting outside frames 1-100');
        end

        function testNoEventInstanceErrors(testCase, Input)
            data = iData(testCase, 100);
            iWriteEvents(testCase.TempFolder, zeros(1, 0), 2);

            ME = iCapture(@() iSplit(testCase, Input, data));

            testCase.assertNotEmpty(ME);
            testCase.verifyEqual(ME.identifier, testCase.ErrId);
            testCase.verifySubstring(ME.message, 'no event instance');
        end

        % ---- interval between events shorter than the baseline -----------

        function testIntervalShorterThanBaselineGivesContiguousTrials(testCase, Input)
            % Onsets 1 s apart with a 2 s baseline, first onset late enough.
            % Characterization: each trial is one interval long and starts
            % 2 s before its onset, so the trials tile the recording from
            % frame 11 and contain no frame after their own onset.
            data = iData(testCase, 100);
            iWriteEvents(testCase.TempFolder, [3 4 5], 2);

            out = iResult(Input, iSplit(testCase, Input, data));

            testCase.verifyEqual(out, iExpected(data, {11:20, 21:30, 31:40}));
        end

        function testUnevenShortIntervalsAreCroppedToTheShortestTrial(testCase, Input)
            % Intervals of 1 s and 2 s with a 2 s baseline: the last two
            % trials are 20 frames long, the first only 10, so every trial
            % is cropped to the shortest.
            data = iData(testCase, 100);
            iWriteEvents(testCase.TempFolder, [3 4 6], 2);

            out = iResult(Input, iSplit(testCase, Input, data));

            testCase.verifyEqual(out, iExpected(data, {11:20, 21:30, 41:50}));
        end

        % ---- trailing edge of the recording -------------------------------

        function testLastTrialRunningPastTheEndIsCropped(testCase, Input)
            % The recording ends at frame 85: the last trial (61:90) loses
            % frames past the end, so every trial is cropped to 25 frames.
            data = iData(testCase, 85);
            iWriteEvents(testCase.TempFolder, [2 5 8], 2);

            out = iResult(Input, iSplit(testCase, Input, data));

            testCase.verifyEqual(out, iExpected(data, {1:25, 31:55, 61:85}));
        end

        function testFilledFirstTrialAndCroppedLastTrialTogether(testCase, Input)
            % Both edges at once: the first trial is filled and the last one
            % runs past the end (recording of 70 frames).
            data = iData(testCase, 70);
            iWriteEvents(testCase.TempFolder, [1 4 7], 2);

            [out, warnId] = iSplitWarn(testCase, Input, data);

            testCase.verifyEqual(warnId, testCase.PadId);
            % Trials start at frames -9, 21, 51; the last one has frames 51:70
            % (20 valid columns of 30), so every trial keeps 20 columns.
            first = cat(3, repmat(mean(data(:, :, 1:10), 3), 1, 1, 10), data(:, :, 1:10));
            expected = cat(4, first, data(:, :, 21:40), data(:, :, 51:70));
            testCase.verifyEqual(iResult(Input, out), expected);
        end

        % ---- both paths agree and leave nothing behind --------------------

        function testFailedSplitLeavesNoOutputFile(testCase)
            % The check runs before the output is created.
            data = iData(testCase, 100);
            iWriteEvents(testCase.TempFolder, [0 3 6], 2);

            iCapture(@() iSplit(testCase, 'dat', data));

            testCase.verifyFalse(isfile(fullfile(testCase.TempFolder, 'dataByEv.dat')));
            testCase.verifyEmpty(dir(fullfile(testCase.TempFolder, '*_writing.dat')));
        end

        function testDatAndUmtErrorsHaveTheSameMessage(testCase)
            data = iData(testCase, 100);
            iWriteEvents(testCase.TempFolder, [0 3 6], 2);

            meDat = iCapture(@() iSplit(testCase, 'dat', data));
            meUmt = iCapture(@() iSplit(testCase, 'umt', data));

            testCase.assertNotEmpty(meDat);
            testCase.assertNotEmpty(meUmt);
            testCase.verifyEqual(meUmt.message, meDat.message);
        end

        function testAggregateFunctionSplitsTheFirstTrialLikeTheSplit(testCase)
            % apply_aggregate_function splits a continuous recording itself
            % (E aggregation): its trials, filled first trial included, must
            % equal split_data_by_event's, and it warns too.
            data = iData(testCase, 100);
            iWriteEvents(testCase.TempFolder, [1 4 7], 2);
            datFile = fullfile(testCase.TempFolder, 'in.dat');
            writeTestDat(datFile, data, testCase.Rate);
            splitTrials = iResult('dat', split_data_by_event(datFile, testCase.TempFolder));

            lastwarn('', '');
            outFile = apply_aggregate_function(datFile, testCase.TempFolder, ...
                'aggregateFcn', 'mean', 'dimensionName', 'E');
            [~, warnId] = lastwarn;

            testCase.verifyEqual(warnId, 'apply_aggregate_function:firstTrialBaselinePadded');
            testCase.verifyEqual(loadData(outFile), mean(splitTrials, 4), 'AbsTol', 1e-4);
        end

        function testDatAndUmtWarningsHaveTheSameMessage(testCase)
            data = iData(testCase, 100);
            iWriteEvents(testCase.TempFolder, [1 4 7], 2);

            [~, ~, msgDat] = iSplitWarn(testCase, 'dat', data);
            [~, ~, msgUmt] = iSplitWarn(testCase, 'umt', data);

            testCase.verifyNotEmpty(msgDat);
            testCase.verifyEqual(msgUmt, msgDat);
        end
    end
end

function data = iData(testCase, nT)
%IDATA Non-square Y-X-T single array whose value is its frame number.
data = repmat(reshape(single(1:nT), 1, 1, nT), testCase.NY, testCase.NX, 1) + ...
    repmat(single(0.01 * (1:testCase.NY)'), 1, testCase.NX, nT) + ...
    repmat(single(0.0001 * (1:testCase.NX)), testCase.NY, 1, nT);
end

function out = iSplit(testCase, kind, data)
%ISPLIT Run split_data_by_event on the data as a .dat file or a UMT struct.
if strcmp(kind, 'dat')
    datFile = fullfile(testCase.TempFolder, 'in.dat');
    writeTestDat(datFile, data, testCase.Rate);
    out = split_data_by_event(datFile, testCase.TempFolder);
else
    umt = genUMTStruct(data, ...
        'kind', 'image', ...
        'entryName', 'main', ...
        'dimNames', {'Y','X','T'}, ...
        'meta', struct('FrameRateHz', testCase.Rate));
    out = split_data_by_event(umt, testCase.TempFolder);
end
end

function [out, warnId, warnMsg] = iSplitWarn(testCase, kind, data)
%ISPLITWARN iSplit that also returns the last warning raised by the call.
lastwarn('', '');
out = iSplit(testCase, kind, data);
[warnMsg, warnId] = lastwarn;
end

function value = iResult(kind, out)
%IRESULT The Y-X-T-E array of a split result, whatever the input kind.
if strcmp(kind, 'dat')
    value = loadData(out);
else
    value = out.data.main.value;
end
value = single(value);
end

function expected = iExpected(data, frameSets)
%IEXPECTED Y-X-T-E array from the data and one frame range per trial,
% cropped to the shortest range.
nFrames = min(cellfun(@numel, frameSets));
expected = zeros(size(data, 1), size(data, 2), nFrames, numel(frameSets), 'single');
for k = 1:numel(frameSets)
    expected(:, :, :, k) = data(:, :, frameSets{k}(1:nFrames));
end
end

function ME = iCapture(fcn)
%ICAPTURE The MException thrown by FCN ([] if it does not throw).
ME = [];
try
    fcn();
catch ME
end
end

function iWriteEvents(folder, onsetsSec, baselineSec, varargin)
%IWRITEEVENTS Write a synthetic events.mat with one condition ("A").
%
% Every event lasts 0.5 s. Without AcqInfos.mat the file is flagged as
% external-signal events so EventsManager loads it from the events alone.
p = inputParser;
addParameter(p, 'Selected', []);
parse(p, varargin{:});

onsetsSec = onsetsSec(:);
n = numel(onsetsSec);
timestamps = reshape([onsetsSec, onsetsSec + 0.5]', [], 1);
state = repmat([true; false], n, 1);
eventID = ones(2 * n, 1, 'uint16');
eventNameList = {'A'};
selectedEvents = true(2 * n, 1);
if ~isempty(p.Results.Selected)
    selectedEvents = repelem(logical(p.Results.Selected(:)), 2);
end
baselinePeriod = baselineSec;
b_hasExternalSignal = true;

save(fullfile(folder, 'events.mat'), 'timestamps', 'state', 'eventID', ...
    'eventNameList', 'selectedEvents', 'baselinePeriod', 'b_hasExternalSignal', '-mat');
end
