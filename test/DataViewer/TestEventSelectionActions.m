classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestEventSelectionActions < matlab.unittest.TestCase
    %TESTEVENTSELECTIONACTIONS Ignore/restore repetitions and save the selection.
    %
    %   DFR-20261006-001. promptIgnoreRepetitions is exercised through its
    %   non-interactive mode (no dialog is opened; the dialog itself is checked
    %   manually). saveEventSelection is DataViewer's explicit save of the event
    %   selection to events.mat.

    properties
        Folder char
        Obj
    end

    methods (TestMethodSetup)
        function createEvents(testCase)
            addpath(genpath(fileparts(fileparts(fileparts(mfilename('fullpath'))))));
            testCase.Folder = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            writeMinimalAcqInfosMat(testCase.Folder, 'FrameRateHz', 60, 'AISampleRate', 10000);
            obj = EventsManager(testCase.Folder, '', 'csv');
            signal = makePulseSignal(22000, [3000 7000 11000], 60, 'Amplitude', 5);
            obj.getTriggersFromSignal(signal, 10000, false);
            obj.saveEvents(testCase.Folder);
            testCase.Obj = obj;
        end
    end

    methods (Test)
        function dialogHelperAppliesOnlyTheDifference(testCase)
            ignoreSecond = promptIgnoreRepetitions([], testCase.Obj, 'extSignal', ...
                'Interactive', false, 'Ignored', 2);
            testCase.verifyEqual(ignoreSecond.severity, 'info', ignoreSecond.message);
            testCase.verifyTrue(ignoreSecond.changed);
            testCase.verifyEqual(testCase.Obj.listRepetitions('extSignal').Ignored, [false; true; false]);
            testCase.verifySubstring(ignoreSecond.message, '1 repetition(s) ignored, 0 included');

            swap = promptIgnoreRepetitions([], testCase.Obj, 'extSignal', ...
                'Interactive', false, 'Ignored', 3);
            testCase.verifyTrue(swap.changed);
            testCase.verifyEqual(testCase.Obj.listRepetitions('extSignal').Ignored, [false; false; true]);
            testCase.verifySubstring(swap.message, '1 repetition(s) ignored, 1 included');

            same = promptIgnoreRepetitions([], testCase.Obj, 'extSignal', ...
                'Interactive', false, 'Ignored', 3);
            testCase.verifyFalse(same.changed);
            testCase.verifySubstring(same.message, 'no repetition changed');

            clearAll = promptIgnoreRepetitions([], testCase.Obj, 'extSignal', ...
                'Interactive', false, 'Ignored', []);
            testCase.verifyTrue(clearAll.changed);
            testCase.verifyTrue(all(testCase.Obj.selectedEvents));
        end

        function ignoringEveryRepetitionSaysTheConditionIsIgnored(testCase)
            result = promptIgnoreRepetitions([], testCase.Obj, 'extSignal', ...
                'Interactive', false, 'Ignored', [1 2 3]);
            testCase.verifyTrue(result.changed);
            testCase.verifySubstring(result.message, 'whole condition is ignored');
            testCase.verifyFalse(any(testCase.Obj.selectedEvents));

            % A wholly ignored condition can be restored one repetition at a time.
            restored = promptIgnoreRepetitions([], testCase.Obj, 'extSignal', ...
                'Interactive', false, 'Ignored', [1 3]);
            testCase.verifyTrue(restored.changed, restored.message);
            testCase.verifyEqual(testCase.Obj.listRepetitions('extSignal').Ignored, [true; false; true]);
        end

        function ignoringEveryRepetitionOfTheOnlyConditionIsRefusedWhenRequired(testCase)
            % The Events Manager passes RequireIncludedEvents: the only condition
            % keeps at least one repetition. Without the option the helper allows it.
            before = testCase.Obj.selectedEvents;
            refused = promptIgnoreRepetitions([], testCase.Obj, 'extSignal', ...
                'Interactive', false, 'Ignored', [1 2 3], 'RequireIncludedEvents', true);
            testCase.verifyFalse(refused.changed);
            testCase.verifyEqual(refused.severity, 'warning');
            testCase.verifySubstring(refused.message, 'At least one repetition must stay included');
            testCase.verifyEqual(testCase.Obj.selectedEvents, before);

            allowed = promptIgnoreRepetitions([], testCase.Obj, 'extSignal', ...
                'Interactive', false, 'Ignored', [1 2], 'RequireIncludedEvents', true);
            testCase.verifyTrue(allowed.changed, allowed.message);
            testCase.verifyEqual(testCase.Obj.listRepetitions('extSignal').Ignored, [true; true; false]);
        end

        function dialogHelperReportsErrorsWithoutThrowing(testCase)
            before = testCase.Obj.selectedEvents;
            unknownRepetition = promptIgnoreRepetitions([], testCase.Obj, 'extSignal', ...
                'Interactive', false, 'Ignored', 9);
            testCase.verifyEqual(unknownRepetition.severity, 'error');
            testCase.verifyFalse(unknownRepetition.changed);

            unknownCondition = promptIgnoreRepetitions([], testCase.Obj, 'Nope', ...
                'Interactive', false, 'Ignored', 1);
            testCase.verifyEqual(unknownCondition.severity, 'error');
            testCase.verifyEqual(testCase.Obj.selectedEvents, before);
        end

        function saveWritesTheSelectionAndKeepsEverythingElse(testCase)
            eventsFile = fullfile(testCase.Folder, 'events.mat');
            before = load(eventsFile);
            testCase.verifyTrue(all(before.selectedEvents));

            testCase.Obj.removeRepetition('extSignal', 2);
            result = saveEventSelection(testCase.Obj, testCase.Folder);

            testCase.verifyTrue(result.changed, result.message);
            testCase.verifyEqual(result.severity, 'info');
            testCase.verifySubstring(result.message, '1 ignored repetition(s)');
            after = load(eventsFile);
            testCase.verifyEqual(after.selectedEvents, testCase.Obj.selectedEvents);
            testCase.verifyFalse(isequal(after.selectedEvents, before.selectedEvents));
            for name = setdiff(fieldnames(before), {'selectedEvents'})'
                testCase.verifyEqual(after.(name{1}), before.(name{1}), ...
                    ['events.mat variable "' name{1} '" must be unchanged.']);
            end

            % A new object (as DataViewer builds on load) sees the saved selection.
            reloaded = EventsManager(testCase.Folder, '', 'csv');
            testCase.verifyEqual(reloaded.listRepetitions('extSignal').Ignored, [false; true; false]);

            again = saveEventSelection(testCase.Obj, testCase.Folder);
            testCase.verifyFalse(again.changed);
            testCase.verifySubstring(again.message, 'already saved');
        end

        function saveReportsMissingFileAndBadInput(testCase)
            emptyFolder = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            noFile = saveEventSelection(testCase.Obj, emptyFolder);
            testCase.verifyFalse(noFile.changed);
            testCase.verifyEqual(noFile.severity, 'warning');
            testCase.verifyFalse(isfile(fullfile(emptyFolder, 'events.mat')));

            noSource = saveEventSelection([], testCase.Folder);
            testCase.verifyEqual(noSource.severity, 'error');
            noFolder = saveEventSelection(testCase.Obj, fullfile(testCase.Folder, 'missing'));
            testCase.verifyEqual(noFolder.severity, 'error');
        end
    end
end
