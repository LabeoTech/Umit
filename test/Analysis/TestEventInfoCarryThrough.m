classdef TestEventInfoCarryThrough < matlab.unittest.TestCase
    %TESTEVENTINFOCARRYTHROUGH Per-instance functions keep the event flags (.dat header Phase 8c).
    %
    %   Functions that keep the E axis must carry the input eventInfo
    %   intact (selected, durationSec, baselinePeriod), so a later reducer
    %   can still exclude the ignored instances. Synthetic event-split UMT.
    %   Only functions that still take UMT input are covered: spatialGaussFilt,
    %   apply_detrend, normalizeLPF and normalizeBSLN take arrays and/or .dat files only, and
    %   an event-split .dat carries no flags (they come from events.mat).

    properties
        Folder char
        UMT struct
    end

    methods (TestMethodSetup)
        function buildInput(testCase)
            testCase.Folder = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            data = 1 + rand(8, 8, 40, 3, 'single');
            umt = genUMTStruct(data, 'kind', 'image', 'entryName', 'main', ...
                'dimNames', {'Y', 'X', 'T', 'E'}, 'meta', struct('FrameRateHz', 10));
            eventInfo = struct('eventID', uint16([1; 2; 1]), 'repetitionIndex', uint16([1; 1; 2]), ...
                'eventName', ["A"; "B"; "A"], 'eventAxisMode', 'instances', ...
                'selected', logical([1; 0; 1]), 'durationSec', [0.5; 0.6; 0.5], ...
                'baselinePeriod', 1);
            testCase.UMT = appendUMTEventInfo(umt, 'eventInfo', eventInfo, 'overwrite', true);
        end
    end

    methods (Test)
        function aggregateOverTKeepsFlags(testCase)
            out = apply_aggregate_function(testCase.UMT, testCase.Folder, 'dimensionName', 'T');
            iVerifyCarried(testCase, out, testCase.UMT.eventInfo);
        end

        function aggregateOverEExcludesIgnored(testCase)
            out = apply_aggregate_function(testCase.UMT, testCase.Folder, 'dimensionName', 'E');
            value = testCase.UMT.data.main.value;
            testCase.verifyEqual(out.data.main.value(:, :, :, 1), ...
                mean(value(:, :, :, [1 3]), 4, 'omitnan'), 'AbsTol', 1e-6);
            testCase.verifyTrue(all(isnan(out.data.main.value(:, :, :, 2)), 'all'), ...
                'condition B: its only instance is ignored');
            testCase.verifyEqual(out.eventInfo.nInstances, [2; 0]);
            testCase.verifyEqual(out.eventInfo.selected, [true; false]);
            testCase.verifyEqual(out.eventInfo.baselinePeriod, 1);
        end
    end
end

function iVerifyCarried(testCase, out, expected)
testCase.assertTrue(isstruct(out) && isfield(out, 'eventInfo'), 'UMT with eventInfo expected');
testCase.verifyEqual(out.eventInfo.selected, expected.selected);
testCase.verifyEqual(out.eventInfo.durationSec, expected.durationSec);
testCase.verifyEqual(out.eventInfo.baselinePeriod, expected.baselinePeriod);
testCase.verifyEqual(double(out.eventInfo.eventID), double(expected.eventID));
end
