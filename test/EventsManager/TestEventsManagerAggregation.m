classdef TestEventsManagerAggregation < matlab.unittest.TestCase
    %TESTEVENTSMANAGERAGGREGATION Centralized aggregation over events (.dat header Phase 8c).
    %
    %   EventsManager.conditionAggregationPlan and reduceByCondition carry
    %   the rules every event-reducing function follows: ignored instances
    %   excluded, conditions in first-appearance order, a NaN slice for a
    %   condition whose instances are all ignored, already aggregated input
    %   refused. Synthetic eventInfo only; no files.

    methods (Test)
        function planExcludesIgnoredAndOrdersByFirstAppearance(testCase)
            % Instances: B1, A1, B2, A2 (A2 ignored), C1 (ignored).
            ev = iInstances([2 1 2 1 3], [1 1 2 2 1], ["B" "A" "B" "A" "C"], ...
                logical([1 1 1 0 0]), [0.5 0.4 0.7 0.4 0.2]);
            plan = EventsManager.conditionAggregationPlan(ev);

            testCase.verifyEqual(plan.conditionID, [2; 1; 3], 'first appearance, not ID order');
            testCase.verifyEqual(plan.conditionName, ["B"; "A"; "C"]);
            testCase.verifyEqual(plan.instanceIdx, {[1 3]; 2; zeros(1, 0)});
            testCase.verifyEqual(plan.nInstances, [2; 1; 0]);
            testCase.verifyEqual(plan.durationSec, [0.6; 0.4; NaN], 'AbsTol', 1e-12);

            out = plan.eventInfoOut;
            testCase.verifyEqual(out.eventAxisMode, 'aggregated_repetitions');
            testCase.verifyEqual(out.repetitionIndex, zeros(3, 1));
            testCase.verifyEqual(out.selected, [true; true; false]);
            testCase.verifyEqual(out.nInstances, [2; 1; 0]);
            testCase.verifyEqual(out.baselinePeriod, 1.5);
            testCase.verifyClass(out.eventID, class(ev.eventID));
        end

        function reduceUsesSelectedInstancesAndNaNForEmptyConditions(testCase)
            ev = iInstances([2 1 2 1 3], [1 1 2 2 1], ["B" "A" "B" "A" "C"], ...
                logical([1 1 1 0 0]), nan(1, 5));
            plan = EventsManager.conditionAggregationPlan(ev);
            data = reshape(single(1:2*2*3*5), 2, 2, 3, 5);

            out = EventsManager.reduceByCondition(data, plan, @(x) mean(x, 4), 4);

            testCase.verifySize(out, [2 2 3 3]);
            testCase.verifyEqual(out(:, :, :, 1), mean(data(:, :, :, [1 3]), 4));
            testCase.verifyEqual(out(:, :, :, 2), data(:, :, :, 2), 'A2 is ignored');
            testCase.verifyTrue(all(isnan(out(:, :, :, 3)), 'all'), 'C is all ignored');
            testCase.verifyClass(out, 'single');
        end

        function slabReductionEqualsWholeArray(testCase)
            ev = iInstances([1 2 1 2], [1 1 2 2], ["A" "B" "A" "B"], true(1, 4), nan(1, 4));
            plan = EventsManager.conditionAggregationPlan(ev);
            data = rand(4, 6, 5, 4, 'single');
            whole = EventsManager.reduceByCondition(data, plan, @(x) median(x, 4), 4);
            slabs = cat(2, ...
                EventsManager.reduceByCondition(data(:, 1:2, :, :), plan, @(x) median(x, 4), 4), ...
                EventsManager.reduceByCondition(data(:, 3:6, :, :), plan, @(x) median(x, 4), 4));
            testCase.verifyEqual(slabs, whole);
        end

        function missingSelectedMeansAllSelected(testCase)
            % eventInfo saved before Phase 8b has no 'selected' field.
            ev = iInstances([1 1 2], [1 2 1], ["A" "A" "B"], true(1, 3), nan(1, 3));
            ev = rmfield(ev, {'selected', 'durationSec'});
            plan = EventsManager.conditionAggregationPlan(ev);
            testCase.verifyEqual(plan.nInstances, [2; 1]);
            testCase.verifyTrue(all(isnan(plan.durationSec)));
        end

        function namesFromEventNameList(testCase)
            % exportEventInfo carries eventNameList instead of eventName.
            ev = struct('eventID', uint16([2; 1]), 'repetitionIndex', uint16([1; 1]), ...
                'eventNameList', {{'A', 'B'}}, 'selected', true(2, 1));
            plan = EventsManager.conditionAggregationPlan(ev);
            testCase.verifyEqual(plan.conditionName, ["B"; "A"]);
        end

        function invalidInputsAreRefused(testCase)
            ev = iInstances([1 2], [1 1], ["A" "B"], true(1, 2), nan(1, 2));
            plan = EventsManager.conditionAggregationPlan(ev);
            testCase.verifyError(@() EventsManager.reduceByCondition(zeros(2, 2, 3, 3), plan, ...
                @(x) mean(x, 4), 4), 'Umitoolbox:EventsManager:eventAxisMismatch');
            planTwo = EventsManager.conditionAggregationPlan( ...
                iInstances([1 1 2], [1 2 1], ["A" "A" "B"], true(1, 3), nan(1, 3)));
            testCase.verifyError(@() EventsManager.reduceByCondition(zeros(2, 2, 3, 3), planTwo, ...
                @(x) x, 4), 'Umitoolbox:EventsManager:invalidReduction');

            aggregated = plan.eventInfoOut;
            testCase.verifyError(@() EventsManager.conditionAggregationPlan(aggregated), ...
                'Umitoolbox:EventsManager:alreadyAggregated');

            allIgnored = iInstances([1 2], [1 1], ["A" "B"], false(1, 2), nan(1, 2));
            testCase.verifyError(@() EventsManager.reduceByCondition(zeros(2, 2, 3, 2), ...
                EventsManager.conditionAggregationPlan(allIgnored), @(x) mean(x, 4), 4), ...
                'Umitoolbox:EventsManager:noSelectedInstances');
        end

        function aggregatedEventInfoIsAValidUMT(testCase)
            ev = iInstances([1 2 1], [1 1 2], ["A" "B" "A"], logical([1 0 1]), [0.5 0.5 0.5]);
            plan = EventsManager.conditionAggregationPlan(ev);
            umt = genUMTStruct(zeros(2, 2, 3, 2, 'single'), 'kind', 'image', ...
                'entryName', 'main', 'dimNames', {'Y', 'X', 'T', 'E'});
            umt = appendUMTEventInfo(umt, 'eventInfo', plan.eventInfoOut, 'overwrite', true);
            validateUMTStruct(umt, 'requireEventInfo', true);
            testCase.verifyEqual(umt.eventInfo.nInstances, [2; 0]);
            testCase.verifyEqual(umt.eventInfo.selected, [true; false]);
            testCase.verifyEqual(umt.eventInfo.baselinePeriod, 1.5);
        end
    end
end

function ev = iInstances(ids, reps, names, selected, durations)
ev = struct();
ev.eventID = uint16(ids(:));
ev.repetitionIndex = reps(:);
ev.eventName = names(:);
ev.eventAxisMode = 'instances';
ev.selected = selected(:);
ev.durationSec = durations(:);
ev.baselinePeriod = 1.5;
end
