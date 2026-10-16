classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestPipelineManagerParameterContainer < matlab.unittest.TestCase
    %TESTPIPELINEMANAGERPARAMETERCONTAINER Embedded and dialog parameter editors.
    %
    %   PipelineManager.setParameters(stepTag, container) builds the editable
    %   parameter controls inside a caller-supplied container (used by the
    %   PipelineManagerTool details panel); without a container it opens a
    %   dedicated uifigure dialog. Both share the same control builders and
    %   the same validate-and-store path.
    %
    %   The fixture pmParamEditor declares one parameter of every editable
    %   kind (numeric with limits, dropdown, multi-select, mixed text/number,
    %   logical). Controls are driven through their callbacks, so no display
    %   interaction is needed.

    properties
        ProjectRoot char
        SaveFolder char
        Figure
        Panel
        AppliedLog cell = {}
    end

    methods (TestClassSetup)
        function resolveProjectRoot(testCase)
            testFolder = fileparts(mfilename('fullpath'));
            testCase.ProjectRoot = extractBefore(testFolder, [filesep 'test']);
        end
    end

    methods (TestMethodSetup)
        function createFixtures(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.SaveFolder = fullfile(fixture.Folder, 'SaveFolder');
            mkdir(testCase.SaveFolder);

            AcqInfoStream = struct('Width', 3, 'Height', 2, 'Length', 4, ...
                'FrameRateHz', 10, 'ExposureMsec', 1);
            AcqInfoStream = appendImportedChannelInfo(AcqInfoStream, ...
                struct('DatFile', 'a.dat', 'Length', 4, 'FrameRateHz', 10));
            save(fullfile(testCase.SaveFolder, 'AcqInfos.mat'), 'AcqInfoStream');
            saveData(fullfile(testCase.SaveFolder, 'a.dat'), ...
                single(reshape(1:24, [2 3 4])), 'DimNames', {'Y', 'X', 'T'}, 'FrameRateHz', 10);

            parserFolder = iCreateFixtureCategory(testCase.ProjectRoot);
            testCase.addTeardown(@() removeAnalysisFixtureCategory(parserFolder, testCase.ProjectRoot));

            testCase.AppliedLog = {};
        end

        function createContainer(testCase)
            testCase.Figure = uifigure('Visible', 'off');
            testCase.addTeardown(@() delete(testCase.Figure));
            testCase.Panel = uipanel(testCase.Figure);
        end
    end

    methods (Test)
        function embeddedEditorBuildsEveryControlKind(testCase)
            [pm, tag] = testCase.newPipeline();
            pm.setParameters(tag, testCase.Panel);

            panel = testCase.Panel;
            testCase.verifyNumElements(findall(panel, 'Type', 'uinumericeditfield'), 1);
            testCase.verifyNumElements(findall(panel, 'Type', 'uidropdown'), 1);
            testCase.verifyNumElements(findall(panel, 'Type', 'uieditfield'), 1);

            checks = findall(panel, 'Type', 'uicheckbox');
            testCase.verifyEqual(sort({checks.Text}), {'', 'a', 'b', 'c'});

            applyBtn = iFindButton(panel, 'Apply');
            testCase.verifyNotEmpty(applyBtn);
            testCase.verifyNotEmpty(iFindButton(panel, 'Reset'));

            statusLabel = iFindStatusLabel(panel);
            testCase.verifyEqual(statusLabel.Text, 'Parameters: not set (defaults will be used)');
        end

        function applyStoresValuesFlagsNodeAndNotifies(testCase)
            [pm, tag] = testCase.newPipeline();
            pm.setParameters(tag, testCase.Panel, 'AppliedFcn', @(t) testCase.logApplied(t));
            panel = testCase.Panel;

            gain = findall(panel, 'Type', 'uinumericeditfield');
            mode = findall(panel, 'Type', 'uidropdown');
            exposure = findall(panel, 'Type', 'uieditfield');
            normalize = iFindCheckbox(panel, '');
            checkB = iFindCheckbox(panel, 'b');
            gain.Value = 4;
            mode.Value = 'slow';
            exposure.Value = '2.5';
            normalize.Value = true;
            checkB.Value = true;

            gain.ValueChangedFcn(gain, []);
            applyBtn = iFindButton(panel, 'Apply');
            testCase.verifyEqual(applyBtn.Text, 'Apply *');

            applyBtn.ButtonPushedFcn(applyBtn, []);

            node = pm.nodes(strcmp({pm.nodes.name}, tag));
            values = iParameterValues(node);
            testCase.verifyEqual(values.Gain, 4);
            testCase.verifyEqual(values.Mode, 'slow');
            testCase.verifyEqual(values.Exposure, '2.5');
            testCase.verifyTrue(values.Normalize);
            testCase.verifyEqual(sort(values.Bands(:)), {'a'; 'b'});
            testCase.verifyTrue(node.runtime.b_paramsSet);
            testCase.verifyEqual(testCase.AppliedLog, {tag});
            testCase.verifyEqual(applyBtn.Text, 'Apply');

            statusLabel = iFindStatusLabel(panel);
            testCase.verifyEqual(statusLabel.Text, 'Parameters: set');
        end

        function invalidValueLeavesNodeUnchanged(testCase)
            [pm, tag] = testCase.newPipeline();
            pm.setParameters(tag, testCase.Panel, 'AppliedFcn', @(t) testCase.logApplied(t));
            panel = testCase.Panel;
            before = pm.nodes(strcmp({pm.nodes.name}, tag));

            gain = findall(panel, 'Type', 'uinumericeditfield');
            exposure = findall(panel, 'Type', 'uieditfield');
            gain.Value = 7;
            exposure.Value = 'not a number';   % neither 'auto' nor numeric

            applyBtn = iFindButton(panel, 'Apply');
            applyBtn.ButtonPushedFcn(applyBtn, []);

            after = pm.nodes(strcmp({pm.nodes.name}, tag));
            testCase.verifyEqual(after.info.parameters, before.info.parameters);
            testCase.verifyFalse(after.runtime.b_paramsSet);
            testCase.verifyEmpty(testCase.AppliedLog);
        end

        function resetRestoresDefaultsInControls(testCase)
            [pm, tag] = testCase.newPipeline();
            pm.setParameters(tag, testCase.Panel);
            panel = testCase.Panel;

            gain = findall(panel, 'Type', 'uinumericeditfield');
            mode = findall(panel, 'Type', 'uidropdown');
            checkC = iFindCheckbox(panel, 'c');
            gain.Value = 9;
            mode.Value = 'slow';
            checkC.Value = true;

            resetBtn = iFindButton(panel, 'Reset');
            resetBtn.ButtonPushedFcn(resetBtn, []);

            testCase.verifyEqual(gain.Value, 1);
            testCase.verifyEqual(mode.Value, 'fast');
            testCase.verifyFalse(checkC.Value);
            testCase.verifyTrue(iFindCheckbox(panel, 'a').Value);
            testCase.verifyEqual(iFindButton(panel, 'Apply').Text, 'Apply *');
        end

        function stepsWithoutParametersShowLabelAndRebuildClearsContainer(testCase)
            [pm, tag] = testCase.newPipeline();
            folderNode = pm.nodes(strcmp({pm.nodes.kind}, 'folder'));
            testCase.verifyNotEmpty(folderNode);

            pm.setParameters(tag, testCase.Panel);
            testCase.verifyNotEmpty(findall(testCase.Panel, 'Type', 'uinumericeditfield'));

            pm.setParameters(folderNode(1).name, testCase.Panel);
            testCase.verifyEmpty(findall(testCase.Panel, 'Type', 'uinumericeditfield'));
            testCase.verifyEmpty(iFindButton(testCase.Panel, 'Apply'));
            label = findall(testCase.Panel, 'Type', 'uilabel');
            testCase.verifyEqual({label.Text}, {'No parameters'});

            pm.setParameters('', testCase.Panel);
            label = findall(testCase.Panel, 'Type', 'uilabel');
            testCase.verifyEqual({label.Text}, {'No step selected'});
        end

        function unsupportedParameterTypeShowsMessageInsteadOfThrowing(testCase)
            % A type the editor cannot render (e.g. struct) must not break
            % selecting the step: the container explains instead.
            pm = PipelineManager(testCase.SaveFolder, '', testCase.ProjectRoot);
            tag = pm.addStep('pmParamUnsupported', 'input', 'a.dat');

            pm.setParameters(tag, testCase.Panel);
            labels = findall(testCase.Panel, 'Type', 'uilabel');
            testCase.verifyTrue(startsWith(labels(1).Text, 'Parameters cannot be edited here'));
            testCase.verifyEmpty(iFindButton(testCase.Panel, 'Apply'));
        end

        function invalidContainerIsRejected(testCase)
            [pm, tag] = testCase.newPipeline();
            deadPanel = uipanel(testCase.Figure);
            delete(deadPanel);
            testCase.verifyError(@() pm.setParameters(tag, deadPanel), ...
                'PipelineManager:setParameters:InvalidContainer');
        end

        function dialogWithoutContainerStillStoresValuesOnOk(testCase)
            [pm, tag] = testCase.newPipeline();
            t = iStartDialogDriver(['Parameters: ' tag], true);
            testCase.addTeardown(@() iStopTimer(t));

            state = pm.setParameters(tag);

            testCase.verifyTrue(state);
            node = pm.nodes(strcmp({pm.nodes.name}, tag));
            values = iParameterValues(node);
            testCase.verifyEqual(values.Gain, 5);
            testCase.verifyTrue(node.runtime.b_paramsSet);
        end

        function dialogClosedWithoutOkChangesNothing(testCase)
            [pm, tag] = testCase.newPipeline();
            t = iStartDialogDriver(['Parameters: ' tag], false);
            testCase.addTeardown(@() iStopTimer(t));

            state = pm.setParameters(tag);

            testCase.verifyFalse(state);
            node = pm.nodes(strcmp({pm.nodes.name}, tag));
            testCase.verifyEqual(iParameterValues(node).Gain, 1);
            testCase.verifyFalse(node.runtime.b_paramsSet);
        end
    end

    methods
        function logApplied(testCase, tag)
            testCase.AppliedLog{end+1} = tag;
        end

        function [pm, tag] = newPipeline(testCase)
            pm = PipelineManager(testCase.SaveFolder, '', testCase.ProjectRoot);
            tag = pm.addStep('pmParamEditor', 'input', 'a.dat');
        end
    end
end

function values = iParameterValues(node)
params = node.info.parameters;
values = struct();
for k = 1:numel(params)
    values.(params(k).name) = params(k).value;
end
end

function btn = iFindButton(container, textStart)
buttons = findall(container, 'Type', 'uibutton');
if isempty(buttons)
    btn = buttons;
    return
end
btn = buttons(startsWith({buttons.Text}, textStart));
end

function lbl = iFindStatusLabel(container)
labels = findall(container, 'Type', 'uilabel');
lbl = labels(startsWith({labels.Text}, 'Parameters:'));
end

function chk = iFindCheckbox(container, label)
checks = findall(container, 'Type', 'uicheckbox');
chk = checks(strcmp({checks.Text}, label));
end

function t = iStartDialogDriver(dialogName, pressOk)
t = timer('ExecutionMode', 'fixedSpacing', 'Period', 0.3, 'StartDelay', 0.3, ...
    'TasksToExecute', 100, ...
    'TimerFcn', @(src, ~) iDriveDialog(src, dialogName, pressOk));
start(t);
end

function iStopTimer(t)
if isvalid(t)
    stop(t);
    delete(t);
end
end

function iDriveDialog(src, dialogName, pressOk)
% Polls while setParameters builds and then waits in uiwait. Once the dialog
% is complete (its OK button exists) edit one value and press OK, or close the
% dialog. The dialog is always released so a failure cannot hang the test.
fig = findall(0, 'Type', 'figure', 'Name', dialogName);
if isempty(fig)
    return
end
fig = fig(1);
ok = findall(fig, 'Type', 'uibutton', 'Text', 'OK');
if isempty(ok)
    if src.TasksExecuted >= src.TasksToExecute
        delete(fig);
    end
    return
end
stop(src);
try
    if pressOk
        gain = findall(fig, 'Type', 'uinumericeditfield');
        gain.Value = 5;
        ok.ButtonPushedFcn(ok, []);
    else
        delete(fig);
    end
catch
    if isvalid(fig)
        delete(fig);
    end
end
end

function parserFolder = iCreateFixtureCategory(projectRoot)
% PipelineManager discovers functions under Analysis/: copy the fixtures into
% a temporary category there (createAnalysisFixtureCategory).
fixtureFolder = fullfile(fileparts(mfilename('fullpath')), 'fixtures');
names = {'pmParamEditor', 'pmParamUnsupported'};
files = struct('name', names, 'source', fullfile(fixtureFolder, strcat(names, '.m')));
parserFolder = createAnalysisFixtureCategory(projectRoot, files);
end
