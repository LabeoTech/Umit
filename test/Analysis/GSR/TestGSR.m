classdef TestGSR < matlab.unittest.TestCase
    % TESTGSR Unit tests for GSR.
    %
    % Required fixture files located in the same folder as this test class:
    %   - green.dat
    %   - AcqInfos.mat
    %   - DataParams.mat
    %
    % Test strategy:
    %   - Use STANDARD MODE as the reference.
    %   - Compare LOW-RAM MODE output to the STANDARD MODE reference.
    %   - If exact equality fails, accept the result only if std(diff) < 1e-4
    %     over finite entries and the NaN pattern matches.

    properties (Access = private)
        TempFolder char = ''
        FixtureFolder char = ''
    end

    methods (TestClassSetup)
        function resolveFixtureFolder(testCase)
            testCase.FixtureFolder = fileparts(mfilename('fullpath'));

            testCase.assertTrue( ...
                isfile(fullfile(testCase.FixtureFolder, 'green.dat')), ...
                'Fixture file "green.dat" not found.');

            testCase.assertTrue( ...
                isfile(fullfile(testCase.FixtureFolder, 'AcqInfos.mat')), ...
                'Fixture file "AcqInfos.mat" not found.');

            testCase.assertTrue( ...
                isfile(fullfile(testCase.FixtureFolder, 'DataParams_original.mat')), ...
                'Fixture file "DataParams_original.mat" not found.');
        end
    end

    methods (TestMethodSetup)
        function createTempFolder(testCase)
            fx = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.TempFolder = fx.Folder;
        end
    end

    methods (TestMethodTeardown)
        function cleanupStrayFixtureFiles(testCase)
            strayFile = fullfile(testCase.FixtureFolder, 'DataParams.mat');
            if isfile(strayFile)
                delete(strayFile);
            end
        end
    end

    methods (Test)
        function testPipelineInfoContract(testCase)
            info = GSR('pipelineInfo');

            testCase.verifyTrue(isstruct(info), ...
                'GSR(''pipelineInfo'') must return a struct.');

            testCase.verifyTrue(isfield(info, 'outputs') && ~isempty(info.outputs), ...
                'pipelineInfo must declare at least one output.');

            outIdx = find(strcmp({info.outputs.name}, 'outData'), 1, 'first');
            testCase.assertNotEmpty(outIdx, ...
                'pipelineInfo must declare the output "outData".');

            testCase.verifyEqual(info.outputs(outIdx).outputMode, 'data');
            testCase.verifyTrue(info.outputs(outIdx).isData);

            paramNames = string.empty(0,1);

            if isfield(info, 'parameters') && ~isempty(info.parameters)
                tmpNames = {info.parameters.name};
                paramNames = [paramNames; string(tmpNames(:))]; %#ok<AGROW>
            end

            if isfield(info, 'arguments') && ~isempty(info.arguments)
                argKinds = {info.arguments.kind};
                argMask = strcmpi(argKinds, 'parameter');
                if any(argMask)
                    tmpNames = {info.arguments(argMask).name};
                    paramNames = [paramNames; string(tmpNames(:))]; %#ok<AGROW>
                end
            end

            paramNames = unique(paramNames, 'stable');

            testCase.verifyTrue(any(strcmpi(paramNames, "UseMask")), ...
                'pipelineInfo must expose "UseMask" as a parameter.');
        end

        function testStandardVsLowRAMWithoutMask(testCase)
            workFolder = localPrepareWorkingFolder(testCase, 'NoMask', false);

            dataStdIn = localReadGreenData(workFolder);
            outStd = GSR(dataStdIn, workFolder, 'UseMask', false);

            outFileLow = GSR('green.dat', workFolder, 'UseMask', false);
            testCase.verifyTrue(isfile(outFileLow), ...
                'LOW-RAM mode did not create the corrected output file.');

            outLow = localReadDatWithAcqInfo(outFileLow, workFolder);

            localVerifyHybridEquality(testCase, outLow, outStd, ...
                'LOW-RAM without mask vs STANDARD without mask');
        end

        function testStandardVsLowRAMWithMask(testCase)
            workFolder = localPrepareWorkingFolder(testCase, 'WithMask', true);

            dataStdIn = localReadGreenData(workFolder);
            outStd = GSR(dataStdIn, workFolder, 'UseMask', true);

            outFileLow = GSR('green.dat', workFolder, 'UseMask', true);
            testCase.verifyTrue(isfile(outFileLow), ...
                'LOW-RAM mode did not create the corrected output file.');

            outLow = localReadDatWithAcqInfo(outFileLow, workFolder);

            localVerifyHybridEquality(testCase, outLow, outStd, ...
                'LOW-RAM with mask vs STANDARD with mask');
        end

        function testMissingDataParamsFallsBackToNoMaskReference(testCase)
            workFolder = localPrepareWorkingFolder(testCase, 'MissingDataParams', false);

            dataStdIn = localReadGreenData(workFolder);
            refStd = GSR(dataStdIn, workFolder, 'UseMask', false);

            [outStdMasked, warnMsgStd, warnIdStd] = localRunWithWarningCapture( ...
                @() GSR(dataStdIn, workFolder, 'UseMask', true));

            localVerifyHybridEquality(testCase, outStdMasked, refStd, ...
                'STANDARD fallback for missing DataParams');

            testCase.verifyTrue(strcmp(warnIdStd, 'Umitoolbox:GSR:MissingDataParams') || ...
                contains(warnMsgStd, 'DataParams.mat not found'), ...
                'Missing DataParams should raise the expected fallback warning.');

            [outFileLow, warnMsgLow, warnIdLow] = localRunWithWarningCapture( ...
                @() GSR('green.dat', workFolder, 'UseMask', true));

            outLowMasked = localReadDatWithAcqInfo(outFileLow, workFolder);

            localVerifyHybridEquality(testCase, outLowMasked, refStd, ...
                'LOW-RAM fallback for missing DataParams');

            testCase.verifyTrue(strcmp(warnIdLow, 'Umitoolbox:GSR:MissingDataParams') || ...
                contains(warnMsgLow, 'DataParams.mat not found'), ...
                'Missing DataParams should raise the expected fallback warning in LOW-RAM mode.');
        end

        function testMissingMaskFieldFallsBackToNoMaskReference(testCase)
            workFolder = localPrepareWorkingFolder(testCase, 'MissingMaskField', true);
            localModifyDataParamsMask(workFolder, 'missingfield');

            dataStdIn = localReadGreenData(workFolder);
            refStd = GSR(dataStdIn, workFolder, 'UseMask', false);

            [outStdMasked, warnMsgStd, warnIdStd] = localRunWithWarningCapture( ...
                @() GSR(dataStdIn, workFolder, 'UseMask', true));

            localVerifyHybridEquality(testCase, outStdMasked, refStd, ...
                'STANDARD fallback for missing mask field');

            testCase.verifyTrue(strcmp(warnIdStd, 'Umitoolbox:GSR:MissingMask') || ...
                contains(warnMsgStd, 'Logical mask was not set'), ...
                'Missing mask field should raise the expected fallback warning.');

            [outFileLow, warnMsgLow, warnIdLow] = localRunWithWarningCapture( ...
                @() GSR('green.dat', workFolder, 'UseMask', true));

            outLowMasked = localReadDatWithAcqInfo(outFileLow, workFolder);

            localVerifyHybridEquality(testCase, outLowMasked, refStd, ...
                'LOW-RAM fallback for missing mask field');

            testCase.verifyTrue(strcmp(warnIdLow, 'Umitoolbox:GSR:MissingMask') || ...
                contains(warnMsgLow, 'Logical mask was not set'), ...
                'Missing mask field should raise the expected fallback warning in LOW-RAM mode.');
        end

        function testEmptyMaskFallsBackToNoMaskReference(testCase)
            workFolder = localPrepareWorkingFolder(testCase, 'EmptyMask', true);
            localModifyDataParamsMask(workFolder, 'empty');

            dataStdIn = localReadGreenData(workFolder);
            refStd = GSR(dataStdIn, workFolder, 'UseMask', false);

            [outStdMasked, warnMsgStd, warnIdStd] = localRunWithWarningCapture( ...
                @() GSR(dataStdIn, workFolder, 'UseMask', true));

            localVerifyHybridEquality(testCase, outStdMasked, refStd, ...
                'STANDARD fallback for empty mask');

            testCase.verifyTrue(strcmp(warnIdStd, 'Umitoolbox:GSR:TrivialMask') || ...
                contains(warnMsgStd, 'empty or all true'), ...
                'Empty mask should raise the expected fallback warning.');

            [outFileLow, warnMsgLow, warnIdLow] = localRunWithWarningCapture( ...
                @() GSR('green.dat', workFolder, 'UseMask', true));

            outLowMasked = localReadDatWithAcqInfo(outFileLow, workFolder);

            localVerifyHybridEquality(testCase, outLowMasked, refStd, ...
                'LOW-RAM fallback for empty mask');

            testCase.verifyTrue(strcmp(warnIdLow, 'Umitoolbox:GSR:TrivialMask') || ...
                contains(warnMsgLow, 'empty or all true'), ...
                'Empty mask should raise the expected fallback warning in LOW-RAM mode.');
        end

        function testAllTrueMaskFallsBackToNoMaskReference(testCase)
            workFolder = localPrepareWorkingFolder(testCase, 'AllTrueMask', true);
            localModifyDataParamsMask(workFolder, 'alltrue');

            dataStdIn = localReadGreenData(workFolder);
            refStd = GSR(dataStdIn, workFolder, 'UseMask', false);

            [outStdMasked, warnMsgStd, warnIdStd] = localRunWithWarningCapture( ...
                @() GSR(dataStdIn, workFolder, 'UseMask', true));

            localVerifyHybridEquality(testCase, outStdMasked, refStd, ...
                'STANDARD fallback for all-true mask');

            testCase.verifyTrue(strcmp(warnIdStd, 'Umitoolbox:GSR:TrivialMask') || ...
                contains(warnMsgStd, 'empty or all true'), ...
                'All-true mask should raise the expected fallback warning.');

            [outFileLow, warnMsgLow, warnIdLow] = localRunWithWarningCapture( ...
                @() GSR('green.dat', workFolder, 'UseMask', true));

            outLowMasked = localReadDatWithAcqInfo(outFileLow, workFolder);

            localVerifyHybridEquality(testCase, outLowMasked, refStd, ...
                'LOW-RAM fallback for all-true mask');

            testCase.verifyTrue(strcmp(warnIdLow, 'Umitoolbox:GSR:TrivialMask') || ...
                contains(warnMsgLow, 'empty or all true'), ...
                'All-true mask should raise the expected fallback warning in LOW-RAM mode.');
        end

        function testWrongSizeMaskErrors(testCase)
            workFolder = localPrepareWorkingFolder(testCase, 'WrongSizeMask', true);
            localModifyDataParamsMask(workFolder, 'wrongsize');

            dataStdIn = localReadGreenData(workFolder);

            testCase.verifyError( ...
                @() GSR(dataStdIn, workFolder, 'UseMask', true), ...
                'Umitoolbox:GSR:InvalidInput');

            testCase.verifyError( ...
                @() GSR('green.dat', workFolder, 'UseMask', true), ...
                'Umitoolbox:GSR:InvalidInput');
        end

        function testEventSplitArrayRegressesEachTrialIndependently(testCase)
            workFolder = localPrepareWorkingFolder(testCase, 'EventSplitArray', false);

            dataYXTE = localBuildTrials(workFolder);
            % A pixel that is invalid (NaN) in the second trial only.
            dataYXTE(5, 5, :, 2) = NaN;

            outYXTE = GSR(dataYXTE, workFolder);

            testCase.verifyEqual(size(outYXTE), size(dataYXTE));
            for iE = 1:size(dataYXTE, 4)
                localVerifyHybridEquality(testCase, outYXTE(:,:,:,iE), ...
                    GSR(dataYXTE(:,:,:,iE), workFolder), ...
                    sprintf('YXTE trial %d vs the same trial regressed alone', iE));
            end
            testCase.verifyTrue(all(isnan(outYXTE(5, 5, :, 2)), 'all'));
            testCase.verifyFalse(any(isnan(outYXTE(5, 5, :, 1)), 'all'));
        end

        function testEventSplitDatMatchesArray(testCase)
            workFolder = localPrepareWorkingFolder(testCase, 'EventSplitDat', true);

            dataYXTE = localBuildTrials(workFolder);
            rate = loadMetaData(fullfile(workFolder, 'green.dat')).frameRateHz;
            if ~isfinite(rate)
                rate = 10;
            end
            inFile = fullfile(workFolder, 'byEvent.dat');
            saveData(inFile, dataYXTE, 'DimNames', {'Y','X','T','E'}, ...
                'FrameRateHz', rate);

            outFile = GSR(inFile, workFolder, 'UseMask', true);

            testCase.verifyEqual(outFile, fullfile(workFolder, 'GSR.dat'));
            testCase.verifyEqual(loadMetaData(outFile).dimNames, {'Y','X','T','E'});
            localVerifyHybridEquality(testCase, single(loadData(outFile)), ...
                GSR(dataYXTE, workFolder, 'UseMask', true), ...
                'YXTE LOW-RAM with mask vs STANDARD with mask');
        end

        function testRejectsUnsupportedDatLayout(testCase)
            workFolder = localPrepareWorkingFolder(testCase, 'BadLayout', false);

            dataYXT = localReadGreenData(workFolder);
            rate = loadMetaData(fullfile(workFolder, 'green.dat')).frameRateHz;
            if ~isfinite(rate)
                rate = 10;
            end
            inFile = fullfile(workFolder, 'yxe.dat');
            saveData(inFile, dataYXT(:,:,1:3), 'DimNames', {'Y','X','E'}, ...
                'FrameRateHz', rate);

            testCase.verifyError(@() GSR(inFile, workFolder), ...
                'Umitoolbox:GSR:unsupportedLayout');
        end

        % =================================================================
        % PipelineManager integration
        %
        % The tests above call GSR directly. Direct calls bypass everything
        % PipelineManager does around a node -- argument marshalling, the
        % array-vs-filename RAM decision, and leaf-output persistence -- so
        % they cannot see defects that only appear in that interaction.
        % GSR's DATA input declares supportsFile=true, so GSR must run under
        % every RAM scenario PipelineManager offers.
        % =================================================================
        function testPMExecutesInAllRamScenarios(testCase)
            for scenario = string(pmRamScenarioList('fileCapable'))
                workFolder = localPreparePMWorkingFolder(testCase, ...
                    char("PM_" + scenario));

                pm = buildPMForScenario(workFolder, 'GSR', char(scenario), ...
                    'Input', 'green.dat');
                pm.executePipeline('PrintSummary', false);

                testCase.verifyTrue( ...
                    isfile(fullfile(workFolder, 'GSR.dat')), ...
                    sprintf('Scenario "%s" did not produce the declared output GSR.dat.', ...
                    scenario));
            end
        end

        function testPMOutputSetIsRamScenarioInvariant(testCase)
            % Which files a node creates is part of its contract. RAM
            % availability is a resource decision, not a data-identity
            % decision, so the produced file set must be identical across
            % scenarios. GSR's low-RAM path writes GSR.dat itself, which can
            % collide with PipelineManager's own leaf-output save of the same
            % declared name and produce an extra renamed duplicate.
            prepareFcn = @() localPreparePMWorkingFolder(testCase, ...
                ['Invariance_' char(matlab.lang.makeValidName( ...
                char(java.util.UUID.randomUUID)))]);

            outputs = pmCollectScenarioOutputs(prepareFcn, 'GSR', ...
                'Input', 'green.dat', 'Extensions', {'.dat'});

            reference = outputs(1).files;
            for iScenario = 2:numel(outputs)
                testCase.verifyEqual(outputs(iScenario).files, reference, ...
                    sprintf(['Scenario "%s" produced %s but scenario "%s" produced %s. ' ...
                    'The set of files a node creates must not depend on the RAM scenario.'], ...
                    outputs(iScenario).scenario, ...
                    localFormatFileList(outputs(iScenario).files), ...
                    outputs(1).scenario, ...
                    localFormatFileList(reference)));
            end
        end

        function testPMDoesNotDuplicateDeclaredOutput(testCase)
            % Companion to the invariance test, stated as an absolute rather
            % than a comparison: a single GSR node must leave exactly one GSR
            % output behind in every scenario, never a GSR.dat plus a
            % collision-renamed GSR_1.dat carrying the same bytes.
            for scenario = string(pmRamScenarioList('fileCapable'))
                workFolder = localPreparePMWorkingFolder(testCase, ...
                    char("NoDup_" + scenario));

                pm = buildPMForScenario(workFolder, 'GSR', char(scenario), ...
                    'Input', 'green.dat');
                pm.executePipeline('PrintSummary', false);

                produced = dir(fullfile(workFolder, 'GSR*.dat'));

                testCase.verifyNumElements(produced, 1, ...
                    sprintf(['Scenario "%s" left %d GSR output file(s) (%s). ' ...
                    'A single node must produce a single declared output.'], ...
                    scenario, numel(produced), ...
                    localFormatFileList({produced.name})));
            end
        end
    end
end

% =========================================================================
% Local helpers
% =========================================================================
function workFolder = localPrepareWorkingFolder(testCase, caseName, includeDataParams)
%LOCALPREPAREWORKINGFOLDER Create a temporary working folder with fixtures.
%
% The original DataParams fixture is never modified in place. When needed,
% DataParams_original.mat is copied into the temporary test folder as
% DataParams.mat and only that copy is modified.

workFolder = fullfile(testCase.TempFolder, caseName);
mkdir(workFolder);

copyfile(fullfile(testCase.FixtureFolder, 'green.dat'), ...
    fullfile(workFolder, 'green.dat'));

copyfile(fullfile(testCase.FixtureFolder, 'AcqInfos.mat'), ...
    fullfile(workFolder, 'AcqInfos.mat'));

if includeDataParams
    copyfile(fullfile(testCase.FixtureFolder, 'DataParams_original.mat'), ...
        fullfile(workFolder, 'DataParams.mat'));
end
end

function workFolder = localPreparePMWorkingFolder(testCase, caseName)
%LOCALPREPAREPMWORKINGFOLDER Working folder accepted by PipelineManager.
%
% PipelineManager.executePipeline refuses any SaveFolder that
% isLegacySchemaFolder reports as legacy. This fixture's AcqInfos.mat has no
% ImportedChannels registry, because the direct-call tests never need one, so
% the copy used for pipeline tests is upgraded in place.

workFolder = localPrepareWorkingFolder(testCase, caseName, true);
ensurePMReadyAcqInfos(workFolder, 'green.dat');
end

function text = localFormatFileList(fileNames)
%LOCALFORMATFILELIST Render a file-name list for assertion messages.

if isempty(fileNames)
    text = '<none>';
    return
end
text = strjoin(cellstr(string(fileNames)), ', ');
end

function dataYXTE = localBuildTrials(workFolder)
%LOCALBUILDTRIALS Two event-like trials: the two halves of green.dat.

dataYXT = localReadGreenData(workFolder);
trialLen = floor(size(dataYXT, 3) / 2);
dataYXTE = cat(4, dataYXT(:,:,1:trialLen), dataYXT(:,:,trialLen + (1:trialLen)));
end

function data = localReadGreenData(workFolder)
%LOCALREADGREENDATA Read green.dat into memory using AcqInfos metadata.

data = localReadDatWithAcqInfo(fullfile(workFolder, 'green.dat'), workFolder);
end

function data = localReadDatWithAcqInfo(datPath, ~)
%LOCALREADDATWITHACQINFO Read a .dat file as single through loadData.
%
% .dat header Phase 4c-2a: GSR writes headered outputs, so the file is
% read through loadData rather than from byte 0 with AcqInfos sizes.

data = single(loadData(datPath));
end

function AcqInfoStream = localLoadAcqInfo(acqInfoPath)
%LOCALLOADACQINFO Load AcqInfoStream from AcqInfos.mat.

S = load(acqInfoPath);

if isfield(S, 'AcqInfoStream')
    AcqInfoStream = S.AcqInfoStream;
elseif isfield(S, 'AcqInfos')
    AcqInfoStream = S.AcqInfos;
else
    fn = fieldnames(S);
    assert(~isempty(fn), 'The MAT file "%s" does not contain any variables.', acqInfoPath);
    AcqInfoStream = S.(fn{1});
end
end

function localVerifyHybridEquality(testCase, actual, reference, contextMsg)
%LOCALVERIFYHYBRIDEQUALITY Compare outputs with exact equality first, then
% fallback to std(diff) < 1e-4 over finite entries.

testCase.verifyEqual(size(actual), size(reference), ...
    sprintf('%s: output size mismatch.', contextMsg));

if isequaln(actual, reference)
    return
end

actual = double(actual);
reference = double(reference);

testCase.verifyTrue(isequaln(isnan(actual), isnan(reference)), ...
    sprintf('%s: NaN pattern mismatch.', contextMsg));

validIdx = isfinite(actual) & isfinite(reference);

if ~any(validIdx(:))
    return
end

diffVals = actual(validIdx) - reference(validIdx);
diffStd = std(diffVals(:));

testCase.verifyTrue(diffStd < 1e-4, ...
    sprintf('%s: exact equality failed and std(diff) = %.6g >= 1e-4.', ...
    contextMsg, diffStd));
end

function [out, warnMsg, warnId] = localRunWithWarningCapture(funHandle)
%LOCALRUNWITHWARNINGCAPTURE Run a function and capture the last warning.

lastwarn('');
out = funHandle();
[warnMsg, warnId] = lastwarn;
lastwarn('');
end

function localModifyDataParamsMask(workFolder, variant)
%LOCALMODIFYDATAPARAMSMASK Modify the stored mask in DataParams.mat.
%
% This helper adapts to common DataParams layouts:
%   - DataParams.logical_mask
%   - DataParams.logicalMask
%   - DataParams.mask
%   - DataParams.mask.logical

dataParamsPath = fullfile(workFolder, 'DataParams.mat');
testStruct = load(dataParamsPath);

fn = fieldnames(testStruct);
assert(~isempty(fn), 'DataParams.mat is empty.');

if isfield(testStruct, 'DataParams')
    varName = 'DataParams';
elseif numel(fn) == 1 && isstruct(testStruct.(fn{1}))
    varName = fn{1};
else
    error('Could not identify the DataParams structure inside DataParams.mat.');
end

dataParams = testStruct.(varName);
AcqInfoStream = localLoadAcqInfo(fullfile(workFolder, 'AcqInfos.mat'));
frameSize = [AcqInfoStream.Height, AcqInfoStream.Width];

style = localDetectMaskStyle(dataParams);

switch lower(variant)
    case 'missingfield'
        dataParams = localRemoveMaskField(dataParams, style);

    case 'empty'
        dataParams = localSetMaskValue(dataParams, style, []);

    case 'alltrue'
        dataParams = localSetMaskValue(dataParams, style, true(frameSize));

    case 'wrongsize'
        badSize = [frameSize(1)+1, frameSize(2)];
        dataParams = localSetMaskValue(dataParams, style, true(badSize));

    otherwise
        error('Unknown DataParams mask variant "%s".', variant);
end

testStruct.(varName) = dataParams;
save(dataParamsPath, '-struct', 'testStruct');
end

function style = localDetectMaskStyle(dataParams)
%LOCALDETECTMASKSTYLE Detect how the logical mask is stored.

if isfield(dataParams, 'logical_mask')
    style = 'logical_mask';
elseif isfield(dataParams, 'logicalMask')
    style = 'logicalMask';
elseif isfield(dataParams, 'mask')
    if isstruct(dataParams.mask) && isfield(dataParams.mask, 'logical')
        style = 'mask.logical';
    else
        style = 'mask';
    end
else
    style = 'logical_mask';
end
end

function dataParams = localRemoveMaskField(dataParams, style)
%LOCALREMOVEMASKFIELD Remove the logical mask field.

switch style
    case 'logical_mask'
        if isfield(dataParams, 'logical_mask')
            dataParams = rmfield(dataParams, 'logical_mask');
        end

    case 'logicalMask'
        if isfield(dataParams, 'logicalMask')
            dataParams = rmfield(dataParams, 'logicalMask');
        end

    case 'mask'
        if isfield(dataParams, 'mask')
            dataParams = rmfield(dataParams, 'mask');
        end

    case 'mask.logical'
        if isfield(dataParams, 'mask') && isstruct(dataParams.mask) && ...
                isfield(dataParams.mask, 'logical')
            dataParams.mask = rmfield(dataParams.mask, 'logical');
        end

    otherwise
        error('Unknown mask style "%s".', style);
end
end

function dataParams = localSetMaskValue(dataParams, style, newMask)
%LOCALSETMASKVALUE Set the logical mask field.

switch style
    case 'logical_mask'
        dataParams.logical_mask = newMask;

    case 'logicalMask'
        dataParams.logicalMask = newMask;

    case 'mask'
        dataParams.mask = newMask;

    case 'mask.logical'
        if ~isfield(dataParams, 'mask') || ~isstruct(dataParams.mask)
            dataParams.mask = struct();
        end
        dataParams.mask.logical = newMask;

    otherwise
        error('Unknown mask style "%s".', style);
end
end