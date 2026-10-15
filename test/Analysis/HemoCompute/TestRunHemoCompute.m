classdef TestRunHemoCompute < matlab.unittest.TestCase
    %TESTRUNHEMOCOMPUTE Unit tests for HemoCompute.
    %
    %   These tests use synthetic AcqInfos-driven datasets to verify wrapper
    %   dispatch, output manifests, and timebase handling when imported
    %   channels have different lengths and frame rates.

    properties
        ProjectRoot
        TempFolder
        RigRoot
    end

    methods (TestMethodSetup)
        function createTempFolder(testCase)
            thisFile = mfilename('fullpath');
            testFolder = fileparts(thisFile);
            projectRoot = extractBefore(testFolder, [filesep 'test']);
            if isempty(projectRoot)
                projectRoot = fileparts(fileparts(testFolder));
            end

            testCase.ProjectRoot = char(projectRoot);
            addpath(genpath(testCase.ProjectRoot));

            testCase.TempFolder = fullfile(tempdir, ...
                ['TestRunHemoCompute_' char(java.util.UUID.randomUUID)]);
            mkdir(testCase.TempFolder);

            iCreateSyntheticHemoDataset(testCase.TempFolder, false);
            % Bound through the folder rig ledger; no library rig needed.
            iBindDatasetToRig(testCase.TempFolder, iConfiguredRigInfo());
            testCase.RigRoot = '';
        end
    end

    methods (TestMethodTeardown)
        function removeTempFolder(testCase)
            if ~isempty(testCase.TempFolder) && isfolder(testCase.TempFolder)
                try %#ok<TRYNC>
                    rmdir(testCase.TempFolder, 's');
                end
            end
            if ~isempty(testCase.RigRoot) && isfolder(testCase.RigRoot)
                rmdir(testCase.RigRoot, 's');
            end
        end
    end

    methods (Test)
        function testRawDatInput(testCase)
            out = HemoCompute(testCase.TempFolder, 'red.dat', ...
                'FilterSet', 'none', ...
                'b_normalize', true, ...
                'Illuminations', {'red','green'});

            testCase.verifyTrue(iscell(out));
            testCase.verifyEqual(numel(out), 2);
            testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, out{1})));
            testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, out{2})));
        end

        function testRepeatedIlluminationRawDatInput(testCase)
            repeatedFolder = fullfile(testCase.TempFolder, 'repeatedTimeline');
            mkdir(repeatedFolder);
            iCreateSyntheticHemoDataset(repeatedFolder, true);
            rigStore = iConfiguredRigInfo();
            iBindDatasetToRig(repeatedFolder, rigStore);

            out = HemoCompute(repeatedFolder, 'green.dat', ...
                'FilterSet', 'none', ...
                'b_normalize', true, ...
                'Illuminations', {'red','green'});

            testCase.verifyEqual(out, {'HbO.dat', 'HbR.dat'});

            outData = loadData(fullfile(repeatedFolder, out{1}));
            greenInfo = loadMetaData(fullfile(repeatedFolder, 'green.dat'));
            redInfo = loadMetaData(fullfile(repeatedFolder, 'red.dat'));

            testCase.verifyGreaterThan(redInfo.frameRateHz, greenInfo.frameRateHz);
            testCase.verifySize(outData, double(greenInfo.dimSizes(:).'));
        end

        function testEmptyIlluminationSelection(testCase)
            testCase.verifyError(@() HemoCompute(testCase.TempFolder, 'red.dat', ...
                'Illuminations', {}), ...
                'Umitoolbox:HemoCompute:InvalidInput');
        end

        function testSingleIlluminationSelectionRejectedByHemoCompute(testCase)
            testCase.verifyError(@() HemoCompute(testCase.TempFolder, 'red.dat', ...
                'Illuminations', {'red'}), ...
                'Umitoolbox:HemoCompute:InvalidInput');
        end

        function testInvalidFilterSet(testCase)
            testCase.verifyError(@() HemoCompute(testCase.TempFolder, 'red.dat', ...
                'FilterSet', 'badFilter'), ...
                'Umitoolbox:UMITRigStore:filterSetNotFound');
        end

        function testForwardPhysiologyParameters(testCase)
            out = HemoCompute(testCase.TempFolder, 'red.dat', ...
                'FilterSet', 'none', ...
                'b_normalize', true, ...
                'Illuminations', {'red','green'}, ...
                'HbT_concentration_uM', 90, ...
                'StO2perc', 55);

            testCase.verifyEqual(out, {'HbO.dat', 'HbR.dat'});
        end

        function testStO2percBoundariesAccepted(testCase)
            for saturation = [0 100]
                out = HemoCompute( ...
                    testCase.TempFolder, 'red.dat', ...
                    'FilterSet', 'none', ...
                    'b_normalize', true, ...
                    'Illuminations', {'red','green'}, ...
                    'StO2perc', saturation);
                testCase.verifyEqual(out, {'HbO.dat', 'HbR.dat'});
            end
        end

        function testStO2percOutsidePercentageRangeRejected(testCase)
            for saturation = [-0.1 100.1]
                testCase.verifyError(@() HemoCompute( ...
                    testCase.TempFolder, 'red.dat', ...
                    'StO2perc', saturation), ...
                    'Umitoolbox:HemoCompute:InvalidStO2perc');
            end
        end

        function testTwoOutputPipelineContract(testCase)
            [hbO, hbR] = HemoCompute( ...
                testCase.TempFolder, 'red.dat', ...
                'FilterSet', 'none', ...
                'Illuminations', {'red','green'});

            testCase.verifyEqual(hbO, 'HbO.dat');
            testCase.verifyEqual(hbR, 'HbR.dat');

            hbOData = loadData(fullfile(testCase.TempFolder, hbO));
            hbRData = loadData(fullfile(testCase.TempFolder, hbR));
            testCase.verifyNotEmpty(hbOData);
            testCase.verifyNotEmpty(hbRData);
            testCase.verifyEqual(size(hbOData), size(hbRData));
        end

        function testOutputsAreYXTFiles(testCase)
            out = HemoCompute(testCase.TempFolder, 'green.dat', ...
                'FilterSet', 'none', 'Illuminations', {'red','green'});

            for k = 1:2
                info = loadMetaData(fullfile(testCase.TempFolder, out{k}));
                testCase.verifyEqual(info.dimNames, {'Y','X','T'});
                testCase.verifyEqual(info.dataClass, 'single');
            end
        end

        function testRejectsUnsupportedInputs(testCase)
            sv = testCase.TempFolder;

            % Arrays are not an input form.
            testCase.verifyError(@() HemoCompute(sv, rand(4, 5, 6, 'single'), ...
                'Illuminations', {'red','green'}), ...
                'Umitoolbox:HemoCompute:UnsupportedInputType');
            testCase.verifyError(@() HemoCompute(sv, 'red.tif', ...
                'Illuminations', {'red','green'}), ...
                'Umitoolbox:HemoCompute:UnsupportedInputFile');
            testCase.verifyError(@() HemoCompute(sv, 'missing.dat', ...
                'Illuminations', {'red','green'}), ...
                'Umitoolbox:HemoCompute:FileNotFound');

            % The wired file must be one of the resolved illumination channels.
            testCase.verifyError(@() HemoCompute(sv, 'yellow.dat', ...
                'Illuminations', {'red','green'}), ...
                'Umitoolbox:HemoCompute:DataInputNotAChannel');
        end

        function testRejectsOtherLayouts(testCase)
            sv = testCase.TempFolder;
            info = loadMetaData(fullfile(sv, 'red.dat'));
            ny = double(info.dimSizes(1));
            nx = double(info.dimSizes(2));

            % The wired file itself is not Y-X-T.
            writeTestDat(fullfile(sv, 'red.dat'), rand(ny, nx, 3, 2, 'single'), 10, 5, ...
                'DimNames', {'Y','X','T','E'});
            testCase.verifyError(@() HemoCompute(sv, 'red.dat', ...
                'Illuminations', {'red','green'}), ...
                'Umitoolbox:HemoCompute:unsupportedLayout');

            % Another resolved channel file is not Y-X-T.
            testCase.verifyError(@() HemoCompute(sv, 'green.dat', ...
                'Illuminations', {'red','green'}), ...
                'Umitoolbox:HemoCompute:unsupportedLayout');
            testCase.verifyFalse(isfile(fullfile(sv, 'HbO.dat')));
        end

        function testPipelineInfo(testCase)
            info = HemoCompute('pipelineInfo');
            testCase.verifyEqual(info.name, 'HemoCompute');
            testCase.verifyFalse(info.legacyOpts);
            testCase.verifyEqual(info.outputs(1).type, {'ImageTimeSeries'});
            testCase.verifyEqual(info.outputs(2).type, {'ImageTimeSeries'});
            testCase.verifyEqual(info.outputs(1).defOutfilename, 'HbO.dat');
            testCase.verifyEqual(info.outputs(2).defOutfilename, 'HbR.dat');
            testCase.verifyEqual({info.outputs.outputMode}, {'data','data'});
            testCase.verifyTrue(all([info.outputs.isData]));
            testCase.verifyEqual({info.parameters.name}, ...
                {'FilterSet','b_normalize','Illuminations', ...
                'HbT_concentration_uM','StO2perc'});

            illuminationParam = info.parameters(strcmp( ...
                {info.parameters.name}, 'Illuminations'));
            testCase.verifyEqual(illuminationParam.type, 'char');
            testCase.verifyEqual(illuminationParam.allowed, ...
                {'red';'green';'yellow'});
            testCase.verifyEqual(illuminationParam.default, ...
                {'red','green','yellow'});
            saveFolderInput = info.inputs(strcmp({info.inputs.name}, 'SaveFolder'));
            testCase.verifyNumElements(saveFolderInput, 1);
            testCase.verifyEqual(saveFolderInput.type, {'SaveFolder'});
            testCase.verifyFalse(saveFolderInput.isData);
            saveFolderArgument = info.arguments( ...
                strcmp({info.arguments.name}, 'SaveFolder'));
            testCase.verifyEqual(saveFolderArgument.kind, 'input');

            dataInput = info.inputs(strcmp({info.inputs.name}, 'data'));
            testCase.verifyNumElements(dataInput, 1);
            testCase.verifyTrue(dataInput.isData);
            testCase.verifyTrue(dataInput.supportsFile);
            testCase.verifyEqual(dataInput.dataMode, 'file');

            saturationParam = info.parameters(strcmp( ...
                {info.parameters.name}, 'StO2perc'));
            testCase.verifyEqual(saturationParam.allowed, [0 100]);
            testCase.verifyTrue(info.supportsPreflight);
        end

        function testPipelineManagerRunsHemoPreflight(testCase)
            pm = PipelineManager({testCase.TempFolder});
            pm.addStep('HemoCompute', 'input', 'red.dat');
            [isValid, report] = pm.validateNodePreflight('throwOnError', false);
            testCase.verifyTrue(isValid, strjoin(report.errors, ' | '));
        end

        function testPipelineManagerExecutesTwoHemoOutputs(testCase)
            pm = PipelineManager({testCase.TempFolder});
            pm.b_skipSteps = false;
            pm.addStep('HemoCompute', 'input', 'red.dat');

            pm.executePipeline();

            testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, 'HbO.dat')));
            testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, 'HbR.dat')));
        end

        function testPipelineManagerAcceptsRamSafeBackingFiles(testCase)
            pm = PipelineManager({testCase.TempFolder});
            pm.ramMode = 'ramsafe';
            pm.b_skipSteps = false;
            pm.addStep('HemoCompute', 'input', 'red.dat');

            pm.executePipeline();

            testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, 'HbO.dat')));
            testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, 'HbR.dat')));
        end

        function testPipelineManagerExecutesAllRamScenarios(testCase)
            rigStore = iConfiguredRigInfo();
            prepareFcn = @() iPreparePMHemoFolder(testCase.TempFolder, rigStore);
            outputs = pmCollectScenarioOutputs(prepareFcn, 'HemoCompute', ...
                'Input', 'red.dat');

            expected = {'HbO.dat', 'HbR.dat'};
            for iScenario = 1:numel(outputs)
                testCase.verifyEqual(outputs(iScenario).files, expected, ...
                    sprintf('HemoCompute output mismatch under %s.', ...
                    outputs(iScenario).scenario));
            end
        end

        function testPipelineManagerReportsHemoPreflightFailure(testCase)
            delete(fullfile(testCase.TempFolder, 'green.dat'));
            pm = PipelineManager({testCase.TempFolder});
            pm.addStep('HemoCompute', 'input', 'red.dat');
            [isValid, report] = pm.validateNodePreflight('throwOnError', false);
            testCase.verifyFalse(isValid);
            testCase.verifyTrue(any(contains(string(report.errors), ...
                'Resolved channel file was not found')));
        end

        function testPipelineManagerInfersMissingImportedChannels(testCase)
            acqPath = fullfile(testCase.TempFolder, 'AcqInfos.mat');
            loaded = load(acqPath, 'AcqInfoStream');
            AcqInfoStream = rmfield(loaded.AcqInfoStream, 'ImportedChannels');
            % The fallback infers channels from the headered files.
            save(acqPath, 'AcqInfoStream');

            pm = PipelineManager({testCase.TempFolder});
            pm.b_skipSteps = false;
            pm.addStep('HemoCompute', 'input', 'red.dat');

            [isValid, report] = pm.validateNodePreflight('throwOnError', false);
            testCase.verifyTrue(isValid, strjoin(report.errors, ' | '));
            pm.executePipeline();
            testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, 'HbO.dat')));
            testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, 'HbR.dat')));

            persisted = load(acqPath, 'AcqInfoStream');
            testCase.verifyEqual(persisted.AcqInfoStream, AcqInfoStream);
        end

        function testPipelineManagerInfersMissingImportedChannelsRamSafe(testCase)
            acqPath = fullfile(testCase.TempFolder, 'AcqInfos.mat');
            loaded = load(acqPath, 'AcqInfoStream');
            AcqInfoStream = rmfield(loaded.AcqInfoStream, 'ImportedChannels');
            % The fallback infers channels from the headered files.
            save(acqPath, 'AcqInfoStream');

            pm = PipelineManager({testCase.TempFolder});
            pm.ramMode = 'ramsafe';
            pm.b_skipSteps = false;
            pm.addStep('HemoCompute', 'input', 'red.dat');

            [isValid, report] = pm.validateNodePreflight('throwOnError', false);
            testCase.verifyTrue(isValid, strjoin(report.errors, ' | '));
            pm.executePipeline();
            testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, 'HbO.dat')));
            testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, 'HbR.dat')));
        end

        function testExplicitImportedChannelsTakePrecedence(testCase)
            loaded = load(fullfile(testCase.TempFolder, 'AcqInfos.mat'), 'AcqInfoStream');
            AcqInfoStream = loaded.AcqInfoStream;
            AcqInfoStream.Illumination1 = struct( ...
                'Color', 'green', 'CamIdx', 2, 'FrameRateHz', 1);

            [resolved, report] = resolveImportedChannelFallback( ...
                AcqInfoStream, testCase.TempFolder, ...
                'RequiredChannels', {'red', 'green'});

            testCase.verifyEqual(report.Source, 'explicit');
            testCase.verifyEqual(resolved.ImportedChannels, AcqInfoStream.ImportedChannels);
        end

        function testPipelineManagerRejectsConflictingLegacyChannelEvidence(testCase)
            acqPath = fullfile(testCase.TempFolder, 'AcqInfos.mat');
            loaded = load(acqPath, 'AcqInfoStream');
            AcqInfoStream = rmfield(loaded.AcqInfoStream, 'ImportedChannels');
            % The fallback infers channels from the headered files.
            AcqInfoStream.Illumination1 = struct('Color', 'red', 'CamIdx', 1);
            AcqInfoStream.Illumination2 = struct('Color', 'red', 'CamIdx', 2);
            save(acqPath, 'AcqInfoStream');

            pm = PipelineManager({testCase.TempFolder});
            pm.addStep('HemoCompute', 'input', 'red.dat');
            [isValid, report] = pm.validateNodePreflight('throwOnError', false);

            testCase.verifyFalse(isValid);
            testCase.verifyTrue(any(contains(string(report.errors), ...
                'conflicting CamIdx values')));
        end

        function testPipelineManagerRejectsIncompleteFallbackEvidence(testCase)
            acqPath = fullfile(testCase.TempFolder, 'AcqInfos.mat');
            loaded = load(acqPath, 'AcqInfoStream');
            AcqInfoStream = rmfield(loaded.AcqInfoStream, 'ImportedChannels');
            % The fallback infers channels from the headered files.
            save(acqPath, 'AcqInfoStream');
            delete(fullfile(testCase.TempFolder, 'green.dat'));

            pm = PipelineManager({testCase.TempFolder});
            pm.addStep('HemoCompute', 'input', 'red.dat');
            [isValid, report] = pm.validateNodePreflight('throwOnError', false);

            testCase.verifyFalse(isValid);
            testCase.verifyTrue(any(contains(string(report.errors), ...
                'Cannot infer ImportedChannels')));
        end

        function testPipelineManagerRejectsAcqInfosBoundChannels(testCase)
            % Old-schema folder whose channel files are headerless without
            % sidecars (described only by AcqInfos.mat): unsupported since
            % .dat header Phase 5b, so the preflight refuses it.
            acqPath = fullfile(testCase.TempFolder, 'AcqInfos.mat');
            loaded = load(acqPath, 'AcqInfoStream');
            AcqInfoStream = rmfield(loaded.AcqInfoStream, 'ImportedChannels');
            iMakeAcqInfosBound(testCase.TempFolder);
            save(acqPath, 'AcqInfoStream');

            pm = PipelineManager({testCase.TempFolder});
            pm.addStep('HemoCompute', 'input', 'red.dat');
            [isValid, report] = pm.validateNodePreflight('throwOnError', false);

            testCase.verifyFalse(isValid);
            testCase.verifyNotEmpty(report.errors);
        end

        function testPipelineSkipsOnlyTheFolderFailingRigPreflight(testCase)
            %TESTPIPELINESKIPSONLYTHEFOLDERFAILINGRIGPREFLIGHT A rig-guard
            %preflight error stops only its folder, which is logged as Failed;
            %the batch still runs. When every folder fails, the run aborts.

            goodFolder = fullfile(testCase.TempFolder, 'boundFolder');
            badFolder = fullfile(testCase.TempFolder, 'unboundFolder');
            mkdir(goodFolder);
            mkdir(badFolder);
            iCreateSyntheticHemoDataset(goodFolder, false);
            iBindDatasetToRig(goodFolder, iConfiguredRigInfo());
            iCreateSyntheticHemoDataset(badFolder, false);

            pm = PipelineManager({goodFolder, badFolder});
            pm.b_skipSteps = false;
            pm.addStep('HemoCompute', 'input', 'red.dat');
            report = pm.diagnosePipeline('verbose', false);
            testCase.verifyTrue(report.isValid, strjoin(report.errors, ' | '));
            testCase.verifyTrue(any(contains(string(report.warnings), 'will be skipped')));

            testCase.verifyWarning(@() pm.executePipeline(), ...
                'PipelineManager:executePipeline:FolderFailed');
            testCase.verifyTrue(isfile(fullfile(goodFolder, 'HbO.dat')));
            testCase.verifyFalse(isfile(fullfile(badFolder, 'HbO.dat')));
            testCase.verifyEqual(string(pm.globalPipeLog.Status), ["Completed"; "Failed"]);
            testCase.verifySubstring(char(pm.globalPipeLog.Messages_short{2}), 'no rig ledger');
            testCase.verifyEqual(UMITRigLedger.read(goodFolder).state, 'BoundCommitted');
            testCase.verifyFalse(UMITRigLedger.exists(badFolder));

            pmAllBad = PipelineManager({badFolder});
            pmAllBad.addStep('HemoCompute', 'input', 'red.dat');
            testCase.verifyError(@() pmAllBad.executePipeline(), ...
                'PipelineManager:validateNodePreflight:Failed');
            testCase.verifyFalse(pmAllBad.diagnosePipeline('verbose', false).isValid);
        end

        function testAcknowledgedConflictForcesRecomputationWhenSkipping(testCase)
            %TESTACKNOWLEDGEDCONFLICTFORCESRECOMPUTATIONWHENSKIPPING Step
            %skipping never reuses a HemoCompute output computed before the
            %folder's latest acknowledged rig conflict; later outputs are reused.

            % HemoCompute -> spatialGaussFilt: re-running an unchanged pipeline
            % reuses the whole branch (HemoCompute is skipped). A lone
            % HemoCompute step is never skipped: its two outputs make the
            % leaf match ambiguous.
            folder = testCase.TempFolder;
            hboPath = fullfile(folder, 'HbO.dat');
            pm = PipelineManager({folder});
            pm.b_skipSteps = true;
            hemoTag = pm.addStep('HemoCompute', 'input', fullfile(folder, 'red.dat'));
            pm.addStep('spatialGaussFilt', 'input', sprintf('%s:HbO', hemoTag));
            runArgs = {'LeafOutputPolicy', 'saveLeaves', 'PrintSummary', false};

            pm.executePipeline(runArgs{:});
            hboBefore = loadData(hboPath);
            stampBefore = iModifiedMs(hboPath);
            pause(0.05);
            pm.executePipeline(runArgs{:});
            testCase.assertEqual(iModifiedMs(hboPath), stampBefore, ...
                'Precondition: an unchanged folder reuses the HemoCompute branch.');

            otherInfo = iConfiguredRigInfo('LED_593nm');
            otherInfo.uuid = 'other-run-test-rig';
            UMITRigLedger.rebind(folder, UMITRigLedger.makeBinding( ...
                UMITRigLedger.sourceRigFromRigInfo(otherInfo), ...
                UMITRigLedger.buildOpticsSnapshot(otherInfo), []), 'Origin', 'test');
            UMITRigLedger.acknowledgeConflict(folder, 'Origin', 'test');
            pause(0.05);

            pm.executePipeline(runArgs{:});
            testCase.verifyNotEqual(iModifiedMs(hboPath), stampBefore);
            testCase.verifyNotEqual(loadData(hboPath), hboBefore, ...
                'HbO must be recomputed with the acknowledged optics.');
            testCase.verifyTrue(UMITRigLedger.read(folder).hadConflict);

            stampAfter = iModifiedMs(hboPath);
            pause(0.05);
            pm.executePipeline(runArgs{:});
            testCase.verifyEqual(iModifiedMs(hboPath), stampAfter, ...
                'Outputs computed after the acknowledgement are reused again.');
        end

        function testUnboundDatasetIsRefusedNotDefaulted(testCase)
            %TESTUNBOUNDDATASETISREFUSEDNOTDEFAULTED No Active-rig fallback:
            %a folder without a rig ledger is an error, even when AcqInfos.mat
            %still names a rig and an Active rig exists, and nothing is written.

            acqPath = fullfile(testCase.TempFolder, 'AcqInfos.mat');
            delete(UMITRigLedger.filePath(testCase.TempFolder));
            loaded = load(acqPath, 'AcqInfoStream');
            AcqInfoStream = loaded.AcqInfoStream;
            AcqInfoStream.rigUUID = 'legacy-rig-uuid';
            AcqInfoStream.rigID = 'LegacyRig';
            saveMatAtomic(acqPath, 'AcqInfoStream', AcqInfoStream);
            acqBefore = dir(acqPath);

            info = HemoCompute('pipelineInfo');
            parameters = info.parameters;
            for iParameter = 1:numel(parameters)
                parameters(iParameter).value = parameters(iParameter).default;
            end
            testCase.verifyError(@() HemoCompute( ...
                'preflight', testCase.TempFolder, parameters), ...
                'Umitoolbox:RigLedger:noLedger');

            pm = PipelineManager({testCase.TempFolder});
            pm.addStep('HemoCompute', 'input', 'red.dat');
            [isValid, report] = pm.validateNodePreflight('throwOnError', false);
            testCase.verifyFalse(isValid);
            testCase.verifyTrue(any(contains(string(report.errors), 'no rig ledger')));

            acqAfter = dir(acqPath);
            testCase.verifyEqual(acqAfter.datenum, acqBefore.datenum);
            testCase.verifyFalse(UMITRigLedger.exists(testCase.TempFolder));
        end
    end
end

function folder = iPreparePMHemoFolder(rootFolder, rigStore)
folder = tempname(rootFolder);
mkdir(folder);
iCreateSyntheticHemoDataset(folder, true);
iBindDatasetToRig(folder, rigStore);
end

function rigInfo = iConfiguredRigInfo(redSpectrumID)
%ICONFIGUREDRIGINFO RigInfo-shaped rig definition over built-in spectra.
if nargin < 1
    redSpectrumID = 'LED_632nm';
end
cameras = struct( ...
    'index', 1, ...
    'displayName', 'Test Camera', ...
    'manufacturer', '', ...
    'model', 'D1024', ...
    'serialNumber', '', ...
    'spectrumID', 'PF1024');
illuminations = struct( ...
    'name', {'red','green','yellow'}, ...
    'displayName', {'Red','Green','Yellow'}, ...
    'manufacturer', {'','',''}, ...
    'model', {'','',''}, ...
    'spectrumID', {redSpectrumID,'LED_521nm','LED_593nm'});
rigInfo = struct('uuid', 'hemo-run-test-rig', 'rigID', 'HemoRunTestRig', ...
    'displayName', 'HemoCompute Test Rig', 'modifiedOn', datetime(2026, 10, 1), ...
    'cameras', cameras, 'illuminations', illuminations);
end

function stamp = iModifiedMs(filePath)
%IMODIFIEDMS File modification time in milliseconds (finer than dir's datenum).
file = java.io.File(filePath);
stamp = file.lastModified();
end

function iBindDatasetToRig(folderPath, rigInfo)
%IBINDDATASETTORIG Bind FOLDERPATH's rig ledger to a RigInfo struct.
UMITRigLedger.bind(folderPath, UMITRigLedger.makeBinding( ...
    UMITRigLedger.sourceRigFromRigInfo(rigInfo), ...
    UMITRigLedger.buildOpticsSnapshot(rigInfo), []), 'Origin', 'TestRunHemoCompute');
end

function iCreateSyntheticHemoDataset(folderPath, repeatedRed)
%ICREATESYNTHETICHEMODATASET Create normalized intrinsic-channel .dat files.

Ny = 8;
Nx = 8;
baseLength = 12;
baseFreq = 10;

if repeatedRed
    redLength = baseLength * 2;
    redFreq = baseFreq * 2;
else
    redLength = baseLength;
    redFreq = baseFreq;
end

AcqInfoStream = struct();
AcqInfoStream.Height = Ny;
AcqInfoStream.Width = Nx;
AcqInfoStream.Length = baseLength;
AcqInfoStream.FrameRateHz = baseFreq;
AcqInfoStream.Camera_Model = 'D1024';
AcqInfoStream.MultiCam = false;
AcqInfoStream.Datatype = 'single';
AcqInfoStream.ExposureMsec = 10;
AcqInfoStream.ImportedChannels = struct( ...
    'DatFile', {}, ...
    'Length', {}, ...
    'FrameRateHz', {}, ...
    'ExposureMsec', {}, ...
    'CamIdx', {});

channels = { ...
    'red.dat', redLength, redFreq, 0.010; ...
    'green.dat', baseLength, baseFreq, 0.015; ...
    'yellow.dat', baseLength, baseFreq, 0.020};

for iChan = 1:size(channels, 1)
    datFile = channels{iChan, 1};
    nFrames = channels{iChan, 2};
    frameRate = channels{iChan, 3};
    amp = channels{iChan, 4};

    data = iMakeNormalizedMovie(Ny, Nx, nFrames, amp);
    iWriteSingleDat(fullfile(folderPath, datFile), data, frameRate, 10);

    AcqInfoStream.ImportedChannels(end+1).DatFile = datFile;
    AcqInfoStream.ImportedChannels(end).Length = nFrames;
    AcqInfoStream.ImportedChannels(end).FrameRateHz = frameRate;
    AcqInfoStream.ImportedChannels(end).ExposureMsec = 10;
    AcqInfoStream.ImportedChannels(end).CamIdx = 1;
end

save(fullfile(folderPath, 'AcqInfos.mat'), 'AcqInfoStream');
end

function data = iMakeNormalizedMovie(Ny, Nx, Nt, amp)
baseMap = single(reshape(linspace(0, 1, Ny * Nx), Ny, Nx));
t = reshape(single(1:Nt), 1, 1, []);
data = single(1 + amp * sin(2*pi*t/max(Nt, 1)) + 0.001 * baseMap);
end

function iWriteSingleDat(filePath, data, frameRateHz, exposureMsec)
% Headered input (.dat header Phase 5a).
writeTestDat(filePath, single(data), frameRateHz, exposureMsec);
end

function iMakeAcqInfosBound(folder)
% Rewrite every .dat as a headerless file without sidecar (described only
% by AcqInfos.mat), which every reader rejects since .dat header Phase 5b.
listing = dir(fullfile(folder, '*.dat'));
for k = 1:numel(listing)
    f = fullfile(folder, listing(k).name);
    values = loadData(f);
    fid = fopen(f, 'w');
    fwrite(fid, values, class(values));
    fclose(fid);
end
end
