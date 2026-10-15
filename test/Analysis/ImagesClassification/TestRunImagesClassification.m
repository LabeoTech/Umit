classdef TestRunImagesClassification < matlab.unittest.TestCase
    %TESTRUNIMAGESCLASSIFICATION Unit tests for ImagesClassification.
    %
    %   This suite verifies wrapper behavior, backup handling, returned
    %   SaveFolder-relative manifests, imported-channel metadata, and
    %   repeated-illumination grouping through a synthetic info.txt mutation.

    properties
        ProjectRoot
        FixtureFolder
        TempSaveFolder
        TempRawFolder
    end

    methods (TestMethodSetup)
        function setup(testCase)
            thisFile = mfilename('fullpath');
            testFolder = fileparts(thisFile);
            projectRoot = extractBefore(testFolder, [filesep 'test']);
            if isempty(projectRoot)
                projectRoot = fileparts(fileparts(testFolder));
            end
            testCase.ProjectRoot = char(projectRoot);
            addpath(genpath(testCase.ProjectRoot));

            testCase.FixtureFolder = fullfile(testCase.ProjectRoot, 'test', 'Analysis', 'TestingData_with_events');
            testCase.TempSaveFolder = fullfile(tempdir, ['TestRunImagesClassification_' char(java.util.UUID.randomUUID)]);
            testCase.TempRawFolder = fullfile(tempdir, ['TestRunImagesClassificationRaw_' char(java.util.UUID.randomUUID)]);
            mkdir(testCase.TempSaveFolder);
            % The fixture is single-camera data: ImagesClassification
            % resolves no rig for it, so no rig-library fixture is needed
            % (rig-ledger redesign, Phase 3).
        end
    end

    methods (TestMethodTeardown)
        function teardown(testCase)
            if isfolder(testCase.TempSaveFolder)
                rmdir(testCase.TempSaveFolder, 's');
            end
            if isfolder(testCase.TempRawFolder)
                rmdir(testCase.TempRawFolder, 's');
            end
        end
    end

    methods (Test)
        function testPipelineManagerExecutesImporter(testCase)
            pm = buildPMForScenario(testCase.TempSaveFolder, ...
                'ImagesClassification', 'auto', ...
                'RawFolder', testCase.FixtureFolder);
            pm.executePipeline('PrintSummary', false);

            testCase.verifyTrue(isfile(fullfile(testCase.TempSaveFolder, ...
                'AcqInfos.mat')));
            files = dir(fullfile(testCase.TempSaveFolder, '*.dat'));
            testCase.verifyNotEmpty(files);
        end

        function testPipelineInfo(testCase)
            %TESTPIPELINEINFO Validate pipelineInfo query.

            info = ImagesClassification('pipelineInfo');
            testCase.verifyTrue(isstruct(info) && isscalar(info));
            testCase.verifyEqual(info.name, 'ImagesClassification');
            testCase.verifyEqual(info.freshSaveFolderRole, ...
                'acquisition-initializer');
            parameterNames = string({info.parameters.name});
            testCase.verifyFalse(any(parameterNames == "ApplyCoregistration"), ...
                'Coregistration is resolved automatically from the Rig backend.');

            % A pipeline-driven run cannot answer an interactive prompt, so it
            % takes this default. It must not be the destructive 'ERASE', which
            % deletes managed files from SaveFolder with no backup (P1-6).
            backupParam = info.parameters(strcmp({info.parameters.name}, 'backupOpts'));
            testCase.verifyEqual(backupParam.default, 'GENBACKUP');
            testCase.verifyEqual(sort(backupParam.allowed), sort({'ERASE','GENBACKUP'}));
        end

        function testReturnsRelativeFileManifest(testCase)
            %TESTRETURNSRELATIVEFILEMANIFEST Return SaveFolder-relative output names.

            [outFile, rigResolution] = ImagesClassification( ...
                testCase.FixtureFolder, ...
                testCase.TempSaveFolder, ...
                'BinningSpatial', 1, ...
                'BinningTemp', 1, ...
                'backupOpts', 'ERASE');

            testCase.verifyTrue(iscell(outFile));
            testCase.verifyTrue(isfile(fullfile(testCase.TempSaveFolder, 'AcqInfos.mat')));
            testCase.verifyEqual(rigResolution.ledgerState, 'NoRig');
            testCase.verifyEmpty(rigResolution.rigUUID);
            testCase.verifyEqual(rigResolution.cameraCoregistration.skipReason, ...
                'singleCameraAcquisition');
            info = load(fullfile(testCase.TempSaveFolder, 'AcqInfos.mat'), 'AcqInfoStream');
            testCase.verifyFalse(any(isfield(info.AcqInfoStream, {'rigUUID', 'rigID'})));
            testCase.verifyEqual(UMITRigLedger.read(testCase.TempSaveFolder).state, 'NoRig');
            testCase.verifyTrue(all(endsWith(outFile, '.dat')), ...
                'ImagesClassification should return imported .dat outputs only.');
            for iFile = 1:numel(outFile)
                testCase.verifyTrue(isfile(fullfile(testCase.TempSaveFolder, outFile{iFile})));
            end
        end

        function testRunsWithoutTformFile(testCase)
            %TESTRUNSWITHOUTTFORMFILE Run without requiring alignment transform files.

            outFile = ImagesClassification( ...
                testCase.FixtureFolder, ...
                testCase.TempSaveFolder, ...
                'BinningSpatial', 1, ...
                'BinningTemp', 1, ...
                'backupOpts', 'ERASE');

            testCase.verifyNotEmpty(outFile);
        end

        function testBackupEraseRemovesManagedFiles(testCase)
            %TESTBACKUPERASEREMOVESMANAGEDFILES Remove old outputs before run.

            dummyFile = fullfile(testCase.TempSaveFolder, 'old_result.dat');
            writeTestDat(dummyFile, single(1), 1);

            testCase.verifyTrue(isfile(dummyFile));

            outFile = ImagesClassification( ...
                testCase.FixtureFolder, ...
                testCase.TempSaveFolder, ...
                'BinningSpatial', 1, ...
                'BinningTemp', 1, ...
                'backupOpts', 'ERASE');

            testCase.verifyFalse(isfile(dummyFile));
            testCase.verifyNotEmpty(outFile);
            testCase.verifyTrue(isfile(fullfile(testCase.TempSaveFolder, 'AcqInfos.mat')));
            testCase.verifyTrue(all(endsWith(outFile, '.dat')));
        end

        function testBackupGenBackupCreatesZip(testCase)
            %TESTBACKUPGENBACKUPCREATESZIP Create backup archive before run.

            dummyFile = fullfile(testCase.TempSaveFolder, 'old_result.dat');
            writeTestDat(dummyFile, single(1), 1);

            outFile = ImagesClassification( ...
                testCase.FixtureFolder, ...
                testCase.TempSaveFolder, ...
                'BinningSpatial', 1, ...
                'BinningTemp', 1, ...
                'backupOpts', 'GENBACKUP');

            zipList = dir(fullfile(testCase.TempSaveFolder, 'bkp_*.zip'));
            testCase.verifyNotEmpty(zipList);
            testCase.verifyFalse(isfile(dummyFile));
            testCase.verifyNotEmpty(outFile);
        end

        function testWrapperCreatesImportedChannels(testCase)
            %TESTWRAPPERCREATESIMPORTEDCHANNELS Validate AcqInfos output organization.

            outFile = ImagesClassification( ...
                testCase.FixtureFolder, ...
                testCase.TempSaveFolder, ...
                'BinningSpatial', 1, ...
                'BinningTemp', 1, ...
                'backupOpts', 'ERASE');

            acqPath = fullfile(testCase.TempSaveFolder, 'AcqInfos.mat');
            testCase.verifyTrue(isfile(acqPath));
            testCase.verifyTrue(all(endsWith(outFile, '.dat')));

            S = load(acqPath, 'AcqInfoStream');
            AcqInfoStream = S.AcqInfoStream;

            testCase.verifyTrue(isfield(AcqInfoStream, 'ImportedChannels'));
            importedChannels = AcqInfoStream.ImportedChannels;
            testCase.verifyNotEmpty(importedChannels);
            testCase.verifyFalse(isfield(importedChannels, 'RepeatCount'));
            testCase.verifyFalse(isfield(importedChannels, 'SequenceIdx'));

            for iChan = 1:numel(importedChannels)
                thisPath = fullfile(testCase.TempSaveFolder, importedChannels(iChan).DatFile);
                testCase.verifyTrue(isfile(thisPath));
                testCase.verifyEqual(double(importedChannels(iChan).CamIdx), 1);
            end
        end

        function testWrapperRepeatingIlluminationManifestAndMetadata(testCase)
            %TESTWRAPPERREPEATINGILLUMINATIONMANIFESTANDMETADATA Validate grouped output.

            [rawFolder, repeatedDatFile] = testCase.makeRepeatedIlluminationFixture();

            outFile = ImagesClassification( ...
                rawFolder, ...
                testCase.TempSaveFolder, ...
                'BinningSpatial', 1, ...
                'BinningTemp', 1, ...
                'backupOpts', 'ERASE');

            for iFile = 1:numel(outFile)
                testCase.verifyTrue(isfile(fullfile(testCase.TempSaveFolder, outFile{iFile})));
            end

            S = load(fullfile(testCase.TempSaveFolder, 'AcqInfos.mat'), 'AcqInfoStream');
            AcqInfoStream = S.AcqInfoStream;
            importedChannels = AcqInfoStream.ImportedChannels;

            testCase.verifyTrue(any(strcmp({importedChannels.DatFile}, repeatedDatFile)), ...
                ['Repeated channel file not found: ' repeatedDatFile]);

            idxRep = find(strcmp({importedChannels.DatFile}, repeatedDatFile), 1, 'first');
            % .dat header Phase 7b: AcqInfos.mat keeps the raw info.txt
            % values and no Length; the base timeline is a non-repeated
            % channel's header.
            rawInfo = ReadInfoFile(rawFolder);
            testCase.verifyEqual(double(AcqInfoStream.FrameRateHz), double(rawInfo.FrameRateHz));
            testCase.verifyFalse(isfield(AcqInfoStream, 'Length'));
            testCase.verifyFalse(isfield(AcqInfoStream, 'Datatype'));
            idxBase = find(~strcmp({importedChannels.DatFile}, repeatedDatFile), 1, 'first');
            baseInfo = loadMetaData(fullfile(testCase.TempSaveFolder, importedChannels(idxBase).DatFile));
            baseLength = datAxisSize(baseInfo, 'T');
            baseFreq = double(baseInfo.frameRateHz);

            testCase.verifyEqual(double(importedChannels(idxRep).Length), 2 * baseLength);
            testCase.verifyEqual(double(importedChannels(idxRep).FrameRateHz), 2 * baseFreq, ...
                'AbsTol', max(1e-9, abs(2 * baseFreq) * 1e-9));

            Info = loadMetaData(fullfile(testCase.TempSaveFolder, repeatedDatFile));
            testCase.verifyEqual(datAxisSize(Info, 'T'), 2 * baseLength);
            testCase.verifyEqual(double(Info.frameRateHz), 2 * baseFreq, ...
                'AbsTol', max(1e-9, abs(2 * baseFreq) * 1e-9));
        end

        function testApplyRigFalseClassifiesOnly(testCase)
            %TESTAPPLYRIGFALSECLASSIFIESONLY Skip Rig binding and DataParams.

            [outFile, rigResolution] = ImagesClassification( ...
                testCase.FixtureFolder, ...
                testCase.TempSaveFolder, ...
                'backupOpts', 'ERASE', ...
                'ApplyRig', false);

            testCase.verifyEmpty(rigResolution);
            testCase.verifyNotEmpty(outFile);
            testCase.verifyTrue(all(endsWith(outFile, '.dat')));
            for iFile = 1:numel(outFile)
                testCase.verifyTrue(isfile(fullfile(testCase.TempSaveFolder, outFile{iFile})));
            end
            S = load(fullfile(testCase.TempSaveFolder, 'AcqInfos.mat'), 'AcqInfoStream');
            testCase.verifyFalse(isfield(S.AcqInfoStream, 'rigUUID'));
            testCase.verifyFalse(isfield(S.AcqInfoStream, 'rigID'));
            testCase.verifyFalse(isfile(fullfile(testCase.TempSaveFolder, 'DataParams.mat')));
            testCase.verifyFalse(UMITRigLedger.exists(testCase.TempSaveFolder));
        end

        function testSingleCameraImportLeavesRigLibraryUntouched(testCase)
            %TESTSINGLECAMERAIMPORTLEAVESRIGLIBRARYUNTOUCHED No rig is
            %resolved, created or activated for single-camera data.

            rigsBefore = UMITRigStore.listRigs();
            ImagesClassification(testCase.FixtureFolder, testCase.TempSaveFolder, ...
                'backupOpts', 'ERASE');
            rigsAfter = UMITRigStore.listRigs();

            testCase.verifyEqual(rigsAfter(:, {'RigUUID', 'Status', 'IsDefault'}), ...
                rigsBefore(:, {'RigUUID', 'Status', 'IsDefault'}));
        end

        function testSingleCameraReimportKeepsBinding(testCase)
            %TESTSINGLECAMERAREIMPORTKEEPSBINDING A re-import into a bound
            %folder writes only a reimport event (decision Q3).

            binding = testCase.makeLedgerBinding();
            UMITRigLedger.bind(testCase.TempSaveFolder, binding, 'Origin', 'test');
            UMITRigLedger.commit(testCase.TempSaveFolder, 'optics', ...
                binding.optics.hash, 'Origin', 'test');

            [~, rigResolution] = ImagesClassification(testCase.FixtureFolder, ...
                testCase.TempSaveFolder, 'backupOpts', 'ERASE');

            ledger = UMITRigLedger.read(testCase.TempSaveFolder);
            types = cellfun(@(e) e.type, ledger.events, 'UniformOutput', false);
            testCase.verifyEqual(types, {'bind', 'commit', 'reimport'});
            testCase.verifyEqual(ledger.state, 'BoundCommitted');
            testCase.verifyEqual(ledger.binding.opticsHash, binding.optics.hash);
            testCase.verifyEqual(rigResolution.ledgerState, 'BoundCommitted');
            testCase.verifyEqual(rigResolution.rigID, 'ClassificationTestRig');
        end

        function testRepeatedSingleCameraImportRecordsReimport(testCase)
            ImagesClassification(testCase.FixtureFolder, testCase.TempSaveFolder, ...
                'backupOpts', 'ERASE');
            ImagesClassification(testCase.FixtureFolder, testCase.TempSaveFolder, ...
                'backupOpts', 'ERASE');

            ledger = UMITRigLedger.read(testCase.TempSaveFolder);
            types = cellfun(@(e) e.type, ledger.events, 'UniformOutput', false);
            testCase.verifyEqual(types, {'noRig', 'reimport'});
            testCase.verifyEqual(ledger.state, 'NoRig');
        end

        function testInvalidLedgerBlocksImportBeforeClassifying(testCase)
            fid = fopen(UMITRigLedger.filePath(testCase.TempSaveFolder), 'w');
            fprintf(fid, '{not json');
            fclose(fid);

            testCase.verifyError(@() ImagesClassification(testCase.FixtureFolder, ...
                testCase.TempSaveFolder, 'backupOpts', 'ERASE'), ...
                'Umitoolbox:RigLedger:invalidLedger');
            testCase.verifyEmpty(dir(fullfile(testCase.TempSaveFolder, '*.dat')));
            testCase.verifyFalse(isfile(fullfile(testCase.TempSaveFolder, 'AcqInfos.mat')));
        end

        function testMissingCameraModelWarnsAndAssumesD1024(testCase)
            %TESTMISSINGCAMERAMODELWARNSANDASSUMESD1024 Warn on the D1024 fallback.

            rawFolder = testCase.TempRawFolder;
            mkdir(rawFolder);
            copyfile(fullfile(testCase.FixtureFolder, '*'), rawFolder);
            infoPath = fullfile(rawFolder, 'info.txt');
            txt = fileread(infoPath);
            testCase.assumeNotEmpty(regexp(txt, '(?m)^Camera Model:', 'once'));
            txt = regexprep(txt, '(?m)^Camera Model:[^\r\n]*(\r?\n)', '');
            fid = fopen(infoPath, 'w');
            testCase.assertNotEqual(fid, -1, 'Failed to open synthetic info.txt for writing.');
            fwrite(fid, txt, 'char');
            fclose(fid);

            testCase.verifyWarning(@() ImagesClassification(rawFolder, ...
                testCase.TempSaveFolder, 'backupOpts', 'ERASE', 'ApplyRig', false), ...
                'Umitoolbox:ImagesClassification:CameraModelMissing');

            S = load(fullfile(testCase.TempSaveFolder, 'AcqInfos.mat'), 'AcqInfoStream');
            testCase.verifyEqual(S.AcqInfoStream.Camera_Model, 'D1024');
        end
    end

    methods (Access = private)
        function [rawFolder, repeatedDatFile] = makeRepeatedIlluminationFixture(testCase)
            %MAKEREPEATEDILLUMINATIONFIXTURE Copy fixture and repeat one illumination.

            rawFolder = testCase.TempRawFolder;
            if isfolder(rawFolder)
                rmdir(rawFolder, 's');
            end
            mkdir(rawFolder);
            copyfile(fullfile(testCase.FixtureFolder, '*'), rawFolder);

            infoPath = fullfile(rawFolder, 'info.txt');
            testCase.assertTrue(isfile(infoPath), 'Fixture folder must contain info.txt.');

            txt = fileread(infoPath);
            expr = '(?m)^Illumination(\d+):\s*([^\r\n]*)';
            tok = regexp(txt, expr, 'tokens');
            testCase.assumeGreaterThanOrEqual(numel(tok), 3, ...
                'Repeating-illumination test requires at least three illumination entries.');

            illumIdx = cellfun(@(x) str2double(x{1}), tok);
            illumValues = cellfun(@(x) strtrim(x{2}), tok, 'UniformOutput', false);
            datTags = cellfun(@(x) testCase.getDatTagFromColor(x), illumValues, 'UniformOutput', false);

            candidateIdx = find(cellfun(@(x) sum(strcmp(datTags, x)) == 1, datTags), 1, 'first');
            testCase.assumeNotEmpty(candidateIdx, ...
                'Repeating-illumination test requires at least one uniquely tagged illumination entry.');

            targetIdx = numel(illumIdx);
            if targetIdx == candidateIdx
                targetIdx = numel(illumIdx) - 1;
            end

            repeatedColor = illumValues{candidateIdx};
            repeatedDatFile = datTags{candidateIdx};

            pattern = sprintf('(?m)^(Illumination%d:\\s*)[^\\r\\n]*', illumIdx(targetIdx));
            txt = regexprep(txt, pattern, ['$1' repeatedColor], 'once');

            fid = fopen(infoPath, 'w');
            testCase.assertNotEqual(fid, -1, 'Failed to open synthetic info.txt for writing.');
            cleaner = onCleanup(@() fclose(fid));
            fwrite(fid, txt, 'char');
        end
    end

    methods (Static, Access = private)
        function binding = makeLedgerBinding()
            %MAKELEDGERBINDING Ledger binding from a RigInfo struct (no library rig).

            rigInfo = struct('uuid', 'classification-test-rig', ...
                'rigID', 'ClassificationTestRig', 'displayName', 'Classification Test Rig', ...
                'modifiedOn', datetime(2026, 10, 1), ...
                'cameras', struct('index', 1, 'spectrumID', 'PF1024', 'model', 'D1024'), ...
                'illuminations', struct('name', {'red', 'green'}, ...
                    'spectrumID', {'LED_632nm', 'LED_521nm'}));
            binding = UMITRigLedger.makeBinding( ...
                UMITRigLedger.sourceRigFromRigInfo(rigInfo), ...
                UMITRigLedger.buildOpticsSnapshot(rigInfo), []);
        end

        function dTag = getDatTagFromColor(colorName)
            %GETDATTAGFROMCOLOR Match ImagesClassification output naming.

            if contains(colorName, {'red', 'green'}, 'IgnoreCase', true)
                dTag = [lower(char(string(colorName))) '.dat'];
            elseif contains(colorName, 'amber', 'IgnoreCase', true)
                dTag = 'yellow.dat';
            elseif contains(colorName, 'fluo', 'IgnoreCase', true)
                waveTag = regexp(colorName, '[0-9]{3}', 'match');
                if ~isempty(waveTag)
                    dTag = ['fluo_' waveTag{1} '.dat'];
                else
                    dTag = 'fluo.dat';
                end
            else
                dTag = 'speckle.dat';
            end
        end
    end
end
