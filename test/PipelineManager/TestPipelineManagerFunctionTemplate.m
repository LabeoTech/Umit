classdef TestPipelineManagerFunctionTemplate < matlab.unittest.TestCase
    %TESTPIPELINEMANAGERFUNCTIONTEMPLATE Validate the analysis templates.

    properties
        ProjectRoot char
        SaveFolder char
        InputData single
    end

    methods (TestClassSetup)
        function resolveProjectRoot(testCase)
            thisFile = mfilename('fullpath');
            testFolder = fileparts(thisFile);
            testCase.ProjectRoot = extractBefore(testFolder, [filesep 'test']);
        end
    end

    methods (TestMethodSetup)
        function createSaveFolder(testCase)
            fixture = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.SaveFolder = fullfile(fixture.Folder, 'SaveFolder');
            mkdir(testCase.SaveFolder);

            testCase.InputData = single(reshape(-12:11, [2 3 4]));

            AcqInfoStream = struct();
            AcqInfoStream.Width = 3;
            AcqInfoStream.Height = 2;
            AcqInfoStream.Length = 4;
            AcqInfoStream.FrameRateHz = 20;
            AcqInfoStream.ExposureMsec = 1;
            AcqInfoStream = appendImportedChannelInfo(AcqInfoStream, struct( ...
                'DatFile', 'input.dat', 'Length', 4, 'FrameRateHz', 20));
            save(fullfile(testCase.SaveFolder, 'AcqInfos.mat'), ...
                'AcqInfoStream');

            saveData(fullfile(testCase.SaveFolder, 'input.dat'), ...
                testCase.InputData, 'DimNames', {'Y', 'X', 'T'}, 'FrameRateHz', 20);
        end
    end

    methods (Test)
        function testOfficialTemplatesRemainHidden(testCase)
            templateNames = testCase.templateNames();
            analysisFolder = fullfile(testCase.ProjectRoot, 'Analysis');

            for iTemplate = 1:numel(templateNames)
                testCase.verifyTrue(isfile(fullfile(analysisFolder, ...
                    [templateNames{iTemplate} '.m'])));
            end

            pm = PipelineManager(testCase.SaveFolder, '', testCase.ProjectRoot);
            discoveredNames = {pm.funcList.name};
            for iTemplate = 1:numel(templateNames)
                testCase.verifyFalse(any(strcmpi(discoveredNames, ...
                    templateNames{iTemplate})));
            end
        end

        function testTemporaryCopiesAreDiscoveredAndParsed(testCase)
            parserFolder = createTemporaryTemplateCategory(testCase.ProjectRoot);
            cleanup = onCleanup(@() cleanupTemplateCategory( ...
                parserFolder, testCase.ProjectRoot));

            pm = PipelineManager(testCase.SaveFolder, '', testCase.ProjectRoot);
            templateNames = testCase.templateNames();

            for iTemplate = 1:numel(templateNames)
                name = templateNames{iTemplate};
                matches = find(strcmpi({pm.funcList.name}, name));
                testCase.verifyNumElements(matches, 1);

                info = pm.getPipelineInfo(name);
                testCase.verifyEqual(info.name, name);
            end
        end

        function testPipelineInfoUsesSaveFolderContract(testCase)
            templateNames = testCase.templateNames();

            for iTemplate = 1:numel(templateNames)
                info = feval(templateNames{iTemplate}, 'pipelineInfo');
                testCase.verifyTrue(all(isfield(info, ...
                    {'inputs','outputs','parameters','arguments'})));

                inputNames = {info.inputs.name};
                outputNames = {info.outputs.name};
                testCase.verifyFalse(any(strcmpi(inputNames, 'metaData')));
                testCase.verifyFalse(any(strcmpi(outputNames, 'metaData')));

                saveFolderInput = info.inputs( ...
                    strcmpi(inputNames, 'SaveFolder'));
                testCase.verifyNumElements(saveFolderInput, 1);
                testCase.verifyEqual(saveFolderInput.type, {'SaveFolder'});
                testCase.verifyFalse(saveFolderInput.isData);

                saveFolderArgument = info.arguments( ...
                    strcmpi({info.arguments.name}, 'SaveFolder'));
                testCase.verifyNumElements(saveFolderArgument, 1);
                testCase.verifyEqual(saveFolderArgument.kind, 'input');
                testCase.verifyEqual(saveFolderArgument.callType, 'positional');
                testCase.verifyFalse(saveFolderArgument.isData);
            end

            numericInfo = funcTemplate('pipelineInfo');
            numericInput = numericInfo.inputs( ...
                strcmp({numericInfo.inputs.name}, 'data'));
            testCase.verifyFalse(numericInput.supportsFile);
            testCase.verifyEqual(numericInput.dataMode, 'ram');

            scale = numericInfo.parameters( ...
                strcmp({numericInfo.parameters.name}, 'Scale'));
            testCase.verifyEqual(scale.default, 1);
            testCase.verifyEqual(scale.allowed, [0 Inf]);

            numericOutput = numericInfo.outputs(1);
            testCase.verifyEqual(numericOutput.outputMode, 'data');
            testCase.verifyEqual(numericOutput.type, {'ImageTimeSeries'});
            testCase.verifyEqual(numericOutput.defOutfilename, ...
                'funcTemplate.dat');
            testCase.verifyEqual(numericOutput.saveFileName, ...
                'funcTemplate.dat');
        end

        function testFreshAcquisitionTemplateRoles(testCase)
            initializerInfo = ...
                funcTemplateAcquisitionInitializer('pipelineInfo');
            companionInfo = ...
                funcTemplateAcquisitionCompanion('pipelineInfo');

            testCase.verifyEqual(initializerInfo.freshSaveFolderRole, ...
                'acquisition-initializer');
            testCase.verifyEqual(companionInfo.freshSaveFolderRole, ...
                'acquisition-companion');

            ordinaryTemplates = { ...
                'funcTemplate', ...
                'funcTemplateUMT', ...
                'funcTemplateFileManifest'};
            for iTemplate = 1:numel(ordinaryTemplates)
                info = feval(ordinaryTemplates{iTemplate}, 'pipelineInfo');
                testCase.verifyEqual(info.freshSaveFolderRole, 'none');
            end
        end

        function testFreshAcquisitionTemplatesPreserveOwnership(testCase)
            fixtureRoot = fileparts(testCase.SaveFolder);
            rawFolder = fullfile(fixtureRoot, 'TemplateRaw');
            freshFolder = fullfile(fixtureRoot, 'TemplateFresh');
            mkdir(rawFolder);
            mkdir(freshFolder);

            imageData = testCase.InputData;
            AcqInfoStream = struct( ...
                'Width', size(imageData, 2), ...
                'Height', size(imageData, 1), ...
                'Length', size(imageData, 3), ...
                'FrameRateHz', 20, ...
                'ExposureMsec', 1, ...
                'Datatype', 'single');
            save(fullfile(rawFolder, 'funcTemplateAcquisition.mat'), ...
                'imageData', 'AcqInfoStream', '-mat');

            initializerFiles = funcTemplateAcquisitionInitializer( ...
                rawFolder, freshFolder);
            testCase.verifyEqual(initializerFiles, ...
                {'funcTemplateImported.dat'});
            testCase.verifyTrue(isfile(fullfile(freshFolder, ...
                'funcTemplateImported.dat')));
            % .dat header Phase 4d: the imported file carries its own rate
            % and exposure.
            importedHeader = readDatHeader(fullfile(freshFolder, 'funcTemplateImported.dat'));
            testCase.verifyEqual(importedHeader.frameRateHz, 20);
            testCase.verifyEqual(importedHeader.exposureMsec, 1);
            testCase.verifyEqual(importedHeader.channelName, 'funcTemplateImported');
            testCase.verifyTrue(isfile(fullfile(freshFolder, ...
                'AcqInfos.mat')));
            [isLegacy, legacyMessage] = isLegacySchemaFolder(freshFolder);
            testCase.verifyFalse(isLegacy, legacyMessage);
            % .dat header Phase 7b: AcqInfos.mat holds the raw acquisition,
            % not the imported Length or Datatype.
            saved = load(fullfile(freshFolder, 'AcqInfos.mat'), 'AcqInfoStream');
            testCase.verifyFalse(isfield(saved.AcqInfoStream, 'Length'));
            testCase.verifyFalse(isfield(saved.AcqInfoStream, 'Datatype'));
            testCase.verifyEqual(saved.AcqInfoStream.FrameRateHz, 20);
            testCase.verifyEqual(saved.AcqInfoStream.ImportedChannels.Length, ...
                size(imageData, 3));

            beforeCompanion = load(fullfile(freshFolder, ...
                'AcqInfos.mat'), 'AcqInfoStream');
            companionFiles = funcTemplateAcquisitionCompanion( ...
                rawFolder, freshFolder);
            afterCompanion = load(fullfile(freshFolder, ...
                'AcqInfos.mat'), 'AcqInfoStream');

            testCase.verifyEqual(companionFiles, ...
                {'funcTemplateCompanion.mat'});
            testCase.verifyTrue(isfile(fullfile(freshFolder, ...
                'funcTemplateCompanion.mat')));
            testCase.verifyEqual(afterCompanion.AcqInfoStream, ...
                beforeCompanion.AcqInfoStream);
        end

        function testNumericTemplateDefaultAndOptions(testCase)
            defaultOut = funcTemplate(testCase.InputData, testCase.SaveFolder);
            testCase.verifyEqual(defaultOut, testCase.InputData);
            testCase.verifyClass(defaultOut, 'single');
            testCase.verifySize(defaultOut, size(testCase.InputData));

            customOut = funcTemplate(testCase.InputData, testCase.SaveFolder, ...
                'Scale', 2, 'ClipNegative', true);
            expected = testCase.InputData .* single(2);
            expected(expected < 0) = 0;

            testCase.verifyEqual(customOut, expected);
            testCase.verifyClass(customOut, 'single');
            testCase.verifySize(customOut, size(testCase.InputData));
        end

        function testUMTTemplateCreatesValidatedOutput(testCase)
            outData = funcTemplateUMT( ...
                testCase.InputData, testCase.SaveFolder, ...
                'EntryName', 'Average');

            validateUMTStruct(outData, 'requireEventInfo', false);
            testCase.verifyEqual(outData.kind, 'image');
            testCase.verifyTrue(isfield(outData.data, 'Average'));

            entry = outData.data.Average;
            testCase.verifyEqual(entry.dimNames, {'Y','X'});
            testCase.verifyEqual(entry.value, ...
                mean(testCase.InputData, 3, 'omitnan'));
            % .dat header Phase 7a: genUMTStruct no longer fills the folder
            % (AcqInfos.mat) rate; this Y-X map has no T axis and no rate.
            testCase.verifyFalse(isfield(entry, 'meta') && isfield(entry.meta, 'FrameRateHz'));

            info = funcTemplateUMT('pipelineInfo');
            output = info.outputs(1);
            testCase.verifyEqual(output.type, {'ProcessedData'});
            testCase.verifyEqual(output.outputMode, 'data');
            testCase.verifyEqual(output.defOutfilename, ...
                'funcTemplateMean.umt');
            testCase.verifyEqual(output.saveFileName, ...
                'funcTemplateMean.umt');
        end

        function testFileManifestMatchesCreatedFiles(testCase)
            outFiles = funcTemplateFileManifest( ...
                testCase.InputData, testCase.SaveFolder);
            expectedFiles = { ...
                'funcTemplateSummary.mat', ...
                'funcTemplateSummary.txt'};
            testCase.verifyEqual(outFiles, expectedFiles);

            for iFile = 1:numel(outFiles)
                [folderPart, ~, ~] = fileparts(outFiles{iFile});
                testCase.verifyEmpty(folderPart);
                testCase.verifyTrue(isfile(fullfile( ...
                    testCase.SaveFolder, outFiles{iFile})));
            end

            loaded = load(fullfile(testCase.SaveFolder, outFiles{1}), ...
                'summary');
            testCase.verifyEqual(loaded.summary.Size, ...
                size(testCase.InputData));
            testCase.verifyEqual(loaded.summary.MatlabClass, 'single');
            testCase.verifyEqual(loaded.summary.Mean, ...
                mean(double(testCase.InputData(:)), 'omitnan'));

            reportText = fileread(fullfile(testCase.SaveFolder, outFiles{2}));
            testCase.verifySubstring(reportText, 'Size: [2 3 4]');
            testCase.verifySubstring(reportText, 'MATLAB class: single');

            info = funcTemplateFileManifest('pipelineInfo');
            output = info.outputs(1);
            testCase.verifyEqual(output.outputMode, 'file');
            testCase.verifyEqual(output.defOutfilename, expectedFiles);
            testCase.verifyFalse(output.isData);
            testCase.verifyEmpty(output.saveFileName);
            testCase.verifyTrue(output.returnsValue);
        end

        function testNonReturningFileEffectMetadataValidation(testCase)
            info = PipelineManager.createPipelineInfo('sideEffectFixture', ...
                'Validate non-returning file-effect metadata.');
            info = PipelineManager.addOutput(info, 'affectedFiles', ...
                'UnknownDataType', 'file', 'Files changed by side effect.', ...
                '*.dat', 1, 'isData', false, 'returnsValue', false);

            testCase.verifyFalse(info.outputs(1).isData);
            testCase.verifyFalse(info.outputs(1).returnsValue);
            testCase.verifyTrue(info.outputs(1).isRequired);
            testCase.verifyEqual(info.outputs(1).defOutfilename, '*.dat');

            info = PipelineManager.addOutput(info, 'optionalReturnedFile', ...
                'UnknownDataType', 'file', 'Optional returned file.', ...
                'optional.txt', 2, 'isData', false, 'isRequired', false);
            testCase.verifyFalse(info.outputs(2).isRequired);
            testCase.verifyTrue(info.outputs(2).returnsValue);

            testCase.verifyError(@() PipelineManager.addOutput(info, ...
                'invalidDataEffect', 'ImageTimeSeries', 'file', '', ...
                '*.dat', 2, 'isData', true, 'returnsValue', false), ...
                'addOutput:InvalidNonReturningOutput');
            testCase.verifyError(@() PipelineManager.addOutput(info, ...
                'invalidRamEffect', 'UnknownDataType', 'data', '', ...
                '', 2, 'isData', false, 'returnsValue', false), ...
                'addOutput:InvalidNonReturningOutput');
            testCase.verifyError(@() PipelineManager.addOutput(info, ...
                'invalidOptionalData', 'ImageTimeSeries', 'data', '', ...
                'optional.dat', 3, 'isData', true, 'isRequired', false), ...
                'addOutput:InvalidOptionalOutput');
            testCase.verifyError(@() PipelineManager.addOutput(info, ...
                'invalidOptionalEffect', 'UnknownDataType', 'file', '', ...
                'optional.txt', 3, 'isData', false, ...
                'returnsValue', false, 'isRequired', false), ...
                'addOutput:InvalidOptionalOutput');
        end

        function testOptionalEmptyFileOutputDoesNotRegisterStaleFile(testCase)
            parserFolder = createTemporaryTemplateCategory(testCase.ProjectRoot);
            cleanup = onCleanup(@() cleanupTemplateCategory( ...
                parserFolder, testCase.ProjectRoot));

            stalePath = fullfile(testCase.SaveFolder, 'optional_result.txt');
            iWriteTextFile(stalePath, 'stale output');

            pm = PipelineManager(testCase.SaveFolder, '', testCase.ProjectRoot);
            pm.b_skipSteps = false;
            pm.addStep('pmOptionalEmptyFileOutput');
            result = pm.executePipeline('PrintSummary', false);

            testCase.verifyEqual(result.status, "completed");
            testCase.verifyFalse(any(strcmpi(result.createdFiles.FileName, ...
                'optional_result.txt')));
            row = result.outputManifest(strcmpi( ...
                result.outputManifest.OutputName, 'optionalFile'), :);
            testCase.verifyNumElements(row.FileExists, 1);
            testCase.verifyFalse(row.IsRequired);
            testCase.verifyFalse(row.FileExists);
            testCase.verifyEqual(row.ActualPersistence, "not_produced");
            testCase.verifyEqual(fileread(stalePath), 'stale output');
        end

        function testOptionalCreatedFileOutputIsRegistered(testCase)
            parserFolder = createTemporaryTemplateCategory(testCase.ProjectRoot);
            cleanup = onCleanup(@() cleanupTemplateCategory( ...
                parserFolder, testCase.ProjectRoot));

            pm = PipelineManager(testCase.SaveFolder, '', testCase.ProjectRoot);
            pm.b_skipSteps = false;
            pm.addStep('pmOptionalCreatedFileOutput');
            result = pm.executePipeline('PrintSummary', false);

            testCase.verifyEqual(result.status, "completed");
            testCase.verifyTrue(any(strcmpi(result.createdFiles.FileName, ...
                'optional_result.txt')));
            row = result.outputManifest(strcmpi( ...
                result.outputManifest.OutputName, 'optionalFile'), :);
            testCase.verifyNumElements(row.FileExists, 1);
            testCase.verifyFalse(row.IsRequired);
            testCase.verifyTrue(row.FileExists);
            testCase.verifyEqual(row.ActualFileName, "optional_result.txt");
        end

        function testMixedReturnedValueAndFileEffectsExecute(testCase)
            parserFolder = createTemporaryTemplateCategory(testCase.ProjectRoot);
            cleanup = onCleanup(@() cleanupTemplateCategory( ...
                parserFolder, testCase.ProjectRoot));

            pm = PipelineManager(testCase.SaveFolder, '', testCase.ProjectRoot);
            pm.b_skipSteps = false;
            pm.addStep('pmMixedReturnedAndFileEffect', ...
                'input', 'input.dat');

            result = pm.executePipeline('PrintSummary', false);

            testCase.verifyEqual(result.status, "completed");
            testCase.verifyEqual(sort(result.createdFiles.FileName), ...
                sort(["effect_one.txt"; "effect_two.log"; "mixedOutput.dat"]));
            testCase.verifyTrue(all(result.createdFiles.FileExists));
            testCase.verifyEqual(loadData(fullfile( ...
                testCase.SaveFolder, 'mixedOutput.dat')), ...
                testCase.InputData + single(1));
        end

        function testMissingDeclaredFileEffectFailsStep(testCase)
            parserFolder = createTemporaryTemplateCategory(testCase.ProjectRoot);
            cleanup = onCleanup(@() cleanupTemplateCategory( ...
                parserFolder, testCase.ProjectRoot));

            pm = PipelineManager(testCase.SaveFolder, '', testCase.ProjectRoot);
            pm.b_skipSteps = false;
            pm.addStep('pmMissingDeclaredFileEffect');

            result = pm.executePipeline('PrintSummary', false);

            testCase.verifyEqual(result.status, "failed");
            testCase.verifySubstring( ...
                string(result.globalPipeLog.Messages_long{1}), ...
                'declared file pattern "never_created_*.txt", but it matched no files');
            testCase.verifyEmpty(dir(fullfile( ...
                testCase.SaveFolder, 'never_created_*.txt')));
        end

        function testInvalidRequiredInputsFailClearly(testCase)
            missingFolder = fullfile(testCase.SaveFolder, 'missing');

            testCase.verifyError(@() funcTemplate( ...
                testCase.InputData(:,:,1), testCase.SaveFolder), ...
                'Umitoolbox:funcTemplate:InvalidData');
            testCase.verifyError(@() funcTemplate( ...
                testCase.InputData, missingFolder), ...
                'Umitoolbox:funcTemplate:InvalidSaveFolder');

            testCase.verifyError(@() funcTemplateUMT( ...
                testCase.InputData, testCase.SaveFolder, ...
                'EntryName', 'not valid'), ...
                'Umitoolbox:funcTemplateUMT:InvalidEntryName');
            testCase.verifyError(@() funcTemplateUMT( ...
                testCase.InputData, missingFolder), ...
                'Umitoolbox:funcTemplateUMT:InvalidSaveFolder');

            testCase.verifyError(@() funcTemplateFileManifest( ...
                1i, testCase.SaveFolder), ...
                'Umitoolbox:funcTemplateFileManifest:InvalidData');
            testCase.verifyError(@() funcTemplateFileManifest( ...
                testCase.InputData, missingFolder), ...
                'Umitoolbox:funcTemplateFileManifest:InvalidSaveFolder');
        end

        function testPipelineManagerInjectsSaveFolder(testCase)
            parserFolder = createTemporaryTemplateCategory(testCase.ProjectRoot);
            cleanup = onCleanup(@() cleanupTemplateCategory( ...
                parserFolder, testCase.ProjectRoot));

            pm = PipelineManager(testCase.SaveFolder, '', testCase.ProjectRoot);
            pm.b_skipSteps = false;
            pm.addStep('funcTemplate', ...
                'input', 'input.dat', ...
                'saveas', 'pipelineTemplate.dat');
            pm.executePipeline();

            outputPath = fullfile(testCase.SaveFolder, ...
                'pipelineTemplate.dat');
            testCase.verifyTrue(isfile(outputPath));
            testCase.verifyEqual(loadData(outputPath), testCase.InputData);
        end

        function testLeafPolicyReusesMatchingPermanentFile(testCase)
            parserFolder = createTemporaryTemplateCategory(testCase.ProjectRoot);
            cleanup = onCleanup(@() cleanupTemplateCategory( ...
                parserFolder, testCase.ProjectRoot));

            pm = PipelineManager(testCase.SaveFolder, '', testCase.ProjectRoot);
            pm.b_skipSteps = false;
            pm.addStep('pmLeafMatchingFile', 'input', 'input.dat');

            result = pm.executePipeline('PrintSummary', false);

            testCase.verifyEqual(result.status, "completed");
            testCase.verifyTrue(isfile(fullfile( ...
                testCase.SaveFolder, 'declared.dat')));
            testCase.verifyFalse(isfile(fullfile( ...
                testCase.SaveFolder, 'declared_1.dat')));
            testCase.verifyEqual(loadData(fullfile( ...
                testCase.SaveFolder, 'declared.dat')), testCase.InputData);
        end

        function testLeafPolicyPreservesForeignCollision(testCase)
            parserFolder = createTemporaryTemplateCategory(testCase.ProjectRoot);
            cleanup = onCleanup(@() cleanupTemplateCategory( ...
                parserFolder, testCase.ProjectRoot));

            foreignData = single(99 .* ones(size(testCase.InputData)));
            saveData(fullfile(testCase.SaveFolder, 'declared.dat'), foreignData, ...
                'DimNames', {'Y', 'X', 'T'}, 'FrameRateHz', 20);

            pm = PipelineManager(testCase.SaveFolder, '', testCase.ProjectRoot);
            pm.b_skipSteps = false;
            pm.addStep('pmLeafMismatchedFile', 'input', 'input.dat');

            result = pm.executePipeline('PrintSummary', false);

            testCase.verifyEqual(result.status, "completed");
            testCase.verifyEqual(loadData(fullfile( ...
                testCase.SaveFolder, 'declared.dat')), foreignData);
            testCase.verifyEqual(loadData(fullfile( ...
                testCase.SaveFolder, 'actual.dat')), testCase.InputData);
            testCase.verifyEqual(loadData(fullfile( ...
                testCase.SaveFolder, 'declared_1.dat')), testCase.InputData);
            testCase.verifyFalse(isfile(fullfile( ...
                testCase.SaveFolder, 'declared_2.dat')));
        end

        function testGlobalLogIncludesExtendedFailureMessage(testCase)
            parserFolder = createTemporaryTemplateCategory(testCase.ProjectRoot);
            cleanup = onCleanup(@() cleanupTemplateCategory( ...
                parserFolder, testCase.ProjectRoot));

            pm = PipelineManager(testCase.SaveFolder, '', testCase.ProjectRoot);
            pm.b_skipSteps = false;
            pm.addStep('funcTemplate', ...
                'input', 'input.dat', ...
                'saveas', 'pipelineTemplate.dat');

            % Leave an unreadable file (headerless, no sidecar: rejected since
            % .dat header Phase 5b) so loading fails during the folder run,
            % where PipelineManager records the MException.
            fid = fopen(fullfile(testCase.SaveFolder, 'input.dat'), 'w');
            testCase.assertNotEqual(fid, -1);
            fileCleanup = onCleanup(@() safeFclose(fid));
            fwrite(fid, single(1), 'single');
            fclose(fid);
            clear fileCleanup

            result = pm.executePipeline('PrintSummary', false);

            testCase.verifyTrue(ismember( ...
                'Messages_long', result.globalPipeLog.Properties.VariableNames));
            shortMessage = string(result.globalPipeLog.Messages_short{1});
            longMessage = string(result.globalPipeLog.Messages_long{1});
            testCase.verifySubstring(longMessage, ...
                'has no header');
            testCase.verifyGreaterThan(strlength(longMessage), ...
                strlength(shortMessage));
        end

        function testTemplatesPassCodeAnalyzer(testCase)
            files = cellfun(@(name) fullfile(testCase.ProjectRoot, ...
                'Analysis', [name '.m']), testCase.templateNames(), ...
                'UniformOutput', false);
            files{end+1} = mfilename('fullpath');

            for iFile = 1:numel(files)
                issues = checkcode(files{iFile}, '-id');
                testCase.verifyEmpty(issues, ...
                    sprintf('Code Analyzer issues in %s.', files{iFile}));
            end
        end
    end

    methods (Static, Access = private)
        function names = templateNames()
            names = { ...
                'funcTemplate', ...
                'funcTemplateUMT', ...
                'funcTemplateFileManifest', ...
                'funcTemplateAcquisitionInitializer', ...
                'funcTemplateAcquisitionCompanion'};
        end
    end
end

function parserFolder = createTemporaryTemplateCategory(projectRoot)
%CREATETEMPORARYTEMPLATECATEGORY Copy templates into a scanned category.

analysisFolder = fullfile(projectRoot, 'Analysis');
parserFolder = tempname(analysisFolder);
mkdir(parserFolder);

templateNames = { ...
    'funcTemplate', ...
    'funcTemplateUMT', ...
    'funcTemplateFileManifest', ...
    'funcTemplateAcquisitionInitializer', ...
    'funcTemplateAcquisitionCompanion'};
for iTemplate = 1:numel(templateNames)
    fileName = [templateNames{iTemplate} '.m'];
    copyfile(fullfile(analysisFolder, fileName), ...
        fullfile(parserFolder, fileName));
end

testFolder = fileparts(mfilename('fullpath'));
fixtureFolder = fullfile(testFolder, 'fixtures');
fixtureNames = { ...
    'pmLeafMatchingFile', ...
    'pmLeafMismatchedFile', ...
    'pmMixedReturnedAndFileEffect', ...
    'pmMissingDeclaredFileEffect', ...
    'pmOptionalEmptyFileOutput', ...
    'pmOptionalCreatedFileOutput'};
for iFixture = 1:numel(fixtureNames)
    fileName = [fixtureNames{iFixture} '.m'];
    copyfile(fullfile(fixtureFolder, fileName), ...
        fullfile(parserFolder, fileName));
end

addpath(parserFolder, '-begin');
clear funcTemplate funcTemplateUMT funcTemplateFileManifest ...
    funcTemplateAcquisitionInitializer funcTemplateAcquisitionCompanion ...
    pmLeafMatchingFile pmLeafMismatchedFile ...
    pmMixedReturnedAndFileEffect pmMissingDeclaredFileEffect ...
    pmOptionalEmptyFileOutput pmOptionalCreatedFileOutput
rehash;
end

function cleanupTemplateCategory(parserFolder, projectRoot)
%CLEANUPTEMPLATECATEGORY Remove the temporary scanned category.

analysisFolder = fullfile(projectRoot, 'Analysis');
expectedPrefix = [analysisFolder filesep];
if ~startsWith(parserFolder, expectedPrefix, 'IgnoreCase', true)
    error('Umitoolbox:TestFunctionTemplate:UnsafeCleanupPath', ...
        'Refusing to remove a folder outside Analysis.');
end

clear funcTemplate funcTemplateUMT funcTemplateFileManifest ...
    funcTemplateAcquisitionInitializer funcTemplateAcquisitionCompanion ...
    pmLeafMatchingFile pmLeafMismatchedFile ...
    pmMixedReturnedAndFileEffect pmMissingDeclaredFileEffect ...
    pmOptionalEmptyFileOutput pmOptionalCreatedFileOutput
if any(strcmp(strsplit(path, pathsep), parserFolder))
    rmpath(parserFolder);
end
if isfolder(parserFolder)
    rmdir(parserFolder, 's');
end
rehash;
end

function iWriteTextFile(filePath, contents)
fid = fopen(filePath, 'w');
assert(fid ~= -1, 'Failed to create test file: %s', filePath);
cleanupObj = onCleanup(@() fclose(fid));
fwrite(fid, contents, 'char');
end
