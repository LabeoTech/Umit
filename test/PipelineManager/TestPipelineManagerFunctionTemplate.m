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
            save(fullfile(testCase.SaveFolder, 'AcqInfos.mat'), ...
                'AcqInfoStream');

            saveData(fullfile(testCase.SaveFolder, 'input.dat'), ...
                testCase.InputData);
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
            testCase.verifyEqual(entry.meta.FrameRateHz, 20);

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
        end

        function testInvalidRequiredInputsFailClearly(testCase)
            missingFolder = fullfile(testCase.SaveFolder, 'missing');

            testCase.verifyError(@() funcTemplate( ...
                testCase.InputData(:,:,1), testCase.SaveFolder), ...
                'umIToolbox:funcTemplate:InvalidData');
            testCase.verifyError(@() funcTemplate( ...
                testCase.InputData, missingFolder), ...
                'umIToolbox:funcTemplate:InvalidSaveFolder');

            testCase.verifyError(@() funcTemplateUMT( ...
                testCase.InputData, testCase.SaveFolder, ...
                'EntryName', 'not valid'), ...
                'umIToolbox:funcTemplateUMT:InvalidEntryName');
            testCase.verifyError(@() funcTemplateUMT( ...
                testCase.InputData, missingFolder), ...
                'umIToolbox:funcTemplateUMT:InvalidSaveFolder');

            testCase.verifyError(@() funcTemplateFileManifest( ...
                1i, testCase.SaveFolder), ...
                'umIToolbox:funcTemplateFileManifest:InvalidData');
            testCase.verifyError(@() funcTemplateFileManifest( ...
                testCase.InputData, missingFolder), ...
                'umIToolbox:funcTemplateFileManifest:InvalidSaveFolder');
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
                'funcTemplateFileManifest'};
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
    'funcTemplateFileManifest'};
for iTemplate = 1:numel(templateNames)
    fileName = [templateNames{iTemplate} '.m'];
    copyfile(fullfile(analysisFolder, fileName), ...
        fullfile(parserFolder, fileName));
end

addpath(parserFolder, '-begin');
clear funcTemplate funcTemplateUMT funcTemplateFileManifest
rehash;
end

function cleanupTemplateCategory(parserFolder, projectRoot)
%CLEANUPTEMPLATECATEGORY Remove the temporary scanned category.

analysisFolder = fullfile(projectRoot, 'Analysis');
expectedPrefix = [analysisFolder filesep];
if ~startsWith(parserFolder, expectedPrefix, 'IgnoreCase', true)
    error('umIToolbox:TestFunctionTemplate:UnsafeCleanupPath', ...
        'Refusing to remove a folder outside Analysis.');
end

clear funcTemplate funcTemplateUMT funcTemplateFileManifest
if any(strcmp(strsplit(path, pathsep), parserFolder))
    rmpath(parserFolder);
end
if isfolder(parserFolder)
    rmdir(parserFolder, 's');
end
rehash;
end
