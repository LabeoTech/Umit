classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestPipelineManagerSingleDatGuard < matlab.unittest.TestCase
    %TESTPIPELINEMANAGERSINGLEDATGUARD PipelineManager enforces single-precision .dat files.
    %
    %   DFR-20261019-006: umIT supports only single-precision .dat files.
    %   PipelineManager refuses a non-single .dat as a file source (per folder:
    %   one bad folder is skipped, the others run) and fails a step whose
    %   function wrote a non-single .dat itself (custom functions).

    properties
        ProjectRoot char
        FixtureFolder char
        SaveFolder char
    end

    methods (TestMethodSetup)
        function createFixtures(testCase)
            addpath(genpath(fileparts(fileparts(fileparts(mfilename('fullpath'))))));
            testFolder = fileparts(mfilename('fullpath'));
            testCase.ProjectRoot = extractBefore(testFolder, [filesep 'test']);

            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.SaveFolder = fullfile(fixture.Folder, 'SaveFolder');
            mkdir(testCase.SaveFolder);
            testCase.prepareFolder(testCase.SaveFolder, 'single');

            testCase.FixtureFolder = createAnalysisFixtureCategory(testCase.ProjectRoot, ...
                [struct('name', 'pmLeafMatchingFile', 'text', '', 'source', ...
                    fullfile(testFolder, 'fixtures', 'pmLeafMatchingFile.m')), ...
                 struct('name', 'pmWriteIntDatFile', 'source', '', ...
                    'text', iWriterText('pmWriteIntDatFile', 'uint16')), ...
                 struct('name', 'pmWriteSingleDatFile', 'source', '', ...
                    'text', iWriterText('pmWriteSingleDatFile', 'single'))]);
            testCase.addTeardown(@() removeAnalysisFixtureCategory( ...
                testCase.FixtureFolder, testCase.ProjectRoot));
        end
    end

    methods (Test)
        function nonSingleFileSourceIsRefusedWhenAdded(testCase)
            folder = testCase.newFolder('int', 'uint16');
            pm = PipelineManager(folder, '', testCase.ProjectRoot);
            testCase.verifyError(@() pm.ensureFileSourceStep('green.dat'), ...
                'Umitoolbox:Dat:nonSingle');
            testCase.verifyEmpty(pm.nodes);
            testCase.verifyError(@() pm.addStep('pmLeafMatchingFile', 'input', 'green.dat'), ...
                'Umitoolbox:Dat:nonSingle');
        end

        function singleSourcesPass(testCase)
            pm = PipelineManager(testCase.SaveFolder, '', testCase.ProjectRoot);
            [stepTag, ~] = pm.ensureFileSourceStep('green.dat');
            testCase.verifyEqual(stepTag, 'green.dat');

            legacyFolder = testCase.newFolder('legacy', 'legacySidecar');
            pmLegacy = PipelineManager(legacyFolder, '', testCase.ProjectRoot);
            [stepTag, ~] = pmLegacy.ensureFileSourceStep('green.dat');
            testCase.verifyEqual(stepTag, 'green.dat');
        end

        function oneNonSingleFolderIsSkippedAndTheOthersRun(testCase)
            badFolder = testCase.newFolder('bad', 'uint16');
            pm = PipelineManager({testCase.SaveFolder, badFolder}, {'', ''}, testCase.ProjectRoot);
            pm.addStep('pmLeafMatchingFile', 'input', fullfile(testCase.SaveFolder, 'green.dat'));

            [~, report] = pm.validateNodePreflight('throwOnError', false);
            testCase.verifyNumElements(report.folderErrors, 1);
            testCase.verifyEqual(report.folderErrors.saveFolder, badFolder);
            testCase.verifySubstring(report.folderErrors.message, 'single-precision');
            testCase.verifySubstring(report.folderErrors.message, 'uint16');

            pm.executePipeline('LeafOutputPolicy', 'saveLeaves', 'PrintSummary', false);
            testCase.verifyTrue(isfile(fullfile(testCase.SaveFolder, 'declared.dat')), ...
                'The single folder must still run.');
            testCase.verifyFalse(isfile(fullfile(badFolder, 'declared.dat')), ...
                'The non-single folder must not run.');
        end

        function customFunctionWritingNonSingleDatFailsItsStep(testCase)
            pm = PipelineManager(testCase.SaveFolder, '', testCase.ProjectRoot);
            pm.addStep('pmWriteIntDatFile');

            result = pm.executePipeline('PrintSummary', false);

            testCase.verifyNotEqual(result.status, "completed");
            testCase.verifySubstring(testCase.allText(result), 'single-precision');
            testCase.verifySubstring(testCase.allText(result), 'uint16');
        end

        function customFunctionWritingSingleDatIsUnaffected(testCase)
            pm = PipelineManager(testCase.SaveFolder, '', testCase.ProjectRoot);
            pm.addStep('pmWriteSingleDatFile');

            result = pm.executePipeline('PrintSummary', false);

            testCase.verifyEqual(result.status, "completed");
            testCase.verifyTrue(isfile(fullfile(testCase.SaveFolder, 'pmWriteSingleDatFile.dat')));
        end
    end

    methods (Access = private)
        function prepareFolder(~, folder, dataClass)
            data = reshape(1:2*3*8, [2 3 8]);
            if strcmp(dataClass, 'legacySidecar')
                writeTestDat(fullfile(folder, 'green.dat'), single(data), 10, ...
                    'Format', 'legacySidecar');
            else
                writeTestDat(fullfile(folder, 'green.dat'), cast(data, dataClass), 10);
            end
            AcqInfoStream = struct('Width', 3, 'Height', 2, 'Length', 8, ...
                'FrameRateHz', 10, 'ExposureMsec', 1);
            AcqInfoStream = appendImportedChannelInfo(AcqInfoStream, ...
                struct('DatFile', 'green.dat', 'Length', 8, 'FrameRateHz', 10));
            save(fullfile(folder, 'AcqInfos.mat'), 'AcqInfoStream');
        end

        function folder = newFolder(testCase, name, dataClass)
            folder = fullfile(fileparts(testCase.SaveFolder), name);
            mkdir(folder);
            testCase.prepareFolder(folder, dataClass);
        end

        function text = allText(~, result)
            text = evalc('disp(result)');
            for fieldName = fieldnames(result)'
                value = result.(fieldName{1});
                if istable(value)
                    text = [text, evalc('disp(value)')]; %#ok<AGROW>
                    for column = value.Properties.VariableNames
                        if isstring(value.(column{1}))
                            text = [text, char(strjoin(value.(column{1}), ' '))]; %#ok<AGROW>
                        end
                    end
                elseif isstring(value) || ischar(value)
                    text = [text, char(strjoin(string(value), ' '))]; %#ok<AGROW>
                end
            end
        end
    end
end

function text = iWriterText(name, dataClass)
%IWRITERTEXT Custom file-output function writing a .dat of the given class.
text = sprintf([ ...
    'function outFile = %s(SaveFolder)\n' ...
    'if nargin == 1 && (ischar(SaveFolder) || isstring(SaveFolder)) && ...\n' ...
    '        strcmpi(strtrim(char(string(SaveFolder))), ''pipelineInfo'')\n' ...
    '    outFile = localPipelineInfo();\n' ...
    '    return\n' ...
    'end\n' ...
    'outFile = ''%s.dat'';\n' ...
    'writeTestDat(fullfile(SaveFolder, outFile), cast(reshape(1:48, [2 3 8]), ''%s''), 10);\n' ...
    'end\n\n' ...
    'function info = localPipelineInfo()\n' ...
    'info = PipelineManager.createPipelineInfo(''%s'', ''Test fixture writing a .dat.'');\n' ...
    'info = PipelineManager.addInput(info, ''SaveFolder'', ''SaveFolder'', ''Output folder.'', ...\n' ...
    '    ''kind'', ''input'', ''position'', 1, ''callType'', ''positional'', ''isData'', false);\n' ...
    'info = PipelineManager.addOutput(info, ''outFile'', ''ImageTimeSeries'', ''file'', ...\n' ...
    '    ''Written .dat.'', ''%s.dat'', 1, ''isData'', false);\n' ...
    'end\n'], name, name, dataClass, name, name);
end
