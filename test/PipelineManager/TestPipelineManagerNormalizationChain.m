classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestPipelineManagerNormalizationChain < matlab.unittest.TestCase
    %TESTPIPELINEMANAGERNORMALIZATIONCHAIN Zero-centered input warning inside a pipeline.
    %
    %   normalizeBSLN followed by normalizeLPF (default Normalize = true) used
    %   to run silently and give spikes of arbitrary size, because the LPF
    %   divides by a baseline that is about 0 after normalizeBSLN. The LPF now
    %   warns (Umitoolbox:normalizeLPF:ZeroCenteredInput) when it runs inside
    %   the pipeline; the run still completes. Synthetic non-square data.

    properties (Constant)
        LPFWarningId = 'Umitoolbox:normalizeLPF:ZeroCenteredInput'
        FrameRateHz = 10
        NumFrames = 1400   % 140 s: longer than the default baseline period (120.5 s)
    end

    properties
        ProjectRoot char
        SaveFolder char
    end

    methods (TestMethodSetup)
        function createFixtures(testCase)
            testFolder = fileparts(mfilename('fullpath'));
            testCase.ProjectRoot = extractBefore(testFolder, [filesep 'test']);

            fixture = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.SaveFolder = fullfile(fixture.Folder, 'SaveFolder');
            mkdir(testCase.SaveFolder);

            AcqInfoStream = struct('Width', 9, 'Height', 6, ...
                'Length', testCase.NumFrames, 'FrameRateHz', testCase.FrameRateHz, ...
                'ExposureMsec', 1);
            AcqInfoStream = appendImportedChannelInfo(AcqInfoStream, struct( ...
                'DatFile', 'input.dat', 'Length', testCase.NumFrames, ...
                'FrameRateHz', testCase.FrameRateHz));
            save(fullfile(testCase.SaveFolder, 'AcqInfos.mat'), 'AcqInfoStream');

            [yy, xx, tt] = ndgrid(1:6, 1:9, 1:testCase.NumFrames);
            raw = single(1000 + 10 * yy + 5 * xx + 20 * sin(2 * pi * tt / 17 + 0.3 * yy));
            saveData(fullfile(testCase.SaveFolder, 'input.dat'), raw, ...
                'DimNames', {'Y', 'X', 'T'}, 'FrameRateHz', testCase.FrameRateHz);
        end
    end

    methods (Test)
        function lpfAfterBaselineNormalizationWarns(testCase)
            pm = PipelineManager(testCase.SaveFolder, '', testCase.ProjectRoot);
            pm.b_skipSteps = false;
            tagBSLN = pm.addStep('normalizeBSLN', 'input', 'input.dat');
            pm.addStep('normalizeLPF', 'input', sprintf('%s:outData', tagBSLN));

            testCase.verifyWarning(@() pm.executePipeline('PrintSummary', false), ...
                testCase.LPFWarningId);
        end

        function lpfOnRawDataDoesNotWarn(testCase)
            pm = PipelineManager(testCase.SaveFolder, '', testCase.ProjectRoot);
            pm.b_skipSteps = false;
            pm.addStep('normalizeLPF', 'input', 'input.dat');

            output = evalc('pm.executePipeline(''PrintSummary'', false);');

            testCase.verifyEmpty(strfind(output, 'looks zero-centered'));
        end

        function chainStillCompletes(testCase)
            % The warning is not an error: the run reports success.
            pm = PipelineManager(testCase.SaveFolder, '', testCase.ProjectRoot);
            pm.b_skipSteps = false;
            tagBSLN = pm.addStep('normalizeBSLN', 'input', 'input.dat');
            pm.addStep('normalizeLPF', 'input', sprintf('%s:outData', tagBSLN));

            warning('off', testCase.LPFWarningId);
            cleanup = onCleanup(@() warning('on', testCase.LPFWarningId));
            result = pm.executePipeline('PrintSummary', false);
            clear cleanup

            testCase.verifyEqual(result.status, "completed");
        end
    end
end
