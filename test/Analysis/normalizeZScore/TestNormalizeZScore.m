classdef TestNormalizeZScore < matlab.unittest.TestCase
    %TESTNORMALIZEZSCORE Unit tests for normalizeZScore.
    %
    % Coverage:
    %   1) pipelineInfo generation
    %   2) Raw 3-D YXT array input
    %   3) Raw .dat filename input
    %   4) Event-split YXTE array: each E slice normalized along T
    %   5) Event-split YXTE .dat input matches the array result
    %   6) Constant-trace handling (std == 0)
    %   7) Rejection of UMT struct and .umt file input
    %   8) Rejection of unsupported .dat layouts
    %
    % Fixture policy:
    %   - Sample files are copied from:
    %         <projectRoot>\test\Analysis\TestingData_with_events
    %   - A fresh temporary SaveFolder is used per test.

    properties
        SaveFolder char
        ProjectRoot char
        SampleDataFolder char
    end

    methods (TestClassSetup)
        function addProjectPathAndLocateFixture(testCase)
            thisFile = mfilename('fullpath');
            testFolder = fileparts(thisFile);
            projectRoot = extractBefore(testFolder, [filesep 'test']);

            if isempty(projectRoot)
                projectRoot = fileparts(fileparts(testFolder));
            end

            testCase.ProjectRoot = char(projectRoot);
            addpath(genpath(testCase.ProjectRoot));

            cfg = [];
            if isappdata(0, 'NormalizeZScoreTestConfig')
                cfg = getappdata(0, 'NormalizeZScoreTestConfig');
            end

            if ~isempty(cfg) && isfield(cfg, 'sampleDataFolder') && ...
                    isfolder(cfg.sampleDataFolder)
                testCase.SampleDataFolder = char(string(cfg.sampleDataFolder));
            else
                testCase.SampleDataFolder = fullfile( ...
                    testCase.ProjectRoot, ...
                    'test', ...
                    'Analysis', ...
                    'TestingData_with_events');
            end

            testCase.assertTrue( ...
                isfile(fullfile(testCase.SampleDataFolder, 'green.dat')), ...
                ['Missing fixture file "green.dat". Put it in: ' ...
                testCase.SampleDataFolder]);

            testCase.assertTrue( ...
                isfile(fullfile(testCase.SampleDataFolder, 'AcqInfos.mat')), ...
                ['Missing fixture file "AcqInfos.mat". Put it in: ' ...
                testCase.SampleDataFolder]);
        end
    end

    methods (TestMethodSetup)
        function createFreshWorkspace(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture
            fx = testCase.applyFixture(TemporaryFolderFixture);
            testCase.SaveFolder = fx.Folder;

            copyfile( ...
                fullfile(testCase.SampleDataFolder, 'green.dat'), ...
                fullfile(testCase.SaveFolder, 'green.dat'));

            copyfile( ...
                fullfile(testCase.SampleDataFolder, 'AcqInfos.mat'), ...
                fullfile(testCase.SaveFolder, 'AcqInfos.mat'));

            deleteIfExists(fullfile(testCase.SaveFolder, 'normZ.dat'));
            deleteIfExists(fullfile(testCase.SaveFolder, 'recording_input.umt'));
        end
    end

    methods (TestMethodTeardown)
        function cleanupCreatedFiles(testCase)
            deleteIfExists(fullfile(testCase.SaveFolder, 'normZ.dat'));
            deleteIfExists(fullfile(testCase.SaveFolder, 'recording_input.umt'));
        end
    end

    methods (Test)

        function testPipelineInfo(testCase)
            info = normalizeZScore('pipelineInfo');

            testCase.verifyTrue(isstruct(info) && isscalar(info));

            reqFields = {'name','description','version','inputs','outputs'};
            testCase.verifyTrue(all(ismember(reqFields, fieldnames(info))), ...
                'pipelineInfo is missing one or more required top-level fields.');

            testCase.verifyEqual(info.name, 'normalizeZScore');
            testCase.verifyEqual(numel(info.inputs), 2);
            testCase.verifyEqual(numel(info.outputs), 1);

            inputNames = {info.inputs.name};
            outputNames = {info.outputs.name};

            testCase.verifyEqual(inputNames, {'data','SaveFolder'});
            testCase.verifyEqual(outputNames, {'outData'});

            dataInput = info.inputs(strcmp(inputNames, 'data'));
            outDataOutput = info.outputs(strcmp(outputNames, 'outData'));

            dataTypes = dataInput.type;
            outTypes = outDataOutput.type;

            if ischar(dataTypes)
                dataTypes = {dataTypes};
            end
            if ischar(outTypes)
                outTypes = {outTypes};
            end

            testCase.verifyTrue(ismember('ImageTimeSeries', dataTypes));
            testCase.verifyTrue(ismember('ProcessedData', dataTypes));
            testCase.verifyTrue(ismember('ImageTimeSeries', outTypes));
            testCase.verifyTrue(ismember('ProcessedData', outTypes));
        end

        function testArrayInput(testCase)
            rawData = single(loadData(fullfile(testCase.SaveFolder, 'green.dat')));
            expected = testCase.expectedZScore(rawData);

            out = normalizeZScore(rawData, testCase.SaveFolder);

            testCase.verifyTrue(isnumeric(out));
            testCase.verifyEqual(size(out), size(rawData));
            testCase.verifyNumericEquivalent(single(out), single(expected));
        end

        function testDatInput(testCase)
            inFile = fullfile(testCase.SaveFolder, 'green.dat');
            rawData = single(loadData(inFile));
            expected = testCase.expectedZScore(rawData);

            out = normalizeZScore(inFile, testCase.SaveFolder);

            testCase.verifyTrue(ischar(out) || (isstring(out) && isscalar(out)));
            outFile = char(string(out));
            testCase.verifyTrue(isfile(outFile));

            actual = single(loadData(outFile));
            testCase.verifyEqual(size(actual), size(rawData));
            testCase.verifyNumericEquivalent(actual, single(expected));
        end

        function testEventSplitArrayNormalizesEachTrialAlongT(testCase)
            trial = testCase.buildTrial();
            % The second trial is an affine copy of the first: per-trial
            % z-scores must therefore be identical. Pooling the statistics
            % over E would not give that.
            rawYXTE = cat(4, trial, 3 .* trial + 7);

            out = normalizeZScore(rawYXTE, testCase.SaveFolder);

            testCase.verifyEqual(size(out), size(rawYXTE));
            testCase.verifyNumericEquivalent(single(out(:,:,:,1)), ...
                single(testCase.expectedZScore(trial)));
            testCase.verifyNumericEquivalent(single(out(:,:,:,2)), ...
                single(out(:,:,:,1)));
        end

        function testEventSplitDatMatchesArray(testCase)
            trial = testCase.buildTrial();
            rawYXTE = cat(4, trial, 3 .* trial + 7, -trial);
            rate = loadMetaData(fullfile(testCase.SaveFolder, 'green.dat')).frameRateHz;
            inFile = fullfile(testCase.SaveFolder, 'yxte_input.dat');
            saveData(inFile, rawYXTE, 'DimNames', {'Y','X','T','E'}, ...
                'FrameRateHz', rate);

            outFile = char(string(normalizeZScore(inFile, testCase.SaveFolder)));

            testCase.verifyEqual(outFile, fullfile(testCase.SaveFolder, 'normZ.dat'));
            testCase.verifyEqual(loadMetaData(outFile).dimNames, {'Y','X','T','E'});
            testCase.verifyNumericEquivalent(single(loadData(outFile)), ...
                single(normalizeZScore(rawYXTE, testCase.SaveFolder)));
        end

        function testConstantTraceArray(testCase)
            data = ones(8, 7, 20, 'single');
            out = normalizeZScore(data, testCase.SaveFolder);

            testCase.verifyTrue(isnumeric(out));
            testCase.verifyEqual(size(out), size(data));
            testCase.verifyEqual(out, zeros(size(data), 'like', out));
        end

        function testRejectsUMTStructAndFile(testCase)
            rawData = single(loadData(fullfile(testCase.SaveFolder, 'green.dat')));
            umt = genUMTStruct(rawData, ...
                'kind', 'image', ...
                'entryName', 'main', ...
                'dimNames', {'Y','X','T'});
            umtFile = fullfile(testCase.SaveFolder, 'recording_input.umt');
            saveData(umtFile, umt);

            testCase.verifyError( ...
                @() normalizeZScore(umt, testCase.SaveFolder), ...
                'normalizeZScore:UnsupportedInputType');
            testCase.verifyError( ...
                @() normalizeZScore(umtFile, testCase.SaveFolder), ...
                'normalizeZScore:UnsupportedInputFile');
        end

        function testRejectsUnsupportedDatLayout(testCase)
            rawData = single(loadData(fullfile(testCase.SaveFolder, 'green.dat')));
            rate = loadMetaData(fullfile(testCase.SaveFolder, 'green.dat')).frameRateHz;
            inFile = fullfile(testCase.SaveFolder, 'yxe_input.dat');
            saveData(inFile, rawData(:,:,1:3), 'DimNames', {'Y','X','E'}, ...
                'FrameRateHz', rate);

            testCase.verifyError( ...
                @() normalizeZScore(inFile, testCase.SaveFolder), ...
                'Umitoolbox:normalizeZScore:unsupportedLayout');
        end

        function testStandardVsLowRAM(testCase)
            rawData = single(loadData(fullfile(testCase.SaveFolder, 'green.dat')));

            outStandard = normalizeZScore(rawData, testCase.SaveFolder);
            outLowRAMFile = normalizeZScore(fullfile(testCase.SaveFolder, 'green.dat'), ...
                testCase.SaveFolder);

            testCase.verifyTrue(ischar(outLowRAMFile) || (isstring(outLowRAMFile) && isscalar(outLowRAMFile)));
            testCase.verifyTrue(isfile(char(string(outLowRAMFile))));

            outLowRAM = single(loadData(char(string(outLowRAMFile))));

            testCase.verifyEqual(size(outStandard), size(outLowRAM));

            if isequaln(outStandard, outLowRAM)
                return
            end

            diffVals = double(outStandard(:)) - double(outLowRAM(:));
            diffVals = diffVals(isfinite(diffVals));

            if isempty(diffVals)
                testCase.verifyEqual(outStandard, outLowRAM);
                return
            end

            testCase.verifyLessThanOrEqual(std(diffVals, 0, 'omitnan'), 1e-4);
        end
    end

    methods (Access = private)

        function expected = expectedZScore(~, dataIn)
            % Reference: z-score every Y x X x T slice (each E slice) on its own.
            origSz = size(dataIn);
            data2D = reshape(single(dataIn), [], origSz(3), prod(origSz(4:end)));

            mu  = mean(data2D, 2, 'omitnan');
            sig = std(data2D, 0, 2, 'omitnan');
            sig(sig == 0) = 1;

            expected = reshape((data2D - mu) ./ sig, origSz);
        end

        function trial = buildTrial(testCase)
            % One short Y x X x T trial taken from the fixture recording.
            rawData = single(loadData(fullfile(testCase.SaveFolder, 'green.dat')));
            trial = rawData(:,:,1:min(30, size(rawData, 3)));
        end

        function verifyNumericEquivalent(testCase, a, b)
            testCase.verifyEqual(size(a), size(b));

            if isequaln(a, b)
                return
            end

            diffVals = double(a(:)) - double(b(:));
            diffVals = diffVals(isfinite(diffVals));

            if isempty(diffVals)
                testCase.verifyEqual(a, b);
                return
            end

            testCase.verifyLessThanOrEqual(std(diffVals, 0, 'omitnan'), 1e-4);
        end

        function verifyStandardVsLowRAM(testCase, actual, expected)
            testCase.verifyEqual(size(actual), size(expected));

            if isequaln(actual, expected)
                return
            end

            diffVals = double(actual(:)) - double(expected(:));
            diffVals = diffVals(isfinite(diffVals));

            if isempty(diffVals)
                testCase.verifyEqual(actual, expected);
                return
            end

            testCase.verifyLessThanOrEqual(std(diffVals, 0, 'omitnan'), 1e-4);
        end
    end
end

function deleteIfExists(filePath)

if isfile(filePath)
    delete(filePath);
end
end

