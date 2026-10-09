classdef TestRunConvertToTiff < matlab.unittest.TestCase
    %TESTRUNCONVERTTOTIFF Unit tests for run_ConvertToTiff.
    %
    % Contract: .dat input only (Y-X-T, Y-X-T-E, or a single Y-X frame);
    % streamed in frame blocks; event-split files get one TIFF per E slice
    % labeled from events.mat.

    methods (Test)
        function testPipelineInfo(testCase)
            info = run_ConvertToTiff('pipelineInfo');
            testCase.verifyTrue(isstruct(info) && isscalar(info));
            testCase.verifyEqual(info.name, 'run_ConvertToTiff');
            testCase.verifyEqual(info.inputs(1).dataMode, 'file');
            testCase.verifyEqual(info.outputs(1).defOutfilename, 'img_out.tif');
        end

        function testDatFileInputUsesInputStemAndKeepsPixels(testCase)
            saveFolder = iTempFolder(testCase);
            data = reshape(single(1:(12*10*4)), 12, 10, 4);
            datFile = fullfile(saveFolder, 'sourceData.dat');
            writeTestDat(datFile, data, 10);

            outFile = run_ConvertToTiff(datFile, saveFolder);

            testCase.verifyEqual(outFile, {'img_sourceData.tif'});
            testCase.verifyEqual(iListExportFiles(saveFolder), {'img_sourceData.tif'});
            tifPath = fullfile(saveFolder, 'img_sourceData.tif');
            testCase.verifyNumElements(imfinfo(tifPath), 4);
            for k = 1:4
                testCase.verifyEqual(single(imread(tifPath, k)), data(:,:,k), ...
                    sprintf('page %d', k));
            end
        end

        function testBareFileNameIsResolvedInSaveFolder(testCase)
            saveFolder = iTempFolder(testCase);
            writeTestDat(fullfile(saveFolder, 'bare.dat'), rand(6, 5, 3, 'single'), 10);

            outFile = run_ConvertToTiff('bare.dat', saveFolder);

            testCase.verifyEqual(outFile, {'img_bare.tif'});
        end

        function testSingleFrameDatIsOnePage(testCase)
            saveFolder = iTempFolder(testCase);
            img = reshape(single(1:30), 6, 5);
            writeTestDat(fullfile(saveFolder, 'frame.dat'), img, 10, 'DimNames', {'Y','X'});

            outFile = run_ConvertToTiff(fullfile(saveFolder, 'frame.dat'), saveFolder);

            testCase.verifyEqual(outFile, {'img_frame.tif'});
            tifPath = fullfile(saveFolder, 'img_frame.tif');
            testCase.verifyNumElements(imfinfo(tifPath), 1);
            testCase.verifyEqual(single(imread(tifPath)), img);
        end

        function testEventSplitDatWritesPerEventTiffAndInfo(testCase)
            saveFolder = iEventsFolder(testCase);   % events: A, B, A (3 instances)
            vals = rand(8, 9, 3, 3, 'single');
            datFile = fullfile(saveFolder, 'trials.dat');
            writeTestDat(datFile, vals, 10, 'DimNames', {'Y','X','T','E'});

            outFile = run_ConvertToTiff(datFile, saveFolder);

            expected = {'img_trials_C1_R1.tif', 'img_trials_C2_R1.tif', ...
                'img_trials_C1_R2.tif', 'img_trials_info.txt'};
            testCase.verifyEqual(outFile, expected);
            testCase.verifyEqual(iListExportFiles(saveFolder), sort(expected));

            % Each TIFF holds the frames of its own E slice.
            for iE = 1:3
                tifPath = fullfile(saveFolder, expected{iE});
                testCase.verifyNumElements(imfinfo(tifPath), 3);
                for k = 1:3
                    testCase.verifyEqual(single(imread(tifPath, k)), vals(:,:,k,iE), ...
                        sprintf('event %d page %d', iE, k));
                end
            end

            txt = fileread(fullfile(saveFolder, 'img_trials_info.txt'));
            testCase.verifySubstring(txt, 'img_trials_C1_R1.tif,A,1');
            testCase.verifySubstring(txt, 'img_trials_C2_R1.tif,B,1');
            testCase.verifySubstring(txt, 'img_trials_C1_R2.tif,A,2');
        end

        function testAggregatedEventSplitUsesConditionLabels(testCase)
            saveFolder = iEventsFolder(testCase);   % 2 conditions: A, B
            datFile = fullfile(saveFolder, 'agg.dat');
            writeTestDat(datFile, rand(8, 9, 3, 2, 'single'), 10, 'DimNames', {'Y','X','T','E'});

            outFile = run_ConvertToTiff(datFile, saveFolder);

            testCase.verifyEqual(outFile, {'img_agg_C1_R0.tif', 'img_agg_C2_R0.tif', 'img_agg_info.txt'});
        end

        function testEventSplitWithoutEventsFileUsesOneCondition(testCase)
            saveFolder = iTempFolder(testCase);
            datFile = fullfile(saveFolder, 'noev.dat');
            writeTestDat(datFile, rand(8, 9, 3, 2, 'single'), 10, 'DimNames', {'Y','X','T','E'});

            outFile = run_ConvertToTiff(datFile, saveFolder);

            testCase.verifyEqual(outFile, ...
                {'img_noev_C1_R1.tif', 'img_noev_C1_R2.tif', 'img_noev_info.txt'});
        end

        function testEventSplitThatDoesNotMatchEventsIsRefused(testCase)
            saveFolder = iEventsFolder(testCase);   % 3 instances, 2 conditions
            datFile = fullfile(saveFolder, 'bad.dat');
            writeTestDat(datFile, rand(8, 9, 3, 4, 'single'), 10, 'DimNames', {'Y','X','T','E'});

            testCase.verifyError(@() run_ConvertToTiff(datFile, saveFolder), ...
                'Umitoolbox:run_ConvertToTiff:EventMappingMismatch');
            testCase.verifyEmpty(iListExportFiles(saveFolder));
        end

        function testForcedMultiBlockStreamingMatchesSingleBlock(testCase)
            % The memory mock forces many frame blocks; the pages written
            % must equal the source frames, for a stack and for event slices.
            testCase.assumeEqual(computer, 'PCWIN64', ...
                'Chunk forcing relies on shadowing the PCWIN64 memory() built-in.');
            saveFolder = iEventsFolder(testCase);
            stack = rand(20, 20, 400, 'single');
            events = rand(20, 20, 300, 3, 'single');
            writeTestDat(fullfile(saveFolder, 'stack.dat'), stack, 10);
            writeTestDat(fullfile(saveFolder, 'split.dat'), events, 10, ...
                'DimNames', {'Y','X','T','E'});

            projectRoot = extractBefore(mfilename('fullpath'), [filesep 'test' filesep]);
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile( ...
                projectRoot, 'test', 'subFunc', 'calculateMaxChunkSize', 'mocks')));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_MODE', 'fixed'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_TOTAL', '10000'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_AVAILABLE', '5000'));
            testCase.assertGreaterThan(calculateMaxChunkSize(20*20*400*4, 2, 0.1), 1, ...
                'The fixture must force more than one block.');

            outStack = {};
            progress = evalc("outStack = run_ConvertToTiff(fullfile(saveFolder, 'stack.dat'), saveFolder);");
            testCase.verifySubstring(progress, 'frames written');
            tifPath = fullfile(saveFolder, outStack{1});
            testCase.assertNumElements(imfinfo(tifPath), 400);
            for k = [1 2 137 399 400]
                testCase.verifyEqual(single(imread(tifPath, k)), stack(:,:,k));
            end

            outSplit = run_ConvertToTiff(fullfile(saveFolder, 'split.dat'), saveFolder);
            tifPath = fullfile(saveFolder, outSplit{3});
            testCase.assertNumElements(imfinfo(tifPath), 300);
            for k = [1 150 300]
                testCase.verifyEqual(single(imread(tifPath, k)), events(:,:,k,3));
            end
        end

        function testRejectsUnsupportedInputs(testCase)
            saveFolder = iTempFolder(testCase);

            % Numeric array
            testCase.verifyError(@() run_ConvertToTiff(rand(6, 5, 3, 'single'), saveFolder), ...
                'Umitoolbox:run_ConvertToTiff:UnsupportedInputType');

            % UMT struct
            umt = genUMTStruct(rand(6, 5, 3, 'single'), 'kind', 'image', ...
                'entryName', 'main', 'dimNames', {'Y','X','T'});
            testCase.verifyError(@() run_ConvertToTiff(umt, saveFolder), ...
                'Umitoolbox:run_ConvertToTiff:UnsupportedInputType');

            % .umt file
            umtFile = fullfile(saveFolder, 'x.umt');
            saveData(umtFile, umt);
            testCase.verifyError(@() run_ConvertToTiff(umtFile, saveFolder), ...
                'Umitoolbox:run_ConvertToTiff:UnsupportedInputFile');

            % Missing file
            testCase.verifyError(@() run_ConvertToTiff('missing.dat', saveFolder), ...
                'Umitoolbox:run_ConvertToTiff:InputFileNotFound');

            % Layout without a T axis but with E
            yxe = fullfile(saveFolder, 'yxe.dat');
            writeTestDat(yxe, rand(6, 5, 3, 'single'), 10, 'DimNames', {'Y','X','E'});
            testCase.verifyError(@() run_ConvertToTiff(yxe, saveFolder), ...
                'Umitoolbox:run_ConvertToTiff:unsupportedLayout');
        end
    end
end

function folder = iTempFolder(testCase)
folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
end

function folder = iEventsFolder(testCase)
%IEVENTSFOLDER Folder with events.mat: conditions A, B, A (3 instances).
folder = iTempFolder(testCase);
writeMinimalAcqInfosMat(folder, 'FrameRateHz', 10, 'AISampleRate', 100);
ev = EventsManager(folder, '', 'csv');
ev.EventFileParseMethod = 'csv';
ev.getTriggersFromSignal(makePulseSignal(400, [50 150 250], 10, 'Amplitude', 5), 100, false);
testCase.assertTrue(ev.readConditionFile( ...
    writeCSVConditionFile(folder, 'labels.csv', {'A', 'B', 'A'})));
ev.saveEvents(folder);
end

function names = iListExportFiles(folderPath)
items = [dir(fullfile(folderPath, '*.tif')); ...
    dir(fullfile(folderPath, '*_info.txt'))];
names = sort({items.name});
end
