classdef TestSplitDataByEvent < matlab.unittest.TestCase
    properties
        ProjectRoot
        SampleDataFolder
        TempFolder
        SourceInfo
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

            cfg = [];
            if isappdata(0, 'SplitDataByEventTestConfig')
                cfg = getappdata(0, 'SplitDataByEventTestConfig');
            end

            if ~isempty(cfg) && isfield(cfg, 'sampleDataFolder') && ...
                    isfolder(cfg.sampleDataFolder)
                testCase.SampleDataFolder = char(string(cfg.sampleDataFolder));
            else
                testCase.SampleDataFolder = fullfile( ...
                    testCase.ProjectRoot, ...
                    'test', 'Analysis', 'TestingData_with_events');
            end

            assert(isfolder(testCase.SampleDataFolder), ...
                'Sample data folder was not found.');

            testCase.TempFolder = fullfile(tempdir, ...
                ['TestSplitDataByEvent_' char(java.util.UUID.randomUUID)]);
            mkdir(testCase.TempFolder);

            iCopyIfExists(fullfile(testCase.SampleDataFolder, 'AcqInfos.mat'), testCase.TempFolder);
            iCopyIfExists(fullfile(testCase.SampleDataFolder, 'events.mat'), testCase.TempFolder);
            iCopyIfExists(fullfile(testCase.SampleDataFolder, 'green.dat'), testCase.TempFolder);

            acqFile = fullfile(testCase.TempFolder, 'AcqInfos.mat');
            assert(isfile(acqFile), ...
                'AcqInfos.mat was not copied into the temporary folder.');

            md = load(acqFile, 'AcqInfoStream');
            testCase.SourceInfo = md.AcqInfoStream;
        end
    end

    methods (TestMethodTeardown)
        function removeTempFolder(testCase)
            if ~isempty(testCase.TempFolder) && isfolder(testCase.TempFolder)
                try %#ok<TRYNC>
                    rmdir(testCase.TempFolder, 's');
                end
            end
        end
    end

    methods (Test)
        function testArrayInputIsRejected(testCase)
            dataYXT = iLoadNumericData(testCase.TempFolder);

            testCase.verifyError(@() split_data_by_event(dataYXT, testCase.TempFolder, ...
                'FrameRateHz', double(testCase.SourceInfo.FrameRateHz)), ...
                'Umitoolbox:split_data_by_event:invalidInputType');
        end

        function testRawDatInput(testCase)
            % A .dat filename runs LOW-RAM MODE: the result is the path of
            % the Y-X-T-E .dat written in SaveFolder.
            outFile = split_data_by_event(fullfile(testCase.TempFolder, 'green.dat'), ...
                testCase.TempFolder);

            testCase.verifyEqual(outFile, fullfile(testCase.TempFolder, 'dataByEv.dat'));
            testCase.verifyTrue(isfile(outFile));
            testCase.verifyEqual(loadMetaData(outFile).dimNames, {'Y','X','T','E'});
        end

        function testLowRamDatMatchesInRamSplit(testCase)
            % The streamed split must equal the in-RAM split of the same
            % data, trial for trial (same frames, same crop to the shortest
            % trial, single class), and leave no scratch file behind.
            datFile = fullfile(testCase.TempFolder, 'green.dat');
            rate = double(testCase.SourceInfo.FrameRateHz);
            expected = iSplitInRam(testCase, iLoadNumericData(testCase.TempFolder), rate);

            outFile = split_data_by_event(datFile, testCase.TempFolder, 'FrameRateHz', rate);

            info = loadMetaData(outFile);
            testCase.verifyEqual(info.dimNames, {'Y','X','T','E'});
            testCase.verifyEqual(info.dimSizes(:).', size(expected));
            testCase.verifyEqual(info.dataClass, 'single');
            testCase.verifyEqual(loadData(outFile), expected);
            testCase.verifyEmpty(dir(fullfile(testCase.TempFolder, '*_writing.dat')));
        end

        function testRawDatUsesTheFileOwnFrameRate(testCase)
            % The file's own rate differs from the folder rate: .dat input
            % must carry the file's rate, not AcqInfoStream's. Since .dat
            % header Phase 5a the file is a legacy (headerless) file whose
            % sidecar declares its rate; the headered case is the next test.
            acqFile = fullfile(testCase.TempFolder, 'AcqInfos.mat');
            datFile = fullfile(testCase.TempFolder, 'green.dat');
            S = load(acqFile, 'AcqInfoStream');
            channelRate = 0.5 * double(S.AcqInfoStream.FrameRateHz);
            writeTestDat(datFile, loadData(datFile), channelRate, 'Format', 'legacySidecar');

            testCase.assertEqual(loadMetaData(datFile).frameRateHz, channelRate, ...
                'Precondition: loadMetaData resolves the file''s own rate.');

            outFile = split_data_by_event(datFile, testCase.TempFolder);
            expected = iSplitInRam(testCase, iLoadNumericData(testCase.TempFolder), channelRate);
            testCase.verifyEqual(loadData(outFile), expected);
        end

        function testHeaderedDatUsesTheHeaderFrameRate(testCase)
            datFile = fullfile(testCase.TempFolder, 'green.dat');
            info = loadMetaData(datFile);
            data = loadData(datFile);
            headerRate = 7.5;
            hdr = struct('dataClass', 'single', 'frameRateHz', headerRate, ...
                'exposureMsec', NaN, 'channelName', 'green', ...
                'dimNames', {{'Y', 'X', 'T'}}, 'dimSizes', info.dimSizes, ...
                'writeComplete', true);
            fid = fopen(datFile, 'w', 'ieee-le');
            fwrite(fid, encodeDatHeader(hdr), 'uint8');
            fwrite(fid, data, 'single');
            fclose(fid);
            testCase.assertTrue(isDatWithHeader(datFile));

            outFile = split_data_by_event(datFile, testCase.TempFolder);
            expected = iSplitInRam(testCase, single(data), headerRate);
            testCase.verifyEqual(loadData(outFile), expected);
        end

        function testUMTStructContinuousInput(testCase)
            dataYXT = iLoadNumericData(testCase.TempFolder);

            umt = genUMTStruct(single(dataYXT), ...
                'kind', 'image', ...
                'entryName', 'main', ...
                'dimNames', {'Y','X','T'}, ...
                'meta', struct('FrameRateHz', double(testCase.SourceInfo.FrameRateHz)));

            % UMT input keeps a UMT output with the shared eventInfo.
            out = split_data_by_event(umt, testCase.TempFolder);

            iVerifyEventSplitUMT(testCase, out, dataYXT);
        end

        function testUMTFileContinuousInput(testCase)
            dataYXT = iLoadNumericData(testCase.TempFolder);

            umt = genUMTStruct(single(dataYXT), ...
                'kind', 'image', ...
                'entryName', 'main', ...
                'dimNames', {'Y','X','T'}, ...
                'meta', struct('FrameRateHz', double(testCase.SourceInfo.FrameRateHz)));
            saveData(fullfile(testCase.TempFolder, 'input.umt'), umt);

            out = split_data_by_event(fullfile(testCase.TempFolder, 'input.umt'), ...
                testCase.TempFolder);

            iVerifyEventSplitUMT(testCase, out, dataYXT);
        end

        function testNonImageRoiUMTIsSplitAlongTime(testCase)
            % Any UMT kind works as long as one entry has a T axis: an
            % ROI x T entry becomes ROI x T x E, keeping kind, entry name,
            % and the ROI labels.
            rate = double(testCase.SourceInfo.FrameRateHz);
            dataYXT = iLoadNumericData(testCase.TempFolder);
            nT = size(dataYXT, 3);
            roiTraces = single(reshape(dataYXT(1:3, 1, :), 3, nT));

            umt = genUMTStruct(roiTraces, ...
                'kind', 'roi', ...
                'entryName', 'traces', ...
                'dimNames', {'ROI','T'}, ...
                'labels', struct('ROI', {{'r1','r2','r3'}}), ...
                'meta', struct('FrameRateHz', rate));

            out = split_data_by_event(umt, testCase.TempFolder);

            validateUMTStruct(out, 'requireEventInfo', true);
            testCase.verifyEqual(lower(char(string(out.kind))), 'roi');
            testCase.verifyEqual(fieldnames(out.data), {'traces'});
            testCase.verifyEqual(out.data.traces.dimNames, {'ROI','T','E'});
            expected = iSplitInRam(testCase, reshape(roiTraces, 3, 1, nT), rate);
            testCase.verifyEqual(out.data.traces.value, ...
                reshape(expected, 3, size(expected, 3), size(expected, 4)));
            testCase.verifyEqual(cellstr(string(out.labels.ROI(:))).', {'r1','r2','r3'});
            testCase.verifyEqual(numel(out.eventInfo.eventID), size(expected, 4));
        end

        function testImageUMTKeepsItsEntryName(testCase)
            rate = double(testCase.SourceInfo.FrameRateHz);
            dataYXT = iLoadNumericData(testCase.TempFolder);
            umt = genUMTStruct(dataYXT, ...
                'kind', 'image', ...
                'entryName', 'green', ...
                'dimNames', {'Y','X','T'}, ...
                'meta', struct('FrameRateHz', rate));

            out = split_data_by_event(umt, testCase.TempFolder);

            testCase.verifyEqual(fieldnames(out.data), {'green'});
            testCase.verifyEqual(out.data.green.dimNames, {'Y','X','T','E'});
        end

        function testUMTWithSeveralTimeEntriesIsRejected(testCase)
            rate = double(testCase.SourceInfo.FrameRateHz);
            dataYXT = iLoadNumericData(testCase.TempFolder);
            umt = genUMTStruct(dataYXT, ...
                'kind', 'image', ...
                'entryName', 'a', ...
                'dimNames', {'Y','X','T'}, ...
                'meta', struct('FrameRateHz', rate));
            umt = genUMTStruct(umt, ...
                'value', dataYXT, ...
                'entryName', 'b', ...
                'dimNames', {'Y','X','T'}, ...
                'meta', struct('FrameRateHz', rate));

            testCase.verifyError(@() split_data_by_event(umt, testCase.TempFolder), ...
                'Umitoolbox:split_data_by_event:multipleCompatibleUMTEntries');
        end

        function testUMTWithoutTimeDimensionIsRejected(testCase)
            umt = genUMTStruct(single(rand(4, 5)), ...
                'kind', 'image', ...
                'entryName', 'main', ...
                'dimNames', {'Y','X'}, ...
                'meta', struct('FrameRateHz', double(testCase.SourceInfo.FrameRateHz)));

            testCase.verifyError(@() split_data_by_event(umt, testCase.TempFolder), ...
                'Umitoolbox:split_data_by_event:noTimeDimension');
        end

        function testIgnoredRepetitionIsSavedAndFlagged(testCase)
            % .dat header Phase 8b/8c: every event instance is saved, ignored
            % ones included; once saved as .dat, the mapping onto events.mat
            % flags them and gives each instance's duration.
            ev = iEv(testCase);
            inst = ev.getEventInstances();
            testCase.assumeGreaterThanOrEqual(numel(inst.eventID), 2);
            cond = char(inst.eventName(1));
            ev.removeRepetition(cond, inst.repetitionIndex(1));
            ev.saveEvents(testCase.TempFolder);

            byEv = split_data_by_event(fullfile(testCase.TempFolder, 'green.dat'), ...
                testCase.TempFolder);
            testCase.verifyEqual(double(loadMetaData(byEv).dimSizes(4)), numel(inst.eventID));
            m = resolveDatEventMapping(loadMetaData(byEv), testCase.TempFolder);
            testCase.verifyEqual(m.status, 'matched');
            expectedSelected = inst.selected;
            expectedSelected(inst.eventID == inst.eventID(1) & ...
                inst.repetitionIndex == inst.repetitionIndex(1)) = false;
            testCase.verifyEqual(m.eventInfo.selected, expectedSelected);
            testCase.verifyEqual(m.eventInfo.durationSec, inst.durationSec, 'AbsTol', 1e-6);
        end

        function testRejectAlreadyEventSplitUMT(testCase)
            dataYXT = iLoadNumericData(testCase.TempFolder);
            dataYXTE = iEv(testCase).splitDataByEvents(dataYXT, 'FrameRateHz', iEv(testCase).AcqInfo.FrameRateHz);

            umt = genUMTStruct(single(dataYXTE), ...
                'kind', 'image', ...
                'entryName', 'main', ...
                'dimNames', {'Y','X','T','E'}, ...
                'SaveFolder', testCase.TempFolder);
            umt = appendUMTEventInfo(umt, ...
                'eventInfo', iEv(testCase).exportEventInfo('FrameRateHz', iEv(testCase).AcqInfo.FrameRateHz), ...
                'overwrite', true);

            testCase.verifyError(@() split_data_by_event(umt, testCase.TempFolder), ...
                'Umitoolbox:split_data_by_event:alreadyEventSplit');
        end

        function testMissingEventsMat(testCase)
            delete(fullfile(testCase.TempFolder, 'events.mat'));

            testCase.verifyError(@() split_data_by_event( ...
                fullfile(testCase.TempFolder, 'green.dat'), testCase.TempFolder), ...
                'Umitoolbox:split_data_by_event:missingEventsFile');
        end

        function testUMTWithoutFrameRateErrors(testCase)
            % A UMT entry without meta.FrameRateHz needs 'FrameRateHz';
            % AcqInfos.mat in the folder is not used for it.
            testCase.assertTrue(isfile(fullfile(testCase.TempFolder, 'AcqInfos.mat')));
            umt = genUMTStruct(iLoadNumericData(testCase.TempFolder), ...
                'kind', 'image', 'entryName', 'main', 'dimNames', {'Y','X','T'});

            testCase.verifyError(@() split_data_by_event(umt, testCase.TempFolder), ...
                'Umitoolbox:split_data_by_event:missingFrameRateHz');
        end

        function testForcedMultiSlabMatchesInRamSplit(testCase)
            % The memory mock forces many X slabs per trial; the streamed
            % split must still equal the in-RAM one.
            testCase.assumeEqual(computer, 'PCWIN64', ...
                'Chunk forcing relies on shadowing the PCWIN64 memory() built-in.');
            datFile = fullfile(testCase.TempFolder, 'green.dat');
            rate = double(testCase.SourceInfo.FrameRateHz);
            dataYXT = iLoadNumericData(testCase.TempFolder);
            expected = iSplitInRam(testCase, dataYXT, rate);

            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile( ...
                testCase.ProjectRoot, 'test', 'subFunc', 'calculateMaxChunkSize', 'mocks')));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_MODE', 'fixed'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_TOTAL', '10000'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_AVAILABLE', '5000'));
            trialBytes = size(dataYXT, 1) * size(dataYXT, 2) * size(expected, 3) * 4;
            testCase.assertGreaterThan(calculateMaxChunkSize(trialBytes, 2, 0.2), 1, ...
                'The fixture must force more than one X slab per trial.');

            outFile = split_data_by_event(datFile, testCase.TempFolder, 'FrameRateHz', rate);

            testCase.verifyEqual(loadData(outFile), expected);
        end

        function testRejectsUnsupportedInputs(testCase)
            sv = testCase.TempFolder;

            testCase.verifyError(@() split_data_by_event('missing.dat', sv), ...
                'Umitoolbox:split_data_by_event:missingInputFile');

            matFile = fullfile(sv, 'x.mat');
            fclose(fopen(matFile, 'w'));
            testCase.verifyError(@() split_data_by_event(matFile, sv), ...
                'Umitoolbox:split_data_by_event:unsupportedInputFile');

            % Already event-split .dat files, and a layout without T.
            layouts = {{'Y','X','T','E'}, [6 5 4 2]; {'Y','X','E'}, [6 5 3]; {'Y','X'}, [6 5]};
            for k = 1:size(layouts, 1)
                f = fullfile(sv, 'layout.dat');
                writeTestDat(f, rand([layouts{k, 2}, 1], 'single'), 10, ...
                    'DimNames', layouts{k, 1});
                testCase.verifyError(@() split_data_by_event(f, sv), ...
                    'Umitoolbox:split_data_by_event:unsupportedLayout');
            end
        end

        function testExplicitFrameRateDifferentFromDatHeaderWarns(testCase)
            % An explicit 'FrameRateHz' that differs from the .dat header
            % warns sourceInfoConflict; the explicit value is used.
            datFile = fullfile(testCase.TempFolder, 'green.dat');
            headerRate = double(loadMetaData(datFile).frameRateHz);
            explicitRate = 0.5 * headerRate;   % events stay inside the recording

            testCase.verifyWarning(@() split_data_by_event(datFile, ...
                testCase.TempFolder, 'FrameRateHz', explicitRate), ...
                'Umitoolbox:split_data_by_event:sourceInfoConflict');
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture( ...
                'Umitoolbox:split_data_by_event:sourceInfoConflict'));
            outFile = split_data_by_event(datFile, testCase.TempFolder, 'FrameRateHz', explicitRate);
            testCase.verifyTrue(isfile(outFile));
        end

        function testPipelineInfo(testCase)
            info = split_data_by_event('pipelineInfo');
            testCase.verifyEqual(info.name, 'split_data_by_event');
            testCase.verifyFalse(info.legacyOpts);
            dataInput = info.inputs(strcmp({info.inputs.name}, 'data'));
            testCase.verifyEqual(dataInput.dataMode, 'file');
            testCase.verifyNumElements(info.outputs, 1);
            testCase.verifyEqual(info.outputs(1).defOutfilename, 'dataByEv.dat');
            testCase.verifyEqual(info.outputs(1).type, {'ProcessedData'});
        end
    end
end

function dataYXT = iLoadNumericData(saveFolder)
loaded = loadData(fullfile(saveFolder, 'green.dat'));
assert(isnumeric(loaded) && ndims(loaded) == 3, ...
    'Fixture data must resolve to numeric YXT data.');
dataYXT = single(loaded);
end

function iVerifyEventSplitUMT(testCase, out, dataYXT)
%IVERIFYEVENTSPLITUMT A UMT output equals the in-RAM split and carries eventInfo.
validateUMTStruct(out, 'requireEventInfo', true);
testCase.verifyEqual(lower(char(string(out.kind))), 'image');
testCase.verifyEqual(out.data.main.dimNames, {'Y','X','T','E'});
expected = iSplitInRam(testCase, dataYXT, double(testCase.SourceInfo.FrameRateHz));
testCase.verifyEqual(out.data.main.value, expected);
testCase.verifyEqual(numel(out.eventInfo.eventID), size(expected, 4));
end

function out = iSplitInRam(testCase, dataYXT, rate)
%ISPLITINRAM Reference split of a Y-X-T array (every instance, ignored included).
out = single(iEv(testCase).splitDataByEvents(dataYXT, 'FrameRateHz', rate, ...
    'IncludeIgnored', true));
end

function iCopyIfExists(srcFile, dstFolder)
if isfile(srcFile)
    copyfile(srcFile, dstFolder);
end
end

function ev = iEv(testCase)
%IEV EventsManager for the test folder (its AcqInfos.mat rate is the data rate).
ev = EventsManager(testCase.TempFolder);
end
