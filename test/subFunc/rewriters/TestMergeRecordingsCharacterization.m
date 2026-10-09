classdef TestMergeRecordingsCharacterization < matlab.unittest.TestCase
    %TESTMERGERECORDINGSCHARACTERIZATION mergeRecordings keeps its merged data and events.
    %
    %   .dat header Phase 4e-3. fixtures/mergeRecordingsReference.mat was
    %   recorded by makeMergeRecordingsReference from the committed,
    %   sidecar-based mergeRecordings (commit 325137f); it is frozen. Each
    %   case runs on legacy sidecar folders and on copies whose green.dat
    %   are headered (no sidecar) and must reproduce the recorded merged
    %   values, events.mat, and copied AcqInfos.mat. The merged file must be
    %   headered and have no metadata .mat next to it. Incompatible inputs
    %   are refused before any output exists.

    properties (TestParameter)
        caseName = {'mergedEvents', 'ignoreEvents', 'trialNamesPermuted', 'missingFileSkipped'}
        inputKind = {'legacy', 'headered'}
        property = {'frameSize', 'dataClass', 'frameRate', 'exposure'}
    end

    properties
        Root
    end

    methods (TestMethodSetup)
        function createRoot(testCase)
            testCase.Root = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
        end
    end

    methods (Test)
        function outputsMatchCommittedReference(testCase, caseName, inputKind)
            [cases, utils] = mergeRecordingsCases();
            c = cases(strcmp({cases.name}, caseName));
            S = load(fullfile(fileparts(mfilename('fullpath')), 'fixtures', ...
                'mergeRecordingsReference.mat'), 'ref');
            ref = S.ref(strcmp({S.ref.name}, caseName));

            utils.build(testCase.Root, strcmp(inputKind, 'headered'));
            c.prep(testCase.Root);
            outFile = utils.outFile(testCase.Root);
            evalc('c.run(testCase.Root, outFile)');

            values = loadData(outFile);
            testCase.verifyEqual(size(values), ref.dimSizes);
            testCase.verifyEqual(utils.hash(values), ref.sha256);
            testCase.verifyEqual(utils.readEvents(fileparts(outFile)), ref.events);
            testCase.verifyEqual(load(fullfile(fileparts(outFile), 'AcqInfos.mat')), ref.acqInfos);

            testCase.assertTrue(isDatWithHeader(outFile), 'the merged file must be headered');
            hdr = readDatHeader(outFile);
            testCase.verifyTrue(hdr.writeComplete);
            testCase.verifyEqual(hdr.dimNames, {'Y', 'X', 'T'});
            testCase.verifyEqual(hdr.dimSizes, ref.dimSizes);
            testCase.verifyEqual(hdr.dataClass, 'single');
            testCase.verifyEqual(hdr.frameRateHz, 10);
            testCase.verifyEqual(hdr.channelName, 'merged');
            testCase.verifyFalse(isfile(strrep(outFile, '.dat', '.mat')), ...
                'no metadata .mat may be written next to the merged file');
        end

        function mergedChannelManifestLengthIsUpdated(testCase)
            % .dat header Phase 7b: the copied AcqInfos.mat entry for the
            % merged file gets the merged frame count; others are kept.
            [~, utils] = mergeRecordingsCases();
            utils.build(testCase.Root, true);
            folders = utils.folders(testCase.Root);
            outFile = utils.outFile(testCase.Root);
            [~, outName, outExt] = fileparts(outFile);
            acqPath = fullfile(folders{end}, 'AcqInfos.mat');
            S = load(acqPath, 'AcqInfoStream');
            AcqInfoStream = rmfield(S.AcqInfoStream, intersect( ...
                fieldnames(S.AcqInfoStream), {'ImportedChannels'}));
            AcqInfoStream = appendImportedChannelInfo(AcqInfoStream, ...
                struct('DatFile', {[outName outExt], 'other.dat'}, ...
                'Length', {1, 7}, 'FrameRateHz', {10, 10}));
            save(acqPath, 'AcqInfoStream');

            evalc('mergeRecordings(outFile, folders, ''green'')');

            S = load(fullfile(fileparts(outFile), 'AcqInfos.mat'), 'AcqInfoStream');
            channels = S.AcqInfoStream.ImportedChannels;
            idx = strcmp({channels.DatFile}, [outName outExt]);
            testCase.verifyEqual(channels(idx).Length, ...
                datAxisSize(loadMetaData(outFile), 'T'));
            testCase.verifyEqual(channels(strcmp({channels.DatFile}, 'other.dat')).Length, 7);
        end

        function incompatibleInputsAreRefused(testCase, property)
            [~, utils] = mergeRecordingsCases();
            utils.build(testCase.Root, true);
            iAlterRecording(fullfile(testCase.Root, 'rec2', 'green.dat'), property);
            outFolder = fileparts(utils.outFile(testCase.Root));
            before = dir(outFolder);

            testCase.verifyError(@() mergeRecordings(utils.outFile(testCase.Root), ...
                utils.folders(testCase.Root), 'green'), ...
                'Umitoolbox:mergeRecordings:incompatibleInputs');
            testCase.verifyFalse(isfile(utils.outFile(testCase.Root)));
            testCase.verifyEqual({dir(outFolder).name}, {before.name});
        end
    end
end

function iAlterRecording(f, property)
% Rewrite a headered recording so one property differs from the others.
info = loadMetaData(f);
values = loadData(f);
hdr = datHeaderFromInfo(info, 'green');
switch property
    case 'frameSize'
        values = values(1:end-1, :, :);
        hdr.dimSizes = size(values);
    case 'dataClass'
        values = uint16(values);
        hdr.dataClass = 'uint16';
    case 'frameRate'
        hdr.frameRateHz = 20;
    case 'exposure'
        hdr.exposureMsec = 7;
end
hdr.writeComplete = true;
fid = fopen(f, 'w', 'ieee-le');
fwrite(fid, encodeDatHeader(hdr), 'uint8');
fwrite(fid, values, hdr.dataClass);
fclose(fid);
end
