classdef TestSameSizeRewritersCharacterization < matlab.unittest.TestCase
    %TESTSAMESIZEREWRITERSCHARACTERIZATION Same-size rewriters keep their values.
    %
    %   .dat header Phase 4e-1. fixtures/sameSizeRewriterReference.mat was
    %   recorded by makeSameSizeRewriterReference from the committed,
    %   headerless applyTform2Cams, correctMotionArtifact, and
    %   applyRegistrationTformOnFolder (commit 7df18d5); it is frozen. Each
    %   case runs on a legacy-sidecar input folder (headerless with sidecars
    %   since .dat header Phase 5a) and on one whose .dat files
    %   are headered, and must reproduce the recorded sizes, value hashes,
    %   and shifts. Outputs must be headered with the input's class, sizes,
    %   frame rate, and exposure, and the output's name as channelName.

    properties (TestParameter)
        % correctMotionArtifact is no longer pinned here: its per-file output
        % and fixed declared name changed on purpose (see TestCorrectMotionArtifact).
        caseName = {'applyTform2Cams_standard', 'applyTform2Cams_RAMsafe', ...
            'applyRegistrationTformOnFolder'}
        inputKind = {'legacySidecar', 'headered'}
    end

    methods (Test)
        function outputsMatchCommittedReference(testCase, caseName, inputKind)
            [cases, utils] = sameSizeRewriterCases();
            c = cases(strcmp({cases.name}, caseName));
            S = load(fullfile(fileparts(mfilename('fullpath')), 'fixtures', ...
                'sameSizeRewriterReference.mat'), 'ref');
            ref = S.ref(strcmp({S.ref.name}, caseName));
            testCase.assertNotEmpty(ref, 'no reference for this case');

            folder = fullfile(tempdir, ['sameSizeChar_' char(java.util.UUID.randomUUID)]);
            cleanupFcn = c.build(folder, strcmp(inputKind, 'headered'));
            testCase.addTeardown(@() iRemove(folder, cleanupFcn));

            inputInfo = containers.Map();
            listing = dir(fullfile(folder, '*.dat'));
            for k = 1:numel(listing)
                [~, info] = evalc('loadMetaData(fullfile(folder, listing(k).name))');
                inputInfo(listing(k).name) = info;
            end

            [~, extra] = evalc('c.run(folder)');
            testCase.verifyEqual(extra, ref.extra, 'extra results (shifts)');

            for j = 1:numel(c.outputs)
                outFile = fullfile(folder, c.outputs{j});
                [~, values] = evalc('loadData(outFile)');
                testCase.verifyEqual(size(values), ref.dimSizes{j}, [c.outputs{j} ' size']);
                testCase.verifyEqual(utils.hash(values), ref.sha256{j}, [c.outputs{j} ' values']);

                if ~ismember(c.outputs{j}, c.rewritten)
                    % Not rewritten (e.g. Camera 1): keeps its input format.
                    testCase.verifyEqual(isDatWithHeader(outFile), strcmp(inputKind, 'headered'), ...
                        [c.outputs{j} ' must be left as it was']);
                    continue
                end
                testCase.assertTrue(isDatWithHeader(outFile), [c.outputs{j} ' must be headered']);
                hdr = readDatHeader(outFile);
                source = inputInfo(strrep(c.outputs{j}, '_MotionCorrected', ''));
                [~, outBase] = fileparts(outFile);
                testCase.verifyTrue(hdr.writeComplete);
                testCase.verifyEqual(hdr.dataClass, source.dataClass);
                testCase.verifyEqual(hdr.dimSizes, source.dimSizes);
                testCase.verifyEqual(hdr.frameRateHz, double(single(source.frameRateHz)));
                testCase.verifyEqual(hdr.exposureMsec, double(single(source.exposureMsec)));
                testCase.verifyEqual(hdr.channelName, outBase);
            end
        end

        function applyTform2CamsIgnoresRawAcqInfosFrameSize(testCase, inputKind)
            % .dat header Phase 7a: after 7b, AcqInfos.mat holds the raw
            % (unbinned, larger) frame size; applyTform2Cams takes the frame
            % size from the files, so the frozen reference still holds.
            caseName = 'applyTform2Cams_standard';
            [cases, utils] = sameSizeRewriterCases();
            c = cases(strcmp({cases.name}, caseName));
            S = load(fullfile(fileparts(mfilename('fullpath')), 'fixtures', ...
                'sameSizeRewriterReference.mat'), 'ref');
            ref = S.ref(strcmp({S.ref.name}, caseName));

            folder = fullfile(tempdir, ['sameSizeChar_' char(java.util.UUID.randomUUID)]);
            cleanupFcn = c.build(folder, strcmp(inputKind, 'headered'));
            testCase.addTeardown(@() iRemove(folder, cleanupFcn));

            acqFile = fullfile(folder, 'AcqInfos.mat');
            A = load(acqFile, 'AcqInfoStream');
            AcqInfoStream = A.AcqInfoStream;
            AcqInfoStream.Height = 2 * AcqInfoStream.Height;
            AcqInfoStream.Width = 2 * AcqInfoStream.Width;
            save(acqFile, 'AcqInfoStream', '-append');

            evalc('c.run(folder)');
            for j = 1:numel(c.outputs)
                [~, values] = evalc('loadData(fullfile(folder, c.outputs{j}))');
                testCase.verifyEqual(utils.hash(values), ref.sha256{j}, [c.outputs{j} ' values']);
            end
        end
    end
end

function iRemove(folder, cleanupFcn)
cleanupFcn();
if isfolder(folder)
    rmdir(folder, 's');
end
end
