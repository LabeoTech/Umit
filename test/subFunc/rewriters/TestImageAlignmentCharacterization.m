classdef TestImageAlignmentCharacterization < matlab.unittest.TestCase
    %TESTIMAGEALIGNMENTCHARACTERIZATION applyImageAlignmentToFolder keeps its values.
    %
    %   .dat header Phase 4e-2. fixtures/imageAlignmentReference.mat was
    %   recorded by makeImageAlignmentReference from the committed,
    %   headerless applyImageAlignmentToFolder (commit 296abb7); it is
    %   frozen. The case runs on a legacy-sidecar folder (headerless with
    %   sidecars since .dat header Phase 5a) and on a headered copy
    %   and must reproduce the recorded sizes and hashes of the .dat and
    %   image .umt outputs. Aligned files must be headered at the new frame
    %   size with the input's class, rate, exposure, and name, and
    %   AcqInfos.mat must be left unchanged. Integer data keep their class.

    properties
        Folder
        CleanupFcn = @() []
    end

    methods (TestMethodTeardown)
        function removeFolder(testCase)
            testCase.CleanupFcn();
            if ~isempty(testCase.Folder) && isfolder(testCase.Folder)
                rmdir(testCase.Folder, 's');
            end
        end
    end

    methods (Test)
        function legacySidecarInputMatchesReference(testCase)
            testCase.verifyAgainstReference('legacySidecar');
        end

        function headeredInputMatchesReference(testCase)
            testCase.verifyAgainstReference('single');
        end

        function integerInputKeepsItsClass(testCase)
            c = imageAlignmentCase();
            testCase.newFolder(c, 'uint16');
            inputs = testCase.readInputs(c);
            view = imref2d(size(c.referenceImage));

            evalc('c.run(testCase.Folder)');

            for j = 1:numel(c.datOutputs)
                f = fullfile(testCase.Folder, c.datOutputs{j});
                [~, out] = evalc('loadData(f)');
                in = inputs{j};
                expected = zeros([size(c.referenceImage), size(in, 3)], 'uint16');
                for t = 1:size(in, 3)
                    expected(:, :, t) = imwarp(in(:, :, t), c.tform, 'OutputView', view, ...
                        'InterpolationMethod', 'linear', 'FillValues', uint16(0));
                end
                testCase.verifyClass(out, 'uint16');
                testCase.verifyEqual(out, expected, c.datOutputs{j});
                testCase.verifyEqual(readDatHeader(f).dataClass, 'uint16');
            end
        end

        function alignedFolderCanBeAlignedAgain(testCase)
            % AcqInfos.mat keeps the raw size, so the second alignment must
            % take the old size from the files themselves.
            c = imageAlignmentCase();
            testCase.newFolder(c, 'legacySidecar');
            evalc('c.run(testCase.Folder)');
            evalc('c.run(testCase.Folder)');

            for j = 1:numel(c.datOutputs)
                hdr = readDatHeader(fullfile(testCase.Folder, c.datOutputs{j}));
                testCase.verifyEqual(hdr.dimSizes(1:2), size(c.referenceImage));
                testCase.verifyTrue(hdr.writeComplete);
            end
        end

        function inconsistentSizesAreRefusedBeforeAnyChange(testCase)
            c = imageAlignmentCase();
            testCase.newFolder(c, 'single');
            % Crop red.dat by one row so the folder's frame sizes differ.
            f = fullfile(testCase.Folder, 'red.dat');
            [~, info] = evalc('loadMetaData(f)');
            [~, values] = evalc('loadData(f)');
            values = values(1:end-1, :, :);
            hdr = datHeaderFromInfo(info, 'red', 'dimSizes', size(values));
            hdr.writeComplete = true;
            fid = fopen(f, 'w', 'ieee-le');
            fwrite(fid, encodeDatHeader(hdr), 'uint8');
            fwrite(fid, values, 'single');
            fclose(fid);
            before = iSnapshot(testCase.Folder);

            testCase.verifyError(@() c.run(testCase.Folder), ...
                'ImageAlignmentTool:InconsistentDatSizes');
            testCase.verifyEqual(iSnapshot(testCase.Folder), before);
        end
    end

    methods (Access = private)
        function newFolder(testCase, c, dataClass)
            testCase.Folder = fullfile(tempdir, ['imageAlignChar_' char(java.util.UUID.randomUUID)]);
            testCase.CleanupFcn = c.build(testCase.Folder, dataClass);
        end

        function inputs = readInputs(testCase, c)
            inputs = cell(1, numel(c.datOutputs));
            for j = 1:numel(c.datOutputs)
                f = fullfile(testCase.Folder, c.datOutputs{j});
                inputs{j} = loadData(f);
            end
        end

        function verifyAgainstReference(testCase, dataClass)
            c = imageAlignmentCase();
            S = load(fullfile(fileparts(mfilename('fullpath')), 'fixtures', ...
                'imageAlignmentReference.mat'), 'ref');
            ref = S.ref;
            testCase.newFolder(c, dataClass);

            sources = cell(1, numel(c.datOutputs));
            for j = 1:numel(c.datOutputs)
                f = fullfile(testCase.Folder, c.datOutputs{j});
                sources{j} = loadMetaData(f);
            end
            acqBefore = load(fullfile(testCase.Folder, 'AcqInfos.mat'));

            [~, report] = evalc('c.run(testCase.Folder)');

            for j = 1:numel(c.datOutputs)
                f = fullfile(testCase.Folder, c.datOutputs{j});
                [~, values] = evalc('loadData(f)');
                testCase.verifyEqual(size(values), ref.datSizes{j}, [c.datOutputs{j} ' size']);
                testCase.verifyEqual(c.hash(values), ref.datSha256{j}, [c.datOutputs{j} ' values']);

                testCase.assertTrue(isDatWithHeader(f), [c.datOutputs{j} ' must be headered']);
                hdr = readDatHeader(f);
                src = sources{j};
                testCase.verifyTrue(hdr.writeComplete);
                testCase.verifyEqual(hdr.dataClass, src.dataClass);
                testCase.verifyEqual(hdr.dimSizes, [size(c.referenceImage), src.dimSizes(3)]);
                testCase.verifyEqual(hdr.frameRateHz, double(single(src.frameRateHz)));
                testCase.verifyEqual(hdr.exposureMsec, double(single(src.exposureMsec)));
                [~, base] = fileparts(f);
                testCase.verifyEqual(hdr.channelName, base);
            end

            U = load(fullfile(testCase.Folder, c.umtOutput), '-mat');
            testCase.verifyEqual(size(U.umt.data.map.value), ref.umtSize);
            testCase.verifyEqual(c.hash(single(U.umt.data.map.value)), ref.umtSha256);

            D = load(fullfile(testCase.Folder, 'DataParams.mat'), 'DataParams');
            testCase.verifyEqual(D.DataParams.view.imageSizeYX, ref.viewImageSizeYX);
            testCase.verifyEqual(load(fullfile(testCase.Folder, 'AcqInfos.mat')), acqBefore, ...
                'AcqInfos.mat must be left unchanged');
            testCase.verifyEqual([acqBefore.AcqInfoStream.Height, acqBefore.AcqInfoStream.Width], ...
                ref.acqHeightWidthBefore);
            testCase.verifyFalse(any(report.metadataFilesUpdated == "AcqInfos.mat"));
            testCase.verifyNotEmpty(dir(fullfile(testCase.Folder, 'manualAlignmentReport_*.mat')));
        end
    end
end

function snapshot = iSnapshot(folder)
% Relative path and checksum of every file under FOLDER.
files = dir(fullfile(folder, '**', '*'));
files = files(~[files.isdir]);
names = strings(numel(files), 1);
sums = strings(numel(files), 1);
for k = 1:numel(files)
    full = fullfile(files(k).folder, files(k).name);
    names(k) = string(extractAfter(full, strlength(folder)));
    sums(k) = string(computeFileChecksum(full));
end
[names, order] = sort(names);
snapshot = table(names, sums(order), 'VariableNames', {'File', 'Checksum'});
end
