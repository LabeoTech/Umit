classdef TestUMTInfoIgnoresAcqInfos < matlab.unittest.TestCase
    %TESTUMTINFOIGNORESACQINFOS .umt Info comes from the file, not AcqInfos.mat (Phase 7a).
    %
    %   AcqInfos.mat describes the raw acquisition. Its Height, Width, and
    %   FrameRateHz must not leak into the Info of a .umt file in the same
    %   folder: size comes from the first entry, and the frame rate from the
    %   entry's meta.FrameRateHz when present.

    properties
        Folder
    end

    methods (TestMethodSetup)
        function createFolder(testCase)
            testCase.Folder = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            AcqInfoStream = struct('Height', 99, 'Width', 88, 'Length', 7, ...
                'FrameRateHz', 1, 'Datatype', 'uint16');
            save(fullfile(testCase.Folder, 'AcqInfos.mat'), 'AcqInfoStream');
        end
    end

    methods (Test)
        function entryDescribesSizeAndRate(testCase)
            umt = genUMTStruct(single(ones(4, 3, 5)), 'kind', 'image', ...
                'dimNames', {'Y', 'X', 'T'}, 'meta', struct('FrameRateHz', 25));
            f = fullfile(testCase.Folder, 'data.umt');
            saveData(f, umt);

            Info = loadMetaData(f);

            testCase.verifyEqual(Info.Height, 4);
            testCase.verifyEqual(Info.Width, 3);
            testCase.verifyEqual(Info.Length, 5);
            testCase.verifyEqual(Info.FrameRateHz, 25);
            testCase.verifyEqual(Info.Freq, 25);
            testCase.verifyEqual(Info.Datatype, 'single');
        end

        function noEntryRateMeansNoRate(testCase)
            umt = genUMTStruct(single(ones(4, 3, 5)), 'kind', 'image', ...
                'dimNames', {'Y', 'X', 'T'});
            f = fullfile(testCase.Folder, 'norate.umt');
            saveData(f, umt);

            Info = loadMetaData(f);

            testCase.verifyFalse(isfield(Info, 'FrameRateHz'), ...
                'the AcqInfos.mat rate must not be used');
            testCase.verifyEqual(Info.Height, 4);
        end

        function umtImageSourceGetsTheEntryRate(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))), ...
                'GUI', 'DataViewer')));
            umt = genUMTStruct(single(ones(4, 3, 5)), 'kind', 'image', ...
                'dimNames', {'Y', 'X', 'T'}, 'meta', struct('FrameRateHz', 25));
            f = fullfile(testCase.Folder, 'viewer.umt');
            saveData(f, umt);

            src = UMTImageSource(f, 'main');
            testCase.verifyEqual(src.FrameRateHz, 25);
        end
    end
end
