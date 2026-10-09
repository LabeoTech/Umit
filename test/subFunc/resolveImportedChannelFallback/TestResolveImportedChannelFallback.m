classdef TestResolveImportedChannelFallback < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addProjectPath(testCase)
            thisFile = mfilename('fullpath');
            projectRoot = extractBefore(fileparts(thisFile), [filesep 'test']);
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                projectRoot, 'IncludingSubfolders', true));
        end
    end

    methods (Test)
        function testBinnedLegacySidecarUsesProcessedGeometry(testCase)
            folder = testCase.createLegacyFolder([6 8], 8, [3 4 8]);
            Info = struct( ...
                'dim_names', {{'Y','X','T'}}, ...
                'datSize', [3 4], ...
                'datLength', 8, ...
                'Freq', 7.5, ...
                'Datatype', 'single');
            save(fullfile(folder, 'red.mat'), 'Info');

            loaded = load(fullfile(folder, 'AcqInfos.mat'));
            [resolved, report] = resolveImportedChannelFallback( ...
                loaded.AcqInfoStream, folder);

            testCase.verifyEqual(report.Source, 'fallback');
            testCase.verifyEqual(resolved.ImportedChannels.Length, 8);
            testCase.verifyEqual(resolved.ImportedChannels.FrameRateHz, 7.5);
            testCase.verifyFalse(isLegacySchemaFolder(folder));
        end

        function testAcqInfosBoundChannelIsRejected(testCase)
            % Headerless channel files without a sidecar (binned or not) are
            % no longer described from AcqInfos.mat (.dat header Phase 5b):
            % the rejection propagates, and the PipelineManager guard treats
            % the folder as legacy.
            for acqYX = {[6 8], [3 4]}
                folder = testCase.createLegacyFolder(acqYX{1}, 8, [3 4 8]);
                loaded = load(fullfile(folder, 'AcqInfos.mat'));

                testCase.verifyError(@() resolveImportedChannelFallback( ...
                    loaded.AcqInfoStream, folder), ...
                    'Umitoolbox:loadMetaData:acqInfosBoundUnsupported');
                testCase.verifyTrue(isLegacySchemaFolder(folder));
            end
        end

        function testHeaderedChannelUsesItsHeader(testCase)
            % A legacy folder (no ImportedChannels) whose channel file is
            % headered, e.g. rewritten by a same-size rewriter: Length and
            % rate come from the header, not from AcqInfos.mat Height/Width.
            folder = testCase.createLegacyFolder([6 8], 99, [3 4 8]);
            delete(fullfile(folder, 'red.dat'));
            writeTestDat(fullfile(folder, 'red.dat'), ...
                single(reshape(1:96, 3, 4, 8)), 7.5);
            loaded = load(fullfile(folder, 'AcqInfos.mat'));

            [resolved, report] = resolveImportedChannelFallback( ...
                loaded.AcqInfoStream, folder);

            testCase.verifyEqual(report.Source, 'fallback');
            testCase.verifyEqual(resolved.ImportedChannels.DatFile, 'red.dat');
            testCase.verifyEqual(resolved.ImportedChannels.Length, 8);
            testCase.verifyEqual(resolved.ImportedChannels.FrameRateHz, 7.5);
            testCase.verifyFalse(isLegacySchemaFolder(folder));
        end

        function testMixedHeaderedAndSidecarChannels(testCase)
            % One channel headered, the other still a legacy sidecar file.
            folder = testCase.createLegacyFolder([3 4], 8, [3 4 8]);
            loaded = load(fullfile(folder, 'AcqInfos.mat'));
            AcqInfoStream = loaded.AcqInfoStream;
            AcqInfoStream.Illumination2 = struct('Color', 'green');
            save(fullfile(folder, 'AcqInfos.mat'), 'AcqInfoStream');
            delete(fullfile(folder, 'red.dat'));
            writeTestDat(fullfile(folder, 'red.dat'), ...
                single(reshape(1:96, 3, 4, 8)), 20);
            writeTestDat(fullfile(folder, 'green.dat'), ...
                single(reshape(1:48, 3, 4, 4)), 10, 'Format', 'legacySidecar');

            resolved = resolveImportedChannelFallback(AcqInfoStream, folder);

            testCase.verifyEqual({resolved.ImportedChannels.DatFile}, {'red.dat', 'green.dat'});
            testCase.verifyEqual([resolved.ImportedChannels.Length], [8 4]);
            testCase.verifyEqual([resolved.ImportedChannels.FrameRateHz], [20 10]);
            testCase.verifyFalse(isLegacySchemaFolder(folder));
        end
    end

    methods (Access = private)
        function folder = createLegacyFolder(testCase, acqYX, acqLength, dataSize)
            fixture = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture);
            folder = fixture.Folder;

            AcqInfoStream = struct( ...
                'Height', acqYX(1), ...
                'Width', acqYX(2), ...
                'Length', acqLength, ...
                'FrameRateHz', 20, ...
                'Datatype', 'single', ...
                'MultiCam', false, ...
                'Illumination1', struct('Color', 'red'));
            save(fullfile(folder, 'AcqInfos.mat'), 'AcqInfoStream');
            data = single(reshape(1:prod(dataSize), dataSize));
            fid = fopen(fullfile(folder, 'red.dat'), 'w');
            testCase.assertGreaterThan(fid, 0);
            cleanup = onCleanup(@() fclose(fid));
            fwrite(fid, data, 'single');
        end
    end
end
