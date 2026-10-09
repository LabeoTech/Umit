classdef TestResolveDataInfoValue < matlab.unittest.TestCase
    %TESTRESOLVEDATAINFOVALUE Precedence explicit > data's own value > error (Phase 6b-2).

    properties
        Folder
        File
    end

    methods (TestMethodSetup)
        function createFile(testCase)
            testCase.Folder = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            testCase.File = fullfile(testCase.Folder, 'data.dat');
            saveData(testCase.File, single(reshape(1:24, 2, 3, 4)), 'DimNames', {'Y', 'X', 'T'}, ...
                'Info', struct('frameRateHz', 25, 'exposureMsec', 3));
            % AcqInfos.mat with another rate must never be used.
            AcqInfoStream = struct('FrameRateHz', 10);
            save(fullfile(testCase.Folder, 'AcqInfos.mat'), 'AcqInfoStream');
        end
    end

    methods (Test)
        function explicitValueWins(testCase)
            v = resolveDataInfoValue('frameRateHz', 30, single(ones(2, 3, 4)), 'fx');
            testCase.verifyEqual(v, 30);
        end

        function headerIsUsedForFileData(testCase)
            testCase.verifyEqual(resolveDataInfoValue('frameRateHz', [], testCase.File, 'fx'), 25);
            testCase.verifyEqual(resolveDataInfoValue('frameRateHz', NaN, testCase.File, 'fx'), 25);
            testCase.verifyEqual(resolveDataInfoValue('exposureMsec', [], testCase.File, 'fx'), 3);
            testCase.verifyEqual(resolveDataInfoValue('dimNames', [], string(testCase.File), 'fx'), ...
                {'Y', 'X', 'T'});
        end

        function explicitConflictWarnsAndWins(testCase)
            v = testCase.verifyWarning(@() resolveDataInfoValue('frameRateHz', 30, testCase.File, 'fx'), ...
                'Umitoolbox:fx:sourceInfoConflict');
            testCase.verifyEqual(v, 30);
            testCase.verifyWarningFree(@() resolveDataInfoValue('frameRateHz', 25, testCase.File, 'fx'));
        end

        function inRamDataWithoutValueErrors(testCase)
            data = single(ones(2, 3, 4));
            try
                resolveDataInfoValue('frameRateHz', [], data, 'fx');
                testCase.verifyFail('expected an error');
            catch ME
                testCase.verifyEqual(ME.identifier, 'Umitoolbox:fx:missingFrameRateHz');
                testCase.verifySubstring(ME.message, '''FrameRateHz''');
                testCase.verifySubstring(ME.message, 'AcqInfos.mat is not used');
            end
            testCase.verifyError(@() resolveDataInfoValue('dimNames', [], data, 'fx'), ...
                'Umitoolbox:fx:missingDimNames');
            testCase.verifyError(@() resolveDataInfoValue('exposureMsec', [], data, 'fx'), ...
                'Umitoolbox:fx:missingExposureMsec');
        end

        function ownValueIsUsedForNonFileData(testCase)
            v = resolveDataInfoValue('frameRateHz', [], struct('kind', 'image'), 'fx', ...
                'OwnValue', 12, 'OwnSource', 'the UMT entry meta');
            testCase.verifyEqual(v, 12);
            v = testCase.verifyWarning(@() resolveDataInfoValue('frameRateHz', 20, [], 'fx', ...
                'OwnValue', 12), 'Umitoolbox:fx:sourceInfoConflict');
            testCase.verifyEqual(v, 20);
            testCase.verifyError(@() resolveDataInfoValue('frameRateHz', [], [], 'fx', 'OwnValue', NaN), ...
                'Umitoolbox:fx:missingFrameRateHz');
        end

        function invalidExplicitValueErrors(testCase)
            testCase.verifyError(@() resolveDataInfoValue('frameRateHz', -5, [], 'fx'), ...
                'Umitoolbox:fx:invalidFrameRateHz');
            testCase.verifyError(@() resolveDataInfoValue('dimNames', {'X', 'Y'}, [], 'fx'), ...
                'Umitoolbox:fx:invalidDimNames');
        end

        function unknownFieldIsRejected(testCase)
            testCase.verifyError(@() resolveDataInfoValue('height', 4, [], 'fx'), ...
                'Umitoolbox:resolveDataInfoValue:invalidInput');
        end
    end
end
