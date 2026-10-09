classdef TestPrepareDualCameraCalibrationData < matlab.unittest.TestCase
    %TESTPREPAREDUALCAMERACALIBRATIONDATA Focused calibration preparation tests.

    properties
        TempRoot
    end

    methods (TestMethodSetup)
        function createTemporaryRoot(testCase)
            testCase.TempRoot = tempname;
            mkdir(testCase.TempRoot);
            testCase.addTeardown(@() testCase.removeFolder(testCase.TempRoot));
        end
    end

    methods (Test)
        function testMultiCameraUsesOwnedTemporarySaveFolder(testCase)
            rawFolder = testCase.createRawFolder();
            invocation = struct('RawFolder', '', 'SaveFolder', '', 'ArgumentCount', 0);
            beforeFolder = '';

            returnedFolder = prepareDualCameraCalibrationData(rawFolder, ...
                'InfoReader', @(~) struct('MultiCam', true), ...
                'ClassificationFcn', @recordClassification, ...
                'BeforeClassificationFcn', @recordBeforeClassification, ...
                'TempRoot', testCase.TempRoot);
            cleanupObject = onCleanup(@() testCase.removeFolder(returnedFolder));

            testCase.verifyEqual(invocation.RawFolder, rawFolder);
            testCase.verifyEqual(invocation.SaveFolder, char(returnedFolder));
            testCase.verifyEqual(invocation.ArgumentCount, 2);
            testCase.verifyEqual(beforeFolder, invocation.SaveFolder);
            testCase.verifyTrue(isfolder(returnedFolder));
            testCase.verifyNotEqual(invocation.SaveFolder, rawFolder);

            function recordClassification(varargin)
                invocation.RawFolder = varargin{1};
                invocation.SaveFolder = varargin{2};
                invocation.ArgumentCount = nargin;
            end
            function recordBeforeClassification(folder)
                beforeFolder = folder;
            end
        end

        function testFalseMultiCameraRejectedBeforeClassification(testCase)
            rawFolder = testCase.createRawFolder();
            wasCalled = false;

            testCase.verifyError(@() prepareDualCameraCalibrationData(rawFolder, ...
                'InfoReader', @(~) struct('MultiCam', false), ...
                'ClassificationFcn', @classify, ...
                'TempRoot', testCase.TempRoot), ...
                'DataViewerCoreg2Cams:NotMultiCameraAcquisition');
            testCase.verifyFalse(wasCalled);
            testCase.verifyEmpty(testCase.childFolders());

            function classify(varargin)
                wasCalled = true;
            end
        end

        function testMissingMultiCameraRejectedBeforeClassification(testCase)
            rawFolder = testCase.createRawFolder();
            wasCalled = false;

            testCase.verifyError(@() prepareDualCameraCalibrationData(rawFolder, ...
                'InfoReader', @(~) struct(), ...
                'ClassificationFcn', @classify, ...
                'TempRoot', testCase.TempRoot), ...
                'DataViewerCoreg2Cams:NotMultiCameraAcquisition');
            testCase.verifyFalse(wasCalled);
            testCase.verifyEmpty(testCase.childFolders());

            function classify(varargin)
                wasCalled = true;
            end
        end

        function testUnreadableMetadataRejectedBeforeClassification(testCase)
            rawFolder = testCase.createRawFolder();
            wasCalled = false;

            testCase.verifyError(@() prepareDualCameraCalibrationData(rawFolder, ...
                'InfoReader', @failRead, ...
                'ClassificationFcn', @classify, ...
                'TempRoot', testCase.TempRoot), ...
                'DataViewerCoreg2Cams:UnreadableAcquisitionMetadata');
            testCase.verifyFalse(wasCalled);
            testCase.verifyEmpty(testCase.childFolders());

            function info = failRead(~)
                info = struct(); %#ok<NASGU>
                error('Test:Unreadable', 'Unreadable metadata.');
            end
            function classify(varargin)
                wasCalled = true;
            end
        end

        function testClassificationFailureRemovesTemporaryFolder(testCase)
            rawFolder = testCase.createRawFolder();

            testCase.verifyError(@() prepareDualCameraCalibrationData(rawFolder, ...
                'InfoReader', @(~) struct('MultiCam', 1), ...
                'ClassificationFcn', @failClassification, ...
                'TempRoot', testCase.TempRoot), 'Test:ClassificationFailed');
            testCase.verifyEmpty(testCase.childFolders());

            function failClassification(varargin)
                error('Test:ClassificationFailed', 'Classification failed.');
            end
        end

        function testDefaultImplementationUsesRequiredInterfaces(testCase)
            % classifyUnregisteredData deliberately calls the legacy,
            % Rig-independent ImagesClassification (not run_ImagesClassification)
            % so calibration source data is never coregistered against an
            % active transform, even when the default Rig has one. See
            % commit a14c260 ("(fix) minor fixes with analysis functions and
            % PipelineManager.").
            source = fileread(which('prepareDualCameraCalibrationData'));
            testCase.verifyNotEmpty(regexp(source, ...
                'InfoReader[^\n]*=\s*@ReadInfoFile', 'once'));
            testCase.verifyNotEmpty(regexp(source, ...
                'ClassificationFcn[^\n]*=\s*@classifyUnregisteredData', 'once'));
            testCase.verifyNotEmpty(regexp(source, ...
                'classifiedFiles\s*=\s*ImagesClassification\(rawFolder,\s*saveFolder,\s*1,\s*1,\s*0\)', ...
                'once'));
        end
    end

    methods (Access = private)
        function rawFolder = createRawFolder(testCase)
            rawFolder = fullfile(testCase.TempRoot, 'raw');
            mkdir(rawFolder);
        end

        function folders = childFolders(testCase)
            entries = dir(testCase.TempRoot);
            entries = entries([entries.isdir]);
            names = string({entries.name});
            folders = names(~ismember(names, [".", "..", "raw"]));
        end
    end

    methods (Static, Access = private)
        function removeFolder(folder)
            if isfolder(folder)
                rmdir(folder, 's');
            end
        end
    end
end
