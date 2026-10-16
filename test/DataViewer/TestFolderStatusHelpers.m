classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestFolderStatusHelpers < matlab.unittest.TestCase
    %TESTFOLDERSTATUSHELPERS DataViewer status back end (project redesign Phase 5).
    %
    %   describeFolderRegistration (alignment status from
    %   DataParams.registration), describeFolderProjectLedger and
    %   formatProjectLedgerHistory. All three are read-only. The DataViewer
    %   labels and dialogs that show them are an app edit (deferred).

    properties
        Folder = ''
    end

    methods (TestMethodSetup)
        function createFolder(testCase)
            addpath(genpath(fileparts(fileparts(fileparts(mfilename('fullpath'))))));
            testCase.Folder = tempname;
            mkdir(testCase.Folder);
            testCase.addTeardown(@() rmdir(testCase.Folder, 's'));
        end
    end

    methods (Test)
        function testRegistrationStates(testCase)
            summary = describeFolderRegistration(testCase.Folder);
            testCase.verifyEqual(summary.state, 'NoDataParams');
            testCase.verifyFalse(summary.isAligned);
            testCase.verifyFalse(isfile(fullfile(testCase.Folder, 'DataParams.mat')), ...
                'Describing never writes.');
            testCase.verifyTrue(any(contains(summary.tooltipLines, 'reused across recordings')), ...
                'project-less note');

            DataParams = createDataParams(testCase.Folder);
            saveDataParams(testCase.Folder, DataParams);
            testCase.verifyEqual(describeFolderRegistration(testCase.Folder).state, 'NoReference');

            [reference, estimate] = testCase.referenceAndEstimate(magic(6), 'Ref A');
            DataParams.registration.reference = reference;
            saveDataParams(testCase.Folder, DataParams);
            summary = describeFolderRegistration(testCase.Folder);
            testCase.verifyEqual(summary.state, 'ReferenceSelected');
            testCase.verifyEqual(summary.referenceName, 'Ref A');
            testCase.verifyEqual(summary.sourceKind, 'file');
            testCase.verifyEqual(summary.shortHash, UMITLedgerCore.shortHash(reference.contentHash));
            testCase.verifyEqual(summary.sourcePath, 'C:\refs\a.tif');

            DataParams.registration.estimate = estimate;
            saveDataParams(testCase.Folder, DataParams);
            testCase.verifyEqual(describeFolderRegistration(testCase.Folder).state, 'Estimated');

            otherReference = testCase.referenceAndEstimate(magic(6)', 'Ref B');
            DataParams.registration.reference = otherReference;
            saveDataParams(testCase.Folder, DataParams);
            summary = describeFolderRegistration(testCase.Folder);
            testCase.verifyEqual(summary.state, 'EstimateStale');
            testCase.verifyTrue(any(contains(summary.tooltipLines, 'Estimate again')));

            DataParams.registration.reference = reference;
            DataParams.view.imageSizeYX = [6 6];
            DataParams.mask.logical = true(6);
            DataParams.registration.apply.isRegistered = true;
            DataParams.registration.apply.appliedReferenceHash = reference.contentHash;
            DataParams.registration.apply.appliedOn = '2026-10-15 10:00:00';
            DataParams.registration.apply.appliedBy = 'test';
            DataParams.registration.apply.confirmationMode = 'pipeline';
            saveDataParams(testCase.Folder, DataParams);
            summary = describeFolderRegistration(testCase.Folder);
            testCase.verifyEqual(summary.state, 'Aligned');
            testCase.verifyTrue(summary.isAligned);
            testCase.verifySubstring(summary.displayText, 'aligned to "Ref A"');
            testCase.verifyEqual(summary.appliedBy, 'test');
            testCase.verifyTrue(any(contains(summary.tooltipLines, 'fresh import')));
        end

        function testRegistrationShowsTheAlignedSourceWarning(testCase)
            DataParams = createDataParams(testCase.Folder);
            reference = testCase.referenceAndEstimate(magic(5), 'From aligned folder');
            reference.source = referenceSourceDescriptor('umitDat', 'D:\other\green.dat');
            reference.source.sourceAligned = true;
            reference.source.warning = 'The source data are already transformed.';
            DataParams.registration.reference = reference;
            saveDataParams(testCase.Folder, DataParams);
            summary = describeFolderRegistration(testCase.Folder);
            testCase.verifyEqual(summary.alignedSourceWarning, reference.source.warning);
            testCase.verifyTrue(any(contains(summary.tooltipLines, 'already transformed')));
        end

        function testInvalidDataParamsIsAStateNotAnError(testCase)
            DataParams = 42;
            save(fullfile(testCase.Folder, 'DataParams.mat'), 'DataParams');
            summary = describeFolderRegistration(testCase.Folder);
            testCase.verifyEqual(summary.state, 'Invalid');
            testCase.verifyNotEmpty(summary.issues);
        end

        function testBoundFolderHasNoProjectNote(testCase)
            UMITProjectLedger.bind(testCase.Folder, testCase.binding(), ...
                testCase.snapshot('P', 'S1', 'Sess1'), 'Origin', 'test');
            summary = describeFolderRegistration(testCase.Folder);
            testCase.verifyTrue(summary.isProjectBound);
            testCase.verifyFalse(any(contains(summary.tooltipLines, 'reused across recordings')));
        end

        function testProjectLedgerStates(testCase)
            summary = describeFolderProjectLedger(testCase.Folder);
            testCase.verifyEqual(summary.state, 'Unknown');
            testCase.verifyFalse(summary.canViewHistory);
            testCase.verifyFalse(UMITProjectLedger.exists(testCase.Folder), 'Describing never writes.');

            binding = testCase.binding();
            UMITProjectLedger.bind(testCase.Folder, binding, ...
                testCase.snapshot('Cortex', 'M01', 'Day1'), 'Origin', 'test');
            summary = describeFolderProjectLedger(testCase.Folder);
            testCase.verifyEqual(summary.state, 'Bound');
            testCase.verifyEqual(summary.binding.sessionUUID, binding.sessionUUID);
            testCase.verifyEqual(summary.displayText, 'Cortex / M01 / Day1');
            testCase.verifyTrue(summary.canViewHistory);

            UMITProjectLedger.unbind(testCase.Folder, 'Reason', 'folderUnbound', 'Origin', 'test');
            summary = describeFolderProjectLedger(testCase.Folder);
            testCase.verifyEqual(summary.state, 'Unbound');
            testCase.verifyEmpty(summary.binding);

            fid = fopen(UMITProjectLedger.filePath(testCase.Folder), 'w');
            fprintf(fid, '{"format": "something else"}');
            fclose(fid);
            summary = describeFolderProjectLedger(testCase.Folder);
            testCase.verifyEqual(summary.state, 'Invalid');
            testCase.verifyFalse(summary.isValid);
        end

        function testSnapshotContentIsDisplayOnly(testCase)
            % Any content is tolerated in the snapshot, including none.
            UMITProjectLedger.bind(testCase.Folder, testCase.binding(), ...
                struct('project', 'not a struct', 'subjectID', {{'a'; 'b'}}), 'Origin', 'test');
            summary = describeFolderProjectLedger(testCase.Folder);
            testCase.verifyEqual(summary.state, 'Bound');
            testCase.verifyEqual(summary.displayText, 'bound');
            lines = formatProjectLedgerHistory(testCase.Folder);
            testCase.verifyTrue(any(contains(lines, 'bound to project')));
        end

        function testProjectLedgerHistory(testCase)
            lines = formatProjectLedgerHistory(testCase.Folder);
            testCase.verifyNumElements(lines, 1);
            testCase.verifySubstring(lines{1}, 'never bound');

            first = testCase.binding();
            second = testCase.binding();
            UMITProjectLedger.bind(testCase.Folder, first, ...
                testCase.snapshot('Cortex', 'M01', 'Day1'), 'Origin', 'test');
            UMITProjectLedger.recordFreshImport(testCase.Folder, 'Origin', 'test');
            UMITProjectLedger.rebind(testCase.Folder, second, ...
                testCase.snapshot('Cortex', 'M02', 'Day1'), 'Origin', 'test');
            UMITProjectLedger.recordOrphanRemoval(testCase.Folder, 'Origin', 'test');
            UMITProjectLedger.bind(testCase.Folder, first, ...
                testCase.snapshot('Cortex', 'M01', 'Day1'), 'Origin', 'test');
            UMITProjectLedger.unbind(testCase.Folder, 'Reason', 'relocated', 'Origin', 'test');

            lines = formatProjectLedgerHistory(testCase.Folder);
            testCase.verifySubstring(lines{2}, 'Unbound');
            body = strjoin(lines, newline);
            testCase.verifySubstring(body, 'bound to project "Cortex", subject "M01", session "Day1"');
            testCase.verifySubstring(body, 'fresh import of the raw data');
            testCase.verifySubstring(body, 'rebound to project "Cortex", subject "M02"');
            testCase.verifySubstring(body, 'its project was not found');
            testCase.verifySubstring(body, 'the session was moved to another folder');
            testCase.verifyNumElements(lines, 4 + 6, 'header lines + one line per event');
        end
    end

    methods (Access = private)
        function [reference, estimate] = referenceAndEstimate(~, image, name)
            [~, reference, estimate] = newRegistrationRecord();
            [reference.contentHash, reference.image] = computeImageContentHash(single(image));
            reference.imageSizeYX = size(image);
            reference.name = name;
            reference.selectedOn = '2026-10-15 09:00:00';
            reference.source = referenceSourceDescriptor('file', 'C:\refs\a.tif');
            estimate.tform = affine2d(eye(3));
            estimate.referenceHash = reference.contentHash;
        end

        function binding = binding(~)
            binding = struct( ...
                'bindingUUID', char(java.util.UUID.randomUUID()), ...
                'projectUUID', char(java.util.UUID.randomUUID()), ...
                'subjectUUID', char(java.util.UUID.randomUUID()), ...
                'sessionUUID', char(java.util.UUID.randomUUID()), ...
                'folderRole', 'processedDataFolder');
        end

        function snapshot = snapshot(~, projectName, subjectID, sessionID)
            snapshot = UMITProjectLedger.makeSnapshot( ...
                struct('projectName', projectName, 'experimenter', 'B. S.'), subjectID, sessionID);
        end
    end
end
