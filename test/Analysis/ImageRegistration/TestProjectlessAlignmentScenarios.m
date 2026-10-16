classdef (SharedTestFixtures = {UMITTemporaryRootFixture('noProjectStore')}) TestProjectlessAlignmentScenarios < matlab.unittest.TestCase
    %TESTPROJECTLESSALIGNMENTSCENARIOS Scenarios 2 and 2b with no project store.
    %
    %   Project redesign Phase 5: alignment works on folders that are not
    %   bound to a project, on a machine without a project store. The class
    %   has a UMIT root of its own (tag 'noProjectStore'); every test checks
    %   that the project store was never even resolved (no 'projects' folder
    %   appears below the root) and that the .dat source folder is left
    %   byte-identical. Scenario 1 (project-sourced references) is covered by
    %   TestCreateRegistrationTform and TestApplyRegistrationTformOnFolder.

    properties
        Root = ''
        Target = ''
        SourceFolder = ''
        SourceDat = ''
        SourceHash = ''
    end

    methods (TestMethodSetup)
        function createFolders(testCase)
            testCase.Root = tempname;
            mkdir(testCase.Root);
            testCase.addTeardown(@() rmdir(testCase.Root, 's'));
            testCase.Target = fullfile(testCase.Root, 'target');
            fixture = setupSyntheticRegistrationFolder(testCase.Target, 'WithStore', false);

            % A UMIT save folder whose green.dat averages to the reference.
            testCase.SourceFolder = fullfile(testCase.Root, 'source');
            mkdir(testCase.SourceFolder);
            writeTestDat(fullfile(testCase.SourceFolder, 'green.dat'), ...
                repmat(fixture.ReferenceImage, 1, 1, 4), 10, 20);
            AcqInfoStream = struct('Width', fixture.Nx, 'Height', fixture.Ny);
            save(fullfile(testCase.SourceFolder, 'AcqInfos.mat'), 'AcqInfoStream');
            testCase.SourceDat = fullfile(testCase.SourceFolder, 'green.dat');
            testCase.SourceHash = computeImageContentHash(fixture.ReferenceImage);
        end
    end

    methods (Test)
        function testScenario2OneOffDatSource(testCase)
            before = testCase.folderChecksums(testCase.SourceFolder);

            createRegistrationTform(testCase.Target, 'ReferenceSource', testCase.SourceDat, ...
                'ShowFigure', false);
            applyRegistrationTformOnFolder(testCase.Target, 'RequireUserConfirmation', false, ...
                'OpenQCFigure', false);

            registration = loadDataParams(testCase.Target).registration;
            testCase.verifyTrue(registration.apply.isRegistered);
            testCase.verifyEqual(registration.apply.appliedReferenceHash, testCase.SourceHash);
            testCase.verifyEqual(registration.reference.source.kind, 'umitDat');
            testCase.verifyEqual(registration.reference.source.path, testCase.SourceDat);
            summary = describeFolderRegistration(testCase.Target);
            testCase.verifyEqual(summary.state, 'Aligned');
            testCase.verifyFalse(summary.isProjectBound);

            % Committed: a new estimate, even against the same source, is refused.
            testCase.verifyError(@() createRegistrationTform(testCase.Target, ...
                'ReferenceSource', testCase.SourceDat, 'ShowFigure', false), ...
                'Umitoolbox:ImageRegistration:alreadyAligned');

            testCase.verifyEqual(testCase.folderChecksums(testCase.SourceFolder), before, ...
                'The .dat source folder is never written.');
            testCase.verifyFalse(UMITProjectLedger.exists(testCase.Target));
            testCase.verifyNoProjectStore();
        end

        function testScenario2bBatchFromOneDatSource(testCase)
            before = testCase.folderChecksums(testCase.SourceFolder);
            targets = {testCase.Target, fullfile(testCase.Root, 'target2'), ...
                fullfile(testCase.Root, 'target3')};
            copyfile(testCase.Target, targets{2});
            copyfile(testCase.Target, targets{3});

            for k = 1:numel(targets)
                createRegistrationTform(targets{k}, 'ReferenceSource', testCase.SourceDat, ...
                    'ShowFigure', false);
                applyRegistrationTformOnFolder(targets{k}, 'RequireUserConfirmation', false, ...
                    'OpenQCFigure', false);
            end

            for k = 1:numel(targets)
                registration = loadDataParams(targets{k}).registration;
                testCase.verifyTrue(registration.apply.isRegistered, targets{k});
                testCase.verifyEqual(registration.apply.appliedReferenceHash, testCase.SourceHash, ...
                    'One .dat source gives one hash across the batch.');
            end
            testCase.verifyEqual(testCase.folderChecksums(testCase.SourceFolder), before);
            testCase.verifyNoProjectStore();
        end

        function testUnboundFolderWithoutSourceIsRefused(testCase)
            testCase.verifyError(@() createRegistrationTform(testCase.Target, 'ShowFigure', false), ...
                'Umitoolbox:ImageRegistration:noReferenceSource');
            testCase.verifyEmpty(loadDataParams(testCase.Target).registration.estimate);
            testCase.verifyNoProjectStore();
        end
    end

    methods (Access = private)
        function verifyNoProjectStore(testCase)
            projectsRoot = fullfile(getenv('UMIT_ROOT'), 'projects');
            testCase.verifyFalse(isfolder(projectsRoot), ...
                'A project-less alignment must not resolve the project store.');
        end

        function sums = folderChecksums(~, folder)
            listing = dir(fullfile(folder, '*'));
            listing = listing(~[listing.isdir]);
            sums = cell2struct(cellfun(@(n) computeFileChecksum(fullfile(folder, n)), ...
                {listing.name}, 'UniformOutput', false), ...
                matlab.lang.makeValidName({listing.name}), 2);
        end
    end
end
