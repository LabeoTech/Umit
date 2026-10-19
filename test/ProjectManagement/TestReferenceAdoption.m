classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestReferenceAdoption < matlab.unittest.TestCase
    %TESTREFERENCEADOPTION Adoption of a folder aligned without a project (Phase 7).
    %
    %   planReferenceAdoption / adoptFolderReference (section 4 Q7): a match
    %   in the subject library is shown only; otherwise the committed
    %   reference can be added. DataParams.registration is never modified.

    properties
        Root = ''
        Folder = ''
        ReferenceImage = []
        Store = []
    end

    methods (TestMethodSetup)
        function alignWithoutProject(testCase)
            testCase.Root = tempname;
            mkdir(testCase.Root);
            testCase.addTeardown(@() rmdir(testCase.Root, 's'));
            testCase.Folder = fullfile(testCase.Root, 'target');
            fixture = setupSyntheticRegistrationFolder(testCase.Folder, 'WithStore', false);
            testCase.ReferenceImage = fixture.ReferenceImage;
            source = fullfile(testCase.Root, 'source');
            mkdir(source);
            writeTestDat(fullfile(source, 'green.dat'), ...
                repmat(fixture.ReferenceImage, 1, 1, 3), 10, 20);
            AcqInfoStream = struct('Width', fixture.Nx, 'Height', fixture.Ny);
            save(fullfile(source, 'AcqInfos.mat'), 'AcqInfoStream');
            createRegistrationTform(testCase.Folder, 'ReferenceSource', ...
                fullfile(source, 'green.dat'), 'ShowFigure', false);
            applyRegistrationTformOnFolder(testCase.Folder, 'RequireUserConfirmation', false, ...
                'OpenQCFigure', false);

            suffix = char(java.util.UUID.randomUUID());
            testCase.Store = UMITProjectStore.create(struct( ...
                'projectName', ['Adoption ' suffix(1:8)], 'description', 'test'));
            projectRoot = testCase.Store.ProjectRoot;
            testCase.addTeardown(@() iRemoveFolder(projectRoot));
            testCase.Store.addSubject(struct('subjectID', 'Mouse_1'));
        end
    end

    methods (Test)
        function testNotBoundAndNotAligned(testCase)
            testCase.verifyEqual(planReferenceAdoption(testCase.Folder).state, 'notBound');

            unaligned = fullfile(testCase.Root, 'unaligned');
            setupSyntheticRegistrationFolder(unaligned, 'WithStore', false);
            testCase.Store.addSession('Mouse_1', struct('sessionID', 'S0', ...
                'processedDataFolder', unaligned));
            testCase.verifyEqual(planReferenceAdoption(unaligned).state, 'notAligned');
        end

        function testBoundFolderWithoutDataParamsIsNotAligned(testCase)
            % DFR-20261016-001: no DataParams.mat is an unaligned folder,
            % not an error, and nothing is written.
            bare = fullfile(testCase.Root, 'noDataParams');
            setupSyntheticRegistrationFolder(bare, 'WithStore', false);
            delete(fullfile(bare, 'DataParams.mat'));
            testCase.Store.addSession('Mouse_1', struct('sessionID', 'S1', ...
                'processedDataFolder', bare));
            before = {dir(bare).name};

            plan = planReferenceAdoption(bare);

            testCase.verifyEqual(plan.state, 'notAligned');
            testCase.verifyEqual(plan.subjectID, 'Mouse_1');
            testCase.verifyFalse(isfile(fullfile(bare, 'DataParams.mat')));
            testCase.verifyEqual({dir(bare).name}, before);
        end

        function testAddTheCommittedReferenceOnce(testCase)
            testCase.bind();
            before = computeFileChecksum(fullfile(testCase.Folder, 'DataParams.mat'));
            committed = loadDataParams(testCase.Folder).registration.apply.appliedReferenceHash;

            plan = planReferenceAdoption(testCase.Folder);
            testCase.verifyEqual(plan.state, 'canAdd');
            testCase.verifyEqual(plan.committedHash, committed);
            testCase.verifyEqual(plan.sourceKind, 'umitDat');

            resourceUUID = adoptFolderReference(testCase.Folder, 'Name', 'Adopted', ...
                'SetActive', true);
            testCase.verifyEqual(computeFileChecksum(fullfile(testCase.Folder, 'DataParams.mat')), ...
                before, 'DataParams.registration is never modified.');
            active = testCase.Store.getActiveImageReference('Mouse_1');
            testCase.verifyEqual(active.uuid, resourceUUID);
            loaded = load(active.absolutePath, 'ImageReference');
            testCase.verifyEqual(loaded.ImageReference.contentHash, committed);
            testCase.verifyEqual(loaded.ImageReference.subjectID, 'Mouse_1');

            plan = planReferenceAdoption(testCase.Folder);
            testCase.verifyEqual(plan.state, 'inLibrary');
            testCase.verifyEqual(plan.match.uuid, resourceUUID);
            testCase.verifyError(@() adoptFolderReference(testCase.Folder), ...
                'Umitoolbox:ProjectManagement:cannotAdoptReference');
        end

        function testExistingLibraryMatchIsShownOnly(testCase)
            ProjectInfo = testCase.Store.getProjectInfo();
            SubjectInfo = testCase.Store.getSubjectInfo('Mouse_1');
            ImageReference = genImageReferenceStruct(testCase.ReferenceImage, ...
                'Name', 'Library copy', 'ProjectUUID', ProjectInfo.projectUUID, ...
                'ProjectName', ProjectInfo.projectName, 'SubjectUUID', SubjectInfo.uuid, ...
                'SubjectID', 'Mouse_1');
            staged = [tempname '.mat'];
            save(staged, 'ImageReference');
            testCase.addTeardown(@() delete(staged));
            libraryUUID = testCase.Store.addImageReference('Mouse_1', staged, ...
                struct('displayName', 'Library copy'));
            % Archive the matching reference behind another active one.
            ImageReference = genImageReferenceStruct(rand(8, 'single'), ...
                'Name', 'Other', 'ProjectUUID', ProjectInfo.projectUUID, ...
                'ProjectName', ProjectInfo.projectName, 'SubjectUUID', SubjectInfo.uuid, ...
                'SubjectID', 'Mouse_1');
            save(staged, 'ImageReference');
            otherUUID = testCase.Store.addImageReference('Mouse_1', staged, struct());
            testCase.Store.setActiveImageReference('Mouse_1', otherUUID);
            testCase.Store.archiveResource(libraryUUID);

            testCase.bind();
            ledgerBefore = fileread(UMITProjectLedger.filePath(testCase.Folder));
            paramsBefore = computeFileChecksum(fullfile(testCase.Folder, 'DataParams.mat'));
            resourcesBefore = testCase.Store.listSubjectResources('Mouse_1', 'Status', {});

            plan = planReferenceAdoption(testCase.Folder);
            testCase.verifyEqual(plan.state, 'inLibrary');
            testCase.verifyEqual(plan.match.uuid, libraryUUID);
            testCase.verifyEqual(plan.match.status, 'archived');
            testCase.verifySubstring(plan.message, 'Library copy');

            testCase.verifyEqual(fileread(UMITProjectLedger.filePath(testCase.Folder)), ledgerBefore);
            testCase.verifyEqual(computeFileChecksum(fullfile(testCase.Folder, 'DataParams.mat')), ...
                paramsBefore);
            testCase.verifyEqual(numel(testCase.Store.listSubjectResources('Mouse_1', 'Status', {})), ...
                numel(resourcesBefore), 'nothing is recorded for a match');
        end
    end

    methods (Access = private)
        function bind(testCase)
            testCase.Store.addSession('Mouse_1', struct('sessionID', 'S1', ...
                'processedDataFolder', testCase.Folder));
        end
    end
end

function iRemoveFolder(folder)
if isfolder(folder)
    rmdir(folder, 's');
end
end
