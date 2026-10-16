classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestOfferReferenceAdoption < matlab.unittest.TestCase
    %TESTOFFERREFERENCEADOPTION DataViewer offer to adopt a folder's reference.
    %
    %   offerReferenceAdoption runs planReferenceAdoption and, for 'canAdd',
    %   adoptFolderReference with the user's choice. Run non-interactively;
    %   the dialogs are checked by hand.

    properties
        Root = ''
        Folder = ''
        Store = []
    end

    methods (TestMethodSetup)
        function alignWithoutProject(testCase)
            testCase.Root = tempname;
            mkdir(testCase.Root);
            testCase.addTeardown(@() rmdir(testCase.Root, 's'));
            testCase.Folder = fullfile(testCase.Root, 'target');
            fixture = setupSyntheticRegistrationFolder(testCase.Folder, 'WithStore', false);
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
                'projectName', ['Offer ' suffix(1:8)], 'description', 'test'));
            projectRoot = testCase.Store.ProjectRoot;
            testCase.addTeardown(@() iRemoveFolder(projectRoot));
            testCase.Store.addSubject(struct('subjectID', 'Mouse_1'));
        end
    end

    methods (Test)
        function testNotBoundOffersNothing(testCase)
            result = offerReferenceAdoption([], testCase.Folder, 'Quiet', true, ...
                'Interactive', false, 'Choice', 'add');
            testCase.verifyEqual(result.state, 'notBound');
            testCase.verifyFalse(result.changed);
            testCase.verifyEmpty(testCase.Store.listSubjectResources('Mouse_1', 'Status', {}));
        end

        function testSkipAddsNothing(testCase)
            testCase.bind();
            result = offerReferenceAdoption([], testCase.Folder, 'Interactive', false);
            testCase.verifyEqual(result.state, 'canAdd');
            testCase.verifyFalse(result.changed);
            testCase.verifySubstring(result.message, 'Add Folder Reference to Library');
            testCase.verifyEmpty(testCase.Store.listSubjectResources('Mouse_1', 'Status', {}));
        end

        function testAddAndMakeActiveThenShownOnly(testCase)
            testCase.bind();
            paramsBefore = computeFileChecksum(fullfile(testCase.Folder, 'DataParams.mat'));

            result = offerReferenceAdoption([], testCase.Folder, 'Interactive', false, ...
                'Choice', 'addActive');
            testCase.verifyTrue(result.changed);
            testCase.verifyEqual(result.severity, 'success');
            testCase.verifySubstring(result.message, 'made active');
            active = testCase.Store.getActiveImageReference('Mouse_1');
            testCase.verifyEqual(active.uuid, result.resourceUUID);
            testCase.verifyEqual(computeFileChecksum(fullfile(testCase.Folder, ...
                'DataParams.mat')), paramsBefore, 'the folder is never changed');

            again = offerReferenceAdoption([], testCase.Folder, 'Interactive', false, ...
                'Choice', 'add');
            testCase.verifyEqual(again.state, 'inLibrary');
            testCase.verifyFalse(again.changed, 'a match is shown only');
            testCase.verifyNumElements(testCase.Store.listSubjectResources('Mouse_1', ...
                'Status', {}), 1);
        end

        function testAddWithoutMakingActive(testCase)
            testCase.bind();
            result = offerReferenceAdoption([], testCase.Folder, 'Interactive', false, ...
                'Choice', 'add');
            testCase.verifyTrue(result.changed);
            testCase.verifyEqual(numel(testCase.Store.listSubjectResources('Mouse_1', ...
                'Status', {})), 1);
            testCase.verifyEqual(planReferenceAdoption(testCase.Folder).state, 'inLibrary');
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
