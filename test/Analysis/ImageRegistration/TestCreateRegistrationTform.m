classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestCreateRegistrationTform < matlab.unittest.TestCase
    %TESTCREATEREGISTRATIONTFORM Unit tests for createRegistrationTform.

    properties
        TempFolder
        Fixture
    end

    methods (TestMethodSetup)
        function setup(testCase)
            testCase.TempFolder = fullfile(tempdir, ['TestCreateRegistrationTform_' char(java.util.UUID.randomUUID)]);
            testCase.Fixture = setupSyntheticRegistrationFolder(testCase.TempFolder);
        end
    end

    methods (TestMethodTeardown)
        function teardown(testCase)
            if isstruct(testCase.Fixture) && ...
                    isfield(testCase.Fixture, 'ProjectRoot') && ...
                    isfolder(testCase.Fixture.ProjectRoot)
                rmdir(testCase.Fixture.ProjectRoot, 's');
            end
            if isstruct(testCase.Fixture) && ...
                    isfield(testCase.Fixture, 'RigRoot') && ...
                    isfolder(testCase.Fixture.RigRoot)
                rmdir(testCase.Fixture.RigRoot, 's');
            end
            if isfolder(testCase.TempFolder)
                rmdir(testCase.TempFolder, 's');
            end
        end
    end

    methods (Test)
        function testPipelineManagerExecutesRegistrationNode(testCase)
            pm = buildPMForScenario(testCase.TempFolder, ...
                'createRegistrationTform', 'auto', ...
                'Parameters', struct('ShowFigure', false));
            result = pm.executePipeline('PrintSummary', false);

            testCase.verifyEqual(result.status, "completed");
            testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, 'DataParams.mat')));
            S = load(fullfile(testCase.TempFolder, 'DataParams.mat'), 'DataParams');
            testCase.verifyNotEmpty(S.DataParams.registration.estimate.tform);
            testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, ...
                char(string(S.DataParams.registration.estimate.qcFigureFile)))));
        end

        function testCreatesTformAndQCArtifacts(testCase)
            outFile = createRegistrationTform(testCase.TempFolder, 'ShowFigure', false);

            testCase.verifyTrue(iscell(outFile));
            testCase.verifyEqual(numel(outFile), 3);
            for i = 1:numel(outFile)
                testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, outFile{i})));
            end

            S = load(fullfile(testCase.TempFolder, 'DataParams.mat'), 'DataParams');
            DP = S.DataParams;
            validateDataParams(DP);
            registration = DP.registration;
            testCase.verifyEqual(registration.version, 2);
            testCase.verifyFalse(registration.apply.isRegistered);
            testCase.verifyFalse(registration.apply.isReviewed);
            testCase.verifyEmpty(registration.apply.appliedReferenceHash);
            estimate = registration.estimate;
            testCase.verifyNotEmpty(estimate.tform);
            testCase.verifyEqual(char(string(estimate.qcFigureFile)), outFile{2});
            testCase.verifyEqual(char(string(estimate.qcPreviewImageFile)), outFile{3});
            testCase.verifyGreaterThan(estimate.qcMetrics.MIAfter, estimate.qcMetrics.MIBefore);
            testCase.verifyEqual(estimate.movingSource, fullfile(testCase.TempFolder, 'green.dat'));
            reference = registration.reference;
            testCase.verifyEqual(reference.image, testCase.Fixture.ReferenceImage);
            % Brief section 1.3: identified by the content hash of the image.
            testCase.verifyEqual(reference.contentHash, ...
                computeImageContentHash(testCase.Fixture.ReferenceImage));
            testCase.verifyEqual(estimate.referenceHash, reference.contentHash);
            testCase.verifyEqual(reference.description, testCase.Fixture.ManagedDescription);
            source = reference.source;
            testCase.verifyEqual(source.kind, 'project');
            testCase.verifyEqual(source.path, testCase.Fixture.ManagedReference.absolutePath);
            testCase.verifyEqual(source.resourceUUID, testCase.Fixture.ImageReferenceUUID);
            testCase.verifyEqual(source.projectUUID, testCase.Fixture.ProjectUUID);
            testCase.verifyEqual(source.subjectID, testCase.Fixture.SubjectID);
            testCase.verifyFalse(any(isfield(registration, ...
                {'referenceFile', 'imageReferenceUUID', 'imageReferenceChecksum', 'resourceUUID'})));
        end

        function testSpecificUseFile(testCase)
            outFile = createRegistrationTform(testCase.TempFolder, 'UseFile', 'red.dat', 'ShowFigure', false);
            testCase.verifyTrue(iscell(outFile));
            S = load(fullfile(testCase.TempFolder, 'DataParams.mat'), 'DataParams');
            testCase.verifyEqual( ...
                S.DataParams.registration.reference.source.path, ...
                testCase.Fixture.ManagedReference.absolutePath);
            testCase.verifyEqual( ...
                S.DataParams.registration.estimate.movingSource, ...
                fullfile(testCase.TempFolder, 'red.dat'));
        end

        function testOtherLayoutsAreRefused(testCase)
            layouts = {{'Y','X','T','E'}, [6 5 4 2]; {'Y','X','E'}, [6 5 3]; {'Y','X'}, [6 5]};
            for k = 1:size(layouts, 1)
                f = fullfile(testCase.TempFolder, 'layout.dat');
                writeTestDat(f, rand([layouts{k, 2}, 1], 'single'), 10, 5, ...
                    'DimNames', layouts{k, 1});
                testCase.verifyError(@() createRegistrationTform(testCase.TempFolder, ...
                    'UseFile', 'layout.dat', 'ShowFigure', false), ...
                    'Umitoolbox:createRegistrationTform:unsupportedLayout');
            end
        end

        function testAutoUsesReferenceSourceInsteadOfDerivedDat(testCase)
            derivedData = testCase.Fixture.RedData;
            testCase.writeSingleDat( ...
                fullfile(testCase.TempFolder, 'GSR.dat'), derivedData, 10, 20);

            createRegistrationTform( ...
                testCase.TempFolder, 'ShowFigure', false);

            DataParams = loadDataParams(testCase.TempFolder);
            testCase.verifyEqual( ...
                DataParams.registration.estimate.movingSource, ...
                fullfile(testCase.TempFolder, 'green.dat'));
        end

        function testAutoWithoutUsableSourceNeverGuessesBetweenChannels(testCase)
            % Audit item 20 / Bruno 2026-10-15: without a usable .dat source,
            % 'auto' takes the only Y-X-T .dat; several are an error.
            for source = {'', 'C:\elsewhere\reference.tif'}
                testCase.updateActiveReferenceSourceFile(source{1});
                testCase.verifyPreflightFailure( ...
                    @() createRegistrationTform( ...
                    testCase.TempFolder, 'ShowFigure', false), ...
                    'Umitoolbox:createRegistrationTform:AutoSourceFileAmbiguous');
            end
        end

        function testAutoWithTifMatOrUmtReferenceUsesTheOnlyYXTDat(testCase)
            % References built from a TIF, MAT or UMT file (audit item 20).
            delete(fullfile(testCase.TempFolder, 'red.dat'));
            writeTestDat(fullfile(testCase.TempFolder, 'map.dat'), ...
                rand(testCase.Fixture.Ny, testCase.Fixture.Nx, 'single'), 10, ...
                'DimNames', {'Y','X'});
            for source = {'reference.tif', 'reference.mat', 'reference.umt'}
                testCase.updateActiveReferenceSourceFile( ...
                    fullfile(testCase.TempFolder, source{1}));

                createRegistrationTform(testCase.TempFolder, 'ShowFigure', false);

                DataParams = loadDataParams(testCase.TempFolder);
                testCase.verifyEqual(DataParams.registration.estimate.movingSource, ...
                    fullfile(testCase.TempFolder, 'green.dat'), source{1});
            end
            delete(fullfile(testCase.TempFolder, 'green.dat'));
            testCase.verifyPreflightFailure( ...
                @() createRegistrationTform(testCase.TempFolder, 'ShowFigure', false), ...
                'Umitoolbox:createRegistrationTform:AutoSourceFileNotFound');
        end

        function testEstimateIsRefusedOnAnAlignedFolder(testCase)
            % Brief section 4 Q1: estimating is refused after apply, and the
            % apply record is never reset (audit item 17).
            createRegistrationTform(testCase.TempFolder, 'ShowFigure', false);
            applyRegistrationTformOnFolder(testCase.TempFolder, ...
                'RequireUserConfirmation', false, 'OpenQCFigure', false);

            testCase.verifyPreflightFailure( ...
                @() createRegistrationTform(testCase.TempFolder, 'ShowFigure', false), ...
                'Umitoolbox:ImageRegistration:alreadyAligned');
            DataParams = loadDataParams(testCase.TempFolder);
            testCase.verifyTrue(DataParams.registration.apply.isRegistered);
            testCase.verifyNotEmpty(DataParams.registration.apply.appliedOn);
        end

        function testDifferentSizeMovingImageIsComposedIntoTransform(testCase)
            % The stored transform maps the folder's own frame into the
            % reference: a 2x larger recording aligns onto the 64x64 reference.
            fx = testCase.Fixture;
            for name = {'green', 'red'}
                data = single(loadData(fullfile(testCase.TempFolder, [name{1} '.dat'])));
                writeTestDat(fullfile(testCase.TempFolder, [name{1} '.dat']), ...
                    imresize(data, 2), 10, 20);
            end
            DataParams = loadDataParams(testCase.TempFolder);
            DataParams.view.imageSizeYX = [];
            DataParams.mask.logical = [];
            saveDataParams(testCase.TempFolder, DataParams);

            createRegistrationTform(testCase.TempFolder, 'ShowFigure', false);
            applyRegistrationTformOnFolder(testCase.TempFolder, ...
                'RequireUserConfirmation', false, 'OpenQCFigure', false);

            aligned = single(loadData(fullfile(testCase.TempFolder, 'green.dat')));
            testCase.verifySize(aligned, [fx.Ny fx.Nx fx.Nt]);
            center = 12:52;
            residual = aligned(center, center, 1) - fx.ReferenceImage(center, center);
            testCase.verifyLessThan(max(abs(residual(:))), 0.1);
        end

        function testDifferentAspectRatioIsRefused(testCase)
            data = single(loadData(fullfile(testCase.TempFolder, 'green.dat')));
            writeTestDat(fullfile(testCase.TempFolder, 'green.dat'), ...
                imresize(data, [64 96]), 10, 20);
            DataParams = loadDataParams(testCase.TempFolder);
            DataParams.view.imageSizeYX = [];
            DataParams.mask.logical = [];
            saveDataParams(testCase.TempFolder, DataParams);

            testCase.verifyPreflightFailure( ...
                @() createRegistrationTform(testCase.TempFolder, 'UseFile', 'green.dat', ...
                'ShowFigure', false), ...
                'Umitoolbox:createRegistrationTform:AnisotropicResize');
        end

        function testNewActiveReferenceSupersedesFolderSnapshot(testCase)
            createRegistrationTform(testCase.TempFolder, 'ShowFigure', false);
            firstHash = loadDataParams(testCase.TempFolder).registration.reference.contentHash;

            secondImage = 1.5 .* testCase.Fixture.ReferenceImage;
            [secondUUID, secondResource] = ...
                testCase.addAndActivateReference( ...
                secondImage, 'Second active reference');

            createRegistrationTform( ...
                testCase.TempFolder, 'ShowFigure', false);

            % Before apply the reference is free: the new active reference
            % replaces the earlier one (brief section 1.5).
            DataParams = loadDataParams(testCase.TempFolder);
            testCase.verifyEqual( ...
                DataParams.registration.reference.image, single(secondImage));
            testCase.verifyEqual(DataParams.registration.reference.contentHash, ...
                computeImageContentHash(secondImage));
            testCase.verifyNotEqual(DataParams.registration.reference.contentHash, firstHash);
            testCase.verifyEqual(DataParams.registration.estimate.referenceHash, ...
                DataParams.registration.reference.contentHash);
            testCase.verifyEqual( ...
                DataParams.registration.reference.source.resourceUUID, secondUUID);
            testCase.verifyEqual( ...
                DataParams.registration.reference.source.path, ...
                secondResource.absolutePath);
        end

        function testNoActiveReferenceDoesNotUseAvailableResource(testCase)
            testCase.clearActiveReferencePointer();

            testCase.verifyPreflightFailure( ...
                @() createRegistrationTform( ...
                testCase.TempFolder, 'ShowFigure', false), ...
                'Umitoolbox:ImageRegistration:referenceNotFound');
        end

        function testMissingManagedReferenceFileDoesNotMutateFolder(testCase)
            delete(testCase.Fixture.ManagedReference.absolutePath);

            testCase.verifyPreflightFailure( ...
                @() createRegistrationTform( ...
                testCase.TempFolder, 'ShowFigure', false), ...
                'Umitoolbox:ImageRegistration:referenceNotFound');
        end

        function testInvalidManagedReferencePayloadDoesNotMutateFolder(testCase)
            managedPath = testCase.Fixture.ManagedReference.absolutePath;
            loaded = load(managedPath, 'ImageReference', '-mat');
            ImageReference = rmfield(loaded.ImageReference, 'image');
            save(managedPath, 'ImageReference', '-mat');

            SubjectInfo = testCase.Fixture.Store.getSubjectInfo( ...
                testCase.Fixture.SubjectID);
            idx = find(strcmp( ...
                {SubjectInfo.resourceRegistry.uuid}, ...
                testCase.Fixture.ImageReferenceUUID), 1, 'first');
            SubjectInfo.resourceRegistry(idx).checksum = ...
                computeFileChecksum(managedPath);
            SubjectInfo.resourceRegistry(idx).modifiedOn = datetime('now');
            SubjectInfo.modifiedOn = datetime('now');
            testCase.saveSubjectInfo(SubjectInfo);

            testCase.verifyPreflightFailure( ...
                @() createRegistrationTform( ...
                testCase.TempFolder, 'ShowFigure', false), ...
                'Umitoolbox:ImageRegistration:invalidSource');
        end

        function testMissingBindingDoesNotMutateFolder(testCase)
            % Without a binding and without ReferenceSource the folder has
            % no reference (Phase 4 resolution order); nothing is written.
            bindingPath = fullfile( ...
                testCase.TempFolder, ...
                UMITProjectStore.getProjectBindingFileName());
            delete(bindingPath);

            testCase.verifyPreflightFailure( ...
                @() createRegistrationTform( ...
                testCase.TempFolder, 'ShowFigure', false), ...
                'Umitoolbox:ImageRegistration:noReferenceSource');
        end

        function testMalformedBindingDoesNotMutateFolder(testCase)
            bindingPath = fullfile( ...
                testCase.TempFolder, ...
                UMITProjectStore.getProjectBindingFileName());
            InvalidBinding = struct('unexpectedField', true);
            save(bindingPath, 'InvalidBinding', '-mat');

            testCase.verifyPreflightFailure( ...
                @() createRegistrationTform( ...
                testCase.TempFolder, 'ShowFigure', false), ...
                ['Umitoolbox:createRegistrationTform:' ...
                 'BindingResolutionFailed']);
        end

        function testInconsistentBindingPathDoesNotMutateFolder(testCase)
            schema = getUMITProjectSchema();
            sessionPath = fullfile( ...
                testCase.Fixture.ProjectRoot, ...
                schema.folders.subjects, ...
                testCase.Fixture.SubjectID, ...
                schema.folders.sessions, ...
                testCase.Fixture.SessionID, ...
                schema.files.sessionMetadata);
            loaded = load(sessionPath, ...
                schema.metadataVariables.session, '-mat');
            SessionInfo = loaded.(schema.metadataVariables.session);
            SessionInfo.processedDataFolder = fullfile( ...
                testCase.TempFolder, 'different_location');
            SessionInfo.modifiedOn = datetime('now');
            saveMatAtomic(sessionPath, ...
                schema.metadataVariables.session, SessionInfo);

            testCase.verifyPreflightFailure( ...
                @() createRegistrationTform( ...
                testCase.TempFolder, 'ShowFigure', false), ...
                ['Umitoolbox:createRegistrationTform:' ...
                 'BindingResolutionFailed']);
        end

        function testEstimateAgainstAThenBThenApplyB(testCase)
            % Brief section 1.5: before apply, selecting and estimating are free.
            [datA, hashA] = testCase.makeDatSource('sourceA', 1);
            [datB, hashB] = testCase.makeDatSource('sourceB', 2);
            testCase.verifyNotEqual(hashA, hashB);

            createRegistrationTform(testCase.TempFolder, 'ReferenceSource', datA, 'ShowFigure', false);
            registration = loadDataParams(testCase.TempFolder).registration;
            testCase.verifyEqual(registration.reference.contentHash, hashA);
            testCase.verifyEqual(registration.reference.source.kind, 'umitDat');
            testCase.verifyEqual(registration.reference.source.path, datA);

            createRegistrationTform(testCase.TempFolder, 'ReferenceSource', datB, 'ShowFigure', false);
            applyRegistrationTformOnFolder(testCase.TempFolder, ...
                'RequireUserConfirmation', false, 'OpenQCFigure', false);

            registration = loadDataParams(testCase.TempFolder).registration;
            testCase.verifyTrue(registration.apply.isRegistered);
            testCase.verifyEqual(registration.apply.appliedReferenceHash, hashB);
        end

        function testApplyAfterTheReferenceChangedIsRefused(testCase)
            [datA, hashA] = testCase.makeDatSource('sourceA', 1);
            [datB, hashB] = testCase.makeDatSource('sourceB', 2);
            createRegistrationTform(testCase.TempFolder, 'ReferenceSource', datA, 'ShowFigure', false);
            selectRegistrationReference(testCase.TempFolder, datB, 'Interactive', false);
            before = testCase.datAndParamsChecksums();

            try
                applyRegistrationTformOnFolder(testCase.TempFolder, ...
                    'RequireUserConfirmation', false, 'OpenQCFigure', false);
                testCase.verifyFail('Expected the apply guard to refuse.');
            catch ME
                testCase.verifyEqual(ME.identifier, ...
                    'Umitoolbox:ImageRegistration:estimateReferenceMismatch');
                testCase.verifySubstring(ME.message, testCase.TempFolder);
                testCase.verifySubstring(ME.message, UMITLedgerCore.shortHash(hashA));
                testCase.verifySubstring(ME.message, UMITLedgerCore.shortHash(hashB));
            end
            testCase.verifyEqual(testCase.datAndParamsChecksums(), before);
        end

        function testDifferentReferenceAfterApplyIsAnError(testCase)
            % Project-sourced and .dat-sourced references follow the same rules.
            [datA, hashA] = testCase.makeDatSource('sourceA', 1);
            datB = testCase.makeDatSource('sourceB', 2);
            createRegistrationTform(testCase.TempFolder, 'ReferenceSource', datA, 'ShowFigure', false);
            applyRegistrationTformOnFolder(testCase.TempFolder, ...
                'RequireUserConfirmation', false, 'OpenQCFigure', false);
            before = testCase.datAndParamsChecksums();

            % The bound project's active reference: refused, both named.
            try
                createRegistrationTform(testCase.TempFolder, 'ShowFigure', false);
                testCase.verifyFail('Expected the committed folder to refuse.');
            catch ME
                testCase.verifyEqual(ME.identifier, 'Umitoolbox:ImageRegistration:alreadyAligned');
                testCase.verifySubstring(ME.message, 'Estimating the image alignment was refused');
                testCase.verifySubstring(ME.message, testCase.TempFolder);
                testCase.verifySubstring(ME.message, 'reference "green.dat" (umitDat, ');
                testCase.verifySubstring(ME.message, UMITLedgerCore.shortHash(hashA));
                testCase.verifySubstring(ME.message, '(project, ');
                testCase.verifySubstring(ME.message, 'fresh import');
            end
            % Selecting: same reference is a no-op, another is refused.
            selectRegistrationReference(testCase.TempFolder, datA, 'Interactive', false);
            testCase.verifyError(@() selectRegistrationReference(testCase.TempFolder, datB, ...
                'Interactive', false), 'Umitoolbox:ImageRegistration:referenceCommitted');
            testCase.verifyError(@() createRegistrationTform(testCase.TempFolder, ...
                'ReferenceSource', datA, 'ShowFigure', false), ...
                'Umitoolbox:ImageRegistration:alreadyAligned');
            testCase.verifyEqual(testCase.datAndParamsChecksums(), before);
        end

        function testBatchFromOneDatSourceSharesOneHash(testCase)
            [datA, hashA] = testCase.makeDatSource('sourceA', 1);
            second = [testCase.TempFolder '_second'];
            copyfile(testCase.TempFolder, second);
            testCase.addTeardown(@() rmdir(second, 's'));

            for folder = {testCase.TempFolder, second}
                createRegistrationTform(folder{1}, 'ReferenceSource', datA, 'ShowFigure', false);
            end

            testCase.verifyEqual(loadDataParams(testCase.TempFolder).registration.reference.contentHash, hashA);
            testCase.verifyEqual(loadDataParams(second).registration.reference.contentHash, hashA);
        end

        function testUnboundFolderNeedsAReferenceSource(testCase)
            % Scenario 2: a one-off .dat source without any project binding.
            unbound = [testCase.TempFolder '_unbound'];
            mkdir(unbound);
            testCase.addTeardown(@() rmdir(unbound, 's'));
            for name = {'green.dat', 'red.dat', 'AcqInfos.mat', 'DataParams.mat'}
                copyfile(fullfile(testCase.TempFolder, name{1}), fullfile(unbound, name{1}));
            end
            testCase.verifyEqual(UMITProjectStore.getProjectBindingState(unbound), 'unbound');

            testCase.verifyError(@() createRegistrationTform(unbound, 'ShowFigure', false), ...
                'Umitoolbox:ImageRegistration:noReferenceSource');

            [datA, hashA] = testCase.makeDatSource('sourceA', 1);
            createRegistrationTform(unbound, 'ReferenceSource', datA, 'ShowFigure', false);
            applyRegistrationTformOnFolder(unbound, 'RequireUserConfirmation', false, ...
                'OpenQCFigure', false);
            registration = loadDataParams(unbound).registration;
            testCase.verifyEqual(registration.apply.appliedReferenceHash, hashA);
            testCase.verifyFalse(UMITProjectLedger.exists(unbound));
        end

        function testVersion1RegistrationIsUpgradedInMemoryOnly(testCase)
            % Existing DataParams.mat files with a version-1 (flat) record
            % load as version 2; the file itself is rewritten only by a writer.
            createRegistrationTform(testCase.TempFolder, 'ShowFigure', false);
            dataParamsPath = fullfile(testCase.TempFolder, 'DataParams.mat');
            v2 = loadDataParams(testCase.TempFolder).registration;
            loaded = load(dataParamsPath, 'DataParams', '-mat');
            DataParams = loaded.DataParams;
            DataParams.registration = struct( ...
                'isRegistered', false, 'isReviewed', false, ...
                'tform', v2.estimate.tform, 'transformType', 'similarity', ...
                'method', 'legacy', 'referenceDescription', 'Legacy ref', ...
                'referenceFile', 'C:\old\ref.mat', ...
                'referenceImage', testCase.Fixture.ReferenceImage, ...
                'imageReferenceUUID', testCase.Fixture.ImageReferenceUUID, ...
                'imageReferenceChecksum', testCase.Fixture.ManagedReference.checksum, ...
                'createdOn', '2026-01-01 10:00:00', 'source', 'green.dat', ...
                'qcStatus', '', 'qcWarning', '', 'qcMetrics', v2.estimate.qcMetrics, ...
                'qcFigureFile', '', 'qcPreviewImageFile', '', 'appliedOn', '', ...
                'appliedBy', '', 'confirmationMode', '', 'notes', '', 'resourceUUID', '');
            save(dataParamsPath, 'DataParams', '-mat');
            before = computeFileChecksum(dataParamsPath);

            upgraded = loadDataParams(testCase.TempFolder).registration;

            testCase.verifyEqual(computeFileChecksum(dataParamsPath), before, ...
                'loading must not write');
            testCase.verifyEqual(upgraded.version, 2);
            testCase.verifyEqual(upgraded.reference.source.kind, 'legacy');
            testCase.verifyEqual(upgraded.reference.source.path, 'C:\old\ref.mat');
            testCase.verifyEqual(upgraded.reference.source.resourceUUID, ...
                testCase.Fixture.ImageReferenceUUID);
            testCase.verifyEqual(upgraded.reference.contentHash, v2.reference.contentHash);
            testCase.verifyEqual(upgraded.estimate.referenceHash, v2.reference.contentHash);
            testCase.verifyEqual(upgraded.estimate.movingSource, 'green.dat');

            % The legacy estimate can still be applied, and is then saved as v2.
            applyRegistrationTformOnFolder(testCase.TempFolder, ...
                'RequireUserConfirmation', false, 'OpenQCFigure', false);
            persisted = load(dataParamsPath, 'DataParams', '-mat');
            testCase.verifyEqual(persisted.DataParams.registration.version, 2);
            testCase.verifyTrue(persisted.DataParams.registration.apply.isRegistered);
        end

        function testPipelineInfoDescribesManagedReferenceRequirement(testCase)
            info = createRegistrationTform('pipelineInfo');

            testCase.verifySubstring(info.description, ...
                'active project-managed Image Reference');
            saveFolderInput = info.inputs(strcmp( ...
                {info.inputs.name}, 'SaveFolder'));
            testCase.verifySubstring(saveFolderInput.description, ...
                'active managed Image Reference');
            useFileInput = info.parameters(strcmp( ...
                {info.parameters.name}, 'UseFile'));
            testCase.verifySubstring(useFileInput.description, ...
                'ImageReference.sourceFile');
            testCase.verifySubstring(useFileInput.description, ...
                'several candidates are an error');
        end
    end

    methods (Access = private)
        function [resourceUUID, resource] = addAndActivateReference( ...
                testCase, imageData, description)
            ProjectInfo = testCase.Fixture.Store.getProjectInfo();
            ImageReference = genImageReferenceStruct( ...
                imageData, ...
                'Name', description, ...
                'ProjectUUID', ProjectInfo.projectUUID, ...
                'ProjectName', ProjectInfo.projectName, ...
                'SubjectUUID', testCase.Fixture.SubjectUUID, ...
                'SubjectID', testCase.Fixture.SubjectID, ...
                'SessionUUID', testCase.Fixture.SessionUUID, ...
                'SessionID', testCase.Fixture.SessionID, ...
                'Description', description, ...
                'SourceFolder', testCase.TempFolder, ...
                'SourceFile', fullfile(testCase.TempFolder, 'green.dat'), ...
                'SourceType', 'synthetic-test', ...
                'CreatedBy', 'TestCreateRegistrationTform');
            sourcePath = fullfile(testCase.TempFolder, ...
                ['reference_' char(java.util.UUID.randomUUID()) '.mat']);
            save(sourcePath, 'ImageReference', '-mat');

            resourceUUID = testCase.Fixture.Store.addImageReference( ...
                testCase.Fixture.SubjectID, sourcePath, struct( ...
                'displayName', description, ...
                'description', description));
            testCase.Fixture.Store.setActiveImageReference( ...
                testCase.Fixture.SubjectID, resourceUUID);
            resource = testCase.Fixture.Store.getResource(resourceUUID);
        end

        function clearActiveReferencePointer(testCase)
            SubjectInfo = testCase.Fixture.Store.getSubjectInfo( ...
                testCase.Fixture.SubjectID);
            idx = find(strcmp( ...
                {SubjectInfo.resourceRegistry.uuid}, ...
                testCase.Fixture.ImageReferenceUUID), 1, 'first');
            SubjectInfo.resourceRegistry(idx).status = 'available';
            SubjectInfo.resourceRegistry(idx).modifiedOn = datetime('now');
            SubjectInfo.activeImageReferenceUUID = '';
            SubjectInfo.modifiedOn = datetime('now');
            testCase.saveSubjectInfo(SubjectInfo);
        end

        function saveSubjectInfo(testCase, SubjectInfo)
            schema = getUMITProjectSchema();
            subjectPath = fullfile( ...
                testCase.Fixture.ProjectRoot, ...
                schema.folders.subjects, ...
                testCase.Fixture.SubjectID, ...
                schema.files.subjectMetadata);
            saveMatAtomic(subjectPath, ...
                schema.metadataVariables.subject, SubjectInfo);
        end

        function updateActiveReferenceSourceFile(testCase, sourceFile)
            managedPath = testCase.Fixture.ManagedReference.absolutePath;
            loaded = load(managedPath, 'ImageReference', '-mat');
            ImageReference = loaded.ImageReference;
            ImageReference.sourceFile = sourceFile;
            save(managedPath, 'ImageReference', '-mat');

            SubjectInfo = testCase.Fixture.Store.getSubjectInfo( ...
                testCase.Fixture.SubjectID);
            idx = find(strcmp( ...
                {SubjectInfo.resourceRegistry.uuid}, ...
                testCase.Fixture.ImageReferenceUUID), 1, 'first');
            SubjectInfo.resourceRegistry(idx).checksum = ...
                computeFileChecksum(managedPath);
            SubjectInfo.resourceRegistry(idx).modifiedOn = datetime('now');
            SubjectInfo.modifiedOn = datetime('now');
            testCase.saveSubjectInfo(SubjectInfo);
        end

        function [datPath, hash] = makeDatSource(testCase, name, scale)
            % A UMIT save folder whose green.dat averages to scale * the
            % fixture reference (a 'umitDat' reference source).
            folder = fullfile(fileparts(testCase.TempFolder), [name '_' char(java.util.UUID.randomUUID())]);
            mkdir(folder);
            testCase.addTeardown(@() rmdir(folder, 's'));
            image = single(scale) .* testCase.Fixture.ReferenceImage;
            writeTestDat(fullfile(folder, 'green.dat'), repmat(image, 1, 1, 3), 10, 20);
            AcqInfoStream = struct('Width', size(image, 2), 'Height', size(image, 1));
            save(fullfile(folder, 'AcqInfos.mat'), 'AcqInfoStream');
            datPath = fullfile(folder, 'green.dat');
            hash = computeImageContentHash(image);
        end

        function sums = datAndParamsChecksums(testCase)
            listing = [dir(fullfile(testCase.TempFolder, '*.dat')); ...
                dir(fullfile(testCase.TempFolder, 'DataParams.mat'))];
            sums = cellfun(@(n) computeFileChecksum(fullfile(testCase.TempFolder, n)), ...
                {listing.name}, 'UniformOutput', false);
        end

        function writeSingleDat(~, filePath, data, frameRateHz, exposureMsec)
            % Headered input (.dat header Phase 5a).
            writeTestDat(filePath, single(data), frameRateHz, exposureMsec);
        end

        function verifyPreflightFailure(testCase, operation, errorID)
            dataParamsPath = fullfile( ...
                testCase.TempFolder, 'DataParams.mat');
            checksumBefore = computeFileChecksum(dataParamsPath);
            qcBefore = dir(fullfile( ...
                testCase.TempFolder, 'registrationQC_*'));

            testCase.verifyError(operation, errorID);

            testCase.verifyEqual( ...
                computeFileChecksum(dataParamsPath), checksumBefore);
            qcAfter = dir(fullfile( ...
                testCase.TempFolder, 'registrationQC_*'));
            testCase.verifyEqual({qcAfter.name}, {qcBefore.name});
        end
    end
end
