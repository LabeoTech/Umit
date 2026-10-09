classdef TestUMITRigOptics < matlab.unittest.TestCase
    %TESTUMITRIGOPTICS Focused schema-v3 optical and preflight tests.

    properties
        TempRoot = ''
        RigRoots = {}
        UserLibraryFiles = {}
        RigActivationFixture = struct()
    end

    methods (TestMethodSetup)
        function createFixture(testCase)
            addpath(genpath(fileparts(fileparts(fileparts(mfilename('fullpath'))))));
            testCase.TempRoot = tempname;
            mkdir(testCase.TempRoot);
            testCase.RigActivationFixture = struct();
        end
    end

    methods (TestMethodTeardown)
        function removeFixture(testCase)
            deactivateRigTemporarily(testCase.RigActivationFixture);
            for iRig = 1:numel(testCase.RigRoots)
                if isfolder(testCase.RigRoots{iRig})
                    rmdir(testCase.RigRoots{iRig}, 's');
                end
            end
            for iFile = 1:numel(testCase.UserLibraryFiles)
                if isfile(testCase.UserLibraryFiles{iFile})
                    delete(testCase.UserLibraryFiles{iFile});
                end
            end
            if isfolder(testCase.TempRoot)
                rmdir(testCase.TempRoot, 's');
            end
        end
    end

    methods (Test)
        function testBuiltInLibraryAndLegacyNumericalRegression(testCase)
            red = UMITRigStore.getSpectrum('illumination', 'LED_632nm');
            green = UMITRigStore.getSpectrum('illumination', 'LED_521nm');
            yellow = UMITRigStore.getSpectrum('illumination', 'LED_593nm');
            camera = UMITRigStore.getSpectrum('camera', 'PF1024');
            filterSet = UMITRigStore.getFilterSet('none');

            testCase.verifySize(red.response, [301 1]);
            testCase.verifyEqual(red.wavelengthNm, (400:700)');
            testCase.verifyGreaterThanOrEqual(min(red.response), 0);
            testCase.verifyLessThanOrEqual(max(red.response), 1);
            testCase.verifyEqual(filterSet.excitationSpectrumID, '');

            optical = struct( ...
                'wavelengthNm', (400:700)', ...
                'activeRows', true(3, 1), ...
                'illuminationResponse', [red.response'; green.response'; yellow.response'], ...
                'cameraResponse', repmat(camera.response', 3, 1), ...
                'excitationResponse', ones(1, 301), ...
                'emissionResponse', ones(1, 301));
            actual = ioi_epsilon_pathlength('Hillman', 100, 60, 40, optical);
            legacy = ioi_epsilon_pathlength( ...
                'Hillman', 100, 60, 40, 'none', 'D1024');
            expected = [607.4858 4677.0937; 3317.8867 3753.0432; 2450.6303 5176.1897];
            testCase.verifyEqual(actual, expected, 'AbsTol', 1e-3);
            testCase.verifyEqual(legacy, expected, 'AbsTol', 1e-3);
            testCase.verifyEqual(actual, legacy, 'AbsTol', 1e-10);
            testCase.verifyError(@() ioi_epsilon_pathlength( ...
                'Hillman', 100, 60, 40, 'none', 'UnknownCamera'), ...
                'Umitoolbox:ioi_epsilon_pathlength:UnknownCamera');
        end

        function testLegacyMigrationProducesExpectedRepertoire(testCase)
            outputRoot = fullfile(testCase.TempRoot, 'migratedLibrary');
            report = migrateLegacyOpticalLibrary(outputRoot, ...
                fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), ...
                'IOIAnalysis', 'SystemOpticalParams'));
            testCase.verifyEqual(report.spectrumCount, 16);
            testCase.verifyEqual(report.filterSetCount, 3);
            values = readmatrix(fullfile(outputRoot, 'spectra', ...
                'illumination', 'LED_632nm.txt'), ...
                'FileType', 'text', 'CommentStyle', '#');
            testCase.verifySize(values, [301 2]);
            testCase.verifyEqual(values(:, 1), (400:700)');
            definition = jsondecode(fileread(fullfile( ...
                outputRoot, 'filterSets', 'GCaMP.json')));
            testCase.verifyEqual(definition.excitationSpectrumID, 'FF01496LP');
            testCase.verifyEqual(definition.emissionSpectrumID, 'FF01496LP');
        end

        function testSpectrumImportResamplesAndNormalizes(testCase)
            source = fullfile(testCase.TempRoot, 'source.csv');
            writematrix([(390:10:710)' linspace(0.1, 2, 33)'], source);
            suffix = strrep(char(java.util.UUID.randomUUID()), '-', '');
            spectrumID = ['UnitTestSpectrum_' suffix(1:8)];

            spectrum = UMITRigStore.importSpectrum( ...
                source, 'camera', spectrumID, struct('displayName', 'Unit Test'));
            testCase.UserLibraryFiles = {spectrum.file};

            testCase.verifySize(spectrum.response, [301 1]);
            testCase.verifyEqual(spectrum.wavelengthNm, (400:700)');
            testCase.verifyEqual(max(spectrum.response), 1, 'AbsTol', eps);
            testCase.verifyEqual(spectrum.metadata.displayName, 'Unit Test');
        end

        function testSpectrumEnumerationIncludesBuiltInAndUserEntries(testCase)
            source = fullfile(testCase.TempRoot, 'listed-source.csv');
            writematrix([(400:700)' ones(301, 1)], source);
            suffix = strrep(char(java.util.UUID.randomUUID()), '-', '');
            spectrumID = ['ListedSpectrum_' suffix(1:8)];
            spectrum = UMITRigStore.importSpectrum(source, 'camera', ...
                spectrumID, struct( ...
                    'displayName', 'Listed Camera', ...
                    'manufacturer', 'Unit Test Maker', ...
                    'model', 'Unit Test Model'));
            testCase.UserLibraryFiles = {spectrum.file};

            cameraSpectra = UMITRigStore.listSpectra('camera');
            builtInRow = strcmp(cameraSpectra.ID, 'PF1024');
            userRow = strcmp(cameraSpectra.ID, spectrumID);
            testCase.verifyTrue(any(builtInRow));
            testCase.verifyTrue(any(userRow));
            testCase.verifyEqual(cameraSpectra.Origin(builtInRow), "builtIn");
            testCase.verifyEqual(cameraSpectra.Origin(userRow), "user");
            testCase.verifyEqual(cameraSpectra.DisplayName(userRow), "Listed Camera");
            testCase.verifyEqual(cameraSpectra.Manufacturer(userRow), "Unit Test Maker");
            testCase.verifyEqual(cameraSpectra.Model(userRow), "Unit Test Model");
            testCase.verifyFalse(any(contains(cameraSpectra.Properties.VariableNames, ...
                {'File', 'Folder', 'Path'}, 'IgnoreCase', true)));

            illuminationSpectra = UMITRigStore.listSpectra('illumination');
            filterSpectra = UMITRigStore.listSpectra('filter');
            testCase.verifyTrue(any(strcmp(illuminationSpectra.ID, 'LED_632nm')));
            testCase.verifyTrue(any(strcmp(filterSpectra.ID, 'FF01496LP')));
        end

        function testFilterSetEnumerationIsStorageIndependent(testCase)
            filterSets = UMITRigStore.listFilterSets();
            gcampRow = strcmp(filterSets.ID, 'GCaMP');

            testCase.verifyTrue(any(gcampRow));
            testCase.verifyEqual(filterSets.Origin(gcampRow), "builtIn");
            testCase.verifyEqual(filterSets.DisplayName(gcampRow), "GCaMP");
            testCase.verifyEqual( ...
                filterSets.ExcitationSpectrumID(gcampRow), "FF01496LP");
            testCase.verifyEqual( ...
                filterSets.EmissionSpectrumID(gcampRow), "FF01496LP");
            testCase.verifyFalse(any(contains(filterSets.Properties.VariableNames, ...
                {'File', 'Folder', 'Path'}, 'IgnoreCase', true)));
        end

        function testUserFilterSetUpdateAndRemoval(testCase)
            suffix = strrep(char(java.util.UUID.randomUUID()), '-', '');
            filterSetID = ['UserFilter_' suffix(1:8)];
            definition = struct( ...
                'id', filterSetID, ...
                'displayName', 'User Filter Set', ...
                'excitationSpectrumID', 'FF01496LP', ...
                'emissionSpectrumID', 'FF01496LP');

            UMITRigStore.importFilterSet(definition);
            imported = UMITRigStore.getFilterSet(filterSetID);
            testCase.UserLibraryFiles{end+1} = imported.file;
            testCase.verifyError(@() UMITRigStore.importFilterSet(definition), ...
                'Umitoolbox:UMITRigStore:filterSetAlreadyExists');

            definition.displayName = 'Updated User Filter Set';
            UMITRigStore.updateFilterSet(filterSetID, definition);
            updated = UMITRigStore.getFilterSet(filterSetID);
            testCase.verifyEqual(updated.displayName, ...
                'Updated User Filter Set');
            testCase.verifyError(@() UMITRigStore.removeFilterSet('GCaMP'), ...
                'Umitoolbox:UMITRigStore:filterSetRemovalFailed');

            UMITRigStore.removeFilterSet(filterSetID);
            testCase.verifyFalse(isfile(imported.file));
            listed = UMITRigStore.listFilterSets();
            testCase.verifyFalse(any(strcmp(listed.ID, filterSetID)));
        end

        function testSpectrumImportRejectsExistingID(testCase)
            source = fullfile(testCase.TempRoot, 'source.csv');
            writematrix([(400:700)' ones(301, 1)], source);
            spectrumID = ['NoOverwrite_' char(java.util.UUID.randomUUID())];
            spectrumID = spectrumID(1:20);
            spectrum = UMITRigStore.importSpectrum(source, 'camera', spectrumID, struct());
            testCase.UserLibraryFiles = {spectrum.file};
            testCase.verifyError(@() UMITRigStore.importSpectrum( ...
                source, 'camera', spectrumID, struct()), ...
                'Umitoolbox:UMITRigStore:spectrumAlreadyExists');
        end

        function testSpectrumImportRejectsMalformedInput(testCase)
            source = fullfile(testCase.TempRoot, 'bad.csv');
            writematrix([(500:10:600)' ones(11, 1)], source);
            testCase.verifyError(@() UMITRigStore.importSpectrum( ...
                source, 'filter', 'BadSpectrum', struct()), ...
                'Umitoolbox:UMITRigStore:invalidSpectrum');
        end

        function testFilterSetRejectsDanglingSpectrumReference(testCase)
            filterSpectrum = UMITRigStore.getSpectrum('filter', 'FF01496LP');
            testCase.verifySize(filterSpectrum.response, [301 1]);

            filterSet = struct( ...
                'id', 'DanglingUnitTest', ...
                'displayName', 'Dangling unit test', ...
                'excitationSpectrumID', 'MissingFilterSpectrum', ...
                'emissionSpectrumID', '');
            testCase.verifyError(@() UMITRigStore.importFilterSet(filterSet), ...
                'Umitoolbox:UMITRigStore:spectrumNotFound');
        end

        function testReferencedSpectrumCannotBeRemoved(testCase)
            source = fullfile(testCase.TempRoot, 'camera.csv');
            writematrix([(400:700)' ones(301, 1)], source);
            suffix = strrep(char(java.util.UUID.randomUUID()), '-', '');
            spectrumID = ['Referenced_' suffix(1:8)];
            spectrum = UMITRigStore.importSpectrum( ...
                source, 'camera', spectrumID, struct());
            testCase.UserLibraryFiles = {spectrum.file};

            store = testCase.createConfiguredRig();
            info = store.getRigInfo();
            info.cameras(1).spectrumID = spectrumID;
            store.setCameras(info.cameras);
            testCase.verifyError(@() UMITRigStore.removeSpectrum( ...
                'camera', spectrumID), 'Umitoolbox:UMITRigStore:spectrumInUse');

            info.cameras(1).spectrumID = 'PF1024';
            store.setCameras(info.cameras);
            UMITRigStore.removeSpectrum('camera', spectrumID);
            testCase.verifyFalse(isfile(spectrum.file));
        end

        function testOpticalResolutionUsesAcquisitionColorAndCameraIndex(testCase)
            store = testCase.createConfiguredRig();
            acq = testCase.makeAcqInfo(store, 'Amber', 2);
            optical = store.resolveOpticalConfiguration( ...
                acq, {'red','amber'}, 'GCaMP', '');

            testCase.verifyEqual([optical.channels.camIdx], [1 2]);
            testCase.verifyEqual({optical.channels.name}, {'red','yellow'});
            testCase.verifyTrue(optical.activeRows(1));
            testCase.verifyTrue(optical.activeRows(3));
            testCase.verifyFalse(optical.activeRows(2));
        end

        function testMissingCameraSpectrumDiagnostic(testCase)
            store = testCase.createConfiguredRig();
            info = store.getRigInfo();
            info.cameras(1).spectrumID = '';
            store.setCameras(info.cameras);
            acq = testCase.makeAcqInfo(store, 'green', 1);
            testCase.verifyError(@() store.resolveOpticalConfiguration( ...
                acq, {'red','green'}, 'none', ''), ...
                'Umitoolbox:UMITRigStore:missingCameraSpectrum');
        end

        function testMissingIlluminationSpectrumDiagnostic(testCase)
            store = testCase.createConfiguredRig();
            info = store.getRigInfo();
            info.illuminations(strcmp({info.illuminations.name}, 'green')).spectrumID = '';
            store.setIlluminations(info.illuminations);
            acq = testCase.makeAcqInfo(store, 'green', 1);
            testCase.verifyError(@() store.resolveOpticalConfiguration( ...
                acq, {'red','green'}, 'none', ''), ...
                'Umitoolbox:UMITRigStore:missingIlluminationSpectrum');
        end

        function testCoregistrationPayloadValidationAndAtomicActivation(testCase)
            store = testCase.createConfiguredRig();
            badFile = fullfile(testCase.TempRoot, 'bad.mat');
            payload = 1;
            save(badFile, 'payload', '-mat');
            testCase.verifyError(@() store.addCameraCoregistration(badFile, struct()), ...
                'Umitoolbox:UMITRigStore:addResourceFailed');

            badInfoFile = fullfile(testCase.TempRoot, 'badInfo.mat');
            tform = affine2d(eye(3));
            tformInfo = repmat(struct('value', 1), 1, 2);
            save(badInfoFile, 'tform', 'tformInfo', '-mat');
            testCase.verifyError(@() store.addCameraCoregistration( ...
                badInfoFile, struct()), ...
                'Umitoolbox:UMITRigStore:addResourceFailed');

            first = testCase.createTform('first.mat', 1);
            second = testCase.createTform('second.mat', 2);
            firstUUID = store.addCameraCoregistration(first, struct());
            secondUUID = store.importAndActivateCameraCoregistration(second, struct());
            testCase.verifyEqual(store.getActiveCameraCoregistration().uuid, secondUUID);
            testCase.verifyEqual(store.getResource(firstUUID).status, 'available');
        end

        function testArchivedRigRejectsMutationUntilRestored(testCase)
            store = testCase.createConfiguredRig();
            replacement = testCase.createConfiguredRig();
            store.archiveRig(replacement.getRigInfo().uuid);
            testCase.verifyError(@() store.updateRigMetadata( ...
                struct('displayName', 'Not allowed')), ...
                'Umitoolbox:UMITRigStore:archivedRigReadOnly');
            store.restoreRig();
            store.updateRigMetadata(struct('displayName', 'Restored'));
            testCase.verifyEqual(store.getRigInfo().displayName, 'Restored');
        end

        function testStoreLevelDefaultInvariant(testCase)
            first = testCase.createConfiguredRig();
            second = testCase.createConfiguredRig();
            testCase.RigActivationFixture = ...
                activateRigTemporarily(first.getRigInfo().uuid);
            UMITRigStore.setDefaultRig(second.getRigInfo().uuid);

            resolved = UMITRigStore.getDefaultRig();
            testCase.verifyEqual(resolved.getRigInfo().uuid, second.getRigInfo().uuid);
            rigs = UMITRigStore.listRigs();
            owned = ismember(rigs.RigUUID, string( ...
                {first.getRigInfo().uuid, second.getRigInfo().uuid}));
            testCase.verifyEqual(sum(rigs.IsDefault(owned)), 1);
        end
    end

    methods (Access = private)
        function store = createConfiguredRig(testCase)
            suffix = strrep(char(java.util.UUID.randomUUID()), '-', '');
            cameras = struct( ...
                'index', {1, 2}, ...
                'displayName', {'Camera 1', 'Camera 2'}, ...
                'manufacturer', {'', ''}, ...
                'model', {'D1024', 'D1024'}, ...
                'serialNumber', {'', ''}, ...
                'spectrumID', {'PF1024', 'PF1024'});
            illuminations = struct( ...
                'name', {'red','green','yellow'}, ...
                'displayName', {'Red','Green','Yellow'}, ...
                'manufacturer', {'','',''}, ...
                'model', {'','',''}, ...
                'spectrumID', {'LED_632nm','LED_521nm','LED_593nm'});
            store = UMITRigStore.create(struct( ...
                'rigID', ['RigOptics_' suffix(1:8)], ...
                'cameras', cameras, ...
                'illuminations', illuminations));
            testCase.RigRoots{end+1} = store.RigRoot;
        end

        function acq = makeAcqInfo(~, store, secondColor, secondCamIdx)
            info = store.getRigInfo();
            acq = struct();
            acq.rigUUID = info.uuid;
            acq.rigID = info.rigID;
            acq.ImportedChannels = struct( ...
                'DatFile', {'red.dat', 'second.dat'}, ...
                'Color', {'Red', secondColor}, ...
                'CamIdx', {1, secondCamIdx});
        end

        function filePath = createTform(testCase, fileName, offset)
            filePath = fullfile(testCase.TempRoot, fileName);
            tform = affine2d([1 0 0; 0 1 0; offset 0 1]);
            tformInfo = struct('CandidateSource', 'unit-test');
            save(filePath, 'tform', 'tformInfo', '-mat');
        end
    end
end
