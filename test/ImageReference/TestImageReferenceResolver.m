classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestImageReferenceResolver < matlab.unittest.TestCase
    %TESTIMAGEREFERENCERESOLVER Reference model (project redesign Phase 2).
    %
    %   resolveImageReferenceSource for the three source kinds, the .dat
    %   source validity rule, the aligned-source warning, read-only
    %   selection, content-hash identity across sources and re-saves, and
    %   standalone reference export/import.

    properties
        Root = ''
        SourceFolder = ''
        Fixture = []
    end

    methods (TestMethodSetup)
        function createFolders(testCase)
            addpath(genpath(fileparts(fileparts(fileparts(mfilename('fullpath'))))));
            testCase.Root = tempname;
            testCase.SourceFolder = fullfile(testCase.Root, 'source');
            mkdir(testCase.SourceFolder);
            frames = single(repmat(magic(8), 1, 1, 6)) + reshape(single(0:5), 1, 1, 6);
            writeTestDat(fullfile(testCase.SourceFolder, 'green.dat'), frames, 10, 5);
            AcqInfoStream = struct('Width', 8, 'Height', 8, 'FrameRateHz', 10);
            save(fullfile(testCase.SourceFolder, 'AcqInfos.mat'), 'AcqInfoStream');
        end
    end

    methods (TestMethodTeardown)
        function removeFolders(testCase)
            if isstruct(testCase.Fixture)
                for name = {'ProjectRoot', 'RigRoot', 'SaveFolder'}
                    if isfield(testCase.Fixture, name{1}) && isfolder(testCase.Fixture.(name{1}))
                        rmdir(testCase.Fixture.(name{1}), 's');
                    end
                end
            end
            if isfolder(testCase.Root)
                rmdir(testCase.Root, 's');
            end
        end
    end

    %% Content hash
    methods (Test)
        function testHashIsContentBasedAndStableAcrossResaves(testCase)
            image = rand(16, 12, 'single');
            hash = computeImageContentHash(image);
            testCase.verifyMatches(hash, '^sha256:[0-9a-f]{64}$');
            testCase.verifyEqual(computeImageContentHash(double(image)), hash, ...
                'same values, other class');
            withNaN = image;
            withNaN(1) = NaN;
            zeroed = image;
            zeroed(1) = 0;
            testCase.verifyEqual(computeImageContentHash(withNaN), computeImageContentHash(zeroed));
            testCase.verifyNotEqual(computeImageContentHash(image.'), hash, 'shape matters');

            % Re-saving a reference never changes its hash (unlike a file checksum).
            ImageReference = genImageReferenceStruct(image, 'Name', 'resave');
            first = fullfile(testCase.Root, 'a.mat');
            second = fullfile(testCase.Root, 'b.mat');
            save(first, 'ImageReference', '-mat');
            pause(1.1);
            save(second, 'ImageReference', '-mat');
            testCase.verifyNotEqual(computeFileChecksum(first), computeFileChecksum(second));
            a = resolveImageReferenceSource(struct('kind', 'file', 'path', first));
            b = resolveImageReferenceSource(struct('kind', 'file', 'path', second));
            testCase.verifyEqual(a.contentHash, hash);
            testCase.verifyEqual(b.contentHash, hash);
        end
    end

    %% Source kinds
    methods (Test)
        function testUmitDatSourceIsReadOnlyAndDescribed(testCase)
            before = testCase.folderChecksums(testCase.SourceFolder);

            reference = resolveImageReferenceSource(struct('kind', 'umitDat', ...
                'saveFolder', testCase.SourceFolder, 'datFile', 'green.dat'));

            testCase.verifyEqual(testCase.folderChecksums(testCase.SourceFolder), before, ...
                'selection must not write to the source folder');
            testCase.verifyEqual(reference.image, mean(single(repmat(magic(8), 1, 1, 6)) + ...
                reshape(single(0:5), 1, 1, 6), 3), 'AbsTol', 1e-5);
            testCase.verifyEqual(reference.contentHash, computeImageContentHash(reference.image));
            source = reference.source;
            testCase.verifyEqual(source.kind, 'umitDat');
            testCase.verifyEqual(source.path, fullfile(testCase.SourceFolder, 'green.dat'));
            testCase.verifyEqual(source.saveFolder, testCase.SourceFolder);
            testCase.verifyEqual(source.extractionRule, 'evenlySpacedMean');
            testCase.verifyEqual(source.extractionRuleVersion, 1);
            testCase.verifyFalse(source.sourceAligned);
            testCase.verifyEmpty(source.warning);
            validateReferenceSourceDescriptor(source);
            % datPath form gives the same result.
            direct = resolveImageReferenceSource(struct('kind', 'umitDat', ...
                'datPath', fullfile(testCase.SourceFolder, 'green.dat')));
            testCase.verifyEqual(direct.contentHash, reference.contentHash);
        end

        function testUmitDatNeedsAReadableAcqInfos(testCase)
            delete(fullfile(testCase.SourceFolder, 'AcqInfos.mat'));
            testCase.verifyError(@() resolveImageReferenceSource(struct('kind', 'umitDat', ...
                'datPath', fullfile(testCase.SourceFolder, 'green.dat'))), ...
                'Umitoolbox:ImageRegistration:invalidDatSource');

            fid = fopen(fullfile(testCase.SourceFolder, 'AcqInfos.mat'), 'w');
            fprintf(fid, 'not a mat file');
            fclose(fid);
            testCase.verifyError(@() resolveImageReferenceSource(struct('kind', 'umitDat', ...
                'datPath', fullfile(testCase.SourceFolder, 'green.dat'))), ...
                'Umitoolbox:ImageRegistration:invalidDatSource');
            testCase.verifyError(@() resolveImageReferenceSource(struct('kind', 'umitDat', ...
                'datPath', fullfile(testCase.SourceFolder, 'missing.dat'))), ...
                'Umitoolbox:ImageRegistration:invalidDatSource');
        end

        function testAlignedSourceIsAllowedWithAWarning(testCase)
            committed = testCase.markSourceFolderAligned();
            before = testCase.folderChecksums(testCase.SourceFolder);

            reference = testCase.verifyWarning(@() resolveImageReferenceSource(struct( ...
                'kind', 'umitDat', 'datPath', fullfile(testCase.SourceFolder, 'green.dat'))), ...
                'Umitoolbox:ImageRegistration:alignedSource');

            testCase.verifyEqual(testCase.folderChecksums(testCase.SourceFolder), before);
            source = reference.source;
            testCase.verifyTrue(source.sourceAligned);
            testCase.verifyEqual(source.sourceCommittedHash, committed.contentHash);
            testCase.verifyEqual(source.sourceCommittedName, committed.name);
            testCase.verifySubstring(source.warning, 'already aligned');
            testCase.verifySubstring(source.warning, 'interpolated a second time');
            testCase.verifySubstring(source.warning, committed.name);
            testCase.verifySubstring(source.warning, UMITLedgerCore.shortHash(committed.contentHash));
        end

        function testFileSourcePlainImages(testCase)
            image = single(magic(10));
            tifPath = fullfile(testCase.Root, 'plain.tif');
            imwrite(uint16(image), tifPath);
            matPath = fullfile(testCase.Root, 'plain.mat');
            img = image;
            other = image + 1;
            save(matPath, 'img', 'other');
            umtPath = fullfile(testCase.Root, 'plain.umt');
            saveData(umtPath, genUMTStruct(image, 'kind', 'image', 'entryName', 'map', ...
                'dimNames', {'Y', 'X'}));

            tif = resolveImageReferenceSource(struct('kind', 'file', 'path', tifPath), ...
                'Interactive', false);
            testCase.verifyEqual(tif.source.fileFormat, 'image');
            testCase.verifyEqual(tif.contentHash, computeImageContentHash(image));

            % Several MAT variables: never a dialog in non-interactive mode.
            testCase.verifyError(@() resolveImageReferenceSource(struct('kind', 'file', ...
                'path', matPath), 'Interactive', false), 'readImageReferenceSource:AmbiguousSource');
            mat = resolveImageReferenceSource(struct('kind', 'file', 'path', matPath, ...
                'variable', 'img'), 'Interactive', false);
            testCase.verifyEqual(mat.source.entry, 'img');
            testCase.verifyEqual(mat.contentHash, tif.contentHash);

            umt = resolveImageReferenceSource(struct('kind', 'file', 'path', umtPath), ...
                'Interactive', false);
            testCase.verifyEqual(umt.contentHash, tif.contentHash);
            testCase.verifyEmpty(umt.source.projectUUID);

            testCase.verifyError(@() resolveImageReferenceSource(struct('kind', 'file', ...
                'path', fullfile(testCase.SourceFolder, 'green.dat'))), ...
                'Umitoolbox:ImageRegistration:invalidSource');
            testCase.verifyError(@() resolveImageReferenceSource(struct('kind', 'other')), ...
                'Umitoolbox:ImageRegistration:invalidSource');
        end

        function testExportedReferenceWithAndWithoutProject(testCase)
            image = single(magic(7));
            free = genImageReferenceStruct(image, 'Name', 'free reference');
            testCase.verifyEqual(free.version, 2);
            freePath = exportImageReferenceFile(free, fullfile(testCase.Root, 'free.mat'));
            resolved = resolveImageReferenceSource(struct('kind', 'file', 'path', freePath));
            testCase.verifyEqual(resolved.source.fileFormat, 'umitReference');
            testCase.verifyEqual(resolved.name, 'free reference');
            testCase.verifyEqual(resolved.contentHash, computeImageContentHash(image));

            projectUUID = lower(char(java.util.UUID.randomUUID()));
            subjectUUID = lower(char(java.util.UUID.randomUUID()));
            owned = genImageReferenceStruct(image, 'Name', 'owned', ...
                'ProjectUUID', projectUUID, 'ProjectName', 'P', ...
                'SubjectUUID', subjectUUID, 'SubjectID', 'S1');
            ownedPath = exportImageReferenceFile(owned, fullfile(testCase.Root, 'owned.mat'));
            resolved = resolveImageReferenceSource(struct('kind', 'file', 'path', ownedPath));
            testCase.verifyEqual(resolved.source.projectUUID, projectUUID);
            testCase.verifyEqual(resolved.source.subjectID, 'S1');
            testCase.verifyEqual(resolved.contentHash, computeImageContentHash(image));

            cleared = exportImageReferenceFile(owned, fullfile(testCase.Root, 'cleared.mat'), ...
                'ClearIdentity', true);
            resolved = resolveImageReferenceSource(struct('kind', 'file', 'path', cleared));
            testCase.verifyEmpty(resolved.source.projectUUID);
            testCase.verifyEqual(resolved.contentHash, computeImageContentHash(image));
        end

        function testVersion1ReferenceIsUpgradedInMemory(testCase)
            image = single(magic(5));
            v2 = genImageReferenceStruct(image, 'Name', 'old', ...
                'ProjectUUID', lower(char(java.util.UUID.randomUUID())), 'ProjectName', 'P', ...
                'SubjectUUID', lower(char(java.util.UUID.randomUUID())), 'SubjectID', 'S1');
            v1 = rmfield(v2, 'contentHash');
            v1.version = 1;
            validateImageReferenceStruct(v1);
            upgraded = normalizeImageReferenceStruct(v1);
            testCase.verifyEqual(upgraded.version, 2);
            testCase.verifyEqual(upgraded.contentHash, v2.contentHash);
            % v1 still requires its identity; v2 does not.
            v1.projectUUID = '';
            testCase.verifyError(@() validateImageReferenceStruct(v1), ...
                'Umitoolbox:validateImageReferenceStruct:invalidInput');
        end
    end

    %% Project kind
    methods (Test)
        function testProjectAndDatSourcesGiveTheSameHash(testCase)
            % The same image from a project library and from a .dat
            % extraction has the same hash (brief section 1.3).
            saveFolder = fullfile(testCase.Root, 'session');
            testCase.Fixture = setupSyntheticRegistrationFolder(saveFolder);
            fx = testCase.Fixture;
            fromDat = resolveImageReferenceSource(struct('kind', 'umitDat', ...
                'datPath', fullfile(testCase.SourceFolder, 'green.dat')));
            ImageReference = genImageReferenceStruct(fromDat.image, 'Name', 'from dat', ...
                'ProjectUUID', fx.ProjectUUID, 'ProjectName', fx.Store.getProjectInfo().projectName, ...
                'SubjectUUID', fx.SubjectUUID, 'SubjectID', fx.SubjectID, ...
                'SourceFile', fullfile(testCase.SourceFolder, 'green.dat'));
            stage = fullfile(testCase.Root, 'stage.mat');
            save(stage, 'ImageReference', '-mat');
            resourceUUID = fx.Store.addImageReference(fx.SubjectID, stage, ...
                struct('displayName', 'from dat', 'description', 'library copy'));
            storeBefore = computeFileChecksum(fx.Store.getResource(resourceUUID).absolutePath);

            fromProject = resolveImageReferenceSource(struct('kind', 'project', ...
                'store', fx.Store, 'subjectID', fx.SubjectID, 'resourceUUID', resourceUUID));

            testCase.verifyEqual(fromProject.contentHash, fromDat.contentHash);
            testCase.verifyEqual(computeFileChecksum(fx.Store.getResource(resourceUUID).absolutePath), ...
                storeBefore, 'reading never writes to the store');
            source = fromProject.source;
            testCase.verifyEqual(source.kind, 'project');
            testCase.verifyEqual(source.resourceUUID, resourceUUID);
            testCase.verifyEqual(source.projectUUID, fx.ProjectUUID);
            testCase.verifyEqual(source.subjectUUID, fx.SubjectUUID);
            testCase.verifyEqual(fromProject.description, 'library copy');
            testCase.verifyEqual(fromProject.sourceFile, fullfile(testCase.SourceFolder, 'green.dat'));

            % Default: the subject's active reference.
            active = resolveImageReferenceSource(struct('kind', 'project', ...
                'store', fx.Store, 'subjectID', fx.SubjectID));
            testCase.verifyEqual(active.source.resourceUUID, fx.ImageReferenceUUID);

            % A library import of an exported file from another subject is refused,
            % and accepted once its identity is cleared.
            foreign = genImageReferenceStruct(single(magic(4)), 'Name', 'foreign', ...
                'ProjectUUID', fx.ProjectUUID, 'ProjectName', 'x', ...
                'SubjectUUID', lower(char(java.util.UUID.randomUUID())), 'SubjectID', 'Other');
            foreignPath = exportImageReferenceFile(foreign, fullfile(testCase.Root, 'foreign.mat'));
            testCase.verifyError(@() importImageReferenceFile(fx.Store, fx.SubjectID, foreignPath), ...
                'Umitoolbox:UMITProjectStore:imageReferenceIdentityMismatch');
            adopted = importImageReferenceFile(fx.Store, fx.SubjectID, foreignPath, ...
                'ClearIdentity', true);
            again = resolveImageReferenceSource(struct('kind', 'project', 'store', fx.Store, ...
                'subjectID', fx.SubjectID, 'resourceUUID', adopted));
            testCase.verifyEqual(again.contentHash, foreign.contentHash);
        end
    end

    %% Helpers
    methods
        function checksums = folderChecksums(~, folder)
            listing = dir(folder);
            listing = listing(~[listing.isdir]);
            checksums = struct('name', {listing.name}, 'sum', '');
            for k = 1:numel(listing)
                checksums(k).sum = computeFileChecksum(fullfile(folder, listing(k).name));
            end
        end

        function reference = markSourceFolderAligned(testCase)
            folder = testCase.SourceFolder;
            DataParams = createDataParams(folder);
            [~, referenceTemplate, estimateTemplate] = newRegistrationRecord();
            reference = referenceTemplate;
            [reference.contentHash, reference.image] = computeImageContentHash(single(magic(8)));
            reference.imageSizeYX = [8 8];
            reference.name = 'committed reference';
            reference.source = referenceSourceDescriptor('file', 'C:\\ref.tif');
            estimate = estimateTemplate;
            estimate.tform = affine2d(eye(3));
            estimate.referenceHash = reference.contentHash;
            DataParams.registration.reference = reference;
            DataParams.registration.estimate = estimate;
            DataParams.registration.apply.isRegistered = true;
            DataParams.registration.apply.appliedReferenceHash = reference.contentHash;
            DataParams.registration.apply.appliedOn = '2026-10-15 10:00:00';
            saveDataParams(folder, DataParams);
        end
    end
end
