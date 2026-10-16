classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestFreshImport < matlab.unittest.TestCase
    %TESTFRESHIMPORT Fresh import of a SaveFolder (project redesign Phase 4).
    %
    %   planFreshImport (deletion scope, brief section 4 Q2), applyFreshImport
    %   (backup, ledger events, deletion) and prepareFreshImport (consent,
    %   refusal before anything is touched). The importer integration is in
    %   TestRunImagesClassification / TestImportFromTif.

    properties
        Folder = ''
    end

    methods (TestMethodSetup)
        function createFolder(testCase)
            addpath(genpath(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))))));
            testCase.Folder = tempname;
            mkdir(testCase.Folder);
        end
    end

    methods (TestMethodTeardown)
        function removeFolder(testCase)
            if isfolder(testCase.Folder)
                rmdir(testCase.Folder, 's');
            end
        end
    end

    methods (Test)
        function testPlanListsKnownArtefactsAndKeepsTheRest(testCase)
            testCase.populateImportedFolder();

            plan = planFreshImport(testCase.Folder);

            testCase.verifyTrue(plan.hasPreviousImport);
            testCase.verifyEqual(sort(plan.toDelete), sort(testCase.derivedNames()));
            testCase.verifyEqual(sort(plan.keptRaw), sort({'img_0001.bin', 'info.txt', ...
                'Comments.txt', 'Snapshot_1.png', 'input_0001.tif'}));
            testCase.verifyEqual(sort(plan.keptLedgers), sort({'UMITRigLedger.json'}));
            testCase.verifyEqual(plan.keptBackups, {'bkp_20260101_000000.zip'});
            testCase.verifyEqual(sort(plan.unknown), sort({'my_notes.docx', ['userStuff' filesep]}));
            testCase.verifyTrue(plan.rigLedgerExists);
            testCase.verifyFalse(plan.projectLedgerBound);
            % Planning never writes.
            testCase.verifyEqual(sort(testCase.entries()), sort(testCase.allNames()));
        end

        function testEmptyFolderIsNotAPreviousImport(testCase)
            testCase.writeFile('img_0001.bin');
            testCase.writeFile('info.txt');
            plan = planFreshImport(testCase.Folder);
            testCase.verifyFalse(plan.hasPreviousImport);
            [report, ~] = prepareFreshImport(testCase.Folder, 'Origin', 'test');
            testCase.verifyEmpty(report);
            testCase.verifyEqual(sort(testCase.entries()), sort({'img_0001.bin', 'info.txt'}));
        end

        function testGenBackupZipsThenDeletesAndLogsTheRigLedger(testCase)
            testCase.populateImportedFolder();
            before = numel(UMITRigLedger.read(testCase.Folder).events);

            report = applyFreshImport(testCase.Folder, 'Backup', 'GENBACKUP', 'Origin', 'test');

            testCase.verifyEqual(sort(report.deleted), sort(testCase.derivedNames()));
            testCase.verifyTrue(isfile(report.backupPath));
            unzipped = fullfile(tempname);
            unzip(report.backupPath, unzipped);
            testCase.addTeardown(@() rmdir(unzipped, 's'));
            listing = dir(unzipped);
            testCase.verifyEqual(sort(setdiff({listing.name}, {'.', '..'})), ...
                sort(testCase.derivedNames()));
            remaining = sort(testCase.entries());
            expected = sort([setdiff(testCase.allNames(), testCase.derivedNames()), ...
                {[report.backupPath(numel(testCase.Folder)+2:end)]}]);
            testCase.verifyEqual(remaining, expected);
            info = UMITRigLedger.read(testCase.Folder);
            testCase.verifyEqual(numel(info.events), before + 1);
            testCase.verifyEqual(info.events{end}.type, 'freshImport');
            testCase.verifyTrue(report.rigEvent);
            testCase.verifyFalse(report.projectEvent);
            testCase.verifyFalse(UMITProjectLedger.exists(testCase.Folder), ...
                'an unbound folder gets no project ledger');
        end

        function testEraseMakesNoZipAndBoundFolderLogsTheProjectLedger(testCase)
            testCase.populateImportedFolder();
            u = @() lower(char(java.util.UUID.randomUUID()));
            UMITProjectLedger.bind(testCase.Folder, struct('bindingUUID', u(), ...
                'projectUUID', u(), 'subjectUUID', u(), 'sessionUUID', u(), ...
                'folderRole', 'processedDataFolder'), struct('note', 'x'), 'Origin', 'test');

            report = applyFreshImport(testCase.Folder, 'Backup', 'ERASE', 'Origin', 'test');

            testCase.verifyEmpty(report.backupPath);
            testCase.verifyEqual(numel(dir(fullfile(testCase.Folder, '*.zip'))), 1, ...
                'only the earlier backup');
            testCase.verifyTrue(report.projectEvent);
            info = UMITProjectLedger.read(testCase.Folder);
            testCase.verifyEqual(info.events{end}.type, 'freshImport');
            testCase.verifyEqual(info.state, 'Bound', 'the binding is unchanged');
        end

        function testBatchDefaultZipsBeforeDeleting(testCase)
            % Pipelines are not refused (decision 2026-10-16): in batch
            % MATLAB an empty backupOpts means GENBACKUP, so the deleted
            % files are kept in a zip. (On a desktop '' opens the dialog.)
            testCase.assumeTrue(batchStartupOptionUsed || ~usejava('desktop'), ...
                'Only meaningful in batch MATLAB.');
            testCase.populateImportedFolder();

            report = prepareFreshImport(testCase.Folder, 'Origin', 'test');

            testCase.verifyEqual(sort(report.deleted), sort(testCase.derivedNames()));
            testCase.verifyNotEmpty(report.backupPath, 'the default zips first');
            testCase.verifyTrue(isfile(report.backupPath));
        end

        function testBackupOptsIsTheOnlyPipelineDecision(testCase)
            % One parameter decides a fresh import (Bruno 2026-10-16): the
            % importers declare backupOpts (GENBACKUP default, or ERASE) and
            % no separate AllowFreshImport.
            for importer = {'ImagesClassification', 'importFromTif'}
                info = feval(importer{1}, 'pipelineInfo');
                names = {info.parameters.name};
                testCase.verifyFalse(ismember('AllowFreshImport', names), importer{1});
                backup = info.parameters(strcmp(names, 'backupOpts'));
                testCase.verifyNumElements(backup, 1, importer{1});
                testCase.verifyEqual(backup.default, 'GENBACKUP', importer{1});
                testCase.verifyEqual(sort(cellstr(backup.allowed)), ...
                    sort({'GENBACKUP', 'ERASE'}), importer{1});
            end
        end

        function testHistoryOutputsAreDeletedButOtherMatFilesAreKept(testCase)
            testCase.writeFile('AcqInfos.mat');
            testCase.writeFile('customResult.mat');
            testCase.writeFile('userData.mat');
            dataHistory = struct('file', {'customResult.mat'}, 'name', {'f'});
            save(fullfile(testCase.Folder, 'dataHistory.mat'), 'dataHistory');

            plan = planFreshImport(testCase.Folder);

            testCase.verifyEqual(sort(plan.toDelete), ...
                sort({'AcqInfos.mat', 'customResult.mat', 'dataHistory.mat'}));
            testCase.verifyEqual(plan.unknown, {'userData.mat'});
        end
    end

    methods
        function populateImportedFolder(testCase)
            for name = testCase.allNames()
                if endsWith(name{1}, filesep)
                    mkdir(fullfile(testCase.Folder, name{1}));
                elseif ~strcmp(name{1}, 'UMITRigLedger.json')
                    testCase.writeFile(name{1});
                end
            end
            UMITRigLedger.writeNoRig(testCase.Folder, 'Origin', 'test');
        end

        function names = derivedNames(~)
            names = {'red.dat', 'green.dat', 'HbO.dat', 'AcqInfos.mat', 'DataParams.mat', ...
                'events.mat', 'dataHistory.mat', 'pipeLog.mat', 'red.mat', 'green_info.mat', ...
                'map.umt', 'cortex.roi', 'registrationQC_2026-10-15_101010.fig', ...
                'registrationQC_2026-10-15_101010.png', 'manualAlignmentReport_20261015_101010.mat', ...
                'ImageAlignmentFailedReport_20261015_101010.txt', 'green_info.txt', 'green.tif', ...
                'red_C1_R1.tif', 'std_speckle.tiff', 'PMTMP_x_y.dat', 'red.dat.aligning.tmp', ...
                'green.dat.aligning.bak', 'pipeErrorLog_20261015.csv', ...
                'red_20261015_ROI_spatialMeasurements.csv'};
        end

        function names = allNames(testCase)
            names = [testCase.derivedNames(), {'img_0001.bin', 'info.txt', 'Comments.txt', ...
                'Snapshot_1.png', 'input_0001.tif', 'UMITRigLedger.json', ...
                'bkp_20260101_000000.zip', 'my_notes.docx', ['userStuff' filesep]}];
        end

        function writeFile(testCase, name)
            fid = fopen(fullfile(testCase.Folder, name), 'w');
            fprintf(fid, 'content of %s', name);
            fclose(fid);
        end

        function names = entries(testCase)
            listing = dir(testCase.Folder);
            listing = listing(~ismember({listing.name}, {'.', '..'}));
            names = {listing.name};
            names([listing.isdir]) = strcat(names([listing.isdir]), filesep);
        end

        function sums = checksums(testCase)
            listing = dir(testCase.Folder);
            listing = listing(~[listing.isdir]);
            sums = cellfun(@(n) computeFileChecksum(fullfile(testCase.Folder, n)), ...
                {listing.name}, 'UniformOutput', false);
        end
    end
end
