classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestExportMatToJson < matlab.unittest.TestCase
    %TESTEXPORTMATTOJSON JSON copies of report and log MAT files.
    %
    %   exportMatToJson turns a report or log MAT file (alignment report, pipeLog,
    %   dataHistory, motion-correction log) into a readable JSON file for users of
    %   the compiled application, which has no MATLAB to open the file with. The
    %   tests build values of the same shapes as the project's writers produce
    %   and read the JSON back with jsondecode. The converter refuses to run
    %   outside the compiled application unless Force is set, which every test
    %   except the gate test does.

    properties
        Folder char
    end

    methods (TestMethodSetup)
        function createFolder(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.Folder = fixture.Folder;
        end
    end

    methods (Test)
        function refusesToRunOutsideTheCompiledApplication(testCase)
            matFile = testCase.writeMat('a.mat', struct('x', 1));

            testCase.verifyError(@() exportMatToJson(matFile, fullfile(testCase.Folder, 'a.json')), ...
                'Umitoolbox:MatToJson:NotDeployed');
            testCase.verifyFalse(isfile(fullfile(testCase.Folder, 'a.json')));
        end

        function anAlignmentReportShapedStructIsWrittenInFull(testCase)
            report = struct();
            report.createdOn = datetime(2026, 11, 25, 10, 30, 15);
            report.finishedOn = NaT;
            report.toolName = 'ImageAlignmentTool';
            report.targetFolder = 'D:\data\Session 01';
            report.oldSizeYX = [64 96];
            report.newSizeYX = [128 192];
            report.finalTform = affinetform2d([1.1 0 3; 0 1.1 -2; 0 0 1]);
            report.datFilesTransformed = ["green.dat"; "red.dat"];
            report.datFilesSkipped = strings(0, 1);
            report.warnings = "One warning with a quote "" and an accent é";
            report.scaleInfo = struct('scale', 1.1, 'isSimilarity', true);
            matFile = testCase.writeMat('manualAlignmentReport_20261125_103015.mat', ...
                struct('report', report));

            result = exportMatToJson(matFile, fullfile(testCase.Folder, 'r.json'), 'Force', true);
            decoded = jsondecode(fileread(result.JsonFile));

            testCase.verifyEqual(decoded.sourceFile, 'manualAlignmentReport_20261125_103015.mat');
            r = decoded.variables.report;
            testCase.verifyTrue(startsWith(r.createdOn, '2026-11-25T10:30:15'));
            testCase.verifyEqual(r.toolName, 'ImageAlignmentTool');
            testCase.verifyEqual(r.oldSizeYX(:).', [64 96]);
            testCase.verifyEqual(r.newSizeYX(:).', [128 192]);
            testCase.verifyEqual(r.datFilesTransformed, {'green.dat'; 'red.dat'});
            testCase.verifyEmpty(r.datFilesSkipped);
            testCase.verifySubstring(r.warnings, 'accent é');
            testCase.verifyEqual(r.scaleInfo.scale, 1.1);
            testCase.verifyTrue(r.scaleInfo.isSimilarity);
            testCase.verifyEqual(r.finalTform.x_class, 'affinetform2d');
            testCase.verifyEmpty(result.SummarisedVariables);
        end

        function aTableLogIsWrittenAsRows(testCase)
            pipeLog = table( ...
                datetime(2026, 11, 25, 9, 0, 0) + minutes([0; 5; 10]), ...
                ["normalizeBSLN_1"; "normalizeLPF_1"; "spatialGaussFilt_1"], ...
                categorical(["completed"; "completed"; "failed"]), ...
                [true; true; false], ...
                {"ok"; ""; "Error: file missing"}, ...
                'VariableNames', {'Time', 'Step', 'Status', 'Success', 'Message'});
            matFile = testCase.writeMat('pipeLog.mat', struct('pipeLog', pipeLog));

            result = exportMatToJson(matFile, fullfile(testCase.Folder, 'p.json'), 'Force', true);
            decoded = jsondecode(fileread(result.JsonFile));

            t = decoded.variables.pipeLog;
            testCase.verifyEqual(t.x_type, 'table');
            testCase.verifyEqual(t.variableNames, {'Time'; 'Step'; 'Status'; 'Success'; 'Message'});
            testCase.verifyNumElements(t.rows, 3);
            testCase.verifyEqual(t.rows(3).Step, 'spatialGaussFilt_1');
            testCase.verifyEqual(t.rows(3).Status, 'failed');
            testCase.verifyFalse(t.rows(3).Success);
            testCase.verifyEqual(t.rows(3).Message, 'Error: file missing');
            testCase.verifyTrue(startsWith(t.rows(2).Time, '2026-11-25T09:05:00'));
        end

        function aStructArrayLogIsWrittenAsAnArrayOfObjects(testCase)
            % dataHistory and MotionCorrectionLog are struct arrays.
            dataHistory = struct('file', {'a.dat'; 'b.dat'; 'c.dat'}, ...
                'info', {struct('step', 'A'); struct('step', 'B'); struct('step', 'C')});
            MotionCorrectionLog = struct('file', {'green.dat'; 'red.dat'}, ...
                'appliedOn', {datetime(2026, 11, 24, 8, 0, 0); datetime(2026, 11, 24, 8, 5, 0)});
            matFile = testCase.writeMat('logs.mat', struct( ...
                'dataHistory', dataHistory, 'MotionCorrectionLog', MotionCorrectionLog));

            result = exportMatToJson(matFile, fullfile(testCase.Folder, 'l.json'), 'Force', true);
            decoded = jsondecode(fileread(result.JsonFile));

            testCase.verifyNumElements(decoded.variables.dataHistory, 3);
            testCase.verifyEqual(decoded.variables.dataHistory(2).file, 'b.dat');
            testCase.verifyEqual(decoded.variables.dataHistory(3).info.step, 'C');
            testCase.verifyNumElements(decoded.variables.MotionCorrectionLog, 2);
            testCase.verifyTrue(startsWith(decoded.variables.MotionCorrectionLog(2).appliedOn, ...
                '2026-11-24T08:05:00'));
        end

        function smallArraysAreWrittenInFullAndLargeOnesAreSummarised(testCase)
            small = magic(4);
            large = reshape(1:90000, [300 300]);
            matFile = testCase.writeMat('arrays.mat', struct('small', small, 'large', large));

            result = exportMatToJson(matFile, fullfile(testCase.Folder, 'x.json'), 'Force', true);
            decoded = jsondecode(fileread(result.JsonFile));

            testCase.verifyEqual(decoded.variables.small, small);
            descriptor = decoded.variables.large;
            testCase.verifyEqual(descriptor.x_type, 'array');
            testCase.verifyEqual(descriptor.size(:).', [300 300]);
            testCase.verifyEqual(descriptor.min, 1);
            testCase.verifyEqual(descriptor.max, 90000);
            testCase.verifyEqual(descriptor.mean, 45000.5);
            testCase.verifySubstring(descriptor.note, 'Full');
            testCase.verifyEqual(result.SummarisedVariables, "large");
            testCase.verifyEqual(decoded.summarisedVariables, {'large'});
        end

        function theFullOptionAndTheLimitOptionChangeWhatIsSummarised(testCase)
            values = (1:50);
            matFile = testCase.writeMat('v.mat', struct('values', values));

            low = exportMatToJson(matFile, fullfile(testCase.Folder, 'low.json'), ...
                'Force', true, 'MaxArrayElements', 10);
            full = exportMatToJson(matFile, fullfile(testCase.Folder, 'full.json'), ...
                'Force', true, 'MaxArrayElements', 10, 'Full', true);

            testCase.verifyEqual(low.SummarisedVariables, "values");
            testCase.verifyEmpty(full.SummarisedVariables);
            decoded = jsondecode(fileread(full.JsonFile));
            testCase.verifyEqual(decoded.variables.values(:).', values);
        end

        function nanAndInfAreWrittenAsNull(testCase)
            matFile = testCase.writeMat('n.mat', struct('a', [1 NaN Inf], 'b', NaN));

            result = exportMatToJson(matFile, fullfile(testCase.Folder, 'n.json'), 'Force', true);
            text = fileread(result.JsonFile);
            decoded = jsondecode(text);

            testCase.verifyEqual(decoded.variables.a(1), 1);
            testCase.verifyTrue(isnan(decoded.variables.a(2)));
            testCase.verifyTrue(isnan(decoded.variables.a(3)));
            testCase.verifyEmpty(decoded.variables.b, 'NaN is written as null.');
            testCase.verifySubstring(text, 'null');
        end

        function functionHandlesMapsAndObjectsAreReadable(testCase)
            map = containers.Map({'alpha', 'beta'}, {1, 'two'});
            handle = @(x) x+1;
            tform = affinetform2d([1 0 5; 0 1 6; 0 0 1]);
            matFile = fullfile(testCase.Folder, 'o.mat');
            save(matFile, 'map', 'handle', 'tform');

            result = exportMatToJson(matFile, fullfile(testCase.Folder, 'o.json'), 'Force', true);
            decoded = jsondecode(fileread(result.JsonFile));

            testCase.verifyEqual(decoded.variables.map.alpha, 1);
            testCase.verifyEqual(decoded.variables.map.beta, 'two');
            testCase.verifySubstring(decoded.variables.handle, 'x+1');
            testCase.verifyEqual(decoded.variables.tform.x_class, 'affinetform2d');
            testCase.verifyEqual(decoded.variables.tform.A, [1 0 5; 0 1 6; 0 0 1]);
        end

        function emptyValuesAndTextAreKept(testCase)
            matFile = testCase.writeMat('e.mat', struct('emptyNumeric', [], 'emptyText', '', ...
                'emptyCell', {{}}, 'emptyStruct', struct([]), 'name', 'ok'));

            result = exportMatToJson(matFile, fullfile(testCase.Folder, 'e.json'), 'Force', true);
            decoded = jsondecode(fileread(result.JsonFile));

            testCase.verifyEmpty(decoded.variables.emptyNumeric);
            testCase.verifyEqual(decoded.variables.emptyText, '');
            testCase.verifyEmpty(decoded.variables.emptyCell);
            testCase.verifyEmpty(decoded.variables.emptyStruct);
            testCase.verifyEqual(decoded.variables.name, 'ok');
        end

        function theJsonKeysAreNotLeftAsPlaceholders(testCase)
            matFile = testCase.writeMat('k.mat', struct('t', table([1; 2], 'VariableNames', {'A'})));

            result = exportMatToJson(matFile, fullfile(testCase.Folder, 'k.json'), 'Force', true);
            text = fileread(result.JsonFile);

            testCase.verifyFalse(contains(text, 'JSONKEY'));
            testCase.verifySubstring(text, '"_type"');
        end

        function theSourceFileIsNeverModifiedOrOverwritten(testCase)
            matFile = testCase.writeMat('s.mat', struct('x', 1));
            before = dir(matFile);
            checksumBefore = computeFileChecksum(matFile);

            testCase.verifyError(@() exportMatToJson(matFile, matFile, 'Force', true), ...
                'Umitoolbox:MatToJson:SameFile');
            exportMatToJson(matFile, fullfile(testCase.Folder, 's.json'), 'Force', true);

            testCase.verifyEqual(computeFileChecksum(matFile), checksumBefore);
            testCase.verifyEqual(dir(matFile).datenum, before.datenum);
        end

        function aMissingFileAndANonMatFileAreReportedClearly(testCase)
            testCase.verifyError(@() exportMatToJson(fullfile(testCase.Folder, 'nope.mat'), ...
                fullfile(testCase.Folder, 'n.json'), 'Force', true), 'Umitoolbox:MatToJson:MissingFile');

            notMat = fullfile(testCase.Folder, 'broken.mat');
            fid = fopen(notMat, 'w');
            fprintf(fid, 'not a mat file');
            fclose(fid);
            testCase.verifyError(@() exportMatToJson(notMat, fullfile(testCase.Folder, 'b.json'), ...
                'Force', true), 'Umitoolbox:MatToJson:NotAMatFile');
            testCase.verifyFalse(isfile(fullfile(testCase.Folder, 'b.json')));
        end
    end

    methods (Access = private)
        function file = writeMat(testCase, name, variables)
            %WRITEMAT Save the fields of a struct as variables of a MAT file.
            file = fullfile(testCase.Folder, name);
            save(file, '-struct', 'variables');
        end
    end
end
