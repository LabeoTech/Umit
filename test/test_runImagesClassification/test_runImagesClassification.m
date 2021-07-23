classdef test_runImagesClassification < matlab.unittest.TestCase
    properties (TestParameter)
        rawFolder = {'\test\test_runImagesClassification\rawData\'};
        expOutFolder = {'\test\test_runImagesClassification\exp_outFiles\'};
        outFolder = {'\test\test_runImagesClassification\outputFolder\'};
        opts = {struct('BinningSpatial', 1, 'BinningTemp', 1, 'b_SubROI', false, 'b_IgnoreStim', false), ...
             struct('BinningSpatial', 1, 'BinningTemp', 1, 'b_SubROI', false, 'b_IgnoreStim', true)}
    end
    
    methods(Test)
        function realSolution(testCase, rawFolder, outFolder, expOutFolder, opts)
            root = getenv('Umitoolbox');  
            run_ImagesClassification(fullfile(root,rawFolder), fullfile(root,outFolder),opts);
            % Test .DAT files:
            dat_files = dir([fullfile(root,outFolder) '*.dat']);
            for i = 1:length(dat_files)
                fid_exp = fopen(fullfile(root, expOutFolder, dat_files(i).name), 'r');
                expFile = fread(fid_exp, 'single=>single');
                fid_out = fopen(fullfile(dat_files(i).folder, dat_files(i).name), 'r');
                outFile = fread(fid_out, 'single=>single');
                testCase.verifyEqual(expFile, outFile);
                fclose(fid_exp);
                fclose(fid_out);
            end
            % Test .MAT files:
            mat_files = dir([fullfile(root,outFolder) '*.mat']);
            for i = 1:length(mat_files)
                Exp = load(fullfile(root, expOutFolder, mat_files(i).name));
                Out = load(fullfile(mat_files(i).folder, mat_files(i).name));
                try
                    Exp = rmfield(Exp,'fileUUID');
                    Out = rmfield(Out, 'fileUUID');
                catch
                end
                testCase.verifyEqual(Exp, Out);
            end
        end
%         function imaginarySolution(testCase)
%             actSolution = quadraticSolver(1,2,10);
%             expSolution = [-1+3i -1-3i];
%             testCase.verifyEqual(actSolution,expSolution)
%         end
    end
    
end