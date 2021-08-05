classdef test_run_HemoCorrection < matlab.unittest.TestCase
    properties (TestParameter)
        rawFolder = {'\test\test_runImagesClassification\exp_outFiles\'};
        expOutFolder = {'\test\test_run_HemoCorrection\exp_outFiles\'};
        outFolder = {'\test\test_run_HemoCorrection\outputFolder\'};
        opts = {struct('Red', true, 'Green', true, 'Amber', true)};
    end
    
    methods(Test)
        function realSolution(testCase, rawFolder, outFolder, expOutFolder, opts)
            root = getenv('Umitoolbox');
            
            % Since the function HemoCorrection works with files in the
            % same folder, copy the raw data from rawFolder to
            % outputFolder:
            dat_files = dir(fullfile(root,rawFolder,'*.dat'));
            mat_files = dir(fullfile(root,rawFolder,'*_info.mat'));
            for i = 1:length(dat_files)
                copyfile(fullfile(dat_files(i).folder, dat_files(i).name), fullfile(root,...
                    outFolder,dat_files(i).name));
                copyfile(fullfile(mat_files(i).folder, mat_files(i).name), fullfile(root,...
                    outFolder,mat_files(i).name));
            end
            
            outFile = run_HemoCorrection(fullfile(root,outFolder,'fluo.dat'),...
                fullfile(root, outFolder),opts);            
            [~,out_filename,~] = fileparts(outFile);
            % Test outFile name:
            testCase.verifyEqual(out_filename, 'hemoCorr_fluo');
            % Test .DAT file:
            fid_exp = fopen(fullfile(root, expOutFolder, 'hemoCorr_fluo.dat'), 'r');
            expFile = fread(fid_exp, 'single=>single');
            fid_out = fopen(fullfile(root, outFolder, 'hemoCorr_fluo.dat'), 'r');
            outFile = fread(fid_out, 'single=>single');
            testCase.verifyEqual(expFile, outFile);
            fclose(fid_exp);
            fclose(fid_out);
            % Test .MAT files:
            Exp = load(fullfile(root, expOutFolder, 'hemoCorr_fluo_info.mat'));
            Out = load(fullfile(root, outFolder, 'hemoCorr_fluo_info.mat'));
            Exp = rmfield(Exp,{'fileUUID', 'datFile'});
            Out = rmfield(Out, {'fileUUID', 'datFile'});
            testCase.verifyEqual(Exp, Out);
            % Remove created files from outFolder:
            files = dir(fullfile(root, outFolder));
            files = files([files.isdir]==0);
            arrayfun(@(x) delete(fullfile(x.folder, x.name)),files);
            disp('Files removed from outFolder!!');
        end
    end
end