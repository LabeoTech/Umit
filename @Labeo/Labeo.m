classdef Labeo < Modality
    % "LabeoTech" imaging class.
    %   This class deals with the imaging data format from LabeoTech.
    %         properties
    %
    %         end
    methods
        function obj = Labeo()
            obj.LastLog.ModalityName = {'Labeo'};
        end
        
        function run_ImagesClassification(obj, FolderName, BinningSpatial, BinningTemp, b_SubROI, b_IgnoreStim)
            state = false;
            try
                ImagesClassification(FolderName, BinningSpatial, BinningTemp, b_SubROI, b_IgnoreStim);
                state = true;
                obj.LastLog.Messages = 'No Issues';
            catch ME
                obj.LastLog.Messages = {getReport(ME)};
            end
            obj.LastLog.Completed = state;
            obj.LastLog.RunDateTime = datetime('now');
        end
        
        function run_Ana_IOI_FullFrame(obj, FolderName, verbose, b_tFilter, b_fitFilter, OStream)
            state = false;
            try
                Ana_IOI_FullFrame(FolderName, verbose, b_tFilter, b_fitFilter, OStream);
                state = true;
                obj.LastLog.Messages = 'No Issues';
            catch ME
                obj.LastLog.Messages = {getReport(ME)};
            end
            obj.LastLog.Completed = state;
            obj.LastLog.RunDateTime = datetime('now');
        end
    end
end

