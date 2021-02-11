classdef (Abstract) Modality < matlab.mixin.Heterogeneous & handle
    % Modality Class.
    % This is a "superclass" for the classes that will handle the
    % different Data types (modalities) from each Acquisition (for info about
    %   Acquisition objects see documentation on the "Acquisition" class.
    
    properties
        ID % Unique Identifier of the object.
        RawFolder % Path of directory containing raw data.
        RawFiles % File(s) containing raw data.
        RecordingSystem % Name of the system used to record the data.
        SampleRateHz % Sampling rate of the recording in Hz.
    end
    properties (SetAccess = {?Protocol})
        SaveFolder % Path of directory containing transformed data.
    end
    properties (SetAccess = {?PipelineManager})
        LastLog % MAT file with a table containing information about the Last Pipeline Operations run by PIPELINEMANAGER.
    end
    methods
        
        function obj = Modality(ID, RawFolder, RawFiles, RecordingSystem, SampleRate, ~)
            % Construct an instance of this class.
            %   The Folder and FileName are defined here.
            %   Folder must be a valid Directory while FileName has to be a
            %   string or a cell array of strings containing valid filenames
            %   that exist in the "Folder" directory.
            if nargin > 0
                obj.ID = ID;
                obj.RawFolder = RawFolder;
                obj.RawFiles = RawFiles;
                obj.RecordingSystem = RecordingSystem;
                obj.SampleRate = SampleRate;
            else
                obj.ID = 'def';
            end            
        end
        
        %%% Property Set Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set.ID(obj, ID)
            % Set function for ID property.
            %   Accepts only non-empty strings.
            mustBeNonzeroLengthText(ID); % Checks if string is empty.
            obj.ID = ID;
        end
        
        function set.RawFolder(obj, RawFolder)
            % Set function of RAWFOLDER property.
            %   This function accepts only existing folders.
            mustBeFolder(RawFolder);
            obj.RawFolder = checkFolder(RawFolder);
        end
        
        function set.RawFiles(obj,RawFiles)
            % Set function of RAWFILES property.
            %   Validates if Files exist in Folder, then sets the RAWFILES
            %   property, otherwise throws an error. Duplicate file names
            %   are ignored.
            
            obj.RawFiles = validateFileName(obj.RawFolder, RawFiles);
        end
        function set.RecordingSystem(obj, RecordingSystem)
            % Set function of RecordingSystem property.
            %   Validates if RECORDINGSYSTEM is a string.
            mustBeText(RecordingSystem)
            obj.RecordingSystem = RecordingSystem;
        end
        
        function set.SampleRateHz(obj, SampleRate)
            % Set function of SampleRate property( in Hz ).
            %   Transforms SAMPLERATE to integer.
            obj.SampleRateHz = int32(SampleRate); % Rounds-up to integer.
        end
        
         function set.SaveFolder(obj, SaveFolder)
            % Set function of SAVEFOLDER property.
            %   This function accepts only existing folders.
            obj.SaveFolder = checkFolder(SaveFolder);
        end
       
    end
    
    methods (Access = protected)
        function delete(obj)
            %             disp(['Modality of type ' class(obj) ' deleted'])
        end
    end
end


% Local Functions
function FileName = validateFileName(Folder, FileName)
% This functions validates if the files in FILENAME exist in FOLDER.
% File(s) not found are removed from the list.

if iscell(FileName)
    % Removes duplicates.
    [~,idx] = ismember(unique(FileName), FileName); % Keeps the first of the list;
    FileName = FileName(idx);
    % Removes non-existant file names from the list "FileName"
    tmp = cellfun(@(x) isfile([Folder x]), FileName, 'UniformOutput', false);
    tmp = cell2mat(tmp);
    if ~all(tmp)
        disp(['The following files were not found in ' Folder ' and were ignored.'])
        disp(FileName(~tmp))
        FileName(~tmp)= [];
    end
elseif ischar(FileName)
    if isfile([Folder FileName])
        return
    else
        disp(['"' FileName '" was not found in ' Folder ' and was ignored']);
        FileName = [];
    end
else
    error('Wrong Data type. FileName must be a String or a cell array containing strings.');
end
end