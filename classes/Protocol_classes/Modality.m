classdef (Abstract) Modality < matlab.mixin.Heterogeneous & handle
    % Modality Class.
    % This is a "superclass" for the classes that will handle the
    % different Data types (modalities) from each Acquisition (for info about
    %   Acquisition objects see documentation on the "Acquisition" class.
    
    properties
        ID % Unique Identifier of the object.
        RecordingSystem char {mustBeNonempty, isscalar} = 'default' % Name of the system used to record the data.
        SampleRateHz double {mustBeNonempty} = 0.0 % Sampling rate of the recording in Hz.        
    end
    properties (SetAccess = {?Protocol, ?PipelineManager, ?Acquisition, ?Subject, ?ObjectListManager})
        LastLog = table.empty % Table containing information about the Last Pipeline Operations run by PIPELINEMANAGER.
        MyParent % Acquisition object that contains MODALITY object.
    end
    properties (Hidden)
        MetaDataFileName char = '' %File containing other information about the recording session. (OBSOLETE. Kept for retrocompatibility).
        RawFiles_FP cell % List of full path of Raw filenames. This property is set/used only by the Protocol Function.
                         % This property will be used to generate just the
                         % list of files (no path) and to create the
                         % "RawFolder" property.
    end
    properties (Dependent, SetAccess = private)
        RawFiles cell % File(s) containing raw data.
        RawFolder char % Path of directory containing raw data.    
        SaveFolder % Path of directory containing transformed data.        
        MetaDataFile % FullPath of file containing other information about the recording session. (OBSOLETE. Kept for retrocompatibility).
    end
    
    methods
        
        function obj = Modality(ID, RawFiles_FP, RecordingSystem, SampleRate)
            % Construct an instance of this class.
            %   The Folder and FileName are defined here.
            %   Folder must be a valid Directory while FileName has to be a
            %   string or a cell array of strings containing valid filenames
            %   that exist in the "Folder" directory.
            if nargin > 0
                obj.ID = ID;                
                obj.RawFiles_FP = RawFiles_FP;
                obj.RecordingSystem = RecordingSystem;
                obj.SampleRateHz = SampleRate;
            else
                obj.ID = 'def';
            end
        end        
        %%% Property Set Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set.ID(obj, ID)
            % Set function for ID property.
            %   Accepts only non-empty strings.
            validateattributes(ID, {'char', 'string'}, {'nonempty'}); % Validates if ID is a non-empty text.
            obj.ID = ID;
        end
                     
        function set.RawFiles_FP(obj,RawFiles)
            % Set function of RAWFILES property.
            %   Validates if Files exist in Folder, then sets the RAWFILES
            %   property, otherwise throws an error. Duplicate file names
            %   are ignored.
            % Input:
            % RawFiles (cell): cell array of strings
            % containing the full path of Raw Files.

            % Control for more than one Raw Folder. All Raw files must be
            % in the same folder:
            [FolderList,~,~] = cellfun(@(x) fileparts(x), RawFiles, 'UniformOutput', false);            
            if numel(unique(FolderList)) > 1
                errID = 'umIToolbox:Modality:TooManyInputs';                
                errMsg = 'Cannot Set RawFiles Property. All files must be in the same folder.';
                error(errID, errMsg);
            end
          
            % Removes duplicates.
            [~,idx] = ismember(unique(RawFiles), RawFiles); % Keeps the first of the list;
            RawFiles = RawFiles(idx);
            % Removes non-existent file names from the list "FileName"
            b_exist = isfile(RawFiles);
            if ( sum(b_exist) < length(b_exist) && sum(b_exist) > 0 )
                disp('The following files were not found and will be ignored:')
                disp(RawFiles(~b_exist)')
                RawFiles(~b_exist)= [];
            elseif ( sum(b_exist) == 0 )
                % If no files exist, Keep only the folder of the first file
                % So the "RawFolder" property can be set.
                RawFiles = {fullfile(fileparts(RawFiles{1}), 'MISSING_FILES')};
            end
            % Set RawFiles property:
            obj.RawFiles_FP = RawFiles;
        end
        
        function set.RecordingSystem(obj, RecordingSystem)
            % Set function of RecordingSystem property.
            %   Validates if RECORDINGSYSTEM is a string.
            validateattributes(RecordingSystem, {'char', 'string'}, {'scalartext'}); % Validates if RecordingSystem is text.
            obj.RecordingSystem = RecordingSystem;
        end
        
        function set.MyParent(obj, MyParent)
            % Set function for MyParent property.
            msgID = 'umIToolbox:Modality:InvalidInput';
            msg = 'Error setting Modality Parent Object. Parent must be Acquitision.';
            assert(isa(MyParent, 'Acquisition'),msgID, msg);
            obj.MyParent = MyParent;
        end
        %%% Validators %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function validate_path(~, input)            
            % This method checks if the input is an existing folder in
            % Matlab's path.
            % Input
            %   input (char): String containing fullpath to a folder.
            errID = 'umIToolbox:Modality:InvalidInput';
            msg = ['The folder "' input '" does not exist or isn''t in Matlab''s Path.'];
            assert(isfolder(input), errID, msg);
        end
        %%% Property Get functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function out = get.RawFolder(obj)
            % Get function of RAWFOLDER property.
            if obj.MyParent.MyParent.MyParent.b_isDummy
                out = obj.MyParent.MyParent.MyParent.MainDir;
                return
            end
            [out,~,~]  = fileparts(obj.RawFiles_FP{1});
        end
        
        function out = get.RawFiles(obj)
            % Get function of RAWFILES property.
            
            [~,out,~]  = cellfun(@(x) fileparts(x), obj.RawFiles_FP, 'UniformOutput',false);
        end
        
        function out = get.SaveFolder(obj)
            % Get function for depentend property SaveFolder.
            if obj.MyParent.MyParent.MyParent.b_isDummy
                out = obj.MyParent.MyParent.MyParent.SaveDir;
                return
            end
            out = fullfile(obj.MyParent.MyParent.MyParent.SaveDir, obj.MyParent.MyParent.ID, obj.MyParent.ID, obj.ID);
            msgID = 'UMIToolbox:FolderNotFound';
            msg = 'Modality SaveFolder doesnt exist.';
            assert(isfolder(out), msgID,msg);
        end
        
        function out = get.MetaDataFile(obj)
            % Get function for depentend property MetaDataFile.
            out = fullfile(obj.RawFolder, obj.MetaDataFileName);
            if isempty(obj.MetaDataFileName)
                out = 'none';
                return
            elseif ~isfile(out)
                out = 'missing';            
            end
        end
    end
    methods (Access = protected)
        
        function delete(obj)
            %             disp(['Modality of type ' class(obj) ' deleted'])
        end             
    end
end