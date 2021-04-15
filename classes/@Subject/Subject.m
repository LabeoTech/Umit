classdef Subject < handle
    % This class creates and manages a list of "Acquisitions" of a Subject.
    %   This class instantiates an object containing information about a
    %   subject (ex. a mouse) in an experiment. It contains the subject's
    %   ID, GroupID (name of experimental group) and an "ObjectListManager"
    %   object containing an array of "Acquisition" objects.
    
    properties
        ID % Subject ID
        Calcium_indicator % Name of the calcium indicator.
        MyParent % Protocol Object.
    end
    properties (SetAccess = {?Protocol, ?PipelineManager, ?ObjectListManager})
        Array % List of Acquisitions.
        GroupID % Experimental Group of Subject.
        LastLog % MAT file with a table containing information about the Last Pipeline Operations run by PIPELINEMANAGER.
    end
    properties (Dependent)
        SaveFolder % Path of directory containing transformed data.
        FilePtr % JSON file containing information of files created using PIPELINEMANAGER.
    end
    methods
        
        function obj = Subject(ID, GroupID, Array, Calcium_indicator)
            % Class Constructor.
            %   This function initiates the Subject class with the
            %   animal's basic information and an Array. Subject's ID
            %   must be provided. If GroupID is empty, a default name "def"
            %   is created. If Array is empty, an empty "ObjectListManager"
            %   is created.
            if nargin > 0
                obj.ID = ID;
                obj.GroupID = GroupID;
                obj.Array = Array;
                obj.Calcium_indicator = Calcium_indicator;
            else
                obj.ID = 'def';
                obj.GroupID = 'def';
                obj.Array = [];
            end
        end
        %%% Property Set Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set.ID(obj, ID)
            % Set function for ID property.
            validateattributes(ID, {'char', 'string'}, {'nonempty'}); % Validates if ID is a non-empty text.
            obj.ID = ID;
        end
        
        function set.GroupID(obj, GroupID)
            % Set function for GroupID property.
            %   Accepts only non-empty strings as input.
            validateattributes(GroupID, {'char', 'string'}, {'scalartext'}); % Validates if GroupID is a non-empty text.
            obj.GroupID = GroupID;
        end
        
        function set.Array(obj, Array)
            % Set function for Array property.
            %   Accepts only a "ObjectListManager" object as input. If
            %   empty, creates an default "ObjectListManager" object.
            
            if isa(Array, 'Acquisition')
                obj.Array.addObj(Array);
            elseif isa(Array, 'ObjectListManager')
                obj.Array = Array;
            else
                obj.Array = ObjectListManager();
            end
        end
        
%         function set.SaveFolder(obj, SaveFolder)
%             % Set function for SaveFolder property.
%             obj.SaveFolder = checkFolder(SaveFolder);
%         end
        
        function set.Calcium_indicator(obj, GECI)
            % Set function for CALCIUM_INDICATOR property.
            validateattributes(GECI, {'char', 'string'}, {'scalartext'}); % Validates if GECI is a non-empty text.
            obj.Calcium_indicator = GECI;
        end
        function set.MyParent(obj, Parent)
            % Set function for MYPARENT property.
            msgID = 'UMIToolbox:InvalidInput';
            msg = 'Parent of Subject object must be a Protocol.';
            assert(isa(Parent,'Protocol'), msgID, msg)
            obj.MyParent = Parent;
        end
        %%% Property Get functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get.SaveFolder(obj)
            % Get function for depentend property SaveFolder.
            out = fullfile(obj.MyParent.SaveDir, obj.ID);
            msgID = 'UMIToolbox:FolderNotFound';
            msg = 'Subject SaveFolder doesnt exist.';
            assert(isfolder(out), msgID,msg);
        end
        function out = get.FilePtr(obj)
            % Get function for depentend property FilePtr.
            out = fullfile(obj.MyParent.SaveDir, obj.ID, 'FilePtr.json');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function createFilePtr(obj)
            % This function creates a JSON file containing basic information from object.
            
            % FilePtr full path:
%             obj.FilePtr = fullfile(obj.SaveFolder, 'FilePtr.json');
            if isfile(obj.FilePtr)
                disp(['Skipped FilePtr creation. File pointer already exists in ' obj.SaveFolder ]);
                return
            end
            A = struct('Type', class(obj), 'ID', obj.ID, 'Files', []);
            txt = jsonencode(A);
            fid = fopen(obj.FilePtr, 'w');
            fprintf(fid, '%s', txt);
            fclose(fid);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function dummyMethodForTesting(obj)
            disp(class(obj))
            disp(['This is a dummy function for testing of SUBJECT ' obj.ID '!'])
        end
    end
    methods (Access = private)
        function delete(obj)
            %             for i = 1:length(obj.Array.ObjList)
            %                 obj.Array.ObjList.delete
            %             end
            %             disp('Subject deleted')
        end
    end
    
end

