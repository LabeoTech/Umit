classdef Subject < handle
    % This class creates and manages a list of "Acquisitions" of a Subject.
    %   This class instantiates an object containing information about a
    %   subject (ex. a mouse) in an experiment. It contains the subject's
    %   ID, GroupID (name of experimental group) and an "ObjectListManager"
    %   object containing an array of "Acquisition" objects.
    
    properties 
        ID % Subject ID (It has it's on set function).
        Sex char {mustBeNonempty} = 'Unknown'% Animal's Gender ID        
        Strain char {mustBeNonempty} = 'Unknown' % Animal's strain (Ex. Mouse strain "C57BL/6J").
        Calcium_indicator char {mustBeNonempty} = 'Unknown'% Name of the calcium indicator.        
    end
    properties (SetAccess = {?Protocol, ?PipelineManager, ?ObjectListManager})
        Array % List of Acquisitions.
        MyParent % Protocol Object.
        GroupID % Experimental Group of Subject.
        LastLog = table.empty % Table containing information about the Last Pipeline Operations run by PIPELINEMANAGER.
    end
    properties (Dependent, SetAccess = private)
        SaveFolder % Path of directory containing transformed data.        
    end
    
    methods
        
        function obj = Subject(ID, GroupID, Array)
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
                obj.Array = ObjectListManager([],obj);
            end
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
            if obj.MyParent.b_isDummy
                out = obj.MyParent.SaveDir;
                return
            end
            out = fullfile(obj.MyParent.SaveDir, obj.ID);
            msgID = 'umIToolbox:FolderNotFound';
            msg = 'Subject SaveFolder doesnt exist.';
            assert(isfolder(out), msgID,msg);
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

