classdef Acquisition < handle
    % This class creates and manages a list of one or more recordings from
    % during an acquisition session.
    %   This class instantiates an object containing information about a
    %   modalities (ex. Brain Imaging, Eye Tracking) in a recording session.
    %   It contains the acquisiion's identifier (ID) an "ObjectListManager"
    %   object containing an array of objects instantiated from a child
    %   class of "Modality".
    
    properties
        ID % Acquisition ID
        Start_datetime % Date and time of the beginning of the acquisition.
    end
    properties (SetAccess = {?Protocol, ?PipelineManager, ?Subject, ?ObjectListManager})
        Array % List of Modalities.
        MyParent % Subject object that contains ACQUISITION object.
        LastLog % MAT file with a table containing information about the Last Pipeline Operations run by PIPELINEMANAGER.
    end
    properties (Dependent)
        SaveFolder % Path of directory containing transformed data.
    end
    
    methods
        
        function obj = Acquisition(ID, Array)
            % Class Constructor.
            %   This function initiates the Acquisition class with the
            %   Acquisition ID an Array. Acquisition's ID
            %   must be provided. If Array is empty, an empty "ObjectListManager"
            %   is created.
            
            if nargin > 0
                obj.ID = ID;
                obj.Array = Array;
            else
                obj.ID = 'def';
                obj.Array = [];
            end
        end
        %%% Property Set Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set.ID(obj, ID)
            % Set function for ID property.
            validateattributes(ID, {'char', 'string'}, {'nonempty'}); % Validates if ID is a non-empty text.
            obj.ID = ID;
        end
        
        function set.Array(obj, Array)
            % Set function for Array property.
            %   Accepts only a "ObjectListManager" object as input. If
            %   empty, creates an default "ObjectListManager" object.
            if isa(Array, 'Modality')
                obj.Array.addObj(Array);
            elseif isa(Array, 'ObjectListManager')
                obj.Array = Array;
            else
                obj.Array = ObjectListManager([],obj);
            end
        end
        
        function set.SaveFolder(obj, SaveFolder)
            % Set function for SAVEFOLDER property and creates LogBook inside SAVEFOLDER.
            obj.SaveFolder = checkFolder(SaveFolder, 'new');
        end
        
        function set.Start_datetime(obj, start_datetime)
            % This function creates a DATETIME object based on Datetime
            % String input from User with the format 'yyyy-MM-dd HH:mm:ss'
            infmt = 'yyyy-MM-dd HH:mm:ss';
            if ~isa(start_datetime, 'datetime') && ~isempty(start_datetime)
                try
                    obj.Start_datetime = datetime(start_datetime, 'InputFormat', infmt);
                catch
                    error(['Wrong DateTime format. It has to be like : ' infmt])
                end
            elseif isempty(start_datetime)
                obj.Start_datetime = [];
            elseif isnat(start_datetime)
                obj.Start_datetime = [];
            else
                validateattributes(start_datetime, {'datetime'}, {'nonempty'});
                obj.Start_datetime = start_datetime;
            end
        end
        
        function set.MyParent(obj, MyParent)
            % Set function for MyParent property.
            msgID = 'umIToolbox:Acquisition:InvalidInput';
            msg = 'Error setting Acquisition Parent Object. Parent must be a Subject.';
            assert(isa(MyParent, 'Subject'), msgID, msg);
            obj.MyParent = MyParent;
        end
        %%% Property Get functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get.SaveFolder(obj)
            % Get function for depentend property SaveFolder.
            out = fullfile(obj.MyParent.MyParent.SaveDir,obj.MyParent.ID, obj.ID);
            msgID = 'umIToolbox:Acquisition:FolderNotFound';
            msg = 'Acquisition SaveFolder doesnt exist.';
            assert(isfolder(out), msgID,msg);
        end
        
        function out = get.Start_datetime(obj)
            % This function transforms datetime objects in strings for
            % display in the GUI.
            if isempty(obj.Start_datetime) || isnat(obj.Start_datetime)
                out = '';
            else
                out = datestr(obj.Start_datetime, 31);
            end
        end
    end
    
    methods (Access = private)
        function delete(obj)
            %             for i = 1:length(obj.Array.ObjList)
            %                 obj.Array.ObjList.delete
            %             end
            %             disp('Acquisition deleted')
        end
    end
end