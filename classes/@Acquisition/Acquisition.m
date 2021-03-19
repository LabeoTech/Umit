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
        Array % List of Modalities.
        Start_datetime % Date and time of the beginning of the acquisition. 
        
    end
    properties (SetAccess = {?Protocol, ?PipelineManager})
        SaveFolder % Path of directory containing transformed data.
        LastLog % MAT file with a table containing information about the Last Pipeline Operations run by PIPELINEMANAGER.
        FilePtr % JSON file containing information of files created using PIPELINEMANAGER.
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
                obj.Array = ObjectListManager();
            end
        end
        
        function set.SaveFolder(obj, SaveFolder)
            % Set function for SAVEFOLDER property and creates LogBook inside SAVEFOLDER.
            obj.SaveFolder = checkFolder(SaveFolder);
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
                
        %%% Property get functions
        function out = get.Start_datetime(obj)
            % This function transforms datetime objects in strings for
            % display in the GUI.
            
            if isempty(obj.Start_datetime) || isnat(obj.Start_datetime)
                out = '';
            else
                out = datestr(obj.Start_datetime, 31);
            end
        end
        
         function createFilePtr(obj)
             % This function creates a JSON file containing basic information from object.
             
             % FilePtr full path:
             obj.FilePtr = fullfile(obj.SaveFolder, 'FilePtr.json');
             A = struct('Type', class(obj), 'ID', obj.ID, 'Files', []);
             txt = jsonencode(A);
             fid = fopen(obj.FilePtr, 'w');
             fprintf(fid, '%s', txt);
             fclose(fid);
         end
        
        %%%%%%%%%%
        function dummyMethodForTesting(obj)
            disp(class(obj))
            disp(['This is a dummy function for testing of Acquisition ' obj.ID '!'])
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