classdef Protocol < handle
    % This class creates and manages a list of "Subject" objects.
    %   This class uses an user-defined function (ProtoFunc) to access the
    %   acquisition information from each subject and assign it to the
    %   pertinent objects. For more information, see documentation of the
    %   following classes: Subject, Acquisition and children of the abstract
    %   class "Modality".
    
    properties 
        MainDir % Folder containing all experiment files.
        SaveDir % Folder to save "Protocol" object and HDF5 files. Default value = current folder.
        ProtoFunc % Function handle of the user-defined OpenProtocol
        % where the Subjects and Acquisition data are created.
        Array % List of Subjects. Default: empty ObjectListManager.
    end
    properties (Access = private)
        RecordingStruct % Recording Structure obtained from ProtoFunc.
    end
    
    methods
        
        function obj = Protocol(MainDir, SaveDir, ProtoFunc, Array)
            % Class constructor.
            %   This function initiates the object "Protocol" with the
            %   properties: MainDir, SaveDir, ProtoFunc and Array. 
            %   All first inputs must be provided. If Array is empty,
            %   the function creates an emtpy Array.
            if nargin > 0
                obj.MainDir = MainDir;
                obj.SaveDir = SaveDir;
                obj.ProtoFunc = ProtoFunc;
                obj.Array = Array;
            else
                obj.Array = ObjectListManager();
            end
        end
        
        %%% Property Set Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set.MainDir(obj, MainDir)
            % Set function for MainDir property.
            %   Accepts only existing Folders as input.
            mustBeFolder(MainDir); % Checks for existing Path..
            obj.MainDir = checkFolder(MainDir);
        end
        
        function set.SaveDir(obj, SaveDir)
            % Set function for SaveDir property.
            %   Accepts only existing Folders as input.
            mustBeFolder(SaveDir); % Checks for existing Path..
            obj.SaveDir = checkFolder(SaveDir);
        end
        
        function set.ProtoFunc(obj, ProtoFunc)
            % Set function for ProtoFunc property.
            %   Accepts only valid function handles.
            mustBeA(ProtoFunc, 'function_handle');
            obj.ProtoFunc = ProtoFunc; % Checks if ProtoFunc is a function handle.
        end
        
        function set.Array(obj, Array)
            % Set function for Array property.
            %   Accepts only a "ObjectListManager" object as input. If
            %   empty, creates an default "ObjectListManager" object.
            if ~isempty(Array)
                mustBeA(Array, 'ObjectListManager'); % Checks if Array is an "ObjectListManager".
                obj.Array = Array;
            else
                obj.Array = ObjectListManager();
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function generateList(obj)
           % This function uses the ProtoFunc to create the lists of
           % Subjects and Acquisitions.
           SubjArray = obj.ProtoFunc(obj);
           obj.Array.addObj(SubjArray);
           disp('List Generated')
        end
        
        function updateList(obj)
            newArray = obj.ProtoFunc(obj);
            % To be continued.
            
        end
        function resetList(obj)
            obj.Array = ObjectListManager;
        end
               
        function delete(obj)
%             for i = 1:length(obj.Array.ObjList)
%                 obj.Array.ObjList.delete
%             end
            disp('Protocol deleted')
        end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%% Extra Methods for Protocol class %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    end
    
end