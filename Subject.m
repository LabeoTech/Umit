classdef Subject < handle
    % This class creates and manages a list of "Acquisitions" of a Subject.
    %   This class instantiates an object containing information about a
    %   subject (ex. a mouse) in an experiment. It contains the subject's
    %   ID, GroupID (name of experimental group) and an "ObjectListManager"
    %   object containing an array of "Acquisition" objects.
    
    properties 
        ID % Subject ID
        GroupID % Experimental Group of Subject.
        Array % List of Acquisitions. 
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
                obj.Array = ObjectListManager();
            end
        end
        %%% Property Set Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set.ID(obj, ID)
            % Set function for ID property.
            mustBeNonzeroLengthText(ID); % Checks if string is empty.
            obj.ID = ID;
        end
        
        function set.GroupID(obj, GroupID)
            % Set function for GroupID property.
            %   Accepts only non-empty strings as input.
            mustBeNonzeroLengthText(GroupID); % Checks if string is empty.
            obj.GroupID = GroupID;
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
       
        function delete(obj)
%             for i = 1:length(obj.Array.ObjList)
%                 obj.Array.ObjList.delete
%             end
            disp('Subject deleted')
        end
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                         %%% Extra Methods for Subject class %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        
        
        
        
    end
end

