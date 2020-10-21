classdef Acquisition < handle
    % This class creates and manages a list of one or more recordings from
    % during an acquisition session.
    %   This class instantiates an object containing information about a
    %   modalities (ex. Brain Imaging, Eye Tracking) in a recording session.
    %   It contains the acquisiion's identifier (ID)an "ObjectListManager"
    %   object containing an array of objects instantiated from a child 
    %   class of "Modality".
    
    properties 
        ID % Acquisition ID
        Array % List of Modalities.
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
                obj.Array = ObjectListManager();
            end
        end
        
         %%% Property Set Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set.ID(obj, ID)
            % Set function for ID property.
            mustBeNonzeroLengthText(ID); % Checks if string is empty.
            obj.ID = ID;
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
%             disp('Acquisition deleted')
        end
        
        
        
        
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     %%% Extra Methods for Acquisition class %%%
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        
        
        
        
    end
end

