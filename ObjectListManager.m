classdef ObjectListManager < handle
    % This is a class created to manage lists of objects.
    %   It creates and manages a vector of objects for use by the classes:
    %   Protocol, Subject and Acquisition.
    %   ObjectListManager only accepts objects with unique "ID" property.
    
    properties
        ObjList(1,:) {mustBeVector} % Array of objects.
    end
    
    methods
        function obj = ObjectListManager(objectList)
            % Class Constructor.
            %   This function initiates the ObjectListManager class by
            %   creating a vector of objects. It checks if the input
            %   (objectList) is an object or an array of objects from
            %   the following classes: Subject, Acquisition or Modality.
            %   If no input argument is provided, it creates an empty
            %   ObjList.
            
            if nargin > 0
                checkIfIsObj(objectList);
                objectList = checkForDuplicates(objectList);
                obj.ObjList = objectList;
            end
        end
        
        function addObj(obj, newObj)
            % Adds a new object to the ObjList list.
            %   This function will accept only valid objects as
            %   input. Otherwise, it will throw an error.
            
            if size(newObj,1) > size(newObj,2)
                % transpose Object list if not organized in columns
                newObj = newObj';
            end
            checkIfIsObj(newObj); % Checks for valid Object types
            newObj = checkForDuplicates(newObj); % Filter out duplicate objects (with the same ID) in newObj
            if ~isempty(obj.ObjList)
                newObj = compareLists(obj, newObj);
            end
            obj.ObjList = [obj.ObjList newObj];
        end
        
        function removeObj(obj,idx)
            % This function removes a specific object from the obj.ObjList.
            %   The input idx is the index of the object(s) in the
            %   obj.ObjList. The index "idx" can be obtained using the
            %   findElement function.
            % To be improved.
            obj.ObjList(idx) = [];
            disp(['Object(s) with index(es) ' num2str(idx) ' successfully removed']);
        end
        
        function indElem = findElement(obj, varargin)
            % This function looks for an specific object inside
            % ObjList.
            %   INDELEM = findElement(CRITERIA#1, PATTERN#1, CRITERIA#2,
            %   PATTERN(s)#2, ...) returns the index(es) of the ObjectList
            %   that have 'criteria#1' AND 'criteria#2' AND ... .           
            %   "Criteria" must be a string with the name of a object's
            %   property. "Pattern" has to be from the data type of the
            %   queried property.
                        
            % Controls for empty variables or odd number of variables.
            if ( isempty(varargin)) || ( mod(length(varargin),2)~=0 ) 
                error('InsufficientInputs')
            end
            
            criteria = varargin(1:2:end); % Keep it as cell array.
            pattern = varargin(2:2:end); 
            propNames = properties(obj.ObjList);
            idx = cellfun(@(x) any(strcmp(propNames,x)), criteria,...
                'UniformOutput', false); idx = [idx{:}];
            if any(~idx)
                disp('The following criteria were not found:')
                disp(criteria(~idx));
            end
            criteria = criteria(idx); % Eliminate criteria inputs not found in the ObjectList.
            pattern = pattern(idx);
            idx = [];
            for i = 1:length(criteria)
                if ischar(pattern{i})
                eval(['idx(i,:) = strcmp({obj.ObjList.' criteria{i} '}, ''' pattern{i} ''');'])
                else
                    % Just in case if the variable is not a string.
                eval(['idx(i,:) = cellfun(@(x) isequal(x, ' pattern{i} '), {obj.ObjList.' criteria{i} '});']);
                end
            end
            indElem = find(all(idx,1)); % Must respect ALL patterns.
            
        end
        
        function eraseObjArray(obj)
            % This function erases the ObjArray property.
            obj.ObjList = [];
            disp('ObjArray successfully erased')
        end
        
        function delete(obj) % To be rechecked.
            
%             for i = 1:length(obj.ObjList)
%                 delete(obj.ObjList(i))
%                 disp(['# ' num2str(i) ' -- Object of type ' class(obj.ObjList(i)) ' deleted from the list'])
%             end
        end
        
        
        
    end
    
end


%%% Local functions %%%

function checkIfIsObj(objectList)
% This function checks if the input (objectList) is an object
% or an array of objects.
%   The function accepts only the objects instantiated from the
%   following classes (or its children): Subject, Acquisition and Modality.

if isa(objectList, 'Subject') || isa(objectList, 'Acquisition')...
        || isa(objectList, 'Modality')
    return;
else
    disp('---------------------------------------------------------------')
    disp(['Object Data Type provided: ' class(objectList)]);
    disp('Please provide valid Data type');
    disp('Valid data types are classes : Subject, Acqusition or Modality')
    disp('---------------------------------------------------------------')
    error('Could not create object Array. Data Type not accepted.')
end
end

function objectList = checkForDuplicates(objectList)
% This function remove duplicate members of objectList.
%   The function looks for unique "ID" property of the object
%   and remove duplicates with the same "ID".

IDlist = {objectList.ID};
uniqList = unique(IDlist);
if numel(uniqList) < numel(IDlist)
    [~,idx] = ismember(uniqList, IDlist); % Keeps the first of the list;
    objectList = objectList(idx);
    disp('More than one object with the same "ID" were found in the input list.');
    disp('Duplicates were ignored.');
end
end

function objectList = compareLists(obj, objectList)
% This function checks if the new object(s) in "objectList" is already
% part of obj.ObjList.
%   The function inputs are an object "ObjectListManager" and a vector
%   of valid objects (see function "checkIfIsObj"). It uses the "ID"
%   property of objects to check for duplicates.

isDup = ismember({objectList.ID}, {obj.ObjList.ID});
if any(isDup)
    objectList = objectList(~isDup);
    disp('The following objects already exist on the ObjList and were ignored.')
    loc = find(isDup); 
    for i = 1:numel(loc)
        disp(obj.ObjList(loc(i))) 
    end
end

end

