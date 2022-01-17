classdef ObjectListManager < handle
    % This is a class created to manage lists of objects.
    %   It creates and manages a vector of objects for use by the classes:
    %   Protocol, Subject and Acquisition.
    %   ObjectListManager only accepts objects with unique "ID" property.
    
    properties
        ObjList(1,:) % Array of objects.
    end
    
    properties (Access = private)
        parentObj % Parent object that contains an ObjectListManager object.
    end
    
    methods
        function obj = ObjectListManager(objectList, parentObj)
            % Class Constructor.
            %   This function initiates the ObjectListManager class by
            %   creating a vector of objects. It checks if the input
            %   (objectList) is an object or an array of objects from
            %   the following classes: Subject, Acquisition or Modality.
                        
            if ~isempty(objectList)
                checkIfIsObj(objectList);
                objectList = checkForDuplicates(objectList);
                obj.ObjList = objectList;
            end
            obj.parentObj = parentObj;
        end 
        function addObj(obj, newObj)
            % Adds a new object to the ObjList list.
            %   This function will accept only valid objects as
            %   input. Otherwise, it will throw an error.
            % Inputs:
            % newObj : scalar or array of valid objects ("Subject",
            % "Acquisition" or "Modality").
            if size(newObj,1) > size(newObj,2)
                % transpose Object list if not organized in columns
                newObj = newObj';
            end
            checkIfIsObj(newObj); % Checks for valid Object types
            newObj = checkForDuplicates(newObj); % Filter out duplicate objects (with the same ID) in newObj
            if ~isempty(obj.ObjList)
                newObj = compareLists(obj, newObj);
            end
            % Set "MyParent" property of newObj:
            for i = 1:length(newObj)
                newObj(i).MyParent = obj.parentObj;
            end 
            obj.ObjList = [obj.ObjList newObj];
            if ~isempty(newObj) && ( isa(newObj, 'Subject') || isa(newObj, 'Acquisition')) 
                % This if statement is to control for objects added
                % manually (one at a time).
                if isempty(newObj(1).Array.ObjList)
                    uiwait(warndlg(['Object of type ' class(newObj) ' is empty. You have to include at least one element inside it before continuing the analysis.'], 'Warning!', 'modal'));
                end
            end
        end
        
        function out = removeObj(obj,idx)
            % This function removes a specific object from the obj.ObjList.
            %   The input idx is the index of the object(s) in the
            %   obj.ObjList. The index "idx" can be obtained using the
            %   findElement function. The function outputs a string array
            %   containing the class and ID of the deleted objects.
            out = arrayfun(@(x) [string(class(x)) string(x.ID) string(x.SaveFolder)], obj.ObjList(idx), 'UniformOutput', false);
            out = [out{:}];
            out = reshape(out, [], length(idx)); out = out';
            obj.ObjList(idx) = [];
            disp(['Object(s) with index(es) ' num2str(idx) ' successfully removed']);
        end
               
        function out = listProp(obj, propName)%#ok
            % LISTPROP lists all values of PROPNAME from objects inside OBJ.OBJLIST.
            idx = arrayfun(@(x) ismember(propName, properties(x)), obj.ObjList);%#ok
            eval(['out = {obj.ObjList(idx).' propName '};']);
        end
                
        function iElem = findElement(obj, PropName, str, varargin)
            % This function looks for an specific object inside
            % ObjList.
            %   INDELEM = findElement(PROPERTY_NAME, STRING, QUERYMETHOD)
            %   returns the index(es) of the ObjectList with 'PROPERTY_NAME'
            %   equals to the STRING. Output is the index of the element in OBJ.OBJLIST
            %   Three query methods are accepted: 'strcmp', 'contains' and
            %   'regexp'.
            p = inputParser;
            validateMethod = @(x) mustBeMember(x, {'contains', 'regexp', 'strcmp'});
            addOptional(p, 'queryMethod', 'contains', validateMethod);
            parse(p,varargin{:})
            queryMethod = p.Results.queryMethod;
            
            if isempty(PropName)
                iElem = 1:length(obj.ObjList);
            else
                
                switch queryMethod
                    case 'strcmp'
                        eval(['iElem = find(strcmp(''' str ''', {obj.ObjList.' PropName '}));'])
                    case 'contains'
                        eval(['iElem = find(contains({obj.ObjList.' PropName '}, ''' str '''));'])
                    case 'regexp'
                        eval(['iElem = regexp({obj.ObjList.' PropName '}, ''' str ''');'])
                        iElem = cellfun(@(x) ~isempty(x), iElem); %#ok
                        iElem = find(iElem);
                end
            end
        end
    end
    methods (Access = private)
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

