classdef (Abstract) Modality < matlab.mixin.Heterogeneous & handle
    % Modality Class.
    % This is a "superclass" for the classes that will handle the
    % different Data types (modalities) from each Acquisition (for info about
    %   Acquisition objects see documentation on the "Acquisition" class.
    
    properties
        ID % Unique Identifier of the object.
        Folder % Path of directory containing raw data.
        FileName % File(s) containing raw data.
    end
    
    methods
        
        function obj = Modality(ID, Folder, FileName)
            % Construct an instance of this class.
            %   The Folder and FileName are defined here.
            %   Folder must be a valid Directory while FileName has to be a
            %   string or a cell array of strings containing valid filenames
            %   that exist in the "Folder" directory.
            if nargin > 0
                obj.ID = ID;
                obj.Folder = Folder;
                obj.FileName = FileName;
            else
                obj.ID = 'def';
            end
        end
        
        %%% Property Set Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set.ID(obj, ID)
            % Set function for ID property.
            %   Accepts only non-empty strings.
            mustBeNonzeroLengthText(ID); % Checks if string is empty.
            obj.ID = [ID '_' class(obj)]; % Adds a '_' + class Name to the ID.
        end
        
        function set.Folder(obj, Folder)
            % Set function of Folder property.
            %   This function accepts only existing folders.
            mustBeFolder(Folder); 
            obj.Folder = checkFolder(Folder);
        end
        
        function set.FileName(obj,FileName)
            % Set function of FileName property.
            %   Validates if Files exist in Folder, then sets the FileName
            %   property, otherwise throws an error. Duplicate file names
            %   are ignored.
            
            obj.FileName = validateFileName(obj.Folder, FileName);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function delete(obj)
%             disp(['Modality of type ' class(obj) ' deleted'])
        end
        
    end
end

% Local Functions
function FileName = validateFileName(Folder, FileName)
% This functions validates if the files in "FileName" exist in "Folder".
% File(s) not found are removed from the list.

if iscell(FileName)
    % Removes duplicates.
    [~,idx] = ismember(unique(FileName), FileName); % Keeps the first of the list;
    FileName = FileName(idx);
    % Removes non-existant file names from the list "FileName"
    tmp = cellfun(@(x) isfile([Folder x]), FileName, 'UniformOutput', false);
    tmp = cell2mat(tmp);
    if ~all(tmp)
        disp(['The following files were not found in ' Folder ' and were ignored.'])
        disp(FileName(~tmp))
        FileName(~tmp)= [];
    end
elseif ischar(FileName)
    if isfile([Folder FileName])
        return
    else
        disp(['"' FileName '" was not found in ' Folder ' and was ignored']);
        FileName = [];
    end
else
    error('Wrong Data type. FileName must be a String or a cell array containing strings.');
end
end
