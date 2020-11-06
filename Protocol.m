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
            SubjArray = obj.ProtoFunc(obj.MainDir);
            obj.Array.addObj(SubjArray);
            disp('List Generated')
        end
        
        function updateList(obj, varargin)
            % This function updates the list of Subjects using
            % obj.ProtoFunc.
            %   The optional input (boolean) allows the user to not discard
            %   (FALSE) elements that were not found during the update. 
            %       *This option does not seem wise since keeping invalid 
            %       paths of filenames may cause problems later on the 
            %       analysis pipeline. "A voir..."
            %   IF discardData == FALSE, a warning message is thrown.
            
            if nargin < 2
                discardData = true;
            else
                discardData = varargin{1};
            end
            newArray = obj.ProtoFunc(obj);
            % 1st, control for new or deleted Subjects:
            iNewSubj = ~ismember({newArray.ID}, {obj.Array.ObjList.ID}); % New subjects in the newArray.
            iMissSubj = ~ismember({obj.Array.ObjList.ID},{newArray.ID});% Subjects from original list that are no longer in the newArray.
            
            % 2nd, control for new or deleted Acquisitions from each Subject: 
            indNwArr = find(~iNewSubj);
            indObj = arrayfun(@(a) find(strcmp({newArray(indNwArr).ID}, a)),...
                {obj.Array.ObjList.ID}, 'UniformOutput', false); indObj = [indObj{:}]; 
            
            iNewAcq = arrayfun(@(a,b) ~ismember({a.Array.ObjList.ID}, {b.Array.ObjList.ID}),...
                newArray(indNwArr), obj.Array.ObjList(indObj), 'UniformOutput', false);
            iMissAcq = arrayfun(@(a,b) ~ismember({b.Array.ObjList.ID}, {a.Array.ObjList.ID}),...
                newArray(indNwArr), obj.Array.ObjList(indObj), 'UniformOutput', false);
            
            %%% Updates the existing list:
            % Add new Acquisitions:
            indSubj = find(cellfun(@(x) any(x), iNewAcq));
            for i = 1:length(indSubj)
                indAcq = find(iNewAcq{indSubj(i)});
                arrayfun(@(x) obj.Array.ObjList(indObj(indSubj(i))).Array.addObj(x.Array.ObjList(indAcq)), ...
                    newArray(indNwArr(indSubj(i))));
            end
            %   Add new Subjects:
            if any(iNewSubj)
                obj.Array.addObj(newArray(iNewSubj));
            end
            
            if discardData
                %   Remove existing Acquisitions:
                indSubj = find(cellfun(@(x) any(x), iMissAcq));
                for i = 1:length(indSubj)
                    indAcq = find(iMissAcq{indSubj(i)});
                    obj.Array.ObjList(indObj(indSubj(i))).Array.removeObj(indAcq);
                end
                %   Remove existing Subjects:
                if any(iMissSubj)
                    obj.Array.removeObj(find(iMissSubj))
                end
            else
                warning('Keeping invalid Paths and/or files may cause problems later on during the analysis')
            end
            disp('update complete!')
        end
      function out = getFilePath(obj)
            % This function gets all the PATHS of the files from the acquisitions.
            out = [];
            acqList = arrayfun(@(x) x.Array.ObjList, obj.Array.ObjList, 'UniformOutput', false);
            for i = 1:length(acqList)
                for j = 1:length(acqList{i})
                    for k = 1:length(acqList{i}(j).Array.ObjList)
                        Folder = acqList{i}(j).Array.ObjList(k).Folder;
                        FileName = acqList{i}(j).Array.ObjList(k).FileName;
                        FullPath = fullfile(Folder, FileName)';
                        out = [out;FullPath];
                    end
                end
            end
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