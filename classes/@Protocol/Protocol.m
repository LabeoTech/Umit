classdef Protocol < handle
    % This class creates and manages a list of "Subject" objects.
    %   This class uses an user-defined function (ProtoFunc) to access the
    %   acquisition information from each subject and assign it to the
    %   pertinent objects. For more information, see documentation of the
    %   following classes: Subject, Acquisition and children of the abstract
    %   class "Modality".
    
    properties
        Name % Title of Project.
        MainDir % Folder containing all Raw recordings.
        SaveDir % Folder to save "Protocol" object and HDF5 files. Default value = current folder.
        ProtoFunc % Function handle of the user-defined OpenProtocol
        % where the Subjects and Acquisition data are created.
        Array % List of Subjects. Default: empty ObjectListManager.
        LogBookFile % MAT file with a table containing information about the Pipeline Operations run by PIPELINEMANAGER.
        garbageList % Table with a list of removed elements
        Idx_Filtered % Structure containing indices of Subject/Acquisition/Modality after using a Query Filter as OBJ.QUERYFILTER
        FilterStruct % Structure containing strings used to filter objects. Used by OBJ.QUERYFILTER.
    end
    properties (SetAccess = {?PipelineManager})
        LastLog % MAT file with a table containing information about the Last Pipeline Operations run by PIPELINEMANAGER.
    end
   
    methods
        
        function obj = Protocol(Name, MainDir, SaveDir, ProtoFunc, Array)
            % Class constructor.
            %   This function initiates the object "Protocol" with the
            %   properties: MainDir, SaveDir, ProtoFunc and Array.
            %   All first inputs must be provided. If Array is empty,
            %   the function creates an emtpy Array.
            if nargin > 0
                obj.Name = Name;
                obj.Array = Array;
                obj.MainDir = MainDir;
                obj.SaveDir = SaveDir;
                obj.ProtoFunc = ProtoFunc;
            else
                obj.Array = ObjectListManager();
            end
            obj.Idx_Filtered = {};
            obj.createLogBookFile
            obj.createFilterStruct
        end
        
        %%% Property Set Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set.Name(obj, Name)
            % Set function for Project Name
            %   Accepts only text.
            mustBeText(Name);
            obj.Name = Name;
        end
        function set.MainDir(obj, MainDir)
            % Set function for MainDir property.
            %   Accepts only existing Folders as input.
            MainDir = checkFolder(MainDir);
            mustBeFolder(MainDir); % Checks for existing Path.
            if isempty(obj.Array.ObjList)
                obj.MainDir = MainDir;
            else
                obj.changeMainDir(MainDir);
                obj.MainDir = MainDir;
            end
        end
        
        function set.SaveDir(obj, SaveDir)
            % Set function for SaveDir property.
            %   Accepts only existing Folders as input. Updates all
            %   SAVEFOLDERS of other Objects contained in OBJ.
            SaveDir = checkFolder(SaveDir);
            mustBeFolder(SaveDir); % Checks for existing Path.
            if isempty(obj.Array.ObjList)
                obj.SaveDir = SaveDir;
            else
                obj.changeSaveDir(SaveDir);
                obj.SaveDir = SaveDir;
            end
        end
        
        function set.ProtoFunc(obj, ProtoFunc)
            % Set function for ProtoFunc property.
            %   Accepts only valid function handles. Returns PROTOFUNC if
            %   empty.
            if ~isempty(ProtoFunc)
                mustBeA(ProtoFunc, 'function_handle');
                obj.ProtoFunc = ProtoFunc; % Checks if ProtoFunc is a function handle.
            end
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
            % Subjects and Acquisitions. Input is an Array of SUBJECT
            % objects.
            SubjArray = obj.ProtoFunc(obj);
            obj.Array.addObj(SubjArray);
            disp('List Generated')
        end
        
        function generateSaveFolders(obj)
            % GENERATESAVEFOLDERS creates directories in OBJ.SAVEDIR
            % containing folders for Subjects, Acquisitions and modalities.
            
            subjList = obj.Array.listProp('ID');
            acqList = arrayfun(@(x) x.Array.listProp('ID'), obj.Array.ObjList, 'UniformOutput', false);
            modList = arrayfun(@(x) arrayfun(@(y) y.Array.listProp('ID'), x.Array.ObjList, 'UniformOutput', false), obj.Array.ObjList, 'UniformOutput', false);
            
            for i = 1:length(subjList)
                tmpS = subjList{i};
                for j = 1:length(acqList{i})
                    tmpA = acqList{i}{j};
                    for k = 1:length(modList{i}{j})
                        tmpM = modList{i}{j}{k};
                        [status, msg, ~] = mkdir(fullfile(obj.SaveDir, tmpS, tmpA, tmpM));
                        if  status
                            obj.Array.ObjList(i).SaveFolder = fullfile(obj.SaveDir, tmpS);
                            obj.Array.ObjList(i).Array.ObjList(j).SaveFolder = fullfile(obj.SaveDir, tmpS, tmpA);
                            obj.Array.ObjList(i).Array.ObjList(j).Array.ObjList(k).SaveFolder = fullfile(obj.SaveDir, tmpS, tmpA, tmpM);
                        else
                            disp(['ERROR trying to create folder in path: ' tmpFolder])
                            disp(msg);
                        end
                    end
                end
            end
            disp('Save Folders generated!');
        end
        function updateList(obj, varargin)
            % This function updates the list of Subjects using
            % obj.ProtoFunc.
            %   The optional input (boolean) allows the user to not discard
            %   (FALSE) elements that were not found during the update.
            %       *This is not a wise option since keeping invalid
            %       paths of filenames may cause problems later on the
            %       analysis pipeline. "A voir..."
            %   IF discardData == FALSE, a warning message is thrown.
            if nargin < 2
                discardData = true;
            else
                discardData = varargin{1};
            end
            newArray = obj.ProtoFunc(obj);
            % Delete objects listed in OBJ.GARBAGELIST from NEWARRAY
            newArray = obj.removeGarbage(newArray);
            
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
            % Add new Folders to OBJ.SAVEDIR
            obj.generateSaveFolders();
            
            if discardData
                %   Remove existing Acquisitions:
                indSubj = find(cellfun(@(x) any(x), iMissAcq));
                for i = 1:length(indSubj)
                    indAcq = find(iMissAcq{indSubj(i)});
                    remInfo = obj.Array.ObjList(indObj(indSubj(i))).Array.removeObj(indAcq);
                    obj.garbageList = [obj.garbageList; remInfo];
                end
                %   Remove existing Subjects:
                if any(iMissSubj)
                    remInfo = obj.Array.removeObj(find(iMissSubj));
                    obj.garbageList = [obj.garbageList; remInfo];
                end
            else
                uiwait(warndlg('Keeping invalid Paths and/or files may cause problems later on during the analysis', 'Warning!', 'modal'));
            end
            uiwait(msgbox('Project update completed!'));
        end
        
        function manualRemoveObj(obj, args)
            % This function manually removes one Subject/Acquisition from
            % Protocol.
            arguments
                obj
                args.SubjectIndex (1,1) = 0;
                args.AcqIndex (1,1) = 0;
            end
            if args.SubjectIndex == 0
                return
            elseif args.AcqIndex == 0
                iRem = obj.Array.removeObj(args.SubjectIndex);
            else
                iRem = obj.Array.ObjList(args.SubjectIndex).Array.removeObj(args.AcqIndex);
            end
            obj.garbageList = [obj.garbageList; iRem];
        end
        
        function newObj = removeGarbage(obj, newObj)
            % This function removes from NEWOBJ elements (Subjects /
            % Acquisitions) listed in OBJ.GARBAGELIST.
            gbList = obj.garbageList;
            for i= 1:size(gbList,1)
                type = gbList(i,1);
                switch type
                    case 'Subject'
                        newObj(strcmp(gbList(i,2),{newObj.ID})) = [];
                    case 'Acquisition'
                        str = erase(gbList(i,3), obj.SaveDir);
                        str = split(str, filesep);
                        subjID = str{1};
                        idxS = find(strcmp(subjID, {newObj.ID}));
                        idx = find(strcmp(gbList(i,2), {newObj(idxS).Array.ObjList.ID}));
                        newObj(idxS).Array.removeObj(idx);
                end
            end
        end
        
        function purgeGarbageList(obj)
            % This function erases the OBJ.GARBAGELIST
            obj.garbageList = [];
            disp('Garbage List was reset!');
        end
        
        function clearLogBook(obj)
            % This function resets the LogBook saved in OBJ.LOGBOOKFILE
            obj.createLogBookFile;
        end
            
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function queryFilter(obj)
            % QUERYFILTER creates a cell array containing 1x3 array
            % of indices of Subject, Acquisition and Modality,
            % respectively. The function uses the property FILTERSTRUCT to
            % filter objects.
            %Reset Idx_Filtered
            obj.Idx_Filtered = [];
            % Filter Subjects
            indS = getIndex(obj, obj, obj.FilterStruct.Subject);
            % Filter Acquisitions
            for i = 1:length(indS)
                targetObj = obj.Array.ObjList(indS(i));
                indA = getIndex(obj, targetObj, obj.FilterStruct.Acquisition);
                % Filter Modality
                for j = 1:length(indA)
                    targetObj = obj.Array.ObjList(indS(i)).Array.ObjList(indA(j));
                    indM = getIndex(obj, targetObj, obj.FilterStruct.Modality);
                    for k = 1:length(indM)
                        obj.Idx_Filtered = [obj.Idx_Filtered ;[indS(i) indA(j) indM(k)]];
                    end
                end
            end
        end
        
        function clearFilterStruct(obj)
           % CLEARFILTERSTRUCT Resets the Filter Parameters to an empty
           % Structure.
            obj.createFilterStruct;
        end
        function LogBook = createEmptyTable(obj)
            % CREATEEMPTYTABLE outputs an empty Table to be filled with the information of pipelines
            % from PIPELINEMANAGER.
            LogBook = table({'None'}, {'None'}, {'None'}, {'None'}, {'None'},  {'None'}, {'None'}, ...
                0,datetime('now'),{'None'}, 'VariableNames', {'Subject', 'Acquisition',...
                'Recording', 'ClassName', 'Job', 'InputFile_UUID', 'InputFile_Path', 'Completed', 'RunDateTime', 'Messages'});
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function s = saveobj(obj)
            s.Name = obj.Name;
            s.MainDir = obj.MainDir;
            s.SaveDir = obj.SaveDir;
            s.ProtoFunc = obj.ProtoFunc;
            s.Array = obj.Array;
            s.LogBookFile = obj.LogBookFile;
            s.garbageList = obj.garbageList;
            s.LastLog = obj.LastLog;
            s.Idx_Filtered = obj.Idx_Filtered;
            s.FilterStruct = obj.FilterStruct;
        end
        
    end
    methods (Access = private)
        
        function createFilterStruct(obj)
            % CREATEFILTERSTRUCT creates an empty structure with the query info.
            Query = struct('PropName', '', 'Expression', '', 'LogicalOperator', '');
            obj.FilterStruct = struct('Subject', Query, 'Acquisition', Query, 'Modality', Query, 'FilterMethod', 'contains'); 
        end
        
        function createLogBookFile(obj)
            % This function creates an empty table LOGBOOK and saves in a .MAT file (LOGBOOKFILE).
            if ~isempty(obj.SaveDir)
                LogBook = obj.createEmptyTable;
                obj.LogBookFile = fullfile(obj.SaveDir, 'LogBook.mat');
                save(obj.LogBookFile, 'LogBook');
                disp(['Emtpy Log Book Created in ' obj.SaveDir]);
            end
        end
        
        
        function changeMainDir(obj, newMainDir)
            % CHANGEMAINDIR changes the main directory where the raw data
            % are located. Paths from raw data inside OBJ will be updated to the new
            % NEWMAINDIR.
            %   This function must be used after the user manually moves the
            %   files from OBJ.MAINDIR to NEWMAINDIR. The paths inside
            %   OBJ.MAINDIR must remain unchanged.
            
            newMainDir = checkFolder(newMainDir);
            if ~isempty(obj.Array.ObjList)
                for i = 1:length(obj.Array.ObjList)
                    tmpS = obj.Array.ObjList(i);
                    for j = 1:length(tmpS.Array.ObjList)
                        tmpA = tmpS.Array.ObjList(j);
                        for k = 1:length(tmpA.Array.ObjList)
                            tmpM = tmpA.Array.ObjList(k);
                            tmpM.RawFolder = strrep(tmpM.RawFolder, obj.MainDir, newMainDir);
                        end
                    end
                end
            end
        end
        
        function changeSaveDir(obj, newSaveDir)
            % CHANGESAVEDIR changes the main directory where the transformed data
            % are located. Paths from transformed data inside OBJ will be updated to the new
            % NEWSAVEDIR.
            %   This function must be used after the user manually moves the
            %   files from OBJ.SAVEDIR to NEWSAVEDIR. The paths inside
            %   OBJ.SAVEDIR must remain unchanged.
            
            newSaveDir = checkFolder(newSaveDir);
            if ~isempty(obj.Array.ObjList)
                for i = 1:length(obj.Array.ObjList)
                    tmpS = obj.Array.ObjList(i);
                    for j = 1:length(tmpS.Array.ObjList)
                        tmpA = tmpS.Array.ObjList(j);
                        for k = 1:length(tmpA.Array.ObjList)
                            tmpM = tmpA.Array.ObjList(k);
                            tmpM.SaveFolder = strrep(tmpM.SaveFolder, obj.SaveDir, newSaveDir);
                        end
                    end
                end
            end
        end
        function out = getIndex(obj, targetObj, FilterStruct)
            % GETINDEX applies a filter (OBJ.FILTERSTRUCT) to TARGETOBJ
            % and outputs a list of indices of filtered elements.
            
            ind = cell(1,length(FilterStruct));
            for i = 1:length(FilterStruct)
                Filt = FilterStruct(i);
                ind{i} = targetObj.Array.findElement(Filt.PropName, Filt.Expression, obj.FilterStruct.FilterMethod);
            end
            if length(FilterStruct) > 1
            for i = 1:length(FilterStruct)-1
                switch FilterStruct(i).LogicalOperator
                    case 'AND'
                        ind{i+1} = intersect(ind{i}, ind{i+1}); % Selects elements common to both arrays.
                    case 'OR'
                        ind{i+1} = setxor(ind{i}, ind{i+1}); % Selects elements that are NOT common to both arrays.
                    case 'NOT'
                        ind{i+1} = setdiff(ind{i}, ind{i+1}); % Selects elements from 1st array that are not in the 2nd array.
                    otherwise
                        ind{i+1} = union(ind{i}, ind{i+1}); % Selects the combined elements from both arrays.
                end         
            end
            end
            out = ind{end};
        end

        
        %         function delete(obj)
        %             disp('Protocol deleted')
        %         end
    end
    methods (Static)
        function obj = loadobj(s)
            if isstruct(s)
                newObj = Protocol;
                newObj.Name = s.Name;
                newObj.MainDir = s.MainDir;
                newObj.SaveDir = s.SaveDir;
                newObj.ProtoFunc = s.ProtoFunc;
                newObj.Array = s.Array;
                newObj.LogBookFile = s.LogBookFile;
                newObj.garbageList = s.garbageList;
                newObj.LastLog = s.LastLog;
                newObj.Idx_Filtered = s.Idx_Filtered;
                newObj.FilterStruct = s.FilterStruct;
                obj = newObj;
            else
                obj = s;
            end
        end
    end
end
