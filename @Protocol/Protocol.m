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
            % Subjects and Acquisitions. Input is an Array of SUBJECT
            % objects.
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
        
        function out = getFilePath(obj, varargin)
            % This function gets the PATHS of the files from the acquisitions.
            if nargin < 2
                FilterExp = createFilterStruct;
            else
                FilterExp = varargin{1};
            end
            out = [];
            indS = queryFilter(obj.Array, FilterExp.Subject);
            for i = 1:numel(indS)
                indA = queryFilter(obj.Array.ObjList(indS(i)).Array, FilterExp.Acquisition);
                for j = 1:numel(indA)
                    indM = queryFilter(obj.Array.ObjList(indS(i)).Array.ObjList(indA(j)).Array, FilterExp.Modality);
                    for k = 1:indM
                        Folder = obj.Array.ObjList(indS(i)).Array.ObjList(indA(j)).Array.ObjList(indM(k)).Folder;
                        FileName = obj.Array.ObjList(indS(i)).Array.ObjList(indA(j)).Array.ObjList(indM(k)).FileName;
                        FullPath = fullfile(Folder, FileName)';
                        out = [out;FullPath];
                    end
                end
            end
        end
        
        function out = getModalityProp(obj, PropName, varargin)
            % This function gets a list of the Properties from a Modality.
            %   PROPNAME can be a string or a cell of strings.
            if nargin < 3
                FilterExp = obj.createFilterStruct;
            else
                FilterExp = varargin{1};
            end
            if ~iscell(PropName)
                PropName = {PropName};
            end
            out = [];
            indS = queryFilter(obj.Array, FilterExp.Subject);
            for i = 1:numel(indS)
                indA = queryFilter(obj.Array.ObjList(indS(i)).Array, FilterExp.Acquisition);
                for j = 1:numel(indA)
                    indM = queryFilter(obj.Array.ObjList(indS(i)).Array.ObjList(indA(j)).Array, FilterExp.Modality);
                    for k = 1:indM
                        for p = 1:numel(PropName)
                            eval(['tmp = {obj.Array.ObjList(indS(i)).Array.ObjList(indA(j)).Array.ObjList(indM(k)).' PropName{p} '};']);
                            out = [out;tmp];
                        end
                    end
                end
            end
        end
        
        function runPreProcessingPipeline(obj, ppPipe, FilterExp)
            % This function runs a Pipeline using a ppPIPE structure.
            %    OPTIONS: query filter(FILTEREXP)
            outputFileName = ['ppPipe_Summary_' num2str(round(posixtime(datetime('now')))) '.csv'];
            filename = fullfile(pwd, outputFileName);
            fid = fopen(filename, 'w'); 
            fprintf(fid,'Subject, Acquisition, FunctionName, Status, RunDateTime\n');
            if isempty(FilterExp)
                FilterExp = obj.createFilterStruct;
            end
            indS = queryFilter(obj.Array, FilterExp.Subject);
            for i = 1:numel(indS)
                tmpS = obj.Array.ObjList(indS(i));
                indA = queryFilter(tmpS.Array, FilterExp.Acquisition);
                for j = 1:numel(indA)
                    tmpA = tmpS.Array.ObjList(indA(j));
                    indM = queryFilter(tmpA.Array, FilterExp.Modality);
                    load(tmpA.LogBookFile);
                    for k = 1:indM
                        pipeSteps = fieldnames(ppPipe);
                        tmpM = tmpA.Array.ObjList(indM(k));
                        for p = 1:numel(pipeSteps)
                            ppStruct = ppPipe.(pipeSteps{p});
                            ppStruct.FuncParams = strrep(ppStruct.FuncParams, 'FOLDER',['''' tmpM.Folder '''']);
                            isLogged = checkInLogBook(LogBook, ppStruct);
                            if ~isLogged
                                %%% RUN THE PIPELINE STEP %%%
                                eval(['tmpM.run_' ppStruct.FuncName '(' ppStruct.FuncParams ')'])
                                %%% LOG NEW ENTRY IN LOGBOOK %%%
                                tmpM.LastLog.FunctionName = {ppStruct.FuncName};
                                tmpM.LastLog.FuncParams = {ppStruct.FuncParams};
                                tmpA.save2LogBook(tmpM.LastLog);
                                failed = ~tmpM.LastLog.Completed;
                                if failed
                                    %%% IF JOB FAILED, SKIP THE REST OF THE
                                    %%% PIPELINE.
                                    fprintf(fid, '%s, %s, %s, %s, %s\n', ...
                                        tmpS.ID, tmpA.ID, ppStruct.FuncName, 'FAILED' ,  tmpM.LastLog.RunDateTime);
                                    disp(repmat('*', 1, 100))
                                    disp(['The function ''' ppStruct.FuncName ''' failed.'])
                                    disp('All dependent functions in the pipeline were aborted.')
                                    disp(['Check the LogBook in ' tmpA.LogBookFile ' and re-run the Pipeline']);
                                    disp(repmat('*', 1, 100))
                                    break
                                end
                                fprintf(fid, '%s, %s, %s, %s, %s\n', ...
                                        tmpS.ID, tmpA.ID, ppStruct.FuncName, 'COMPLETED' ,  tmpM.LastLog.RunDateTime);
                            else
                                disp(repmat('*', 1, 100))
                                disp(['The function ''' ppStruct.FuncName ''' was previously run on the Data and will be skipped.'])
                                disp('Trying to run the next step on the Pipeline...')
                                disp(repmat('*', 1, 100))
                            end
                        end
                    end
                end
            end
            fclose(fid);
        end
 
        function FilterExp = createFilterStruct(~)
            % Creates an empty structure with the query info used in some
            % "get" functions.
            Query = struct('PropName', [], 'Expression', [], 'logicalOperator', []);
            FilterExp = struct('Subject', Query, 'Acquisition', Query, 'Modality', Query);
        end
        function delete
            disp('Protocol deleted')
        end
    end
end

% Local functions

function isLogged = checkInLogBook(LogBook, pipeStruct)
% This function checks if the job in the pipeline (PIPESTRUCT) has already
% been processed.
isLogged = false;
if isempty(LogBook)
    return
end
a = strcmp(pipeStruct.ClassName, LogBook.ModalityName);
b = strcmp(pipeStruct.FuncName, LogBook.FunctionName);
FP_pipe = regexprep(pipeStruct.FuncParams, '\s', '');
FP_logBook = regexprep(LogBook.FuncParams, '\s', '');
c = strcmp(FP_pipe, FP_logBook);
d = LogBook.Completed == 1;

if sum(a & b & c & d)
    isLogged = true;
end
end

function ind = queryFilter(objectArray, FilterExp)
% This function applies a filter (FILTEREXP) to the OBJECTARRAY and outputs the
% indices of filtered elements.

if isempty(FilterExp.PropName)|| isempty(FilterExp)
    ind = 1:length(objectArray.ObjList);
    return
end
idx = false(length(FilterExp), length(objectArray.ObjList));
for i = 1:length(FilterExp)
    PropName = strtrim(FilterExp.PropName);
    exp = FilterExp.Expression;
    idx(i,:) = objectArray.findElement(PropName, exp);
end

if ~isempty(FilterExp(1).logicalOperator)
    tmp = idx(1,:);
    for i = 1:length(FilterExp) - 1
        logicOp = strtrim(FilterExp(i).logicalOperator);
        if strcmpi(logicOp, 'AND') || strcmp(logicOp, '&')
            tmp = ( tmp & idx(i+1,:) );
        elseif strcmpi(logicOp, 'OR') || strcmp(logicOp, '|')
            tmp = ( tmp | idx(i+1,:) );
        else
            warning(['Invalid logical operator : ' FilterExp(i).logicalOperator '. Used an AND(&) instead. Check query parameters and try again.']);
            tmp = ( tmp & idx(i+1,:) );
        end
    end
    ind = find(tmp);
else
    ind = find(idx);
end
end
