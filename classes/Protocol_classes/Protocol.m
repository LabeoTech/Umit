classdef Protocol < handle
    % This class creates and manages a list of "Subject" objects.
    %   This class uses an user-defined function (ProtoFunc) to access the
    %   acquisition information from each subject and assign it to the
    %   pertinent objects. For more information, see documentation of the
    %   following classes: Subject, Acquisition and children of the abstract
    %   class "Modality".
    
    properties
        Name char % Title of Project.
        MainDir  % List of folders containing all Raw recordings.
        SaveDir char % Folder to save "Protocol" object and HDF5 files. Default value = current folder.
        ProtoFunc % Function handle of the user-defined OpenProtocol
        % where the Subjects and Acquisition data are created.
        Array % List of Subjects. Default: empty ObjectListManager.
        garbageList % Table with a list of removed elements
        Idx_Filtered  = []% Array containing indices of Subject/Acquisition/Modality after using a Query Filter as OBJ.QUERYFILTER
        FilterStruct % Structure containing strings used to filter objects. Used by OBJ.QUERYFILTER.
    end
    properties (SetAccess = {?PipelineManager})
        LastLog % MAT file with a table containing information about the Last Pipeline Operations run by PIPELINEMANAGER.
    end
    properties (SetAccess = private)
        b_isDummy  = false % Used by DataViewer as standalone ONLY!!
    end
    properties (Dependent)
        LogBookFile char % MAT file with a table containing information about the Pipeline Operations run by PIPELINEMANAGER.
    end
    methods
        function obj = Protocol(Name, MainDir, SaveDir, ProtoFunc, Array, varargin)
            % Class constructor.
            %   This function initiates the object "Protocol" with the
            %   properties: MainDir, SaveDir, ProtoFunc and Array.
            %   All first inputs must be provided. If Array is empty,
            %   the function creates an emtpy Array.
            
            if nargin > 0                
                if ~isempty([varargin{:}])
                    b_isDummy = varargin{:};
                end
            else
                % Create dummy protocol object, if no input is provided
                % (used for testing and dev.);
                mod = FluorescenceImaging('DummyMod',{''},'Dummy',1);
                acq = Acquisition('DummyAcq',[]);acq.Array.addObj(mod);
                sbj = Subject('DummySbj','def',''); sbj.Array.addObj(acq);
                Array = sbj;
                Name = 'DummyProtocol';
                MainDir = pwd;
                SaveDir = pwd;
                ProtoFunc = @(x) x;
                b_isDummy = true;                
            end
            obj.Name = Name;
            obj.Array = Array;
            obj.MainDir = MainDir;
            obj.SaveDir = SaveDir;
            obj.ProtoFunc = ProtoFunc;
            obj.b_isDummy = b_isDummy;
            
            %
            if ~obj.b_isDummy
                obj.createLogBookFile
            end
            obj.createFilterStruct
        end
        %%% Property Set Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set.Name(obj, Name)
            % Set function for Project Name
            %   Accepts only text.
            validateattributes(Name,{'char'}, {'nonempty'}, 'set.Name');
            obj.Name = Name;
        end
        function set.MainDir(obj, MainDir)
            % Set function for MainDir property.
            %   Accepts only existing Folders as input.       
            validateattributes(MainDir,{'char', 'cell'}, {'nonempty'}, 'set.MainDir');
            if ischar(MainDir)
                MainDir = {MainDir};
            end            
            MainDir = checkFolder(MainDir, 'raw');
            
            if isempty(obj.MainDir)
                obj.MainDir = MainDir;
            else                
                obj.changeMainDir(MainDir);
                obj.MainDir = unique(MainDir,'stable');
            end
            
        end
        function set.SaveDir(obj, SaveDir)
            % Set function for SaveDir property.
            %   Accepts only existing Folders as input. Updates all
            %   SAVEFOLDERS of other Objects contained in OBJ.
            SaveDir = checkFolder(SaveDir, 'new');
            obj.validate_path(SaveDir); % Checks for existing Path.
            obj.SaveDir = SaveDir;
        end
        function set.ProtoFunc(obj, ProtoFunc)
            % Set function for ProtoFunc property.
            %   Accepts only valid function handles. Returns PROTOFUNC if
            %   empty.
            if ~isempty(ProtoFunc)
                % Check if ProtoFunc is a function handle:
                validateattributes(ProtoFunc, {'function_handle'}, {'nonempty'}, 'set.ProtoFunc');
                obj.ProtoFunc = ProtoFunc;
            end
        end
        function set.Array(obj, Array)
            % Set function for Array property.
            %   Accepts "ObjectListManager" or "Subject" objects as input. If
            %   empty, creates an default "ObjectListManager" object.
            obj.Array = ObjectListManager([],obj);
            if isa(Array, 'Subject')
                obj.Array.addObj(Array);
            elseif isa(Array, 'ObjectListManager')
                obj.Array = Array;
                %             else
                %                 obj.Array = ObjectListManager([],obj);
            end
        end
        %%% Property Get functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get.LogBookFile(obj)
            % Get function for getting the full path of Log Book .MAT file.
            % It is dependent from the property "SaveFolder":
            out = fullfile(obj.SaveDir, 'LogBook.mat');
        end
        %%% Validators %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function validate_path(~, input)
            msgID = 'umIToolbox:Protocol:InvalidInput';
            str = strrep(input, filesep, repmat(filesep,1,2));
            msg = ['The input "' str '" is not a valid folder or it is not in MATLAB''s path:'];
            assert(isfolder(input),msgID, msg)
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function generateList(obj)
            % This function uses the ProtoFunc to create the lists of
            % Subjects and Acquisitions. Input is an Array of SUBJECT
            % objects.
            
            sList = [];
            for ii = 1:length(obj.MainDir)
                sList = [sList,obj.ProtoFunc(obj.MainDir{ii})];%#ok
            end
            obj.Array.addObj(obj.mergeSubjectList(sList));
            disp('List Generated')
        end
        function generateSaveFolders(obj)
            % GENERATESAVEFOLDERS creates directories in OBJ.SAVEDIR
            % containing folders for Subjects, Acquisitions and modalities.
            subjList = obj.Array.listProp('ID');
            acqList = arrayfun(@(x) x.Array.listProp('ID'), obj.Array.ObjList, 'UniformOutput', false);
            modList = arrayfun(@(x) arrayfun(@(y) y.Array.listProp('ID'),...
                x.Array.ObjList, 'UniformOutput', false), obj.Array.ObjList, 'UniformOutput', false);
            
            for i = 1:length(subjList)
                tmpS = subjList{i};
                indx = obj.Array.findElement('ID', tmpS, 'strcmp');
                tmpS_obj = obj.Array.ObjList(indx);
                [~,~] = mkdir(fullfile(obj.SaveDir, tmpS_obj.ID));
                
                for j = 1:length(acqList{i})
                    tmpA = acqList{i}{j};
                    indx = tmpS_obj.Array.findElement('ID', tmpA, 'strcmp');
                    tmpA_obj = tmpS_obj.Array.ObjList(indx);
                    [~,~] = mkdir(fullfile(obj.SaveDir, tmpS_obj.ID, tmpA_obj.ID));
                    for k = 1:length(modList{i}{j})
                        tmpM = modList{i}{j}{k};
                        [~,~] = mkdir(fullfile(obj.SaveDir, tmpS, tmpA, tmpM));
                    end
                end
            end
            disp('Save Folders generated!');
        end
        function [objArray, idxStruc,msg] = lookForNewData(obj, b_verbose)
            % This function re-runs the protocol function and compares the
            % newly created array with the current object array in
            % "protocol".
            % Input:
            %   b_verbose(bool): If TRUE, outputs a formatted text with a
            %   summary with the number of new/missing items in "protocol".
            % Outputs:
            %   objArray(Subject): Subject array after running the protocol
            %       fcn.
            %   idxStruc (struct): list of boolean indices for new/missing
            %       items in protocol.
            %   msg (char): If "b_verbose" is TRUE, contains a text with
            %       total number of new/missing items in "protocol".
            
%             iNewSubj = {};
%             iMissSubj = {};
            iNewAcq = {};
            iMissAcq = {};
            
            idxStruc = struct('iNewSubj',{},'iMissSubj', {}, 'iNewAcq', {},'iMissAcq',{});
            objArray = [];
            msg = '';
            for ii = 1:length(obj.MainDir)
                addpath(obj.MainDir{ii}); % To ensure that the proto function is found.
                objArray = [objArray, obj.ProtoFunc(obj.MainDir{ii})];%#ok
            end
            objArray = obj.mergeSubjectList(objArray);
            if isempty(objArray)                
                warning('The protocol function returned an empty output! Is the raw folder empty?')
                return
            end
            
            % Delete objects listed in OBJ.GARBAGELIST from NEWARRAY
            objArray = obj.removeGarbage(objArray);
            
            % 1st, control for new or deleted Subjects:
            iNewSubj = ~ismember({objArray.ID}, {obj.Array.ObjList.ID}); % New subjects in the newArray.
            iMissSubj = ~ismember({obj.Array.ObjList.ID},{objArray.ID});% Subjects from original list that are no longer in the newArray.
            
            % 2nd, control for new or deleted Acquisitions from each Subject:
            indNwArr = find(~iNewSubj);
            if indNwArr
                % Locate the sujects from "newArray" in "obj":
                [~,indObj] = ismember({objArray(indNwArr).ID},{obj.Array.ObjList.ID});
                % Get new Acquisitions from existing Subjects:
                iNewAcq = arrayfun(@(a,b) ~ismember({a.Array.ObjList.ID}, {b.Array.ObjList.ID}),...
                    objArray(indNwArr), obj.Array.ObjList(indObj), 'UniformOutput', false);
                % Get missing Acquisitions from existing Subjects:
                iMissAcq = arrayfun(@(a,b) ~ismember({b.Array.ObjList.ID}, {a.Array.ObjList.ID}),...
                    objArray(indNwArr), obj.Array.ObjList(indObj), 'UniformOutput', false);
                % Exclude "virtual" acquisitions from the missing list:
                indSubj = find(cellfun(@any, iMissAcq));
                for ii = 1:length(indSubj)
                    % REMOVE virtual ACQS:
                    iMissAcq{indSubj(ii)} = iMissAcq{indSubj(ii)} & ~[obj.Array.ObjList(indObj(indSubj(ii))).Array.ObjList.b_isVirtual];
                end
            end
            % Create message with total number of new/missing items:
            if b_verbose
                msg = sprintf(['List of New/missing Items:\nNewSubjects: %d\n'...
                    'Missing Subjects: %d\nNewAcquisitions: %d\nMissingAcquisitions: %d'],...
                    sum(iNewSubj), sum(iMissSubj), sum([iNewAcq{:}]),sum([iMissAcq{:}]));
            else
                msg = '';
            end
            % Populate output structure:
            idxStruc(1).iNewSubj = iNewSubj;
            idxStruc(1).iMissSubj = iMissSubj;
            idxStruc(1).iNewAcq = iNewAcq;
            idxStruc(1).iMissAcq = iMissAcq;
        end
        function updateList(obj, varargin)
            % This function updates the list of Subjects using
            % obj.ProtoFunc.
            % Input:
            %   b_discardData(bool) : (Optional) If TRUE, discards elements
            %   that were not found in the raw folder (Default).
            %   If FALSE, elements that were not found during the update
            %   will be kept.
            %       *This is not a wise option since keeping invalid
            %       paths of filenames may cause problems later on the
            %       analysis pipeline. *Consider to eliminate this option in
            %       a later Release of umIT if not pertinent.
            
            if nargin < 2
                discardData = true;
            else
                discardData = varargin{1};
            end
            % Look for new data:
            [newArray,idxInfo] = obj.lookForNewData(false);
            if isempty(newArray)
                return
            end
            %%% Update the existing list:
            % Add new Acquisitions:
            indSubj = find(cellfun(@any, idxInfo.iNewAcq));
            for ii = 1:length(indSubj)
                indAcq = find(idxInfo.iNewAcq{indSubj(ii)});
                % Find subject in "protocol":
                currSubj = obj.Array.ObjList(indSubj(ii));
                % Find subject in new array:
                indSubjNew = find(strcmp({newArray.ID},currSubj.ID));
                % Add new acquisition(s) to the existing Subject:
                arrayfun(@(x) currSubj.Array.addObj(x.Array.ObjList(indAcq)), ...
                    newArray(indSubjNew));%#ok
            end
            %   Add new Subjects:
            if any(idxInfo.iNewSubj)
                obj.Array.addObj(newArray(idxInfo.iNewSubj));
            end
            % Add new Folders to OBJ.SAVEDIR
            obj.generateSaveFolders();
            if discardData
                %   Remove existing Acquisitions:
                indSubj = find(cellfun(@any, idxInfo.iMissAcq));
                for ii = 1:length(indSubj)
                    indAcq = find(idxInfo.iMissAcq{indSubj(ii)});
                    remInfo = obj.Array.ObjList(indSubj(ii)).Array.removeObj(indAcq);%#ok
                    obj.garbageList = [obj.garbageList; remInfo];
                end
                %   Remove existing Subjects:
                if any(idxInfo.iMissSubj)
                    remInfo = obj.Array.removeObj(find(idxInfo.iMissSubj));%#ok
                    obj.garbageList = [obj.garbageList; remInfo];
                end
                
                
                % Update element's properties, if changes were made in the
                % protocol function:
                for ii = 1:length(newArray)
                    % Compare Subjects' properties:
                    idxS = strcmp(newArray(ii).ID, obj.Array.listProp('ID'));
                    propInfo = metaclass(newArray(ii));
                    settableProps = protocolCanSet(propInfo.PropertyList);
                    % Compare settable properties between the old and new
                    % arrays:
                    b_HasChanged = cellfun(@(x) ~isequaln(newArray(ii).(x), obj.Array.ObjList(idxS).(x)), settableProps);
                    if any(b_HasChanged)
                        cellfun(@(x) evalin('caller',['obj.Array.ObjList(' num2str(find(idxS)) ').' x ' = newArray(' num2str(ii) ').' x ';']), ...
                            settableProps(b_HasChanged))
                    end
                    % Compare Acquisitions:
                    for jj = 1:length(newArray(ii).Array.ObjList)
                        idxA = strcmp(newArray(ii).Array.ObjList(jj).ID, obj.Array.ObjList(idxS).Array.listProp('ID'));
                        propInfo = metaclass(newArray(ii).Array.ObjList(jj));
                        settableProps = protocolCanSet(propInfo.PropertyList);
                        % Compare settable properties between the old and new
                        % arrays:
                        b_HasChanged = cellfun(@(x) ~isequaln(newArray(ii).Array.ObjList(jj).(x), obj.Array.ObjList(idxS).Array.ObjList(idxA).(x)), settableProps);
                        if any(b_HasChanged)
                            cellfun(@(x) evalin('caller',['obj.Array.ObjList(' num2str(find(idxS)) ').Array.ObjList(' num2str(find(idxA)) ').' x...
                                ' = newArray(' num2str(ii) ').Array.ObjList(' num2str(jj) ').' x ';']), ...
                                settableProps(b_HasChanged))
                        end
                        % Compare modalities:
                        for kk = 1:length(newArray(ii).Array.ObjList(jj).Array.ObjList)
                            idxM = strcmp(newArray(ii).Array.ObjList(jj).Array.ObjList(kk).ID, obj.Array.ObjList(idxS).Array.ObjList(idxA).Array.listProp('ID'));
                            propInfo = metaclass(newArray(ii).Array.ObjList(jj).Array.ObjList(kk));
                            settableProps = protocolCanSet(propInfo.PropertyList);
                            % Compare settable properties between the old and new
                            % arrays:
                            b_HasChanged = cellfun(@(x) ~isequaln(newArray(ii).Array.ObjList(jj).Array.ObjList(kk).(x), obj.Array.ObjList(idxS).Array.ObjList(idxA).Array.ObjList(idxM).(x)), settableProps);
                            if any(b_HasChanged)
                                cellfun(@(x) evalin('caller',['obj.Array.ObjList(' num2str(find(idxS)) ').Array.ObjList(' num2str(find(idxA)) ').Array.ObjList(' num2str(find(idxM)) ').' x...
                                    ' = newArray(' num2str(ii) ').Array.ObjList(' num2str(jj) ').Array.ObjList(' num2str(kk) ').' x ';']), ...
                                    settableProps(b_HasChanged))
                            end
                        end
                    end
                end
            end
            disp('Project update completed!');
            % Local functions
            function propNames = protocolCanSet(propList)
                b_pass = false(size(propList));
                propNames = {propList.Name}';
                for ind = 1:length(propList)
                    % Remove MyParent From Settable Array:
                    if ismember(lower(propList(ind).Name), {'myparent', 'array', 'lastlog'})
                        continue
                    end
                    
                    if ischar(propList(ind).SetAccess)
                        b_pass(ind) = strcmpi(propList(ind).SetAccess, 'public');
                    else
                        b_pass(ind) = any(cellfun(@(x) strcmpi(x.Name, 'protocol'), propList(ind).SetAccess));
                    end
                end
                propNames = propNames(b_pass);
            end
            
        end
        function updateMainDir(obj)
            % This method updates the RawFolder and RawFiles properties of
            % all modalities when the user moved the raw data between
            % MainDir folders.
            
            obj.changeMainDir(obj.MainDir);
        end
        function out = getSelectedItems(obj)
            % This method generates a list of full IDs (Subject --
            % Acquisition -- Modality) from the selected data stored in
            % "Idx_filtered" property.
            
            items = obj.extractFilteredObjects(3);
            out = cellfun(@(x) strjoin({x.MyParent.MyParent.ID, x.MyParent.ID, x.ID},' -- '), items, 'UniformOutput',false);
        end
        function [modHandle, AcqHandle] = manualAddModality(obj, modClass, modID, AcqID, subjHandle)
            % MANUALADDMODALITY creates a new Acquisition inside a "Subject" Object.
            % and adds the modality object of class "modClass" as it's child. The
            % Acquisition ID is stored in "acqID".
            % The Subject object handle is "subjHandle".
            % The output are the handles for the new modality and
            % acquisition.
            %         **Use the outputs to edit the elements' properties.**
            disp('Manually adding object to Parent...');
            modHandle = [];
            AcqHandle = [];
            %%% Validation:
            % Check if the modality exists:
            errID = 'umIToolbox:Protocol:WrongInput';
            if ~exist(modClass, 'class')
                error(errID, ['The modality ' modClass 'does not exist!']);
            end
            if ~isa(subjHandle, 'Subject')
                error(errID,'The subject handle is invalid!');
            end
            % Check if the acquisition was added to the subject:
            idxS = obj.Array.findElement('ID',subjHandle.ID);
            if any(strcmp(AcqID, obj.Array.ObjList(idxS).Array.listProp('ID')))
                warning('Acquisition already exists in the selected Subject!')
                return
            end
            % Create Acquisition:
            AcqHandle = Acquisition();
            AcqHandle.ID = AcqID;
            % Create modality:
            eval(['modHandle = ' modClass '();'])
            modHandle.ID = modID;
            modHandle.RawFiles_FP = {'MISSING'}; % Create dummy raw file name in new modality.
            % Add the modality to the Acquisition:
            AcqHandle.Array.addObj(modHandle)
            % Add the acquisition to the Subject:
            subjHandle.Array.addObj(AcqHandle);
            % Create Save folder for the new acquisition:
            obj.generateSaveFolders;
            % Mark Acquisition as "virtual":
            AcqHandle.b_isVirtual = true;
            disp('Done')
        end
        function manualRemoveObj(obj, SubjectIndex, varargin)
            % This function manually removes one Subject/Acquisition from
            % Protocol.
            % Inputs:
            %   SubjectIndex (int): Index of the Subject object in Protocol
            %   AcqIndex (int): Index of the Acquisition object in Protocol
            %   b_delFolder(bool): If TRUE, delete SaveFolder of the object
            %   to be deleted and its contents.
            p = inputParser;
            addOptional(p,'SubjectIndex', 0, @isnumeric);
            addOptional(p,'AcqIndex', 0, @isnumeric);
            addOptional(p, 'b_delFolder', true, @islogical);
            parse(p, SubjectIndex, varargin{:})
            if p.Results.SubjectIndex == 0
                return
            elseif p.Results.AcqIndex == 0
                iRem = obj.Array.removeObj(p.Results.SubjectIndex);
            else
                iRem = obj.Array.ObjList(p.Results.SubjectIndex).Array.removeObj(p.Results.AcqIndex);
            end
            if isempty(iRem)
                return
            end
            if p.Results.b_delFolder && iRem(3) ~= ""
                answer = questdlg(['Delete Folder "' iRem(3) '" and it''s contents?'],...
                    'Delete Save Folder?', 'Yes', 'No', 'No');
                if strcmp(answer, 'Yes')
                    rmdir(iRem(3), 's') % Delete save folder and its contents.
                end
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
                        idxS = strcmp(gbList(i,2),{newObj.ID});
                        if isempty(idxS)
                            continue
                        end
                        newObj(idxS) = [];
                    case 'Acquisition'
                        str = split(gbList(i,3), filesep); % Get subject ID from parsing the folder path.
                        if str == ""
                            % skip if the folder doesn't exist. This should
                            % not be reached.
                            continue
                        end
                        subjID = str{end-1};
                        idxS = strcmp(subjID, {newObj.ID});
                        if isempty(idxS)
                            continue
                        end
                        idxAcq = strcmp(gbList(i,2), {newObj(idxS).Array.ObjList.ID});
                        newObj(idxS).Array.ObjList(idxAcq) = [];
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
        function out = extractFilteredObjects(obj, lvl)
            %EXTRACTFILTEREDOBJECTS generates a cell array containing all
            % objects listed in OBJ.IDX_FILTERED at the level LVL.
            % Input:
            % lvl:
            % 1 - Subject
            % 2 - Acquisition
            % 3 - Modality
            out = {};
            switch lvl
                case 1
                    indx = unique(obj.Idx_Filtered(:,1));
                    out = arrayfun(@(x) obj.Array.ObjList(x), indx, 'UniformOutput', false);
                case 2
                    indx = unique(obj.Idx_Filtered(:,[1 2]), 'rows');
                    out = arrayfun(@(x,y) obj.Array.ObjList(x).Array.ObjList(y),...
                        indx(:,1), indx(:,2), 'UniformOutput', false);
                case 3
                    indx = unique(obj.Idx_Filtered, 'rows');
                    out = arrayfun(@(x,y,z) obj.Array.ObjList(x).Array.ObjList(y).Array.ObjList(z),...
                        indx(:,1), indx(:,2), indx(:,3), 'UniformOutput', false);
            end
        end
        function addTextEvent(obj, Text, dateAndTime, varargin)
            % This method adds an event to a "Text_events.mat" file.
            % This method is a wrapper of the function "create_Text_eventsFile.m".
            
            % Inputs:
            %   Text(char): String that will be added to "Text_events.mat"
            %       file.
            %   dateAndTime(datetime): timestamp associated with the Text.
            %   flag (char): (Optional) {'add', 'overwrite'}. Use 'add'
            %   (default) to add a Text to the file OR 'overwrite' to
            %   ovewrite the content of the "Text_events.mat" file with a
            %   new Text.
            
            default_flag = 'add';
            p = inputParser;
            validFlag = @(x) ischar(x) && ismember(x, {'add', 'overwrite'});
            validText = @(x) ischar(x) || (iscell(x) && (ischar([x{:}])));
            addRequired(p, 'Text', validText);
            addRequired(p, 'dateAndTime', @isdatetime);
            addOptional(p,'flag', default_flag, validFlag);
            parse(p, Text, dateAndTime, varargin{:});
            eventID = p.Results.Text;
            timestamps = p.Results.dateAndTime;
            flag = p.Results.flag;
            
            txtEv_file = fullfile(obj.SaveDir, 'Text_events.mat');
            if ~isfile(txtEv_file) && strcmp(flag, 'add')
                flag = 'overwrite';
                answer = questdlg(['Text Event file not found in ' obj.SaveDir '. Create a new file?'], ...
                    'Create Text Event File', 'Yes', 'Cancel', 'Cancel');
                if strcmp(answer, 'Cancel') || isempty(answer)
                    disp('Operation cancelled by user')
                    return
                end
            elseif isfile(txtEv_file) && strcmp(flag, 'overwrite')
                answer = questdlg('Text Event file will be overwritten. Proceed?', ...
                    'Create Text Event File', 'Yes', 'No. Add a new entry instead', 'Cancel', 'Cancel');
                switch answer
                    case {'Cancel', 0}
                        disp('Operation cancelled by user')
                        return
                    case 'No. Add a new entry'
                        flag = 'add';
                end
            end
            
            switch flag
                case 'add'
                    f = '-a';
                otherwise
                    f = '-w';
            end
            create_Text_eventsFile(obj.SaveDir, eventID, timestamps, f)
        end
        function queryFilter(obj)
            % QUERYFILTER creates a cell array containing 1x3 array
            % of indices of Subject, Acquisition and Modality,
            % respectively. The function uses the property FILTERSTRUCT to
            % filter objects.
            
            new_Idx_Filtered = [];
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
                        new_Idx_Filtered = [new_Idx_Filtered ;[indS(i) indA(j) indM(k)]];%#ok
                    end
                end
            end
            
            if isempty(new_Idx_Filtered)
                return
            end
            % Remove objects that are not in the pre-filtered array:
            if ~isempty(obj.Idx_Filtered)
                new_Idx_Filtered = new_Idx_Filtered(ismember(new_Idx_Filtered,obj.Idx_Filtered,'rows'),:);
            end
            obj.Idx_Filtered = new_Idx_Filtered;
        end
        function clearFilterStruct(obj)
            % CLEARFILTERSTRUCT erases the list of filtered objects and
            % resets the Filter Parameters to an empty structure.
            obj.Idx_Filtered = [];
            obj.createFilterStruct;
        end
        function LogBook = createEmptyTable(~)
            % CREATEEMPTYTABLE outputs an empty Table to be filled with the information of pipelines
            % from PIPELINEMANAGER.
            LogBook = table({'None'}, {'None'}, {'None'}, {'None'}, {'None'},...
                {'None'}, {'None'}, 0,datetime('now'),{'None'}, {'None'},...
                'VariableNames', {'Subject', 'Acquisition', 'Recording', 'ClassName','TaskName',...
                'Job', 'InputFile_Path', 'Completed', 'RunDateTime', 'Messages', 'Messages_short'});
        end
        function varargout = save(obj, saveFolder)
            % Saves the protocol object in the SaveFolder as the project's
            % name. It handles retrocompatibility.
            if ~exist('saveFolder','var')
                saveFolder = obj.SaveDir;
            end
            
            filename = fullfile(saveFolder, obj.Name);
            if isfile([filename '.prt']) || (~isfile([filename '.prt']) && ~isfile([filename '.mat']))
                filename = [filename '.prt'];
            else
                % For retrocompatibility
                filename = [filename '.mat'];
            end
            protocol = obj;
            save(filename, 'protocol','-mat');     
            disp(['Protocol saved as: "' filename '"']);
            if nargout
                varargout = {filename};
            end
        end
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
            s.b_isDummy = obj.b_isDummy;
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
                save(obj.LogBookFile, 'LogBook');
                disp(['Emtpy Log Book Created in ' obj.SaveDir]);
            end
        end
        function changeMainDir(obj, newMainDir)
            % CHANGEMAINDIR changes the main directory where the raw data
            % are located. Paths from raw data inside OBJ will be updated to the new
            % NEWMAINDIR.
            %   This function must be used only AFTER the user manually moves the
            %   files from OBJ.MAINDIR to NEWMAINDIR. The paths inside
            %   OBJ.MAINDIR must remain unchanged.            
            newMainDir = checkFolder(newMainDir);
            curr_idxFiltered = obj.Idx_Filtered;
            obj.clearFilterStruct;obj.queryFilter;
            modItems = obj.extractFilteredObjects(3);
            for ii = 1:length(modItems)
                tmpM = modItems{ii};
                thisFolder = obj.MainDir(cellfun(@(x) startsWith(tmpM.RawFolder,x), obj.MainDir));
                if isempty(thisFolder)
                    continue
                end
                thisFolder = thisFolder{:};                
                for jj = 1:length(newMainDir)
                    newDir = newMainDir{jj};
                    oldPath = tmpM.RawFiles_FP{1};
                    newPath = strrep(oldPath, thisFolder,newMainDir{jj});
                    if isfile(newPath) && ~strcmp(oldPath, newPath)
                        tmpM.RawFiles_FP = cellfun(@(x) strrep(x, thisFolder,newDir),tmpM.RawFiles_FP, 'UniformOutput',false);
                        continue
                    end
                end
            end
            obj.Idx_Filtered = curr_idxFiltered;      
        end
        function sList = mergeSubjectList(~, sList)
            % This function is a helper function for all public methods
            % that add new subjects to the protocol using the protocol
            % function. In brief, it merges the acquisitions of duplicate
            % subjects in the list "sList". The presence of duplicates may
            % arise from the existence of different acquisitions of the
            % same subject across different raw folders (i.e. MainDir(s)).
            sIDlist = unique({sList.ID},'stable');
            idxDelete = false(size(sList));
            for ii = 1:length(sIDlist)
                % Merge acquisitions of a single subject:
                sIdx = find(strcmp(sIDlist{ii},{sList.ID}));
                for jj = 2:length(sIdx)
                    sList(sIdx(1)).Array.addObj(sList(sIdx(jj)).Array.ObjList);
                end
                idxDelete(sIdx(2:end)) = true;
            end
            sList(idxDelete) = []; % Remove duplicates
        end
        function out = getIndex(obj, targetObj, FilterStruct)
            % GETINDEX applies a filter (OBJ.FILTERSTRUCT) to TARGETOBJ
            % and outputs a list of indices of filtered elements.
            % Input:
            %   targetObj : Protocol Object of class Subject, Acquisition
            %       or Modality.
            %   FilterStruct(struct): Protocol's filter structure
            %   containing info for filtering objects inside "targetObj".
            
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
    end
    methods (Static)
        function obj = loadobj(s)
            if isstruct(s)
                % Update Save Directory based on the current location
                % of the protocol file:
                if isfile(which([s.Name, '.prt']))
                    ext = '.prt';
                    s.SaveDir = fileparts(which([s.Name, ext]));
                elseif isfile(which([s.Name, '.mat']))
                    % For retrocompatibility
                    ext = '.mat';
                    s.SaveDir = fileparts(which([s.Name, ext]));
                end
                newObj = Protocol;
                newObj.Name = s.Name;
                % Check MainDir and SaveDir existance:
                errID = 'umIToolbox:Protocol:InvalidInput';
                errMsg = ' is not an existing folder!';
                assert(isfolder(s.SaveDir),errID, [strrep(s.MainDir,filesep, repmat(filesep,1,2)), errMsg]);
                %%%
                newObj.MainDir = s.MainDir;
                newObj.SaveDir = s.SaveDir;
                newObj.ProtoFunc = s.ProtoFunc;
                newObj.Array = s.Array;
                newObj.garbageList = s.garbageList;
                newObj.LastLog = s.LastLog;
                newObj.Idx_Filtered = s.Idx_Filtered;
                newObj.FilterStruct = s.FilterStruct;
                if isfield(s,'b_isDummy')
                    % For retrocompatibility
                    newObj.b_isDummy = s.b_isDummy;
                end
                obj = newObj;
            else
                obj = s;
            end
            % Add Main and Save directories to Matlab's path (this may be SLOW for directories with a lot of folders...):
            cellfun(@(x) addpath(genpath(x)), obj.MainDir);
            addpath(genpath(obj.SaveDir));
            % Rebuild handle of "MyParent" property of elements from
            % Protocol:
            obj.Array.parentObj = obj; % Update handle in ObjectListManager.
            % Update handles in Subjects, Acquisitions and Modalities:
            for i = 1:numel(obj.Array.ObjList)
                obj.Array.ObjList(i).MyParent = obj;
                for j = 1:numel(obj.Array.ObjList(i).Array.ObjList)
                    obj.Array.ObjList(i).Array.ObjList(j).MyParent = obj.Array.ObjList(i);
                    for k = 1:numel(obj.Array.ObjList(i).Array.ObjList(j).Array.ObjList)
                        obj.Array.ObjList(i).Array.ObjList(j).Array.ObjList(k).MyParent = obj.Array.ObjList(i).Array.ObjList(j);
                    end
                end
            end
        end
    end
end
