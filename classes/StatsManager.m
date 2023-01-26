classdef StatsManager < handle
    % STATSMANAGER manages data formatted to statistical analysis (see
    % documentation to create functions for statistics). This class
    % provides several methods to output summary statistics, perform
    % hypothesis tests, normality tests and other options.
    % Inputs:
    %  list_of_objs : cell array of elements from "Protocol" object.
    %  stats_filename: filename of .MAT file saved using "save2Mat.m"
    %  function.
    
    properties
        list_of_objs   % Cell array of elements from "Protocol" object.
        obs_list       % current list of observations from stats_filename.
        list_of_groups % Cell array of group names of objects from list_of_objs.
        list_of_events % current list of event Names from all objects.
        % !!If the input data has no events, the list will contain the string"NoEvents" !!
        list_of_subjects % Current list of Subject Names.
        stats_filename % Filename of .MAT file saved using "save2Mat.m"
        % function.
        MfileArr       % Array of MATFILE objects.
        time_resolution = 'none' % Time resolution for grouping each observation
        % Options: "none", "minute", "hour, "day", "week", "month".
        %%%%% STATS PARAMETERS %%%%%
        pAlpha = 0.05 % "Significance" threshold (Alpha) for statistical analysis.
    end
    properties (SetAccess = private)
        %         timestamp_list % List of timestamps associated with each object in list_of_objs.
        b_hasStatsToolbox  % True, if Matlab contains the Statistics and Machine learning toolbox.
        inputFeatures % Structure containing some information about the input data. This will be used by plotting tools and the umIToolbox app.
        dataArr  % structure with data extracted from "stats_data" that can be averaged using the method "setAcquisitionRange".
        flatData % structure derived from "dataArr" used by some plotting tools and stats functions.
        obs_list_original % Cell Array of observations from stats_filename when this object was created.
        list_of_events_original% Cell array of event Names from all objects when this object was created.
        list_of_subjects_original % Cell array of Subject Names from all objects when this object was created.
        curr_test = ''; % Name of the current test from the list_of_tests.
        indepVars % Independent variables used to perform statistical comparisons.
        splitVar % Data dimension to use to split the data into subsets for statistical analysis.
        results_stats % structure containing the results of hypothesis tests between groups.
    end
    properties (Access = private)
        stats_data  = {} % cell array containing all data and metaData created.
        avg_stats_data = {} % cell array similar to stats_data containing the average values of acquisitions. (See method averageData).
        headers = {} % cell array with the _stats_data column names as keys and indices as values.
        
        % Structure with information about all possible independent variables (see prop "indepVar") for statistics:
        indepVarInfo = struct('Name',{'Group','Subject','ROI','Acquisition', 'Event'}, ...
            'fieldName',{'groupID','SubjectID','ObsID','AcquisitionID','EventID'},...
            'indexName',{'gIndx','sIndx','rIndx','aIndx','eIndx'},...
            'b_isRepeatedMeasure',{false,false,false,true,true});       
    end
    
    methods
        function obj = StatsManager(list_of_objs, obs_list, list_of_groups, stats_filename)
            % Class constructor.
            %   This function initiates the object "StatsManager" with the
            %   properties: list_of_objs, obs_list, list_of_groups and stats_filename.
            %   All inputs must be provided.
            obj.list_of_objs = list_of_objs;
            obj.obs_list = obs_list; obj.obs_list_original = obs_list;
            obj.list_of_groups = list_of_groups;
            obj.stats_filename = stats_filename;
            % Check if inputs are correct and Generate list of MatFile
            % handles:
            obj.validateObject;
            obj.createDataArray;
            obj.getEventList; % Get list of events and store it to "list_of_events" property.
            obj.list_of_events_original = obj.list_of_events;
            obj.gen_dataArr;
            obj.list_of_subjects = unique({obj.dataArr.SubjectID});
            obj.list_of_subjects_original = obj.list_of_subjects;
            obj.validateData; % Identifies the type of input data and other properties.
            
            % Check if Matlab's Stats. toolbox exist:
            a = ver;
            if any(strcmp('Statistics and Machine Learning Toolbox', {a.Name}))
                obj.b_hasStatsToolbox = true;
                obj.setStatsVariables; % Update list of independent variable(s) used in statistical comparisons.
            else
                obj.b_hasStatsToolbox = false;
            end
        end
        %%% Property Set Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set.list_of_objs(obj,list_of_objs)
            % Set function for list_of_objs.
            % It checks if the input is a cell array of valid objects.
            types = cellfun(@(x) isa(x, 'Protocol') || isa(x, 'Subject') ||...
                isa(x, 'Acquisition') || isa(x, 'Modality'), list_of_objs);
            errID = 'Umitoolbox:StatsManager:InvalidInput';
            if sum(types) == 0
                errMsg = ['Invalid list of objects. Objects must be a cell array ' ...
                    'containing Protocol/Subject/Acquisition/Modality objects'];
                error(errID, errMsg);
            elseif ~iscell(list_of_objs) || isempty(list_of_objs)
                errMsg = ['Wrong input type. List of Objects must be a non-empty cell array.' ...
                    'containing Protocol/Subject/Acquisition/Modality objects'];
                error(errID, errMsg);
            else
                obj.list_of_objs = list_of_objs;
            end
        end
        
        function set.obs_list(obj,obs_list)
            % Set function for obs_list.
            % It checks if obs_list is a cell array of CHAR.
            if iscell(obs_list) && ischar([obs_list{:}])
                obj.obs_list = obs_list;
            else
                errID = 'umIToolbox:StatsManager:WrongInput';
                errMsg = 'List of observations must be a non-empty cell array of characters';
                error(errID, errMsg);
            end
        end
        
        function set.list_of_groups(obj,list_of_groups)
            % Set function for list_of_groups.
            % It checks if list_of_groups is a cell array of CHAR.
            if iscell(list_of_groups) && ischar([list_of_groups{:}])
                obj.list_of_groups = list_of_groups;
            else
                errID = 'umIToolbox:StatsManager:WrongInput';
                errMsg = 'List of groups must be a non-empty cell array of characters';
                error(errID, errMsg);
            end
        end
        
        function set.stats_filename(obj,stats_filename)
            % Set function for stats_filename.
            % It checks if stats_filename is a .MAT file.
            errID = 'umIToolbox:StatsManager:InvalidInput';
            errMsg = 'Wrong input. Stats file must be a .MAT file';
            assert(isa(stats_filename, 'char') & endsWith(stats_filename, '.mat'),...
                errID, errMsg);
            obj.stats_filename = stats_filename;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [out,uniqLabels] = createTable(obj,varargin)
            % This function creates a table or an array of tables
            % containing ROI data from each observation.
            % Input:
            % tableType (char) =  ['raw' (default), 'summary']
            % Output:
            %   out (table): table containing data and
            %       metadata from all observations in obj.stats_data structure.
            %   uniqLabels(cell array of char): list of all labels found in
            %       the "stats_data" structure.
            
            p = inputParser;
            addRequired(p, 'obj');
            addOptional(p,'tableType', 'raw', @(x) ismember(x, {'raw', 'summary'}));
            parse(p,obj, varargin{:});
            obj = p.Results.obj;
            tableType = p.Results.tableType;
            %
            disp('Creating table...')
            % Get all unique labels from "dataArr" structure. This will
            % be used by "getObsData" method to put NaNs on missing
            % data (i.e. labels that are missing for a given recording).
            if size(obj.dataArr(1).labels,1) < size(obj.dataArr(1).labels,2)
                labelList = [obj.dataArr.labels]';
            else
                labelList = vertcat(obj.dataArr.labels); labelList(cellfun(@isempty, labelList)) = []; % remove empty labels.
            end
            uniqLabels = unique(labelList, 'stable'); % Sort unique labels in alphabetical order.
            % Get observations's labels from the stats_data structure:
            if strcmp(tableType, 'raw')
                tableArr = cellfun(@(x) obj.getObsData(x, uniqLabels), obj.obs_list,'UniformOutput',false);
                out = vertcat(tableArr{:});
            else
                disp('No!')
            end
            disp('Table created!')
        end
        
        function exportToCSV(obj,filename)
            % This function creates a .CSV file containing all data created
            % by the method "createTable".
            % Input:
            % filename (char): valid path for a .CSV file.
            if ~obj.inputFeatures.b_isExportable
                error(['Data type ' obj.inputFeatures.dataType ' cannot be exported to CSV'])
            end
            [data, labels] = obj.createTable;
            % Transform data table to cell array to avoid invalid variable names from labels:
            varNames = data.Properties.VariableNames;
            varNames(7:end) = labels;
            data = table2cell(data);
            data = vertcat(varNames, data);
            % Write to .csv file:
            disp('Writing table to .CSV file...')
            writecell(data,filename);
            msgbox(['Data saved to file : ' filename], 'to CSV');
        end
        
        function out = getAcqIndexList(obj,type)
            % GETACQINDEXLIST provides the list of the acquisition indices
            % that are available for merging by "setAcquisitionRange"
            % method OR the original list of acquisition indices.
            % Input:
            %   type(char): 'original', 'available', 'current'
            % Output:
            %   out (int or cell): list of acquisition indices as a numeric
            %   array for input type 'original' OR 'available' or as a
            %   cell array for type 'current'. The order of the cell array
            %   corresponds to the index of the "indx_avg_data" in
            %   stats_data.
            
            if ~ismember(lower(type), {'original', 'available', 'current'})
                error('Input should be either "original", "available" or "current".')
            end
            
            switch lower(type)
                case 'original'
                    out =  unique([obj.stats_data{:, obj.hMap('AcquisitionIndx')}]);
                case 'available'
                    idx = [obj.stats_data{:,obj.hMap('indx_avg_data')}]' == 0;
                    out =  unique([obj.stats_data{idx, obj.hMap('AcquisitionIndx')}]);
                case 'current'
                    indxAvg = [obj.stats_data{:,obj.hMap('indx_avg_data')}];
                    acqIndx = [obj.stats_data{:,obj.hMap('AcquisitionIndx')}];
                    b_isCurrent = indxAvg ~=0;
                    if ~any(b_isCurrent)
                        out = {};
                        return
                    end
                    for ii = 1:sum(unique(indxAvg) ~= 0)
                        out{ii} = sort(acqIndx(indxAvg == ii));%#ok
                    end
            end
        end
        
        function setAcquisitionRange(obj,indxRange)
            % SETACQUISITIONRANGE regroups a range of acquisitions to be
            % merged using the method "averageData".
            % Input:
            %   indxRange(1x2 int): low and high indices of the acquisition
            %   to be merged.
            
            % Block function usage ONLY for equal number of acquisitions:
            assert(obj.inputFeatures.b_hasSameAcqN & obj.inputFeatures.b_AcqHasSameDimSize,...
                'Grouping acquisitions is available for data with equal number of acquisitions and with the same data size per Subject!')
            %%%%% Input validation: %%%%%
            indxRange = round(indxRange);
            % Get available acquisition ranges:
            av_range = [min([obj.stats_data{[obj.stats_data{:, obj.hMap('indx_avg_data')}] == 0,obj.hMap('AcquisitionIndx')}]), ...
                max([obj.stats_data{[obj.stats_data{:, obj.hMap('indx_avg_data')}] == 0,obj.hMap('AcquisitionIndx')}])];
            assert(~isempty(av_range), 'All acquisitions already grouped! Reset the average data array and try again!')
            assert(isequal(indxRange, sort(indxRange)), 'Invalid index range!Please provide the low and high acquisition index in ascending order.');
            assert(all(isnumeric(indxRange)), 'Invalid input! Index range must be numeric')
            assert(numel(indxRange) == 2, 'Invalid input! Number of indices must be equal to 2!')
            assert(min(indxRange >= min(av_range)) & max(indxRange) <= max(av_range),...
                ['Index range must be between ' num2str(av_range(1)) ' and ' num2str(av_range(2)) '!'])
            %%%%%%%%%%
            % Update indices in the "stats_data" array:
            nextIndx = max([obj.stats_data{:,obj.hMap('indx_avg_data')}])+1;
            idx = ismember([obj.stats_data{:,obj.hMap('AcquisitionIndx')}], [indxRange(1):indxRange(2)]);
            obj.stats_data(idx,obj.hMap('indx_avg_data')) = repmat({nextIndx},sum(idx),1);
            disp(['Acquisitions ' num2str(indxRange(1)) ' to ' num2str(indxRange(2))...
                ' added to average data array as group No. #' num2str(nextIndx) '.']);
            obj.gen_dataArr;
        end
        
        function setStatsVariables(obj,varargin)
            % This method sets the "indepVars" property based on the
            % available data in "dataArr".
            % Inputs:
            %   varNames (cell): name(s) of the 2 data dimension to set as
            %       Independent variable(s).
            %   groupVar (char): name of the dimension to split the data.
            %   For instance, if the groupVar is "Observation", the statistical
            %   analysis will be performed for each one of them separately.
            
            p = inputParser;
            addRequired(p, 'obj');
            name_validation = @(x) all(ismember(x, {obj.indepVarInfo.Name}));
            addOptional(p, 'varNames', {'Acquisition','Group'},@(x) name_validation(x) & numel(x) == 2);
            addOptional(p, 'groupVar', 'ROI',@(x) name_validation(x));
            parse(p, obj, varargin{:});
            varNames = p.Results.varNames;
            groupVar = p.Results.groupVar;
            % Further validation:
            if any(strcmpi(groupVar, varNames))
                error('Grouping variable cannot be an independent variable');
            end
            % Use this only for "scalar" data (FOR NOW):
            if ~strcmpi(obj.inputFeatures.dataType, 'scalar')
                warning(['Statistical comparisons not available for data of type ' obj.inputFeatures.dataType]);
                return
            end
            disp(repmat('-',1,100));
            fprintf('Current independent variable(s) are "%s" and "%s".\n', varNames{:});
            [~,lb] = ismember(varNames, {obj.indepVarInfo.Name});
            obj.indepVars = obj.indepVarInfo(lb);
            fprintf('Statistical analysis will be applied to each "%s" item.\n', groupVar);
            disp(repmat('-',1,100));
            obj.splitVar = obj.indepVarInfo(strcmpi(groupVar, {obj.indepVarInfo.Name}));
            disp('Done');
        end
        
        function out = getStatsVariables(obj)
            % This function gives the list of stats variables in the
            % following order:
            % 1 - independent variable 1
            % 2 - independent variable 2
            % 3 - split variable
            out = {obj.indepVars.Name, obj.splitVar.Name};
        end
        
        function resetAvgIndex(obj)
            % RESETAVGINDEX sets all acquisition indices in the average
            % data to zero.
            obj.stats_data(:,obj.hMap('indx_avg_data')) = repmat({0},size(obj.stats_data,1),1);
            disp('Acquisition average group indices reset!')
            obj.gen_dataArr;
        end
        
        function gen_dataArr(obj)
            % GEN_DATAARR updates the values of dataArr property. This
            % function is used by the public methods 'setAcquisitionRange'
            % and 'resetAvgIndex'.
            
            % "dataArr" is a structure containing the data and meta data
            % stored in this class. The data will average any acquisition with
            % the value in the column "indx_avg_data" higher than zero.
            
            
            disp('updating data Array...')
            obj.dataArr = struct.empty(0,1);
            idx = [obj.stats_data{:,obj.hMap('indx_avg_data')}]' == 0;
            % Populate structure with all acquisitions that will not be
            % averaged:
            indxZero = find(idx);
            myHeaders = setdiff(obj.headers, {'indx_avg_data'}, 'stable'); % Remove non-pertinent columns.
            for i = 1:length(indxZero)
                for j = 1:length(myHeaders)
                    obj.dataArr(i).(myHeaders{j}) = obj.stats_data{indxZero(i),obj.hMap(myHeaders{j})};
                end
            end
            % Average acquisitions and append to non-averaged data:
            gNames = unique(obj.list_of_groups);
            tmp = {};
            for iG = 1:numel(gNames)
                idxG = strcmp(obj.stats_data(:, obj.hMap('groupID')), gNames(iG));
                sNames = unique(obj.stats_data(idxG, obj.hMap('SubjectID')));
                for iS = 1:numel(sNames)
                    idxS = strcmp(obj.stats_data(:,obj.hMap('SubjectID')), sNames(iS));
                    indxAcq = [obj.stats_data{idxG & idxS, obj.hMap('indx_avg_data')}];
                    nMerge = setdiff(unique(indxAcq), 0);
                    for iAcq = 1:length(nMerge)
                        idxA = [obj.stats_data{:, obj.hMap('indx_avg_data')}]' == nMerge(iAcq);
                        tmp = [tmp; {averageData(obj, obj.stats_data(idxG & idxS & idxA,:))}];%#ok
                    end
                end
            end
            if isempty(obj.dataArr)
                obj.dataArr = horzcat(tmp{:});
            else
                obj.dataArr = horzcat(obj.dataArr,tmp{:});
            end
            % Remap acquisition indices of non-averaged data to have a
            % continuous range of acquisition indices.
            for iG = 1:numel(gNames)
                idxG = strcmp({obj.dataArr.groupID}, gNames{iG});
                sNames = unique({obj.dataArr(idxG).SubjectID});
                for iS = 1:numel(sNames)
                    idxS = strcmp({obj.dataArr.SubjectID}, sNames{iS});
                    acqList = sort([obj.dataArr(idxG & idxS).AcquisitionIndx]);
                    newList = [min(acqList):numel(acqList)];
                    mapIndx = containers.Map(acqList, newList);
                    indxAcq = find(idxG & idxS);
                    for iA = 1:length(indxAcq)
                        obj.dataArr(indxAcq(iA)).AcquisitionIndx = mapIndx(obj.dataArr(indxAcq(iA)).AcquisitionIndx);
                    end
                end
            end
            % Add some indices to obj.dataArr:
            idxErase = false(size(obj.dataArr));
            for ind = 1:length(obj.dataArr)
                obj.dataArr(ind).gIndx = find(strcmp(obj.dataArr(ind).groupID, unique(obj.list_of_groups)));
                obj.dataArr(ind).sIndx = find(strcmp(obj.dataArr(ind).SubjectID, unique({obj.dataArr.SubjectID})));
                [idxObs, obj.dataArr(ind).rIndx] = ismember(obj.dataArr(ind).observationID, obj.obs_list);
                if obj.dataArr(ind).b_hasEvents
                    % Flag if events exist:
                    [idxE,obj.dataArr(ind).eIndx] = ismember(obj.dataArr(ind).MatFile.eventNameList, obj.list_of_events);
                    evID = obj.dataArr(ind).MatFile.eventID;
                    evNames = obj.dataArr(ind).MatFile.eventNameList;
                    dimE = find(strcmpi(obj.dataArr(ind).MatFile.dim_names, 'E'));
                    if ( numel(evID) > numel(evNames) )
                        % Automatically average any event repetitions in
                        % data:
                        disp('Averaging events ...');
                        for jnd = 1:numel(obj.dataArr(ind).data)
                            obj.dataArr(ind).data{jnd} = averageEvents(obj.dataArr(ind).data{jnd});
                        end
                    end
                else
                    idxE = true;
                    obj.dataArr(ind).eIndx = 0;
                end
                % Manage missing events:
                if ( obj.dataArr(ind).b_hasEvents )
                    if  any(idxE)
                        indxE = find(strcmp(obj.dataArr(ind).MatFile.dim_names,'E'));
                        permE = [indxE, setdiff(1:numel(obj.dataArr(ind).MatFile.dim_names), indxE)];
                        
                        % Remove data corresponding to missing events:
                        for kk = 1:length(obj.dataArr(ind).data)
                            tmp = permute(obj.dataArr(ind).data{kk},permE);
                            sz = size(tmp);
                            tmp = reshape(tmp,sz(1),[]);
                            tmp(~idxE,:) = []; % Erase missing events
                            sz(1) = sum(idxE);
                            tmp = reshape(tmp,sz);
                            obj.dataArr(ind).data{kk} = ipermute(tmp,permE);
                        end
                        % Update meta data:
                        obj.dataArr(ind).eIndx(~idxE) = [];
                        obj.dataArr(ind).labels(~idxE) = [];
                        clear tmp sz indxE permE
                    else
                        idxErase(ind) = true; % Erase everything.
                    end
                end
                % Remove data corresponding to missing observations:
                if any(idxObs)
                    obj.dataArr(ind).data(~idxObs) = [];
                    obj.dataArr(ind).dataSize(~idxObs) = [];
                    obj.dataArr(ind).rIndx(~idxObs) = [];
                    obj.dataArr(ind).observationID(~idxObs) = [];
                else
                    idxErase(ind) = true;
                end
            end
            % Erase empty items from datArr:
            obj.dataArr(idxErase) = [];
            
            disp('Done');
            % Local functions:
            function out = averageData(obj, dataIn)
                % AVERAGEDATA calculates the average of all data (from "stats_data")
                % set with indices greater than zero in the column
                % "indx_avg_data".
                % Output:
                %   out (struct): structure containing the average of  the tagged acquisitions.
                % Instantiate output variable:
                out = struct();
                % Average data
                cols2copy = {'groupID','SubjectID','ModalityID', 'labels',...
                    'dataSize', 'MatFile', 'b_hasEvents'};
                for ii = 1:length(cols2copy)
                    out.(cols2copy{ii}) = dataIn{1,obj.hMap(cols2copy{ii})};
                end
                % Set as baseline if one of the elements in dataIn is a
                % baseline:
                out.b_isBaseline = any([dataIn{:,obj.hMap('b_isBaseline')}]);
                % Use the earliest recording start datetime:
                [~,k] = min(datetime(string(dataIn(:,obj.hMap('RecStartDateTime')))));
                out.('RecStartDateTime') = dataIn{k,obj.hMap('RecStartDateTime')};
                % Average the data per observation:
                obsList = unique(vertcat(dataIn{:,obj.hMap('observationID')}), 'stable');
                %
                dimCat = numel(dataIn{1,obj.hMap('MatFile')}.dim_names) + 1;
                avg = {};
                currObs = {};
                for iOb = 1:length(obsList)
                    % Loop across each observation and average the data:
                    idxOb = cellfun(@(x) ismember(x, obsList(iOb)),...
                        dataIn(:,obj.hMap('observationID')), 'UniformOutput',false);
                    dat = cellfun(@(x,y) x(y),dataIn(:,obj.hMap('data')), idxOb, 'UniformOutput',false);
                    idxEmpty = cellfun(@isempty,dat);
                    if all(idxEmpty)
                        continue
                    else
                        dat = [dat{~idxEmpty}]';
                        % Calculate average:
                        avg = [avg; {mean(cat(dimCat,dat{:}),dimCat,'omitnan')}];%#ok
                        currObs = [currObs; obsList(iOb)];%#ok
                    end
                end
                out.data = avg;
                out.observationID = currObs;
                out.AcquisitionID = ['AverageAcq_' num2str(dataIn{1,obj.hMap('indx_avg_data')})];
                out.AcquisitionIndx = min([dataIn{:,obj.hMap('AcquisitionIndx')}]);
                out.dataFile = '';
            end
            function out = averageEvents(dataIn)
                % AVERAGEEVENTS averages the event repetitions in "dataIn".
                newOrder = [dimE, setdiff(1:ndims(dataIn),dimE)];
                dataIn = permute(dataIn, newOrder);
                datsz = size(dataIn);
                dataIn = reshape(dataIn, datsz(1), []);
                out = zeros([numel(evNames), size(dataIn,2)], 'single');
                for jj = 1:numel(evNames)
                    out(jj,:) = mean(dataIn(evID == jj,:), 1,'omitnan');
                end
                out = reshape(out, [numel(evNames), datsz(2:end)]);
                out = ipermute(out, newOrder);
            end
        end
        
        function flatDataArr(obj)
            % This function flattens the data Array into one vector for
            % scalar data and a X by T matrix for time-series data. The
            % dimensions indices are also provided to retrieve the meta
            % data information.
            
            if ( any(strcmpi(obj.inputFeatures.dataType,{'map','matrix'})) )
                error(['Operation aborted. Plot option not available for data type ' ...
                    obj.inputFeatures.dataType])
            end
            % Repackage data for plotting:
            disp('repackaging data...')
            % Check for time-series:
            %             [hasT,indxT] = ismember('T', obj.inputFeatures.dim_names);
            indxT = find(strcmpi('t',obj.inputFeatures.dim_names));
            dataSize = arrayfun(@(x) cellfun(@(y) size(y), x.data, 'UniformOutput',false), obj.dataArr, 'UniformOutput', false)';
            dataSize = vertcat(dataSize{:}); dataSize = vertcat(dataSize{:});
            if ( indxT )
                Xsz = max(dataSize(indxT,:));
            else
                Xsz = 1; % For scalar data.
            end
            % Preallocate output array with maximum possible data size:
            nRec = numel(obj.dataArr);
            nEv = numel(obj.list_of_events);
            nROI = numel(obj.obs_list);
            datOut = nan(prod([nRec,nEv,nROI]), Xsz, 'single');
            % Populate output data structure:
            obj.flatData = struct('data',[],'gIndx', [],'sIndx', [],'rIndx',[],'eIndx',[],'aIndx',[],...
                'groupID',{},'SubjectID',{},'ObsID',{},'EventID',{},'AcquisitionID', {},...
                'AcquisitionNames',{},'RecStartDateTime',{});
            % *Here, "AcquisitionID" are just generic strings (e.g. Acq#1,
            % Acq#2...) while "AcquisitionNames" are the original IDs of
            % the recordings.
            indxEv = find(strcmpi('E', obj.inputFeatures.dim_names));
            cnt = 1;
            for ii = 1:length(obj.dataArr)
                % Concatenate the data:
                dat = vertcat(obj.dataArr(ii).data{:});
                if ( indxEv )
                    datSz1 = size(obj.dataArr(ii).data{1},indxEv) * numel(obj.dataArr(ii).data);
                    dat = reshape(permute(dat,[indxEv, setdiff(1:ndims(dat),indxEv)]),datSz1,[]); % Reshape data;
                else
                    datSz1 =  length(obj.dataArr(ii).data);
                end
                
                if ( indxT )
                    datSz2 = size(obj.dataArr(ii).data{1},indxT);
                else
                    datSz2 = 1;
                end
                % Add data to output array:
                datOut(cnt:cnt+datSz1-1,1:datSz2) = dat;
                % Get meta data:
                RecStartDateTime = repmat({obj.dataArr(ii).RecStartDateTime},datSz1,1);
                AcqNames = repmat({obj.dataArr(ii).AcquisitionID},datSz1,1);
                aIndx = repmat(obj.dataArr(ii).AcquisitionIndx, datSz1,1);
                gIndx = repmat(obj.dataArr(ii).gIndx,datSz1,1);
                sIndx = repmat(obj.dataArr(ii).sIndx,datSz1,1);
                rIndx = repelem(obj.dataArr(ii).rIndx, datSz1/length(obj.dataArr(ii).rIndx),1);
                eIndx = repmat(obj.dataArr(ii).eIndx, datSz1/length(obj.dataArr(ii).eIndx),1);
                % Add meta data to output structure:
                obj.flatData(1).gIndx = [obj.flatData.gIndx; gIndx];
                obj.flatData.sIndx = [obj.flatData.sIndx; sIndx];
                obj.flatData.rIndx = [obj.flatData.rIndx; rIndx];
                obj.flatData.eIndx = [obj.flatData.eIndx; eIndx];
                obj.flatData.aIndx = [obj.flatData.aIndx; aIndx];
                obj.flatData.AcquisitionNames = [obj.flatData.AcquisitionNames; AcqNames];
                obj.flatData.RecStartDateTime = [obj.flatData.RecStartDateTime; RecStartDateTime];
                cnt = cnt + datSz1;
            end
            obj.flatData.data = datOut(~all(isnan(datOut),2),:);
            obj.flatData.groupID = unique(obj.list_of_groups);
            obj.flatData.SubjectID = obj.list_of_subjects;
            obj.flatData.ObsID = obj.obs_list;
            obj.flatData.EventID = obj.list_of_events;
            obj.flatData.AcquisitionID = arrayfun(@(x) ['Acq#' num2str(x)],...
                unique(obj.flatData.aIndx), 'UniformOutput',false)';
            disp('Done.')
        end
        
        function resetLists(obj)
            % RESETLISTS will reset the observation and event lists to the
            % original state when this object was created. Additionally, the
            % data array "dataArr" will be updated.
            
            obj.obs_list = obj.obs_list_original;
            obj.list_of_events = obj.list_of_events_original;
            % Update data array:
            obj.gen_dataArr;
        end
        %%% Stats functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function varargout = runStats(obj)
            % This function checks if the data can be used for statistical
            % analysis (using the available tests) and performs the
            % statistical analysis on the data.
            % IMPORTANT: !! currently, this method is available only for
            % "scalar" data!!
            % Output:
            %   statsInfo (optional): structure with the results of
            %   hypothesis tests.
            
            obj.results_stats = [];            
            % Validation:
            if ~obj.b_hasStatsToolbox
                varargout{1} = warning('Operation Aborted! "Statistics and Machine learning" toolbox is necessary!');
                return
            end
            if ~strcmpi(obj.inputFeatures.dataType, 'scalar')
                varargout{1} = warning('Operation Aborted! Statistical Comparison available ONLY for "scalar" data types!');
                return
            end
            % Use the independent variables list ("indepVars") and the data
            % from "dataArr" to decide which statistical test to use.
            obj.chooseStatsTest;
            if strcmpi(obj.curr_test, 'unknown')
                varargout{1} = warning('Operation Aborted! Failed to Identify Hypothesis test based on the data.');
                return
            end
            % Validate assumptions to use parametric tests;
            if ~obj.checkNormality && any(contains(obj.curr_test, {'Two','Repeated'}, 'IgnoreCase',true))
                % If the data is not normal for ANOVA test, abort!
                varargout{1} = warning('The data does not follow a normal distribution. Statistical tests will not be applied!');
                return
            elseif ~obj.checkNormality
                % If is not normal for Paired and Unpaired 2-sampled data,
                % switch comparison to non-parametric versions:
                ws = warning('The data does not follow a normal distribution. Non-parametric tests will be performed.');
            end
            if obj.checkHomogeneity && any(contains(obj.curr_test, {'Two','Repeated'}, 'IgnoreCase',true))
                ws = warning('The data does not contain homogeneous variances. Interpret results with caution!');               
            elseif ~obj.checkHomogeneity
                ws = warning('The data does not contain homogeneous variances. Non-parametric tests will be performed. Interpret results with caution!');         
            end
            obj.runHypothesisTest;
            if nargout
                varargout{1} = obj.results_stats;
                varargout{2} = ws;
            end                        
        end 
        
        function genStatsReport(obj,filename)
           % This method creates a .TXT file with the information from the
           % "results_stats" property.
           
           if ~endsWith(filename,'.txt','IgnoreCase',true)
               filename = [filename '.txt'];
           end
           
           if isempty(obj.results_stats)
               warning('Statistical results not found! Run statistical comparisons in order to generate a report');
               return
           end
           
           div = repmat('-',1,80);
           % Add header with some general information about the data and
           % the statistical comparisons performed:           
           str = sprintf('-----------------------Statistical Hypothesis test report-----------------------\nRun date: %s\n\n',datestr(datetime('now')));
                      
           str = [str, sprintf(['The data was grouped by %d %s(s) vs %d %s(s) and split by a total of %d %s(s).\n\n'...
               'The hypothesis test executed was "%s"\n%s\n'],...
               obj.indepVars(1).len, obj.indepVars(1).Name,obj.indepVars(2).len, obj.indepVars(2).Name, ...
               length(obj.results_stats), obj.splitVar.Name, obj.curr_test, div)];                      
           % Write the results of each test.            
           hypMap = containers.Map([true,false],{'Yes','No'});
           for ii = 1:length(obj.results_stats)
               % Add name of Split Variable:
               str = [str, sprintf('Subgroup name: %s\n',obj.results_stats(ii).subGroupName)];
               % Add outputs of statistical tests:
               info = obj.results_stats(ii).GroupComparison;
               switch obj.curr_test                    
                    case {'PairedTest','UnpairedTest'}
                        str = [str, sprintf('Comparison performed: "%s"\nRejected null hypothesis: "%s"\np-value: %0.4f\nStats:\n',...
                            info.Name, hypMap(info.p < obj.pAlpha),info.p)];  
                            % Add "stats" fields:
                            fn = fieldnames(info.stats);                            
                            for kk = 1:length(fn)
                                str = [str,sprintf('\t%s: %0.4f\n',fn{kk},info.stats.(fn{kk}))];
                            end                                                                                                                    
                    case {'OneWayANOVA', 'TwoWayANOVA'}                        
                        str = [str, sprintf('Comparison performed: "%s"\nStats:\n',info.Name)];
                        str = [str, obj.tab2char(info.ANOVAtab)];
                        if ~isempty(info.postHocTab)
                            str = [str,sprintf('Post-hoc tests:\nTest name: "dunn-sidak"\nP-Value table:\n')];
                            str = [str,obj.tab2char(info.postHocTab)];
                        end                     
                    case {'OneWayRepeatedMeasures', 'TwoWayRepeatedMeasures'}                        
                        str = [str, sprintf('Comparison performed: "%s"\nStats:\n',info.Name)];
                        txt = obj.tab2char(info.ANOVAtab);
                        str = [str, txt];
                        if ~isempty(info.postHocTab)
                            str = [str,sprintf('Post-hoc tests:\nTest name: "tukey-kramer"\n%"s"\n', info.postHocTab.Name)];
                            str = [str,obj.tab2char(info.postHocTab.p)];
                        end             
                    otherwise
                        error('unknown hypothesis test!')                                      
               end
               str = [str, sprintf('%s\n',div)];
           end           
           
           % Save text to file:
           fid = fopen(filename,'w');
           fprintf(fid,'%s',str);           
           fclose(fid);
           disp(['Report saved to file: ' filename]);                                 
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    methods(Access = private)
        %%% Validator Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function validateObject(obj)
            % VALIDATEOBJECT checks if all objects in obj.list_of_objs
            % have obj.stats_filename in their respective SaveFolders.
            % Also, it checks if the sizes of list_of_objs and list_of_groups are equal.
            
            % If True, it creates an array of MAT-file objects and
            % saves to obj.MatFileList.
            % Else, it displays a warning and removes the object with the
            % missing obj.stats_filename from obj.list_of_objs.
            
            % Checks if list_of_objs have the same length as list_of_groups
            if size(obj.list_of_objs) ~= size(obj.list_of_groups)
                errID = 'umIToolbox:StatsManager:IncompatibleInputSizes';
                errMsg = ['The size of the list of objects is different from ' ...
                    'the size of list of groups'];
                error(errID, errMsg);
            end
            % Find Stats file of each object in list_of_objs:
            b_remObj = false(1,length(obj.list_of_objs));
            matHandleArr = cell(size(b_remObj));
            disp('Opening stats files...')
            for i = 1:length(obj.list_of_objs)
                elem = obj.list_of_objs{i};
                try
                    matHandleArr{i} = matfile(fullfile(elem.SaveFolder, obj.stats_filename));
                catch
                    warning(['Cannot find file "' obj.stats_filename '" in ' elem.SaveFolder ...
                        '. This object will be removed from the statistical group.'])
                    b_remObj(i) = true;
                end
            end
            disp('Done.')
            % Save MatFile handles:
            if all(b_remObj)
                % If none of the objects have the stats file, throw an
                % error:
                errID = 'umIToolbox:StatsManager:MissingData';
                errMsg = ['None of the objects listed contain "' obj.stats_filename '"'];
                error(errID, errMsg);
            else
                % Update lists:
                obj.MfileArr = matHandleArr(~b_remObj);
                obj.list_of_objs = obj.list_of_objs(~b_remObj);
                obj.list_of_groups = obj.list_of_groups(~b_remObj);
            end
        end
        
        function validateData(obj)
            % VALIDATEDATA performs a series of checks on the input data to
            % ensure that the StatsManager class works properly.
            % The input data must follow the constraints below:
            %   1) There must be an equal number of acquisitions per
            %   subject, per group.
            %   2) All input data must have the same dimension names.
            %   3) All input data must have the same dimension sizes.
            % In addition, it classifies the input data in one of the
            % following categories: 'scalar', 'time-vector','map', 'matrix' or 'unknown'.
            % And if it is "exportable" as a .CSV file. To be exportable,
            % the data must be scalar or be a 1xN array per observation.
            
            % Instantiate "inputFeatures" structure:
            obj.inputFeatures = struct('b_hasSameAcqN',false,'b_hasSameDimNames',false,...
                'b_AcqHasSameDimSize',false, 'b_SubjHasSameDimSize', false,  'b_isExportable',false,...
                'b_hasSameObs',false,'b_hasValidEvent',false, 'dataType', 'unknown', 'dim_names', {'unknown'});
            
            % Check #1 - are there is an equal number of acquisitions per
            % subject per group?
            gNames = unique(obj.list_of_groups);
            sNames = unique(obj.stats_data(:,obj.hMap('SubjectID')));
            acqIndx = {};
            for iG = 1:numel(gNames)
                idxG = strcmp(obj.stats_data(:,obj.hMap('groupID')), gNames{iG});
                for iS = 1:numel(sNames)
                    idxS = strcmp(obj.stats_data(:,obj.hMap('SubjectID')), sNames{iS});
                    acqIndx = [acqIndx; {sort([obj.stats_data{idxG & idxS, obj.hMap('AcquisitionIndx')}])}]; %#ok.
                end
            end
            acqIndx(cellfun(@isempty, acqIndx)) = [];
            obj.inputFeatures.b_hasSameAcqN = all(cellfun(@(x) all(isequal(x, acqIndx{1})), acqIndx));
            % Check #2 - Do the input data have the same dimensions?
            dim_names = arrayfun(@(x) x.MatFile.dim_names, obj.dataArr, 'UniformOutput',false);
            if isscalar(dim_names)
                obj.inputFeatures.b_hasSameDimNames = true;
            else
                obj.inputFeatures.b_hasSameDimNames = isequaln(dim_names{:});
            end
            % Check #3 - All observations exist across all data?
            obj.inputFeatures.b_hasSameObs = all(arrayfun(@(x) all(ismember(obj.obs_list, x.observationID)), obj.dataArr));
            % Check #4 - All data have "E"vents. If so, check if all necessary event info is present.
            if any(strcmpi(dim_names{1}, 'E'))
                obj.inputFeatures.b_hasValidEvent = all([obj.dataArr.b_hasEvents]);
            end
            % Check #5 - Do the input data have the same size across acquisitions per subject?
            %%% TO DO %%%
            sNames = unique(obj.stats_data(:, obj.hMap('SubjectID')));
            b_equalAcqSize = false(size(sNames));
            for ii = 1:numel(sNames)
                idxS = strcmp(obj.stats_data(:,obj.hMap('SubjectID')), sNames{ii});
                dim_sizes = obj.stats_data(idxS,obj.hMap('dataSize')); dim_sizes = vertcat(dim_sizes{:});
                if isscalar(dim_sizes)
                    b_equalAcqSize(ii) = true;
                else
                    b_equalAcqSize(ii) = isequaln(dim_sizes{:});
                end
            end
            if all(b_equalAcqSize)
                obj.inputFeatures.b_AcqHasSameDimSize = true;
            end
            % Check #6 - Do the input data have the same size across all
            % Subjects?
            if obj.inputFeatures.b_AcqHasSameDimSize
                dim_sizes = obj.stats_data(:,obj.hMap('dataSize')); dim_sizes = vertcat(dim_sizes{:});
                if isscalar(dim_sizes)
                    obj.inputFeatures.b_SubjHasSameDimSize = true;
                else
                    obj.inputFeatures.b_SubjHasSameDimSize = isequaln(dim_sizes{:});
                    obj.inputFeatures.b_SubjHasSameDimSize;
                end
            end
            % Check if the data is exportable to a .CSV file:
            % Here, two criteria are applied: 1) the data must be a 1xN
            % vector or a matrix; 2) the labels must be the same across
            obj.inputFeatures.b_isExportable = all(cellfun(@(x) isequaln(prod(x{1}),max(x{1})),...
                obj.stats_data(:,obj.hMap('dataSize'))));
            % Set input data type:
            if ~obj.inputFeatures.b_hasSameDimNames
                % If input data is heterogeneous, set data type as
                % "unknown".
                return
            else
                obj.inputFeatures.dim_names = obj.stats_data{1,obj.hMap('MatFile')}.dim_names;
            end
            % For homogeneous data, guess type of data from first element
            % of "stats_data" array:
            if sum(strcmp(obj.inputFeatures.dim_names, 'O')) == 2
                obj.inputFeatures.dataType = 'matrix'; % Correlation Matrix with dimensions {'O','O'};
            elseif ( isscalar(obj.stats_data{1,obj.hMap('data')}{1}) ||...
                    ( isequaln(prod(obj.stats_data{1, obj.hMap('dataSize')}{1}),max(obj.stats_data{1, obj.hMap('dataSize')}{1})) && ...
                    numel(setdiff(obj.inputFeatures.dim_names, {'E'})) == 1 ) )
                obj.inputFeatures.dataType = 'scalar'; % Single value per observation.
            elseif all(ismember(obj.inputFeatures.dim_names, {'O', 'T'})) || all(ismember(obj.inputFeatures.dim_names, {'O', 'T', 'E'}))
                obj.inputFeatures.dataType = 'time-vector';
            elseif all(ismember(obj.inputFeatures.dim_names, {'Y', 'X', 'O', 'E'}))
                obj.inputFeatures.dataType = 'map'; % Map with dimensions {'Y','X'} per observations.
                %             elseif all(ismember(obj.inputFeatures.dim_names, {'Y', 'X', 'T','O'}))
                %                 obj.inputFeatures.dataType = 'time-series';
                %             elseif all(ismember(obj.inputFeatures.dim_names, {'Y', 'X', 'T', 'E','O'}))
                %                 obj.inputFeatures.dataType = 'time-series-by-event';
                %             elseif all(ismember(obj.inputFeatures.dim_names, {'Y', 'X', 'E','O'}))
                %                 obj.inputFeatures.dataType = 'map-by-event';
            else
                obj.inputFeatures.dataType = 'unknown';
            end
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function createDataArray(obj)
            % CREATEDATAARRAY generates a cell array containing all
            % data and metadata provided as inputs to STATSMANAGER class.
            
            disp('Creating Data array. Please wait...')
            obj.headers = {'groupID', 'SubjectID', 'AcquisitionID', 'ModalityID',...
                'RecStartDateTime', 'MatFile', 'dataFile','labels','observationID',...
                'data', 'dataSize','AcquisitionIndx', 'b_isBaseline', 'b_hasEvents', 'indx_avg_data'};
            obj.stats_data = cell(numel(obj.list_of_objs),length(obj.headers));
            for i = 1:numel(obj.list_of_objs)
                obj.stats_data{i,1} = obj.list_of_groups{i}; % Group ID
                obj.stats_data{i,2} = obj.getElementInfo(obj.list_of_objs{i},'Subject', 'ID'); % Subject ID
                obj.stats_data{i,3} = obj.getElementInfo(obj.list_of_objs{i},'Acquisition', 'ID'); % Acquisition ID
                obj.stats_data{i,4} = obj.getElementInfo(obj.list_of_objs{i},'Modality', 'ID'); % modality ID
                obj.stats_data{i,5} = obj.getElementInfo(obj.list_of_objs{i},...
                    'Acquisition', 'Start_datetime'); % Acquisition start timestamp
                obj.stats_data{i,6} = obj.MfileArr{i}; % matfile handle
                obj.stats_data{i,7} = obj.MfileArr{i}.Properties.Source; % data file path
                obj.stats_data{i,8} = obj.MfileArr{i}.label; % data labels
                indx = find(ismember(obj.MfileArr{i}.obsID, obj.obs_list)); % Get observations listed in "obs_list" prop.
                obsID = obj.MfileArr{i}.obsID;
                data = obj.MfileArr{i}.data;
                obj.stats_data{i,9} = obsID(indx,:); % observation ID
                obj.stats_data{i,10}= data(indx,:); % observation data
                obj.stats_data{i,11}= cellfun(@(x) size(squeeze(x)), ...
                    obj.stats_data{i,10}, 'UniformOutput', false);% size of observation data
                obj.stats_data{i,obj.hMap('b_hasEvents')} = ( isprop(obj.MfileArr{i}, 'eventID') && isprop(obj.MfileArr{i}, 'eventNameList') );
            end
            
            % Create Relative time per subject's acquisitions per group:
            gNames= unique(obj.list_of_groups);
            for iG = 1:numel(gNames)
                idxG = strcmp(obj.stats_data(:,obj.hMap('groupID')),gNames{iG});
                subjs = unique(obj.stats_data(:,obj.hMap('SubjectID')));
                for iS = 1:numel(subjs)
                    indxS = find(strcmp(subjs{iS}, obj.stats_data(:,obj.hMap('SubjectID'))) & idxG);
                    acq_time_list = datetime(string(obj.stats_data(indxS,obj.hMap('RecStartDateTime'))));
                    [~,tm_idx] = sort(acq_time_list);
                    [~,rel_time]= sort(tm_idx);
                    obj.stats_data(indxS,obj.hMap('AcquisitionIndx')) = arrayfun(@(x) x, rel_time, 'UniformOutput', false);
                    obj.stats_data(indxS,obj.hMap('b_isBaseline')) = num2cell([obj.stats_data{indxS,obj.hMap('AcquisitionIndx')}] == 1);
                end
            end
            obj.stats_data(:,obj.hMap('indx_avg_data')) = num2cell(zeros(size(obj.stats_data,1),1, 'single'));
            disp('Done!')
        end
        
        function out = getElementInfo(~,elem,className,propName)
            % getElementInfo retrieves the property PROPNAME from an
            % object ELEM of class CLASSNAME or from its Parents.
            
            % Find the object of class CLASSNAME:
            b_isClass = false;
            while ~b_isClass
                if isa(elem, className)|| isa(elem, 'Protocol')
                    b_isClass = true;
                else
                    elem = elem.MyParent;
                end
            end
            % Get PROPNAME VALUE
            if isprop(elem, propName)
                out = elem.(propName);
            else
                out = [];
            end
        end
        
        function out = getObsData(obj,obs_ID,labelList)
            % This function retrieves the observation data from all
            % recordings containing the observation "obs_ID".
            % Input
            %   obs_ID (char) = Name of observation.
            %   labelList(cell array of char): list of all possible labels in
            %       "dataArr".
            % Outputs:
            %   out (table) : table containing the recording
            %       information and numerical data of the observation
            %       "obs_ID" for all recordings.
            
            % Find indexes of the observation inside dataArr:
            indx_obs = zeros(1,length(obj.dataArr));
            for i = 1:length(indx_obs)
                indx = find(strcmp(obs_ID, obj.dataArr(i).observationID));
                if isempty(indx)
                    continue
                end
                indx_obs(i) = indx;
            end
            % Get info only from data containing the observation:
            stats_info = obj.dataArr(indx_obs ~= 0);
            indx_obs = indx_obs(indx_obs~=0);
            % Get observation's data:
            % Preallocate data with NaNs based on the length of unique
            % labels:
            data = nan(length(stats_info),length(labelList));
            for i = 1:length(stats_info)
                [~,locB] = ismember(stats_info(i).labels,labelList);
                data(i,locB) = stats_info(i).data{indx_obs(i)};
            end
            % Prepare data to be added to "out" table:
            data = num2cell(data);
            % Prepare metaData to be added to "out" table:
            % Build table:
            out = table('Size', [length(stats_info), 6 + size(data,2)],...
                'VariableTypes',...
                [{'cellstr', 'cellstr', 'cellstr', 'cellstr', 'datetime', 'cellstr'}...
                repmat({obj.MfileArr{1}.Datatype},1,size(data,2))]);
            out(:,1:5) = [{stats_info.groupID}', {stats_info.SubjectID}', {stats_info.AcquisitionID}', ...
                {stats_info.ModalityID}' arrayfun(@(x) datetime(x.RecStartDateTime), stats_info, 'UniformOutput',false)'];
            out(:,6) = repmat({obs_ID},length(stats_info),1);
            out(:,7:end) = data;
            % Set table header:
            generic_label = arrayfun(@(x) ['Val_' num2str(x)], 1:size(data,2), 'UniformOutput',false);
            out.Properties.VariableNames = [{'GroupName','Subject', 'Acquisition', 'Modality',...
                'RecordingTime', 'ObsID'}, generic_label];
        end
        
        function getEventList(obj)
            % GETEVEENTLIST creates a list of unique values from "eventNameList" variable
            % inside the datas' meta data file. If no events exist, the list
            % will contain the string "NoEvents".
            
            if isempty(obj.stats_data)
                error('Failed to get event list.The property "stats_data" is empty!')
            end
            idxEv = obj.stats_data{:,obj.hMap('b_hasEvents')};
            if ~any(idxEv)
                obj.list_of_events = {'NoEvents'};
                return
            end
            all_evNames = cellfun(@(x) x.eventNameList, obj.stats_data(idxEv,obj.hMap('MatFile')), 'UniformOutput', false)';
            obj.list_of_events = unique(vertcat(all_evNames{:}));
            
            
        end
        
        function out = hMap(obj,colName)
            % HMAP gives the index of the "stats_data" column with
            % the header "colName".
            if ischar(colName)
                colName = {colName};
            end
            [~, out]= ismember(lower(colName), lower(obj.headers));
        end
        
        %%%%%%%%%%%%  Auxiliary Stats functions %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function out = calculateVariation(~,varType,data,dim)
            % CALCULATEVARIATION calculates one of the following variation
            % measures on "data":
            %   'STD': standard deviation
            %   'SEM': standard error of the mean
            %   'CI': 95% Confidence Interval.
            % Inputs:
            %   varType (char): name of the variation measure (see above).
            %   data (numerical array): data to calculate variation.
            %   dim (int scalar): dimension to calculate the variation.
            % !! Note #1: NaNs are excluded in all calculations.
            % Output:
            %   out(numerical array): variation of the input data. For the
            %   CI measure, the upper bound is provided.
            
            if strcmpi(varType, 'ci') && ~obj.b_hasStatsToolbox
                error('Cannot calculate CI. Statistics and Machine Learning Toolbox is needed!')
            end
            
            switch upper(varType)
                case 'STD'
                    out = std(data,0,dim, 'omitnan');
                case 'SEM'
                    out = std(data,0,dim,'omitnan')./sqtr(size(data,dim));
                case 'CI'
                    out = (tinv(1-0.025,size(data,dim)-1)).*(std(data,0,dim,'omitnan')./sqtr(size(data,dim)));
                otherwise
                    error('Unknown variation measure!')
            end
        end
        
        function chooseStatsTest(obj)
            % This method uses the current independent variables and the
            % size of the data in "dataArr" to decide which statistical test
            % to apply.
            
            obj.flatDataArr; % Prepare data for stats.
            disp('Choosing statistical test...')
            % Get first independent variable:
            obj.indepVars(1).len = length(unique(obj.flatData.(obj.indepVars(1).fieldName)));
            % Secondary independent variable:
            obj.indepVars(2).len = length(unique(obj.flatData.(obj.indepVars(2).fieldName)));
            % Decide which test to use:
            % One group: abort:
            disp(repmat('-',1,100));
            fprintf('The data is grouped by %d "%s"(s) and %d "%s"(s) across %d "%s"(s).\n',...
                obj.indepVars(1).len, obj.indepVars(1).Name, obj.indepVars(2).len, obj.indepVars(2).Name,...
                numel(unique(obj.flatData.(obj.splitVar.indexName))), obj.splitVar.Name)
            if obj.indepVars(1).len == 1 && obj.indepVars(2).len == 1
                disp('No stats for one group with one independent sample')
                disp(repmat('-',1,100));
                return
            end
            if ( obj.indepVars(1).len == 2 && obj.indepVars(2).len == 1 && obj.indepVars(1).b_isRepeatedMeasure ) ||...
                    ( obj.indepVars(1).len == 1 && obj.indepVars(2).len == 2 && obj.indepVars(2).b_isRepeatedMeasure )
                % Paired two-sided test                
                obj.curr_test = 'PairedTest';
            elseif ( obj.indepVars(1).len == 2 && obj.indepVars(2).len == 1 && ~obj.indepVars(1).b_isRepeatedMeasure ) ||...
                    ( obj.indepVars(1).len == 1 && obj.indepVars(2).len == 2 && ~obj.indepVars(2).b_isRepeatedMeasure )
                % Unpaired two-sided test                
                obj.curr_test = 'UnpairedTest';
            elseif ( obj.indepVars(1).len == 1 && obj.indepVars(2).len > 2 && ~obj.indepVars(2).b_isRepeatedMeasure ) ||...
                    ( obj.indepVars(1).len > 2 && obj.indepVars(2).len == 1 && ~obj.indepVars(1).b_isRepeatedMeasure )
                % 3+ independent (unpaired) samples test                
                obj.curr_test = 'OneWayANOVA';
            elseif ( all(~[obj.indepVars.b_isRepeatedMeasure]) ) && ( ( obj.indepVars(1).len > 1 && obj.indepVars(2).len > 2 ) ||...
                    ( obj.indepVars(1).len > 2 && obj.indepVars(2).len > 1 ) )
                % 3+ independent (unpaired) samples test                
                obj.curr_test = 'TwoWayANOVA';
            elseif ( obj.indepVars(1).len == 1 && obj.indepVars(2).len > 2 && obj.indepVars(2).b_isRepeatedMeasure && ~obj.indepVars(1).b_isRepeatedMeasure ) ||...
                    ( obj.indepVars(1).len > 2 && obj.indepVars(2).len == 1 && obj.indepVars(1).b_isRepeatedMeasure && ~obj.indepVars(2).b_isRepeatedMeasure )
                % 3+ dependent (repeated measures) samples test                
                obj.curr_test = 'OneWayRepeatedMeasures';
            elseif ( (obj.indepVars(1).len > 1 && obj.indepVars(2).len >= 2) && obj.indepVars(1).b_isRepeatedMeasure ) || ...
                    ( (obj.indepVars(1).len >= 2 && obj.indepVars(2).len > 1) && obj.indepVars(2).b_isRepeatedMeasure )                    
                % 3+ dependent (repeated measures) with 2+ experimental groups
                obj.curr_test = 'TwoWayRepeatedMeasures';
            else
                obj.curr_test = 'Unknown';
            end
            % Update "flatData" with the list of names from stats grouping variables:
            
            % Get the list of unique values for each independent variable:
            obj.flatData.indPNames = cell(2,1);
            for ii = 1:2
                % Get the names of the variable keeping its original
                % order:
                names = unique(obj.flatData.(obj.indepVars(ii).fieldName), 'stable');
                if size(names,1) > size(names,2)
                    names = names';
                end
                obj.flatData.indPNames{ii} = names;
            end
            % Get the list of names of the split variable:
            obj.flatData.splitNames = obj.flatData.(obj.splitVar.fieldName);
            % Display messages:
            fprintf('Hypothesis test chosen: %s.\n', obj.curr_test);
            disp(repmat('-',1,100));
            disp('Done')
        end
        
        function runHypothesisTest(obj)
            % This method performs the hypothesis test between groups using
            % the information created by the method "chooseStatsTest".
            % The results will be stored in the property "results_stats".
            % The available tests are listed in the property
            % "list_of_tests".
            %   **If a test is added or removed from this
            %   `method, please update the "list_of_tests" array**.
            
            disp('Performing hypothesis tests...')
            % For each subset of groups in "splitVar", perform the
            % statistical comparison:
            indxG = unique(obj.flatData.(obj.splitVar.indexName));
            % Instantiate stats structure:
            statsInfo = repmat(struct('subGroupName','','GroupComparison', struct('Name','None', 'h',[],'p',[],'stats',[], 'ANOVAtab',[],...
                'postHocTab',[])),1,numel(indxG));
            
            for ii = 1:numel(indxG)
                statsInfo(ii).subGroupName = obj.flatData.(obj.splitVar.fieldName){ii};
                % Perform stats in data subset:
                idxSplit = obj.flatData.(obj.splitVar.indexName) == indxG(ii);
                
                % Perform stats:
                switch obj.curr_test
                    case 'PairedTest'
                        % Paired test with 2 groups:
                        idxIndep = [obj.indepVars.len] > 1;
                        dv = dummyvar(obj.flatData.(obj.indepVars(idxIndep).indexName));
                        % Segregate data into two subsamples:
                        x = obj.flatData.data(idxSplit & dv(:,1));
                        y = obj.flatData.data(idxSplit & ~dv(:,1));
                        if ( isempty(x) || isempty(y) ) || ( ~isequal(length(x),length(y)) )
                            % Skip if one of the subsamples are missing:
                            warning(['Skipped Paired test because of missing data on comparison :' statsInfo(ii).subGroupName ' --> ' statsInfo(ii).GroupComparison.Name]);
                            continue
                        end
                        % Pair data across the independent variable with
                        % repeated measures making sure that the other
                        % variables are ordered:
                        otherVars = obj.indepVarInfo(~strcmpi(obj.indepVars(idxIndep).Name, {obj.indepVarInfo.Name}));
                        indxList = arrayfun(@(x) obj.flatData.(x.indexName), otherVars, 'UniformOutput', false);
                        indxList = horzcat(indxList{:});
                        indxX = indxList(idxSplit & dv(:,1),:); indxY = indxList(idxSplit & ~dv(:,1),:);
                        for jj = 1:size(indxX,2)
                            % For X:
                            [~, indxSort] = sort(indxX(:,jj));
                            x = x(indxSort);
                            indxX = indxX(indxSort,:);
                            % For Y:
                            [~, indxSort] = sort(indxY(:,jj));
                            y = y(indxSort);
                            indxY = indxY(indxSort,:);
                        end
                        clear indxSort indxX indxY otherVars
                        
                        % Create ID for comparison:
                        statsInfo(ii).GroupComparison.Name = strjoin(obj.flatData.(obj.indepVars(idxIndep).fieldName), '_vs_');                        
                        % Perform statistical comparison:
                        if obj.checkNormality
                            % Parametric T-Test:
                            [statsInfo(ii).GroupComparison.h,statsInfo(ii).GroupComparison.p,~,statsInfo(ii).GroupComparison.stats] = ttest(x,y,'Alpha',obj.pAlpha);
                        else
                            % Non-parametric Wilcoxon rank-sum test:
                            [statsInfo(ii).GroupComparison.p,statsInfo(ii).GroupComparison.h,statsInfo(ii).GroupComparison.stats] = signrank(x,y,'alpha',obj.pAlpha);
                        end
                        
                    case 'UnpairedTest'
                        % Unpaired test with 2 groups:
                        % Segregate data into two subsamples:
                        idxIndep = [obj.indepVars.len] > 1;
                        dv = dummyvar(obj.flatData.(obj.indepVars(idxIndep).indexName));
                        x = obj.flatData.data(idxSplit & dv(:,1));
                        y = obj.flatData.data(idxSplit & ~dv(:,1));
                        % Create ID for comparison:
                        statsInfo(ii).GroupComparison.Name = strjoin(obj.flatData.(obj.indepVars(idxIndep).fieldName), '_vs_');
                        if isempty(x) || isempty(y)
                            % Skip if one of the subsamples are missing:
                            warning(['Skipped UnPaired Wilcoxon test because of missing data on comparison :' statsInfo(ii).subGroupName ' --> ' statsInfo(ii).GroupComparison.Name]);
                            continue
                        end
                        % Compensate for unbalanced variances in T test:
                        varMap = containers.Map([true,false], {'equal', 'unequal'});
                        % Perform statistical comparison:
                        if obj.checkNormality
                            % Parametric T-Test:
                            [statsInfo(ii).GroupComparison.h,statsInfo(ii).GroupComparison.p,~,statsInfo(ii).GroupComparison.stats] = ttest2(x,y,'Alpha',obj.pAlpha,'Vartype',varMap(obj.checkHomogeneity));
                        else
                            % Non-parametric Wilcoxon rank-sum test:
                            [statsInfo(ii).GroupComparison.p,statsInfo(ii).GroupComparison.h,statsInfo(ii).GroupComparison.stats] = ranksum(x,y,'Alpha',obj.pAlpha);
                        end
                    case 'OneWayANOVA'                        
                        idxHas2PlusItems = [obj.indepVars.len] > 1;
                        % Segregate data into 3+ subsamples:
                        gVar = obj.flatData.(obj.indepVars(idxHas2PlusItems).indexName);
                        gNames = obj.flatData.(obj.indepVars(idxHas2PlusItems).fieldName);
                        gMap = containers.Map(unique(gVar),gNames);
                        gVar = arrayfun(@(x) gMap(x),gVar, 'UniformOutput',false);
                        x = obj.flatData.data(idxSplit);
                        % Create ID for comparison:
                        statsInfo(ii).GroupComparison.Name = strjoin(obj.flatData.(obj.indepVars(idxHas2PlusItems).fieldName), '_vs_');
                        if obj.checkHomogeneity
                            % Perform one-way ANOVA:
                            [statsInfo(ii).GroupComparison.p,statsInfo(ii).GroupComparison.ANOVAtab, statsInfo(ii).GroupComparison.stats] = anova1(x, gVar(idxSplit), 'off');
                        else
                            % Perform Kruskal-Wallis:
                            [statsInfo(ii).GroupComparison.p,statsInfo(ii).GroupComparison.ANOVAtab, statsInfo(ii).GroupComparison.stats] = kruskalwallis(x, gVar(idxSplit), 'off');
                        end                        
                        % Perform postHoc test when a difference was detected in the ANOVA test:
                        if statsInfo(ii).GroupComparison.p < obj.pAlpha
                            % Prepare identifiers for pair-wise group comparison:                            
                            Names = cell(length(gNames),1);                            
                            c = multcompare(statsInfo(ii).GroupComparison.stats,'Alpha', obj.pAlpha,'CType', 'dunn-sidak', 'Display','off');                            
                            for kk = 1:size(c,1)
                                Names{kk}= strjoin([gNames(c(kk,1)), gNames(c(kk,2))], '_vs_');                                
                            end
                            statsInfo(ii).GroupComparison.postHocTab = table(Names,c(:,end),'VariableNames',{'ComparisonID','p_value'});                            
                        end
                    case 'TwoWayANOVA'                        
                        % Create groups for anovan function
                        indxVar1 = obj.flatData.(obj.indepVars(1).indexName)(idxSplit);
                        g1 = obj.flatData.(obj.indepVars(1).fieldName)(indxVar1);
                        indxVar2 = obj.flatData.(obj.indepVars(2).indexName)(idxSplit);
                        g2 = obj.flatData.(obj.indepVars(2).fieldName)(indxVar2);
                        if size(g1)~=size(g2)
                            g1 = g1';
                        end
                        data = obj.flatData.data(idxSplit);
                        % Create GENERIC ID for comparison:
                        statsInfo(ii).GroupComparison.Name = strjoin({obj.indepVars.Name},'_vs_');
                        % Perform ANOVA:
                        [statsInfo(ii).GroupComparison.p,statsInfo(ii).GroupComparison.ANOVAtab, statsInfo(ii).GroupComparison.stats] = ...
                            anovan(data,{g1,g2}, 'model','interaction', 'varnames',{obj.indepVars.Name}, 'display','off');                                                                                        
                        % Perform postHoc test when a difference was detected in the ANOVA test:
                        % Create table with postHoc info:                        
                        if any(statsInfo(ii).GroupComparison.p < obj.pAlpha)
                            % Prepare identifiers for pair-wise group comparison:
                            g1Names = unique(g1); g2Names = unique(g2);
                            [m,n] = meshgrid([1:length(g1Names)],[1:length(g2Names)]); gPairs = [m(:),n(:)];
                            gPairIDList = cell(size(gPairs,1),1);
                            for kk = 1:size(gPairs,1)
                                gPairIDList{kk} = strjoin([g1Names(gPairs(kk,1)), g2Names(gPairs(kk,2))],',_');
                            end
                            c = multcompare(statsInfo(ii).GroupComparison.stats,'Alpha', obj.pAlpha,'CType', 'dunn-sidak', 'Display','off', 'Dimension',[1 2]);                            
                            Names = cell(size(c,1),1);
                            for kk = 1:size(c,1)
                                Names{kk} = strjoin([gPairIDList(c(kk,1)), gPairIDList(c(kk,2))], '_vs_');                                
                            end
                            statsInfo(ii).GroupComparison.postHocTab = table(Names,c(:,end),'VariableNames',{'ComparisonID','p_value'});                            
                        end                                               
                    case {'OneWayRepeatedMeasures', 'TwoWayRepeatedMeasures'}                        
                        % Get independent variable without repeated
                        % measures:
                        idxPred = obj.flatData.(obj.indepVars([obj.indepVars.b_isRepeatedMeasure] == 0).indexName)(idxSplit);   
                        predNames = obj.flatData.(obj.indepVars([obj.indepVars.b_isRepeatedMeasure] == 0).fieldName)(idxPred);   
                        % Get Within-Subject indices and names:
                        idxWithin = obj.flatData.(obj.indepVars([obj.indepVars.b_isRepeatedMeasure] == 1).indexName)(idxSplit);                       
                        withinList = unique(idxWithin);
                        withinVarName = obj.indepVars([obj.indepVars.b_isRepeatedMeasure] == 1).Name;
                        betweenVarName = obj.indepVars([obj.indepVars.b_isRepeatedMeasure] == 0).Name;
                        withinNames = matlab.lang.makeValidName(obj.flatData.(obj.indepVars([obj.indepVars.b_isRepeatedMeasure] == 1).fieldName)(withinList));
                        % Create GENERIC ID for comparison:
                        statsInfo(ii).GroupComparison.Name = [betweenVarName ' by ' withinVarName];
                        % Get list of subject indices:
                        idxSubj = obj.flatData.sIndx(idxSplit);                        
                        subjList = unique(idxSubj);
                        % Create tables for RepeatedMeasures Model fit:
                        dataIn = obj.flatData.data(idxSplit);
                        dataOut = nan(numel(subjList),numel(withinList));                                                
                        predVar = cell(size(dataOut,1),1);
                        % 
                        for jj = 1:length(subjList)
                            idx1 = idxSubj == subjList(jj);                            
                            for kk = 1:length(withinList)
                                idx2 = idxWithin == withinList(kk);
                                predVar(jj) = predNames(idx1 & idx2);
                                dataOut(jj,kk) = dataIn(idx1 & idx2);
                            end                            
                        end
                        dataOut = num2cell(dataOut); 
                        t = cell2table([predVar,dataOut],'VariableNames',[{betweenVarName}, withinNames]);
                        
                        if contains(obj.curr_test, 'oneway','IgnoreCase',true)
                            notation = [withinNames{1} '-' withinNames{end} ' ~ 1'];
                            postHocVar = withinVarName;
                        else
                            notation = [withinNames{1} '-' withinNames{end} ' ~ ' t.Properties.VariableNames{1}];
                            postHocVar = betweenVarName;
                        end
                        rm = fitrm(t,notation,'WithinDesign',table(withinList,'VariableNames',{withinVarName}));
                        statsInfo(ii).GroupComparison.ANOVAtab = ranova(rm);
                        pv_col = startsWith(statsInfo(ii).GroupComparison.ANOVAtab.Properties.VariableNames, 'pValue','IgnoreCase',true);
                        pvals = cellfun(@(x) (x < obj.pAlpha), table2cell(statsInfo(ii).GroupComparison.ANOVAtab(:, pv_col)));
                        if any(pvals,'all')
                            % Post-hoc test:
                            statsInfo(ii).GroupComparison.postHocTab.Name = ['Comparison by ' t.Properties.VariableNames{1}];
                            statsInfo(ii).GroupComparison.postHocTab.p = multcompare(rm,postHocVar);
                        end
                    otherwise
                        error('unknown hypothesis test!')
                end
                obj.results_stats = statsInfo;  
            end
            disp('Done')
        end
        
        function status = checkNormality(obj)
            % Check for normality
            % Inputs:
            % idx (array | matrix): array of stats. groups indices to test.
            % CheckName (str): {'normality', 'homoskedacity','sphericity'}
            indx = bin2dec([dec2bin(obj.flatData.(obj.splitVar.indexName)),...
                dec2bin(obj.flatData.(obj.indepVars(1).indexName)), ...
                dec2bin(obj.flatData.(obj.indepVars(2).indexName))]);
            indxList = unique(indx);
            status = false(size(indxList,1),1);
            for ii = 1:size(indxList,1)
                data = obj.flatData.data(indx == indxList(ii,:));
                if numel(data) >= 4
                    status(ii) = lillietest(data);
                end
            end
            status = all(~status);
        end
        
        function status = checkHomogeneity(obj)
            % Check for homogeneity of variances for ANOVA
            % Here, use the 4:1 ratio as a maximum disparity of variances.
            % For Repeated measures, sphericity WILL NOT BE TESTED HERE!
            status = true;
            if any(strcmpi(obj.curr_test, 'PairedTest'))
                return
            end
            
            idxG = obj.flatData.(obj.splitVar.indexName);gList = unique(idxG);
            status = false(length(gList),1);
            if strcmpi(obj.curr_test, 'OneWayANOVA')
                % For one way ANOVA:
                idxVar1 = obj.flatData.(obj.indepVars([obj.indepVars.len] == 1).indexName);varList1 = unique(idxVar1);
                idxVar2 = obj.flatData.(obj.indepVars([obj.indepVars.len]>1).indexName);varList2 = unique(idxVar2);
            elseif strcmpi(obj.curr_test, 'TwoWayANOVA')
                % For two way ANOVA:
                % Use the order of variables in "indepVars".
                idxVar1 = obj.flatData.(obj.indepVars(1).indexName); varList1 = unique(idxVar1);
                idxVar2 = obj.flatData.(obj.indepVars(2).indexName); varList2 = unique(idxVar2);
            else
                idxVar1 = obj.flatData.(obj.indepVars(~[obj.indepVars.b_isRepeatedMeasure]).indexName);varList1 = unique(idxVar1);
                idxVar2 = obj.flatData.(obj.indepVars([obj.indepVars.b_isRepeatedMeasure]).indexName);varList2 = unique(idxVar2);
            end
            for ii = 1:length(gList)
                for jj = 1:length(varList1)
                    groupVars = zeros(length(varList2),1);
                    for kk = 1:length(varList2)
                        groupVars(kk) = var(obj.flatData.data(idxG == gList(ii) & idxVar1 == varList1(jj) & idxVar2 == varList2(kk)));
                    end
                    if min(groupVars)/max(groupVars) > 0.25
                        status(ii) = true;
                    end
                end
            end
            status = all(status);
        end
        
        function out = tab2char(obj,tab)
            % Transforms a table as text.           
            % Parse table:
            if iscell(tab)
                % Deal with info already in "cells"                 
                txt = cellfun(@num2str,tab, 'UniformOutput',false);
            else                
                % Transform table to cell array with strings:
                txt = table2cell(tab);
                txt = cellfun(@num2str,txt,'UniformOutput',false);
                % Append headers and row names:
                txt = [tab.Properties.RowNames,txt];
                if ~isempty(tab.Properties.RowNames)
                    txt = vertcat([{''}, tab.Properties.VariableNames], txt);
                else
                    txt = vertcat(tab.Properties.VariableNames, txt);
                end
            end
            % Get largest number of charaters:
            charLen = max(cellfun(@length,txt),[],'all');
            txt = cellfun(@(x) pad(x,charLen,'right'),txt,'UniformOutput',false);
            % Build table:
            divV = ' | ';
            divH = repmat('-',1,(charLen + length(divV))*size(txt,2));
            
            out = '';
            for ii = 1:size(txt,1)
                out = [out,sprintf('%s\n',strjoin(txt(ii,:),divV))];
                if ii == 1 || ii == size(txt,1)
                    out =[out,sprintf('%s\n',divH)];
                end
            end            
        end
    end
end











































