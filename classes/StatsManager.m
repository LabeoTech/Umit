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
        obs_list       % Cell Array of observations from stats_filename.
        list_of_groups % Cell array of group names of objects from list_of_objs.
        list_of_events % Cell array of event Names from all objects.
        % !!If the input data has no events, the list will contain the string"NoEvents" !!
        stats_filename % Filename of .MAT file saved using "save2Mat.m"
        % function.
        MfileArr       % Array of MATFILE objects.
        time_resolution = 'none' % Time resolution for grouping each observation
        % Options: "none", "minute", "hour, "day", "week", "month".
    end
    properties (SetAccess = private)
        %         timestamp_list % List of timestamps associated with each object in list_of_objs.
        b_hasStatsToolbox  % True, if Matlab contains the Statistics and Machine learning toolbox.
        inputFeatures % Structure containing some information about the input data. This will be used by plotting tools and the umIToolbox app.
        dataArr  % structure with data extracted from "stats_data" that can be averaged using the method "setAcquisitionRange".
    end
    properties (Access = private)
        stats_data  = {} % cell array containing all data and metaData created.
        avg_stats_data = {} % cell array similar to stats_data containing the average values of acquisitions. (See method averageData).
        headers = {} % cell array with the _stats_data column names as keys and indices as values.
    end
    
    methods
        function obj = StatsManager(list_of_objs, obs_list, list_of_groups, stats_filename)
            % Class constructor.
            %   This function initiates the object "StatsManager" with the
            %   properties: list_of_objs, obs_list, list_of_groups and stats_filename.
            %   All inputs must be provided.
            obj.list_of_objs = list_of_objs;
            obj.obs_list = obs_list;
            obj.list_of_groups = list_of_groups;
            obj.stats_filename = stats_filename;
            % Check if inputs are correct and Generate list of MatFile
            % handles:
            obj.validateObject;
            obj.createDataArray;
            obj.getEventList; % Get list of events and store it to "list_of_events" property.
            obj.update_dataArr;
            obj.validateData; % Identifies the type of input data and other properties.
            
            % Check if Matlab's Stats. toolbox exist:
            a = ver;
            if any(strcmp('Statistics and Machine Learning Toolbox', {a.Name}))
                obj.b_hasStatsToolbox = true;
            else
                obj.b_hasStatsToolbox = false;
            end
        end
        %%% Property Set Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set.list_of_objs(obj, list_of_objs)
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
        
        function set.obs_list(obj, obs_list)
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
        
        function set.list_of_groups(obj, list_of_groups)
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
        
        function set.stats_filename(obj, stats_filename)
            % Set function for stats_filename.
            % It checks if stats_filename is a .MAT file.
            errID = 'umIToolbox:StatsManager:InvalidInput';
            errMsg = 'Wrong input. Stats file must be a .MAT file';
            assert(isa(stats_filename, 'char') & endsWith(stats_filename, '.mat'),...
                errID, errMsg);
            obj.stats_filename = stats_filename;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [out,uniqLabels] = createTable(obj, varargin)
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
        
        function exportToCSV(obj, filename)
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
        
        function out = getAcqIndexList(obj, type)
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
                error('Input should be either "original" or "available".')
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
        
        function setAcquisitionRange(obj, indxRange)
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
            obj.update_dataArr;
        end
        
        function resetAvgIndex(obj)
            % RESETAVGINDEX sets all acquisition indices in the average
            % data to zero.
            obj.stats_data(:,obj.hMap('indx_avg_data')) = repmat({0},size(obj.stats_data,1),1);
            disp('Acquisition average group indices reset!')
            obj.update_dataArr;
        end
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
        
        function out = getElementInfo(~,elem,className, propName)
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
        
        function out = getObsData(obj, obs_ID, labelList)
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
        
        function out = hMap(obj, colName)
            % HMAP gives the index of the "stats_data" column with
            % the header "colName".
            if ischar(colName)
                colName = {colName};
            end
            [~, out]= ismember(lower(colName), lower(obj.headers));
        end
        
        function update_dataArr(obj)
            % UPDATE_DATAARR updates the values of dataArr property. This
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
            disp('Done');
            % Add some indices to obj.dataArr:
            for ind = 1:length(obj.dataArr)
                obj.dataArr(ind).gIndx = find(strcmp(obj.dataArr(ind).groupID, unique(obj.list_of_groups)));
                obj.dataArr(ind).sIndx = find(strcmp(obj.dataArr(ind).SubjectID, unique({obj.dataArr.SubjectID})));
                [~, obj.dataArr(ind).rIndx] = ismember(obj.dataArr(ind).observationID, obj.obs_list);
                if obj.dataArr(ind).b_hasEvents
                    % Flag if events exist:
                    [~,obj.dataArr(ind).eIndx] = ismember(obj.dataArr(ind).MatFile.eventNameList, obj.list_of_events);
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
                    obj.dataArr(ind).eIndx = 0;
                end
            end
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
                    'dataSize', 'MatFile'};
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
        
        %%%%%%%%%%%%  Auxiliary Stats functions %%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = calculateVariation(~, varType,data, dim)
            % CALCULATEVARIATION calculates one of the following variation
            % measures on "data":
            %   'STD': standard deviation
            %   'SEM': standard error of the mean
            %   'CI': 95% Confidence Interval.
            % Inputs:
            %   varType (char): name of the variation measure (see above).
            %   data (numerical array): data to calculate variation.
            %   dim (int scalar): number of the dimension to calculate the
            %   variation.
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
        
        %%%%%%%%%%%%%%%%%%
    end
end
