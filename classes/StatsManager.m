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
        stats_filename % Filename of .MAT file saved using "save2Mat.m"
                       % function. 
        MfileArr       % Array of MATFILE objects.
        time_resolution = 'none' % Time resolution for grouping each observation 
                       % Options: "none", "minute", "hour, "day", "week", "month".                  
    end
    properties (SetAccess = private)
        timestamp_list % List of timestamps associated with each object in list_of_objs.
%         stats_data     = struct() % Structure containing all data and metadata created.
        stats_data  = {} % cell array containing all data and metaData created.
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
            obj.validateStatsInputs;            
            obj.createDataArray;
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
                errID = 'Umitoolbox:StatsManager:WrongInput';
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
                errID = 'Umitoolbox:StatsManager:WrongInput';
                errMsg = 'List of groups must be a non-empty cell array of characters';
                error(errID, errMsg);
            end
        end
            
        function set.stats_filename(obj, stats_filename)
            % Set function for stats_filename.
            % It checks if stats_filename is a .MAT file.
            errID = 'Umitoolbox:StatsManager:InvalidInput';
            errMsg = 'Wrong input. Stats file must be a .MAT file';
            assert(isa(stats_filename, 'char') & endsWith(stats_filename, '.mat'),...
                errID, errMsg);
            obj.stats_filename = stats_filename;
        end
        
        function set.time_resolution(obj, time_res)
            % Set function for time_resolution.
            % It checks if time_res is a valid string.
            errID = 'Umitoolbox:StatsManager:InvalidInput';
            errMsg = ['Invalid time resolution. Valid options are: ' ...
                '{"none", "second", "minute", "hour", "day", "week", "month"}'];
            assert(ismember(time_res, {'none', 'minute', 'hour', 'day', ...
                'week', 'month'}), errID, errMsg);
            % Set time_res:
            obj.time_resolution = time_res;
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
            
            % Get all unique labels from "stats_data" structure. This will
            % be used by "getObsData" method to put NaNs on missing
            % data (i.e. labels that are missing for a given recording).
            
            labels = arrayfun(@(x) x.MatFile.label, obj.stats_data, 'UniformOutput',false);
            labels = labels'; 
            uniqLabels = unique(vertcat(labels{:})); % Sort unique labels in alphabetical order.
                        
            % Get observations's labels from the stats_data structure:
                
            if strcmp(tableType, 'raw')                                
                tableArr = cellfun(@(x) obj.getObsData(x, uniqLabels), obj.obs_list,'UniformOutput',false);
                out = vertcat(tableArr{:});
            else
                disp('No!')
            end            
            disp('Table created!')
        end
        
        function out = packageData(obj)
            % This function creates a structure with all
            % "raw" stats data.
            % Output:
            %   out(struct): structure containing all data and metadata stored
            %   in obj.stats_data.
            
            % Get metaData from Matfile:
            fields = {'groupID', 'SubjectID', 'AcquisitionID', 'ModalityID',...
                'RecStartDateTime', 'MatFile', 'dataFile','labels','observationID',...
                'data','AcquisitionIndx'};
            out = cell2struct(obj.stats_data, fields,2);
            out(1).metaData = [];
            h = waitbar(0,'Compiling stats data...');
            for i = 1:length(out)
                metaData_fn = properties(out(i).MatFile);
                metaData_fn = setdiff(metaData_fn, {'Properties','data',...
                    'obsID','label', 'datFile','datLength','datSize'});
                for j = 1:numel(metaData_fn)
                    out(i).metaData.(metaData_fn{j}) = out(i).MatFile.(metaData_fn{j});
                end
                % Pack observation info:
                out(i).observations = struct('ID','', 'data',[]);
                for j = 1:numel(out(i).observationID)
                    out(i).observations(j).ID = out(i).observationID{j};
                    out(i).observations(j).data = out(i).data{j};
                end
                waitbar(i/length(out),h);
            end
            out = rmfield(out, {'MatFile', 'observationID', 'data'});            
            waitbar(1,'Done!');
            pause(1);
            close(h);
        end
        function exportToCSV(obj, filename)
            % This function creates a .CSV file containing all data created
            % by the method "createTable".
            % Input:
            % filename (char): valid path for a .CSV file.
            
            [data, labels] = obj.createTable;             
            data.Properties.VariableNames(7:end) = labels;
            writetable(data,filename);
            msgbox(['Data saved to file : ' filename], 'to CSV');
        end
    end
   
    methods(Access = private)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function setTimeStamps(obj)
            % SETTIMESTAMPS sets the resolution for grouping data in
            % time based on the time_resolution property. New datetime values will 
            % be shifted to the start of the period 
            % (e.g., setTimeStamps("13-Jun-2021", 'month') = 01-Jun-2021);
            % Obs: "none" defaults to a time resolution of seconds.
            
            % Get Acquisitions' timestamps:
            tmstmp_list = NaT(length(obj.list_of_objs),1);
            for i=1:numel(obj.list_of_objs)
                elem = obj.list_of_objs{i};
                if isa(elem, 'Modality')
                    elem = elem.MyParent;
                end
                tmstmp_list(i) = elem.Start_datetime;
            end
            % Set resolution of tmpstmp_list:
            if strcmp(obj.time_resolution, 'none')
                obj.timestamp_list = dateshift(tmstmp_list, 'start', 'second');
            else
                obj.timestamp_list = dateshift(tmstmp_list, 'start', obj.time_resolution);
            end
        end
        %%% Validator Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function validateStatsInputs(obj)
            % VALIDATESTATSINPUTS checks if all objects in obj.list_of_objs
            % have obj.stats_filename in their respective SaveFolders.
            % Also, it checks if the sizes of list_of_objs and list_of_groups are equal.
            
            % If True, it creates an array of MAT-file objects and
            % saves to obj.MatFileList. 
            % Else, it displays a warning and removes the object with the 
            % missing obj.stats_filename from obj.list_of_objs.
            
            % Checks if list_of_objs have the same length as list_of_groups
            if size(obj.list_of_objs) ~= size(obj.list_of_groups)
                errID = 'Umitoolbox:StatsManager:IncompatibleInputSizes';
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
                errID = 'Umitoolbox:StatsManager:MissingData';
                errMsg = ['None of the objects listed contain "' obj.stats_filename '"'];
                error(errID, errMsg);
            else
                % Update lists:
                obj.MfileArr = matHandleArr(~b_remObj);
                obj.list_of_objs = obj.list_of_objs(~b_remObj);
                obj.list_of_groups = obj.list_of_groups(~b_remObj);
            end
        end
               
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function createDataArray(obj)
            % CREATEDATAArray generates a cell array containing all
            % data and metadata provided as inputs to STATSMANAGER class.
            
            disp('Creating Data array. Please wait...')
            
            obj.stats_data = cell(numel(obj.list_of_objs),11);
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
                obj.stats_data{i,9} = obj.stats_data{i,6}.obsID; % observation ID
                obj.stats_data{i,10}= obj.stats_data{i,6}.data; % observation data
            end
            
            % Create Relative time per subject's acquisitions:
            subjs = unique(obj.stats_data(:,2));
            for i = 1:numel(subjs)
                indx = find(strcmp(subjs{i}, obj.stats_data(:,2)));
                acq_time_list = datetime(vertcat(obj.stats_data{indx,5}));
                [~,tm_idx] = sort(acq_time_list);
                [~,rel_time]= sort(tm_idx);
               obj.stats_data(indx,11) = arrayfun(@(x) x, rel_time, 'UniformOutput', false);
            end
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
            %       "stats_data".
            % Outputs:
            %   out (table) : table containing the recording
            %       information and numerical data of the observation
            %       "obs_ID" for all recordings.
                        
            % Filter out recordings without the observation:
            idx = arrayfun(@(x) any(strcmp(obs_ID, {x.observations.ID})), obj.stats_data);            
            statsStruct = obj.stats_data(idx);
            
            indx_obs = arrayfun(@(x) find(strcmp(obs_ID, {x.observations.ID})), statsStruct);            
                        
            % Retrieve recordings' info and store in a cell array:            
            info = arrayfun(@(x,y) {x.subjID, x.acqID, x.modID,...
            x.groupName, datetime(x.acqTimeStamp), obs_ID},...
            statsStruct,indx_obs, 'UniformOutput', false);           
            info = info'; 
            info = vertcat(info{:});
            
            % Get observation's data:
            
            % Preallocate data with NaNs based on the length of unique
            % labels:
            data = nan(size(info,1),length(labelList));
            for i = 1:length(statsStruct)                
                [~,locB] = ismember(statsStruct(i).observations(indx_obs(i)).label,labelList);
                data(i,locB) = statsStruct(i).observations(indx_obs(i)).data;
            end
            % Prepare data to be added to "out" table:
            data = num2cell(data);
            
            % build table:
            out = table('Size', [size(info,1) size(info,2) + size(data,2)],...
                'VariableTypes',...
                [{'cellstr', 'cellstr', 'cellstr', 'cellstr', 'datetime', 'cellstr'}...
                 repmat({'single'},1,size(data,2))]);
            out(:,1:6) = info;
            out(:,7:end) = data;
            % Set table header:
            generic_label = arrayfun(@(x) ['Val_' num2str(x)], 1:size(data,2), 'UniformOutput',false);
            out.Properties.VariableNames = [{'Subject', 'Acquisition', 'Modality', 'GroupName',...
                    'RecordingTime', 'ObsID'}, generic_label];            
        end
        %%%%%%%%%%%%%%%%%%
        
    end
            
            
        
        
        
        
        
   
    
    
    
    
    
    
    
    
end
