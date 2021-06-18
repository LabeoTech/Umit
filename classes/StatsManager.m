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
        stats_data     % Structure containing all data and metadata created.
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
            obj.createDataStruct;
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
            assert(isa(stats_filename, 'char') & endsWith(stats_filename, '.mat'), errID, errMsg);
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
            % VALIDATESTATSINPUTSchecks if all objects in obj.list_of_objs
            % have obj.stats_filename in their respective FilePtr.json
            % files and if the sizes of list_of_objs and list_of_groups are equal.
            
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
            for i = 1:length(obj.list_of_objs)
                elem = obj.list_of_objs{i};
                filePtr = elem.FilePtr;
                txt = fileread(filePtr);
                a = jsondecode(txt);
                if isempty(a.Files)
                    b_remObj(i) = true;
                else
                    filenames = {a.Files.Name};
                    idx_file = strcmp(filenames, obj.stats_filename);
                    if sum(idx_file)== 0
                        b_remObj(i) = true;
                    else
                        filePath = fullfile(a.Files(idx_file).Folder, a.Files(idx_file).Name);
                        filePath = tokenizePath(filePath, elem, 'detokenize');
                        h = matfile(filePath);
                        matHandleArr{i} = h;
                    end
                end
            end
            
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
                
        function createDataStruct(obj)
            % CREATEDATASTRUCT generates a structure containing all
            % data and metadata provided as inputs to STATSMANAGER class.
            
            stats_data = struct();
            groupNames = unique(obj.list_of_groups);
            for i = 1:numel(groupNames)
                idx_group = strcmp(obj.list_of_groups, groupNames{i});
                stats_data(i).groupName = groupNames{i};
                stats_data(i).data = getElementInfo(obj,idx_group);
                
                
            end
            
            function out = getElementInfo(obj,idx)
                % getElementInfo gathers all data and metadata of all
                % elements from a group and saves it in a structure.
                obj_list = obj.list_of_objs(idx);
                matfile_list = obj.MfileArr(idx);
                
                
                
                
                
                
            end
                
            
            
        end
            
            
        
        
        
        
        
    end
    
    
    
    
    
    
    
    
end