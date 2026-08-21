classdef DataViewerPreferences
%DATAVIEWERPREFERENCES Persist and validate user-level DataViewer settings.
%
%   DataViewer settings are stored under the dataViewer field in the shared
%   umIT application preferences file. Theme remains the shared root theme
%   used by all umIT apps. Preferences never write dataset, Project, Rig,
%   or SaveFolder metadata.

    methods (Static)
        function preferences = defaults()
            %DEFAULTS Return safe DataViewer preference defaults.

            preferences = struct( ...
                'theme', 'light', ...
                'defaultColormap', 'parula', ...
                'reopenLastFile', false, ...
                'rememberLastSaveFolder', false, ...
                'defaultFolder', '', ...
                'lastFile', '', ...
                'lastSaveFolder', '');
        end

        function preferences = load(varargin)
            %LOAD Read and normalize persisted DataViewer preferences.
            %
            %   preferences = DataViewerPreferences.load()
            %   preferences = DataViewerPreferences.load( ...
            %       'AvailableColormaps', keys, ...
            %       'PreferencesFolder', folder)

            options = DataViewerPreferences.parseOptions(varargin{:});
            rootPreferences = readUmitPreferences( ...
                options.PreferencesFolder, true);

            saved = struct();
            if isfield(rootPreferences, 'dataViewer')
                saved = rootPreferences.dataViewer;
            end
            if isfield(rootPreferences, 'theme')
                saved.theme = rootPreferences.theme;
            end

            preferences = DataViewerPreferences.normalize( ...
                saved, options.AvailableColormaps);
        end

        function [preferencesFile, preferences] = save(preferences, varargin)
            %SAVE Validate and atomically persist DataViewer preferences.

            options = DataViewerPreferences.parseOptions(varargin{:});
            preferences = DataViewerPreferences.normalize( ...
                preferences, options.AvailableColormaps);

            rootPreferences = readUmitPreferences( ...
                options.PreferencesFolder, true);
            rootPreferences.schemaVersion = 1;
            rootPreferences.theme = preferences.theme;
            rootPreferences.dataViewer = rmfield(preferences, 'theme');
            preferencesFile = writeUmitPreferences( ...
                rootPreferences, options.PreferencesFolder);
        end

        function [preferencesFile, preferences] = reset(varargin)
            %RESET Restore DataViewer defaults while preserving other settings.

            options = DataViewerPreferences.parseOptions(varargin{:});
            [preferencesFile, preferences] = DataViewerPreferences.save( ...
                DataViewerPreferences.defaults(), ...
                'AvailableColormaps', options.AvailableColormaps, ...
                'PreferencesFolder', options.PreferencesFolder);
        end

        function [preferencesFile, preferences] = recordSuccessfulOpen( ...
                filePath, varargin)
            %RECORDSUCCESSFULOPEN Remember a successfully opened data file.

            filePath = DataViewerPreferences.normalizeText(filePath, '');
            if ~DataViewerPreferences.isSupportedDataFile(filePath) || ...
                    ~isfile(filePath)
                error('Umitoolbox:DataViewerPreferences:InvalidLastFile', ...
                    ['The last-opened file must be an existing supported ' ...
                    'DataViewer file.']);
            end

            options = DataViewerPreferences.parseOptions(varargin{:});
            preferences = DataViewerPreferences.load( ...
                'AvailableColormaps', options.AvailableColormaps, ...
                'PreferencesFolder', options.PreferencesFolder);
            if preferences.reopenLastFile
                preferences.lastFile = filePath;
            end

            saveFolder = fileparts(filePath);
            if preferences.rememberLastSaveFolder && isfolder(saveFolder)
                preferences.lastSaveFolder = saveFolder;
            end

            [preferencesFile, preferences] = DataViewerPreferences.save( ...
                preferences, ...
                'AvailableColormaps', options.AvailableColormaps, ...
                'PreferencesFolder', options.PreferencesFolder);
        end

        function filePath = resolveStartupFile(explicitFilePath, preferences)
            %RESOLVESTARTUPFILE Choose an explicit or remembered startup file.

            explicitFilePath = DataViewerPreferences.normalizeText( ...
                explicitFilePath, '');
            if ~isempty(explicitFilePath)
                filePath = explicitFilePath;
                return
            end

            preferences = DataViewerPreferences.normalize(preferences, {});
            rememberedFile = preferences.lastFile;
            if preferences.reopenLastFile && ...
                    DataViewerPreferences.isSupportedDataFile(rememberedFile) && ...
                    isfile(rememberedFile)
                filePath = rememberedFile;
            else
                filePath = '';
            end
        end

        function folder = resolveDialogFolder(fallbackFolder, preferences)
            %RESOLVEDIALOGFOLDER Apply remembered/default/fallback precedence.

            preferences = DataViewerPreferences.normalize(preferences, {});

            if preferences.rememberLastSaveFolder && ...
                    isfolder(preferences.lastSaveFolder)
                folder = preferences.lastSaveFolder;
                return
            end

            if isfolder(preferences.defaultFolder)
                folder = preferences.defaultFolder;
                return
            end

            fallbackFolder = DataViewerPreferences.normalizeText( ...
                fallbackFolder, '');
            if isfolder(fallbackFolder)
                folder = fallbackFolder;
            else
                folder = pwd;
            end
        end

        function tf = isSupportedDataFile(filePath)
            %ISSUPPORTEDDATAFILE True for DataViewer .dat and .umt files.

            filePath = DataViewerPreferences.normalizeText(filePath, '');
            [~, ~, extension] = fileparts(filePath);
            tf = ismember(lower(extension), {'.dat', '.umt'});
        end
    end

    methods (Static, Access = private)
        function preferences = normalize(value, availableColormaps)
            preferences = DataViewerPreferences.defaults();
            if ~(isstruct(value) && isscalar(value))
                value = struct();
            end

            preferences.theme = ...
                DataViewerPreferences.normalizeFieldTheme( ...
                value, 'theme', preferences.theme);
            preferences.defaultColormap = ...
                DataViewerPreferences.normalizeFieldText( ...
                value, 'defaultColormap', preferences.defaultColormap);
            preferences.reopenLastFile = ...
                DataViewerPreferences.normalizeFieldLogical( ...
                value, 'reopenLastFile', preferences.reopenLastFile);
            preferences.rememberLastSaveFolder = ...
                DataViewerPreferences.normalizeFieldLogical( ...
                value, 'rememberLastSaveFolder', ...
                preferences.rememberLastSaveFolder);
            preferences.defaultFolder = ...
                DataViewerPreferences.normalizeFieldText( ...
                value, 'defaultFolder', preferences.defaultFolder);
            preferences.lastFile = ...
                DataViewerPreferences.normalizeFieldText( ...
                value, 'lastFile', preferences.lastFile);
            preferences.lastSaveFolder = ...
                DataViewerPreferences.normalizeFieldText( ...
                value, 'lastSaveFolder', preferences.lastSaveFolder);

            availableColormaps = DataViewerPreferences.normalizeColormapList( ...
                availableColormaps);
            if isempty(availableColormaps)
                if ~isvarname(preferences.defaultColormap)
                    preferences.defaultColormap = 'parula';
                end
                return
            end

            if ismember('parula', availableColormaps)
                fallbackColormap = 'parula';
            else
                fallbackColormap = availableColormaps{1};
            end

            if ~ismember(preferences.defaultColormap, availableColormaps)
                preferences.defaultColormap = fallbackColormap;
            end
        end

        function value = normalizeFieldText(source, fieldName, fallback)
            if isfield(source, fieldName)
                value = DataViewerPreferences.normalizeText( ...
                    source.(fieldName), fallback);
            else
                value = fallback;
            end
        end

        function value = normalizeFieldLogical(source, fieldName, fallback)
            if isfield(source, fieldName) && ...
                    islogical(source.(fieldName)) && ...
                    isscalar(source.(fieldName))
                value = source.(fieldName);
            else
                value = fallback;
            end
        end

        function value = normalizeFieldTheme(source, fieldName, fallback)
            value = DataViewerPreferences.normalizeFieldText( ...
                source, fieldName, fallback);
            value = lower(value);
            if ~ismember(value, {'light', 'dark'})
                value = fallback;
            end
        end

        function value = normalizeText(value, fallback)
            if ischar(value) && (isrow(value) || isempty(value))
                value = strtrim(value);
            elseif isstring(value) && isscalar(value) && ~ismissing(value)
                value = strtrim(char(value));
            else
                value = fallback;
            end
        end

        function colormaps = normalizeColormapList(value)
            if isempty(value)
                colormaps = {};
                return
            end

            if ischar(value)
                value = {value};
            elseif isstring(value)
                value = cellstr(value(:));
            end

            if ~iscell(value)
                colormaps = {};
                return
            end

            colormaps = cell(0, 1);
            for idx = 1:numel(value)
                item = DataViewerPreferences.normalizeText(value{idx}, '');
                if ~isempty(item) && isvarname(item)
                    colormaps{end + 1, 1} = item; %#ok<AGROW>
                end
            end
            colormaps = unique(colormaps, 'stable');
        end

        function options = parseOptions(varargin)
            parser = inputParser;
            parser.FunctionName = 'DataViewerPreferences';
            addParameter(parser, 'AvailableColormaps', {}, ...
                @(value) ischar(value) || isstring(value) || iscell(value));
            addParameter(parser, 'PreferencesFolder', '', ...
                @(value) ischar(value) || ...
                (isstring(value) && isscalar(value)));
            parse(parser, varargin{:});

            options = parser.Results;
            options.PreferencesFolder = ...
                DataViewerPreferences.normalizeText( ...
                options.PreferencesFolder, '');
        end
    end
end
