classdef DataViewer_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                  matlab.ui.Figure
        FileMenu                  matlab.ui.container.Menu
        OpenMenu                  matlab.ui.container.Menu
        PreviewRawMenu            matlab.ui.container.Menu
        SaveMenu                  matlab.ui.container.Menu
        ImportMenu                matlab.ui.container.Menu
        frombinMenu               matlab.ui.container.Menu
        fromtifMenu               matlab.ui.container.Menu
        ExportMenu                matlab.ui.container.Menu
        totifMenu                 matlab.ui.container.Menu
        EditmetadataMenu          matlab.ui.container.Menu
        SetRawFolderMenu          matlab.ui.container.Menu
        HelpMenu_2                matlab.ui.container.Menu
        PreferencesMenu           matlab.ui.container.Menu
        OnlinedocumentationMenu   matlab.ui.container.Menu
        GridLayoutMain            matlab.ui.container.GridLayout
        GridLayoutSub1            matlab.ui.container.GridLayout
        GridLayoutRight           matlab.ui.container.GridLayout
        EventSwitchGrid           matlab.ui.container.GridLayout
        Switch                    matlab.ui.control.Switch
        EventsLabel               matlab.ui.control.Label
        EventsGrid                matlab.ui.container.GridLayout
        RestoreButton             matlab.ui.control.Button
        DeleteRepetitionButton    matlab.ui.control.Button
        DeleteConditionButton     matlab.ui.control.Button
        RepetitionDropDown        matlab.ui.control.DropDown
        RepetitionDropDownLabel   matlab.ui.control.Label
        ConditionDropDown         matlab.ui.control.DropDown
        ConditionDropDownLabel    matlab.ui.control.Label
        UITable                   matlab.ui.control.Table
        ROItableLabel             matlab.ui.control.Label
        TemporalProfileLabel      matlab.ui.control.Label
        PlotAxes                  matlab.ui.control.UIAxes
        GridLayoutLeft            matlab.ui.container.GridLayout
        DataFolderContextLabel    matlab.ui.control.Label
        StatusLabel               matlab.ui.control.Label
        ImageViewLabel            matlab.ui.control.Label
        TopGrid                   matlab.ui.container.GridLayout
        NextFrameButton           matlab.ui.control.Button
        Slider                    matlab.ui.control.Slider
        MovieSpeedLabel           matlab.ui.control.Label
        MovieSpeedDropDown        matlab.ui.control.DropDown
        PreviousFrameButton       matlab.ui.control.Button
        PlayMovieButton           matlab.ui.control.Button
        BottomGrid                matlab.ui.container.GridLayout
        ImageStatusLabel          matlab.ui.control.Label
        SetClipButton             matlab.ui.control.Button
        ClipSliderLabel           matlab.ui.control.Label
        HidecrosshairCheckBox     matlab.ui.control.CheckBox
        Panel                     matlab.ui.container.Panel
        ClipSliderRange           matlab.ui.control.RangeSlider
        AutoButton                matlab.ui.control.Button
        InvertCheckBox            matlab.ui.control.CheckBox
        ColormapDropDown          matlab.ui.control.DropDown
        ColormapDropDownLabel     matlab.ui.control.Label
        ImageAxes                 matlab.ui.control.UIAxes
        TabGroup                  matlab.ui.container.TabGroup
        UtilitiesTab              matlab.ui.container.Tab
        GridUtilitiesTab          matlab.ui.container.GridLayout
        AlignmentToolLabel        matlab.ui.control.Label
        ImageAlignButton          matlab.ui.control.Button
        OiS2CamAlignmentLabel     matlab.ui.control.Label
        OiSDUalCamCoregButton     matlab.ui.control.Button
        EventsManagerLabel        matlab.ui.control.Label
        DataHistoryLabel          matlab.ui.control.Label
        ImageReferenceLabel       matlab.ui.control.Label
        LogicalMaskLabel          matlab.ui.control.Label
        ViewCalibrationLabel      matlab.ui.control.Label
        EventsManagerButton       matlab.ui.control.Button
        DataHistoryButton         matlab.ui.control.Button
        ImageRefButton            matlab.ui.control.Button
        LogicalMaskButton         matlab.ui.control.Button
        ViewCalibButton           matlab.ui.control.Button
        PipelineTab               matlab.ui.container.Tab
        GridPipelineTab           matlab.ui.container.GridLayout
        LatestsummaryLabel        matlab.ui.control.Label
        CurrentpipelinefileLabel  matlab.ui.control.Label
        LastrunoutputsLabel       matlab.ui.control.Label
        PipelineStatusLabel       matlab.ui.control.Label
        PipelineManagerLabel      matlab.ui.control.Label
        PipeLauncherButton        matlab.ui.control.Button
        ROIsTab                   matlab.ui.container.Tab
        GridROIsTab               matlab.ui.container.GridLayout
        Sep                       matlab.ui.container.Panel
        AllenAtlasLabel           matlab.ui.control.Label
        ThresholdLabel            matlab.ui.control.Label
        ExportROIdataLabel        matlab.ui.control.Label
        ImportLabel               matlab.ui.control.Label
        ROIAllenButton            matlab.ui.control.Button
        ROIbyTresholdButton       matlab.ui.control.Button
        ExportROIDataButton       matlab.ui.control.Button
        ImportROIButton           matlab.ui.control.Button
        LoadLabel                 matlab.ui.control.Label
        LoadROIButton             matlab.ui.control.Button
        SaveLabel                 matlab.ui.control.Label
        DeleteLabel               matlab.ui.control.Label
        PolygonLabel              matlab.ui.control.Label
        EllipeLabel               matlab.ui.control.Label
        RectangleLabel            matlab.ui.control.Label
        SaveROIButton             matlab.ui.control.Button
        DeleteROIButton           matlab.ui.control.Button
        DrawPolygonButton         matlab.ui.control.Button
        DrawEllipseButton         matlab.ui.control.Button
        DrawRectangleButton       matlab.ui.control.Button
    end


    properties (Access = private)
        % =====================================================================
        % Data backend and file identity
        % =====================================================================
        DataSource = []
        CurrentFile char = ''
        CurrentEntry char = ''
        DataRawFolder char = ''
        DataRawFolderSaveFolder char = ''
        RawFolderInvalidWarningSaveFolder char = ''

        % =====================================================================
        % PipelineManagerTool integration
        % =====================================================================
        PipelineManagerObj = []
        PipelineManagerToolApp = []
        PipelineDataViewerFileSourceNodeID double = NaN
        LastPipelineResult = []
        PipelineRunHistory = cell(0,1)
        PipelineTemporaryFiles table = table( ...
            strings(0,1), ...
            datetime.empty(0,1), ...
            strings(0,1), ...
            'VariableNames', {'FilePath','CreatedOn','Source'})

        % =====================================================================
        % Interaction and viewer state
        % =====================================================================
        InteractionMode char = 'idle'
        CurrentFrame double = 1
        CrosshairXY double = [1 1]   % [X Y]

        % =====================================================================
        % Event/view state
        % =====================================================================
        ViewMode char = 'normal'     % 'normal' or 'event'
        EventSource = []             % EventsManager for .dat; [] for .umt
        EventInfo struct = struct()  % Normalized event info for GUI use
        CurrentCondition char = ''
        CurrentRepetition char = 'AVERAGE'
        EventFrameMatrix double = []
        EventConditionIDList = []
        EventRepetitionList = []
        EventConditionColors double = zeros(0, 3)
        bWarnedMissingUMTBaseline logical = false

        % =====================================================================
        % Display/cache state and graphics
        % =====================================================================
        ShowCacheRectangle logical = true
        OriginalClipSliderLimits double = [0 1]
        ColormapLibrary struct = struct()

        ImageHandle = matlab.graphics.primitive.Image.empty
        CrosshairHandles = gobjects(0)
        CacheRectHandle = gobjects(0)
        CrossTraceHandle = gobjects(0)
        CrossTraceSEMHandle = gobjects(0)
        TimeBarHandle = gobjects(0)
        EventPatchHandles = gobjects(0)
        PlotAxesYLimListener = []

        % =====================================================================
        % Image context menu
        % =====================================================================
        ImageContextMenu matlab.ui.container.ContextMenu
        UpdateCacheMenu matlab.ui.container.Menu
        ToggleCacheLockMenu matlab.ui.container.Menu
        ToggleCacheAreaMenu matlab.ui.container.Menu

        % =====================================================================
        % Spatial metadata and origin graphics
        % =====================================================================
        DataParams struct = struct()
        OriginCrosshairHandles = gobjects(0)

        % =====================================================================
        % Movie playback
        % =====================================================================
        MovieTimer timer = timer.empty
        MovieStartFrame double = 1
        MovieStartTic = []
        PlayIconFile char = 'icon_rightTriangle_Black.png'
        StopIconFile char = 'icon_stopButton_Black.png'

        % =====================================================================
        % ROI core model and selection
        % =====================================================================
        ROIList struct = struct([])
        SelectedROIID double = NaN
        ROIColorList double = distinguishable_colors(50, {'w','b'})
        ROINextColorIndex double = 1
        ROISelectionOrder double = []      % order-dependent operations, e.g. subtraction

        LastImageClickTic = []
        LastImageClickXY double = [NaN NaN]

        % =====================================================================
        % ROI edit runtime
        % =====================================================================
        GroupROIEditState struct = struct()
        GroupEditBoxHandle = gobjects(1)
        GroupEditPreviewHandles = gobjects(0)

        % =====================================================================
        % ROI delete context menu
        % =====================================================================
        DeleteROIContextMenu matlab.ui.container.ContextMenu
        DeleteSelectedROIsMenu matlab.ui.container.Menu
        DeleteAllROIsMenu matlab.ui.container.Menu

        % =====================================================================
        % Group ROI edit / Boolean Ops context menus
        % =====================================================================
        GroupROIEditContextMenu matlab.ui.container.ContextMenu
        GroupROIDeleteMenu matlab.ui.container.Menu
        GroupROICancelMenu matlab.ui.container.Menu

        GroupROIBooleanContextMenu matlab.ui.container.ContextMenu
        GroupROIBooleanOpsMenu matlab.ui.container.Menu
        GroupROIIntersectMenu matlab.ui.container.Menu
        GroupROIMergeMenu matlab.ui.container.Menu
        GroupROIXORMenu matlab.ui.container.Menu
        GroupROISubtractMenu matlab.ui.container.Menu

        % =====================================================================
        % Allen Atlas ROI placement runtime
        % =====================================================================
        AllenAtlasROICreatorApp = []
        AtlasROIPlacementState struct = struct()
        AtlasROIPlacementBoxHandle = gobjects(1)
        AtlasROIPlacementRegionHandles = gobjects(0)
        AtlasROIPlacementBoundaryHandles = gobjects(0)
        AtlasROIPlacementBregmaHandle = gobjects(1)
        AtlasROIPlacementLambdaHandle = gobjects(1)
        AtlasROIPlacementContextMenu matlab.ui.container.ContextMenu
        AtlasROIPlacementConfirmMenu matlab.ui.container.Menu
        AtlasROIPlacementCancelMenu matlab.ui.container.Menu

        % =====================================================================
        % Logical mask runtime state
        % =====================================================================
        LogicalMaskOverlayHandle = matlab.graphics.primitive.Image.empty
        LogicalMaskDraftHandles = gobjects(0)
        LogicalMaskDraftListeners cell = {}
        LogicalMaskPreviousImageAxesContextMenu = []
        LogicalMaskPreviousImageHandleContextMenu = []
        LogicalMaskPreviousKeyFcn = []

        LogicalMaskContextMenu matlab.ui.container.ContextMenu
        LogicalMaskAddPolygonMenu matlab.ui.container.Menu
        LogicalMaskConfirmMenu matlab.ui.container.Menu
        LogicalMaskCancelMenu matlab.ui.container.Menu

        % =====================================================================
        % Logical mask button context menu
        % =====================================================================
        LogicalMaskButtonContextMenu matlab.ui.container.ContextMenu
        CreateManualLogicalMaskMenu matlab.ui.container.Menu
        CreateAutomaticLogicalMaskMenu matlab.ui.container.Menu
        ResetLogicalMaskMenu matlab.ui.container.Menu

        % =====================================================================
        % Image reference button context menu
        % =====================================================================
        SetReferenceButtonContextMenu matlab.ui.container.ContextMenu
        PreviewImageReferenceMenu matlab.ui.container.Menu


    end

    methods (Access = private)
        % =====================================================================
        % Runtime setup and dialogs
        % =====================================================================
        function initializeViewerGraphics(app)
            %INITIALIZEVIEWERGRAPHICS Create persistent graphics objects and static styles.
            %
            %   This method runs once during startup. It creates placeholder graphics
            %   objects and sets non-data-dependent visual properties. Data-dependent
            %   properties such as image size, axes limits, CLim, tick labels, and cache
            %   rectangle position are configured later by setAxes.

            cla(app.ImageAxes);
            cla(app.PlotAxes);

            % -------------------------------------------------------------------------
            % Image placeholder
            % -------------------------------------------------------------------------
            app.ImageHandle = imagesc(app.ImageAxes, zeros(100, 100, 'single'));
            app.ImageHandle.ButtonDownFcn = @(src, evt) app.onImageClicked(src, evt);
            app.ImageAxes.ButtonDownFcn = @(src, evt) app.onImageClicked(src, evt);
            app.ImageHandle.PickableParts = 'all';
            app.ImageHandle.HitTest = 'on';

            axis(app.ImageAxes, 'image');
            app.ImageAxes.Box = 'on';
            app.ImageAxes.YDir = 'reverse';
            app.ImageAxes.TickDir = 'out';
            app.ImageAxes.Title.Interpreter = 'none';

            % Cache menu is created now but attached only when a partial cache is used.
            app.ImageAxes.ContextMenu = [];
            app.ImageHandle.ContextMenu = [];

            % -------------------------------------------------------------------------
            % Spatial origin crosshair. Hidden unless DataParams origin differs from [1 1].
            % -------------------------------------------------------------------------
            app.OriginCrosshairHandles = gobjects(2, 1);

            app.OriginCrosshairHandles(1) = line(app.ImageAxes, [nan nan], [nan nan], ...
                'Color', [0 0 0], ...
                'LineWidth', 0.5, ...
                'HitTest', 'off', ...
                'PickableParts', 'none', ...
                'HandleVisibility', 'off', ...
                'Visible', 'off');

            app.OriginCrosshairHandles(2) = line(app.ImageAxes, [nan nan], [nan nan], ...
                'Color', [0 0 0], ...
                'LineWidth', 0.5, ...
                'HitTest', 'off', ...
                'PickableParts', 'none', ...
                'HandleVisibility', 'off', ...
                'Visible', 'off');

            % -------------------------------------------------------------------------
            % Active white crosshair
            % -------------------------------------------------------------------------
            app.CrosshairHandles = gobjects(2, 1);

            app.CrosshairHandles(1) = line(app.ImageAxes, [1 1], [1 100], ...
                'Color', [1 1 1], ...
                'LineWidth', 1, ...
                'HitTest', 'off', ...
                'PickableParts', 'none', ...
                'HandleVisibility', 'off');

            app.CrosshairHandles(2) = line(app.ImageAxes, [1 100], [1 1], ...
                'Color', [1 1 1], ...
                'LineWidth', 1, ...
                'HitTest', 'off', ...
                'PickableParts', 'none', ...
                'HandleVisibility', 'off');

            % -------------------------------------------------------------------------
            % Cache rectangle
            % -------------------------------------------------------------------------
            app.CacheRectHandle = rectangle(app.ImageAxes, ...
                'Position', [1 1 1 1], ...
                'EdgeColor', [1 1 1], ...
                'LineStyle', '--', ...
                'LineWidth', 1, ...
                'Visible', 'off', ...
                'HitTest', 'off', ...
                'PickableParts', 'none', ...
                'HandleVisibility', 'off');

            % -------------------------------------------------------------------------
            % Temporal plot handles
            % -------------------------------------------------------------------------
            hold(app.PlotAxes, 'on');

            % Event patches are created on demand by refreshEventPatches.
            app.EventPatchHandles = gobjects(0);

            app.CrossTraceSEMHandle = patch(app.PlotAxes, ...
                'XData', nan, ...
                'YData', nan, ...
                'FaceColor', [0 0 0], ...
                'FaceAlpha', 0.18, ...
                'EdgeColor', 'none', ...
                'Visible', 'off', ...
                'HitTest', 'off', ...
                'PickableParts', 'none', ...
                'HandleVisibility', 'off');

            app.CrossTraceHandle = plot(app.PlotAxes, nan, nan, ...
                'LineWidth', 1, ...
                'HandleVisibility', 'off', ...
                'Color', [.5 .5 .5]);

            app.TimeBarHandle = line(app.PlotAxes, [nan nan], [nan nan], ...
                'Color', [0 0 0], ...
                'LineStyle', '-', ...
                'LineWidth', 1.5, ...
                'Visible', 'off', ...
                'HitTest', 'off', ...
                'PickableParts', 'none', ...
                'HandleVisibility', 'off');

            hold(app.PlotAxes, 'off');

            app.PlotAxes.Box = 'on';
            app.PlotAxes.TickDir = 'out';

            xlabel(app.PlotAxes, 'Time (s)');
            ylabel(app.PlotAxes, 'Amplitude');

            % Keep helper graphics aligned with PlotAxes vertical zoom/pan/restore.
            try
                app.PlotAxesYLimListener = addlistener( ...
                    app.PlotAxes, ...
                    'YLim', ...
                    'PostSet', ...
                    @(~, ~) app.onPlotAxesYLimChanged());
            catch
                app.PlotAxesYLimListener = [];
            end

            app.createImageContextMenu();
            app.setStatusMessage('Ready.');

        end

        function createImageContextMenu(app)
            %CREATEIMAGECONTEXTMENU Create right-click menu for image/cache actions.
            %
            %   The menu is created once during startup, but it is only
            %   attached to the image axis when the active backend exposes a
            %   partial temporal cache.

            app.ImageContextMenu = uicontextmenu(app.UIFigure);

            app.UpdateCacheMenu = uimenu(app.ImageContextMenu);
            app.UpdateCacheMenu.Text = 'Update cache around crosshair';
            app.UpdateCacheMenu.MenuSelectedFcn = @(src, evt) app.onUpdateCacheAroundCrosshair();

            app.ToggleCacheLockMenu = uimenu(app.ImageContextMenu);
            app.ToggleCacheLockMenu.Text = 'Lock cache';
            app.ToggleCacheLockMenu.MenuSelectedFcn = @(src, evt) app.onToggleCacheLock();

            app.ToggleCacheAreaMenu = uimenu(app.ImageContextMenu);
            app.ToggleCacheAreaMenu.Text = 'Hide cache area';
            app.ToggleCacheAreaMenu.MenuSelectedFcn = @(src, evt) app.onToggleCacheRectangle();

            app.ImageAxes.ContextMenu = [];
            if ~isempty(app.ImageHandle) && isvalid(app.ImageHandle)
                app.ImageHandle.ContextMenu = [];
            end

        end

        function selectedEntry = selectUMTEntryDialog(app, entrySummary)
            %SELECTUMTENTRYDIALOG Modal dialog for choosing a UMT data entry.

            selectedEntry = '';

            if isempty(entrySummary) || height(entrySummary) == 0
                return
            end

            selectedRow = 1;

            dlg = uifigure( ...
                'Name', 'Select UMT entry', ...
                'WindowStyle', 'modal', ...
                'Position', [100 100 560 320], ...
                'Visible','off',...
                'CloseRequestFcn', @onCancel);

            grid = uigridlayout(dlg);
            grid.RowHeight = {25, '1x', 25, 35};
            grid.ColumnWidth = {'1x'};
            grid.Padding = [12 12 12 12];

            titleLabel = uilabel(grid);
            titleLabel.Text = 'Select the UMT data entry to display:';
            titleLabel.FontWeight = 'bold';
            titleLabel.Layout.Row = 1;
            titleLabel.Layout.Column = 1;

            entryTable = uitable(grid);
            entryTable.Data = entrySummary;
            entryTable.ColumnEditable = false(1, width(entrySummary));
            entryTable.RowName = {};
            entryTable.Layout.Row = 2;
            entryTable.Layout.Column = 1;
            entryTable.SelectionChangedFcn = @onSelectionChanged;

            % Prefer whole-row selection when available.
            if isprop(entryTable, 'SelectionType')
                entryTable.SelectionType = 'row';
            end

            if isprop(entryTable, 'Multiselect')
                entryTable.Multiselect = 'off';
            end

            % Initialize visible selection on first row.
            try
                if isprop(entryTable, 'SelectionType') && strcmpi(entryTable.SelectionType, 'row')
                    entryTable.Selection = 1;
                else
                    entryTable.Selection = [1 1];
                end
            catch
            end

            statusLabel = uilabel(grid);
            statusLabel.Text = '';
            statusLabel.FontColor = [0.65 0 0];
            statusLabel.Layout.Row = 3;
            statusLabel.Layout.Column = 1;

            buttonGrid = uigridlayout(grid);
            buttonGrid.RowHeight = {'1x'};
            buttonGrid.ColumnWidth = {'1x', '1x'};
            buttonGrid.Padding = [0 0 0 0];
            buttonGrid.Layout.Row = 4;
            buttonGrid.Layout.Column = 1;

            okButton = uibutton(buttonGrid, 'push');
            okButton.Text = 'OK';
            okButton.Layout.Row = 1;
            okButton.Layout.Column = 1;
            okButton.ButtonPushedFcn = @onOK;

            cancelButton = uibutton(buttonGrid, 'push');
            cancelButton.Text = 'Cancel';
            cancelButton.Layout.Row = 1;
            cancelButton.Layout.Column = 2;
            cancelButton.ButtonPushedFcn = @onCancel;

            if exist('placeAppInsideCaller', 'file') == 2
                placeAppInsideCaller(app, dlg, 'center');
            end
            dlg.Visible = 'on';

            uiwait(dlg);

            if isvalid(dlg)
                delete(dlg);
            end

            function onSelectionChanged(~, event)
                selectedRow = iGetSelectedRow(event);

                if isempty(selectedRow) || selectedRow < 1 || selectedRow > height(entrySummary)
                    selectedRow = [];
                    return
                end

                % Force whole-row highlight when possible.
                try
                    if isprop(entryTable, 'SelectionType') && strcmpi(entryTable.SelectionType, 'row')
                        entryTable.Selection = selectedRow;
                    else
                        entryTable.Selection = [selectedRow, 1];
                    end
                catch
                end

                statusLabel.Text = '';
            end

            function rowIdx = iGetSelectedRow(event)
                rowIdx = [];

                if isprop(event, 'Selection') && ~isempty(event.Selection)
                    selection = event.Selection;
                elseif isprop(entryTable, 'Selection') && ~isempty(entryTable.Selection)
                    selection = entryTable.Selection;
                else
                    return
                end

                if isvector(selection)
                    rowIdx = selection(1);
                else
                    rowIdx = selection(1, 1);
                end

                rowIdx = round(double(rowIdx));
            end

            function onOK(~, ~)
                if isempty(selectedRow) || selectedRow < 1 || selectedRow > height(entrySummary)
                    statusLabel.Text = 'Select one entry.';
                    return
                end

                selectedEntry = char(entrySummary.Name(selectedRow));
                uiresume(dlg);
            end

            function onCancel(~, ~)
                selectedEntry = '';
                uiresume(dlg);
            end

        end

        % =====================================================================
        % Data loading and metadata
        % =====================================================================

        function openDataFromDialog(app, startFolder)
            %OPENDATAFROMDIALOG Open .dat or .umt image data.
            %
            %   openDataFromDialog(app) opens the file picker from the current data
            %   folder when available, otherwise from pwd.
            %
            %   openDataFromDialog(app, startFolder) opens the file picker from the
            %   supplied folder. This is useful after data import, where the imported
            %   files are expected in the selected save folder.

            if nargin < 2 || isempty(startFolder)
                if isempty(app.CurrentFile)
                    startFolder = pwd;
                else
                    startFolder = fileparts(app.CurrentFile);
                end
            end

            startFolder = char(string(startFolder));

            if ~isfolder(startFolder)
                if isempty(app.CurrentFile)
                    startFolder = pwd;
                else
                    startFolder = fileparts(app.CurrentFile);
                end
            end

            [file, folder] = uigetfile( ...
                {'*.dat;*.umt', 'UMIT image data (*.dat, *.umt)'; ...
                '*.dat', 'Raw DAT image data (*.dat)'; ...
                '*.umt', 'UMT image data (*.umt)'}, ...
                'Open image data', ...
                startFolder);

            if isequal(file, 0)
                return
            end

            filePath = fullfile(folder, file);

            % Opening another file is one of the lifecycle boundaries where
            % temporary PipelineManager outputs should be cleaned. Preserve the file
            % that is about to be opened, because it may itself be a temporary output.
            canContinue = app.cleanupPipelineTemporaryFiles( ...
                'PreserveFile', filePath, ...
                'AskForCurrentTemporaryFile', true);

            if ~canContinue
                app.setStatusMessage('Open file cancelled during temporary-file cleanup.');
                return
            end

            bLoaded = app.loadDataSource(filePath);

            if bLoaded
                app.refreshViewerAfterLoad();
            end

            % Focus back on the app.
            figure(app.UIFigure);
        end

        function bLoaded = loadDataSource(app, filePath)
            %LOADDATASOURCE Create the correct backend for a selected image file.

            bLoaded = false;

            if ~isfile(filePath)
                error('DataViewer:fileNotFound', ...
                    'File not found: "%s".', filePath);
            end

            previousSource = app.DataSource;
            previousFile = app.CurrentFile;
            previousEntry = app.CurrentEntry;
            previousEventSource = app.EventSource;
            previousEventInfo = app.EventInfo;
            previousViewMode = app.ViewMode;

            app.setInteractionMode('loading');

            try
                [~, ~, ext] = fileparts(filePath);
                ext = lower(ext);

                app.setStatusMessage(sprintf('Opening "%s"...', filePath));

                switch ext
                    case '.dat'
                        newSource = DatImageSource(filePath);
                        newEntry = '';

                    case '.umt'
                        entrySummary = UMTImageSource.inspectEntries(filePath);

                        if height(entrySummary) == 1
                            selectedEntry = char(entrySummary.Name(1));
                        else
                            selectedEntry = app.selectUMTEntryDialog(entrySummary);

                            if isempty(selectedEntry)
                                app.setStatusMessage('UMT loading cancelled.');
                                app.setInteractionMode('idle');
                                return
                            end
                        end

                        newSource = UMTImageSource(filePath, selectedEntry);
                        newEntry = selectedEntry;

                    otherwise
                        error('DataViewer:unsupportedFileType', ...
                            'Unsupported file type: "%s".', ext);
                end

                app.resetPipelineContextForNewSaveFolder(previousFile, filePath);

                app.DataSource = newSource;
                app.CurrentFile = filePath;
                app.CurrentEntry = newEntry;
                app.CurrentFrame = 1;
                app.bWarnedMissingUMTBaseline = false;

                app.loadDataParamsForCurrentFile();
                app.initializeEventsForLoadedData(ext);

                sz = app.getDataSize();
                app.CrosshairXY = [round(sz(2) / 2), round(sz(1) / 2)];

                if isempty(newEntry)
                    app.setStatusMessage(sprintf('Loaded "%s".', filePath));
                else
                    app.setStatusMessage(sprintf('Loaded "%s" | Entry: %s.', ...
                        filePath, newEntry));
                end

                bLoaded = true;
                app.setInteractionMode('idle');
                app.updateDataFolderContextLabel();
                app.updatePipelineTabState();

            catch ME
                app.DataSource = previousSource;
                app.CurrentFile = previousFile;
                app.CurrentEntry = previousEntry;
                app.EventSource = previousEventSource;
                app.EventInfo = previousEventInfo;
                app.ViewMode = previousViewMode;

                app.setInteractionMode('idle');
                app.updateDataFolderContextLabel();
                app.updatePipelineTabState();
                rethrow(ME)
            end
        end

        function loadDataParamsForCurrentFile(app)
            %LOADDATAPARAMSFORCURRENTFILE Load folder-global DataParams.mat.

            app.DataParams = struct();

            if isempty(app.CurrentFile)
                app.ensureDataParamsFolderFields();
                app.updateImageStatusLabel();
                return
            end

            folderPath = fileparts(app.CurrentFile);
            dataParamsPath = fullfile(folderPath, 'DataParams.mat');

            if ~isfile(dataParamsPath)
                app.setStatusMessage('DataParams.mat not found. Using pixel coordinates and full logical mask.');
                app.ensureDataParamsViewFields();
                app.ensureDataParamsMaskFields();
                app.ensureDataParamsFolderFields();
                app.updateImageStatusLabel();
                return
            end

            try
                % loadDataParams normalizes older DataParams schemas before validation.
                app.DataParams = loadDataParams(folderPath);
                app.ensureDataParamsViewFields();
                app.ensureDataParamsMaskFields();
                app.ensureDataParamsFolderFields();

            catch ME
                app.DataParams = struct();
                app.ensureDataParamsViewFields();
                app.ensureDataParamsMaskFields();
                app.ensureDataParamsFolderFields();

                warnStatus = warning();
                warning('off','backtrace');
                warning('DataViewer:invalidDataParams', ...
                    'DataParams.mat could not be loaded or validated.\n%s', ME.message);
                warning(warnStatus);
                app.setStatusMessage('Invalid DataParams.mat. Using pixel coordinates and full logical mask.');
            end

            app.updateImageStatusLabel();

        end

        function Info = getSourceInfo(app)
            %GETSOURCEINFO Return metadata from the active backend.

            Info = struct();

            if ~app.hasData()
                return
            end

            if ismethod(app.DataSource, 'getInfo')
                Info = app.DataSource.getInfo();
            end

        end

        function sz = getSourceDataSize(app)
            %GETSOURCEDATASIZE Return normalized backend size [Y X T E].

            rawSize = app.DataSource.getSize();

            if numel(rawSize) < 3
                rawSize(3) = 1;
            end

            if numel(rawSize) < 4
                rawSize(4) = 1;
            end

            sz = double(rawSize(1:4));

        end

        function sourceType = getSourceType(app)
            %GETSOURCETYPE Return current backend source type: 'dat', 'umt', or ''.

            sourceType = '';

            if isempty(app.CurrentFile)
                return
            end

            [~, ~, ext] = fileparts(app.CurrentFile);
            switch lower(ext)
                case '.dat'
                    sourceType = 'dat';
                case '.umt'
                    sourceType = 'umt';
            end

        end

        function tf = hasData(app)
            %HASDATA Return true when a data source is loaded.

            tf = ~isempty(app.DataSource) && isvalid(app.DataSource);

        end

        function sz = getDataSize(app)
            %GETDATASIZE Return active display size [Y X T E].

            rawSize = app.getSourceDataSize();
            sz = rawSize;

            if strcmpi(app.ViewMode, 'event') && strcmpi(app.getSourceType(), 'dat') && ...
                    ~isempty(app.EventFrameMatrix)
                sz(3) = size(app.EventFrameMatrix, 2);
                sz(4) = max(1, size(app.EventFrameMatrix, 1));
            end

        end

        function tf = hasValidFrameRate(app)
            %HASVALIDFRAMERATE Return true when metadata has a usable frame rate.

            frameRateHz = app.getFrameRateHz();
            tf = ~isempty(frameRateHz) && isfinite(frameRateHz) && frameRateHz > 0;

        end

        function frameRateHz = getFrameRateHz(app)
            %GETFRAMERATEHZ Return frame rate from active backend metadata.

            frameRateHz = [];

            Info = app.getSourceInfo();

            if isempty(Info) || ~isstruct(Info)
                return
            end

            if isfield(Info, 'FrameRateHz') && ~isempty(Info.FrameRateHz)
                frameRateHz = double(Info.FrameRateHz);
            elseif isfield(Info, 'Freq') && ~isempty(Info.Freq)
                frameRateHz = double(Info.Freq);
            end

            if ~isscalar(frameRateHz) || ~isfinite(frameRateHz) || frameRateHz <= 0
                frameRateHz = [];
            end

        end

        % =====================================================================
        % Event metadata and selection state
        % =====================================================================

        function initializeEventsForLoadedData(app, fileExt)
            %INITIALIZEEVENTSFORLOADEDDATA Load/normalize event metadata for current data.

            app.EventSource = [];
            app.EventInfo = struct();
            app.CurrentCondition = '';
            app.CurrentRepetition = 'AVERAGE';
            app.EventFrameMatrix = [];
            app.EventConditionIDList = [];
            app.EventRepetitionList = [];
            app.EventConditionColors = zeros(0, 3);
            app.deleteEventPatches();

            if ~app.hasData()
                app.ViewMode = 'normal';
                app.populateEventControls();
                return
            end

            switch lower(fileExt)
                case '.dat'
                    app.ViewMode = 'normal';
                    saveFolder = fileparts(app.CurrentFile);
                    eventsPath = fullfile(saveFolder, 'events.mat');

                    if isfile(eventsPath)
                        try
                            app.EventSource = EventsManager(saveFolder);
                            app.EventInfo = app.normalizeDatEventInfo();
                        catch ME
                            app.EventSource = [];
                            app.EventInfo = struct();
                            warning('DataViewer:EventsLoadFailed', ...
                                'Failed to load events.mat.\n%s', ME.message);
                            app.setStatusMessage('events.mat could not be loaded. Events disabled.');
                        end
                    else
                        app.setStatusMessage('No events.mat found. Use Event Manager utility to create one.');
                    end

                case '.umt'
                    app.EventInfo = app.normalizeUMTEventInfo();

                    sourceSize = app.getSourceDataSize();
                    if numel(sourceSize) >= 4 && sourceSize(4) > 1
                        app.ViewMode = 'event';
                    else
                        app.ViewMode = 'normal';
                    end
            end

            if app.hasNormalizedEvents()
                nCond = numel(app.EventInfo.EventNames);
                if nCond > 0
                    app.EventConditionColors = lines(nCond);
                    app.CurrentCondition = app.EventInfo.EventNames{1};
                end
            end

            app.syncSwitchToViewMode();
            app.populateEventControls();
            app.rebuildEventFrameMatrix();

        end

        function tf = hasNormalizedEvents(app)
            %HASNORMALIZEDEVENTS True when normalized EventInfo is usable.

            tf = isstruct(app.EventInfo) && isscalar(app.EventInfo) && ...
                isfield(app.EventInfo, 'HasEvents') && app.EventInfo.HasEvents;

        end

        function ev = normalizeDatEventInfo(app)
            %NORMALIZEDATEVENTINFO Normalize EventsManager metadata for GUI use.

            ev = struct();
            ev.SourceType = 'dat';
            ev.EventAxisMode = 'instances';
            ev.HasEvents = false;
            ev.EventIDs = [];
            ev.RepetitionIndex = [];
            ev.EventNames = {};
            ev.EventNamePerEvent = {};
            ev.HasRepetitions = false;
            ev.BaselinePeriod = [];

            if isempty(app.EventSource) || isempty(app.EventSource.eventID)
                return
            end

            EM = app.EventSource;

            if isempty(EM.selectedEvents) || ~islogical(EM.selectedEvents) || ...
                    ~isequal(size(EM.selectedEvents), size(EM.eventID))
                EM.selectedEvents = true(size(EM.eventID));
            end

            onMask = EM.state & EM.selectedEvents;

            if ~any(onMask)
                return
            end

            eventIDs = double(EM.eventID(onMask));
            repIdx = double(EM.repetitionID(onMask));

            validIDs = unique(eventIDs, 'stable');
            eventNames = EM.eventNameList(validIDs);
            eventNamePerEvent = EM.eventNameList(eventIDs);

            ev.SourceType = 'dat';
            ev.EventAxisMode = 'instances';
            ev.HasEvents = true;
            ev.EventIDs = eventIDs(:);
            ev.RepetitionIndex = repIdx(:);
            ev.EventNames = eventNames(:).';
            ev.EventNamePerEvent = eventNamePerEvent(:);
            ev.HasRepetitions = any(repIdx > 0);
            ev.BaselinePeriod = double(EM.baselinePeriod);

        end

        function ev = normalizeUMTEventInfo(app)
            %NORMALIZEUMTEVENTINFO Normalize embedded UMT eventInfo for GUI use.

            ev = struct();
            ev.SourceType = 'umt';
            ev.EventAxisMode = '';
            ev.HasEvents = false;
            ev.EventIDs = [];
            ev.RepetitionIndex = [];
            ev.EventNames = {};
            ev.EventNamePerEvent = {};
            ev.HasRepetitions = false;
            ev.BaselinePeriod = [];

            if ~app.hasData() || ~ismethod(app.DataSource, 'getEventInfo')
                return
            end

            info = app.DataSource.getEventInfo();

            if ~isstruct(info) || ~isscalar(info) || isempty(fieldnames(info))
                return
            end

            requiredFields = {'eventID','repetitionIndex','eventName','eventAxisMode'};
            if ~all(isfield(info, requiredFields))
                return
            end

            eventIDs = double(info.eventID(:));
            repIdx = double(info.repetitionIndex(:));
            eventNamesPerEvent = cellstr(string(info.eventName(:)));
            axisMode = lower(char(string(info.eventAxisMode)));

            uniqueIDs = unique(eventIDs, 'stable');
            uniqueNames = cell(1, numel(uniqueIDs));
            for iID = 1:numel(uniqueIDs)
                idx = find(eventIDs == uniqueIDs(iID), 1, 'first');
                uniqueNames{iID} = eventNamesPerEvent{idx};
            end

            ev.SourceType = 'umt';
            ev.EventAxisMode = axisMode;
            ev.HasEvents = true;
            ev.EventIDs = eventIDs(:);
            ev.RepetitionIndex = repIdx(:);
            ev.EventNames = uniqueNames;
            ev.EventNamePerEvent = eventNamesPerEvent(:);
            ev.HasRepetitions = strcmp(axisMode, 'instances') && any(repIdx > 0);

            ev.BaselinePeriod = app.getUMTBaselinePeriod();

        end

        function baselinePeriod = getUMTBaselinePeriod(app)
            %GETUMTBASELINEPERIOD Return baselinePeriod from UMT if available.

            baselinePeriod = [];

            if ~app.hasData()
                return
            end

            % Future-compatible probes. Current UMT schema may not yet allow this field.
            try
                if isprop(app.DataSource, 'UMT') && isfield(app.DataSource.UMT, 'baselinePeriod')
                    baselinePeriod = double(app.DataSource.UMT.baselinePeriod);
                    return
                end
            catch
            end

            try
                eventInfo = app.DataSource.getEventInfo();
                if isstruct(eventInfo) && isfield(eventInfo, 'baselinePeriod')
                    baselinePeriod = double(eventInfo.baselinePeriod);
                    return
                end
            catch
            end

            Info = app.getSourceInfo();
            if isstruct(Info) && isfield(Info, 'baselinePeriod') && ~isempty(Info.baselinePeriod)
                baselinePeriod = double(Info.baselinePeriod);
            end

        end

        function baselinePeriod = getActiveBaselinePeriod(app)
            %GETACTIVEBASELINEPERIOD Return baseline period for current event source.

            baselinePeriod = [];

            if strcmpi(app.getSourceType(), 'dat') && ~isempty(app.EventSource)
                baselinePeriod = double(app.EventSource.baselinePeriod);
                return
            end

            if app.hasNormalizedEvents() && isfield(app.EventInfo, 'BaselinePeriod') && ...
                    ~isempty(app.EventInfo.BaselinePeriod)
                baselinePeriod = double(app.EventInfo.BaselinePeriod);
                return
            end

            if strcmpi(app.getSourceType(), 'umt') && strcmpi(app.ViewMode, 'event') && ...
                    ~app.bWarnedMissingUMTBaseline
                app.bWarnedMissingUMTBaseline = true;

                try
                    uialert(app.UIFigure, ...
                        ['baselinePeriod was not found in the UMT metadata. ' ...
                        'Time will be shown from frame 1 instead of event-aligned.'], ...
                        'Missing baselinePeriod', ...
                        'Icon', 'warning');
                catch
                    warning('DataViewer:MissingUMTBaselinePeriod', ...
                        'baselinePeriod was not found in the UMT metadata.');
                end
            end

        end

        function tf = hasDatEventSource(app)
            %HASDATEVENTSOURCE True when a .dat EventsManager object is loaded.

            tf = strcmpi(app.getSourceType(), 'dat') && ...
                ~isempty(app.EventSource) && ...
                isvalid(app.EventSource) && ...
                isprop(app.EventSource, 'eventID') && ...
                ~isempty(app.EventSource.eventID);

        end

        function tf = hasIgnoredDatEvents(app)
            %HASIGNOREDDATEVENTS True when any .dat event is ignored in EventsManager.

            tf = false;

            if ~app.hasDatEventSource()
                return
            end

            if ~isprop(app.EventSource, 'selectedEvents') || ...
                    isempty(app.EventSource.selectedEvents)
                return
            end

            tf = any(~app.EventSource.selectedEvents(:));

        end

        function syncSwitchToViewMode(app)
            %SYNCSWITCHTOVIEWMODE Update Normal/Event switch value from app.ViewMode.

            if isempty(app.Switch) || ~isvalid(app.Switch)
                return
            end

            switch lower(app.ViewMode)
                case 'event'
                    app.Switch.Value = 'Event mode';
                otherwise
                    app.Switch.Value = 'Normal mode';
            end

        end

        function populateEventControls(app)
            %POPULATEEVENTCONTROLS Populate condition/repetition dropdowns.

            if isempty(app.ConditionDropDown) || ~isvalid(app.ConditionDropDown) || ...
                    isempty(app.RepetitionDropDown) || ~isvalid(app.RepetitionDropDown)
                return
            end

            if ~app.hasNormalizedEvents()
                app.ConditionDropDown.Items = {'No events'};
                app.ConditionDropDown.ItemsData = {''};
                app.ConditionDropDown.Value = '';

                app.RepetitionDropDown.Items = {'AVERAGE'};
                app.RepetitionDropDown.ItemsData = {'AVERAGE'};
                app.RepetitionDropDown.Value = 'AVERAGE';
                return
            end

            conditionNames = app.EventInfo.EventNames;
            if isempty(conditionNames)
                conditionNames = {'No events'};
            end

            app.ConditionDropDown.Items = conditionNames;
            app.ConditionDropDown.ItemsData = conditionNames;

            if isempty(app.CurrentCondition) || ~ismember(app.CurrentCondition, conditionNames)
                app.CurrentCondition = conditionNames{1};
            end
            app.ConditionDropDown.Value = app.CurrentCondition;

            app.populateRepetitionDropDownForCurrentCondition();

        end

        function populateRepetitionDropDownForCurrentCondition(app)
            %POPULATEREPETITIONDROPDOWNFORCURRENTCONDITION Populate repetitions.

            if ~app.hasNormalizedEvents() || isempty(app.CurrentCondition)
                app.RepetitionDropDown.Items = {'AVERAGE'};
                app.RepetitionDropDown.ItemsData = {'AVERAGE'};
                app.RepetitionDropDown.Value = 'AVERAGE';
                app.CurrentRepetition = 'AVERAGE';
                return
            end

            condID = app.getConditionID(app.CurrentCondition);

            switch app.EventInfo.SourceType
                case 'umt'
                    if strcmpi(app.EventInfo.EventAxisMode, 'aggregated_repetitions')
                        app.RepetitionDropDown.Items = {'AVERAGE'};
                        app.RepetitionDropDown.ItemsData = {'AVERAGE'};
                        app.RepetitionDropDown.Value = 'AVERAGE';
                        app.CurrentRepetition = 'AVERAGE';
                        return
                    end
            end

            idxCond = app.EventInfo.EventIDs == condID;
            repList = unique(app.EventInfo.RepetitionIndex(idxCond), 'stable');
            repList = repList(repList > 0);

            items = [{'AVERAGE'}, cellstr(string(repList(:).'))];
            itemsData = items;

            app.RepetitionDropDown.Items = items;
            app.RepetitionDropDown.ItemsData = itemsData;

            if isempty(app.CurrentRepetition) || ~ismember(app.CurrentRepetition, itemsData)
                app.CurrentRepetition = 'AVERAGE';
            end

            app.RepetitionDropDown.Value = app.CurrentRepetition;

        end

        function condID = getConditionID(app, conditionName)
            %GETCONDITIONID Return numeric condition ID for one condition name.

            condID = [];

            if ~app.hasNormalizedEvents() || isempty(conditionName)
                return
            end

            idx = find(strcmp(app.EventInfo.EventNames, conditionName), 1, 'first');
            if isempty(idx)
                return
            end

            % EventNames are built from stable unique EventIDs.
            uniqueIDs = unique(app.EventInfo.EventIDs, 'stable');
            condID = uniqueIDs(idx);

        end

        function rebuildEventFrameMatrix(app)
            %REBUILDEVENTFRAMEMATRIX Build .dat event frame matrix for current selection.

            app.EventFrameMatrix = [];
            app.EventConditionIDList = [];
            app.EventRepetitionList = [];

            if ~strcmpi(app.getSourceType(), 'dat') || ~strcmpi(app.ViewMode, 'event') || ...
                    isempty(app.EventSource) || isempty(app.CurrentCondition)
                return
            end

            sourceSize = app.getSourceDataSize();
            datLen = sourceSize(3);

            if strcmpi(app.CurrentRepetition, 'AVERAGE')
                [frMat, conditionIDlist, repetitionList] = ...
                    app.EventSource.getFrameMatrix(datLen, app.CurrentCondition);

                firstInvalidCol = find(any(isnan(frMat), 1), 1, 'first');
                if ~isempty(firstInvalidCol)
                    frMat = frMat(:, 1:firstInvalidCol-1);
                end
            else
                repIdx = str2double(app.CurrentRepetition);
                [frMat, conditionIDlist, repetitionList] = ...
                    app.EventSource.getFrameMatrix(datLen, app.CurrentCondition, repIdx);
            end

            app.EventFrameMatrix = frMat;
            app.EventConditionIDList = conditionIDlist;
            app.EventRepetitionList = repetitionList;

        end

        function eIdxList = getSelectedUMTEventIndices(app)
            %GETSELECTEDUMTEVENTINDICES Return UMT E indices for current selection.

            eIdxList = [];

            if ~app.hasNormalizedEvents() || ~strcmpi(app.EventInfo.SourceType, 'umt')
                return
            end

            condID = app.getConditionID(app.CurrentCondition);
            if isempty(condID)
                return
            end

            idx = find(app.EventInfo.EventIDs == condID);

            if strcmpi(app.EventInfo.EventAxisMode, 'aggregated_repetitions')
                eIdxList = idx(:).';
                return
            end

            if strcmpi(app.CurrentRepetition, 'AVERAGE')
                eIdxList = idx(:).';
            else
                repIdx = str2double(app.CurrentRepetition);
                idx = idx(app.EventInfo.RepetitionIndex(idx) == repIdx);
                eIdxList = idx(:).';
            end

        end

        function refreshDatEventsAfterEdit(app)
            %REFRESHDATEVENTSAFTEREDIT Refresh GUI after EventsManager selection edits.
            %
            %   This method is used after removeCondition, removeRepetition, or
            %   clearIgnoredEvents. It rebuilds the normalized event metadata from the
            %   EventsManager state and refreshes all event-dependent display elements.
            %
            %   If all events are ignored, the viewer automatically returns to normal
            %   mode and keeps Restore available so the user can recover event display.

            if ~app.hasDatEventSource()
                return
            end

            app.EventInfo = app.normalizeDatEventInfo();

            bNoEventsAvailable = ~app.hasNormalizedEvents();
            bHasIgnoredEvents = app.hasIgnoredDatEvents();

            if bNoEventsAvailable
                app.ViewMode = 'normal';
                app.CurrentCondition = '';
                app.CurrentRepetition = 'AVERAGE';
                app.EventFrameMatrix = [];
                app.EventConditionIDList = [];
                app.EventRepetitionList = [];

                app.syncSwitchToViewMode();
                app.populateEventControls();

                app.CurrentFrame = 1;

                % Force refresh from the normal YXT source.
                app.refreshFrameControls();
                app.refreshImageFrame();
                title(app.ImageAxes, app.getImageTitle());
                app.refreshTemporalProfile();
                app.refreshEventPatches();

                % Deleting the last condition changes the ROI trace frame basis from
                % event-aligned frames back to normal source frames.
                app.refreshROITraces();
                app.updateROIStatsForCurrentFrame();

                if bHasIgnoredEvents
                    app.setStatusMessage('All events are ignored. Click Restore to show events again.');

                    try
                        uialert(app.UIFigure, ...
                            ['All events are currently ignored, so the viewer was switched ' ...
                            'back to normal mode. Click Restore in the Events panel to ' ...
                            'show the events again.'], ...
                            'All events ignored', ...
                            'Icon', 'warning');
                    catch
                        warning('DataViewer:AllEventsIgnored', ...
                            ['All events are currently ignored. Click Restore in the Events ' ...
                            'panel to show the events again.']);
                    end
                end

                app.updateGUIEnabledState();
                return
            end

            conditionNames = app.EventInfo.EventNames;

            if isempty(app.CurrentCondition) || ~ismember(app.CurrentCondition, conditionNames)
                app.CurrentCondition = conditionNames{1};
            end

            app.CurrentRepetition = 'AVERAGE';

            app.syncSwitchToViewMode();
            app.populateEventControls();

            app.CurrentFrame = 1;
            app.rebuildEventFrameMatrix();
            app.refreshViewerForModeChange();

            % Explicit final state update. refreshViewerForModeChange already calls this,
            % but keeping it here makes event-edit button refresh unambiguous.
            app.updateGUIEnabledState();

        end

        % =====================================================================
        % Axes setup and spatial calibration
        % =====================================================================

        function refreshViewerAfterLoad(app)
            %REFRESHVIEWERAFTERLOAD Initialize all viewer panels after loading data.

            if ~app.hasData()
                return
            end

            app.setAxes();
            app.refreshCrosshairGraphics();
            app.refreshFrameControls();
            app.refreshTemporalProfile();
            app.refreshCacheOverlay();
            app.refreshCacheMenuState();
            app.refreshEventPatches();
            app.updateImageStatusLabel();

            % Synchronize SaveFolder-level RawFolder context from DataParams.mat after
            % the viewer is ready. This may show a one-time warning if the stored path
            % is invalid, for example after an external-drive letter changed.
            app.synchronizeRawFolderFromDataParams();

            % ROIs intentionally persist across opened data. If existing ROI masks match
            % the new image size, their traces and current-frame spatial statistics are
            % recomputed for the new source.
            app.refreshROITraces();
            app.updateROIStatsForCurrentFrame();

            app.updateDataFolderContextLabel();
            app.updatePipelineTabState();
            app.updateGUIEnabledState();
        end


        function setAxes(app)
            %SETAXES Configure data-dependent axes properties after loading a dataset.
            %
            %   This method runs after a new data source is loaded. It does not recreate
            %   graphics objects. It updates loaded-data-dependent properties such as
            %   CData, axes limits, CLim, spatial ticks, origin marker, cache overlay,
            %   temporal-axis labeling, and cache context-menu attachment.

            if ~app.hasData()
                return
            end

            frame = app.getCurrentFrame();
            sz = app.getDataSize();

            % -------------------------------------------------------------------------
            % Image content and image axes
            % -------------------------------------------------------------------------
            set(app.ImageHandle, 'CData', frame);

            app.ImageAxes.XLim = [0.5, sz(2) + 0.5];
            app.ImageAxes.YLim = [0.5, sz(1) + 0.5];

            title(app.ImageAxes, app.getImageTitle());

            % Initial CLim for the loaded dataset/frame. Later frame updates must not
            % reset CLim.
            finiteFrame = frame(isfinite(frame));

            if isempty(finiteFrame)
                clim = [0 1];
            else
                clim = double([min(finiteFrame(:)), max(finiteFrame(:))]);

                if clim(1) == clim(2)
                    clim = clim + [-1, 1];
                end
            end

            % Data loading is one of the allowed moments to reset slider limits.
            app.setDisplayCLim(clim, true);

            % Store data-loading limits as persistent restore target.
            app.OriginalClipSliderLimits = clim;

            % Image tick labels and origin marker from DataParams.
            app.setImageAxisSpatialCalibration();

            % -------------------------------------------------------------------------
            % Crosshair and overlays
            % -------------------------------------------------------------------------
            app.refreshLogicalMaskOverlay();
            app.refreshCrosshairGraphics();
            app.refreshCacheOverlay();

            % -------------------------------------------------------------------------
            % Temporal plot axis
            % -------------------------------------------------------------------------
            app.PlotAxes.Box = 'on';
            app.PlotAxes.TickDir = 'out';
            app.PlotAxes.YLimMode = 'manual';
            ylabel(app.PlotAxes, 'Amplitude');

            if app.hasValidFrameRate()
                xlabel(app.PlotAxes, 'Time (s)');
            else
                xlabel(app.PlotAxes, 'Frame');
            end

            % Trace data are updated separately so cache misses/status can be handled by
            % refreshTemporalProfile.
            app.hideTemporalTrace();

            % -------------------------------------------------------------------------
            % Cache context menu
            % -------------------------------------------------------------------------
            app.configureCacheContextMenu();

        end

        function refreshFrameControls(app)
            %REFRESHFRAMECONTROLS Synchronize frame slider values safely.
            %
            %   This method updates only frame-control values/properties that depend on
            %   the active display size. It does not enable or disable GUI elements.

            if ~app.hasData()
                return
            end

            if isempty(app.Slider) || ~isvalid(app.Slider)
                return
            end

            sz = app.getDataSize();
            nFrames = max(1, round(sz(3)));

            if nFrames <= 1
                app.Slider.Limits = [1 2];
                app.Slider.Value = 1;
                app.Slider.MajorTicks = 1;
            else
                app.Slider.Limits = [1 nFrames];
                app.Slider.Value = min(max(app.CurrentFrame, 1), nFrames);
                app.Slider.MajorTicks = unique(round(linspace(1, nFrames, min(5, nFrames))));
            end

        end

        function setImageAxisSpatialCalibration(app)
            %SETIMAGEAXISSPATIALCALIBRATION Set image tick labels from DataParams.
            %
            %   Image axes remain internally in pixel coordinates. Tick labels are shown
            %   relative to DataParams.view.origin_xy_px. If pixelSize_px_per_mm is
            %   empty, labels are shown as pixel offsets. If positive, labels are shown
            %   in mm.

            if ~app.hasData()
                return
            end

            sz = app.getDataSize();
            Ny = sz(1);
            Nx = sz(2);

            % Defaults.
            originXY = [1 1];
            pxPerMm = [];

            if ~isempty(app.DataParams) && isfield(app.DataParams, 'view')
                viewInfo = app.DataParams.view;

                if isfield(viewInfo, 'origin_xy_px') && ~isempty(viewInfo.origin_xy_px)
                    candidateOrigin = double(viewInfo.origin_xy_px(:).');

                    if numel(candidateOrigin) == 2 && all(isfinite(candidateOrigin)) && ...
                            candidateOrigin(1) >= 1 && candidateOrigin(1) <= Nx && ...
                            candidateOrigin(2) >= 1 && candidateOrigin(2) <= Ny
                        originXY = candidateOrigin;
                    else
                        app.setStatusMessage('Invalid DataParams origin. Using [1 1].');
                    end
                end

                if isfield(viewInfo, 'imageSizeYX') && ~isempty(viewInfo.imageSizeYX)
                    dataParamsSize = double(viewInfo.imageSizeYX(:).');
                    if numel(dataParamsSize) == 2 && ~isequal(dataParamsSize, [Ny Nx])
                        app.setStatusMessage('DataParams image size does not match loaded data. Using current image size.');
                    end
                end

                if isfield(viewInfo, 'pixelSize_px_per_mm') && ...
                        ~isempty(viewInfo.pixelSize_px_per_mm)

                    candidatePxPerMm = double(viewInfo.pixelSize_px_per_mm(:).');

                    if isscalar(candidatePxPerMm) && isfinite(candidatePxPerMm) && candidatePxPerMm > 0
                        pxPerMm = candidatePxPerMm;
                    else
                        app.setStatusMessage('Invalid pixel ratio. Showing pixel offsets.');
                    end
                end
            end

            xTicks = app.makePixelTicks(Nx);
            yTicks = app.makePixelTicks(Ny);

            if isempty(pxPerMm)
                % Pixel-offset mode: still center labels at origin.
                xValues = xTicks - originXY(1);
                yValues = yTicks - originXY(2);

                app.ImageAxes.XTick = xTicks;
                app.ImageAxes.YTick = yTicks;
                app.ImageAxes.XTickLabel = app.formatSpatialTickLabels(xValues);
                app.ImageAxes.YTickLabel = app.formatSpatialTickLabels(yValues);

                xlabel(app.ImageAxes, 'X (px)');
                ylabel(app.ImageAxes, 'Y (px)');
            else
                % Calibrated mm mode.
                xValues = (xTicks - originXY(1)) ./ pxPerMm;
                yValues = (yTicks - originXY(2)) ./ pxPerMm;

                app.ImageAxes.XTick = xTicks;
                app.ImageAxes.YTick = yTicks;
                app.ImageAxes.XTickLabel = app.formatSpatialTickLabels(xValues);
                app.ImageAxes.YTickLabel = app.formatSpatialTickLabels(yValues);

                xlabel(app.ImageAxes, 'X (mm)');
                ylabel(app.ImageAxes, 'Y (mm)');
            end

            app.updateOriginCrosshair(originXY);

        end

        function setImageAxisPixelTicks(app, Ny, Nx)
            %SETIMAGEAXISPIXELTICKS Set default pixel-based image ticks.

            xTicks = app.makePixelTicks(Nx);
            yTicks = app.makePixelTicks(Ny);

            app.ImageAxes.XTick = xTicks;
            app.ImageAxes.YTick = yTicks;
            app.ImageAxes.XTickLabel = app.formatSpatialTickLabels(xTicks);
            app.ImageAxes.YTickLabel = app.formatSpatialTickLabels(yTicks);

            xlabel(app.ImageAxes, 'X (px)');
            ylabel(app.ImageAxes, 'Y (px)');

        end

        function ticks = makePixelTicks(app, nPixels) %#ok<INUSD>
            %MAKEPIXELTICKS Create a small set of readable pixel tick positions.

            nPixels = max(1, round(nPixels));
            nTicks = min(6, nPixels);

            ticks = unique(round(linspace(1, nPixels, nTicks)));

        end

        function labels = formatSpatialTickLabels(app, values) %#ok<INUSD>
            %FORMATSPATIALTICKLABELS Convert numeric tick values to compact labels.

            values = double(values(:));
            labels = cell(numel(values), 1);

            for iVal = 1:numel(values)
                if abs(values(iVal)) >= 100
                    labels{iVal} = sprintf('%.0f', values(iVal));
                elseif abs(values(iVal)) >= 10
                    labels{iVal} = sprintf('%.1f', values(iVal));
                else
                    labels{iVal} = sprintf('%.2f', values(iVal));
                end
            end

        end

        function updateOriginCrosshair(app, originXY)
            %UPDATEORIGINCROSSHAIR Show or hide the spatial 0,0 origin marker.

            if isempty(app.OriginCrosshairHandles) || any(~isvalid(app.OriginCrosshairHandles))
                return
            end

            if isempty(originXY) || numel(originXY) ~= 2 || any(~isfinite(originXY))
                set(app.OriginCrosshairHandles, 'Visible', 'off');
                return
            end

            % Do not show the origin marker for the default top-left pixel origin.
            if isequal(round(originXY), [1 1])
                set(app.OriginCrosshairHandles, 'Visible', 'off');
                return
            end

            sz = app.getDataSize();
            Ny = sz(1);
            Nx = sz(2);

            x0 = originXY(1);
            y0 = originXY(2);

            set(app.OriginCrosshairHandles(1), ...
                'XData', [x0 x0], ...
                'YData', [0.5 Ny + 0.5], ...
                'Visible', 'on');

            set(app.OriginCrosshairHandles(2), ...
                'XData', [0.5 Nx + 0.5], ...
                'YData', [y0 y0], ...
                'Visible', 'on');

        end

        % =====================================================================
        % Image, temporal trace, and event patch refresh
        % =====================================================================

        function refreshImageFrame(app)
            %REFRESHIMAGEFRAME Update displayed image content only.
            %
            %   This method is used for frame navigation/movie playback. It must not
            %   recreate graphics handles or reset axes properties such as CLim,
            %   colormap, zoom, spatial calibration, or cache overlays.

            if ~app.hasData()
                return
            end

            frame = app.getCurrentFrame();
            set(app.ImageHandle, 'CData', frame);

        end

        function refreshCrosshairGraphics(app)
            %REFRESHCROSSHAIRGRAPHICS Update crosshair line positions.

            if isempty(app.CrosshairHandles) || any(~isvalid(app.CrosshairHandles))
                return
            end

            sz = app.getDataSize();
            x = app.CrosshairXY(1);
            y = app.CrosshairXY(2);

            set(app.CrosshairHandles(1), ...
                'XData', [x x], ...
                'YData', [1 sz(1)]);

            set(app.CrosshairHandles(2), ...
                'XData', [1 sz(2)], ...
                'YData', [y y]);

        end

        function refreshTemporalProfile(app)
            %REFRESHTEMPORALPROFILE Update the pixel temporal trace.
            %
            %   PlotAxes.YLim is intentionally not updated here. The temporal plot
            %   y-limits are bound to ImageAxes.CLim and are controlled only by
            %   setDisplayCLim.

            if ~app.hasData()
                return
            end

            usesCache = app.sourceUsesTemporalCache();

            if usesCache
                x = app.CrosshairXY(1);
                y = app.CrosshairXY(2);
                wasInsideCache = app.DataSource.isInsideCache(y, x);

                if ~wasInsideCache && ~app.isCacheLocked()
                    app.setStatusMessage('Updating temporal cache...');
                    drawnow limitrate
                end
            else
                wasInsideCache = true;
            end

            [trace, status, semTrace, xData] = app.getCurrentPixelTrace();

            switch status
                case 'outside_locked_cache'
                    app.hideTemporalTrace();
                    app.setStatusMessage('Selected pixel is outside the locked temporal cache.');
                    return

                case {'ok', 'cache_rebuilt'}
                    % Continue.

                otherwise
                    app.hideTemporalTrace();
                    app.setStatusMessage(sprintf('No temporal trace available. Status: %s', status));
                    return
            end

            if isempty(trace) || isscalar(trace)
                app.hideTemporalTrace();
                app.setStatusMessage('No temporal trace available for selected pixel.');
                return
            end

            set(app.CrossTraceHandle, ...
                'XData', xData, ...
                'YData', trace(:), ...
                'Visible', 'on');

            app.updateCrosshairSEMPatch(xData, trace(:), semTrace(:));

            % Respect the current crosshair visibility setting.
            if isvalid(app.HidecrosshairCheckBox) && app.HidecrosshairCheckBox.Value
                app.CrossTraceHandle.Visible = 'off';
                if ~isempty(app.CrossTraceSEMHandle) && isvalid(app.CrossTraceSEMHandle)
                    app.CrossTraceSEMHandle.Visible = 'off';
                end
            end

            if numel(xData) > 1
                app.PlotAxes.XLim = [xData(1) xData(end)];
            else
                app.PlotAxes.XLim = [0 1];
            end

            app.PlotAxes.YLimMode = 'manual';
            app.refreshTimeBar();
            app.refreshEventPatches();

            if usesCache && (~wasInsideCache || strcmp(status, 'cache_rebuilt'))
                app.refreshCacheOverlay();
                app.setStatusMessage('Temporal cache updated.');
            end

        end

        function updateCrosshairSEMPatch(app, xData, meanTrace, semTrace)
            %UPDATECROSSHAIRSEMPATCH Update SEM patch behind crosshair trace.

            if isempty(app.CrossTraceSEMHandle) || ~isvalid(app.CrossTraceSEMHandle)
                return
            end

            if isempty(semTrace) || numel(semTrace) ~= numel(meanTrace) || ...
                    ~any(isfinite(semTrace))
                set(app.CrossTraceSEMHandle, 'XData', nan, 'YData', nan, 'Visible', 'off');
                return
            end

            upperTrace = meanTrace(:).' + semTrace(:).';
            lowerTrace = meanTrace(:).' - semTrace(:).';

            patchX = [xData(:).', fliplr(xData(:).')];
            patchY = [upperTrace, fliplr(lowerTrace)];

            set(app.CrossTraceSEMHandle, ...
                'XData', patchX, ...
                'YData', patchY, ...
                'FaceColor', app.CrossTraceHandle.Color, ...
                'Visible', 'on');

            uistack(app.CrossTraceSEMHandle, 'bottom');
            uistack(app.CrossTraceHandle, 'top');

        end

        function hideTemporalTrace(app)
            %HIDETEMPORALTRACE Hide crosshair temporal trace, SEM patch, and frame marker.

            if ~isempty(app.CrossTraceHandle) && isvalid(app.CrossTraceHandle)
                app.CrossTraceHandle.XData = nan;
                app.CrossTraceHandle.YData = nan;
                app.CrossTraceHandle.Visible = 'off';
            end

            if ~isempty(app.CrossTraceSEMHandle) && isvalid(app.CrossTraceSEMHandle)
                app.CrossTraceSEMHandle.XData = nan;
                app.CrossTraceSEMHandle.YData = nan;
                app.CrossTraceSEMHandle.Visible = 'off';
            end

            if ~isempty(app.TimeBarHandle) && isvalid(app.TimeBarHandle)
                app.TimeBarHandle.Visible = 'off';
            end

        end

        function refreshTimeBar(app, frameIdx)
            %REFRESHTIMEBAR Update current-frame marker on temporal plot.

            if nargin < 2 || isempty(frameIdx)
                frameIdx = app.CurrentFrame;
            end

            if isempty(app.TimeBarHandle) || ~isvalid(app.TimeBarHandle)
                return
            end

            yl = app.PlotAxes.YLim;
            xNow = app.getCurrentFrameTime(frameIdx);

            set(app.TimeBarHandle, ...
                'XData', [xNow xNow], ...
                'YData', yl, ...
                'Visible', 'on');

            uistack(app.TimeBarHandle, 'top');

        end

        function xData = getTraceTimeVector(app, nSamples)
            %GETTRACETIMEVECTOR Return trace XData in seconds or frames.

            nSamples = max(1, round(nSamples));
            frameRateHz = app.getFrameRateHz();

            if isempty(frameRateHz)
                xData = 1:nSamples;
                return
            end

            if strcmpi(app.ViewMode, 'event')
                baselinePeriod = app.getActiveBaselinePeriod();

                if ~isempty(baselinePeriod)
                    eventOnsetFrame = round(baselinePeriod * frameRateHz) + 1;
                    xData = ((1:nSamples) - eventOnsetFrame) ./ frameRateHz;
                else
                    xData = (0:nSamples-1) ./ frameRateHz;
                end
            else
                xData = (0:nSamples-1) ./ frameRateHz;
            end

        end

        function xNow = getCurrentFrameTime(app, frameIdx)
            %GETCURRENTFRAMETIME Return frame location on plot X axis.

            if nargin < 2 || isempty(frameIdx)
                frameIdx = app.CurrentFrame;
            end

            xData = app.getTraceTimeVector(max(1, frameIdx));
            xNow = xData(min(max(round(frameIdx), 1), numel(xData)));

        end

        function onPlotAxesYLimChanged(app)
            %ONPLOTAXESYLIMCHANGED Keep PlotAxes helper graphics aligned to Y limits.

            if isempty(app.PlotAxes) || ~isvalid(app.PlotAxes)
                return
            end

            yLim = app.PlotAxes.YLim;

            if numel(yLim) ~= 2 || any(~isfinite(yLim)) || yLim(1) == yLim(2)
                return
            end

            app.updateEventPatchYExtents(yLim);

            if app.isUsableGraphicsHandle(app.TimeBarHandle) && isprop(app.TimeBarHandle, 'YData')
                try
                    xData = app.TimeBarHandle.XData;

                    % Only resize an already valid time bar. Do not revive a hidden or
                    % uninitialized marker during axis interactions.
                    if ~isempty(xData) && any(isfinite(xData(:)))
                        app.TimeBarHandle.YData = yLim;
                    end
                catch
                end
            end

            try
                app.stackEventPatchesBottom();

                if app.isUsableGraphicsHandle(app.TimeBarHandle)
                    uistack(app.TimeBarHandle, 'top');
                end
            catch
            end

        end

        function refreshEventPatches(app)
            %REFRESHEVENTPATCHES Draw event patches in the temporal plot.

            app.deleteEventPatches();

            if ~app.hasNormalizedEvents() || ~app.hasData()
                legend(app.PlotAxes, 'off');
                return
            end

            switch app.ViewMode
                case 'normal'
                    app.drawNormalModeEventPatches();
                case 'event'
                    app.drawEventModeEventPatch();
            end

            if ~isempty(app.EventPatchHandles) && any(isvalid(app.EventPatchHandles))
                try
                    legend(app.PlotAxes, 'show', 'Location', 'northeast', 'Interpreter', 'none');
                catch
                end
            else
                legend(app.PlotAxes, 'off');
            end

            % Ensure helper graphics remain above event patches.
            try
                if ~isempty(app.CrossTraceSEMHandle) && isvalid(app.CrossTraceSEMHandle)
                    uistack(app.CrossTraceSEMHandle, 'top');
                end
                if ~isempty(app.CrossTraceHandle) && isvalid(app.CrossTraceHandle)
                    uistack(app.CrossTraceHandle, 'top');
                end
                if ~isempty(app.TimeBarHandle) && isvalid(app.TimeBarHandle)
                    uistack(app.TimeBarHandle, 'top');
                end
            catch
            end

        end

        function drawNormalModeEventPatches(app)
            %DRAWNORMALMODEEVENTPATCHES Draw continuous-time event patches.

            if ~strcmpi(app.getSourceType(), 'dat') || isempty(app.EventSource)
                return
            end

            EM = app.EventSource;
            if isempty(EM.timestamps) || isempty(EM.eventID) || isempty(EM.state)
                return
            end

            yLim = app.PlotAxes.YLim;
            conditionIDs = unique(double(EM.eventID(EM.state == 1)), 'stable');
            shownCond = false(size(conditionIDs));
            patchList = gobjects(0);

            for iCond = 1:numel(conditionIDs)
                condID = conditionIDs(iCond);

                onMask = EM.state == 1 & double(EM.eventID) == condID;
                offMask = EM.state == 0 & double(EM.eventID) == condID;

                if ~isempty(EM.selectedEvents) && islogical(EM.selectedEvents) && ...
                        isequal(size(EM.selectedEvents), size(EM.eventID))
                    onMask = onMask & EM.selectedEvents;
                    offMask = offMask & EM.selectedEvents;
                end

                onTimes = double(EM.timestamps(onMask));
                offTimes = double(EM.timestamps(offMask));
                repOn = double(EM.repetitionID(onMask));
                repOff = double(EM.repetitionID(offMask));

                for iOn = 1:numel(onTimes)
                    idxOff = find(repOff == repOn(iOn), 1, 'first');
                    if isempty(idxOff)
                        continue
                    end

                    x0 = onTimes(iOn);
                    x1 = offTimes(idxOff);

                    if ~isfinite(x0) || ~isfinite(x1) || x1 <= x0
                        continue
                    end

                    color = app.getEventColor(condID);
                    displayName = app.getEventNameFromID(condID);

                    h = patch(app.PlotAxes, ...
                        [x0 x1 x1 x0], ...
                        [yLim(1) yLim(1) yLim(2) yLim(2)], ...
                        color, ...
                        'FaceAlpha', 0.16, ...
                        'EdgeColor', 'none', ...
                        'HitTest', 'off', ...
                        'PickableParts', 'none', ...
                        'DisplayName', displayName);

                    if shownCond(iCond)
                        h.HandleVisibility = 'off';
                    else
                        h.HandleVisibility = 'on';
                        shownCond(iCond) = true;
                    end

                    patchList(end+1) = h; %#ok<AGROW>
                end
            end

            app.EventPatchHandles = patchList;
            app.stackEventPatchesBottom();

        end

        function drawEventModeEventPatch(app)
            %DRAWEVENTMODEEVENTPATCH Draw selected event patch/onset marker.

            yLim = app.PlotAxes.YLim;
            patchList = gobjects(0);

            switch app.getSourceType()
                case 'dat'
                    duration = app.getSelectedDatEventDuration();
                    if isempty(duration) || ~isfinite(duration) || duration <= 0
                        return
                    end
                    x0 = 0;
                    x1 = duration;
                    displayName = app.CurrentCondition;

                case 'umt'
                    baselinePeriod = app.getActiveBaselinePeriod();
                    if isempty(baselinePeriod)
                        return
                    end

                    frameRateHz = app.getFrameRateHz();
                    if isempty(frameRateHz)
                        halfWidth = 0.01;
                    else
                        halfWidth = 0.5 ./ frameRateHz;
                    end

                    x0 = -halfWidth;
                    x1 = halfWidth;
                    displayName = app.CurrentCondition;

                otherwise
                    return
            end

            h = patch(app.PlotAxes, ...
                [x0 x1 x1 x0], ...
                [yLim(1) yLim(1) yLim(2) yLim(2)], ...
                [0 1 1], ...
                'FaceAlpha', 0.20, ...
                'EdgeColor', 'none', ...
                'HitTest', 'off', ...
                'PickableParts', 'none', ...
                'DisplayName', displayName, ...
                'HandleVisibility', 'on');

            patchList(end+1) = h;
            app.EventPatchHandles = patchList;
            app.stackEventPatchesBottom();

        end

        function deleteEventPatches(app)
            %DELETEEVENTPATCHES Delete current temporal event patches.

            if ~isempty(app.EventPatchHandles)
                for iPatch = 1:numel(app.EventPatchHandles)
                    try
                        if isvalid(app.EventPatchHandles(iPatch))
                            delete(app.EventPatchHandles(iPatch));
                        end
                    catch
                    end
                end
            end

            app.EventPatchHandles = gobjects(0);
            legend(app.PlotAxes, 'off');

        end

        function stackEventPatchesBottom(app)
            %STACKEVENTPATCHESBOTTOM Place event patches behind all other plot objects.

            if isempty(app.EventPatchHandles)
                return
            end

            for iPatch = 1:numel(app.EventPatchHandles)
                try
                    if isvalid(app.EventPatchHandles(iPatch))
                        uistack(app.EventPatchHandles(iPatch), 'bottom');
                    end
                catch
                end
            end

        end

        function updateEventPatchYExtents(app, yLim)
            %UPDATEEVENTPATCHYEXTENTS Resize event patches to current PlotAxes YLim.
            %
            %   updateEventPatchYExtents(app) uses app.PlotAxes.YLim.
            %   updateEventPatchYExtents(app, yLim) uses the supplied limits.

            if nargin < 2 || isempty(yLim)
                if isempty(app.PlotAxes) || ~isvalid(app.PlotAxes)
                    return
                end

                yLim = app.PlotAxes.YLim;
            end

            if numel(yLim) ~= 2 || any(~isfinite(yLim)) || yLim(1) == yLim(2)
                return
            end

            if isempty(app.EventPatchHandles)
                return
            end

            for iPatch = 1:numel(app.EventPatchHandles)
                h = app.EventPatchHandles(iPatch);

                if ~app.isUsableGraphicsHandle(h) || ~isprop(h, 'YData')
                    continue
                end

                try
                    yData = h.YData;

                    if numel(yData) == 4
                        h.YData = [yLim(1) yLim(1) yLim(2) yLim(2)];
                    end
                catch
                end
            end

        end

        function color = getEventColor(app, condID)
            %GETEVENTCOLOR Return color for one condition ID.

            uniqueIDs = unique(app.EventInfo.EventIDs, 'stable');
            idx = find(uniqueIDs == condID, 1, 'first');

            if isempty(idx) || idx > size(app.EventConditionColors, 1)
                color = [0 1 1];
            else
                color = app.EventConditionColors(idx, :);
            end

        end

        function name = getEventNameFromID(app, condID)
            %GETEVENTNAMEFROMID Return condition name for one condition ID.

            name = sprintf('Event %g', condID);

            if ~app.hasNormalizedEvents()
                return
            end

            uniqueIDs = unique(app.EventInfo.EventIDs, 'stable');
            idx = find(uniqueIDs == condID, 1, 'first');

            if ~isempty(idx) && idx <= numel(app.EventInfo.EventNames)
                name = app.EventInfo.EventNames{idx};
            end

        end

        function duration = getSelectedDatEventDuration(app)
            %GETSELECTEDDATEVENTDURATION Return event duration for current .dat selection.

            duration = [];

            if isempty(app.EventSource) || isempty(app.CurrentCondition)
                return
            end

            EM = app.EventSource;
            condID = app.getConditionID(app.CurrentCondition);

            if isempty(condID)
                return
            end

            onMask = EM.state == 1 & double(EM.eventID) == condID;
            offMask = EM.state == 0 & double(EM.eventID) == condID;

            if ~isempty(EM.selectedEvents) && islogical(EM.selectedEvents) && ...
                    isequal(size(EM.selectedEvents), size(EM.eventID))
                onMask = onMask & EM.selectedEvents;
                offMask = offMask & EM.selectedEvents;
            end

            onTimes = double(EM.timestamps(onMask));
            offTimes = double(EM.timestamps(offMask));
            repOn = double(EM.repetitionID(onMask));
            repOff = double(EM.repetitionID(offMask));

            durations = [];

            for iOn = 1:numel(onTimes)
                if ~strcmpi(app.CurrentRepetition, 'AVERAGE')
                    repWanted = str2double(app.CurrentRepetition);
                    if repOn(iOn) ~= repWanted
                        continue
                    end
                end

                idxOff = find(repOff == repOn(iOn), 1, 'first');
                if isempty(idxOff)
                    continue
                end

                thisDuration = offTimes(idxOff) - onTimes(iOn);
                if isfinite(thisDuration) && thisDuration > 0
                    durations(end+1) = thisDuration; %#ok<AGROW>
                end
            end

            if isempty(durations)
                return
            end

            if strcmpi(app.CurrentRepetition, 'AVERAGE')
                duration = mean(durations, 'omitnan');
            else
                duration = durations(1);
            end

        end

        function titleText = getImageTitle(app)
            %GETIMAGETITLE Build image axis title from current state.

            if isempty(app.CurrentFile)
                titleText = 'Image';
                return
            end

            [~, name, ext] = fileparts(app.CurrentFile);
            fileText = [name, ext];

            entryText = '';
            if ~isempty(app.CurrentEntry)
                entryText = sprintf(' | Entry: %s', app.CurrentEntry);
            end

            modeText = '';
            if strcmpi(app.ViewMode, 'event')
                modeText = ' | Event mode';
            end

            sz = app.getDataSize();

            if sz(3) > 1
                titleText = sprintf('%s%s%s | Frame %d/%d', ...
                    fileText, entryText, modeText, app.CurrentFrame, sz(3));
            else
                titleText = sprintf('%s%s%s', fileText, entryText, modeText);
            end

        end

        % =====================================================================
        % Data access
        % =====================================================================

        function frame = getCurrentFrame(app)
            %GETCURRENTFRAME Read current display frame from active backend/view mode.

            sz = app.getDataSize();
            app.CurrentFrame = min(max(round(app.CurrentFrame), 1), max(1, sz(3)));

            if ~strcmpi(app.ViewMode, 'event')
                frame = app.DataSource.getFrame(app.CurrentFrame);
                frame = squeeze(frame);
                return
            end

            switch app.getSourceType()
                case 'dat'
                    frame = app.getCurrentDatEventFrame();

                case 'umt'
                    frame = app.getCurrentUMTEventFrame();

                otherwise
                    frame = app.DataSource.getFrame(app.CurrentFrame);
            end

            frame = squeeze(frame);

        end

        function frame = getCurrentDatEventFrame(app)
            %GETCURRENTDATEVENTFRAME Return current .dat event-mode frame.

            sourceSize = app.getSourceDataSize();
            Ny = sourceSize(1);
            Nx = sourceSize(2);

            if isempty(app.EventFrameMatrix)
                frame = nan(Ny, Nx, 'single');
                return
            end

            colIdx = min(max(app.CurrentFrame, 1), size(app.EventFrameMatrix, 2));
            sourceFrames = app.EventFrameMatrix(:, colIdx);
            sourceFrames = sourceFrames(isfinite(sourceFrames));

            if isempty(sourceFrames)
                frame = nan(Ny, Nx, 'single');
                return
            end

            if ~strcmpi(app.CurrentRepetition, 'AVERAGE')
                frame = app.DataSource.getFrame(sourceFrames(1));
                return
            end

            n = 0;
            meanFrame = [];

            for iFrame = 1:numel(sourceFrames)
                thisFrame = single(app.DataSource.getFrame(sourceFrames(iFrame)));

                if isempty(meanFrame)
                    meanFrame = zeros(size(thisFrame), 'single');
                end

                n = n + 1;
                meanFrame = meanFrame + (thisFrame - meanFrame) ./ n;
            end

            frame = meanFrame;

        end

        function frame = getCurrentUMTEventFrame(app)
            %GETCURRENTUMTEVENTFRAME Return current UMT event-mode frame.

            eIdxList = app.getSelectedUMTEventIndices();

            if isempty(eIdxList)
                frame = app.DataSource.getFrame(app.CurrentFrame, 1);
                return
            end

            if ~strcmpi(app.CurrentRepetition, 'AVERAGE') || numel(eIdxList) == 1
                frame = app.DataSource.getFrame(app.CurrentFrame, eIdxList(1));
                return
            end

            n = 0;
            meanFrame = [];

            for iE = 1:numel(eIdxList)
                thisFrame = single(app.DataSource.getFrame(app.CurrentFrame, eIdxList(iE)));

                if isempty(meanFrame)
                    meanFrame = zeros(size(thisFrame), 'single');
                end

                n = n + 1;
                meanFrame = meanFrame + (thisFrame - meanFrame) ./ n;
            end

            frame = meanFrame;

        end

        function [trace, status, semTrace, xData] = getCurrentPixelTrace(app)
            %GETCURRENTPIXELTRACE Read current crosshair temporal trace.

            status = 'ok';
            semTrace = [];
            xData = [];

            sz = app.getDataSize();

            x = app.CrosshairXY(1);
            y = app.CrosshairXY(2);

            x = min(max(round(x), 1), sz(2));
            y = min(max(round(y), 1), sz(1));

            if strcmpi(app.ViewMode, 'event')
                switch app.getSourceType()
                    case 'dat'
                        [trace, status, semTrace] = app.getDatEventPixelTrace(y, x);
                    case 'umt'
                        [trace, status, semTrace] = app.getUMTEventPixelTrace(y, x);
                    otherwise
                        [trace, status] = app.DataSource.getPixelTrace(y, x);
                end
            else
                [trace, status] = app.DataSource.getPixelTrace(y, x);
            end

            trace = trace(:);
            semTrace = semTrace(:);
            xData = app.getTraceTimeVector(numel(trace));

        end

        function [trace, status, semTrace] = getDatEventPixelTrace(app, y, x)
            %GETDATEVENTPIXELTRACE Return .dat event-mode pixel trace.

            semTrace = [];

            [fullTrace, status] = app.DataSource.getPixelTrace(y, x);
            fullTrace = fullTrace(:);

            if ~strcmp(status, 'ok') && ~strcmp(status, 'cache_rebuilt')
                trace = [];
                return
            end

            if isempty(app.EventFrameMatrix)
                trace = [];
                status = 'no_event_frames';
                return
            end

            frMat = app.EventFrameMatrix;

            if ~strcmpi(app.CurrentRepetition, 'AVERAGE')
                frameIdx = frMat(1, :);
                valid = isfinite(frameIdx);
                frameIdx = frameIdx(valid);
                trace = fullTrace(frameIdx);
                return
            end

            [trace, semTrace] = app.computeWelfordTraceFromFrameMatrix(fullTrace, frMat);

        end

        function [trace, status, semTrace] = getUMTEventPixelTrace(app, y, x)
            %GETUMTEVENTPIXELTRACE Return UMT event-mode pixel trace.

            status = 'ok';
            semTrace = [];
            eIdxList = app.getSelectedUMTEventIndices();

            if isempty(eIdxList)
                [trace, status] = app.DataSource.getPixelTrace(y, x, 1);
                return
            end

            if ~strcmpi(app.CurrentRepetition, 'AVERAGE') || numel(eIdxList) == 1
                [trace, status] = app.DataSource.getPixelTrace(y, x, eIdxList(1));
                return
            end

            n = [];
            meanTrace = [];
            M2 = [];

            for iE = 1:numel(eIdxList)
                thisTrace = app.DataSource.getPixelTrace(y, x, eIdxList(iE));
                thisTrace = double(thisTrace(:));

                if isempty(meanTrace)
                    meanTrace = zeros(size(thisTrace));
                    M2 = zeros(size(thisTrace));
                    n = zeros(size(thisTrace));
                end

                valid = isfinite(thisTrace);
                n(valid) = n(valid) + 1;

                delta = thisTrace(valid) - meanTrace(valid);
                meanTrace(valid) = meanTrace(valid) + delta ./ n(valid);
                delta2 = thisTrace(valid) - meanTrace(valid);
                M2(valid) = M2(valid) + delta .* delta2;
            end

            trace = meanTrace;
            semTrace = zeros(size(meanTrace));
            validSEM = n > 1;
            semTrace(validSEM) = sqrt((M2(validSEM) ./ (n(validSEM) - 1)) ./ n(validSEM));
            semTrace(n <= 1) = 0;

        end

        function [meanTrace, semTrace] = computeWelfordTraceFromFrameMatrix(app, fullTrace, frMat)
            %COMPUTEWELFORDTRACEFROMFRAMEMATRIX Mean/SEM across event repetitions.

            if isempty(frMat)
                meanTrace = [];
                semTrace = [];
                return
            end

            % AVERAGE mode always crops to the shortest valid trial.
            firstInvalidCol = find(any(isnan(frMat), 1), 1, 'first');
            if ~isempty(firstInvalidCol)
                frMat = frMat(:, 1:firstInvalidCol-1);
            end

            if isempty(frMat) || size(frMat, 2) == 0
                meanTrace = [];
                semTrace = [];
                return
            end

            nFrames = size(frMat, 2);
            n = zeros(1, nFrames);
            meanTrace = zeros(1, nFrames);
            M2 = zeros(1, nFrames);

            for iRep = 1:size(frMat, 1)
                frameIdx = frMat(iRep, :);
                valid = isfinite(frameIdx);

                x = nan(1, nFrames);
                x(valid) = double(fullTrace(frameIdx(valid)));

                valid = isfinite(x);
                n(valid) = n(valid) + 1;

                delta = x(valid) - meanTrace(valid);
                meanTrace(valid) = meanTrace(valid) + delta ./ n(valid);
                delta2 = x(valid) - meanTrace(valid);
                M2(valid) = M2(valid) + delta .* delta2;
            end

            semTrace = zeros(size(meanTrace));
            validSEM = n > 1;
            semTrace(validSEM) = sqrt((M2(validSEM) ./ (n(validSEM) - 1)) ./ n(validSEM));
            semTrace(n <= 1) = 0;

            meanTrace = meanTrace(:);
            semTrace = semTrace(:);

        end

        % =====================================================================
        % Cache and image interaction
        % =====================================================================

        function refreshCacheOverlay(app)
            %REFRESHCACHEOVERLAY Update cache rectangle overlay.

            if isempty(app.CacheRectHandle) || ~isvalid(app.CacheRectHandle)
                return
            end

            if ~app.hasData() || ~app.ShowCacheRectangle || ~app.sourceUsesTemporalCache()
                app.CacheRectHandle.Visible = 'off';
                return
            end

            rectPos = app.DataSource.getCacheRectangle();

            if isempty(rectPos)
                app.CacheRectHandle.Visible = 'off';
                return
            end

            set(app.CacheRectHandle, ...
                'Position', rectPos, ...
                'Visible', 'on');

        end

        function refreshCacheMenuState(app)
            %REFRESHCACHEMENUSTATE Update context-menu text for cache state.

            app.configureCacheContextMenu();

            if isempty(app.ToggleCacheLockMenu) || ~isvalid(app.ToggleCacheLockMenu)
                return
            end

            if ~app.sourceUsesTemporalCache()
                app.UpdateCacheMenu.Enable = 'off';
                app.ToggleCacheLockMenu.Enable = 'off';
                app.ToggleCacheAreaMenu.Enable = 'off';
                return
            end

            app.UpdateCacheMenu.Enable = 'on';
            app.ToggleCacheLockMenu.Enable = 'on';
            app.ToggleCacheAreaMenu.Enable = 'on';

            if app.isCacheLocked()
                app.ToggleCacheLockMenu.Text = 'Unlock cache';
            else
                app.ToggleCacheLockMenu.Text = 'Lock cache';
            end

            if app.ShowCacheRectangle
                app.ToggleCacheAreaMenu.Text = 'Hide cache area';
            else
                app.ToggleCacheAreaMenu.Text = 'Show cache area';
            end

        end

        function configureCacheContextMenu(app)
            %CONFIGURECACHECONTEXTMENU Attach cache context menu only when useful.

            if ~app.hasData() || ~app.sourceUsesTemporalCache()
                app.ImageAxes.ContextMenu = [];

                if ~isempty(app.ImageHandle) && isvalid(app.ImageHandle)
                    app.ImageHandle.ContextMenu = [];
                end

                return
            end

            app.ImageAxes.ContextMenu = app.ImageContextMenu;

            if ~isempty(app.ImageHandle) && isvalid(app.ImageHandle)
                app.ImageHandle.ContextMenu = app.ImageContextMenu;
            end

        end

        function tf = sourceUsesTemporalCache(app)
            %SOURCEUSESTEMPORALCACHE Return true only for partial temporal-cache mode.

            tf = false;

            if ~app.hasData()
                return
            end

            if ismethod(app.DataSource, 'hasPartialTemporalCache')
                tf = app.DataSource.hasPartialTemporalCache();
                return
            end

            % Fallback for older backend objects: infer from cache rectangle.
            if ~ismethod(app.DataSource, 'getCacheRectangle')
                return
            end

            rectPos = app.DataSource.getCacheRectangle();

            if isempty(rectPos)
                return
            end

            sz = app.getDataSize();

            coversFullImage = ...
                rectPos(1) <= 0.5 && ...
                rectPos(2) <= 0.5 && ...
                rectPos(1) + rectPos(3) >= sz(2) + 0.5 && ...
                rectPos(2) + rectPos(4) >= sz(1) + 0.5;

            tf = ~coversFullImage;

        end

        function tf = isCacheLocked(app)
            %ISCACHELOCKED Return true when backend cache is locked.

            tf = false;

            if ~app.hasData()
                return
            end

            if ismethod(app.DataSource, 'isCacheLocked')
                tf = app.DataSource.isCacheLocked();
            end

        end

        function [x, y] = clampImageCoordinates(app, x, y)
            %CLAMPIMAGECOORDINATES Clamp coordinates to image bounds.

            sz = app.getDataSize();

            x = min(max(round(x), 1), sz(2));
            y = min(max(round(y), 1), sz(1));

        end

        function onImageClicked(app, ~, event)
            %ONIMAGECLICKED Handle image/axis click interaction.
            %
            %   Single-click updates the inspected pixel and temporal profile. Double-click
            %   inside an ROI still starts edit mode as a shortcut, but table Select is now
            %   the primary edit-mode controller.

            if ~app.hasData()
                return
            end

            if isprop(event, 'Button') && event.Button ~= 1
                return
            end

            coords = [];

            try
                if isprop(event, 'IntersectionPoint') && ~isempty(event.IntersectionPoint)
                    coords = round(event.IntersectionPoint(1:2));
                end
            catch
                coords = [];
            end

            if isempty(coords)
                try
                    currentPoint = app.ImageAxes.CurrentPoint;
                    coords = round(currentPoint(1, 1:2));
                catch
                    return
                end
            end

            x = coords(1);
            y = coords(2);

            [x, y] = app.clampImageCoordinates(x, y);

            bDoubleClick = false;

            try
                if isprop(event, 'SelectionType') && strcmpi(event.SelectionType, 'double')
                    bDoubleClick = true;
                end
            catch
                bDoubleClick = false;
            end

            if ~bDoubleClick && ~isempty(app.LastImageClickTic)
                try
                    elapsedSec = toc(app.LastImageClickTic);
                    samePixel = all(abs([x y] - app.LastImageClickXY) <= 3);
                    bDoubleClick = isfinite(elapsedSec) && elapsedSec <= 0.35 && samePixel;
                catch
                    bDoubleClick = false;
                end
            end

            app.LastImageClickTic = tic;
            app.LastImageClickXY = [x y];

            if bDoubleClick
                roiID = app.findVisibleROIAtPixel(x, y);

                if ~isempty(roiID)
                    app.LastImageClickTic = [];
                    app.LastImageClickXY = [NaN NaN];

                    selectedIDs = app.getSelectedROIIDList();

                    if isempty(selectedIDs)
                        app.clearROISelectionState();
                        roiIdx = app.getROIIndexByID(roiID);
                        if ~isempty(roiIdx)
                            app.setROISelectedStateByIndex(roiIdx, true);
                            app.refreshROITable();
                            app.handleROISelectionStateChanged();
                        end
                    elseif numel(selectedIDs) > 1 && any(selectedIDs == roiID)
                        app.startGroupROIEditByIDs(selectedIDs);
                    elseif any(selectedIDs == roiID)
                        app.handleROISelectionStateChanged();
                    else
                        app.clearROISelectionState();
                        roiIdx = app.getROIIndexByID(roiID);
                        if ~isempty(roiIdx)
                            app.setROISelectedStateByIndex(roiIdx, true);
                            app.refreshROITable();
                            app.handleROISelectionStateChanged();
                        end
                    end
                    return
                end
            end

            app.CrosshairXY = [x y];

            app.refreshCrosshairGraphics();
            app.refreshTemporalProfile();

        end

        function startGroupROIEditByIDs(app, roiIDList)
            %STARTGROUPROIEDITBYIDS Edit selected ROIs through one bounding box.

            roiIDList = unique(double(roiIDList(:).'), 'stable');
            roiIDList = roiIDList(isfinite(roiIDList));

            if numel(roiIDList) < 2
                if numel(roiIDList) == 1
                    app.handleROISelectionStateChanged();
                end
                return
            end

            if ~app.hasData()
                return
            end

            if ~strcmp(app.InteractionMode, 'idle')
                app.stopActiveROIEditForSelectionChange();
            end

            roiIdxList = [];
            for iID = 1:numel(roiIDList)
                roiIdx = app.getROIIndexByID(roiIDList(iID));
                if ~isempty(roiIdx)
                    roiIdxList(end+1) = roiIdx; %#ok<AGROW>
                end
            end

            if numel(roiIdxList) < 2
                app.setStatusMessage('Select at least two valid ROIs for group editing.');
                return
            end

            [state, ok, msg] = app.buildGroupROIEditState(roiIdxList);
            if ~ok
                app.setStatusMessage(msg);
                return
            end

            app.deleteGroupEditRuntimeGraphics();

            previousKeyFcn = app.UIFigure.WindowKeyPressFcn;
            state.previousKeyFcn = previousKeyFcn;
            state.previousMode = app.InteractionMode;

            for k = 1:numel(roiIdxList)
                roiIdx = roiIdxList(k);
                if isfield(app.ROIList(roiIdx).runtime, 'ROIHandle') && ...
                        app.isUsableGraphicsHandle(app.ROIList(roiIdx).runtime.ROIHandle)
                    try
                        delete(app.ROIList(roiIdx).runtime.ROIHandle);
                    catch
                    end
                end
                app.ROIList(roiIdx).runtime.ROIHandle = gobjects(1);
                app.setROISelectedStateByIndex(roiIdx, true);
            end

            bbox = state.originalBBox;

            hold(app.ImageAxes, 'on');

            % Group bounding boxes can be translated/resized beyond image limits.
            % Commit-time rasterization determines whether any selected ROI remains valid.
            hBox = drawrectangle(app.ImageAxes, ...
                'Position', bbox, ...
                'Color', [1 1 1], ...
                'FaceAlpha', 0, ...
                'LineWidth', 1.5, ...
                'Rotatable', true, ...
                'DrawingArea', 'unlimited', ...
                'InteractionsAllowed', 'all');
            app.setROIObjectPropertyIfAvailable(hBox, 'LineStyle', '--');
            app.updateGroupEditBoxSizeLabel(hBox);

            listeners = {};
            listeners{end+1} = addlistener(hBox, 'MovingROI', ...
                @(src, evt) app.onGroupEditBoxMoving(src, evt));
            listeners{end+1} = addlistener(hBox, 'ROIMoved', ...
                @(src, evt) app.onGroupEditBoxMoved(src, evt));
            listeners{end+1} = addlistener(hBox, 'ROIClicked', ...
                @(src, evt) app.onGroupEditBoxClicked(src, evt));

            try
                listeners{end+1} = addlistener(hBox, 'DeletingROI', ...
                    @(src, evt) app.onGroupEditBoxDeleting(src, evt));
            catch
            end

            state.listeners = listeners;
            app.GroupROIEditState = state;
            app.GroupEditBoxHandle = hBox;
            % Attach the full dynamic group-edit menu: Delete/Cancel, separator, Boolean Ops.
            app.createGroupROIEditContextMenu();
            app.UIFigure.WindowKeyPressFcn = @(src, evt) app.onActiveGroupROIEditKeyPress(src, evt);

            app.createGroupEditPreviewGraphics();
            app.updateGroupEditPreview();

            app.setInteractionMode('editingROI');
            app.refreshROITable();
            app.setStatusMessage('Group editing ROIs. Move/resize/rotate the box; double-click or press Enter to confirm; uncheck all selected ROIs to exit.');

        end

        function updateGroupEditBoxSizeLabel(app, hBox)
            %UPDATEGROUPEDITBOXSIZELABEL Display current group box size while resizing.

            if ~app.isUsableGraphicsHandle(hBox)
                return
            end

            try
                pos = double(hBox.Position(:).');
            catch
                return
            end

            if numel(pos) < 4 || any(~isfinite(pos(3:4)))
                return
            end

            [scaleFactor, unitText] = app.getROIAxisUnitScale();
            sizeXY = abs(pos(3:4)) .* scaleFactor;

            labelText = sprintf('Group | %.3g x %.3g %s', sizeXY(1), sizeXY(2), unitText);

            app.setROIObjectPropertyIfAvailable(hBox, 'Label', labelText);
            app.setROIObjectPropertyIfAvailable(hBox, 'LabelVisible', 'on');

        end

        function updateGroupEditPreview(app)
            %UPDATEGROUPEDITPREVIEW Apply current group box transform to previews.
            %
            %   Primitive ROIs are previewed using the same approximation that will be
            %   committed. This avoids showing a sheared rectangle/ellipse that later
            %   snaps to an editable primitive after the user clicks it.

            if isempty(fieldnames(app.GroupROIEditState)) || ...
                    ~app.isUsableGraphicsHandle(app.GroupEditBoxHandle)
                return
            end

            transformInfo = app.getCurrentGroupEditTransform();
            roiIdxList = app.GroupROIEditState.roiIdxList;

            for k = 1:numel(roiIdxList)
                if k > numel(app.GroupEditPreviewHandles) || ...
                        ~app.isUsableGraphicsHandle(app.GroupEditPreviewHandles(k))
                    continue
                end

                roiIdx = roiIdxList(k);
                roiType = lower(char(string(app.ROIList(roiIdx).type)));

                verticesXY = app.GroupROIEditState.originalVertices{k};
                verticesAffine = app.transformVerticesFromGroupEdit(verticesXY, transformInfo);

                [verticesPreview, ~] = app.approximateGroupEditedROIPrimitive(roiType, verticesAffine);

                try
                    app.GroupEditPreviewHandles(k).XData = verticesPreview(:, 1);
                    app.GroupEditPreviewHandles(k).YData = verticesPreview(:, 2);
                catch
                end
            end

        end

        function transformInfo = getCurrentGroupEditTransform(app)
            %GETCURRENTGROUPEDITTRANSFORM Return transform from original group box.

            state = app.GroupROIEditState;
            hBox = app.GroupEditBoxHandle;

            pos = double(hBox.Position(:).');
            newSize = max(abs(pos(3:4)), eps);
            newCenter = [pos(1) + pos(3)/2, pos(2) + pos(4)/2];
            angleDeg = app.getNumericROIObjectProperty(hBox, 'RotationAngle', 0);

            transformInfo = struct();
            transformInfo.originalCenter = state.originalCenter;
            transformInfo.newCenter = newCenter;
            transformInfo.scaleXY = newSize ./ max(state.originalSize, eps);
            transformInfo.angleDeg = angleDeg;

        end

        function verticesOut = transformVerticesFromGroupEdit(app, verticesIn, transformInfo) %#ok<INUSL>
            %TRANSFORMVERTICESFROMGROUPEDIT Apply group edit affine transform.

            verticesIn = double(verticesIn);

            localXY = verticesIn - transformInfo.originalCenter;
            localXY(:, 1) = localXY(:, 1) .* transformInfo.scaleXY(1);
            localXY(:, 2) = localXY(:, 2) .* transformInfo.scaleXY(2);

            % Match the image-coordinate rotation convention used for ROI storage.
            angleRad = -deg2rad(transformInfo.angleDeg);
            R = [cos(angleRad), -sin(angleRad); sin(angleRad), cos(angleRad)];

            verticesOut = localXY * R.' + transformInfo.newCenter;

        end

        function commitGroupROIEdit(app)
            %COMMITGROUPROIEDIT Store transformed geometry for all selected ROIs.
            %
            %   The group edit box is intentionally unbounded. At commit time, each
            %   transformed ROI is rasterized into image coordinates and clipped to the
            %   active logical mask. ROIs that no longer contain valid pixels are deleted
            %   only after user confirmation. Split clipped ROIs are preserved as
            %   <ROINAME>_part<N> polygon ROIs.

            if isempty(fieldnames(app.GroupROIEditState)) || ...
                    ~isfield(app.GroupROIEditState, 'roiIdxList')
                return
            end

            transformInfo = app.getCurrentGroupEditTransform();
            roiIdxList = app.GroupROIEditState.roiIdxList;
            editedIDs = app.GroupROIEditState.roiIDList;

            nROI = numel(roiIdxList);
            editResults = repmat(struct( ...
                'roiIdx', [], ...
                'roiID', [], ...
                'action', '', ...
                'type', '', ...
                'vertices', [], ...
                'mask', [], ...
                'pgon', [], ...
                'params', [], ...
                'componentMasks', {{}}, ...
                'baseName', '', ...
                'bUsePartNames', false), 1, nROI);

            invalidROIIDs = [];
            invalidROINames = {};
            nClipped = 0;
            nSplitOriginalROIs = 0;

            % -------------------------------------------------------------------------
            % First pass: compute and validate all transformed/clipped outputs.
            % -------------------------------------------------------------------------
            for k = 1:nROI
                roiIdx = roiIdxList(k);

                if roiIdx < 1 || roiIdx > numel(app.ROIList)
                    app.setStatusMessage('Group ROI edit cancelled: selected ROI index is invalid.');
                    app.cancelGroupROIEdit(false);
                    return
                end

                roiID = app.ROIList(roiIdx).ID;
                roiType = lower(char(string(app.ROIList(roiIdx).type)));
                baseName = char(string(app.ROIList(roiIdx).name));

                verticesXY = app.GroupROIEditState.originalVertices{k};
                verticesAffine = app.cleanROIVertices( ...
                    app.transformVerticesFromGroupEdit(verticesXY, transformInfo));

                [verticesStored, paramsStored] = app.approximateGroupEditedROIPrimitive( ...
                    roiType, verticesAffine);

                verticesStored = app.cleanROIVertices(verticesStored);

                editResults(k).roiIdx = roiIdx;
                editResults(k).roiID = roiID;
                editResults(k).baseName = baseName;

                if size(verticesStored, 1) < 3
                    editResults(k).action = 'delete';
                    invalidROIIDs(end+1) = roiID; %#ok<AGROW>
                    invalidROINames{end+1} = baseName; %#ok<AGROW>
                    continue
                end

                rawMask = app.createMaskFromVertices(verticesStored);
                if isempty(rawMask) || ~any(rawMask(:))
                    editResults(k).action = 'delete';
                    invalidROIIDs(end+1) = roiID; %#ok<AGROW>
                    invalidROINames{end+1} = baseName; %#ok<AGROW>
                    continue
                end

                mask = app.clipROIMaskToActiveLogicalMask(rawMask);
                if isempty(mask) || ~any(mask(:))
                    editResults(k).action = 'delete';
                    invalidROIIDs(end+1) = roiID; %#ok<AGROW>
                    invalidROINames{end+1} = baseName; %#ok<AGROW>
                    continue
                end

                if ~isequal(mask, rawMask)
                    nClipped = nClipped + 1;

                    componentMasks = app.maskToConnectedComponentMasks(mask);
                    if isempty(componentMasks)
                        editResults(k).action = 'delete';
                        invalidROIIDs(end+1) = roiID; %#ok<AGROW>
                        invalidROINames{end+1} = baseName; %#ok<AGROW>
                        continue
                    end

                    editResults(k).action = 'replace';
                    editResults(k).componentMasks = componentMasks;
                    editResults(k).bUsePartNames = numel(componentMasks) > 1;

                    if numel(componentMasks) > 1
                        nSplitOriginalROIs = nSplitOriginalROIs + 1;
                    end

                    % Validate that every component can become a polygon before changing ROIList.
                    for iComp = 1:numel(componentMasks)
                        verticesComp = app.maskToSimplifiedPolygonVertices(componentMasks{iComp});
                        verticesComp = app.cleanROIVertices(verticesComp);

                        if size(verticesComp, 1) < 3
                            app.setStatusMessage('Group ROI edit cancelled: one split ROI component has invalid geometry.');
                            app.cancelGroupROIEdit(false);
                            return
                        end

                        pgonComp = polyshape(verticesComp(:, 1), verticesComp(:, 2), 'Simplify', true);
                        if isempty(pgonComp.Vertices)
                            app.setStatusMessage('Group ROI edit cancelled: one split ROI component produced an invalid polyshape.');
                            app.cancelGroupROIEdit(false);
                            return
                        end
                    end

                    continue
                end

                % No mask clipping: keep the transformed primitive/polygon approximation.
                pgon = polyshape(verticesStored(:, 1), verticesStored(:, 2), 'Simplify', true);
                if isempty(pgon.Vertices)
                    editResults(k).action = 'delete';
                    invalidROIIDs(end+1) = roiID; %#ok<AGROW>
                    invalidROINames{end+1} = baseName; %#ok<AGROW>
                    continue
                end

                editResults(k).action = 'update';
                editResults(k).type = roiType;
                editResults(k).vertices = verticesStored;
                editResults(k).mask = rawMask;
                editResults(k).pgon = pgon;
                editResults(k).params = paramsStored;
            end

            invalidROIIDs = unique(double(invalidROIIDs(:).'), 'stable');
            invalidROIIDs = invalidROIIDs(isfinite(invalidROIIDs));

            if ~isempty(invalidROIIDs)
                if numel(invalidROIIDs) == 1
                    promptText = sprintf( ...
                        ['The group edit moved ROI "%s" completely outside the image ' ...
                        'or active logical mask. Delete this invalid ROI?'], ...
                        invalidROINames{1});
                else
                    promptText = sprintf( ...
                        ['The group edit moved %d ROIs completely outside the image ' ...
                        'or active logical mask. Delete these invalid ROIs?'], ...
                        numel(invalidROIIDs));
                end

                if ~app.confirmROIDeletion(promptText)
                    app.cancelGroupROIEdit(false);
                    app.setStatusMessage('Group ROI edit cancelled. Invalid ROI deletion was not confirmed.');
                    return
                end
            end

            % -------------------------------------------------------------------------
            % Second pass: remove temporary edit graphics and apply data-model changes.
            % Apply from high index to low index to keep original indices valid when split
            % replacements insert additional ROIList entries.
            % -------------------------------------------------------------------------
            app.deleteGroupEditRuntimeGraphics();

            validResults = editResults(~strcmp({editResults.action}, 'delete'));
            nUpdatedOutputROIs = 0;
            nReplacementFailures = 0;

            if ~isempty(validResults)
                [~, orderDesc] = sort([validResults.roiIdx], 'descend');

                for iOrder = 1:numel(orderDesc)
                    r = validResults(orderDesc(iOrder));
                    roiIdx = r.roiIdx;

                    if roiIdx < 1 || roiIdx > numel(app.ROIList)
                        nReplacementFailures = nReplacementFailures + 1;
                        continue
                    end

                    switch r.action
                        case 'replace'
                            nParts = app.replaceROIByMaskComponents( ...
                                roiIdx, ...
                                r.componentMasks, ...
                                r.baseName, ...
                                r.bUsePartNames);

                            if nParts > 0
                                nUpdatedOutputROIs = nUpdatedOutputROIs + nParts;
                            else
                                nReplacementFailures = nReplacementFailures + 1;
                            end

                        case 'update'
                            app.ROIList(roiIdx).type = r.type;
                            app.ROIList(roiIdx).geometry.polyshape = r.pgon;
                            app.ROIList(roiIdx).geometry.verticesXY_px = r.vertices;
                            app.ROIList(roiIdx).geometry.ROIType = r.type;
                            app.ROIList(roiIdx).geometry.ROIParameters = r.params;

                            app.ROIList(roiIdx).mask = r.mask;
                            app.ROIList(roiIdx).stats = app.computeROIStatsFromMask(r.mask);
                            app.ROIList(roiIdx).modifiedOn = datetime('now');

                            app.setROISelectedStateByIndex(roiIdx, false);

                            app.ROIList(roiIdx).runtime.ROIHandle = ...
                                app.createStaticROIOverlayFromROI(app.ROIList(roiIdx));

                            nUpdatedOutputROIs = nUpdatedOutputROIs + 1;
                    end
                end
            end

            if ~isempty(invalidROIIDs)
                app.deleteROIsByIDList(invalidROIIDs, false);
            end

            app.clearROISelectionState();
            app.setSelectedROI(NaN);
            app.setInteractionMode('idle');

            app.refreshROITable();
            app.refreshROITraces();
            app.refreshEventPatches();
            app.stackROITraceGraphics();
            app.updateGUIEnabledState();

            nDeletedInvalid = numel(invalidROIIDs);

            if nReplacementFailures > 0
                app.setStatusMessage(sprintf( ...
                    'Updated group ROI edit with %d replacement failure(s). Review ROI table.', ...
                    nReplacementFailures));
            elseif nDeletedInvalid > 0
                app.setStatusMessage(sprintf( ...
                    'Updated group ROI edit. Deleted %d invalid ROI(s); kept %d valid output ROI(s).', ...
                    nDeletedInvalid, nUpdatedOutputROIs));
            elseif nSplitOriginalROIs > 0
                app.setStatusMessage(sprintf( ...
                    ['Updated %d selected ROIs into %d output ROIs. ' ...
                    '%d ROI(s) were clipped and %d ROI(s) were split by the logical mask.'], ...
                    numel(editedIDs), nUpdatedOutputROIs, nClipped, nSplitOriginalROIs));
            elseif nClipped > 0
                app.setStatusMessage(sprintf( ...
                    'Updated %d ROIs. %d ROI(s) were clipped to the logical mask.', ...
                    numel(editedIDs), nClipped));
            else
                app.setStatusMessage(sprintf('Updated %d ROIs.', numel(editedIDs)));
            end

        end

        function [verticesOut, params] = approximateGroupEditedROIPrimitive(app, roiType, verticesXY)
            %APPROXIMATEGROUPEDITEDROIPRIMITIVE Return editable primitive geometry.
            %
            %   Group resizing/rotation can shear rectangle and ellipse vertices. MATLAB
            %   editable rectangle/ellipse ROI objects cannot represent shear. For those
            %   ROI types, this method converts the affine-transformed vertices to the
            %   closest editable primitive before committing or previewing the result.
            %
            %   Polygon ROIs keep the exact transformed vertices.

            roiType = lower(char(string(roiType)));
            verticesXY = app.cleanROIVertices(verticesXY);

            params = struct();
            params.ROIType = roiType;
            params.GroupTransformed = true;
            params.PrimitiveApproximation = false;

            verticesOut = verticesXY;

            if isempty(verticesXY) || size(verticesXY, 2) ~= 2 || size(verticesXY, 1) < 3
                verticesOut = zeros(0, 2);
                return
            end

            switch roiType
                case 'rectangle'
                    if size(verticesXY, 1) < 4
                        return
                    end

                    % Use the transformed quadrilateral edges to build the closest rotated
                    % rectangle. This intentionally discards shear.
                    v = verticesXY(1:4, :);
                    centerXY = mean(v, 1);

                    edge12 = v(2, :) - v(1, :);
                    edge23 = v(3, :) - v(2, :);
                    edge34 = v(3, :) - v(4, :);
                    edge41 = v(4, :) - v(1, :);

                    widthValue = mean([norm(edge12), norm(edge34)]);
                    heightValue = mean([norm(edge23), norm(edge41)]);

                    if ~isfinite(widthValue) || ~isfinite(heightValue) || ...
                            widthValue <= 0 || heightValue <= 0
                        verticesOut = zeros(0, 2);
                        return
                    end

                    if norm(edge12) > eps
                        axisX = edge12 ./ norm(edge12);
                    else
                        axisX = [1 0];
                    end

                    % Stored ROI vertices use the image-coordinate rotation convention.
                    storageAngleDeg = atan2d(axisX(2), axisX(1));
                    rotationAngleDeg = -storageAngleDeg;

                    pos = [ ...
                        centerXY(1) - widthValue/2, ...
                        centerXY(2) - heightValue/2, ...
                        widthValue, ...
                        heightValue];

                    baseVertices = [ ...
                        pos(1),             pos(2); ...
                        pos(1) + pos(3),    pos(2); ...
                        pos(1) + pos(3),    pos(2) + pos(4); ...
                        pos(1),             pos(2) + pos(4)];

                    verticesOut = iRotateForROIStorage(baseVertices, centerXY, rotationAngleDeg);

                    params.Position = pos;
                    params.RotationAngle = rotationAngleDeg;
                    params.Vertices = verticesOut;
                    params.PrimitiveApproximation = true;

                case 'ellipse'
                    centerXY = mean(verticesXY, 1);
                    centered = verticesXY - centerXY;

                    if size(centered, 1) < 3
                        verticesOut = zeros(0, 2);
                        return
                    end

                    C = (centered.' * centered) ./ max(1, size(centered, 1) - 1);
                    [V, D] = eig(C);
                    [~, order] = sort(diag(D), 'descend');
                    V = V(:, order);

                    axisX = V(:, 1).';
                    axisY = V(:, 2).';

                    projX = centered * axisX.';
                    projY = centered * axisY.';

                    semiAxes = [range(projX), range(projY)] ./ 2;

                    if any(~isfinite(semiAxes)) || any(semiAxes <= 0)
                        verticesOut = zeros(0, 2);
                        return
                    end

                    storageAngleDeg = atan2d(axisX(2), axisX(1));
                    rotationAngleDeg = -storageAngleDeg;

                    theta = linspace(0, 2*pi, 129).';
                    theta(end) = [];

                    baseVertices = [ ...
                        centerXY(1) + semiAxes(1) .* cos(theta), ...
                        centerXY(2) + semiAxes(2) .* sin(theta)];

                    verticesOut = iRotateForROIStorage(baseVertices, centerXY, rotationAngleDeg);

                    params.Center = centerXY;
                    params.SemiAxes = semiAxes;
                    params.RotationAngle = rotationAngleDeg;
                    params.Position = [ ...
                        centerXY(1) - semiAxes(1), ...
                        centerXY(2) - semiAxes(2), ...
                        2 * semiAxes(1), ...
                        2 * semiAxes(2)];
                    params.Vertices = verticesOut;
                    params.PrimitiveApproximation = true;

                otherwise
                    % Polygon ROIs preserve the exact group affine transform.
                    params.Vertices = verticesXY;
                    params.Position = verticesXY;
                    verticesOut = verticesXY;
            end

            function verticesRot = iRotateForROIStorage(verticesIn, centerIn, rotationAngleDegIn)
                if isempty(verticesIn) || rotationAngleDegIn == 0
                    verticesRot = verticesIn;
                    return
                end

                angleRad = -deg2rad(rotationAngleDegIn);
                R = [cos(angleRad), -sin(angleRad); sin(angleRad), cos(angleRad)];

                verticesRot = (verticesIn - centerIn) * R.' + centerIn;
            end

        end

        function cancelGroupROIEdit(app, bSilent)
            %CANCELGROUPROIEDIT Cancel group edit and restore static overlays.

            if nargin < 2
                bSilent = false;
            end

            roiIdxList = [];
            if isfield(app.GroupROIEditState, 'roiIdxList')
                roiIdxList = app.GroupROIEditState.roiIdxList;
            end

            previousMode = 'idle';
            if isfield(app.GroupROIEditState, 'previousMode') && ...
                    ~isempty(app.GroupROIEditState.previousMode)
                previousMode = app.GroupROIEditState.previousMode;
            end

            app.deleteGroupEditRuntimeGraphics();

            for k = 1:numel(roiIdxList)
                roiIdx = roiIdxList(k);
                if roiIdx >= 1 && roiIdx <= numel(app.ROIList)
                    app.setROISelectedStateByIndex(roiIdx, false);
                    if ~isfield(app.ROIList(roiIdx).runtime, 'ROIHandle') || ...
                            ~app.isUsableGraphicsHandle(app.ROIList(roiIdx).runtime.ROIHandle)
                        app.ROIList(roiIdx).runtime.ROIHandle = app.createStaticROIOverlayFromROI(app.ROIList(roiIdx));
                    end
                end
            end

            app.setSelectedROI(NaN);

            try
                app.setInteractionMode(previousMode);
            catch
                app.setInteractionMode('idle');
            end

            app.refreshROITable();
            app.updateGUIEnabledState();

            if ~bSilent
                app.setStatusMessage('Group ROI edit cancelled.');
            end

        end

        function deleteGroupEditRuntimeGraphics(app)
            %DELETEGROUPEDITRUNTIMEGRAPHICS Delete temporary group edit handles.
            app.deleteGroupROIEditContextMenu();
            try
                if isfield(app.GroupROIEditState, 'listeners') && ~isempty(app.GroupROIEditState.listeners)
                    listeners = app.GroupROIEditState.listeners;
                    for iListener = 1:numel(listeners)
                        try
                            if ~isempty(listeners{iListener}) && isvalid(listeners{iListener})
                                delete(listeners{iListener});
                            end
                        catch
                        end
                    end
                end
            catch
            end

            try
                if isfield(app.GroupROIEditState, 'previousKeyFcn') && isvalid(app.UIFigure)
                    app.UIFigure.WindowKeyPressFcn = app.GroupROIEditState.previousKeyFcn;
                end
            catch
            end

            try
                if app.isUsableGraphicsHandle(app.GroupEditBoxHandle)
                    delete(app.GroupEditBoxHandle);
                end
            catch
            end

            app.GroupEditBoxHandle = gobjects(1);
            app.deleteGroupEditPreviewGraphics();
            app.GroupROIEditState = struct();

        end

        function deleteGroupEditPreviewGraphics(app)
            %DELETEGROUPEDITPREVIEWGRAPHICS Delete temporary group preview patches.

            if ~isempty(app.GroupEditPreviewHandles)
                for iHandle = 1:numel(app.GroupEditPreviewHandles)
                    try
                        if app.isUsableGraphicsHandle(app.GroupEditPreviewHandles(iHandle))
                            delete(app.GroupEditPreviewHandles(iHandle));
                        end
                    catch
                    end
                end
            end

            app.GroupEditPreviewHandles = gobjects(0);

        end

        function onGroupEditBoxMoving(app, src, ~)
            %ONGROUPEDITBOXMOVING Update group ROI preview while moving/resizing box.

            app.updateGroupEditBoxSizeLabel(src);
            app.updateGroupEditPreview();

        end

        function onGroupEditBoxMoved(app, src, ~)
            %ONGROUPEDITBOXMOVED Update group ROI preview after moving/resizing box.

            app.updateGroupEditBoxSizeLabel(src);
            app.updateGroupEditPreview();

        end

        function onGroupEditBoxClicked(app, ~, evt)
            %ONGROUPEDITBOXCLICKED Confirm group edit on double-click.

            try
                if isprop(evt, 'SelectionType') && strcmpi(evt.SelectionType, 'double')
                    app.commitGroupROIEdit();
                end
            catch
            end

        end

        function onGroupEditBoxDeleting(app, ~, ~)
            %ONGROUPEDITBOXDELETING Cancel group edit if box is deleted.

            app.cancelGroupROIEdit(false);

        end

        function onActiveGroupROIEditKeyPress(app, ~, evt)
            %ONACTIVEGROUPROIEDITKEYPRESS Handle group edit keyboard shortcuts.

            switch lower(evt.Key)
                case {'return', 'enter'}
                    app.commitGroupROIEdit();

                case 'escape'
                    app.cancelGroupROIEdit(false);
            end

        end

        function mask = createMaskFromVertices(app, verticesXY)
            %CREATEMASKFROMVERTICES Rasterize polygon vertices to current image size.

            mask = [];

            if isempty(verticesXY) || size(verticesXY, 2) ~= 2 || ~app.hasData()
                return
            end

            sz = app.getDataSize();
            Ny = sz(1);
            Nx = sz(2);

            x = verticesXY(:, 1);
            y = verticesXY(:, 2);

            try
                mask = poly2mask(x, y, Ny, Nx);
            catch
                [X, Y] = meshgrid(1:Nx, 1:Ny);
                mask = inpolygon(X, Y, x, y);
            end

            mask = logical(mask);

        end

        function createGroupEditPreviewGraphics(app)
            %CREATEGROUPEDITPREVIEWGRAPHICS Create temporary preview patches.

            app.deleteGroupEditPreviewGraphics();

            if isempty(fieldnames(app.GroupROIEditState)) || ...
                    ~isfield(app.GroupROIEditState, 'roiIdxList')
                return
            end

            roiIdxList = app.GroupROIEditState.roiIdxList;
            previewHandles = gobjects(0);

            hold(app.ImageAxes, 'on');

            for k = 1:numel(roiIdxList)
                roiIdx = roiIdxList(k);
                verticesXY = app.GroupROIEditState.originalVertices{k};

                roiColor = [1 0 0];
                if isfield(app.ROIList(roiIdx), 'color') && numel(app.ROIList(roiIdx).color) == 3
                    roiColor = min(max(double(app.ROIList(roiIdx).color(:).'), 0), 1);
                end

                h = patch(app.ImageAxes, ...
                    'XData', verticesXY(:, 1), ...
                    'YData', verticesXY(:, 2), ...
                    'FaceColor', roiColor, ...
                    'FaceAlpha', 0.18, ...
                    'EdgeColor', roiColor, ...
                    'LineWidth', 1.5, ...
                    'HitTest', 'off', ...
                    'PickableParts', 'none', ...
                    'HandleVisibility', 'off');

                previewHandles(end+1) = h; %#ok<AGROW>
            end

            app.GroupEditPreviewHandles = previewHandles;

        end

        function [state, ok, msg] = buildGroupROIEditState(app, roiIdxList)
            %BUILDGROUPROIEDITSTATE Collect original geometry for group ROI edit.

            state = struct();
            ok = false;
            msg = '';

            allVertices = zeros(0, 2);
            originalVertices = cell(1, numel(roiIdxList));

            for k = 1:numel(roiIdxList)
                roiIdx = roiIdxList(k);

                if ~isfield(app.ROIList(roiIdx), 'geometry') || ...
                        ~isfield(app.ROIList(roiIdx).geometry, 'verticesXY_px') || ...
                        isempty(app.ROIList(roiIdx).geometry.verticesXY_px)
                    msg = sprintf('ROI "%s" has no editable vertices.', app.ROIList(roiIdx).name);
                    return
                end

                verticesXY = double(app.ROIList(roiIdx).geometry.verticesXY_px);
                verticesXY = app.cleanROIVertices(verticesXY);

                if size(verticesXY, 1) < 3
                    msg = sprintf('ROI "%s" has invalid geometry.', app.ROIList(roiIdx).name);
                    return
                end

                originalVertices{k} = verticesXY;
                allVertices = [allVertices; verticesXY]; %#ok<AGROW>
            end

            minXY = min(allVertices, [], 1);
            maxXY = max(allVertices, [], 1);
            boxSize = maxXY - minXY;

            if any(~isfinite(boxSize)) || any(boxSize <= 0)
                msg = 'Selected ROIs have an invalid group bounding box.';
                return
            end

            state.roiIdxList = roiIdxList(:).';
            state.roiIDList = [app.ROIList(roiIdxList).ID];
            state.originalVertices = originalVertices;
            state.originalBBox = [minXY(1), minXY(2), boxSize(1), boxSize(2)];
            state.originalCenter = minXY + boxSize ./ 2;
            state.originalSize = boxSize;
            state.listeners = {};

            ok = true;

        end

        function roiID = findVisibleROIAtPixel(app, x, y)
            %FINDVISIBLEROIATPIXEL Return topmost visible ROI containing one pixel.

            roiID = [];

            if isempty(app.ROIList)
                return
            end

            sz = app.getDataSize();
            Ny = sz(1);
            Nx = sz(2);

            x = min(max(round(x), 1), Nx);
            y = min(max(round(y), 1), Ny);

            for iROI = numel(app.ROIList):-1:1
                ROI = app.ROIList(iROI);

                if isfield(ROI, 'runtime') && isfield(ROI.runtime, 'visible') && ...
                        ~ROI.runtime.visible
                    continue
                end

                bInside = false;

                if isfield(ROI, 'mask') && ~isempty(ROI.mask)
                    mask = logical(ROI.mask);
                    if isequal(size(mask), [Ny Nx])
                        bInside = mask(y, x);
                    end
                end

                if ~bInside && isfield(ROI, 'geometry') && ...
                        isfield(ROI.geometry, 'polyshape') && ~isempty(ROI.geometry.polyshape)
                    try
                        bInside = isinterior(ROI.geometry.polyshape, x, y);
                    catch
                        bInside = false;
                    end
                end

                if bInside
                    if isfield(ROI, 'ID') && ~isempty(ROI.ID) && isfinite(ROI.ID)
                        roiID = ROI.ID;
                    else
                        roiID = iROI;
                    end
                    return
                end
            end

        end

        function onUpdateCacheAroundCrosshair(app)
            %ONUPDATECACHEAROUNDCROSSHAIR Manually rebuild cache around selected pixel.

            if ~app.hasData() || ~app.sourceUsesTemporalCache()
                return
            end

            x = app.CrosshairXY(1);
            y = app.CrosshairXY(2);

            app.setStatusMessage('Updating temporal cache around crosshair...');
            drawnow limitrate

            app.DataSource.updateCacheAround(y, x);

            app.refreshCacheOverlay();
            app.refreshTemporalProfile();

            app.setStatusMessage('Temporal cache updated around crosshair.');

        end

        function onToggleCacheLock(app)
            %ONTOGGLECACHELOCK Toggle temporal cache lock state.

            if ~app.hasData() || ~app.sourceUsesTemporalCache()
                return
            end

            newState = ~app.isCacheLocked();
            app.DataSource.setCacheLocked(newState);

            app.refreshCacheMenuState();

            if newState
                app.setStatusMessage('Temporal cache locked.');
            else
                app.setStatusMessage('Temporal cache unlocked.');
            end

            app.refreshTemporalProfile();

        end

        function onToggleCacheRectangle(app)
            %ONTOGGLECACHERECTANGLE Show or hide cache overlay.

            if ~app.hasData() || ~app.sourceUsesTemporalCache()
                return
            end

            app.ShowCacheRectangle = ~app.ShowCacheRectangle;

            app.refreshCacheMenuState();
            app.refreshCacheOverlay();

        end

        % =====================================================================
        % Frame navigation and movie playback
        % =====================================================================

        function setCurrentFrame(app, frameIdx)
            %SETCURRENTFRAME Set current frame and refresh frame-dependent views.
            %
            %   This method intentionally does not refresh the temporal trace. The
            %   trace represents the full time profile of the selected pixel, so frame
            %   navigation only needs to update the image content and timebar. ROI
            %   spatial statistics are refreshed from the newly displayed frame, except
            %   during movie playback where table updates would slow rendering.

            if ~app.hasData()
                return
            end

            sz = app.getDataSize();
            app.CurrentFrame = min(max(round(frameIdx), 1), sz(3));

            app.refreshImageFrame();
            title(app.ImageAxes, app.getImageTitle());
            app.refreshTimeBar();
            app.refreshFrameControls();

            if ~strcmp(app.InteractionMode, 'playingMovie')
                app.updateROIStatsForCurrentFrame();
            end

        end

        function goToPreviousFrame(app)
            %GOTOPREVIOUSFRAME Move to previous frame.

            app.setCurrentFrame(app.CurrentFrame - 1);

        end

        function goToNextFrame(app)
            %GOTONEXTFRAME Move to next frame.

            app.setCurrentFrame(app.CurrentFrame + 1);

        end

        function toggleMoviePlayback(app)
            %TOGGLEMOVIEPLAYBACK Start or stop movie playback.

            if ~app.hasData()
                return
            end

            if strcmp(app.InteractionMode, 'playingMovie')
                app.stopMoviePlayback();
                return
            end

            sz = app.getDataSize();
            if sz(3) <= 1
                app.setStatusMessage('Movie playback is unavailable for single-frame data.');
                return
            end

            speedFactor = app.getMovieSpeedFactor();
            frameRateHz = app.getFrameRateHz();

            if isempty(frameRateHz)
                frameRateHz = 10;
                app.setStatusMessage('Frame-rate metadata not found. Using 10 fps for playback.');
            end

            targetFPS = max(eps, frameRateHz * speedFactor);
            maxDisplayFPS = 50;
            timerPeriod = 1 / min(targetFPS, maxDisplayFPS);

            if isempty(app.MovieTimer) || ~isvalid(app.MovieTimer)
                app.MovieTimer = timer( ...
                    'ExecutionMode', 'fixedRate', ...
                    'BusyMode', 'drop', ...
                    'Name', 'DataViewerMovieTimer', ...
                    'TimerFcn', @(~, ~) app.movieTimerTick());
            end

            if strcmpi(app.MovieTimer.Running, 'on')
                stop(app.MovieTimer);
            end

            app.MovieTimer.Period = timerPeriod;
            app.MovieStartFrame = app.CurrentFrame;
            app.MovieStartTic = tic;
            app.setInteractionMode('playingMovie');

            app.PlayMovieButton.Icon = app.StopIconFile;
            app.updateMovieSpeedLabel();

            start(app.MovieTimer);
            app.setStatusMessage('Movie playback started.');

        end

        function stopMoviePlayback(app)
            %STOPMOVIEPLAYBACK Stop movie playback and restore play icon.

            if ~isempty(app.MovieTimer) && isvalid(app.MovieTimer) && ...
                    strcmpi(app.MovieTimer.Running, 'on')
                stop(app.MovieTimer);
            end

            app.setInteractionMode('idle');

            if isprop(app, 'PlayMovieButton') && ~isempty(app.PlayMovieButton) && ...
                    isvalid(app.PlayMovieButton)
                app.PlayMovieButton.Icon = app.PlayIconFile;
            end

            % Movie ticks intentionally skip ROI table refresh for speed. Update the
            % table once when playback stops so it reflects the final displayed frame.
            app.updateROIStatsForCurrentFrame();

            app.setStatusMessage('Movie playback stopped.');

        end

        function movieTimerTick(app)
            %MOVIETIMERTICK Advance movie using elapsed time and drop-frame timing.
            %
            %   The timer does not advance by exactly one frame per callback. Instead,
            %   it computes the frame that should be visible from elapsed wall-clock
            %   time. If MATLAB cannot render every frame, skipped frames are dropped
            %   and playback timing remains close to the requested movie speed.

            if ~strcmp(app.InteractionMode, 'playingMovie') || ~app.hasData()
                return
            end

            sz = app.getDataSize();
            nFrames = sz(3);

            if nFrames <= 1
                app.stopMoviePlayback();
                return
            end

            frameRateHz = app.getFrameRateHz();
            if isempty(frameRateHz)
                frameRateHz = 10;
            end

            targetFPS = max(eps, frameRateHz * app.getMovieSpeedFactor());
            framesElapsed = floor(toc(app.MovieStartTic) * targetFPS);

            newFrame = mod(app.MovieStartFrame - 1 + framesElapsed, nFrames) + 1;

            if newFrame == app.CurrentFrame
                return
            end

            app.CurrentFrame = newFrame;

            app.refreshImageFrame();
            title(app.ImageAxes, app.getImageTitle());
            app.refreshTimeBar();
            app.refreshFrameControls();

            drawnow limitrate

        end

        function speedFactor = getMovieSpeedFactor(app)
            %GETMOVIESPEEDFACTOR Parse playback speed factor from dropdown value.

            speedFactor = 1;

            if isempty(app.MovieSpeedDropDown) || ~isvalid(app.MovieSpeedDropDown)
                return
            end

            rawValue = char(string(app.MovieSpeedDropDown.Value));
            tok = regexp(rawValue, '([0-9]*\.?[0-9]+)', 'tokens', 'once');

            if isempty(tok)
                return
            end

            parsedValue = str2double(tok{1});

            if isfinite(parsedValue) && parsedValue > 0
                speedFactor = parsedValue;
            end

        end

        function updateMovieSpeedLabel(app)
            %UPDATEMOVIESPEEDLABEL Show target movie progression speed in fps.

            if ~isprop(app, 'MovieSpeedLabel') || isempty(app.MovieSpeedLabel) || ...
                    ~isvalid(app.MovieSpeedLabel)
                return
            end

            frameRateHz = app.getFrameRateHz();
            if isempty(frameRateHz)
                frameRateHz = 10;
            end

            targetFPS = frameRateHz * app.getMovieSpeedFactor();

            if targetFPS >= 100
                labelText = sprintf('%.0f fps', targetFPS);
            elseif targetFPS >= 10
                labelText = sprintf('%.1f fps', targetFPS);
            else
                labelText = sprintf('%.2f fps', targetFPS);
            end

            app.MovieSpeedLabel.Text = labelText;

        end

        function cleanupAppResources(app)
            %CLEANUPAPPRESOURCES Stop and delete runtime resources safely.
            %
            %   This method is intentionally idempotent. It can be called more than once
            %   without error.

            if ~isempty(app.MovieTimer) && isvalid(app.MovieTimer)
                try
                    if strcmp(app.MovieTimer.Running, 'on')
                        stop(app.MovieTimer);
                    end
                catch
                    % Timer may already be partially destroyed.
                end

                try
                    delete(app.MovieTimer);
                catch
                    % Timer may already be deleted.
                end
            end

            app.MovieTimer = timer.empty;

            try
                if ~isempty(app.PlotAxesYLimListener) && isvalid(app.PlotAxesYLimListener)
                    delete(app.PlotAxesYLimListener);
                end
            catch
            end

            app.PlotAxesYLimListener = [];

            try
                app.deleteGroupEditRuntimeGraphics();
            catch
            end

            try
                if isfinite(app.SelectedROIID)
                    roiIdx = app.getROIIndexByID(app.SelectedROIID);
                    if ~isempty(roiIdx)
                        app.cancelROIEditingByIndex(roiIdx, true);
                    end
                end
            catch
            end

            app.InteractionMode = 'idle';

        end

        function populateColormapDropdown(app)
            %POPULATECOLORMAPDROPDOWN Populate dropdown with MATLAB and seaborn maps.
            %
            %   Colormap resources are expected under:
            %       <app folder>/colormaps/seaborn_color_palettes.mat
            %       <app folder>/colormaps/icons

            app.ColormapLibrary = struct();

            appFolder = fileparts(mfilename('fullpath'));
            cmapFolder = fullfile(appFolder, 'colormaps');
            iconFolder = fullfile(cmapFolder, 'icons');
            seabornFile = fullfile(cmapFolder, 'seaborn_color_palettes.mat');

            if ~isfolder(iconFolder)
                mkdir(iconFolder);
            end

            items = {};
            itemKeys = {};
            iconFiles = {};

            % -------------------------------------------------------------------------
            % MATLAB colormaps
            % -------------------------------------------------------------------------
            matlabMaps = {'parula', 'turbo', 'nebula', 'jet', 'hsv', 'hot', 'gray'};
            nColors = 64;
            iconHeight = 18;

            for iMap = 1:numel(matlabMaps)
                mapName = matlabMaps{iMap};

                try
                    cmap = feval(mapName, nColors);
                catch
                    % Some MATLAB versions may not have all maps.
                    continue
                end

                key = matlab.lang.makeValidName(mapName);
                app.ColormapLibrary.(key) = cmap;

                iconFile = fullfile(iconFolder, [key, '.png']);
                if ~isfile(iconFile)
                    img = repmat(reshape(cmap, [1, nColors, 3]), iconHeight, 1, 1);
                    imwrite(img, iconFile);
                end

                items{end+1} = mapName; %#ok<AGROW>
                itemKeys{end+1} = key; %#ok<AGROW>
                iconFiles{end+1} = iconFile; %#ok<AGROW>
            end

            % -------------------------------------------------------------------------
            % Seaborn/matplotlib colormaps exported from Python
            % -------------------------------------------------------------------------
            if isfile(seabornFile)
                S = load(seabornFile, '-mat');

                if isfield(S, 'cmapNames')
                    cmapNames = cellstr(string(S.cmapNames(:)));
                else
                    cmapNames = {};
                end

                if isfield(S, 'cmapDisplayNames')
                    displayNames = cellstr(string(S.cmapDisplayNames(:)));
                else
                    displayNames = cmapNames;
                end

                for iMap = 1:numel(cmapNames)
                    key = matlab.lang.makeValidName(cmapNames{iMap});

                    if ~isfield(S, key)
                        continue
                    end

                    cmap = double(S.(key));

                    if size(cmap, 2) ~= 3 || size(cmap, 1) < 2
                        continue
                    end

                    cmap = min(max(cmap, 0), 1);
                    app.ColormapLibrary.(key) = cmap;

                    iconFile = fullfile(iconFolder, [key, '.png']);

                    items{end+1} = displayNames{iMap}; %#ok<AGROW>
                    itemKeys{end+1} = key; %#ok<AGROW>
                    iconFiles{end+1} = iconFile; %#ok<AGROW>
                end
            else
                app.setStatusMessage('Seaborn colormap file not found. MATLAB colormaps loaded only.');
            end

            if isempty(items)
                items = {'gray'};
                itemKeys = {'gray'};
                app.ColormapLibrary.gray = gray(nColors);
                iconFiles = {''};
            end

            app.ColormapDropDown.Items = items;
            app.ColormapDropDown.ItemsData = itemKeys;

            if ismember('parula', itemKeys)
                app.ColormapDropDown.Value = 'parula';
            else
                app.ColormapDropDown.Value = itemKeys{1};
            end

            % Add item icons when available.
            for iItem = 1:numel(iconFiles)
                if isfile(iconFiles{iItem})
                    styleObj = uistyle( ...
                        'Icon', iconFiles{iItem}, ...
                        'IconAlignment', 'rightmargin');
                    addStyle(app.ColormapDropDown, styleObj, 'item', iItem);
                end
            end

        end

        function applySelectedColormap(app)
            %APPLYSELECTEDCOLORMAP Apply selected colormap to the image axes.

            if isempty(app.ColormapLibrary) || ~isstruct(app.ColormapLibrary)
                return
            end

            key = app.ColormapDropDown.Value;
            key = matlab.lang.makeValidName(char(string(key)));

            if ~isfield(app.ColormapLibrary, key)
                return
            end

            cmap = app.ColormapLibrary.(key);

            if app.InvertCheckBox.Value
                cmap = flipud(cmap);
            end

            app.ImageAxes.Colormap = cmap;

        end

        function setDisplayCLim(app, clim, updateSliderLimits, updateImageCLimOnly)
            %SETDISPLAYCLIM Set image CLim and bind temporal plot YLim to it.
            %
            %   setDisplayCLim(app, clim)
            %   setDisplayCLim(app, clim, updateSliderLimits)
            %
            %   updateSliderLimits should be true only during data loading,
            %   confirmed manual Set Clip, or Auto when the computed limits fall
            %   outside the current slider range.

            arguments
                app
                clim
                updateSliderLimits  = false
                updateImageCLimOnly = false
            end


            clim = double(clim(:).');

            if numel(clim) ~= 2 || any(~isfinite(clim)) || clim(1) >= clim(2)
                error('DataViewer:InvalidCLim', ...
                    'CLim must be finite and satisfy min < max.');
            end

            app.ImageAxes.CLim = clim;
            if updateImageCLimOnly
                return
            end
            app.PlotAxes.YLimMode = 'manual';
            app.PlotAxes.YLim = clim;
            app.updateEventPatchYExtents();
            if isempty(app.ClipSliderRange) || ~isvalid(app.ClipSliderRange)
                return
            end

            if updateSliderLimits
                currentLimits = double(app.ClipSliderRange.Limits);
                currentValue  = double(app.ClipSliderRange.Value);

                % Avoid transient RangeSlider errors when the current Value is
                % outside the target Limits.
                tempLimits = [ ...
                    min([currentLimits(1), clim(1), currentValue(1)]), ...
                    max([currentLimits(2), clim(2), currentValue(2)])];

                if tempLimits(1) == tempLimits(2)
                    tempLimits = tempLimits + [-1, 1];
                end

                app.ClipSliderRange.Limits = tempLimits;
                app.ClipSliderRange.Value  = clim;

                return
            end

            sliderLimits = double(app.ClipSliderRange.Limits);

            if clim(1) >= sliderLimits(1) && clim(2) <= sliderLimits(2)
                app.ClipSliderRange.Value = clim;
            end

        end

        function setAutoCLimFromCurrentFrame(app)
            %SETAUTOCLIMFROMCURRENTFRAME Set CLim from center 75 percent of image.

            if isempty(app.ImageHandle) || ~isvalid(app.ImageHandle)
                return
            end

            frame = app.ImageHandle.CData;

            if isempty(frame) || ~isnumeric(frame)
                return
            end

            [Ny, Nx] = size(frame);

            y1 = max(1, floor(0.125 * Ny) + 1);
            y2 = min(Ny, ceil(0.875 * Ny));
            x1 = max(1, floor(0.125 * Nx) + 1);
            x2 = min(Nx, ceil(0.875 * Nx));

            centerFrame = frame(y1:y2, x1:x2);
            centerFrame = centerFrame(isfinite(centerFrame));

            if isempty(centerFrame)
                app.setStatusMessage('Auto clip failed: no finite pixels in center region.');
                return
            end

            clim = double([min(centerFrame(:)), max(centerFrame(:))]);

            if clim(1) == clim(2)
                clim = clim + [-1, 1];
            end

            sliderLimits = app.ClipSliderRange.Limits;
            bOutsideSliderLimits = clim(1) < sliderLimits(1) || clim(2) > sliderLimits(2);

            % Auto normally preserves slider limits, but expands them if needed to avoid
            % RangeSlider assignment errors.
            app.setDisplayCLim(clim, bOutsideSliderLimits);

        end

        function setCrosshairVisibility(app, tf)
            %SETCROSSHAIRVISIBILITY Show or hide crosshair and pixel trace.
            %
            %   This does not control TimeBarHandle, OriginCrosshairHandles, or
            %   CacheRectHandle.

            tf = logical(tf);

            if ~isempty(app.CrosshairHandles) && all(isvalid(app.CrosshairHandles))
                if tf
                    set(app.CrosshairHandles, 'Visible', 'on');
                else
                    set(app.CrosshairHandles, 'Visible', 'off');
                end
            end

            if isempty(app.CrossTraceHandle) || ~isvalid(app.CrossTraceHandle)
                return
            end

            if ~tf
                if strcmpi(app.CrossTraceHandle.Visible, 'on')
                    app.CrossTraceHandle.Visible = 'off';
                end

                if ~isempty(app.CrossTraceSEMHandle) && isvalid(app.CrossTraceSEMHandle) && ...
                        strcmpi(app.CrossTraceSEMHandle.Visible, 'on')
                    app.CrossTraceSEMHandle.Visible = 'off';
                end
                return
            end

            % Do not revive a trace that was hidden because it has no valid data.
            xData = app.CrossTraceHandle.XData;
            yData = app.CrossTraceHandle.YData;

            hasValidTrace = ~isempty(xData) && ~isempty(yData) && ...
                any(isfinite(xData(:))) && any(isfinite(yData(:)));

            if hasValidTrace
                app.CrossTraceHandle.Visible = 'on';
            end

            if isempty(app.CrossTraceSEMHandle) || ~isvalid(app.CrossTraceSEMHandle)
                return
            end

            patchX = app.CrossTraceSEMHandle.XData;
            patchY = app.CrossTraceSEMHandle.YData;
            hasValidSEM = ~isempty(patchX) && ~isempty(patchY) && ...
                any(isfinite(patchX(:))) && any(isfinite(patchY(:)));

            if hasValidTrace && hasValidSEM
                app.CrossTraceSEMHandle.Visible = 'on';
            end

        end

        % =====================================================================
        % GUI state management
        % =====================================================================

        function setStatusMessage(app, msg)
            %SETSTATUSMESSAGE Update status/message label.

            msg = char(string(msg));

            if ~isempty(app.StatusLabel) && isvalid(app.StatusLabel)
                app.StatusLabel.Text = msg;
                fprintf('[DataViewer] %s\n', msg);
            end

            % drawnow limitrate

        end

        function setInteractionMode(app, mode)
            %SETINTERACTIONMODE Set app interaction mode and refresh GUI availability.

            mode = lower(char(string(mode)));

            validModes = { ...
                'idle', ...
                'loading', ...
                'playingmovie', ...
                'settingclip', ...
                'settingviewcalibration', ...
                'drawingroi', ...
                'editingroi', ...
                'drawinglogicalmask', ...
                'runningpipeline', ...
                'editingevents'};

            if ~ismember(mode, validModes)
                error('DataViewer:InvalidInteractionMode', ...
                    'Invalid interaction mode: "%s".', mode);
            end

            % Keep canonical mixed-case values for readability.
            switch mode
                case 'playingmovie'
                    app.InteractionMode = 'playingMovie';
                case 'settingclip'
                    app.InteractionMode = 'settingClip';
                case 'settingviewcalibration'
                    app.InteractionMode = 'settingViewCalibration';
                case 'drawingroi'
                    app.InteractionMode = 'drawingROI';
                case 'editingroi'
                    app.InteractionMode = 'editingROI';
                case 'drawinglogicalmask'
                    app.InteractionMode = 'drawingLogicalMask';
                case 'runningpipeline'
                    app.InteractionMode = 'runningPipeline';
                case 'editingevents'
                    app.InteractionMode = 'editingEvents';
                otherwise
                    app.InteractionMode = mode;
            end

            app.updateGUIEnabledState();

        end

        function caps = getDataCapabilities(app)
            %GETDATACAPABILITIES Return factual capabilities of the loaded data.

            caps = struct();
            caps.hasData = app.hasData();
            caps.hasTemporalDimension = false;
            caps.hasEventDimension = false;
            caps.hasEventMetadata = false;
            caps.hasROIs = false;
            caps.hasPartialTemporalCache = false;

            if ~caps.hasData
                return
            end

            sourceSize = app.getSourceDataSize();
            if numel(sourceSize) >= 4
                caps.hasEventDimension = sourceSize(4) > 1;
            end

            sz = app.getDataSize();
            if numel(sz) >= 3
                caps.hasTemporalDimension = sz(3) > 1;
            end

            if numel(sz) >= 4
                caps.hasEventDimension = sz(4) > 1;
            end

            % Use normalized event metadata as the GUI source of truth.
            caps.hasEventMetadata = app.hasNormalizedEvents();

            if ismethod(app.DataSource, 'hasPartialTemporalCache')
                caps.hasPartialTemporalCache = app.DataSource.hasPartialTemporalCache();
            else
                caps.hasPartialTemporalCache = app.sourceUsesTemporalCache();
            end

            % ROI subsystem.
            caps.hasROIs = ~isempty(app.ROIList);

        end

        function updateGUIEnabledState(app)
            %UPDATEGUIENABLEDSTATE Enable/disable GUI elements from data state and mode.

            caps = app.getDataCapabilities();

            mode = char(string(app.InteractionMode));
            isIdle = strcmp(mode, 'idle');
            isLoading = strcmp(mode, 'loading');
            isPlaying = strcmp(mode, 'playingMovie');
            isSettingClip = strcmp(mode, 'settingClip');
            isSettingViewCalibration = strcmp(mode, 'settingViewCalibration');
            isDrawingLogicalMask = strcmp(mode, 'drawingLogicalMask');

            sourceType = app.getSourceType();
            isDAT = strcmpi(sourceType, 'dat');
            isUMT = strcmpi(sourceType, 'umt');

            % -------------------------------------------------------------------------
            % Menus
            % -------------------------------------------------------------------------
            if ~isempty(app.OpenMenu) && isvalid(app.OpenMenu)
                app.setComponentEnabled(app.OpenMenu, ...
                    ~isLoading && ~isPlaying && ~isSettingClip && ~isSettingViewCalibration && ~isDrawingLogicalMask);
            end

            if ~isempty(app.SaveMenu) && isvalid(app.SaveMenu)
                app.setComponentEnabled(app.SaveMenu, caps.hasData && isIdle);
            end

            if ~isempty(app.ExportMenu) && isvalid(app.ExportMenu)
                app.setComponentEnabled(app.ExportMenu, caps.hasData && isIdle);
            end

            if ~isempty(app.ImportMenu) && isvalid(app.ImportMenu)
                app.setComponentEnabled(app.ImportMenu, ...
                    ~isLoading && ~isPlaying && ~isSettingClip && ~isSettingViewCalibration && ~isDrawingLogicalMask);
            end

            if ~isempty(app.ROIbyTresholdButton) && isvalid(app.ROIbyTresholdButton)
                app.setComponentEnabled(app.ROIbyTresholdButton, caps.hasData && isIdle);
            end

            if ~isempty(app.ROIAllenButton) && isvalid(app.ROIAllenButton)
                app.setComponentEnabled(app.ROIAllenButton, caps.hasData && isIdle);
            end

            % -------------------------------------------------------------------------
            % TopGrid: frame/movie controls
            % -------------------------------------------------------------------------
            enableFrameControls = caps.hasData && caps.hasTemporalDimension && isIdle;
            enableMovieControl = caps.hasData && caps.hasTemporalDimension && ...
                ~isLoading && ~isSettingClip && ~isSettingViewCalibration && ~isDrawingLogicalMask;

            if ~isempty(app.PreviousFrameButton) && isvalid(app.PreviousFrameButton)
                app.setComponentEnabled(app.PreviousFrameButton, enableFrameControls);
            end

            if ~isempty(app.NextFrameButton) && isvalid(app.NextFrameButton)
                app.setComponentEnabled(app.NextFrameButton, enableFrameControls);
            end

            if ~isempty(app.Slider) && isvalid(app.Slider)
                app.setComponentEnabled(app.Slider, enableFrameControls);
            end

            if ~isempty(app.PlayMovieButton) && isvalid(app.PlayMovieButton)
                app.setComponentEnabled(app.PlayMovieButton, enableMovieControl || isPlaying);
            end

            if ~isempty(app.MovieSpeedDropDown) && isvalid(app.MovieSpeedDropDown)
                app.setComponentEnabled(app.MovieSpeedDropDown, enableMovieControl || isPlaying);
            end

            % -------------------------------------------------------------------------
            % BottomGrid: image display controls
            % -------------------------------------------------------------------------
            enableDisplayControls = caps.hasData && ~isLoading && ~isSettingViewCalibration && ~isDrawingLogicalMask;

            if ~isempty(app.BottomGrid) && isvalid(app.BottomGrid)
                app.setContainerEnabled(app.BottomGrid, enableDisplayControls);
            end

            if ~isempty(app.SetClipButton) && isvalid(app.SetClipButton)
                app.setComponentEnabled(app.SetClipButton, enableDisplayControls && ~isPlaying);
            end

            if ~isempty(app.ClipSliderRange) && isvalid(app.ClipSliderRange)
                app.setComponentEnabled(app.ClipSliderRange, enableDisplayControls);
            end

            % -------------------------------------------------------------------------
            % Utilities-tab launchers
            % -------------------------------------------------------------------------
            if ~isempty(app.ViewCalibButton) && isvalid(app.ViewCalibButton)
                app.setComponentEnabled(app.ViewCalibButton, caps.hasData && isIdle);
            end

            if ~isempty(app.LogicalMaskButton) && isvalid(app.LogicalMaskButton)
                app.setComponentEnabled(app.LogicalMaskButton, caps.hasData && isIdle);
            end

            if ~isempty(app.ImageRefButton) && isvalid(app.ImageRefButton)
                app.setComponentEnabled(app.ImageRefButton, isIdle);
            end

            if ~isempty(app.DataHistoryButton) && isvalid(app.DataHistoryButton)
                app.setComponentEnabled(app.DataHistoryButton, caps.hasData && isIdle);
            end

            % Events Manager should remain reachable for manual event loading/editing.
            % RawFolder-dependent controls should be handled inside DataViewer_EventsManager.
            if ~isempty(app.EventsManagerButton) && isvalid(app.EventsManagerButton)
                app.setComponentEnabled(app.EventsManagerButton, ...
                    caps.hasData && isIdle && strcmpi(app.getSourceType(), 'dat'));
            end

            % Dual-camera coregistration genuinely requires valid raw-folder context
            % and MultiCam metadata.
            if ~isempty(app.OiSDUalCamCoregButton) && isvalid(app.OiSDUalCamCoregButton)
                canOpenDualCamCoreg = false;
                if isIdle
                    [canOpenDualCamCoreg, ~, ~] = app.canOpenOiSDualCamCoreg();
                end
                app.setComponentEnabled(app.OiSDUalCamCoregButton, canOpenDualCamCoreg);
            end

            app.refreshLogicalMaskButtonContextMenuState();
            app.refreshSetReferenceButtonContextMenuState();

            % -------------------------------------------------------------------------
            % Pipeline Tab
            % -------------------------------------------------------------------------
            if ~isempty(app.PipeLauncherButton) && isvalid(app.PipeLauncherButton)
                app.setComponentEnabled(app.PipeLauncherButton, caps.hasData && isIdle);
            end

            % -------------------------------------------------------------------------
            % Temporal profile region
            % -------------------------------------------------------------------------
            enableTemporalProfile = caps.hasData && caps.hasTemporalDimension;

            if ~enableTemporalProfile
                if ~isempty(app.CrossTraceHandle) && isvalid(app.CrossTraceHandle)
                    app.CrossTraceHandle.Visible = 'off';
                end

                if ~isempty(app.CrossTraceSEMHandle) && isvalid(app.CrossTraceSEMHandle)
                    app.CrossTraceSEMHandle.Visible = 'off';
                end

                if ~isempty(app.TimeBarHandle) && isvalid(app.TimeBarHandle)
                    app.TimeBarHandle.Visible = 'off';
                end

                if ~isempty(app.TemporalProfileLabel) && isvalid(app.TemporalProfileLabel)
                    if caps.hasData
                        app.TemporalProfileLabel.Text = 'Temporal Profile (not available)';
                    else
                        app.TemporalProfileLabel.Text = 'Temporal Profile';
                    end
                end
            else
                if ~isempty(app.TemporalProfileLabel) && isvalid(app.TemporalProfileLabel)
                    app.TemporalProfileLabel.Text = 'Temporal Profile';
                end
            end

            % -------------------------------------------------------------------------
            % ROI table region
            % -------------------------------------------------------------------------
            enableROITable = caps.hasData && caps.hasROIs && ...
                (isIdle || strcmp(mode, 'editingROI'));

            if ~isempty(app.UITable) && isvalid(app.UITable)
                app.setComponentEnabled(app.UITable, enableROITable);
            end

            if ~isempty(app.DeleteROIButton) && isvalid(app.DeleteROIButton)
                enableDeleteROI = caps.hasData && caps.hasROIs && ...
                    (isIdle || strcmp(mode, 'editingROI'));
                app.setComponentEnabled(app.DeleteROIButton, enableDeleteROI);
            end

            if ~isempty(app.SaveROIButton) && isvalid(app.SaveROIButton)
                app.setComponentEnabled(app.SaveROIButton, caps.hasData && caps.hasROIs && isIdle);
            end

            if ~isempty(app.LoadROIButton) && isvalid(app.LoadROIButton)
                app.setComponentEnabled(app.LoadROIButton, caps.hasData && isIdle);
            end

            if ~isempty(app.ImportROIButton) && isvalid(app.ImportROIButton)
                app.setComponentEnabled(app.ImportROIButton, caps.hasData && isIdle);
            end

            if ~isempty(app.ExportROIDataButton) && isvalid(app.ExportROIDataButton)
                app.setComponentEnabled(app.ExportROIDataButton, caps.hasData && caps.hasROIs && isIdle);
            end

            app.refreshDeleteROIContextMenuState();

            % Logical-mask drawing launcher.
            if ~isempty(app.LogicalMaskButton) && isvalid(app.LogicalMaskButton)
                app.setComponentEnabled(app.LogicalMaskButton, caps.hasData && isIdle);
            end
            app.refreshLogicalMaskButtonContextMenuState();

            % -------------------------------------------------------------------------
            % Events switch/control state
            % -------------------------------------------------------------------------
            datEventSourceLoaded = app.hasDatEventSource();
            hasIgnoredDatEvents = app.hasIgnoredDatEvents();

            eventsAvailable = caps.hasData && app.hasNormalizedEvents();

            eventsFromDAT = datEventSourceLoaded || ...
                (eventsAvailable && isfield(app.EventInfo, 'SourceType') && ...
                strcmpi(app.EventInfo.SourceType, 'dat'));

            eventsFromUMT = eventsAvailable && ...
                isfield(app.EventInfo, 'SourceType') && ...
                strcmpi(app.EventInfo.SourceType, 'umt');

            isEventMode = strcmpi(app.ViewMode, 'event');

            eventsPanelReachable = eventsAvailable || datEventSourceLoaded || hasIgnoredDatEvents;

            if ~isempty(app.EventSwitchGrid) && isvalid(app.EventSwitchGrid)
                app.setContainerEnabled(app.EventSwitchGrid, eventsPanelReachable && isIdle);
            end

            if ~isempty(app.Switch) && isvalid(app.Switch)
                app.setComponentEnabled(app.Switch, datEventSourceLoaded && isIdle);
            end

            enableEventsGrid = isIdle && ( ...
                (eventsAvailable && isEventMode) || ...
                (datEventSourceLoaded && hasIgnoredDatEvents));

            if ~isempty(app.EventsGrid) && isvalid(app.EventsGrid)
                app.setContainerEnabled(app.EventsGrid, enableEventsGrid);
            end

            if ~isempty(app.ConditionDropDown) && isvalid(app.ConditionDropDown)
                app.setComponentEnabled(app.ConditionDropDown, eventsAvailable && isEventMode && isIdle);
            end

            if ~isempty(app.RepetitionDropDown) && isvalid(app.RepetitionDropDown)
                enableRep = eventsAvailable && isEventMode && isIdle;

                if eventsFromUMT && isfield(app.EventInfo, 'EventAxisMode') && ...
                        strcmpi(app.EventInfo.EventAxisMode, 'aggregated_repetitions')
                    enableRep = false;
                end

                app.setComponentEnabled(app.RepetitionDropDown, enableRep);
            end

            enableDatEventDelete = eventsAvailable && eventsFromDAT && isEventMode && isIdle;

            if ~isempty(app.DeleteConditionButton) && isvalid(app.DeleteConditionButton)
                app.setComponentEnabled(app.DeleteConditionButton, enableDatEventDelete);
            end

            if ~isempty(app.DeleteRepetitionButton) && isvalid(app.DeleteRepetitionButton)
                enableDeleteRep = enableDatEventDelete;

                if ~isempty(app.RepetitionDropDown) && isvalid(app.RepetitionDropDown)
                    repValue = char(string(app.RepetitionDropDown.Value));
                    enableDeleteRep = enableDeleteRep && ~strcmpi(repValue, 'AVERAGE');
                end

                app.setComponentEnabled(app.DeleteRepetitionButton, enableDeleteRep);
            end

            if ~isempty(app.RestoreButton) && isvalid(app.RestoreButton)
                enableRestore = datEventSourceLoaded && isIdle && hasIgnoredDatEvents;
                app.setComponentEnabled(app.RestoreButton, enableRestore);
            end

            if ~isempty(app.EventsLabel) && isvalid(app.EventsLabel)
                if ~caps.hasData
                    app.EventsLabel.Text = 'Events';
                elseif hasIgnoredDatEvents && ~eventsAvailable
                    app.EventsLabel.Text = 'Events (all ignored)';
                elseif ~eventsPanelReachable
                    app.EventsLabel.Text = 'Events (not available)';
                elseif eventsFromUMT
                    app.EventsLabel.Text = 'Events (UMT embedded)';
                elseif eventsFromDAT
                    app.EventsLabel.Text = 'Events';
                elseif isUMT
                    app.EventsLabel.Text = 'Events';
                elseif isDAT
                    app.EventsLabel.Text = 'Events';
                else
                    app.EventsLabel.Text = 'Events';
                end
            end

            drawnow limitrate

        end

        function setContainerEnabled(app, containerObj, tf, varargin)
            %SETCONTAINERENABLED Recursively enable/disable children of a container.
            %
            %   setContainerEnabled(app, containerObj, tf)
            %   setContainerEnabled(app, containerObj, tf, 'exclude', {obj1, obj2})
            %
            %   Components tagged as 'AlwaysEnabled' are skipped.

            p = inputParser;
            p.FunctionName = 'setContainerEnabled';

            addRequired(p, 'containerObj');
            addRequired(p, 'tf', @(x) islogical(x) && isscalar(x));
            addParameter(p, 'exclude', {}, @(x) iscell(x) || isobject(x));

            parse(p, containerObj, tf, varargin{:});

            excludeList = p.Results.exclude;
            if ~iscell(excludeList)
                excludeList = {excludeList};
            end

            if isempty(containerObj) || ~isvalid(containerObj)
                return
            end

            if ~isprop(containerObj, 'Children')
                app.setComponentEnabled(containerObj, tf);
                return
            end

            children = containerObj.Children;

            for iChild = 1:numel(children)
                child = children(iChild);

                if app.isExcludedComponent(child, excludeList)
                    continue
                end

                if isprop(child, 'Tag') && strcmpi(char(string(child.Tag)), 'AlwaysEnabled')
                    continue
                end

                if isprop(child, 'Children') && ~isempty(child.Children)
                    app.setContainerEnabled(child, tf, 'exclude', excludeList);
                end

                app.setComponentEnabled(child, tf);
            end

        end

        function setComponentEnabled(~, componentObj, tf)
            %SETCOMPONENTENABLED Safely set Enable on one component if supported.

            if isempty(componentObj) || ~isvalid(componentObj)
                return
            end

            if isprop(componentObj, 'Tag') && strcmpi(char(string(componentObj.Tag)), 'AlwaysEnabled')
                return
            end

            if ~isprop(componentObj, 'Enable')
                return
            end

            if tf
                componentObj.Enable = 'on';
            else
                componentObj.Enable = 'off';
            end

        end

        function tf = isExcludedComponent(app, componentObj, excludeList) %#ok<INUSL>
            %ISEXCLUDEDCOMPONENT Return true when componentObj is in excludeList.

            tf = false;

            if isempty(excludeList)
                return
            end

            for iEx = 1:numel(excludeList)
                ex = excludeList{iEx};

                if isempty(ex)
                    continue
                end

                for j = 1:numel(ex)
                    try
                        if isequal(componentObj, ex(j))
                            tf = true;
                            return
                        end
                    catch
                    end
                end
            end

        end

        function refreshViewerForModeChange(app)
            %REFRESHVIEWERFORMODECHANGE Refresh display after mode/event selection changes.

            app.refreshFrameControls();
            app.refreshImageFrame();
            title(app.ImageAxes, app.getImageTitle());
            app.refreshTemporalProfile();
            app.refreshEventPatches();

            % Mode, condition, and repetition changes alter both the frame index basis
            % for ROI traces and the currently displayed frame used for ROI statistics.
            app.refreshROITraces();
            app.updateROIStatsForCurrentFrame();

            app.updateGUIEnabledState();

        end

        function openViewCalibrationDialog(app)
            %OPENVIEWCALIBRATIONDIALOG Open modal utility for spatial view calibration.

            if ~app.hasData()
                app.setStatusMessage('Load image data before opening View Calibration.');
                return
            end

            originalDataParams = app.DataParams;
            [pixelRatio, originXY] = app.getCurrentViewCalibration();

            if isempty(pixelRatio)
                pixelRatioText = '0';
            else
                pixelRatioText = sprintf('%.8g', pixelRatio);
            end

            sz = app.getDataSize();
            Ny = sz(1);
            Nx = sz(2);

            previousMode = app.InteractionMode;
            app.setInteractionMode('settingViewCalibration');
            cleanupMode = onCleanup(@() app.setInteractionMode(previousMode));

            dlg = uifigure( ...
                'Name', 'View Calibration', ...
                'WindowStyle', 'modal', ...
                'Position', [100 100 400 260], ...
                'Visible', 'off', ...
                'CloseRequestFcn', @onCancel);

            grid = uigridlayout(dlg);
            grid.RowHeight = {28, 28, 28, 28, 28, '1x', 36};
            grid.ColumnWidth = {120, '1x', 95};
            grid.Padding = [12 12 12 12];

            pixelLabel = uilabel(grid);
            pixelLabel.Text = 'Pixel ratio (px/mm)';
            pixelLabel.Tooltip = {'Pixel ratio in px/mm. Use 0 or [] to disable mm calibration.'};
            pixelLabel.Layout.Row = 1;
            pixelLabel.Layout.Column = 1;

            pixelField = uieditfield(grid, 'text');
            pixelField.Value = pixelRatioText;
            pixelField.Tooltip = {'Enter a positive scalar in px/mm. Use 0 or [] for no pixel ratio.'};
            pixelField.Layout.Row = 1;
            pixelField.Layout.Column = [2 3];

            xLabel = uilabel(grid);
            xLabel.Text = 'Origin X (px)';
            xLabel.Layout.Row = 2;
            xLabel.Layout.Column = 1;

            xField = uieditfield(grid, 'numeric');
            xField.Value = originXY(1);
            xField.Limits = [1 Nx];
            xField.Tooltip = {'Origin X coordinate in image pixels.'};
            xField.Layout.Row = 2;
            xField.Layout.Column = [2 3];

            yLabel = uilabel(grid);
            yLabel.Text = 'Origin Y (px)';
            yLabel.Layout.Row = 3;
            yLabel.Layout.Column = 1;

            yField = uieditfield(grid, 'numeric');
            yField.Value = originXY(2);
            yField.Limits = [1 Ny];
            yField.Tooltip = {'Origin Y coordinate in image pixels.'};
            yField.Layout.Row = 3;
            yField.Layout.Column = [2 3];

            interactiveButton = uibutton(grid, 'push');
            interactiveButton.Text = 'Set origin interactively';
            interactiveButton.Layout.Row = 4;
            interactiveButton.Layout.Column = [1 2];
            interactiveButton.ButtonPushedFcn = @onInteractiveOrigin;

            resetButton = uibutton(grid, 'push');
            resetButton.Text = 'Reset origin';
            resetButton.Layout.Row = 4;
            resetButton.Layout.Column = 3;
            resetButton.ButtonPushedFcn = @onResetOrigin;

            infoLabel = uilabel(grid);
            infoLabel.Text = sprintf('Image size: X = %d px, Y = %d px', Nx, Ny);
            infoLabel.FontAngle = 'italic';
            infoLabel.Layout.Row = 5;
            infoLabel.Layout.Column = [1 3];

            statusLabel = uilabel(grid);
            statusLabel.Text = '';
            statusLabel.FontColor = [0.65 0 0];
            statusLabel.Layout.Row = 6;
            statusLabel.Layout.Column = [1 3];

            buttonGrid = uigridlayout(grid);
            buttonGrid.RowHeight = {'1x'};
            buttonGrid.ColumnWidth = {'1x', '1x', '1x'};
            buttonGrid.Padding = [0 0 0 0];
            buttonGrid.Layout.Row = 7;
            buttonGrid.Layout.Column = [1 3];

            previewButton = uibutton(buttonGrid, 'push');
            previewButton.Text = 'Preview';
            previewButton.Layout.Row = 1;
            previewButton.Layout.Column = 1;
            previewButton.ButtonPushedFcn = @onPreview;

            applyButton = uibutton(buttonGrid, 'push');
            applyButton.Text = 'Apply';
            applyButton.Layout.Row = 1;
            applyButton.Layout.Column = 2;
            applyButton.ButtonPushedFcn = @onApply;

            cancelButton = uibutton(buttonGrid, 'push');
            cancelButton.Text = 'Cancel';
            cancelButton.Layout.Row = 1;
            cancelButton.Layout.Column = 3;
            cancelButton.ButtonPushedFcn = @onCancel;

            placeAppInsideCaller(app, dlg, 'center');
            dlg.Visible = 'on';

            uiwait(dlg);
            if isvalid(dlg)
                delete(dlg);
            end

            function [tf, parsedPixelRatio, parsedOriginXY] = validateDialogValues()
                tf = false;
                parsedPixelRatio = [];
                parsedOriginXY = [nan nan];

                rawPixelRatio = strtrim(char(string(pixelField.Value)));

                if isempty(rawPixelRatio) || strcmp(rawPixelRatio, '[]')
                    parsedPixelRatio = [];
                else
                    val = str2double(rawPixelRatio);

                    if ~isfinite(val) || val < 0
                        statusLabel.Text = 'Pixel ratio must be a non-negative scalar, 0, or [].';
                        return
                    end

                    if val == 0
                        parsedPixelRatio = [];
                    else
                        parsedPixelRatio = val;
                    end
                end

                parsedOriginXY = double([xField.Value, yField.Value]);

                if numel(parsedOriginXY) ~= 2 || any(~isfinite(parsedOriginXY))
                    statusLabel.Text = 'Origin coordinates must be finite numeric values.';
                    return
                end

                if parsedOriginXY(1) < 1 || parsedOriginXY(1) > Nx || ...
                        parsedOriginXY(2) < 1 || parsedOriginXY(2) > Ny
                    statusLabel.Text = sprintf( ...
                        'Origin must be inside image bounds: X [1 %d], Y [1 %d].', Nx, Ny);
                    return
                end

                statusLabel.Text = '';
                tf = true;
            end

            function onPreview(~, ~)
                [tf, parsedPixelRatio, parsedOriginXY] = validateDialogValues();
                if ~tf
                    return
                end

                app.previewViewCalibration(parsedPixelRatio, parsedOriginXY);
                app.setStatusMessage('View calibration preview updated.');
            end

            function onApply(~, ~)
                [tf, parsedPixelRatio, parsedOriginXY] = validateDialogValues();
                if ~tf
                    return
                end

                applyButton.Enable = 'off';
                drawnow limitrate

                try
                    app.applyViewCalibration(parsedPixelRatio, parsedOriginXY);
                    app.setStatusMessage('View calibration saved.');
                    uiresume(dlg);

                catch ME
                    applyButton.Enable = 'on';
                    statusLabel.Text = ME.message;
                    app.setStatusMessage(sprintf('View calibration failed: %s', ME.message));

                    fprintf(2, '\nView Calibration apply failed\n');
                    fprintf(2, '=============================\n');
                    fprintf(2, '%s\n', getReport(ME, 'extended', 'hyperlinks', 'on'));
                end
            end

            function onCancel(~, ~)
                app.DataParams = originalDataParams;
                app.ensureDataParamsViewFields();
                app.setImageAxisSpatialCalibration();
                app.updateImageStatusLabel();
                uiresume(dlg);
            end

            function onResetOrigin(~, ~)
                xField.Value = 1;
                yField.Value = 1;
                statusLabel.Text = '';

                try
                    app.previewViewCalibration(app.parsePixelRatioText(pixelField.Value), [1 1]);
                catch ME
                    statusLabel.Text = ME.message;
                end
            end

            function onInteractiveOrigin(~, ~)
                % If the current origin is default [1 1], start at image center.
                % Otherwise, start at the current origin.
                currentOrigin = double([xField.Value, yField.Value]);

                if isequal(round(currentOrigin), [1 1])
                    startXY = [round(Nx / 2), round(Ny / 2)];
                else
                    startXY = currentOrigin;
                end

                startXY(1) = min(max(round(startXY(1)), 1), Nx);
                startXY(2) = min(max(round(startXY(2)), 1), Ny);

                dlg.Visible = 'off';
                drawnow limitrate

                [newOriginXY, bConfirmed] = app.startInteractiveOriginSelection(startXY);

                if isvalid(dlg)
                    placeAppInsideCaller(app, dlg, 'center');
                    dlg.Visible = 'on';
                end

                if bConfirmed
                    xField.Value = newOriginXY(1);
                    yField.Value = newOriginXY(2);
                    statusLabel.Text = '';

                    try
                        app.previewViewCalibration(app.parsePixelRatioText(pixelField.Value), newOriginXY);
                    catch ME
                        statusLabel.Text = ME.message;
                    end
                else
                    statusLabel.Text = 'Interactive origin selection cancelled.';
                end
            end
        end

        function pixelRatio = parsePixelRatioText(app, txt) %#ok<INUSL>
            %PARSEPIXELRATIOTEXT Parse user-facing pixel-ratio text.
            %
            %   Returns [] for empty, [], or 0. Returns a positive scalar otherwise.

            txt = strtrim(char(string(txt)));

            if isempty(txt) || strcmp(txt, '[]')
                pixelRatio = [];
                return
            end

            pixelRatio = str2double(txt);

            if ~isfinite(pixelRatio) || pixelRatio < 0
                error('DataViewer:InvalidPixelRatio', ...
                    'Pixel ratio must be a non-negative scalar, 0, or [].');
            end

            if pixelRatio == 0
                pixelRatio = [];
            end

        end

        function [pixelRatio, originXY] = getCurrentViewCalibration(app)
            %GETCURRENTVIEWCALIBRATION Return current pixel ratio and origin.

            app.ensureDataParamsViewFields();

            pixelRatio = app.DataParams.view.pixelSize_px_per_mm;
            if isempty(pixelRatio)
                pixelRatio = [];
            else
                pixelRatio = double(pixelRatio(1));
                if ~isfinite(pixelRatio) || pixelRatio <= 0
                    pixelRatio = [];
                end
            end

            originXY = double(app.DataParams.view.origin_xy_px(:).');
            if numel(originXY) ~= 2 || any(~isfinite(originXY))
                originXY = [1 1];
            end

            if app.hasData()
                sz = app.getDataSize();
                originXY(1) = min(max(originXY(1), 1), sz(2));
                originXY(2) = min(max(originXY(2), 1), sz(1));
            end

        end

        function ensureDataParamsViewFields(app)
            %ENSUREDATAPARAMSVIEWFIELDS Ensure DataParams.view has calibration fields.

            if isempty(app.DataParams) || ~isstruct(app.DataParams)
                app.DataParams = struct();
            end

            if ~isfield(app.DataParams, 'view') || ~isstruct(app.DataParams.view)
                app.DataParams.view = struct();
            end

            if ~isfield(app.DataParams.view, 'pixelSize_px_per_mm') || ...
                    isempty(app.DataParams.view.pixelSize_px_per_mm)
                app.DataParams.view.pixelSize_px_per_mm = [];
            else
                pxPerMm = double(app.DataParams.view.pixelSize_px_per_mm);
                pxPerMm = pxPerMm(:).';

                if isempty(pxPerMm) || any(~isfinite(pxPerMm)) || pxPerMm(1) <= 0
                    app.DataParams.view.pixelSize_px_per_mm = [];
                else
                    app.DataParams.view.pixelSize_px_per_mm = pxPerMm(1);
                end
            end

            if ~isfield(app.DataParams.view, 'origin_xy_px') || ...
                    isempty(app.DataParams.view.origin_xy_px)
                app.DataParams.view.origin_xy_px = [1 1];
            end

            originXY = double(app.DataParams.view.origin_xy_px(:).');
            if numel(originXY) ~= 2 || any(~isfinite(originXY))
                originXY = [1 1];
            end

            if app.hasData()
                sz = app.getDataSize();
                originXY(1) = min(max(originXY(1), 1), sz(2));
                originXY(2) = min(max(originXY(2), 1), sz(1));
            end

            app.DataParams.view.origin_xy_px = originXY;

            if app.hasData()
                sz = app.getDataSize();
                app.DataParams.view.imageSizeYX = double(sz(1:2));
            end

        end

        function previewViewCalibration(app, pixelRatio, originXY)
            %PREVIEWVIEWCALIBRATION Apply calibration to viewer without saving.

            app.setViewCalibrationInDataParams(pixelRatio, originXY);
            app.setImageAxisSpatialCalibration();
            app.updateImageStatusLabel();

        end

        function applyViewCalibration(app, pixelRatio, originXY)
            %APPLYVIEWCALIBRATION Save view calibration and refresh image axes.
            %
            %   pixelRatio:
            %       [] or positive scalar px/mm.
            %
            %   originXY:
            %       [x y] origin coordinates in image pixels.

            if ~app.hasData()
                error('DataViewer:NoDataLoaded', ...
                    'Load image data before applying view calibration.');
            end

            pixelRatio = app.normalizePixelRatioValue(pixelRatio);
            originXY = app.validateOriginXY(originXY);

            folderPath = app.getCurrentDataFolder();

            % Ensure a valid folder-level DataParams.mat exists before using
            % updateDataParam, which expects intermediate fields to exist.
            app.DataParams = app.ensureDataParamsFileForCurrentFolder();
            app.ensureDataParamsViewFields();

            sz = app.getDataSize();
            imageSizeYX = double(sz(1:2));

            % Persist using the existing folder-centric update function.
            % Validate only after the final related update to avoid transient
            % intermediate validation issues.
            DataParams = updateDataParam(folderPath, ...
                'view.imageSizeYX', imageSizeYX, ...
                'validateAfterSet', false);

            DataParams = updateDataParam(folderPath, ...
                'view.origin_xy_px', originXY, ...
                'validateAfterSet', false);

            DataParams = updateDataParam(folderPath, ...
                'view.pixelSize_px_per_mm', pixelRatio, ...
                'validateAfterSet', true);

            % Critical: update runtime state from the persisted/validated struct.
            app.DataParams = DataParams;
            app.ensureDataParamsViewFields();

            app.setImageAxisSpatialCalibration();
            app.updateImageStatusLabel();
            app.refreshImageFrame();
        end

        function originXY = validateOriginXY(app, originXY)
            %VALIDATEORIGINXY Validate origin [x y] against current image size.

            if ~app.hasData()
                error('DataViewer:NoDataLoaded', ...
                    'Load image data before setting the origin.');
            end

            sz = app.getDataSize();
            Ny = sz(1);
            Nx = sz(2);

            originXY = double(originXY(:).');

            if numel(originXY) ~= 2 || any(~isfinite(originXY))
                error('DataViewer:InvalidOrigin', ...
                    'Origin must be a finite [x y] coordinate.');
            end

            if originXY(1) < 1 || originXY(1) > Nx || ...
                    originXY(2) < 1 || originXY(2) > Ny
                error('DataViewer:OriginOutOfBounds', ...
                    'Origin must be inside image bounds: X [1 %d], Y [1 %d].', Nx, Ny);
            end

            originXY = round(originXY);

        end

        function pixelRatio = normalizePixelRatioValue(app, pixelRatio) %#ok<INUSD>
            %NORMALIZEPIXELRATIOVALUE Normalize user-facing pixel ratio.
            %
            %   User-facing 0 or [] means no pixel ratio. Backend storage uses [].

            if isempty(pixelRatio)
                pixelRatio = [];
                return
            end

            if ischar(pixelRatio) || isstring(pixelRatio)
                txt = strtrim(char(string(pixelRatio)));

                if isempty(txt) || strcmp(txt, '[]') || strcmp(txt, '0')
                    pixelRatio = [];
                    return
                end

                pixelRatio = str2double(txt);
            end

            pixelRatio = double(pixelRatio);

            if ~isscalar(pixelRatio) || ~isfinite(pixelRatio) || pixelRatio < 0
                error('DataViewer:InvalidPixelRatio', ...
                    'Pixel ratio must be [], 0, or a positive scalar px/mm.');
            end

            if pixelRatio == 0
                pixelRatio = [];
            end

        end

        function setViewCalibrationInDataParams(app, pixelRatio, originXY)
            %SETVIEWCALIBRATIONINDATAPARAMS Update app.DataParams.view calibration fields.

            app.ensureDataParamsViewFields();

            if isempty(pixelRatio) || pixelRatio == 0
                app.DataParams.view.pixelSize_px_per_mm = [];
            else
                validateattributes(pixelRatio, {'numeric'}, ...
                    {'scalar', 'finite', 'positive'}, ...
                    'setViewCalibrationInDataParams', 'pixelRatio');
                app.DataParams.view.pixelSize_px_per_mm = double(pixelRatio);
            end

            originXY = double(originXY(:).');

            if numel(originXY) ~= 2 || any(~isfinite(originXY))
                error('DataViewer:InvalidOrigin', ...
                    'Origin must be a finite [X Y] coordinate pair.');
            end

            if app.hasData()
                sz = app.getDataSize();
                if originXY(1) < 1 || originXY(1) > sz(2) || ...
                        originXY(2) < 1 || originXY(2) > sz(1)
                    error('DataViewer:InvalidOrigin', ...
                        'Origin must be inside image bounds: X [1 %d], Y [1 %d].', ...
                        sz(2), sz(1));
                end
            end

            app.DataParams.view.origin_xy_px = originXY;

        end

        function saveDataParamsForCurrentFile(app)
            %SAVEDATAPARAMSFORCURRENTFILE Validate and save current app.DataParams.

            folderPath = app.getCurrentDataFolder();

            if isempty(app.DataParams) || ~isstruct(app.DataParams)
                app.DataParams = app.ensureDataParamsFileForCurrentFolder();
            end

            app.ensureDataParamsViewFields();
            app.ensureDataParamsMaskFields();
            app.ensureDataParamsFolderFields();

            DataParams = app.DataParams;
            saveDataParams(folderPath, DataParams);

            % Reload saved version so app.DataParams includes lastModified update.
            app.DataParams = loadDataParams(folderPath);
            app.ensureDataParamsViewFields();
            app.ensureDataParamsMaskFields();
            app.ensureDataParamsFolderFields();

        end

        function updateImageStatusLabel(app)
            %UPDATEIMAGESTATUSLABEL Show view calibration summary in ImageStatusLabel.

            if isempty(app.ImageStatusLabel) || ~isvalid(app.ImageStatusLabel)
                return
            end

            if ~app.hasData()
                app.ImageStatusLabel.Text = 'pixel ratio (px/mm): none | origin (x,y): [ ]';
                return
            end

            [pixelRatio, originXY] = app.getCurrentViewCalibration();

            if isempty(pixelRatio)
                pixelText = 'none';
            else
                pixelText = sprintf('%.6g', pixelRatio);
            end

            app.ImageStatusLabel.Text = sprintf( ...
                'pixel ratio (px/mm): %s | origin (x,y): [%.3g %.3g]', ...
                pixelText, originXY(1), originXY(2));

        end

        function [originXY, bConfirmed] = startInteractiveOriginSelection(app, startOriginXY)
            %STARTINTERACTIVEORIGINSELECTION Select image origin with an interactive point.
            %
            %   User confirms by double-clicking the point or pressing Enter. Escape
            %   cancels selection. The point is green when idle and red while moving.

            originXY = double(startOriginXY(:).');
            bConfirmed = false;

            if ~app.hasData()
                return
            end

            sz = app.getDataSize();
            Ny = sz(1);
            Nx = sz(2);

            if numel(originXY) ~= 2 || any(~isfinite(originXY))
                originXY = [round(Nx / 2), round(Ny / 2)];
            end

            originXY(1) = min(max(round(originXY(1)), 1), Nx);
            originXY(2) = min(max(round(originXY(2)), 1), Ny);

            previousMode = app.InteractionMode;
            previousKeyFcn = app.UIFigure.WindowKeyPressFcn;

            hPoint = [];
            movingListener = [];
            movedListener = [];
            clickedListener = [];

            bDone = false;

            try
                app.setInteractionMode('settingViewCalibration');
                app.setStatusMessage('Move origin point. Double-click or press Enter to confirm; press Escape to cancel.');

                drawnow limitrate
                figure(app.UIFigure);

                hold(app.ImageAxes, 'on');

                hPoint = drawpoint(app.ImageAxes, ...
                    'Position', originXY, ...
                    'Color', [0 0.8 0], ...
                    'Label', sprintf('X %.0f, Y %.0f', originXY(1), originXY(2)), ...
                    'InteractionsAllowed', 'translate');

                movingListener = addlistener(hPoint, 'MovingROI', @onMovingPoint);
                movedListener = addlistener(hPoint, 'ROIMoved', @onMovedPoint);
                clickedListener = addlistener(hPoint, 'ROIClicked', @onClickedPoint);

                app.UIFigure.WindowKeyPressFcn = @onKeyPress;

                while ~bDone && isvalid(app.UIFigure) && ~isempty(hPoint) && isvalid(hPoint)
                    drawnow
                    pause(0.02)
                end

                if bConfirmed && ~isempty(hPoint) && isvalid(hPoint)
                    pos = double(hPoint.Position(:).');
                    originXY = [round(pos(1)), round(pos(2))];
                    originXY(1) = min(max(originXY(1), 1), Nx);
                    originXY(2) = min(max(originXY(2), 1), Ny);
                end

            catch ME
                cleanupInteractiveOrigin();
                rethrow(ME)
            end

            cleanupInteractiveOrigin();

            function onMovingPoint(src, evt)
                posNow = clampPointPosition(evt.CurrentPosition);

                try
                    src.Color = [0.9 0 0];
                    src.Label = sprintf('X %.0f, Y %.0f', posNow(1), posNow(2));
                catch
                end
            end

            function onMovedPoint(src, evt)
                posNow = clampPointPosition(evt.CurrentPosition);

                try
                    src.Position = posNow;
                    src.Color = [0 0.8 0];
                    src.Label = sprintf('X %.0f, Y %.0f', posNow(1), posNow(2));
                catch
                end
            end

            function onClickedPoint(~, evt)
                if isprop(evt, 'SelectionType') && strcmpi(evt.SelectionType, 'double')
                    bConfirmed = true;
                    bDone = true;
                end
            end

            function onKeyPress(~, evt)
                switch lower(evt.Key)
                    case {'return', 'enter'}
                        bConfirmed = true;
                        bDone = true;

                    case 'escape'
                        bConfirmed = false;
                        bDone = true;
                end
            end

            function pos = clampPointPosition(pos)
                pos = double(pos(:).');
                pos(1) = min(max(pos(1), 1), Nx);
                pos(2) = min(max(pos(2), 1), Ny);
            end

            function cleanupInteractiveOrigin()
                try
                    if ~isempty(movingListener) && isvalid(movingListener)
                        delete(movingListener);
                    end
                catch
                end

                try
                    if ~isempty(movedListener) && isvalid(movedListener)
                        delete(movedListener);
                    end
                catch
                end

                try
                    if ~isempty(clickedListener) && isvalid(clickedListener)
                        delete(clickedListener);
                    end
                catch
                end

                try
                    if ~isempty(hPoint) && isvalid(hPoint)
                        delete(hPoint);
                    end
                catch
                end

                try
                    if isvalid(app.UIFigure)
                        app.UIFigure.WindowKeyPressFcn = previousKeyFcn;

                    end
                catch
                end

                try
                    app.setInteractionMode(previousMode);
                catch
                    app.InteractionMode = previousMode;
                end
            end

        end

        function folderPath = getCurrentDataFolder(app)
            %GETCURRENTDATAFOLDER Return folder for the currently loaded data file.

            if isempty(app.CurrentFile)
                error('DataViewer:NoCurrentFile', ...
                    'No current file is loaded.');
            end

            folderPath = fileparts(app.CurrentFile);

            if ~isfolder(folderPath)
                error('DataViewer:InvalidCurrentFolder', ...
                    'Current data folder does not exist: %s', folderPath);
            end

        end

        function DataParams = ensureDataParamsFileForCurrentFolder(app)
            %ENSUREDATAPARAMSFILEFORCURRENTFOLDER Ensure DataParams.mat exists.

            folderPath = app.getCurrentDataFolder();
            dataParamsPath = fullfile(folderPath, 'DataParams.mat');

            if ~isfile(dataParamsPath)
                DataParams = createDataParams(folderPath);
            else
                DataParams = loadDataParams(folderPath);
            end

            app.DataParams = DataParams;
            app.ensureDataParamsViewFields();
            app.ensureDataParamsMaskFields();
            app.ensureDataParamsFolderFields();
            DataParams = app.DataParams;

        end



        % =========================================================================
        % ROI table setup and callbacks
        % =========================================================================

        function setupROITable(app)
            %SETUPROITABLE Configure ROI table columns and callbacks.
            %
            %   The Select checkbox is the edit-mode controller:
            %       - one selected ROI starts single-ROI edit mode
            %       - two or more selected ROIs start/update group edit mode
            %       - no selected ROIs exits edit mode

            if isempty(app.UITable) || ~isvalid(app.UITable)
                return
            end

            columnNames = app.getROITableColumnNames();

            app.UITable.ColumnName = columnNames;
            app.UITable.RowName = {};

            app.UITable.ColumnEditable = [ ...
                true, ...   % Visible
                true, ...   % Select
                true, ...   % Name
                false, ...  % Type
                false, ...  % Color
                false(1, numel(columnNames) - 5)];

            app.UITable.ColumnFormat = { ...
                'logical', ...
                'logical', ...
                'char', ...
                'char', ...
                'char', ...
                'numeric', ...
                'numeric', ...
                'numeric', ...
                'numeric', ...
                'numeric', ...
                'numeric', ...
                'numeric', ...
                'numeric', ...
                'numeric', ...
                'numeric'};

            try
                if isprop(app.UITable, 'Multiselect')
                    app.UITable.Multiselect = 'on';
                end
            catch
            end

            app.UITable.CellEditCallback = @(src, event) app.onROITableCellEdit(event);
            app.UITable.CellSelectionCallback = @(src, event) app.onROITableCellSelection(event);
            app.createDeleteROIButtonContextMenu();
            app.refreshROITable();

        end

        function columnNames = getROITableColumnNames(app)
            %GETROITABLECOLUMNNAMES Return ROI table columns with coordinate units.

            if app.roiTableUsesMM()
                centroidX = 'Centroid X (mm)';
                centroidY = 'Centroid Y (mm)';
                distance  = 'Distance (mm)';
                areaText  = 'Area (mm²)';
            else
                centroidX = 'Centroid X (px)';
                centroidY = 'Centroid Y (px)';
                distance  = 'Distance (px)';
                areaText  = 'Area (px²)';
            end

            columnNames = { ...
                'Visible', ...
                'Select', ...
                'Name', ...
                'Type', ...
                'Color', ...
                centroidX, ...
                centroidY, ...
                distance, ...
                areaText, ...
                'N pixels', ...
                'Mean', ...
                'Std', ...
                'Median', ...
                'Min', ...
                'Max'};

        end

        function tf = roiTableUsesMM(app)
            %ROITABLEUSESMM True when ROI coordinate columns should display mm units.

            tf = false;

            if isempty(app.DataParams) || ~isfield(app.DataParams, 'view')
                return
            end

            viewInfo = app.DataParams.view;

            if ~isfield(viewInfo, 'pixelSize_px_per_mm') || isempty(viewInfo.pixelSize_px_per_mm)
                return
            end

            pxPerMm = double(viewInfo.pixelSize_px_per_mm);

            tf = isscalar(pxPerMm) && isfinite(pxPerMm) && pxPerMm > 0;

        end

        function refreshROITable(app)
            %REFRESHROITABLE Rebuild ROI table from app.ROIList.
            %
            %   The table Data is stored as a cell array instead of a MATLAB table to
            %   avoid fragile generated variable names from display headers containing
            %   spaces, parentheses, and special unit symbols.

            if isempty(app.UITable) || ~isvalid(app.UITable)
                return
            end

            columnNames = app.getROITableColumnNames();
            app.UITable.ColumnName = columnNames;

            nCol = numel(columnNames);

            if isempty(app.ROIList)
                app.UITable.Data = cell(0, nCol);
                app.UITable.UserData = struct('RowROIIDs', []);
                app.refreshROITableColorStyles();
                return
            end

            nROI = numel(app.ROIList);
            tableData = cell(nROI, nCol);
            rowROIIDs = nan(nROI, 1);

            useMM = app.roiTableUsesMM();

            for iROI = 1:nROI
                ROI = app.ROIList(iROI);

                if isfield(ROI, 'ID') && ~isempty(ROI.ID) && isfinite(ROI.ID)
                    rowROIIDs(iROI) = ROI.ID;
                else
                    rowROIIDs(iROI) = iROI;
                end

                visibleValue = true;
                if isfield(ROI, 'runtime') && isfield(ROI.runtime, 'visible')
                    visibleValue = logical(ROI.runtime.visible);
                end

                selectValue = false;
                if isfield(ROI, 'runtime') && isfield(ROI.runtime, 'selected')
                    selectValue = logical(ROI.runtime.selected);
                end

                roiName = '';
                if isfield(ROI, 'name') && ~isempty(ROI.name)
                    roiName = char(string(ROI.name));
                end

                roiType = '';
                if isfield(ROI, 'type') && ~isempty(ROI.type)
                    roiType = char(string(ROI.type));
                end

                centroidX = NaN;
                centroidY = NaN;
                distanceValue = NaN;
                areaValue = NaN;
                nPixels = NaN;
                meanValue = NaN;
                stdValue = NaN;
                medianValue = NaN;
                minValue = NaN;
                maxValue = NaN;

                if isfield(ROI, 'stats') && isstruct(ROI.stats)
                    stats = ROI.stats;

                    if useMM
                        if isfield(stats, 'centroidXY_mm') && numel(stats.centroidXY_mm) == 2
                            centroidX = stats.centroidXY_mm(1);
                            centroidY = stats.centroidXY_mm(2);
                        end
                        if isfield(stats, 'distanceFromOrigin_mm') && ~isempty(stats.distanceFromOrigin_mm)
                            distanceValue = stats.distanceFromOrigin_mm;
                        end
                        if isfield(stats, 'areaMM2') && ~isempty(stats.areaMM2)
                            areaValue = stats.areaMM2;
                        end
                    else
                        if isfield(stats, 'centroidXY_px') && numel(stats.centroidXY_px) == 2
                            centroidX = stats.centroidXY_px(1);
                            centroidY = stats.centroidXY_px(2);
                        end
                        if isfield(stats, 'distanceFromOrigin_px') && ~isempty(stats.distanceFromOrigin_px)
                            distanceValue = stats.distanceFromOrigin_px;
                        end
                        if isfield(stats, 'areaPx2') && ~isempty(stats.areaPx2)
                            areaValue = stats.areaPx2;
                        end
                    end

                    if isfield(stats, 'NPixels') && ~isempty(stats.NPixels)
                        nPixels = stats.NPixels;
                    end
                    if isfield(stats, 'spatialMean') && ~isempty(stats.spatialMean)
                        meanValue = stats.spatialMean;
                    end
                    if isfield(stats, 'spatialStd') && ~isempty(stats.spatialStd)
                        stdValue = stats.spatialStd;
                    end
                    if isfield(stats, 'spatialMedian') && ~isempty(stats.spatialMedian)
                        medianValue = stats.spatialMedian;
                    end
                    if isfield(stats, 'spatialMin') && ~isempty(stats.spatialMin)
                        minValue = stats.spatialMin;
                    end
                    if isfield(stats, 'spatialMax') && ~isempty(stats.spatialMax)
                        maxValue = stats.spatialMax;
                    end
                end

                tableData(iROI, :) = { ...
                    visibleValue, ...
                    selectValue, ...
                    roiName, ...
                    roiType, ...
                    '', ...
                    centroidX, ...
                    centroidY, ...
                    distanceValue, ...
                    areaValue, ...
                    nPixels, ...
                    meanValue, ...
                    stdValue, ...
                    medianValue, ...
                    minValue, ...
                    maxValue};
            end

            app.UITable.Data = tableData;
            app.UITable.UserData = struct('RowROIIDs', rowROIIDs);

            app.refreshROITableColorStyles();

        end

        function refreshROITableColorStyles(app)
            %REFRESHROITABLECOLORSTYLES Color the ROI table Color column.

            if isempty(app.UITable) || ~isvalid(app.UITable)
                return
            end

            try
                removeStyle(app.UITable);
            catch
            end

            if isempty(app.ROIList)
                return
            end

            colorCol = app.getROITableColumnIndex('Color');
            if isempty(colorCol)
                return
            end

            nROI = numel(app.ROIList);

            for iROI = 1:nROI
                thisColor = [1 1 1];

                if isfield(app.ROIList(iROI), 'color') && numel(app.ROIList(iROI).color) == 3
                    thisColor = double(app.ROIList(iROI).color(:).');
                    thisColor = min(max(thisColor, 0), 1);
                end

                try
                    styleObj = uistyle( ...
                        'BackgroundColor', thisColor, ...
                        'FontColor', thisColor);
                    addStyle(app.UITable, styleObj, 'cell', [iROI, colorCol]);
                catch
                end
            end

        end

        function colIdx = getROITableColumnIndex(app, columnName)
            %GETROITABLECOLUMNINDEX Return table column index by visible column name.

            colIdx = [];

            if isempty(app.UITable) || ~isvalid(app.UITable)
                return
            end

            names = string(app.UITable.ColumnName);
            idx = find(strcmpi(names, string(columnName)), 1, 'first');

            if ~isempty(idx)
                colIdx = idx;
            end

        end

        function onROITableCellSelection(app, event)
            %ONROITABLECELLSELECTION Handle non-edit table clicks.
            %
            %   ROI edit selection is controlled by the explicit Select checkbox.
            %   Plain row/cell selection does not enter edit mode and does not update
            %   app.SelectedROIID, which avoids stale checked Select cells after table
            %   refreshes.

            if isempty(event.Indices)
                return
            end

            rowIdx = event.Indices(1);
            colIdx = event.Indices(2);

            roiID = app.getROIIDFromTableRow(rowIdx);
            if isempty(roiID)
                return
            end

            colorCol = app.getROITableColumnIndex('Color');
            if ~isempty(colorCol) && colIdx == colorCol
                app.pickROIColorFromTableRow(rowIdx);
            end

        end

        function onROITableCellEdit(app, event)
            %ONROITABLECELLEDIT Handle editable ROI table fields.

            if isempty(event.Indices)
                return
            end

            rowIdx = event.Indices(1);
            colIdx = event.Indices(2);

            roiID = app.getROIIDFromTableRow(rowIdx);
            if isempty(roiID)
                app.refreshROITable();
                return
            end

            roiIdx = app.getROIIndexByID(roiID);
            if isempty(roiIdx)
                app.refreshROITable();
                return
            end

            colNames = string(app.UITable.ColumnName);
            colName = char(colNames(colIdx));

            switch colName
                case 'Visible'
                    app.setROIVisibilityByIndex(roiIdx, logical(event.NewData));

                case 'Select'
                    % Preserve checkbox state and read all selected rows from the table. Do
                    % not rebuild the table here; rebuilding can interrupt sequential
                    % multi-ROI checkbox edits in some MATLAB releases.
                    try
                        tableData = app.UITable.Data;
                        if iscell(tableData) && rowIdx <= size(tableData, 1) && colIdx <= size(tableData, 2)
                            tableData{rowIdx, colIdx} = logical(event.NewData);
                            app.UITable.Data = tableData;
                        end
                    catch
                    end

                    app.updateROISelectionOrder(roiID, logical(event.NewData));

                    app.syncROISelectionStateFromTable();
                    app.handleROISelectionStateChanged();
                    return

                case 'Name'
                    requestedName = strtrim(char(string(event.NewData)));

                    if isempty(requestedName)
                        app.setStatusMessage('ROI name cannot be empty.');
                        app.refreshROITable();
                        return
                    end

                    uniqueName = app.makeUniqueROIName(requestedName, roiID);

                    if ~strcmp(uniqueName, requestedName)
                        app.setStatusMessage(sprintf( ...
                            'ROI name "%s" already exists. Renamed to "%s".', ...
                            requestedName, uniqueName));
                    end

                    app.ROIList(roiIdx).name = uniqueName;
                    app.ROIList(roiIdx).modifiedOn = datetime('now');

                    app.updateROITraceDisplayName(roiIdx);

                otherwise
                    app.refreshROITable();
                    return
            end

            app.refreshROITable();
            app.updateGUIEnabledState();

        end

        function syncROISelectionStateFromTable(app)
            %SYNCROISELECTIONSTATEFROMTABLE Copy Select checkbox state into ROIList.

            if isempty(app.ROIList) || isempty(app.UITable) || ~isvalid(app.UITable)
                return
            end

            selectCol = app.getROITableColumnIndex('Select');
            if isempty(selectCol)
                return
            end

            if isempty(app.UITable.UserData) || ~isstruct(app.UITable.UserData) || ...
                    ~isfield(app.UITable.UserData, 'RowROIIDs')
                return
            end

            tableData = app.UITable.Data;
            rowROIIDs = app.UITable.UserData.RowROIIDs;

            if isempty(tableData) || ~iscell(tableData)
                return
            end

            nRows = min(size(tableData, 1), numel(rowROIIDs));

            for iRow = 1:nRows
                roiID = rowROIIDs(iRow);
                roiIdx = app.getROIIndexByID(roiID);

                if isempty(roiIdx)
                    continue
                end

                selectedValue = false;

                try
                    selectedValue = logical(tableData{iRow, selectCol});
                catch
                    selectedValue = false;
                end

                app.setROISelectedStateByIndex(roiIdx, selectedValue);
            end

        end

        function handleROISelectionStateChanged(app)
            %HANDLEROISELECTIONSTATECHANGED Route Select-checkbox state to edit mode.
            %
            %   Expected behavior:
            %       1 selected ROI  -> single-ROI edit mode
            %       2+ selected ROI -> group-ROI edit mode
            %       0 selected ROI  -> exit edit mode

            selectedIDs = app.getSelectedROIIDList();

            if isempty(selectedIDs)
                app.stopActiveROIEditForSelectionChange();
                app.setSelectedROI(NaN);
                app.updateGUIEnabledState();
                app.setStatusMessage('ROI edit mode exited.');
                return
            end

            if numel(selectedIDs) == 1
                roiID = selectedIDs(1);
                roiIdx = app.getROIIndexByID(roiID);

                if isempty(roiIdx)
                    app.updateGUIEnabledState();
                    return
                end

                activeSingleID = app.getActiveSingleEditROIID();
                activeGroupIDs = app.getActiveGroupEditROIIDList();

                if isequal(activeSingleID, roiID) && isempty(activeGroupIDs)
                    app.setSelectedROI(roiID);
                    app.updateGUIEnabledState();
                    return
                end

                app.stopActiveROIEditForSelectionChange();
                app.setSelectedROI(roiID);
                app.beginROIEditingByIndex(roiIdx);
                return
            end

            selectedIDs = unique(double(selectedIDs(:).'), 'stable');
            activeGroupIDs = app.getActiveGroupEditROIIDList();

            if isequal(activeGroupIDs, selectedIDs)
                app.updateGUIEnabledState();
                return
            end

            app.stopActiveROIEditForSelectionChange();
            app.setSelectedROI(selectedIDs(end));
            app.startGroupROIEditByIDs(selectedIDs);

        end

        function roiID = getActiveSingleEditROIID(app)
            %GETACTIVESINGLEEDITROIID Return ROI ID currently edited as a single ROI.

            roiID = [];

            if isempty(app.ROIList)
                return
            end

            for iROI = 1:numel(app.ROIList)
                if isfield(app.ROIList(iROI), 'runtime') && ...
                        isfield(app.ROIList(iROI).runtime, 'editHandle') && ...
                        app.isUsableGraphicsHandle(app.ROIList(iROI).runtime.editHandle)

                    if isfield(app.ROIList(iROI), 'ID') && isfinite(app.ROIList(iROI).ID)
                        roiID = app.ROIList(iROI).ID;
                    end
                    return
                end
            end

        end

        function roiIDList = getActiveGroupEditROIIDList(app)
            %GETACTIVEGROUPEDITROIIDLIST Return ROI IDs currently edited as a group.

            roiIDList = [];

            if isempty(fieldnames(app.GroupROIEditState)) || ...
                    ~isfield(app.GroupROIEditState, 'roiIDList')
                return
            end

            roiIDList = double(app.GroupROIEditState.roiIDList(:).');
            roiIDList = roiIDList(isfinite(roiIDList));

        end

        function stopActiveROIEditForSelectionChange(app)
            %STOPACTIVEROIEDITFORSELECTIONCHANGE Stop edit graphics without clearing Select.
            %
            %   This is used when the table Select column changes. It exits the current
            %   single/group edit visualization, restores static overlays, and leaves
            %   ROI.runtime.selected untouched so the new selection state can drive the next
            %   edit mode.

            % Stop active group edit runtime.
            groupRoiIdxList = [];
            if isfield(app.GroupROIEditState, 'roiIdxList')
                groupRoiIdxList = app.GroupROIEditState.roiIdxList;
            end

            if ~isempty(groupRoiIdxList)
                app.deleteGroupEditRuntimeGraphics();

                for k = 1:numel(groupRoiIdxList)
                    roiIdx = groupRoiIdxList(k);
                    if roiIdx >= 1 && roiIdx <= numel(app.ROIList)
                        if ~isfield(app.ROIList(roiIdx).runtime, 'ROIHandle') || ...
                                ~app.isUsableGraphicsHandle(app.ROIList(roiIdx).runtime.ROIHandle)
                            app.ROIList(roiIdx).runtime.ROIHandle = ...
                                app.createStaticROIOverlayFromROI(app.ROIList(roiIdx));
                        end
                    end
                end
            end

            % Stop active single-ROI edit runtime.
            for iROI = 1:numel(app.ROIList)
                if isfield(app.ROIList(iROI), 'runtime') && ...
                        isfield(app.ROIList(iROI).runtime, 'editHandle') && ...
                        app.isUsableGraphicsHandle(app.ROIList(iROI).runtime.editHandle)

                    app.cleanupROIEditRuntimeByIndex(iROI);

                    if ~isfield(app.ROIList(iROI).runtime, 'ROIHandle') || ...
                            ~app.isUsableGraphicsHandle(app.ROIList(iROI).runtime.ROIHandle)
                        app.ROIList(iROI).runtime.ROIHandle = ...
                            app.createStaticROIOverlayFromROI(app.ROIList(iROI));
                    end
                end
            end

            try
                app.setInteractionMode('idle');
            catch
                app.InteractionMode = 'idle';
                app.updateGUIEnabledState();
            end

        end

        function clearROISelectionState(app, roiIDList)
            %CLEARROISELECTIONSTATE Clear temporary ROI table-selection state.

            if isempty(app.ROIList)
                app.SelectedROIID = NaN;
                return
            end

            if nargin < 2 || isempty(roiIDList)
                roiIDList = [app.ROIList.ID];
            end

            roiIDList = double(roiIDList(:).');

            for iROI = 1:numel(app.ROIList)
                if isfield(app.ROIList(iROI), 'ID') && any(app.ROIList(iROI).ID == roiIDList)
                    app.setROISelectedStateByIndex(iROI, false);
                end
            end

            selectedIDs = app.getSelectedROIIDList();
            if isempty(selectedIDs)
                app.SelectedROIID = NaN;
            elseif ~any(selectedIDs == app.SelectedROIID)
                app.SelectedROIID = selectedIDs(end);
            end

            if ~isempty(app.UITable) && isvalid(app.UITable)
                selectCol = app.getROITableColumnIndex('Select');

                if ~isempty(selectCol) && iscell(app.UITable.Data) && ...
                        ~isempty(app.UITable.UserData) && isstruct(app.UITable.UserData) && ...
                        isfield(app.UITable.UserData, 'RowROIIDs')

                    tableData = app.UITable.Data;
                    rowROIIDs = app.UITable.UserData.RowROIIDs;
                    nRows = min(size(tableData, 1), numel(rowROIIDs));

                    for iRow = 1:nRows
                        if any(rowROIIDs(iRow) == roiIDList)
                            tableData{iRow, selectCol} = false;
                        end
                    end

                    app.UITable.Data = tableData;
                end
            end

            if nargin < 2 || isempty(roiIDList)
                app.ROISelectionOrder = [];
            else
                app.ROISelectionOrder = app.ROISelectionOrder(~ismember(app.ROISelectionOrder, roiIDList));
            end
        end

        function roiID = getROIIDFromTableRow(app, rowIdx)
            %GETROIIDFROMTABLEROW Return ROI runtime ID mapped to one table row.

            roiID = [];

            if isempty(app.UITable) || ~isvalid(app.UITable)
                return
            end

            if isempty(app.UITable.UserData) || ~isstruct(app.UITable.UserData) || ...
                    ~isfield(app.UITable.UserData, 'RowROIIDs')
                return
            end

            rowROIIDs = app.UITable.UserData.RowROIIDs;

            if rowIdx < 1 || rowIdx > numel(rowROIIDs)
                return
            end

            roiID = rowROIIDs(rowIdx);

            if ~isfinite(roiID)
                roiID = [];
            end

        end

        function idx = getROIIndexByID(app, roiID)
            %GETROIINDEXBYID Return ROIList index for a runtime ROI ID.

            idx = [];

            if isempty(app.ROIList)
                return
            end

            roiIDs = [app.ROIList.ID];
            idx = find(roiIDs == roiID, 1, 'first');

        end



        function setSelectedROI(app, roiID)
            %SETSELECTEDROI Store selected ROI ID.

            if isempty(roiID) || ~isfinite(roiID)
                app.SelectedROIID = NaN;
            else
                app.SelectedROIID = roiID;
            end

        end

        function pickROIColorFromTableRow(app, rowIdx)
            %PICKROICOLORFROMTABLEROW Open color picker for clicked ROI color cell.

            roiID = app.getROIIDFromTableRow(rowIdx);
            if isempty(roiID)
                return
            end

            roiIdx = app.getROIIndexByID(roiID);
            if isempty(roiIdx)
                return
            end

            oldColor = [1 0 0];

            if isfield(app.ROIList(roiIdx), 'color') && numel(app.ROIList(roiIdx).color) == 3
                oldColor = double(app.ROIList(roiIdx).color(:).');
            end

            newColor = uisetcolor(oldColor, sprintf('Select color for %s', app.ROIList(roiIdx).name));

            if ~(isnumeric(newColor) && isequal(size(newColor), [1 3]) && all(isfinite(newColor)))
                return
            end

            newColor = min(max(double(newColor), 0), 1);

            app.setROIColorByIndex(roiIdx, newColor);
            app.refreshROITable();
            app.updateGUIEnabledState();

        end

        function setROIVisibilityByIndex(app, roiIdx, tf)
            %SETROIVISIBILITYBYINDEX Apply ROI visibility to runtime graphics.
            %
            %   Visibility is a graphics-only operation. ROI traces are not recomputed
            %   here because all ROI traces are computed by refreshROITraces regardless
            %   of visibility state.

            if isempty(roiIdx) || roiIdx < 1 || roiIdx > numel(app.ROIList)
                return
            end

            tf = logical(tf);
            app.ROIList(roiIdx).runtime.visible = tf;

            if tf
                visText = 'on';
            else
                visText = 'off';
            end

            runtimeFields = {'ROIHandle', 'traceLine', 'tracePatch'};

            for iField = 1:numel(runtimeFields)
                fieldName = runtimeFields{iField};

                if ~isfield(app.ROIList(roiIdx).runtime, fieldName)
                    continue
                end

                h = app.ROIList(roiIdx).runtime.(fieldName);

                if app.isUsableGraphicsHandle(h) && isprop(h, 'Visible')
                    try
                        h.Visible = visText;
                    catch
                    end
                end
            end

        end

        function setROIColorByIndex(app, roiIdx, newColor)
            %SETROICOLORBYINDEX Apply ROI color to data model and runtime graphics.

            if isempty(roiIdx) || roiIdx < 1 || roiIdx > numel(app.ROIList)
                return
            end

            newColor = double(newColor(:).');
            newColor = min(max(newColor, 0), 1);

            app.ROIList(roiIdx).color = newColor;
            app.ROIList(roiIdx).modifiedOn = datetime('now');

            runtime = app.ROIList(roiIdx).runtime;

            if isfield(runtime, 'ROIHandle')
                app.applyColorToGraphicsHandle(runtime.ROIHandle, newColor);
            end

            if isfield(runtime, 'traceLine')
                app.applyColorToGraphicsHandle(runtime.traceLine, newColor);
            end

            if isfield(runtime, 'tracePatch')
                app.applyColorToGraphicsHandle(runtime.tracePatch, newColor);
            end

        end

        function applyColorToGraphicsHandle(app, h, newColor)
            %APPLYCOLORTOGRAPHICSHANDLE Apply color to common graphics handle types.

            if ~app.isUsableGraphicsHandle(h)
                return
            end

            try
                if isprop(h, 'Color')
                    if strcmp(h.Color,'none');return;end
                    h.Color = newColor;
                end
            catch
            end

            try
                if isprop(h, 'FaceColor')
                    if strcmp(h.FaceColor,'none');return;end
                    h.FaceColor = newColor;
                end
            catch
            end

            try
                if isprop(h, 'EdgeColor')
                    if strcmp(h.EdgeColor,'none');return;end
                    h.EdgeColor = newColor;
                end
            catch
            end

        end

        function updateROITraceDisplayName(app, roiIdx)
            %UPDATEROITRACEDISPLAYNAME Update ROI trace label after rename.

            if isempty(roiIdx) || roiIdx < 1 || roiIdx > numel(app.ROIList)
                return
            end

            if ~isfield(app.ROIList(roiIdx), 'runtime') || ...
                    ~isfield(app.ROIList(roiIdx).runtime, 'traceLine')
                return
            end

            h = app.ROIList(roiIdx).runtime.traceLine;

            if app.isUsableGraphicsHandle(h) && isprop(h, 'DisplayName')
                try
                    h.DisplayName = app.ROIList(roiIdx).name;
                catch
                end
            end

        end

        function uniqueName = makeUniqueROIName(app, requestedName, excludeID)
            %MAKEUNIQUEROINAME Return unique ROI name in current app session.

            if nargin < 3
                excludeID = [];
            end

            existingNames = strings(0, 1);

            if ~isempty(app.ROIList)
                for iROI = 1:numel(app.ROIList)
                    if ~isempty(excludeID) && isfield(app.ROIList(iROI), 'ID') && ...
                            app.ROIList(iROI).ID == excludeID
                        continue
                    end

                    if isfield(app.ROIList(iROI), 'name') && ~isempty(app.ROIList(iROI).name)
                        existingNames(end+1, 1) = string(app.ROIList(iROI).name); %#ok<AGROW>
                    end
                end
            end

            uniqueName = app.makeUniqueNameAgainstList(requestedName, existingNames);

        end

        function createNewROI(app, ROIType)
            %CREATENEWROI Interactively draw and register a new ROI.
            %
            %   createNewROI(app, ROIType)
            %
            %   ROIType must be one of:
            %       'rectangle'
            %       'ellipse'
            %       'polygon'
            %
            %   Temporary image.roi objects are used for drawing/editing. Once the
            %   user confirms the ROI, it is converted to:
            %       - polyshape for static geometry
            %       - logical mask for pixel membership
            %       - static overlay for display
            %
            %   Confirmation:
            %       - double-click ROI body
            %       - press Enter
            %
            %   Rename while drawing:
            %       - single-click ROI label when MATLAB exposes label hit info
            %       - press N as reliable fallback
            %
            %   Cancel:
            %       - press Escape
            %       - delete ROI through default context menu

            ROIType = lower(char(string(ROIType)));

            if ~ismember(ROIType, {'rectangle', 'ellipse', 'polygon'})
                error('DataViewer:InvalidROIType', ...
                    'ROIType must be rectangle, ellipse, or polygon.');
            end

            if ~app.hasData()
                app.setStatusMessage('Load image data before creating ROIs.');
                return
            end

            sz = app.getDataSize();
            Ny = sz(1);
            Nx = sz(2);

            % Drawing is intentionally unbounded. Commit-time rasterization determines
            % whether any part of the ROI falls inside the image/logical mask.
            drawingArea = 'unlimited';

            previousMode = app.InteractionMode;
            previousKeyFcn = app.UIFigure.WindowKeyPressFcn;


            hROI = [];
            movingListener = [];
            movedListener = [];
            clickedListener = [];
            deletingListener = [];

            bDone = false;
            bConfirmed = false;
            bCancelled = false;
            bCleanedUp = false;

            roiColor = getNextROIColor(false);  % Advance only after confirmed creation.
            roiID = getNextROIID();
            roiName = app.makeUniqueROIName(sprintf('ROI_%d', roiID), []);

            lastStableSize = [];

            try
                app.setInteractionMode('drawingROI');
                app.UIFigure.WindowKeyPressFcn = @onKeyPress;

                app.setStatusMessage(sprintf( ...
                    ['Draw %s ROI. Double-click ROI or press Enter to confirm. ' ...
                    'Press N to rename. Press Escape to cancel.'], ROIType));

                hold(app.ImageAxes, 'on');

                switch ROIType
                    case 'rectangle'
                        hROI = drawrectangle(app.ImageAxes, ...
                            'Color', roiColor, ...
                            'FaceAlpha', 0.2, ...
                            'LineWidth', 1.5, ...
                            'DrawingArea', drawingArea, ...
                            'InteractionsAllowed', 'all');

                    case 'ellipse'
                        hROI = drawellipse(app.ImageAxes, ...
                            'Color', roiColor, ...
                            'FaceAlpha', 0.2, ...
                            'LineWidth', 1.5, ...
                            'DrawingArea', drawingArea, ...
                            'InteractionsAllowed', 'all');

                    case 'polygon'
                        hROI = drawpolygon(app.ImageAxes, ...
                            'Color', roiColor, ...
                            'FaceAlpha', 0.2, ...
                            'LineWidth', 1.5, ...
                            'DrawingArea', drawingArea, ...
                            'InteractionsAllowed', 'all');
                end

                if isempty(hROI) || ~isvalid(hROI)
                    cleanupTemporaryROI();
                    app.setStatusMessage('ROI creation cancelled.');
                    return
                end

                configureInteractiveROI(hROI);
                restoreROILabel();
                % Preserve MATLAB's built-in ROI context menu and replace only the built-in
                % "Delete <shape>" item with an explicit cancel action.
                app.modifyROICreationContextMenu(hROI, @(src, evt) cancelROIFromContextMenu());
                lastStableSize = getUnrotatedBoundingBoxSize(hROI);

                movingListener = addlistener(hROI, 'MovingROI', @onMovingROI);
                movedListener = addlistener(hROI, 'ROIMoved', @onMovedROI);
                clickedListener = addlistener(hROI, 'ROIClicked', @onROIClicked);

                try
                    deletingListener = addlistener(hROI, 'DeletingROI', @onDeletingROI);
                catch
                    deletingListener = [];
                end

                while ~bDone && isvalid(app.UIFigure) && ~isempty(hROI) && isvalid(hROI)
                    drawnow
                    pause(0.02)
                end

                if ~isvalid(app.UIFigure)
                    cleanupTemporaryROI();
                    return
                end

                if isempty(hROI) || ~isvalid(hROI)
                    bCancelled = true;
                end

                if bCancelled || ~bConfirmed
                    cleanupTemporaryROI();
                    app.setStatusMessage('ROI creation cancelled.');
                    return
                end

                % -----------------------------------------------------------------
                % Convert interactive ROI into persistent/static ROI representation.
                % Apply the active logical mask immediately so stored ROIs never
                % include pixels outside the user-defined processing mask.
                % -----------------------------------------------------------------
                rawMask = app.createMaskFromROIObject(hROI);

                if ~isequal(size(rawMask), [Ny Nx])
                    cleanupTemporaryROI();
                    error('DataViewer:InvalidROIMaskSize', ...
                        'ROI mask size does not match image size.');
                end

                rawMask = logical(rawMask);
                mask = app.clipROIMaskToActiveLogicalMask(rawMask);

                if ~any(mask(:))
                    cleanupTemporaryROI();
                    app.setStatusMessage('ROI creation cancelled: ROI is outside the logical mask.');
                    return
                end

                bWasClippedByMask = ~isequal(mask, rawMask);

                if bWasClippedByMask
                    % Clipped rectangles/ellipses are no longer true primitives.
                    % If clipping splits the ROI into multiple islands, preserve all islands as
                    % independent polygon ROIs named <ROINAME>_part<N>.
                    componentMasks = app.maskToConnectedComponentMasks(mask);

                    if isempty(componentMasks)
                        cleanupTemporaryROI();
                        app.setStatusMessage('ROI creation cancelled: ROI is outside the logical mask.');
                        return
                    end

                    if numel(componentMasks) > 1
                        cleanupTemporaryROI();

                        nAdded = app.addPolygonROIsFromMaskComponents( ...
                            componentMasks, roiName, roiColor, true);

                        if nAdded < 1
                            app.setStatusMessage('ROI creation cancelled: split ROI components were invalid.');
                            return
                        end

                        advanceNextROIColor();

                        app.setSelectedROI(NaN);
                        app.refreshROITable();
                        app.refreshROITraces();
                        app.refreshEventPatches();
                        app.stackROITraceGraphics();
                        app.updateGUIEnabledState();

                        app.setStatusMessage(sprintf('Created split ROI "%s" as %d parts.', roiName, nAdded));
                        return
                    end

                    mask = componentMasks{1};
                    ROIType = 'polygon';
                    verticesXY = app.maskToSimplifiedPolygonVertices(mask);
                    roiParameters = app.makePolygonROIParametersFromVertices(verticesXY);
                else
                    verticesXY = getVerticesFromROIObject(hROI, ROIType);
                    roiParameters = getROIParametersFromObject(hROI, ROIType);
                end

                verticesXY = app.cleanROIVertices(verticesXY);

                if size(verticesXY, 1) < 3
                    cleanupTemporaryROI();
                    app.setStatusMessage('ROI creation cancelled: invalid ROI geometry.');
                    return
                end

                pgon = polyshape(verticesXY(:, 1), verticesXY(:, 2), 'Simplify', true);

                if isempty(pgon.Vertices)
                    cleanupTemporaryROI();
                    app.setStatusMessage('ROI creation cancelled: invalid polyshape.');
                    return
                end
                ROI = struct();

                ROI.name = roiName;
                ROI.type = ROIType;
                ROI.DOC = datetime('now');
                ROI.modifiedOn = ROI.DOC;
                ROI.color = roiColor;
                ROI.notes = '';

                ROI.geometry = struct();
                ROI.geometry.polyshape = pgon;
                ROI.geometry.verticesXY_px = verticesXY;
                ROI.geometry.ROIType = ROIType;
                ROI.geometry.ROIParameters = roiParameters;
                ROI.mask = mask;
                ROI.stats = app.computeROIStatsFromMask(mask);

                ROI.ID = roiID;
                ROI.runtime = struct();
                ROI.runtime.visible = true;
                ROI.runtime.ROIHandle = gobjects(1);
                ROI.runtime.editHandle = gobjects(1);
                ROI.runtime.traceLine = gobjects(1);
                ROI.runtime.tracePatch = gobjects(1);

                ROI.runtime.trace = struct();
                ROI.runtime.trace.XData = [];
                ROI.runtime.trace.Mean = [];
                ROI.runtime.trace.SEM = [];
                ROI.runtime.trace.Mode = '';
                ROI.runtime.trace.Status = 'not computed';

                cleanupTemporaryROI();

                ROI.runtime.ROIHandle = createStaticROIOverlay(ROI);

                if isempty(app.ROIList)
                    app.ROIList = ROI;
                else
                    app.ROIList(end+1) = ROI;
                end

                advanceNextROIColor();

                % Newly created ROIs are not automatically selected. The ROI table Select
                % checkbox is the source of truth for edit/delete intent.
                app.setROISelectedStateByIndex(numel(app.ROIList), false);
                app.setSelectedROI(NaN);

                app.refreshROITable();
                app.refreshROITraces();

                % Creating ROI trace graphics can disturb the PlotAxes legend state.
                % Rebuild event patches/legend explicitly so the legend is not lost
                % until the next image click.
                app.refreshEventPatches();
                app.stackROITraceGraphics();

                app.updateGUIEnabledState();

                drawnow
                app.setStatusMessage(sprintf('Created ROI "%s".', roiName));

            catch ME
                cleanupTemporaryROI();
                app.setStatusMessage(sprintf('ROI creation failed: %s', ME.message));
            end

            function configureInteractiveROI(h)
                %CONFIGUREINTERACTIVEROI Apply visual and interaction defaults.

                setPropertyIfAvailable(h, 'LineStyle', '-');
                setPropertyIfAvailable(h, 'FaceAlpha', 0.2);
                setPropertyIfAvailable(h, 'Color', roiColor);
                setPropertyIfAvailable(h, 'LineWidth', 1.5);
                setPropertyIfAvailable(h, 'Rotatable', true);

                try
                    h.DrawingArea = drawingArea;
                catch
                end

                try
                    h.Label = roiName;
                catch
                end

                try
                    if isprop(h, 'LabelVisible')
                        h.LabelVisible = 'on';
                    end
                catch
                end
            end

            function setPropertyIfAvailable(h, propName, propValue)
                try
                    if isprop(h, propName)
                        h.(propName) = propValue;
                    end
                catch
                end
            end

            function restoreROILabel()
                try
                    if ~isempty(hROI) && isvalid(hROI)
                        hROI.Label = roiName;
                        if isprop(hROI, 'LabelVisible')
                            hROI.LabelVisible = 'on';
                        end
                    end
                catch
                end
            end

            function onMovingROI(src, ~)
                currentSize = getUnrotatedBoundingBoxSize(src);

                if isempty(currentSize) || isempty(lastStableSize)
                    restoreROILabel();
                    return
                end

                isResize = any(abs(currentSize - lastStableSize) > 0.25);

                if isResize
                    [scaleFactor, unitText] = getAxisUnitScale();
                    displaySize = currentSize .* scaleFactor;

                    try
                        src.Label = sprintf('%s | %.3g x %.3g %s', ...
                            roiName, displaySize(1), displaySize(2), unitText);
                    catch
                    end
                else
                    restoreROILabel();
                end
            end

            function onMovedROI(~, ~)
                if ~isempty(hROI) && isvalid(hROI)
                    lastStableSize = getUnrotatedBoundingBoxSize(hROI);
                end

                restoreROILabel();
            end

            function onROIClicked(src, evt)
                %ONROICLICKED Rename label on single click, confirm ROI on double-click body.

                if isROILabelClick(evt)
                    renameROIWhileDrawing(src);
                    return
                end

                if isprop(evt, 'SelectionType') && strcmpi(evt.SelectionType, 'double')
                    bConfirmed = true;
                    bDone = true;
                end
            end

            function onDeletingROI(~, ~)
                bCancelled = true;
                bDone = true;
            end

            function cancelROIFromContextMenu()
                bConfirmed = false;
                bCancelled = true;
                bDone = true;
            end
            function onKeyPress(~, evt)
                switch lower(evt.Key)
                    case {'return', 'enter'}
                        bConfirmed = true;
                        bDone = true;

                    case 'escape'
                        bConfirmed = false;
                        bCancelled = true;
                        bDone = true;

                    case 'n'
                        if ~isempty(hROI) && isvalid(hROI)
                            renameROIWhileDrawing(hROI);
                        end
                end
            end

            function renameROIWhileDrawing(src)
                newName = promptForROIName(roiName);

                if isempty(newName)
                    restoreROILabel();
                    return
                end

                roiName = app.makeUniqueROIName(newName, []);

                try
                    src.Label = roiName;
                    if isprop(src, 'LabelVisible')
                        src.LabelVisible = 'on';
                    end
                catch
                end
            end

            function roiColorOut = getNextROIColor(bAdvance)
                %GETNEXTROICOLOR Return next color from app.ROIColorList.

                if nargin < 1
                    bAdvance = false;
                end

                if ~isempty(app.ROIColorList) && ...
                        isnumeric(app.ROIColorList) && size(app.ROIColorList, 2) == 3
                    colorList = double(app.ROIColorList);
                else
                    colorList = lines(50);
                end

                colorList = min(max(colorList, 0), 1);

                if ~isempty(app.ROINextColorIndex) && ...
                        isfinite(app.ROINextColorIndex) && app.ROINextColorIndex >= 1
                    idxColor = round(app.ROINextColorIndex);
                else
                    idxColor = 1;
                end

                idxColor = mod(idxColor - 1, size(colorList, 1)) + 1;
                roiColorOut = colorList(idxColor, :);

                if bAdvance
                    advanceNextROIColor();
                end
            end

            function advanceNextROIColor()
                %ADVANCENEXTROICOLOR Advance app ROI color index after confirmed creation.

                if isempty(app.ROIColorList) || ...
                        ~isnumeric(app.ROIColorList) || size(app.ROIColorList, 2) ~= 3
                    nColors = 50;
                else
                    nColors = size(app.ROIColorList, 1);
                end

                if isempty(app.ROINextColorIndex) || ~isfinite(app.ROINextColorIndex)
                    app.ROINextColorIndex = 1;
                end

                app.ROINextColorIndex = round(app.ROINextColorIndex) + 1;

                if app.ROINextColorIndex > nColors
                    app.ROINextColorIndex = 1;
                end
            end


            function tf = isROILabelClick(evt)
                %ISROILABELCLICK Best-effort detection of ROI label click.

                tf = false;

                candidateParts = {'SelectedPart', 'HitPart', 'ClickedPart', 'Part'};

                for iPart = 1:numel(candidateParts)
                    propName = candidateParts{iPart};

                    try
                        if isprop(evt, propName)
                            value = lower(char(string(evt.(propName))));
                            tf = contains(value, 'label');
                            if tf
                                return
                            end
                        end
                    catch
                    end
                end

                try
                    if isprop(evt, 'HitObject') && ~isempty(evt.HitObject)
                        h = evt.HitObject;
                        if isprop(h, 'Type')
                            tf = contains(lower(char(string(h.Type))), 'text');
                        end
                    end
                catch
                end
            end

            function newName = promptForROIName(currentName)
                %PROMPTFORROINAME Modal dialog for renaming ROI while drawing.

                newName = '';

                dlg = uifigure( ...
                    'Name', 'Rename ROI', ...
                    'WindowStyle', 'modal', ...
                    'Position', [100 100 320 120], ...
                    'Visible', 'off', ...
                    'CloseRequestFcn', @onCancel);

                grid = uigridlayout(dlg);
                grid.RowHeight = {28, '1x', 32};
                grid.ColumnWidth = {70, '1x'};
                grid.Padding = [12 12 12 12];

                nameLabel = uilabel(grid);
                nameLabel.Text = 'Name';
                nameLabel.Layout.Row = 1;
                nameLabel.Layout.Column = 1;

                nameField = uieditfield(grid, 'text');
                nameField.Value = currentName;
                nameField.Layout.Row = 1;
                nameField.Layout.Column = 2;

                statusLabel = uilabel(grid);
                statusLabel.Text = '';
                statusLabel.FontColor = [0.65 0 0];
                statusLabel.Layout.Row = 2;
                statusLabel.Layout.Column = [1 2];

                buttonGrid = uigridlayout(grid);
                buttonGrid.RowHeight = {'1x'};
                buttonGrid.ColumnWidth = {'1x', '1x'};
                buttonGrid.Padding = [0 0 0 0];
                buttonGrid.Layout.Row = 3;
                buttonGrid.Layout.Column = [1 2];

                okButton = uibutton(buttonGrid, 'push');
                okButton.Text = 'OK';
                okButton.Layout.Row = 1;
                okButton.Layout.Column = 1;
                okButton.ButtonPushedFcn = @onOK;

                cancelButton = uibutton(buttonGrid, 'push');
                cancelButton.Text = 'Cancel';
                cancelButton.Layout.Row = 1;
                cancelButton.Layout.Column = 2;
                cancelButton.ButtonPushedFcn = @onCancel;

                if exist('placeAppInsideCaller', 'file') == 2
                    placeAppInsideCaller(app, dlg, 'center');
                end

                dlg.Visible = 'on';
                drawnow
                uiwait(dlg);

                if isvalid(dlg)
                    delete(dlg);
                end

                function onOK(~, ~)
                    requestedName = strtrim(char(string(nameField.Value)));

                    if isempty(requestedName)
                        statusLabel.Text = 'ROI name cannot be empty.';
                        return
                    end

                    newName = requestedName;
                    uiresume(dlg);
                end

                function onCancel(~, ~)
                    newName = '';
                    uiresume(dlg);
                end
            end

            function roiID = getNextROIID()
                %GETNEXTROIID Return next runtime ROI ID.

                if isempty(app.ROIList)
                    roiID = 1;
                    return
                end

                existingIDs = [app.ROIList.ID];
                existingIDs = existingIDs(isfinite(existingIDs));

                if isempty(existingIDs)
                    roiID = 1;
                else
                    roiID = max(existingIDs) + 1;
                end
            end

            function verticesXY = getVerticesFromROIObject(h, typeName)
                %GETVERTICESFROMROIOBJECT Extract/synthesize ROI vertices in image XY px.

                verticesXY = [];

                switch lower(typeName)
                    case 'rectangle'
                        if isprop(h, 'Position') && ~isempty(h.Position)
                            pos = double(h.Position);
                            verticesXY = rectanglePositionToVertices(pos);

                            angleDeg = getNumericProperty(h, 'RotationAngle', 0);
                            verticesXY = rotateVertices(verticesXY, ...
                                [pos(1) + pos(3)/2, pos(2) + pos(4)/2], angleDeg);
                        end

                    case 'ellipse'
                        verticesXY = ellipseToVertices(h, 128);

                    case 'polygon'
                        if isprop(h, 'Position') && ~isempty(h.Position)
                            verticesXY = double(h.Position);
                        elseif isprop(h, 'Vertices') && ~isempty(h.Vertices)
                            verticesXY = double(h.Vertices);
                        end
                end

                verticesXY = cleanVertices(verticesXY);
            end

            function verticesXY = rectanglePositionToVertices(pos)
                x = pos(1);
                y = pos(2);
                w = pos(3);
                h = pos(4);

                verticesXY = [ ...
                    x,     y; ...
                    x + w, y; ...
                    x + w, y + h; ...
                    x,     y + h];
            end

            function verticesXY = ellipseToVertices(h, nPts)
                theta = linspace(0, 2*pi, nPts + 1).';
                theta(end) = [];

                if isprop(h, 'Center') && isprop(h, 'SemiAxes') && ...
                        ~isempty(h.Center) && ~isempty(h.SemiAxes)
                    centerXY = double(h.Center(:).');
                    semiAxes = double(h.SemiAxes(:).');
                elseif isprop(h, 'Position') && ~isempty(h.Position)
                    pos = double(h.Position);
                    centerXY = [pos(1) + pos(3)/2, pos(2) + pos(4)/2];
                    semiAxes = [pos(3)/2, pos(4)/2];
                else
                    verticesXY = [];
                    return
                end

                verticesXY = [semiAxes(1) .* cos(theta), semiAxes(2) .* sin(theta)];

                angleDeg = getNumericProperty(h, 'RotationAngle', 0);
                verticesXY = rotateVertices(verticesXY + centerXY, centerXY, angleDeg);
            end

            function value = getNumericProperty(h, propName, defaultValue)
                value = defaultValue;

                try
                    if isprop(h, propName) && ~isempty(h.(propName))
                        candidate = double(h.(propName));
                        if isscalar(candidate) && isfinite(candidate)
                            value = candidate;
                        end
                    end
                catch
                    value = defaultValue;
                end
            end

            function verticesOut = rotateVertices(verticesIn, centerXY, angleDeg)
                if isempty(verticesIn) || angleDeg == 0
                    verticesOut = verticesIn;
                    return
                end

                % images.roi RotationAngle is expressed in image coordinates, where
                % Y increases downward. polyshape/plot use the axis coordinate system,
                % so the sign must be inverted to avoid mirrored rotation.
                angleRad = -deg2rad(angleDeg);
                R = [cos(angleRad), -sin(angleRad); sin(angleRad), cos(angleRad)];

                verticesOut = (verticesIn - centerXY) * R.' + centerXY;
            end

            function verticesXY = cleanVertices(verticesXY)
                if isempty(verticesXY) || size(verticesXY, 2) ~= 2
                    verticesXY = zeros(0, 2);
                    return
                end

                verticesXY = double(verticesXY);
                verticesXY = verticesXY(all(isfinite(verticesXY), 2), :);

                if size(verticesXY, 1) > 1
                    duplicate = [false; all(abs(diff(verticesXY, 1, 1)) < eps, 2)];
                    verticesXY(duplicate, :) = [];
                end
            end

            function params = getROIParametersFromObject(h, typeName)
                %GETROIPARAMETERSFROMOBJECT Store parameters needed for reconstruction.

                params = struct();
                params.ROIType = typeName;

                propList = { ...
                    'Position', ...
                    'Center', ...
                    'SemiAxes', ...
                    'RotationAngle', ...
                    'Vertices'};

                for iProp = 1:numel(propList)
                    propName = propList{iProp};

                    try
                        if isprop(h, propName) && ~isempty(h.(propName))
                            params.(propName) = h.(propName);
                        end
                    catch
                    end
                end
            end

            function sizeXY = getUnrotatedBoundingBoxSize(h)
                %GETUNROTATEDBOUNDINGBOXSIZE Estimate unrotated bounding-box size.

                sizeXY = [];

                try
                    if isprop(h, 'SemiAxes') && ~isempty(h.SemiAxes)
                        semiAxes = double(h.SemiAxes(:).');
                        if numel(semiAxes) >= 2
                            sizeXY = 2 .* semiAxes(1:2);
                            return
                        end
                    end
                catch
                end

                try
                    if isprop(h, 'Position') && ~isempty(h.Position)
                        pos = double(h.Position);
                        if numel(pos) >= 4
                            sizeXY = abs(pos(3:4));
                            return
                        end
                    end
                catch
                end

                try
                    if isprop(h, 'Vertices') && ~isempty(h.Vertices)
                        verts = double(h.Vertices);
                    elseif isprop(h, 'Position') && size(h.Position, 2) == 2
                        verts = double(h.Position);
                    else
                        verts = [];
                    end

                    if ~isempty(verts)
                        sizeXY = [range(verts(:, 1)), range(verts(:, 2))];
                    end
                catch
                end
            end

            function hStatic = createStaticROIOverlay(ROI)
                %CREATESTATICROIOVERLAY Draw passive persistent/static ROI overlay.

                hStatic = app.createStaticROIOverlayFromROI(ROI);
            end

            function onStaticROIClick(roiID, evt)
                app.setSelectedROI(roiID);

                if isprop(evt, 'SelectionType') && strcmpi(evt.SelectionType, 'double')
                    app.setStatusMessage('ROI editing is not implemented yet.');
                end
            end

            function [scaleFactor, unitText] = getAxisUnitScale()
                %GETAXISUNITSCALE Return size conversion for ROI dimension label.

                scaleFactor = 1;
                unitText = 'px';

                if isempty(app.DataParams) || ~isfield(app.DataParams, 'view')
                    return
                end

                viewInfo = app.DataParams.view;

                if ~isfield(viewInfo, 'pixelSize_px_per_mm') || isempty(viewInfo.pixelSize_px_per_mm)
                    return
                end

                pxPerMm = double(viewInfo.pixelSize_px_per_mm);

                if isscalar(pxPerMm) && isfinite(pxPerMm) && pxPerMm > 0
                    scaleFactor = 1 ./ pxPerMm;
                    unitText = 'mm';
                end
            end

            function cleanupTemporaryROI()
                %CLEANUPTEMPORARYROI Delete listeners/temp ROI and restore app state.

                if bCleanedUp
                    return
                end

                bCleanedUp = true;

                try
                    if ~isempty(movingListener) && isvalid(movingListener)
                        delete(movingListener);
                    end
                catch
                end

                try
                    if ~isempty(movedListener) && isvalid(movedListener)
                        delete(movedListener);
                    end
                catch
                end

                try
                    if ~isempty(clickedListener) && isvalid(clickedListener)
                        delete(clickedListener);
                    end
                catch
                end

                try
                    if ~isempty(deletingListener) && isvalid(deletingListener)
                        delete(deletingListener);
                    end
                catch
                end

                try
                    if ~isempty(hROI) && isvalid(hROI)
                        delete(hROI);
                    end
                catch
                end

                try
                    if isvalid(app.UIFigure)
                        app.UIFigure.WindowKeyPressFcn = previousKeyFcn;
                    end
                catch
                end

                try
                    app.setInteractionMode(previousMode);
                catch
                    app.InteractionMode = previousMode;
                end

            end

        end

        function selectROIForEditingByID(app, roiID)
            %SELECTROIFOREDITINGBYID Enter interactive edit mode for one ROI.
            %
            %   Table selection alone does not call this method. Single-ROI edit is
            %   triggered from a double-click inside the ROI or another explicit edit
            %   command.

            if isempty(roiID) || ~isfinite(roiID)
                return
            end

            if ~app.hasData()
                return
            end

            if strcmp(app.InteractionMode, 'playingMovie')
                app.setStatusMessage('Stop movie playback before editing ROIs.');
                app.refreshROITable();
                return
            end

            if strcmp(app.InteractionMode, 'editingROI')
                app.setStatusMessage('Finish the current ROI edit before starting another edit.');
                app.refreshROITable();
                return
            end

            roiIdx = app.getROIIndexByID(roiID);
            if isempty(roiIdx)
                app.refreshROITable();
                return
            end

            app.clearROISelectionState();
            app.setROISelectedStateByIndex(roiIdx, true);
            app.setSelectedROI(roiID);
            app.beginROIEditingByIndex(roiIdx);

        end

        function roiIDList = getSelectedROIIDList(app)
            %GETSELECTEDROIIDLIST Return selected ROI IDs in user selection order.
            %
            %   ROISelectionOrder stores the chronological checkbox order. Any selected ROI
            %   missing from that order is appended in ROIList order as a safety fallback.

            roiIDList = [];

            if isempty(app.ROIList)
                app.ROISelectionOrder = [];
                return
            end

            selectedIDs = [];

            for iROI = 1:numel(app.ROIList)
                if isfield(app.ROIList(iROI), 'runtime') && ...
                        isfield(app.ROIList(iROI).runtime, 'selected') && ...
                        app.ROIList(iROI).runtime.selected

                    selectedIDs(end+1) = app.ROIList(iROI).ID; %#ok<AGROW>
                end
            end

            selectedIDs = unique(double(selectedIDs(:).'), 'stable');
            selectedIDs = selectedIDs(isfinite(selectedIDs));

            if isempty(selectedIDs)
                app.ROISelectionOrder = [];
                return
            end

            app.ROISelectionOrder = double(app.ROISelectionOrder(:).');
            app.ROISelectionOrder = app.ROISelectionOrder(isfinite(app.ROISelectionOrder));

            orderedSelected = app.ROISelectionOrder(ismember(app.ROISelectionOrder, selectedIDs));
            missingSelected = selectedIDs(~ismember(selectedIDs, orderedSelected));

            roiIDList = [orderedSelected, missingSelected];

            % Keep runtime order clean.
            app.ROISelectionOrder = roiIDList;

        end

        function setROISelectedStateByIndex(app, roiIdx, tf)
            %SETROISELECTEDSTATEBYINDEX Set temporary table-selection state.

            if isempty(roiIdx) || roiIdx < 1 || roiIdx > numel(app.ROIList)
                return
            end

            if ~isfield(app.ROIList(roiIdx), 'runtime') || isempty(app.ROIList(roiIdx).runtime)
                app.ROIList(roiIdx).runtime = struct();
            end

            app.ROIList(roiIdx).runtime.selected = logical(tf);

        end

        function beginROIEditingByIndex(app, roiIdx)
            %BEGINROIEDITINGBYINDEX Replace static ROI overlay with editable ROI object.
            %
            %   Edit mode deletes the stored polyshape overlay, rebuilds an images.roi
            %   object from the stored ROI parameters, and attaches editing callbacks
            %   equivalent to the creation workflow. Confirm with Enter or double-click;
            %   cancel with Escape or by deleting the temporary ROI object.

            if isempty(roiIdx) || roiIdx < 1 || roiIdx > numel(app.ROIList)
                return
            end

            if ~app.hasData()
                return
            end

            ROI = app.ROIList(roiIdx);
            roiID = ROI.ID;

            if isfield(ROI, 'runtime') && isfield(ROI.runtime, 'editHandle') && ...
                    app.isUsableGraphicsHandle(ROI.runtime.editHandle)
                app.refreshROITable();
                return
            end

            if ~isfield(app.ROIList(roiIdx), 'runtime') || isempty(app.ROIList(roiIdx).runtime)
                app.ROIList(roiIdx).runtime = struct();
            end

            previousMode = app.InteractionMode;
            previousKeyFcn = app.UIFigure.WindowKeyPressFcn;

            % Delete the passive polyshape/patch overlay. It will be recreated on
            % confirm or cancel.
            if isfield(app.ROIList(roiIdx).runtime, 'ROIHandle') && ...
                    app.isUsableGraphicsHandle(app.ROIList(roiIdx).runtime.ROIHandle)
                try
                    delete(app.ROIList(roiIdx).runtime.ROIHandle);
                catch
                end
            end
            app.ROIList(roiIdx).runtime.ROIHandle = gobjects(1);

            try
                hEdit = app.createEditableROIObjectFromStoredROI(ROI);
            catch ME
                app.ROIList(roiIdx).runtime.ROIHandle = app.createStaticROIOverlayFromROI(ROI);
                app.setStatusMessage(sprintf('ROI edit failed: %s', ME.message));
                app.refreshROITable();
                return
            end

            if isempty(hEdit) || ~isvalid(hEdit)
                app.ROIList(roiIdx).runtime.ROIHandle = app.createStaticROIOverlayFromROI(ROI);
                app.setStatusMessage('ROI edit failed: could not create editable ROI object.');
                app.refreshROITable();
                return
            end

            app.ROIList(roiIdx).runtime.editHandle = hEdit;
            app.createEditableROIContextMenu(hEdit, roiID);
            app.ROIList(roiIdx).runtime.editName = ROI.name;
            app.ROIList(roiIdx).runtime.editOriginalName = ROI.name;
            app.ROIList(roiIdx).runtime.editPreviousMode = previousMode;
            app.ROIList(roiIdx).runtime.editPreviousKeyFcn = previousKeyFcn;

            listeners = {};
            listeners{end+1} = addlistener(hEdit, 'MovingROI', ...
                @(src, evt) app.onEditableROIMoving(src, evt, roiID));
            listeners{end+1} = addlistener(hEdit, 'ROIMoved', ...
                @(src, evt) app.onEditableROIMoved(src, evt, roiID));
            listeners{end+1} = addlistener(hEdit, 'ROIClicked', ...
                @(src, evt) app.onEditableROIClicked(src, evt, roiID));

            try
                listeners{end+1} = addlistener(hEdit, 'DeletingROI', ...
                    @(src, evt) app.onEditableROIDeleting(src, evt, roiID));
            catch
            end

            app.ROIList(roiIdx).runtime.editListeners = listeners;
            app.UIFigure.WindowKeyPressFcn = @(src, evt) app.onActiveROIEditKeyPress(src, evt, roiID);

            app.setInteractionMode('editingROI');
            app.refreshROITable();
            app.setStatusMessage(sprintf( ...
                'Editing ROI "%s". Double-click or press Enter to confirm; press N to rename; press Escape to cancel.', ...
                ROI.name));

        end

        function hEdit = createEditableROIObjectFromStoredROI(app, ROI)
            %CREATEEDITABLEROIOBJECTFROMSTOREDROI Rebuild an images.roi object.
            %
            %   hEdit = createEditableROIObjectFromStoredROI(app, ROI) recreates the
            %   editable images.roi object from the stored ROI geometry. The object is
            %   initialized in-place; it must never start a new interactive drawing
            %   operation during edit mode.
            %
            %   Important:
            %       Group-transformed rectangles and ellipses must remain editable as
            %       rectangles and ellipses. Do not coerce them to polygon just because
            %       ROIParameters contains Vertices.

            hEdit = [];

            if ~isfield(ROI, 'type') || isempty(ROI.type)
                error('DataViewer:InvalidROIType', 'ROI type is missing.');
            end

            roiType = lower(char(string(ROI.type)));

            sz = app.getDataSize();
            Nx = sz(2);
            Ny = sz(1);
            % Editing is intentionally unbounded. Commit-time rasterization determines
            % whether any part of the ROI remains valid.

            drawingArea = 'unlimited';

            roiColor = [1 0 0];
            if isfield(ROI, 'color') && numel(ROI.color) == 3
                roiColor = min(max(double(ROI.color(:).'), 0), 1);
            end

            params = struct();
            if isfield(ROI, 'geometry') && isfield(ROI.geometry, 'ROIParameters') && ...
                    isstruct(ROI.geometry.ROIParameters)
                params = ROI.geometry.ROIParameters;
            end

            hold(app.ImageAxes, 'on');

            switch roiType
                case 'rectangle'
                    pos = [];

                    if isfield(params, 'Position') && ~isempty(params.Position)
                        pos = double(params.Position(:).');
                    elseif isfield(ROI, 'geometry') && isfield(ROI.geometry, 'verticesXY_px') && ...
                            ~isempty(ROI.geometry.verticesXY_px)
                        verticesXY = double(ROI.geometry.verticesXY_px);
                        verticesXY = verticesXY(all(isfinite(verticesXY), 2), :);

                        if size(verticesXY, 1) >= 4
                            minXY = min(verticesXY(:, 1:2), [], 1);
                            maxXY = max(verticesXY(:, 1:2), [], 1);
                            boxSize = maxXY - minXY;
                            pos = [minXY(1), minXY(2), boxSize(1), boxSize(2)];
                        end
                    end

                    if numel(pos) < 4 || any(~isfinite(pos(1:4))) || any(pos(3:4) <= 0)
                        error('DataViewer:InvalidROIParameters', ...
                            'Rectangle ROI is missing valid Position parameters.');
                    end

                    hEdit = drawrectangle(app.ImageAxes, ...
                        'Position', pos(1:4), ...
                        'Color', roiColor, ...
                        'FaceAlpha', 0.2, ...
                        'LineWidth', 1.5, ...
                        'DrawingArea', drawingArea, ...
                        'InteractionsAllowed', 'all');

                    if isfield(params, 'RotationAngle') && ~isempty(params.RotationAngle)
                        app.setROIObjectPropertyIfAvailable(hEdit, 'RotationAngle', params.RotationAngle);
                    end

                case 'ellipse'
                    centerXY = [];
                    semiAxes = [];
                    angleDeg = 0;

                    if isfield(params, 'RotationAngle') && ~isempty(params.RotationAngle)
                        angleCandidate = double(params.RotationAngle);
                        if isscalar(angleCandidate) && isfinite(angleCandidate)
                            angleDeg = angleCandidate;
                        end
                    end

                    if isfield(params, 'Center') && isfield(params, 'SemiAxes') && ...
                            ~isempty(params.Center) && ~isempty(params.SemiAxes)

                        centerXY = double(params.Center(:).');
                        semiAxes = double(params.SemiAxes(:).');

                    elseif isfield(params, 'Position') && ~isempty(params.Position)
                        pos = double(params.Position(:).');

                        if numel(pos) >= 4
                            centerXY = [pos(1) + pos(3)/2, pos(2) + pos(4)/2];
                            semiAxes = abs([pos(3), pos(4)]) ./ 2;
                        end

                    elseif isfield(ROI, 'geometry') && isfield(ROI.geometry, 'verticesXY_px') && ...
                            ~isempty(ROI.geometry.verticesXY_px)

                        verticesXY = double(ROI.geometry.verticesXY_px);
                        verticesXY = verticesXY(all(isfinite(verticesXY), 2), :);

                        if size(verticesXY, 1) >= 3
                            minXY = min(verticesXY(:, 1:2), [], 1);
                            maxXY = max(verticesXY(:, 1:2), [], 1);

                            centerXY = mean([minXY; maxXY], 1);
                            semiAxes = max((maxXY - minXY) ./ 2, 0.5);
                        end
                    end

                    if numel(centerXY) ~= 2 || numel(semiAxes) ~= 2 || ...
                            any(~isfinite(centerXY)) || any(~isfinite(semiAxes)) || ...
                            any(semiAxes <= 0)
                        error('DataViewer:InvalidROIParameters', ...
                            'Ellipse ROI is missing valid Center/SemiAxes parameters.');
                    end

                    centerXY(1) = min(max(centerXY(1), 1), Nx);
                    centerXY(2) = min(max(centerXY(2), 1), Ny);
                    semiAxes = max(double(semiAxes(:).'), 0.5);

                    try
                        hEdit = drawellipse(app.ImageAxes, ...
                            'Center', centerXY, ...
                            'SemiAxes', semiAxes, ...
                            'RotationAngle', angleDeg, ...
                            'Color', roiColor, ...
                            'FaceAlpha', 0.2, ...
                            'LineWidth', 1.5, ...
                            'DrawingArea', drawingArea, ...
                            'InteractionsAllowed', 'all');
                    catch
                        pos = [ ...
                            centerXY(1) - semiAxes(1), ...
                            centerXY(2) - semiAxes(2), ...
                            2 * semiAxes(1), ...
                            2 * semiAxes(2)];

                        hEdit = drawellipse(app.ImageAxes, ...
                            'Position', pos, ...
                            'Color', roiColor, ...
                            'FaceAlpha', 0.2, ...
                            'LineWidth', 1.5, ...
                            'DrawingArea', drawingArea, ...
                            'InteractionsAllowed', 'all');

                        app.setROIObjectPropertyIfAvailable(hEdit, 'RotationAngle', angleDeg);
                    end

                case 'polygon'
                    verticesXY = [];

                    if isfield(params, 'Position') && ~isempty(params.Position) && ...
                            size(params.Position, 2) == 2
                        verticesXY = double(params.Position);
                    elseif isfield(params, 'Vertices') && ~isempty(params.Vertices) && ...
                            size(params.Vertices, 2) == 2
                        verticesXY = double(params.Vertices);
                    elseif isfield(ROI, 'geometry') && isfield(ROI.geometry, 'verticesXY_px')
                        verticesXY = double(ROI.geometry.verticesXY_px);
                    end

                    verticesXY = verticesXY(all(isfinite(verticesXY), 2), :);

                    if isempty(verticesXY) || size(verticesXY, 1) < 3 || size(verticesXY, 2) ~= 2
                        error('DataViewer:InvalidROIParameters', ...
                            'Polygon ROI is missing valid vertices.');
                    end

                    hEdit = drawpolygon(app.ImageAxes, ...
                        'Position', verticesXY, ...
                        'Color', roiColor, ...
                        'FaceAlpha', 0.2, ...
                        'LineWidth', 1.5, ...
                        'DrawingArea', drawingArea, ...
                        'InteractionsAllowed', 'all');

                otherwise
                    error('DataViewer:InvalidROIType', ...
                        'Unsupported ROI type: %s.', roiType);
            end

            app.setROIObjectPropertyIfAvailable(hEdit, 'Rotatable', true);
            app.setROIObjectPropertyIfAvailable(hEdit, 'LineStyle', '-');
            app.setROIObjectPropertyIfAvailable(hEdit, 'Label', char(string(ROI.name)));
            app.setROIObjectPropertyIfAvailable(hEdit, 'LabelVisible', 'on');

        end

        function commitROIEditingByIndex(app, roiIdx)
            %COMMITROIEDITINGBYINDEX Store edited ROI geometry and rebuild static overlay.
            %
            %   The editable ROI object may be moved outside the image because edit
            %   interactions are intentionally unbounded. At commit time, the ROI is
            %   rasterized into image coordinates and clipped to the active logical mask.
            %   If no valid pixels remain, the user is asked whether the now-invalid ROI
            %   should be deleted.

            if isempty(roiIdx) || roiIdx < 1 || roiIdx > numel(app.ROIList)
                return
            end

            if ~isfield(app.ROIList(roiIdx), 'runtime') || ...
                    ~isfield(app.ROIList(roiIdx).runtime, 'editHandle') || ...
                    ~app.isUsableGraphicsHandle(app.ROIList(roiIdx).runtime.editHandle)
                return
            end

            hEdit = app.ROIList(roiIdx).runtime.editHandle;
            roiType = lower(char(string(app.ROIList(roiIdx).type)));

            sz = app.getDataSize();
            Ny = sz(1);
            Nx = sz(2);

            rawMask = app.createMaskFromROIObject(hEdit);
            if ~isequal(size(rawMask), [Ny Nx])
                app.setStatusMessage('ROI edit cancelled: mask size mismatch.');
                app.cancelROIEditingByIndex(roiIdx, true);
                return
            end

            rawMask = logical(rawMask);
            mask = app.clipROIMaskToActiveLogicalMask(rawMask);

            editName = app.ROIList(roiIdx).name;
            if isfield(app.ROIList(roiIdx).runtime, 'editName') && ...
                    ~isempty(app.ROIList(roiIdx).runtime.editName)
                editName = char(string(app.ROIList(roiIdx).runtime.editName));
            end
            editName = app.makeUniqueROIName(editName, app.ROIList(roiIdx).ID);

            if ~any(mask(:))
                promptText = sprintf( ...
                    ['The edited ROI "%s" does not include any valid image pixels ' ...
                    'inside the active logical mask. Delete this ROI?'], ...
                    editName);

                if app.confirmROIDeletion(promptText)
                    roiID = app.ROIList(roiIdx).ID;
                    app.deleteROIsByIDList(roiID, true);
                    app.setStatusMessage(sprintf('Deleted invalid ROI "%s".', editName));
                else
                    app.cancelROIEditingByIndex(roiIdx, false);
                    app.setStatusMessage(sprintf('ROI edit cancelled. Kept original ROI "%s".', editName));
                end
                return
            end

            bWasClipped = ~isequal(mask, rawMask);

            if bWasClipped
                componentMasks = app.maskToConnectedComponentMasks(mask);

                if isempty(componentMasks)
                    promptText = sprintf( ...
                        ['The edited ROI "%s" has no valid mask components after clipping. ' ...
                        'Delete this ROI?'], editName);

                    if app.confirmROIDeletion(promptText)
                        roiID = app.ROIList(roiIdx).ID;
                        app.deleteROIsByIDList(roiID, true);
                        app.setStatusMessage(sprintf('Deleted invalid ROI "%s".', editName));
                    else
                        app.cancelROIEditingByIndex(roiIdx, false);
                        app.setStatusMessage(sprintf('ROI edit cancelled. Kept original ROI "%s".', editName));
                    end
                    return
                end

                if numel(componentMasks) > 1
                    % One edited ROI became multiple disconnected valid regions. Replace the
                    % original ROI with one polygon ROI per component.
                    app.cleanupROIEditRuntimeByIndex(roiIdx);

                    nParts = app.replaceROIByMaskComponents(roiIdx, componentMasks, editName, true);

                    app.SelectedROIID = NaN;
                    app.setInteractionMode('idle');
                    app.refreshROITable();
                    app.refreshROITraces();
                    app.refreshEventPatches();
                    app.stackROITraceGraphics();
                    app.updateGUIEnabledState();

                    if nParts < 1
                        app.setStatusMessage('ROI edit cancelled: split ROI components were invalid.');
                    else
                        app.setStatusMessage(sprintf('Updated ROI "%s" as %d mask-clipped parts.', editName, nParts));
                    end
                    return
                end

                % Single connected clipped region: store as one polygon under the edited name.
                mask = componentMasks{1};
                roiType = 'polygon';
                verticesXY = app.maskToSimplifiedPolygonVertices(mask);
                roiParams = app.makePolygonROIParametersFromVertices(verticesXY);
            else
                verticesXY = app.getVerticesFromROIObjectForStorage(hEdit, roiType);
                roiParams = app.getROIParametersFromObjectForStorage(hEdit, roiType);
            end

            verticesXY = app.cleanROIVertices(verticesXY);
            if size(verticesXY, 1) < 3
                app.setStatusMessage('ROI edit cancelled: invalid ROI geometry.');
                app.cancelROIEditingByIndex(roiIdx, true);
                return
            end

            pgon = polyshape(verticesXY(:, 1), verticesXY(:, 2), 'Simplify', true);
            if isempty(pgon.Vertices)
                app.setStatusMessage('ROI edit cancelled: invalid polyshape.');
                app.cancelROIEditingByIndex(roiIdx, true);
                return
            end

            app.ROIList(roiIdx).name = editName;
            app.ROIList(roiIdx).type = roiType;
            app.ROIList(roiIdx).geometry.polyshape = pgon;
            app.ROIList(roiIdx).geometry.verticesXY_px = verticesXY;
            app.ROIList(roiIdx).geometry.ROIType = roiType;
            app.ROIList(roiIdx).geometry.ROIParameters = roiParams;
            app.ROIList(roiIdx).mask = mask;
            app.ROIList(roiIdx).stats = app.computeROIStatsFromMask(mask);
            app.ROIList(roiIdx).modifiedOn = datetime('now');

            app.cleanupROIEditRuntimeByIndex(roiIdx);
            app.ROIList(roiIdx).runtime.ROIHandle = app.createStaticROIOverlayFromROI(app.ROIList(roiIdx));
            app.setROISelectedStateByIndex(roiIdx, false);
            app.setSelectedROI(NaN);

            app.setInteractionMode('idle');
            app.refreshROITable();
            app.refreshROITraces();
            app.refreshEventPatches();
            app.stackROITraceGraphics();
            app.updateGUIEnabledState();

            if bWasClipped
                app.setStatusMessage(sprintf('Updated ROI "%s" and clipped it to the logical mask.', app.ROIList(roiIdx).name));
            else
                app.setStatusMessage(sprintf('Updated ROI "%s".', app.ROIList(roiIdx).name));
            end

        end

        function cancelROIEditingByIndex(app, roiIdx, bSilent)
            %CANCELROIEDITINGBYINDEX Cancel ROI edit and restore static overlay.

            if nargin < 3
                bSilent = false;
            end

            if isempty(roiIdx) || roiIdx < 1 || roiIdx > numel(app.ROIList)
                return
            end

            previousMode = 'idle';
            if isfield(app.ROIList(roiIdx), 'runtime') && ...
                    isfield(app.ROIList(roiIdx).runtime, 'editPreviousMode') && ...
                    ~isempty(app.ROIList(roiIdx).runtime.editPreviousMode)
                previousMode = app.ROIList(roiIdx).runtime.editPreviousMode;
            end

            app.cleanupROIEditRuntimeByIndex(roiIdx);

            if ~isfield(app.ROIList(roiIdx).runtime, 'ROIHandle') || ...
                    ~app.isUsableGraphicsHandle(app.ROIList(roiIdx).runtime.ROIHandle)
                app.ROIList(roiIdx).runtime.ROIHandle = app.createStaticROIOverlayFromROI(app.ROIList(roiIdx));
            end

            app.setROISelectedStateByIndex(roiIdx, false);
            app.setSelectedROI(NaN);

            try
                app.setInteractionMode(previousMode);
            catch
                app.setInteractionMode('idle');
            end

            app.refreshROITable();
            app.updateGUIEnabledState();

            if ~bSilent
                app.setStatusMessage(sprintf('ROI edit cancelled for "%s".', app.ROIList(roiIdx).name));
            end

        end

        function cleanupROIEditRuntimeByIndex(app, roiIdx)
            %CLEANUPROIEDITRUNTIMEBYINDEX Delete temporary edit ROI and listeners.

            if isempty(app.ROIList) || roiIdx < 1 || roiIdx > numel(app.ROIList)
                return
            end

            if ~isfield(app.ROIList(roiIdx), 'runtime')
                return
            end

            runtime = app.ROIList(roiIdx).runtime;
            try
                if isfield(runtime, 'editContextMenu') && ...
                        ~isempty(runtime.editContextMenu) && ...
                        isvalid(runtime.editContextMenu)
                    delete(runtime.editContextMenu);
                end
            catch
            end

            try
                if isfield(runtime, 'editPreviousKeyFcn') && isvalid(app.UIFigure)
                    app.UIFigure.WindowKeyPressFcn = runtime.editPreviousKeyFcn;
                end
            catch
            end

            if isfield(runtime, 'editListeners') && ~isempty(runtime.editListeners)
                listeners = runtime.editListeners;
                for iListener = 1:numel(listeners)
                    try
                        if ~isempty(listeners{iListener}) && isvalid(listeners{iListener})
                            delete(listeners{iListener});
                        end
                    catch
                    end
                end
            end

            if isfield(runtime, 'editHandle') && app.isUsableGraphicsHandle(runtime.editHandle)
                try
                    delete(runtime.editHandle);
                catch
                end
            end

            cleanupFields = { ...
                'editHandle', ...
                'editListeners', ...
                'editName', ...
                'editOriginalName', ...
                'editPreviousMode', ...
                'editPreviousKeyFcn'};

            for iField = 1:numel(cleanupFields)
                fieldName = cleanupFields{iField};
                if isfield(app.ROIList(roiIdx).runtime, fieldName)
                    app.ROIList(roiIdx).runtime = rmfield(app.ROIList(roiIdx).runtime, fieldName);
                end
            end

            app.ROIList(roiIdx).runtime.editHandle = gobjects(1);
            app.ROIList(roiIdx).runtime.editContextMenu = matlab.ui.container.ContextMenu.empty;
        end

        function onEditableROIMoving(app, src, evt, roiID)
            %ONEDITABLEROIMOVING Update temporary label while editing ROI.

            roiIdx = app.getROIIndexByID(roiID);
            if isempty(roiIdx)
                return
            end

            editName = app.getCurrentROIEditName(roiIdx);
            currentSize = app.getUnrotatedROIObjectSize(src);

            if isempty(currentSize)
                labelText = editName;
            else
                [scaleFactor, unitText] = app.getROIAxisUnitScale();
                displaySize = currentSize .* scaleFactor;
                labelText = sprintf('%s | %.3g x %.3g %s', ...
                    editName, displaySize(1), displaySize(2), unitText);
            end

            app.setROIObjectPropertyIfAvailable(src, 'Label', labelText);

            if nargin >= 3 && isprop(evt, 'CurrentPosition')
                % Keep callback signature explicit. Position clamping is handled by
                % the ROI DrawingArea.
            end

        end

        function onEditableROIMoved(app, src, evt, roiID)
            %ONEDITABLEROIMOVED Restore label after moving/resizing an edited ROI.

            roiIdx = app.getROIIndexByID(roiID);
            if isempty(roiIdx)
                return
            end

            app.setROIObjectPropertyIfAvailable(src, 'Label', app.getCurrentROIEditName(roiIdx));

            if nargin >= 3 && isprop(evt, 'CurrentPosition')
                % Position is accepted as delivered by the images.roi object.
            end

        end

        function onEditableROIClicked(app, src, evt, roiID)
            %ONEDITABLEROICLICKED Handle rename and confirmation during ROI edit.

            roiIdx = app.getROIIndexByID(roiID);
            if isempty(roiIdx)
                return
            end

            if app.isROILabelClickEvent(evt)
                app.renameROIWhileEditingByIndex(roiIdx, src);
                return
            end

            if isprop(evt, 'SelectionType') && strcmpi(evt.SelectionType, 'double')
                app.commitROIEditingByIndex(roiIdx);
            end

        end

        function onEditableROIDeleting(app, ~, ~, roiID)
            %ONEDITABLEROIDELETING Cancel edit when temporary ROI is deleted by user.

            roiIdx = app.getROIIndexByID(roiID);
            if isempty(roiIdx)
                return
            end

            app.cancelROIEditingByIndex(roiIdx, false);

        end

        function onActiveROIEditKeyPress(app, ~, evt, roiID)
            %ONACTIVEROIEDITKEYPRESS Handle keyboard shortcuts for active ROI edit.

            roiIdx = app.getROIIndexByID(roiID);
            if isempty(roiIdx)
                return
            end

            switch lower(evt.Key)
                case {'return', 'enter'}
                    app.commitROIEditingByIndex(roiIdx);

                case 'escape'
                    app.cancelROIEditingByIndex(roiIdx, false);

                case 'n'
                    if isfield(app.ROIList(roiIdx).runtime, 'editHandle') && ...
                            app.isUsableGraphicsHandle(app.ROIList(roiIdx).runtime.editHandle)
                        app.renameROIWhileEditingByIndex( ...
                            roiIdx, app.ROIList(roiIdx).runtime.editHandle);
                    end
            end

        end

        function renameROIWhileEditingByIndex(app, roiIdx, hEdit)
            %RENAMEROIWHILEEDITINGBYINDEX Rename active ROI edit label.

            if isempty(roiIdx) || roiIdx < 1 || roiIdx > numel(app.ROIList)
                return
            end

            currentName = app.getCurrentROIEditName(roiIdx);
            newName = app.promptForROINameDialog(currentName);

            if isempty(newName)
                app.setROIObjectPropertyIfAvailable(hEdit, 'Label', currentName);
                return
            end

            newName = app.makeUniqueROIName(newName, app.ROIList(roiIdx).ID);
            app.ROIList(roiIdx).runtime.editName = newName;
            app.setROIObjectPropertyIfAvailable(hEdit, 'Label', newName);
            app.setROIObjectPropertyIfAvailable(hEdit, 'LabelVisible', 'on');

        end

        function editName = getCurrentROIEditName(app, roiIdx)
            %GETCURRENTROIEDITNAME Return temporary edit name for one ROI.

            editName = app.ROIList(roiIdx).name;

            if isfield(app.ROIList(roiIdx), 'runtime') && ...
                    isfield(app.ROIList(roiIdx).runtime, 'editName') && ...
                    ~isempty(app.ROIList(roiIdx).runtime.editName)
                editName = char(string(app.ROIList(roiIdx).runtime.editName));
            end

        end

        function hStatic = createStaticROIOverlayFromROI(app, ROI)
            %CREATESTATICROIOVERLAYFROMROI Draw clickable persistent/static ROI overlay.
            %
            %   The static overlay is the stored/display representation. Clicking it
            %   enters ROI edit mode by replacing this overlay with an images.roi object.

            hStatic = gobjects(1);

            roiColor = [1 0 0];
            if isfield(ROI, 'color') && numel(ROI.color) == 3
                roiColor = min(max(double(ROI.color(:).'), 0), 1);
            end

            roiID = [];
            if isfield(ROI, 'ID') && ~isempty(ROI.ID)
                roiID = ROI.ID;
            end

            try
                hStatic = plot(app.ImageAxes, ROI.geometry.polyshape);
                app.configureStaticROIOverlayHandle(hStatic, roiColor, roiID);
            catch
                if ~isfield(ROI, 'geometry') || ~isfield(ROI.geometry, 'verticesXY_px')
                    return
                end

                verts = ROI.geometry.verticesXY_px;
                hStatic = patch(app.ImageAxes, ...
                    'XData', verts(:, 1), ...
                    'YData', verts(:, 2), ...
                    'FaceColor', roiColor, ...
                    'FaceAlpha', 0.2, ...
                    'EdgeColor', roiColor, ...
                    'LineStyle', '-', ...
                    'LineWidth', 1.5, ...
                    'HitTest', 'on', ...
                    'PickableParts', 'all', ...
                    'HandleVisibility', 'off');

                app.configureStaticROIOverlayHandle(hStatic, roiColor, roiID);
            end

        end

        function configureStaticROIOverlayHandle(app, hStatic, roiColor, roiID) %#ok<INUSD>
            %CONFIGURESTATICROIOVERLAYHANDLE Apply visual settings to static ROI overlay.
            %
            %   Static ROI overlays are intentionally non-interactive. Click handling is
            %   routed through ImageAxes/ImageHandle so single-click can inspect pixels
            %   and double-click inside an ROI can enter single or group edit mode.

            if ~app.isUsableGraphicsHandle(hStatic)
                return
            end

            try
                if isprop(hStatic, 'FaceColor')
                    hStatic.FaceColor = roiColor;
                end
                if isprop(hStatic, 'FaceAlpha')
                    hStatic.FaceAlpha = 0.2;
                end
                if isprop(hStatic, 'EdgeColor')
                    hStatic.EdgeColor = roiColor;
                end
                if isprop(hStatic, 'LineStyle')
                    hStatic.LineStyle = '-';
                end
                if isprop(hStatic, 'LineWidth')
                    hStatic.LineWidth = 1.5;
                end
                if isprop(hStatic, 'HitTest')
                    hStatic.HitTest = 'off';
                end
                if isprop(hStatic, 'PickableParts')
                    hStatic.PickableParts = 'none';
                end
                if isprop(hStatic, 'HandleVisibility')
                    hStatic.HandleVisibility = 'off';
                end
                if isprop(hStatic, 'ButtonDownFcn')
                    hStatic.ButtonDownFcn = [];
                end
            catch
            end

        end

        function mask = createMaskFromROIObject(app, hROI)
            %CREATEMASKFROMROIOBJECT Create logical image-sized ROI mask.

            try
                mask = createMask(hROI, app.ImageHandle);
            catch
                mask = createMask(hROI);
            end

            mask = logical(mask);

        end

        function verticesXY = getVerticesFromROIObjectForStorage(app, hROI, roiType)
            %GETVERTICESFROMROIOBJECTFORSTORAGE Extract/synthesize ROI vertices.

            verticesXY = [];
            roiType = lower(char(string(roiType)));

            switch roiType
                case 'rectangle'
                    if isprop(hROI, 'Position') && ~isempty(hROI.Position)
                        pos = double(hROI.Position);
                        verticesXY = [ ...
                            pos(1),          pos(2); ...
                            pos(1) + pos(3), pos(2); ...
                            pos(1) + pos(3), pos(2) + pos(4); ...
                            pos(1),          pos(2) + pos(4)];

                        angleDeg = app.getNumericROIObjectProperty(hROI, 'RotationAngle', 0);
                        verticesXY = app.rotateROIVerticesForStorage(verticesXY, ...
                            [pos(1) + pos(3)/2, pos(2) + pos(4)/2], angleDeg);
                    end

                case 'ellipse'
                    theta = linspace(0, 2*pi, 129).';
                    theta(end) = [];

                    if isprop(hROI, 'Center') && isprop(hROI, 'SemiAxes') && ...
                            ~isempty(hROI.Center) && ~isempty(hROI.SemiAxes)
                        centerXY = double(hROI.Center(:).');
                        semiAxes = double(hROI.SemiAxes(:).');
                    elseif isprop(hROI, 'Position') && ~isempty(hROI.Position)
                        pos = double(hROI.Position);
                        centerXY = [pos(1) + pos(3)/2, pos(2) + pos(4)/2];
                        semiAxes = [pos(3)/2, pos(4)/2];
                    else
                        verticesXY = [];
                        return
                    end

                    verticesXY = [semiAxes(1) .* cos(theta), semiAxes(2) .* sin(theta)];
                    angleDeg = app.getNumericROIObjectProperty(hROI, 'RotationAngle', 0);
                    verticesXY = app.rotateROIVerticesForStorage(verticesXY + centerXY, centerXY, angleDeg);

                case 'polygon'
                    if isprop(hROI, 'Position') && ~isempty(hROI.Position)
                        verticesXY = double(hROI.Position);
                    elseif isprop(hROI, 'Vertices') && ~isempty(hROI.Vertices)
                        verticesXY = double(hROI.Vertices);
                    end
            end

            verticesXY = app.cleanROIVertices(verticesXY);

        end

        function params = getROIParametersFromObjectForStorage(app, hROI, roiType) %#ok<INUSL>
            %GETROIPARAMETERSFROMOBJECTFORSTORAGE Store ROI reconstruction parameters.

            params = struct();
            params.ROIType = roiType;

            propList = { ...
                'Position', ...
                'Center', ...
                'SemiAxes', ...
                'RotationAngle', ...
                'Vertices'};

            for iProp = 1:numel(propList)
                propName = propList{iProp};

                try
                    if isprop(hROI, propName) && ~isempty(hROI.(propName))
                        params.(propName) = hROI.(propName);
                    end
                catch
                end
            end

        end

        function verticesOut = rotateROIVerticesForStorage(app, verticesIn, centerXY, angleDeg) %#ok<INUSL>
            %ROTATEROIVERTICESFORSTORAGE Rotate ROI vertices for polyshape storage.

            if isempty(verticesIn) || angleDeg == 0
                verticesOut = verticesIn;
                return
            end

            % images.roi RotationAngle is expressed in image coordinates, where Y
            % increases downward. Invert the sign so the stored polyshape matches the
            % displayed ROI geometry in ImageAxes.
            angleRad = -deg2rad(angleDeg);
            R = [cos(angleRad), -sin(angleRad); sin(angleRad), cos(angleRad)];

            verticesOut = (verticesIn - centerXY) * R.' + centerXY;

        end

        function verticesXY = cleanROIVertices(app, verticesXY) %#ok<INUSL>
            %CLEANROIVERTICES Remove invalid and duplicate consecutive vertices.

            if isempty(verticesXY) || size(verticesXY, 2) ~= 2
                verticesXY = zeros(0, 2);
                return
            end

            verticesXY = double(verticesXY);
            verticesXY = verticesXY(all(isfinite(verticesXY), 2), :);

            if size(verticesXY, 1) > 1
                duplicate = [false; all(abs(diff(verticesXY, 1, 1)) < eps, 2)];
                verticesXY(duplicate, :) = [];
            end

        end

        function value = getNumericROIObjectProperty(app, hROI, propName, defaultValue) %#ok<INUSL>
            %GETNUMERICROIOBJECTPROPERTY Read scalar numeric ROI property safely.

            value = defaultValue;

            try
                if isprop(hROI, propName) && ~isempty(hROI.(propName))
                    candidate = double(hROI.(propName));
                    if isscalar(candidate) && isfinite(candidate)
                        value = candidate;
                    end
                end
            catch
                value = defaultValue;
            end

        end

        function sizeXY = getUnrotatedROIObjectSize(app, hROI) %#ok<INUSL>
            %GETUNROTATEDROIOBJECTSIZE Estimate unrotated ROI size in pixels.

            sizeXY = [];

            try
                if isprop(hROI, 'SemiAxes') && ~isempty(hROI.SemiAxes)
                    semiAxes = double(hROI.SemiAxes(:).');
                    if numel(semiAxes) >= 2
                        sizeXY = 2 .* semiAxes(1:2);
                        return
                    end
                end
            catch
            end

            try
                if isprop(hROI, 'Position') && ~isempty(hROI.Position)
                    pos = double(hROI.Position);
                    if numel(pos) >= 4
                        sizeXY = abs(pos(3:4));
                        return
                    end
                end
            catch
            end

            try
                if isprop(hROI, 'Vertices') && ~isempty(hROI.Vertices)
                    verts = double(hROI.Vertices);
                elseif isprop(hROI, 'Position') && size(hROI.Position, 2) == 2
                    verts = double(hROI.Position);
                else
                    verts = [];
                end

                if ~isempty(verts)
                    sizeXY = [range(verts(:, 1)), range(verts(:, 2))];
                end
            catch
            end

        end

        function [scaleFactor, unitText] = getROIAxisUnitScale(app)
            %GETROIAXISUNITSCALE Return size conversion for ROI dimension label.

            scaleFactor = 1;
            unitText = 'px';

            if isempty(app.DataParams) || ~isfield(app.DataParams, 'view')
                return
            end

            viewInfo = app.DataParams.view;

            if ~isfield(viewInfo, 'pixelSize_px_per_mm') || isempty(viewInfo.pixelSize_px_per_mm)
                return
            end

            pxPerMm = double(viewInfo.pixelSize_px_per_mm);

            if isscalar(pxPerMm) && isfinite(pxPerMm) && pxPerMm > 0
                scaleFactor = 1 ./ pxPerMm;
                unitText = 'mm';
            end

        end

        function setROIObjectPropertyIfAvailable(app, hROI, propName, propValue) %#ok<INUSL>
            %SETROIOBJECTPROPERTYIFAVAILABLE Set ROI object property if it exists.

            try
                if ~isempty(hROI) && isvalid(hROI) && isprop(hROI, propName)
                    hROI.(propName) = propValue;
                end
            catch
            end

        end

        function tf = isROILabelClickEvent(app, evt) %#ok<INUSL>
            %ISROILABELCLICKEVENT Best-effort detection of ROI label click.

            tf = false;

            candidateParts = {'SelectedPart', 'HitPart', 'ClickedPart', 'Part'};

            for iPart = 1:numel(candidateParts)
                propName = candidateParts{iPart};

                try
                    if isprop(evt, propName)
                        value = lower(char(string(evt.(propName))));
                        tf = contains(value, 'label');
                        if tf
                            return
                        end
                    end
                catch
                end
            end

            try
                if isprop(evt, 'HitObject') && ~isempty(evt.HitObject)
                    h = evt.HitObject;
                    if isprop(h, 'Type')
                        tf = contains(lower(char(string(h.Type))), 'text');
                    end
                end
            catch
            end

        end

        function newName = promptForROINameDialog(app, currentName)
            %PROMPTFORROINAMEDIALOG Modal dialog for renaming an ROI.

            newName = '';

            dlg = uifigure( ...
                'Name', 'Rename ROI', ...
                'WindowStyle', 'modal', ...
                'Position', [100 100 320 120], ...
                'Visible', 'off', ...
                'CloseRequestFcn', @onCancel);

            grid = uigridlayout(dlg);
            grid.RowHeight = {28, '1x', 32};
            grid.ColumnWidth = {70, '1x'};
            grid.Padding = [12 12 12 12];

            nameLabel = uilabel(grid);
            nameLabel.Text = 'Name';
            nameLabel.Layout.Row = 1;
            nameLabel.Layout.Column = 1;

            nameField = uieditfield(grid, 'text');
            nameField.Value = currentName;
            nameField.Layout.Row = 1;
            nameField.Layout.Column = 2;

            statusLabel = uilabel(grid);
            statusLabel.Text = '';
            statusLabel.FontColor = [0.65 0 0];
            statusLabel.Layout.Row = 2;
            statusLabel.Layout.Column = [1 2];

            buttonGrid = uigridlayout(grid);
            buttonGrid.RowHeight = {'1x'};
            buttonGrid.ColumnWidth = {'1x', '1x'};
            buttonGrid.Padding = [0 0 0 0];
            buttonGrid.Layout.Row = 3;
            buttonGrid.Layout.Column = [1 2];

            okButton = uibutton(buttonGrid, 'push');
            okButton.Text = 'OK';
            okButton.Layout.Row = 1;
            okButton.Layout.Column = 1;
            okButton.ButtonPushedFcn = @onOK;

            cancelButton = uibutton(buttonGrid, 'push');
            cancelButton.Text = 'Cancel';
            cancelButton.Layout.Row = 1;
            cancelButton.Layout.Column = 2;
            cancelButton.ButtonPushedFcn = @onCancel;

            if exist('placeAppInsideCaller', 'file') == 2
                placeAppInsideCaller(app, dlg, 'center');
            end

            dlg.Visible = 'on';
            drawnow
            uiwait(dlg);

            if isvalid(dlg)
                delete(dlg);
            end

            function onOK(~, ~)
                requestedName = strtrim(char(string(nameField.Value)));

                if isempty(requestedName)
                    statusLabel.Text = 'ROI name cannot be empty.';
                    return
                end

                newName = requestedName;
                uiresume(dlg);
            end

            function onCancel(~, ~)
                newName = '';
                uiresume(dlg);
            end

        end

        function updateROIStatsForCurrentFrame(app)
            %UPDATEROISTATSFORCURRENTFRAME Refresh ROI statistics from displayed frame.
            %
            %   This updates only frame-dependent spatial statistics and the ROI table.
            %   ROI geometry, masks, and temporal traces are left unchanged. During movie
            %   playback this method should be skipped and called once after playback
            %   stops to avoid expensive table updates on every timer tick.

            if ~app.hasData() || isempty(app.ROIList)
                return
            end

            frame = [];

            if app.isUsableGraphicsHandle(app.ImageHandle) && isprop(app.ImageHandle, 'CData')
                frame = app.ImageHandle.CData;
            end

            if isempty(frame) || ~isnumeric(frame)
                frame = app.getCurrentFrame();
            end

            frame = squeeze(frame);

            if isempty(frame) || ~ismatrix(frame)
                return
            end

            frameSize = size(frame);

            for iROI = 1:numel(app.ROIList)
                if ~isfield(app.ROIList(iROI), 'mask') || isempty(app.ROIList(iROI).mask)
                    continue
                end

                mask = logical(app.ROIList(iROI).mask);

                if ~isequal(size(mask), frameSize)
                    continue
                end

                app.ROIList(iROI).stats = app.computeROIStatsFromMask(mask, frame);
                app.ROIList(iROI).modifiedOn = datetime('now');
            end

            app.refreshROITable();

        end

        function stats = computeROIStatsFromMask(app, mask, frame)
            %COMPUTEROISTATSFROMMASK Compute ROI geometry and current-frame statistics.
            %
            %   stats = computeROIStatsFromMask(app, mask) computes ROI statistics from
            %   the frame currently shown in ImageAxes.
            %
            %   stats = computeROIStatsFromMask(app, mask, frame) computes the same
            %   statistics from the supplied 2-D frame. This form avoids rereading data
            %   when a frame was just displayed.

            if nargin < 3 || isempty(frame)
                frame = [];
                if app.isUsableGraphicsHandle(app.ImageHandle) && isprop(app.ImageHandle, 'CData')
                    frame = app.ImageHandle.CData;
                end
            end

            mask = logical(mask);
            stats = struct();

            stats.computedOn = datetime('now');
            stats.NPixels = double(nnz(mask));
            stats.areaPx2 = double(nnz(mask));
            stats.areaMM2 = [];

            stats.centroidXY_px = [];
            stats.centroidXY_mm = [];

            stats.distanceFromOrigin_px = [];
            stats.distanceFromOrigin_mm = [];

            stats.spatialMean = [];
            stats.spatialStd = [];
            stats.spatialMedian = [];
            stats.spatialMin = [];
            stats.spatialMax = [];

            if stats.NPixels == 0
                return
            end

            [yIdx, xIdx] = find(mask);
            centroidXY = [mean(xIdx), mean(yIdx)];
            stats.centroidXY_px = centroidXY;

            originXY = [1 1];
            pxPerMm = [];

            if ~isempty(app.DataParams) && isfield(app.DataParams, 'view')
                viewInfo = app.DataParams.view;

                if isfield(viewInfo, 'origin_xy_px') && numel(viewInfo.origin_xy_px) == 2
                    originCandidate = double(viewInfo.origin_xy_px(:).');
                    if all(isfinite(originCandidate))
                        originXY = originCandidate;
                    end
                end

                if isfield(viewInfo, 'pixelSize_px_per_mm') && ...
                        ~isempty(viewInfo.pixelSize_px_per_mm)
                    pxCandidate = double(viewInfo.pixelSize_px_per_mm);
                    if isscalar(pxCandidate) && isfinite(pxCandidate) && pxCandidate > 0
                        pxPerMm = pxCandidate;
                    end
                end
            end

            deltaXY_px = centroidXY - originXY;
            stats.distanceFromOrigin_px = hypot(deltaXY_px(1), deltaXY_px(2));

            if ~isempty(pxPerMm)
                stats.centroidXY_mm = deltaXY_px ./ pxPerMm;
                stats.distanceFromOrigin_mm = stats.distanceFromOrigin_px ./ pxPerMm;
                stats.areaMM2 = stats.areaPx2 ./ (pxPerMm .^ 2);
            end

            frame = squeeze(frame);

            if isempty(frame) || ~isnumeric(frame) || ~isequal(size(frame), size(mask))
                return
            end

            roiPixels = double(frame(mask));
            roiPixels = roiPixels(isfinite(roiPixels));

            if isempty(roiPixels)
                stats.spatialMean = NaN;
                stats.spatialStd = NaN;
                stats.spatialMedian = NaN;
                stats.spatialMin = NaN;
                stats.spatialMax = NaN;
                return
            end

            stats.spatialMean = mean(roiPixels);
            stats.spatialStd = std(roiPixels);
            stats.spatialMedian = median(roiPixels);
            stats.spatialMin = min(roiPixels);
            stats.spatialMax = max(roiPixels);

        end

        function refreshROITraces(app)
            %REFRESHROITRACES Compute and display temporal traces for all ROIs.
            %
            %   ROI traces are runtime-only display data. This method computes traces
            %   for every ROI in app.ROIList, including hidden ROIs, so toggling ROI
            %   visibility later is only a graphics operation and does not trigger
            %   expensive recomputation.

            if ~app.hasData() || isempty(app.ROIList)
                return
            end

            caps = app.getDataCapabilities();

            if ~caps.hasTemporalDimension
                for iROI = 1:numel(app.ROIList)
                    app.hideROITraceGraphicsByIndex(iROI);
                    app.ROIList(iROI).runtime.trace.Status = 'no temporal dimension';
                end
                return
            end

            [roiMasks, roiIdxList] = app.collectROIMasksForTrace();

            if isempty(roiIdxList)
                return
            end

            try
                [traceMatrix, traceMode, readMode] = app.computeROITraceMatrixForCurrentState(roiMasks);

                app.updateROITraceGraphicsFromMatrix( ...
                    roiIdxList, ...
                    traceMatrix, ...
                    traceMode, ...
                    readMode);

            catch ME
                for k = 1:numel(roiIdxList)
                    roiIdx = roiIdxList(k);
                    app.hideROITraceGraphicsByIndex(roiIdx);

                    if isfield(app.ROIList(roiIdx), 'runtime') && ...
                            isfield(app.ROIList(roiIdx).runtime, 'trace')
                        app.ROIList(roiIdx).runtime.trace.Status = 'error';
                    end
                end

                app.setStatusMessage(sprintf('ROI trace update failed: %s', ME.message));
            end

        end

        function [roiMasks, roiIdxList] = collectROIMasksForTrace(app)
            %COLLECTROIMASKSFORTRACE Return valid ROI masks matching current image size.

            roiMasks = {};
            roiIdxList = [];

            if isempty(app.ROIList)
                return
            end

            sz = app.getDataSize();
            expectedSize = double(sz(1:2));

            for iROI = 1:numel(app.ROIList)
                if ~isfield(app.ROIList(iROI), 'mask') || isempty(app.ROIList(iROI).mask)
                    app.hideROITraceGraphicsByIndex(iROI);
                    app.ROIList(iROI).runtime.trace.Status = 'missing mask';
                    continue
                end

                mask = logical(app.ROIList(iROI).mask);

                if ~isequal(double(size(mask)), expectedSize)
                    app.hideROITraceGraphicsByIndex(iROI);
                    app.ROIList(iROI).runtime.trace.Status = 'mask size mismatch';
                    continue
                end

                if ~any(mask(:))
                    app.hideROITraceGraphicsByIndex(iROI);
                    app.ROIList(iROI).runtime.trace.Status = 'empty mask';
                    continue
                end

                roiMasks{end+1} = mask; %#ok<AGROW>
                roiIdxList(end+1) = iROI; %#ok<AGROW>
            end

        end

        function [traceMatrix, traceMode, readMode] = computeROITraceMatrixForCurrentState(app, roiMasks)
            %COMPUTEROITRACEMATRIXFORCURRENTSTATE Compute ROI traces for the current view.
            %
            %   Output shape:
            %       Normal mode or single event repetition:
            %           traceMatrix = nROI x nFrames
            %
            %       Event average mode:
            %           traceMatrix = nROI x nTrials x nFrames

            traceMatrix = [];
            traceMode = 'normal';
            readMode = 'generic';

            sourceType = app.getSourceType();

            switch lower(sourceType)
                case 'dat'
                    [traceMatrix, traceMode, readMode] = app.computeDatROITraceMatrix(roiMasks);

                otherwise
                    [traceMatrix, traceMode, readMode] = app.computeGenericROITraceMatrix(roiMasks);
            end

        end

        function [traceMatrix, traceMode, readMode] = computeDatROITraceMatrix(app, roiMasks)
            %COMPUTEDATROITRACEMATRIX Compute ROI traces using DatImageSource.

            traceMatrix = [];
            readMode = 'generic';

            if ~ismethod(app.DataSource, 'getROIMeanTraceMatrix')
                [traceMatrix, traceMode, readMode] = app.computeGenericROITraceMatrix(roiMasks);
                return
            end

            if strcmpi(app.ViewMode, 'event')
                if isempty(app.EventFrameMatrix)
                    error('DataViewer:MissingEventFrameMatrix', ...
                        'Event frame matrix is empty.');
                end

                if strcmpi(app.CurrentRepetition, 'AVERAGE')
                    frameIdx = app.EventFrameMatrix;
                    traceMode = 'event_average';
                else
                    frameIdx = app.EventFrameMatrix(1, :);
                    traceMode = 'event';
                end
            else
                sourceSize = app.getSourceDataSize();
                frameIdx = 1:sourceSize(3);
                traceMode = 'normal';
            end

            [traceMatrix, readMode] = app.DataSource.getROIMeanTraceMatrix(roiMasks, frameIdx);

        end

        function [traceMatrix, traceMode, readMode] = computeGenericROITraceMatrix(app, roiMasks)
            %COMPUTEGENERICROITRACEMATRIX Compute ROI traces using DataSource.getFrame.
            %
            %   This is a correctness fallback for non-DAT sources. DAT sources should
            %   normally use DatImageSource.getROIMeanTraceMatrix for performance.

            readMode = 'generic_getFrame';

            if strcmpi(app.ViewMode, 'event') && strcmpi(app.getSourceType(), 'umt')
                eIdxList = app.getSelectedUMTEventIndices();

                if isempty(eIdxList)
                    eIdxList = 1;
                end

                sz = app.getDataSize();
                nFrames = sz(3);

                if strcmpi(app.CurrentRepetition, 'AVERAGE') && numel(eIdxList) > 1
                    traceMode = 'event_average';
                    traceMatrix = nan(numel(roiMasks), numel(eIdxList), nFrames);

                    for iE = 1:numel(eIdxList)
                        for tIdx = 1:nFrames
                            frame = app.DataSource.getFrame(tIdx, eIdxList(iE));
                            traceMatrix(:, iE, tIdx) = app.computeROIMeansFromFrame(roiMasks, frame);
                        end
                    end
                else
                    traceMode = 'event';
                    traceMatrix = nan(numel(roiMasks), nFrames);
                    eIdx = eIdxList(1);

                    for tIdx = 1:nFrames
                        frame = app.DataSource.getFrame(tIdx, eIdx);
                        traceMatrix(:, tIdx) = app.computeROIMeansFromFrame(roiMasks, frame);
                    end
                end

                return
            end

            sourceSize = app.getSourceDataSize();
            nFrames = sourceSize(3);

            traceMode = 'normal';
            traceMatrix = nan(numel(roiMasks), nFrames);

            for tIdx = 1:nFrames
                frame = app.DataSource.getFrame(tIdx);
                traceMatrix(:, tIdx) = app.computeROIMeansFromFrame(roiMasks, frame);
            end

        end

        function roiMeans = computeROIMeansFromFrame(app, roiMasks, frame) %#ok<INUSL>
            %COMPUTEROIMEANSFROMFRAME Compute all ROI spatial means from one image frame.

            frame = squeeze(frame);
            roiMeans = nan(numel(roiMasks), 1);

            for iROI = 1:numel(roiMasks)
                mask = roiMasks{iROI};

                if ~isequal(size(mask), size(frame))
                    continue
                end

                pix = double(frame(mask));
                pix = pix(isfinite(pix));

                if isempty(pix)
                    roiMeans(iROI) = NaN;
                else
                    roiMeans(iROI) = mean(pix);
                end
            end

        end

        function updateROITraceGraphicsFromMatrix(app, roiIdxList, traceMatrix, traceMode, readMode)
            %UPDATEROITRACEGRAPHICSFROMMATRIX Store and plot ROI traces.

            if isempty(roiIdxList) || isempty(traceMatrix)
                return
            end

            switch lower(traceMode)
                case 'event_average'
                    nSamples = size(traceMatrix, 3);
                otherwise
                    nSamples = size(traceMatrix, 2);
            end

            xData = app.getTraceTimeVector(nSamples);

            for k = 1:numel(roiIdxList)
                roiIdx = roiIdxList(k);

                switch lower(traceMode)
                    case 'event_average'
                        trialData = squeeze(traceMatrix(k, :, :)); % nTrials x nFrames

                        if isvector(trialData)
                            trialData = reshape(trialData, 1, []);
                        end

                        [meanTrace, semTrace] = app.computeMeanAndSEMFromTrials(trialData);

                    otherwise
                        meanTrace = squeeze(traceMatrix(k, :));
                        meanTrace = meanTrace(:);
                        semTrace = [];
                end

                app.ROIList(roiIdx).runtime.trace.XData = xData(:);
                app.ROIList(roiIdx).runtime.trace.Mean = meanTrace(:);
                app.ROIList(roiIdx).runtime.trace.SEM = semTrace(:);
                app.ROIList(roiIdx).runtime.trace.Mode = traceMode;
                app.ROIList(roiIdx).runtime.trace.Status = sprintf('ok:%s', readMode);

                app.updateSingleROITraceGraphics(roiIdx);
            end

            app.stackROITraceGraphics();

        end

        function [meanTrace, semTrace] = computeMeanAndSEMFromTrials(app, trialData) %#ok<INUSL>
            %COMPUTEMEANANDSEMFROMTRIALS Compute mean/SEM across trial rows.

            trialData = double(trialData);

            nValid = sum(isfinite(trialData), 1);

            meanTrace = nan(1, size(trialData, 2));
            semTrace = nan(1, size(trialData, 2));

            for iFrame = 1:size(trialData, 2)
                vals = trialData(:, iFrame);
                vals = vals(isfinite(vals));

                if isempty(vals)
                    continue
                end

                meanTrace(iFrame) = mean(vals);

                if numel(vals) > 1
                    semTrace(iFrame) = std(vals, 0) ./ sqrt(numel(vals));
                else
                    semTrace(iFrame) = 0;
                end
            end

            semTrace(nValid == 0) = NaN;

            meanTrace = meanTrace(:);
            semTrace = semTrace(:);

        end

        function updateSingleROITraceGraphics(app, roiIdx)
            %UPDATESINGLEROITRACEGRAPHICS Create/update one ROI trace and SEM patch.

            if isempty(roiIdx) || roiIdx < 1 || roiIdx > numel(app.ROIList)
                return
            end

            ROI = app.ROIList(roiIdx);

            if ~isfield(ROI, 'runtime') || ~isfield(ROI.runtime, 'trace')
                return
            end

            trace = ROI.runtime.trace;

            xData = trace.XData(:);
            meanTrace = trace.Mean(:);
            semTrace = trace.SEM(:);

            hasTrace = ~isempty(xData) && ~isempty(meanTrace) && ...
                numel(xData) == numel(meanTrace) && ...
                any(isfinite(xData)) && any(isfinite(meanTrace));

            if ~hasTrace
                app.hideROITraceGraphicsByIndex(roiIdx);
                return
            end

            visibleState = app.getROIVisibleState(roiIdx);

            if ~isfield(app.ROIList(roiIdx).runtime, 'traceLine') || ...
                    ~app.isUsableGraphicsHandle(app.ROIList(roiIdx).runtime.traceLine)
                app.ROIList(roiIdx).runtime.traceLine = plot(app.PlotAxes, nan, nan, ...
                    'LineWidth', 1.25, ...
                    'Color', ROI.color, ...
                    'HandleVisibility', 'off', ...
                    'HitTest', 'off', ...
                    'PickableParts', 'none');
            end

            hLine = app.ROIList(roiIdx).runtime.traceLine;

            set(hLine, ...
                'XData', xData, ...
                'YData', meanTrace, ...
                'Color', ROI.color, ...
                'Visible', visibleState);

            hasSEM = ~isempty(semTrace) && numel(semTrace) == numel(meanTrace) && ...
                any(isfinite(semTrace));

            if hasSEM
                upperTrace = meanTrace(:).' + semTrace(:).';
                lowerTrace = meanTrace(:).' - semTrace(:).';

                patchX = [xData(:).', fliplr(xData(:).')];
                patchY = [upperTrace, fliplr(lowerTrace)];

                if ~isfield(app.ROIList(roiIdx).runtime, 'tracePatch') || ...
                        ~app.isUsableGraphicsHandle(app.ROIList(roiIdx).runtime.tracePatch)
                    app.ROIList(roiIdx).runtime.tracePatch = patch(app.PlotAxes, ...
                        nan, nan, ROI.color, ...
                        'FaceAlpha', 0.16, ...
                        'EdgeColor', 'none', ...
                        'HandleVisibility', 'off', ...
                        'HitTest', 'off', ...
                        'PickableParts', 'none');
                end

                hPatch = app.ROIList(roiIdx).runtime.tracePatch;

                set(hPatch, ...
                    'XData', patchX, ...
                    'YData', patchY, ...
                    'FaceColor', ROI.color, ...
                    'EdgeColor', 'none', ...
                    'Visible', visibleState);
            else
                if isfield(ROI.runtime, 'tracePatch') && ...
                        app.isUsableGraphicsHandle(ROI.runtime.tracePatch) && ...
                        isprop(ROI.runtime.tracePatch, 'Visible')
                    ROI.runtime.tracePatch.Visible = 'off';
                end
            end

        end

        function visibleState = getROIVisibleState(app, roiIdx)
            %GETROIVISIBLESTATE Return 'on' or 'off' for ROI runtime visibility.

            visibleState = 'on';

            if roiIdx < 1 || roiIdx > numel(app.ROIList)
                visibleState = 'off';
                return
            end

            if isfield(app.ROIList(roiIdx), 'runtime') && ...
                    isfield(app.ROIList(roiIdx).runtime, 'visible') && ...
                    ~app.ROIList(roiIdx).runtime.visible
                visibleState = 'off';
            end

        end

        function hideROITraceGraphicsByIndex(app, roiIdx)
            %HIDEROITRACEGRAPHICSBYINDEX Hide one ROI temporal trace and patch.

            if isempty(app.ROIList) || roiIdx < 1 || roiIdx > numel(app.ROIList)
                return
            end

            if ~isfield(app.ROIList(roiIdx), 'runtime')
                return
            end

            runtimeFields = {'traceLine', 'tracePatch'};

            for iField = 1:numel(runtimeFields)
                fieldName = runtimeFields{iField};

                if ~isfield(app.ROIList(roiIdx).runtime, fieldName)
                    continue
                end

                h = app.ROIList(roiIdx).runtime.(fieldName);

                if app.isUsableGraphicsHandle(h) && isprop(h, 'Visible')
                    try
                        h.Visible = 'off';
                    catch
                    end
                end
            end

        end

        function deleteROIByID(app, roiID)
            %DELETEROIBYID Delete one ROI by runtime ROI ID.
            %
            %   This wrapper keeps backward compatibility with older single-ROI delete
            %   calls while routing deletion through the batch-safe implementation.

            if isempty(roiID) || ~isfinite(roiID)
                return
            end

            app.deleteROIsByIDList(roiID, true);

        end

        function deleteSelectedROIs(app)
            %DELETESELECTEDROIS Delete ROI(s) selected in the ROI table.
            %
            %   The table Select column is the primary deletion target. If nothing is
            %   checked, SelectedROIID is used as a fallback when available.

            if isempty(app.ROIList)
                app.setStatusMessage('No ROIs available to delete.');
                app.updateGUIEnabledState();
                return
            end

            roiIDList = app.getSelectedROIIDList();

            % Do not fall back to SelectedROIID. The Select checkbox is the only delete
            % target for the Delete ROI button. This prevents deleting a newly created
            % ROI that is internally selected but not checked in the table.

            roiIDList = unique(double(roiIDList(:).'), 'stable');
            roiIDList = roiIDList(isfinite(roiIDList));

            if isempty(roiIDList)
                app.setStatusMessage('Select one or more ROIs to delete.');
                return
            end

            validIDMask = false(size(roiIDList));
            roiNames = cell(size(roiIDList));

            for iID = 1:numel(roiIDList)
                roiIdx = app.getROIIndexByID(roiIDList(iID));

                if ~isempty(roiIdx)
                    validIDMask(iID) = true;
                    roiNames{iID} = char(string(app.ROIList(roiIdx).name));
                end
            end

            roiIDList = roiIDList(validIDMask);
            roiNames = roiNames(validIDMask);

            if isempty(roiIDList)
                app.setStatusMessage('Selected ROI(s) were not found.');
                app.refreshROITable();
                app.updateGUIEnabledState();
                return
            end

            if numel(roiIDList) == 1
                promptText = sprintf('Delete ROI "%s"?', roiNames{1});
            else
                promptText = sprintf('Delete %d selected ROIs?', numel(roiIDList));
            end

            if ~app.confirmROIDeletion(promptText)
                app.setStatusMessage('ROI deletion cancelled.');
                return
            end

            app.deleteROIsByIDList(roiIDList, true);

        end

        function bConfirmed = confirmROIDeletion(app, promptText)
            %CONFIRMROIDeleTION Ask user to confirm ROI deletion.

            bConfirmed = false;

            try
                choice = uiconfirm(app.UIFigure, ...
                    promptText, ...
                    'Delete ROI', ...
                    'Options', {'Delete', 'Cancel'}, ...
                    'DefaultOption', 'Cancel', ...
                    'CancelOption', 'Cancel', ...
                    'Icon', 'warning');

                bConfirmed = strcmp(choice, 'Delete');
                return
            catch
            end

            try
                choice = questdlg(promptText, ...
                    'Delete ROI', ...
                    'Delete', ...
                    'Cancel', ...
                    'Cancel');

                bConfirmed = strcmp(choice, 'Delete');
            catch
                bConfirmed = false;
            end

        end

        function deleteROIsByIDList(app, roiIDList, bRefreshAfterDelete)
            %DELETEROISBYIDLIST Delete one or more ROIs by runtime ROI IDs.
            %
            %   deleteROIsByIDList(app, roiIDList, bRefreshAfterDelete) removes all ROI
            %   graphics, edit/runtime handles, table-selection state, trace graphics, and
            %   ROIList entries. Indices are deleted in descending order to avoid index
            %   shifts.

            if nargin < 3
                bRefreshAfterDelete = true;
            end

            if isempty(app.ROIList) || isempty(roiIDList)
                return
            end

            roiIDList = unique(double(roiIDList(:).'), 'stable');
            roiIDList = roiIDList(isfinite(roiIDList));

            if isempty(roiIDList)
                return
            end

            % Exit active edit runtime first. This prevents stale images.roi or group edit
            % graphics from surviving after the data model entry is removed.
            try
                app.deleteGroupEditRuntimeGraphics();
            catch
            end

            for iROI = 1:numel(app.ROIList)
                try
                    if isfield(app.ROIList(iROI), 'runtime') && ...
                            isfield(app.ROIList(iROI).runtime, 'editHandle') && ...
                            app.isUsableGraphicsHandle(app.ROIList(iROI).runtime.editHandle)
                        app.cleanupROIEditRuntimeByIndex(iROI);
                    end
                catch
                end
            end

            roiIdxList = [];

            for iID = 1:numel(roiIDList)
                roiIdx = app.getROIIndexByID(roiIDList(iID));

                if ~isempty(roiIdx)
                    roiIdxList(end+1) = roiIdx; %#ok<AGROW>
                end
            end

            roiIdxList = unique(roiIdxList, 'stable');

            if isempty(roiIdxList)
                if bRefreshAfterDelete
                    app.refreshROITable();
                    app.updateGUIEnabledState();
                end
                return
            end

            roiIdxList = sort(roiIdxList, 'descend');
            nDeleted = numel(roiIdxList);

            for iIdx = 1:numel(roiIdxList)
                roiIdx = roiIdxList(iIdx);

                if roiIdx < 1 || roiIdx > numel(app.ROIList)
                    continue
                end

                app.deleteROIGraphicsByIndex(roiIdx);
                app.ROIList(roiIdx) = [];
            end

            app.SelectedROIID = NaN;

            if ~isempty(app.ROIList)
                app.clearROISelectionState();

                try
                    validIDs = [app.ROIList.ID];

                    if ~isempty(validIDs)
                        app.SelectedROIID = NaN;
                    end
                catch
                    app.SelectedROIID = NaN;
                end
            end

            if bRefreshAfterDelete
                app.refreshROITable();
                app.refreshROITraces();

                % Rebuild event patch legend/stacking because deleting ROI traces can disturb
                % PlotAxes object order and legend state.
                app.refreshEventPatches();
                app.stackROITraceGraphics();

                app.setInteractionMode('idle');
                app.updateGUIEnabledState();

                if nDeleted == 1
                    app.setStatusMessage('Deleted 1 ROI.');
                else
                    app.setStatusMessage(sprintf('Deleted %d ROIs.', nDeleted));
                end
            end

        end

        function createDeleteROIButtonContextMenu(app)
            %CREATEDELETEROIBUTTONCONTEXTMENU Create Delete ROI button context menu.
            %
            %   Left-clicking DeleteROIButton keeps the standard behavior of deleting
            %   selected ROI(s). Right-clicking the button exposes explicit delete actions,
            %   including Delete all ROIs.

            if isempty(app.DeleteROIButton) || ~isvalid(app.DeleteROIButton)
                return
            end

            if ~isempty(app.DeleteROIContextMenu) && isvalid(app.DeleteROIContextMenu)
                try
                    app.DeleteROIButton.ContextMenu = app.DeleteROIContextMenu;
                catch
                end

                app.refreshDeleteROIContextMenuState();
                return
            end

            app.DeleteROIContextMenu = uicontextmenu(app.UIFigure);

            app.DeleteSelectedROIsMenu = uimenu(app.DeleteROIContextMenu);
            app.DeleteSelectedROIsMenu.Text = 'Delete selected ROI(s)';
            app.DeleteSelectedROIsMenu.MenuSelectedFcn = @(src, evt) app.deleteSelectedROIs();

            app.DeleteAllROIsMenu = uimenu(app.DeleteROIContextMenu);
            app.DeleteAllROIsMenu.Text = 'Delete all ROIs';
            app.DeleteAllROIsMenu.Separator = 'on';
            app.DeleteAllROIsMenu.MenuSelectedFcn = @(src, evt) app.deleteAllROIs();

            try
                app.DeleteROIButton.ContextMenu = app.DeleteROIContextMenu;
            catch
            end

            app.refreshDeleteROIContextMenuState();

        end

        function refreshDeleteROIContextMenuState(app)
            %REFRESHDELETEROICONTEXTMENUSTATE Enable/disable ROI delete menu items.

            hasROIs = ~isempty(app.ROIList);
            hasSelectedROIs = false;

            if hasROIs
                selectedIDs = app.getSelectedROIIDList();
                hasSelectedROIs = ~isempty(selectedIDs);
            end

            if ~isempty(app.DeleteSelectedROIsMenu) && isvalid(app.DeleteSelectedROIsMenu)
                if hasSelectedROIs
                    app.DeleteSelectedROIsMenu.Enable = 'on';
                else
                    app.DeleteSelectedROIsMenu.Enable = 'off';
                end
            end

            if ~isempty(app.DeleteAllROIsMenu) && isvalid(app.DeleteAllROIsMenu)
                if hasROIs
                    app.DeleteAllROIsMenu.Enable = 'on';
                else
                    app.DeleteAllROIsMenu.Enable = 'off';
                end
            end

        end


        function deleteAllROIs(app)
            %DELETEALLROIS Delete every ROI after user confirmation.

            if isempty(app.ROIList)
                app.setStatusMessage('No ROIs available to delete.');
                app.updateGUIEnabledState();
                return
            end

            nROI = numel(app.ROIList);

            if nROI == 1
                promptText = sprintf('Delete the only ROI "%s"?', app.ROIList(1).name);
            else
                promptText = sprintf('Delete all %d ROIs? This cannot be undone.', nROI);
            end

            if ~app.confirmROIDeletion(promptText)
                app.setStatusMessage('Delete all ROIs cancelled.');
                return
            end

            roiIDList = [app.ROIList.ID];
            app.deleteROIsByIDList(roiIDList, true);

            app.SelectedROIID = NaN;
            app.refreshDeleteROIContextMenuState();

            if nROI == 1
                app.setStatusMessage('Deleted all ROIs.');
            else
                app.setStatusMessage(sprintf('Deleted all %d ROIs.', nROI));
            end

        end

        function deleteROIGraphicsByIndex(app, roiIdx)
            %DELETEROIGRAPHICSBYINDEX Delete all graphics handles for one ROI.

            if isempty(app.ROIList) || roiIdx < 1 || roiIdx > numel(app.ROIList)
                return
            end

            if ~isfield(app.ROIList(roiIdx), 'runtime')
                return
            end

            app.cleanupROIEditRuntimeByIndex(roiIdx);

            runtimeFields = {'ROIHandle', 'editHandle', 'traceLine', 'tracePatch'};

            for iField = 1:numel(runtimeFields)
                fieldName = runtimeFields{iField};

                if ~isfield(app.ROIList(roiIdx).runtime, fieldName)
                    continue
                end

                h = app.ROIList(roiIdx).runtime.(fieldName);

                if app.isUsableGraphicsHandle(h)
                    try
                        delete(h);
                    catch
                    end
                end
            end

        end

        function stackROITraceGraphics(app)
            %STACKROITRACEGRAPHICS Keep event patches behind traces and timebar on top.

            try
                app.stackEventPatchesBottom();
            catch
            end

            for iROI = 1:numel(app.ROIList)
                if ~isfield(app.ROIList(iROI), 'runtime')
                    continue
                end

                try
                    hPatch = app.ROIList(iROI).runtime.tracePatch;
                    if app.isUsableGraphicsHandle(hPatch)
                        uistack(hPatch, 'top');
                    end
                catch
                end

                try
                    hLine = app.ROIList(iROI).runtime.traceLine;
                    if app.isUsableGraphicsHandle(hLine)
                        uistack(hLine, 'top');
                    end
                catch
                end
            end

            try
                if app.isUsableGraphicsHandle(app.CrossTraceSEMHandle)
                    uistack(app.CrossTraceSEMHandle, 'top');
                end

                if app.isUsableGraphicsHandle(app.CrossTraceHandle)
                    uistack(app.CrossTraceHandle, 'top');
                end

                if app.isUsableGraphicsHandle(app.TimeBarHandle)
                    uistack(app.TimeBarHandle, 'top');
                end
            catch
            end

        end

        function tf = isUsableGraphicsHandle(app, h) %#ok<INUSL>
            %ISUSABLEGRAPHICSHANDLE True for non-placeholder valid graphics handles.
            %
            %   A freshly initialized gobjects(1) value can be a
            %   GraphicsPlaceholder. It cannot safely receive graphics property
            %   updates such as XData, YData, Visible, or Color.

            tf = false;

            if isempty(h)
                return
            end

            try
                if isa(h, 'matlab.graphics.GraphicsPlaceholder')
                    return
                end
            catch
            end

            try
                tf = all(isvalid(h));
            catch
                tf = false;
            end
        end

        function [position, angleDeg] = rectangleParametersFromTransformedVertices(app, verticesXY) %#ok<INUSL>
            %RECTANGLEPARAMETERSFROMTRANSFORMEDVERTICES Approximate editable rectangle.
            %
            %   verticesXY is expected to be ordered around the transformed rectangle.
            %   The returned Position is the unrotated box used by drawrectangle, and
            %   RotationAngle follows the images.roi sign convention used elsewhere in
            %   this app.

            verticesXY = double(verticesXY);

            if size(verticesXY, 1) < 4 || size(verticesXY, 2) ~= 2
                minXY = min(verticesXY, [], 1);
                maxXY = max(verticesXY, [], 1);
                sizeXY = max(maxXY - minXY, 1);
                centerXY = minXY + sizeXY ./ 2;
                position = [centerXY(1) - sizeXY(1)/2, centerXY(2) - sizeXY(2)/2, sizeXY];
                angleDeg = 0;
                return
            end

            v = verticesXY(1:4, :);
            centerXY = mean(v, 1);

            edge1 = v(2, :) - v(1, :);
            edge2 = v(3, :) - v(2, :);

            widthValue = max(norm(edge1), 1);
            heightValue = max(norm(edge2), 1);

            % Stored vertices are generated with angleRad = -deg2rad(RotationAngle).
            % Invert that convention when reconstructing the images.roi angle.
            angleDeg = -rad2deg(atan2(edge1(2), edge1(1)));

            position = [centerXY(1) - widthValue/2, ...
                centerXY(2) - heightValue/2, ...
                widthValue, ...
                heightValue];

        end

        function [centerXY, semiAxes, angleDeg] = ellipseParametersFromTransformedVertices(app, verticesXY) %#ok<INUSL>
            %ELLIPSEPARAMETERSFROMTRANSFORMEDVERTICES Approximate editable ellipse.
            %
            %   The transformed ellipse outline is summarized with PCA. For uniformly
            %   sampled ellipse perimeter points, the covariance eigenvalues are
            %   approximately semiAxis^2 / 2.

            verticesXY = double(verticesXY);
            verticesXY = verticesXY(all(isfinite(verticesXY), 2), :);

            if size(verticesXY, 1) < 3 || size(verticesXY, 2) ~= 2
                centerXY = [mean(verticesXY(:, 1), 'omitnan'), mean(verticesXY(:, 2), 'omitnan')];
                if any(~isfinite(centerXY))
                    centerXY = [1 1];
                end
                semiAxes = [1 1];
                angleDeg = 0;
                return
            end

            centerXY = mean(verticesXY, 1);
            centeredXY = verticesXY - centerXY;

            C = (centeredXY.' * centeredXY) ./ max(size(centeredXY, 1), 1);
            [V, D] = eig(C);
            eigValues = diag(D);
            [eigValues, order] = sort(eigValues, 'descend');
            V = V(:, order);

            semiAxes = sqrt(max(2 .* eigValues(:).', eps));
            semiAxes = max(semiAxes(1:2), 0.5);

            majorVector = V(:, 1).';

            % Sign convention inverse of rotateVertices(..., RotationAngle), where
            % angleRad = -deg2rad(RotationAngle).
            angleDeg = -rad2deg(atan2(majorVector(2), majorVector(1)));

            if semiAxes(2) > semiAxes(1)
                semiAxes = fliplr(semiAxes);
                angleDeg = angleDeg + 90;
            end

            angleDeg = mod(angleDeg + 180, 360) - 180;

        end

        function ensureDataParamsMaskFields(app)
            %ENSUREDATAPARAMSMASKFIELDS Ensure app.DataParams has valid mask fields.
            %
            %   Empty DataParams.mask.logical is allowed and means no user-set mask. The
            %   runtime interpretation is a full TRUE mask with current image size.

            if isempty(app.DataParams) || ~isstruct(app.DataParams)
                app.DataParams = struct();
            end

            if ~isfield(app.DataParams, 'mask') || ~isstruct(app.DataParams.mask) || ...
                    ~isscalar(app.DataParams.mask)
                app.DataParams.mask = struct();
            end

            reqFields = {'logical','name','description','space','createdOn','source'};
            defaultValues = {[], '', '', 'native', [], ''};

            for iField = 1:numel(reqFields)
                if ~isfield(app.DataParams.mask, reqFields{iField})
                    app.DataParams.mask.(reqFields{iField}) = defaultValues{iField};
                end
            end

            if app.hasData()
                sz = app.getDataSize();
                imageSizeYX = double(sz(1:2));

                app.ensureDataParamsViewFields();
                app.DataParams.view.imageSizeYX = imageSizeYX;

                if ~isempty(app.DataParams.mask.logical)
                    mask = logical(app.DataParams.mask.logical);

                    if ~isequal(size(mask), imageSizeYX)
                        warning('DataViewer:LogicalMaskSizeMismatch', ...
                            'DataParams.mask.logical size does not match current image. Using full-field mask at runtime.');
                        app.DataParams.mask.logical = [];
                        app.DataParams.mask.name = '';
                        app.DataParams.mask.description = '';
                        app.DataParams.mask.createdOn = [];
                        app.DataParams.mask.source = '';
                    else
                        app.DataParams.mask.logical = mask;
                    end
                end
            end

        end

        function mask = getActiveLogicalMask(app)
            %GETACTIVELOGICALMASK Return user mask or full TRUE mask.

            if ~app.hasData()
                mask = [];
                return
            end

            sz = app.getDataSize();
            Ny = sz(1);
            Nx = sz(2);

            mask = [];

            try
                if isfield(app.DataParams, 'mask') && isfield(app.DataParams.mask, 'logical')
                    mask = app.DataParams.mask.logical;
                end
            catch
                mask = [];
            end

            if isempty(mask) || ~islogical(mask) || ~isequal(size(mask), [Ny Nx])
                mask = true(Ny, Nx);
            else
                mask = logical(mask);
            end

        end

        function tf = hasUserLogicalMask(app)
            %HASUSERLOGICALMASK True when DataParams contains a user-set logical mask.
            %
            %   A user-set mask may still be all TRUE. Runtime algorithms interpret an
            %   empty mask as full-field inclusion, so the distinction here is whether
            %   DataParams.mask.logical is populated and size-compatible.

            tf = false;

            if ~app.hasData()
                return
            end

            try
                mask = app.DataParams.mask.logical;
                sz = app.getDataSize();
                tf = islogical(mask) && isequal(size(mask), sz(1:2));
            catch
                tf = false;
            end

        end

        function clippedMask = clipROIMaskToActiveLogicalMask(app, roiMask)
            %CLIPROIMASKTOACTIVELOGICALMASK Intersect one ROI mask with active mask.

            if isempty(roiMask)
                clippedMask = [];
                return
            end

            roiMask = logical(roiMask);
            activeMask = app.getActiveLogicalMask();

            if isempty(activeMask) || ~isequal(size(activeMask), size(roiMask))
                clippedMask = roiMask;
                return
            end

            clippedMask = roiMask & activeMask;

        end

        function refreshLogicalMaskOverlay(app)
            %REFRESHLOGICALMASKOVERLAY Dim pixels outside the active logical mask.

            if ~app.hasData() || isempty(app.ImageAxes) || ~isvalid(app.ImageAxes)
                return
            end

            activeMask = app.getActiveLogicalMask();
            if isempty(activeMask)
                return
            end

            [Ny, Nx] = size(activeMask);
            outsideMask = ~activeMask;

            if isempty(app.LogicalMaskOverlayHandle) || ~isvalid(app.LogicalMaskOverlayHandle)
                hold(app.ImageAxes, 'on');
                app.LogicalMaskOverlayHandle = image(app.ImageAxes, ...
                    'CData', zeros(Ny, Nx, 3), ...
                    'AlphaData', 0.2 * double(outsideMask), ...
                    'XData', [1 Nx], ...
                    'YData', [1 Ny], ...
                    'HitTest', 'off', ...
                    'PickableParts', 'none', ...
                    'HandleVisibility', 'off');
            else
                app.LogicalMaskOverlayHandle.CData = zeros(Ny, Nx, 3);
                app.LogicalMaskOverlayHandle.AlphaData = 0.2 * double(outsideMask);
                app.LogicalMaskOverlayHandle.XData = [1 Nx];
                app.LogicalMaskOverlayHandle.YData = [1 Ny];
                app.LogicalMaskOverlayHandle.Visible = 'on';
            end

            % Keep the mask above the image but below crosshair, cache, and ROI overlays.
            try
                uistack(app.LogicalMaskOverlayHandle, 'top');
                if app.isUsableGraphicsHandle(app.CacheRectHandle)
                    uistack(app.CacheRectHandle, 'top');
                end
                if ~isempty(app.CrosshairHandles) && all(isvalid(app.CrosshairHandles))
                    uistack(app.CrosshairHandles, 'top');
                end
                if ~isempty(app.ROIList)
                    for iROI = 1:numel(app.ROIList)
                        if isfield(app.ROIList(iROI), 'runtime') && ...
                                isfield(app.ROIList(iROI).runtime, 'ROIHandle') && ...
                                app.isUsableGraphicsHandle(app.ROIList(iROI).runtime.ROIHandle)
                            uistack(app.ROIList(iROI).runtime.ROIHandle, 'top');
                        end
                    end
                end
            catch
            end

        end

        function saveLogicalMaskToDataParams(app, mask, maskName, descriptionText)
            %SAVELOGICALMASKTODATAPARAMS Persist logical mask into DataParams.mat.

            if ~app.hasData()
                error('DataViewer:NoDataLoaded', 'Load image data before saving a logical mask.');
            end

            folderPath = app.getCurrentDataFolder();

            if isempty(app.DataParams) || ~isstruct(app.DataParams) || ...
                    ~isfield(app.DataParams, 'schemaVersion')
                app.DataParams = app.ensureDataParamsFileForCurrentFolder();
            end

            app.ensureDataParamsViewFields();
            app.ensureDataParamsMaskFields();

            sz = app.getDataSize();
            app.DataParams.view.imageSizeYX = double(sz(1:2));

            if isempty(mask)
                app.DataParams.mask.logical = [];
                app.DataParams.mask.name = '';
                app.DataParams.mask.description = '';
                app.DataParams.mask.space = 'native';
                app.DataParams.mask.createdOn = [];
                app.DataParams.mask.source = '';
            else
                mask = logical(mask);
                if ~isequal(size(mask), sz(1:2))
                    error('DataViewer:LogicalMaskSizeMismatch', ...
                        'Logical mask size does not match current image size.');
                end

                if nargin < 3 || isempty(maskName)
                    maskName = 'user logical mask';
                end
                if nargin < 4 || isempty(descriptionText)
                    descriptionText = 'User-drawn inclusion mask from DataViewer.';
                end

                app.DataParams.mask.logical = mask;
                app.DataParams.mask.name = char(string(maskName));
                app.DataParams.mask.description = char(string(descriptionText));
                app.DataParams.mask.space = 'native';
                app.DataParams.mask.createdOn = datetime('now');
                app.DataParams.mask.source = 'DataViewer';
            end

            saveDataParams(folderPath, app.DataParams);

            S = load(fullfile(folderPath, 'DataParams.mat'), 'DataParams');
            if isfield(S, 'DataParams')
                app.DataParams = S.DataParams;
            end

            app.ensureDataParamsMaskFields();

        end

        function createLogicalMaskButtonContextMenu(app)
            %CREATELOGICALMASKBUTTONCONTEXTMENU Create Logical Mask button menu.

            if isempty(app.LogicalMaskButton) || ~isvalid(app.LogicalMaskButton)
                return
            end

            bReuseMenu = ...
                ~isempty(app.LogicalMaskButtonContextMenu) && ...
                isvalid(app.LogicalMaskButtonContextMenu) && ...
                ~isempty(app.CreateManualLogicalMaskMenu) && ...
                isvalid(app.CreateManualLogicalMaskMenu) && ...
                ~isempty(app.CreateAutomaticLogicalMaskMenu) && ...
                isvalid(app.CreateAutomaticLogicalMaskMenu) && ...
                ~isempty(app.ResetLogicalMaskMenu) && ...
                isvalid(app.ResetLogicalMaskMenu);

            if bReuseMenu
                try
                    app.LogicalMaskButton.ContextMenu = app.LogicalMaskButtonContextMenu;
                catch
                end
                app.refreshLogicalMaskButtonContextMenuState();
                return
            end

            try
                if ~isempty(app.LogicalMaskButtonContextMenu) && isvalid(app.LogicalMaskButtonContextMenu)
                    delete(app.LogicalMaskButtonContextMenu);
                end
            catch
            end

            app.LogicalMaskButtonContextMenu = uicontextmenu(app.UIFigure);

            app.CreateManualLogicalMaskMenu = uimenu(app.LogicalMaskButtonContextMenu);
            app.CreateManualLogicalMaskMenu.Text = 'Create Manual Mask';
            app.CreateManualLogicalMaskMenu.MenuSelectedFcn = @(src, evt) app.startLogicalMaskDrawing();

            app.CreateAutomaticLogicalMaskMenu = uimenu(app.LogicalMaskButtonContextMenu);
            app.CreateAutomaticLogicalMaskMenu.Text = 'Create Automatic Mask';
            app.CreateAutomaticLogicalMaskMenu.MenuSelectedFcn = @(src, evt) app.startAutomaticLogicalMaskDrawing();

            app.ResetLogicalMaskMenu = uimenu(app.LogicalMaskButtonContextMenu);
            app.ResetLogicalMaskMenu.Text = 'Clear Existing Mask';
            app.ResetLogicalMaskMenu.Separator = 'on';
            app.ResetLogicalMaskMenu.MenuSelectedFcn = @(src, evt) app.resetLogicalMask();

            try
                app.LogicalMaskButton.ContextMenu = app.LogicalMaskButtonContextMenu;
            catch
            end

            app.refreshLogicalMaskButtonContextMenuState();

        end

        function refreshLogicalMaskButtonContextMenuState(app)
            %REFRESHLOGICALMASKBUTTONCONTEXTMENUSTATE Enable Logical Mask menu items.

            hasDataIdle = app.hasData() && strcmp(app.InteractionMode, 'idle');

            if ~isempty(app.CreateManualLogicalMaskMenu) && isvalid(app.CreateManualLogicalMaskMenu)
                if hasDataIdle
                    app.CreateManualLogicalMaskMenu.Enable = 'on';
                else
                    app.CreateManualLogicalMaskMenu.Enable = 'off';
                end
            end

            if ~isempty(app.CreateAutomaticLogicalMaskMenu) && isvalid(app.CreateAutomaticLogicalMaskMenu)
                if hasDataIdle
                    app.CreateAutomaticLogicalMaskMenu.Enable = 'on';
                else
                    app.CreateAutomaticLogicalMaskMenu.Enable = 'off';
                end
            end

            if ~isempty(app.ResetLogicalMaskMenu) && isvalid(app.ResetLogicalMaskMenu)
                if hasDataIdle && app.hasUserLogicalMask()
                    app.ResetLogicalMaskMenu.Enable = 'on';
                else
                    app.ResetLogicalMaskMenu.Enable = 'off';
                end
            end

        end

        function resetLogicalMask(app)
            %RESETLOGICALMASK Clear user mask and return to full-field runtime mask.

            if ~app.hasData()
                app.setStatusMessage('Load image data before resetting the logical mask.');
                return
            end

            if ~app.hasUserLogicalMask()
                app.setStatusMessage('Logical mask is already full-field.');
                return
            end

            try
                choice = uiconfirm(app.UIFigure, ...
                    'Reset the logical mask to full image inclusion?', ...
                    'Reset logical mask', ...
                    'Options', {'Reset Mask', 'Cancel'}, ...
                    'DefaultOption', 'Cancel', ...
                    'CancelOption', 'Cancel', ...
                    'Icon', 'warning');
                if ~strcmp(choice, 'Reset Mask')
                    app.setStatusMessage('Logical mask reset cancelled.');
                    return
                end
            catch
                % Continue without modal confirmation when uiconfirm is unavailable.
            end

            app.saveLogicalMaskToDataParams([], '', '');
            app.refreshLogicalMaskOverlay();
            app.refreshLogicalMaskButtonContextMenuState();
            app.setStatusMessage('Logical mask reset to full image inclusion.');

        end

        function startLogicalMaskDrawing(app)
            %STARTLOGICALMASKDRAWING Draw one or more polygons to define logical mask.

            app.startLogicalMaskDrawingFromPolygons( ...
                {}, ...
                'user logical mask', ...
                'User-drawn inclusion mask from DataViewer.');

        end

        function startAutomaticLogicalMaskDrawing(app)
            %STARTAUTOMATICLOGICALMASKDRAWING Initialize editable mask polygons automatically.

            if ~app.hasData()
                app.setStatusMessage('Load image data before creating a logical mask.');
                return
            end

            if app.hasUserLogicalMask()
                try
                    choice = uiconfirm(app.UIFigure, ...
                        ['An existing logical mask is already saved. Automatic mask creation ' ...
                        'will replace it if you confirm the edited draft. Continue?'], ...
                        'Existing logical mask', ...
                        'Options', {'Continue', 'Cancel'}, ...
                        'DefaultOption', 'Cancel', ...
                        'CancelOption', 'Cancel', ...
                        'Icon', 'warning');

                    if ~strcmp(choice, 'Continue')
                        app.setStatusMessage('Automatic logical mask creation cancelled.');
                        return
                    end
                catch
                    % Continue without modal confirmation when uiconfirm is unavailable.
                end
            end

            try
                app.setStatusMessage('Computing temporal average for automatic logical mask...');
                drawnow limitrate

                avgImg = app.computeDisplayedTemporalAverage();
                autoMask = app.createAutomaticLogicalMaskFromAverage(avgImg);

                if isempty(autoMask) || ~any(autoMask(:))
                    app.alertAutomaticLogicalMaskFailed();
                    return
                end

                maxRegions = 6;
                polygonList = app.logicalMaskToDraftPolygons(autoMask, maxRegions);

                if isempty(polygonList)
                    app.alertAutomaticLogicalMaskFailed();
                    return
                end

                app.startLogicalMaskDrawingFromPolygons( ...
                    polygonList, ...
                    'automatic logical mask', ...
                    'Automatic threshold-based inclusion mask from DataViewer.');

            catch ME
                app.setStatusMessage(sprintf('Automatic logical mask creation failed: %s', ME.message));

                try
                    uialert(app.UIFigure, ...
                        ['Automatic logical mask creation failed. Create the logical mask manually.' newline newline ME.message], ...
                        'Automatic logical mask failed', ...
                        'Icon', 'warning');
                catch
                end
            end

        end

        function startLogicalMaskDrawingFromPolygons(app, initialPolygons, maskName, descriptionText)
            %STARTLOGICALMASKDRAWINGFROMPOLYGONS Draw/refine logical-mask polygons.
            %
            %   initialPolygons is a cell array of [x y] vertex arrays. When empty, the
            %   method starts the standard manual drawpolygon flow.

            if nargin < 2 || isempty(initialPolygons)
                initialPolygons = {};
            end

            if nargin < 3 || isempty(maskName)
                maskName = 'user logical mask';
            end

            if nargin < 4 || isempty(descriptionText)
                descriptionText = 'User-drawn inclusion mask from DataViewer.';
            end

            if ~iscell(initialPolygons)
                initialPolygons = {initialPolygons};
            end

            if ~app.hasData()
                app.setStatusMessage('Load image data before creating a logical mask.');
                return
            end

            if strcmp(app.InteractionMode, 'playingMovie')
                app.stopMoviePlayback();
            end

            if ~strcmp(app.InteractionMode, 'idle')
                app.stopActiveROIEditForSelectionChange();
            end

            sz = app.getDataSize();
            Ny = sz(1);
            Nx = sz(2);
            drawingArea = [0.5, 0.5, Nx, Ny];

            previousMode = app.InteractionMode;
            app.LogicalMaskPreviousKeyFcn = app.UIFigure.WindowKeyPressFcn;
            app.LogicalMaskPreviousImageAxesContextMenu = app.ImageAxes.ContextMenu;
            app.LogicalMaskPreviousImageHandleContextMenu = app.ImageHandle.ContextMenu;

            bDone = false;
            bConfirmed = false;

            app.deleteLogicalMaskDraftRuntime();
            app.createLogicalMaskDrawingContextMenu();

            % Bind this drawing session's nested callbacks to the context menu.
            app.LogicalMaskAddPolygonMenu.MenuSelectedFcn = @(src, evt) addDraftPolygon([]);
            app.LogicalMaskConfirmMenu.MenuSelectedFcn = @(src, evt) confirmFromContextMenu();
            app.LogicalMaskCancelMenu.MenuSelectedFcn = @(src, evt) cancelFromContextMenu();

            try
                app.setInteractionMode('drawingLogicalMask');
                app.UIFigure.WindowKeyPressFcn = @onKeyPress;

                % Attach drawing context menu only when it belongs to the same figure as the
                % graphics object. This avoids stale-menu parent mismatch after reset/reuse.
                app.safeSetGraphicsContextMenu(app.ImageAxes, app.LogicalMaskContextMenu);
                app.safeSetGraphicsContextMenu(app.ImageHandle, app.LogicalMaskContextMenu);

                if isempty(initialPolygons)
                    app.setStatusMessage(['Draw logical mask polygon. Right-click to add another polygon. ' ...
                        'Double-click a polygon or press Enter to confirm. Press Escape to cancel.']);

                    addDraftPolygon([]);
                else
                    app.setStatusMessage(['Automatic logical mask draft created. Refine polygon(s), right-click to add another polygon, ' ...
                        'double-click a polygon or press Enter to confirm. Press Escape to cancel.']);

                    for iPoly = 1:numel(initialPolygons)
                        addDraftPolygon(initialPolygons{iPoly});
                    end
                end

                while ~bDone && isvalid(app.UIFigure)
                    drawnow
                    pause(0.02)
                end

                if bConfirmed
                    confirmMaskDrawing();
                else
                    app.setStatusMessage('Logical mask creation cancelled.');
                end

            catch ME
                app.setStatusMessage(sprintf('Logical mask creation failed: %s', ME.message));
            end

            cleanupMaskDrawingRuntime();

            try
                app.setInteractionMode(previousMode);
            catch
                app.setInteractionMode('idle');
            end

            function addDraftPolygon(positionXY)
                if ~isvalid(app.UIFigure)
                    return
                end

                hold(app.ImageAxes, 'on');

                if nargin < 1 || isempty(positionXY)
                    hPoly = drawpolygon(app.ImageAxes, ...
                        'Color', [1 1 0], ...
                        'FaceAlpha', 0.12, ...
                        'LineWidth', 1.5, ...
                        'DrawingArea', drawingArea, ...
                        'InteractionsAllowed', 'all');
                else
                    positionXY = double(positionXY);

                    if size(positionXY, 2) ~= 2 || size(positionXY, 1) < 3 || ...
                            any(~isfinite(positionXY(:)))
                        return
                    end

                    hPoly = drawpolygon(app.ImageAxes, ...
                        'Position', positionXY, ...
                        'Color', [1 1 0], ...
                        'FaceAlpha', 0.12, ...
                        'LineWidth', 1.5, ...
                        'DrawingArea', drawingArea, ...
                        'InteractionsAllowed', 'all');
                end

                if isempty(hPoly) || ~isvalid(hPoly)
                    return
                end

                app.setROIObjectPropertyIfAvailable(hPoly, 'Label', sprintf('Mask polygon %d', numel(app.LogicalMaskDraftHandles) + 1));
                app.setROIObjectPropertyIfAvailable(hPoly, 'LabelVisible', 'on');

                app.safeSetGraphicsContextMenu(hPoly, app.LogicalMaskContextMenu);

                app.LogicalMaskDraftHandles(end+1) = hPoly;

                listeners = {};
                listeners{end+1} = addlistener(hPoly, 'ROIClicked', @onPolygonClicked);
                try
                    listeners{end+1} = addlistener(hPoly, 'DeletingROI', @onPolygonDeleting);
                catch
                end
                app.LogicalMaskDraftListeners{end+1} = listeners;
            end

            function onPolygonClicked(~, evt)
                try
                    if isprop(evt, 'SelectionType') && strcmpi(evt.SelectionType, 'double')
                        bConfirmed = true;
                        bDone = true;
                    end
                catch
                end
            end

            function onPolygonDeleting(~, ~)
                drawnow limitrate
                validMask = arrayfun(@(h) app.isUsableGraphicsHandle(h), app.LogicalMaskDraftHandles);
                app.LogicalMaskDraftHandles = app.LogicalMaskDraftHandles(validMask);
            end

            function onKeyPress(~, evt)
                switch lower(evt.Key)
                    case {'return', 'enter'}
                        bConfirmed = true;
                        bDone = true;
                    case 'escape'
                        bConfirmed = false;
                        bDone = true;
                end
            end

            function confirmMaskDrawing()
                newMask = false(Ny, Nx);

                for iMask = 1:numel(app.LogicalMaskDraftHandles)
                    h = app.LogicalMaskDraftHandles(iMask);
                    if ~app.isUsableGraphicsHandle(h)
                        continue
                    end

                    try
                        subMask = createMask(h, app.ImageHandle);
                    catch
                        subMask = createMask(h);
                    end

                    if isequal(size(subMask), [Ny Nx])
                        newMask = newMask | logical(subMask);
                    end
                end

                if ~any(newMask(:))
                    app.setStatusMessage('Logical mask creation cancelled: empty mask.');
                    return
                end

                emptyROIIDs = app.getROIIDsThatWouldBeEmptyAfterLogicalMask(newMask);
                if ~isempty(emptyROIIDs)
                    if numel(emptyROIIDs) == 1
                        promptText = '1 ROI falls fully outside the new logical mask and will be deleted. Continue?';
                    else
                        promptText = sprintf('%d ROIs fall fully outside the new logical mask and will be deleted. Continue?', numel(emptyROIIDs));
                    end

                    try
                        choice = uiconfirm(app.UIFigure, promptText, ...
                            'Logical mask clips ROIs', ...
                            'Options', {'Continue', 'Cancel'}, ...
                            'DefaultOption', 'Cancel', ...
                            'CancelOption', 'Cancel', ...
                            'Icon', 'warning');
                        if ~strcmp(choice, 'Continue')
                            app.setStatusMessage('Logical mask creation cancelled.');
                            return
                        end
                    catch
                    end
                end

                app.saveLogicalMaskToDataParams(newMask, maskName, descriptionText);
                app.refreshLogicalMaskOverlay();
                [nUpdated, nDeleted] = app.clipExistingROIsToActiveLogicalMask(false);

                app.refreshROITable();
                app.refreshROITraces();
                app.refreshEventPatches();
                app.stackROITraceGraphics();
                app.updateGUIEnabledState();

                if nUpdated > 0 || nDeleted > 0
                    app.setStatusMessage(sprintf('Logical mask saved. Updated %d ROI(s); deleted %d ROI(s) outside mask.', nUpdated, nDeleted));
                else
                    app.setStatusMessage('Logical mask saved.');
                end
            end

            function cleanupMaskDrawingRuntime()
                try
                    app.UIFigure.WindowKeyPressFcn = app.LogicalMaskPreviousKeyFcn;
                catch
                end
                app.restoreGraphicsContextMenu(app.ImageAxes, app.LogicalMaskPreviousImageAxesContextMenu);
                app.restoreGraphicsContextMenu(app.ImageHandle, app.LogicalMaskPreviousImageHandleContextMenu);

                app.deleteLogicalMaskDraftRuntime();

                app.LogicalMaskPreviousImageAxesContextMenu = [];
                app.LogicalMaskPreviousImageHandleContextMenu = [];
                app.LogicalMaskPreviousKeyFcn = [];
            end

            function confirmFromContextMenu()
                bConfirmed = true;
                bDone = true;
            end

            function cancelFromContextMenu()
                bConfirmed = false;
                bDone = true;
            end

        end

        function avgImg = computeDisplayedTemporalAverage(app)
            %COMPUTEDISPLAYEDTEMPORALAVERAGE Average the currently displayed frame stack.
            %
            %   This method uses getCurrentFrame so it respects the active viewer mode,
            %   including event-mode display, selected condition/repetition, and UMT entry.

            if ~app.hasData()
                error('DataViewer:NoDataLoaded', ...
                    'Load image data before computing an automatic logical mask.');
            end

            sz = app.getDataSize();
            Ny = sz(1);
            Nx = sz(2);
            nFrames = max(1, round(sz(3)));

            sumImg = zeros(Ny, Nx, 'double');
            countImg = zeros(Ny, Nx, 'double');

            originalFrame = app.CurrentFrame;
            cleanupFrame = onCleanup(@() restoreOriginalFrame());

            for iFrame = 1:nFrames
                app.CurrentFrame = iFrame;
                frame = double(app.getCurrentFrame());

                if ~isequal(size(frame), [Ny Nx])
                    frame = squeeze(frame);
                end

                if ~isequal(size(frame), [Ny Nx])
                    error('DataViewer:InvalidDisplayedFrameSize', ...
                        'Displayed frame size does not match the active image size.');
                end

                valid = isfinite(frame);
                sumImg(valid) = sumImg(valid) + frame(valid);
                countImg(valid) = countImg(valid) + 1;

                if mod(iFrame, 50) == 0 || iFrame == nFrames
                    app.setStatusMessage(sprintf( ...
                        'Computing temporal average for automatic logical mask... %d/%d frames', ...
                        iFrame, nFrames));
                    drawnow limitrate
                end
            end

            avgImg = nan(Ny, Nx);
            validCount = countImg > 0;
            avgImg(validCount) = sumImg(validCount) ./ countImg(validCount);

            function restoreOriginalFrame()
                try
                    app.CurrentFrame = min(max(round(originalFrame), 1), nFrames);
                    app.refreshImageFrame();
                    title(app.ImageAxes, app.getImageTitle());
                    app.refreshTimeBar();
                    app.refreshFrameControls();
                catch
                end
            end

        end

        function mask = createAutomaticLogicalMaskFromAverage(app, avgImg)
            %CREATEAUTOMATICLOGICALMASKFROMAVERAGE Threshold and clean temporal average.
            %
            %   The initial threshold keeps the highest 10 percent of finite average
            %   intensities. Morphology uses image-size-relative parameters.

            mask = [];

            if isempty(avgImg) || ~ismatrix(avgImg)
                return
            end

            avgImg = double(avgImg);
            finiteMask = isfinite(avgImg);

            finiteValues = avgImg(finiteMask);
            if isempty(finiteValues)
                return
            end

            finiteValues = sort(finiteValues(:), 'ascend');
            nFinite = numel(finiteValues);

            if nFinite < 2 || finiteValues(1) == finiteValues(end)
                return
            end

            topFraction = 0.10;
            thresholdIdx = max(1, min(nFinite, ceil((1 - topFraction) * nFinite)));
            thresholdValue = finiteValues(thresholdIdx);

            mask = finiteMask & avgImg >= thresholdValue;

            if ~any(mask(:))
                mask = [];
                return
            end

            [Ny, Nx] = size(mask);
            imageArea = Ny * Nx;
            minDim = max(1, min(Ny, Nx));

            minRegionArea = max(1, round(0.001 * imageArea));
            closeRadius = max(1, round(0.010 * minDim));
            openRadius = max(1, round(0.003 * minDim));

            try
                closeSE = strel('disk', closeRadius, 0);
            catch
                closeSE = strel('disk', closeRadius);
            end

            try
                openSE = strel('disk', openRadius, 0);
            catch
                openSE = strel('disk', openRadius);
            end

            % Merge nearby high-intensity fragments into broader regions.
            mask = imclose(mask, closeSE);

            % Remove small holes inside candidate regions.
            mask = imfill(mask, 'holes');

            % Remove small isolated regions without imposing an upper size limit.
            mask = bwareaopen(mask, minRegionArea);

            if ~any(mask(:))
                mask = [];
                return
            end

            % Light opening removes very thin accidental bridges/spurs. Fall back to
            % the pre-open mask if this is too aggressive for the current image.
            maskBeforeOpen = mask;
            maskAfterOpen = imopen(mask, openSE);
            maskAfterOpen = imfill(maskAfterOpen, 'holes');
            maskAfterOpen = bwareaopen(maskAfterOpen, minRegionArea);

            if any(maskAfterOpen(:))
                mask = maskAfterOpen;
            else
                mask = maskBeforeOpen;
            end

            mask = logical(mask);

        end

        function polygonList = logicalMaskToDraftPolygons(app, mask, maxRegions)
            %LOGICALMASKTODRAFTPOLYGONS Convert mask components to editable polygons.
            %
            %   Keeps up to maxRegions largest regions.

            if nargin < 3 || isempty(maxRegions)
                maxRegions = 6;
            end

            polygonList = {};

            if isempty(mask) || ~islogical(mask) || ~any(mask(:))
                return
            end

            maxRegions = max(1, round(maxRegions));

            CC = bwconncomp(mask);
            if CC.NumObjects == 0
                return
            end

            regionAreas = cellfun(@numel, CC.PixelIdxList);
            [~, order] = sort(regionAreas, 'descend');
            order = order(1:min(maxRegions, numel(order)));

            [Ny, Nx] = size(mask);


            for iRegion = 1:numel(order)
                compMask = false(Ny, Nx);
                compMask(CC.PixelIdxList{order(iRegion)}) = true;

                boundaries = bwboundaries(compMask, 'noholes');
                if isempty(boundaries)
                    continue
                end

                boundaryLengths = cellfun(@(b) size(b, 1), boundaries);
                [~, iLongest] = max(boundaryLengths);

                boundaryYX = boundaries{iLongest};

                if size(boundaryYX, 1) < 3
                    continue
                end

                verticesXY = [boundaryYX(:, 2), boundaryYX(:, 1)];
                verticesXY = double(verticesXY);

                if size(verticesXY, 1) > 2 && isequal(verticesXY(1, :), verticesXY(end, :))
                    verticesXY(end, :) = [];
                end

                verticesXY(:, 1) = min(max(verticesXY(:, 1), 1), Nx);
                verticesXY(:, 2) = min(max(verticesXY(:, 2), 1), Ny);

                % Simplify only if available. This keeps the draft interactive on large
                % masks while preserving the overall periphery.
                % reducepoly uses a normalized tolerance. Keep it conservative so the
                % automatic mask follows the detected boundary instead of collapsing into a
                % coarse polygon.
                maxDraftVertices = 1200;

                if size(verticesXY, 1) > maxDraftVertices
                    simplifyTol = max(0.001, max(0.00005, 1 / (4 * size(verticesXY, 1))));

                    try
                        verticesReduced = reducepoly(verticesXY, simplifyTol);

                        % Accept simplification only if it remains a meaningful polygon.
                        if size(verticesReduced, 1) >= 12
                            verticesXY = verticesReduced;
                        end
                    catch
                    end
                end

                verticesXY = app.cleanROIVertices(verticesXY);

                if size(verticesXY, 1) < 3
                    continue
                end

                polygonList{end+1} = verticesXY; %#ok<AGROW>
            end

        end

        function alertAutomaticLogicalMaskFailed(app)
            %ALERTAUTOMATICLOGICALMASKFAILED Tell user to create mask manually.

            msg = ['Automatic thresholding did not identify a usable logical mask. ' ...
                'Create the logical mask manually.'];

            app.setStatusMessage(msg);

            try
                uialert(app.UIFigure, ...
                    msg, ...
                    'Automatic logical mask failed', ...
                    'Icon', 'warning');
            catch
            end

        end

        function createLogicalMaskDrawingContextMenu(app)
            %CREATELOGICALMASKDRAWINGCONTEXTMENU Create right-click menu while drawing mask.
            %
            %   The context menu must belong to the same figure as ImageAxes/ImageHandle.
            %   MATLAB raises an error if a UIContextMenu from one figure is assigned to
            %   an object parented to another figure. Recreate the menu when stale.

            parentFig = [];

            try
                parentFig = ancestor(app.ImageAxes, 'figure');
            catch
                parentFig = [];
            end

            if isempty(parentFig) || ~isvalid(parentFig)
                parentFig = app.UIFigure;
            end

            bReuseMenu = ...
                ~isempty(app.LogicalMaskContextMenu) && ...
                isvalid(app.LogicalMaskContextMenu);

            if bReuseMenu
                try
                    bReuseMenu = isequal(ancestor(app.LogicalMaskContextMenu, 'figure'), parentFig);
                catch
                    bReuseMenu = false;
                end
            end

            if bReuseMenu
                return
            end

            try
                if ~isempty(app.LogicalMaskContextMenu) && isvalid(app.LogicalMaskContextMenu)
                    delete(app.LogicalMaskContextMenu);
                end
            catch
            end

            app.LogicalMaskContextMenu = uicontextmenu(parentFig);

            app.LogicalMaskAddPolygonMenu = uimenu(app.LogicalMaskContextMenu);
            app.LogicalMaskAddPolygonMenu.Text = 'Add mask polygon';

            app.LogicalMaskConfirmMenu = uimenu(app.LogicalMaskContextMenu);
            app.LogicalMaskConfirmMenu.Text = 'Confirm logical mask';
            app.LogicalMaskConfirmMenu.Separator = 'on';

            app.LogicalMaskCancelMenu = uimenu(app.LogicalMaskContextMenu);
            app.LogicalMaskCancelMenu.Text = 'Cancel mask drawing';

        end

        function deleteLogicalMaskDraftRuntime(app)
            %DELETELOGICALMASKDRAFTRUNTIME Delete temporary mask drawing handles/listeners.

            for iCell = 1:numel(app.LogicalMaskDraftListeners)
                listeners = app.LogicalMaskDraftListeners{iCell};
                for iListener = 1:numel(listeners)
                    try
                        if ~isempty(listeners{iListener}) && isvalid(listeners{iListener})
                            delete(listeners{iListener});
                        end
                    catch
                    end
                end
            end
            app.LogicalMaskDraftListeners = {};

            if ~isempty(app.LogicalMaskDraftHandles)
                for iHandle = 1:numel(app.LogicalMaskDraftHandles)
                    try
                        if app.isUsableGraphicsHandle(app.LogicalMaskDraftHandles(iHandle))
                            delete(app.LogicalMaskDraftHandles(iHandle));
                        end
                    catch
                    end
                end
            end

            app.LogicalMaskDraftHandles = gobjects(0);

        end

        function roiIDList = getROIIDsThatWouldBeEmptyAfterLogicalMask(app, candidateMask)
            %GETROIIDTHATWOULDBEEMPTYAFTERLOGICALMASK Find ROIs erased by candidate mask.

            roiIDList = [];

            if isempty(app.ROIList) || isempty(candidateMask)
                return
            end

            candidateMask = logical(candidateMask);

            for iROI = 1:numel(app.ROIList)
                if ~isfield(app.ROIList(iROI), 'mask') || isempty(app.ROIList(iROI).mask)
                    continue
                end

                roiMask = logical(app.ROIList(iROI).mask);
                if ~isequal(size(roiMask), size(candidateMask))
                    continue
                end

                if ~any(roiMask(:) & candidateMask(:))
                    roiIDList(end+1) = app.ROIList(iROI).ID; %#ok<AGROW>
                end
            end

        end

        function [nUpdated, nDeleted] = clipExistingROIsToActiveLogicalMask(app, bAskBeforeDelete)
            %CLIPEXISTINGROISTOACTIVELOGICALMASK Clip current ROIs to active mask.
            %
            %   Split intersections are preserved as multiple polygon ROIs. If one ROI is
            %   split into several connected components, the original ROI is replaced by:
            %
            %       <OriginalName>_part1
            %       <OriginalName>_part2
            %       ...
            %
            %   The replacement is performed from the end of ROIList to the beginning so
            %   ROI indices remain valid while entries are deleted/inserted.

            if nargin < 2
                bAskBeforeDelete = true;
            end

            nUpdated = 0;
            nDeleted = 0;

            if isempty(app.ROIList) || ~app.hasData()
                return
            end

            activeMask = app.getActiveLogicalMask();
            if isempty(activeMask)
                return
            end

            emptyIDs = app.getROIIDsThatWouldBeEmptyAfterLogicalMask(activeMask);

            if ~isempty(emptyIDs) && bAskBeforeDelete
                if numel(emptyIDs) == 1
                    promptText = '1 ROI falls fully outside the logical mask and will be deleted. Continue?';
                else
                    promptText = sprintf('%d ROIs fall fully outside the logical mask and will be deleted. Continue?', numel(emptyIDs));
                end

                try
                    choice = uiconfirm(app.UIFigure, promptText, ...
                        'Logical mask clips ROIs', ...
                        'Options', {'Continue', 'Cancel'}, ...
                        'DefaultOption', 'Cancel', ...
                        'CancelOption', 'Cancel', ...
                        'Icon', 'warning');

                    if ~strcmp(choice, 'Continue')
                        return
                    end
                catch
                end
            end

            for iROI = numel(app.ROIList):-1:1
                if ~isfield(app.ROIList(iROI), 'mask') || isempty(app.ROIList(iROI).mask)
                    continue
                end

                oldMask = logical(app.ROIList(iROI).mask);
                if ~isequal(size(oldMask), size(activeMask))
                    continue
                end

                clippedMask = oldMask & activeMask;

                if ~any(clippedMask(:))
                    app.deleteROIGraphicsByIndex(iROI);
                    app.ROIList(iROI) = [];
                    nDeleted = nDeleted + 1;
                    continue
                end

                if isequal(clippedMask, oldMask)
                    continue
                end

                componentMasks = app.maskToConnectedComponentMasks(clippedMask);

                if isempty(componentMasks)
                    app.deleteROIGraphicsByIndex(iROI);
                    app.ROIList(iROI) = [];
                    nDeleted = nDeleted + 1;
                    continue
                end

                baseName = char(string(app.ROIList(iROI).name));
                bUsePartNames = numel(componentMasks) > 1;

                nParts = app.replaceROIByMaskComponents(iROI, componentMasks, baseName, bUsePartNames);

                if nParts > 0
                    nUpdated = nUpdated + nParts;
                else
                    % All candidate components failed polygon conversion. The original ROI was
                    % removed by replaceROIByMaskComponents, so count it as deleted.
                    nDeleted = nDeleted + 1;
                end
            end

            app.SelectedROIID = NaN;
            app.clearROISelectionState();

        end

        function verticesXY = maskToSimplifiedPolygonVertices(app, mask)
            %MASKTOSIMPLIFIEDPOLYGONVERTICES Convert one mask component to polygon vertices.
            %
            %   This helper preserves the current strict mask-preserving simplification
            %   strategy. If a multi-component mask is supplied by mistake, the largest
            %   component is used for backward compatibility. Use
            %   maskToConnectedComponentMasks before this function when all split regions
            %   must be preserved.

            verticesXY = zeros(0, 2);

            if isempty(mask) || ~any(mask(:))
                return
            end

            mask = logical(mask);

            componentMasks = app.maskToConnectedComponentMasks(mask);
            if isempty(componentMasks)
                return
            end

            keepMask = componentMasks{1};

            B = bwboundaries(keepMask, 'noholes');
            if isempty(B)
                return
            end

            % Choose the longest boundary from this connected component.
            [~, idxBoundary] = max(cellfun(@(x) size(x, 1), B));
            boundaryRC = double(B{idxBoundary});       % [row col]
            xy = [boundaryRC(:, 2), boundaryRC(:, 1)]; % [x y]

            if size(xy, 1) < 3
                return
            end

            boundaryMask = bwperim(keepMask);
            maxAllowedDiff = max(1, round(0.005 * nnz(keepMask)));
            maxAllowedDiff = max(maxAllowedDiff, round(0.25 * nnz(boundaryMask)));

            % reducepoly accepts tolerances in [0 1].
            candidateTol = [1 0.75 0.5 0.25 0.1 0.05 0.025 0.01 0.005 0];

            bestXY = xy;

            for iTol = 1:numel(candidateTol)
                tol = candidateTol(iTol);

                if tol > 0
                    xyCandidate = reducepoly(xy, tol);
                else
                    xyCandidate = xy;
                end

                xyCandidate = app.cleanROIVertices(xyCandidate);
                if size(xyCandidate, 1) < 3
                    continue
                end

                testMask = poly2mask(xyCandidate(:, 1), xyCandidate(:, 2), size(mask, 1), size(mask, 2));
                nDiff = nnz(xor(testMask, keepMask));

                if nDiff <= maxAllowedDiff
                    bestXY = xyCandidate;
                    break
                end
            end

            verticesXY = app.cleanROIVertices(bestXY);

        end

        function componentMasks = maskToConnectedComponentMasks(app, mask) %#ok<INUSD>
            %MASKTOCONNECTEDCOMPONENTMASKS Split a logical mask into component masks.
            %
            %   Components are sorted from largest to smallest so _part1 corresponds to
            %   the dominant remaining region.

            componentMasks = {};

            if isempty(mask) || ~any(mask(:))
                return
            end

            mask = logical(mask);

            CC = bwconncomp(mask, 8);
            if CC.NumObjects < 1
                return
            end

            nPix = cellfun(@numel, CC.PixelIdxList);
            [~, order] = sort(nPix, 'descend');

            componentMasks = cell(1, numel(order));

            for iComp = 1:numel(order)
                thisMask = false(size(mask));
                thisMask(CC.PixelIdxList{order(iComp)}) = true;
                componentMasks{iComp} = thisMask;
            end

        end

        function roiID = getNextAvailableROIID(app)
            %GETNEXTAVAILABLEROIID Return next available runtime ROI ID.

            if isempty(app.ROIList)
                roiID = 1;
                return
            end

            roiIDs = [app.ROIList.ID];
            roiIDs = roiIDs(isfinite(roiIDs));

            if isempty(roiIDs)
                roiID = 1;
            else
                roiID = max(roiIDs) + 1;
            end

        end

        function ROI = makePolygonROIFromMaskComponent(app, templateROI, componentMask, roiName, roiID, bCreateOverlay)
            %MAKEPOLYGONROIFROMMASKCOMPONENT Build one polygon ROI from a mask component.
            %
            %   The template ROI supplies color, notes, visibility and trace structure
            %   defaults. Geometry and mask fields are rebuilt from componentMask.
            %   bCreateOverlay controls whether the static display handle is created
            %   immediately. Load/import paths pass false and rebuild overlays after
            %   ROIList assignment.

            if nargin < 6 || isempty(bCreateOverlay)
                bCreateOverlay = true;
            end

            if isempty(componentMask) || ~any(componentMask(:))
                error('DataViewer:EmptyROIComponent', ...
                    'Cannot create ROI from an empty mask component.');
            end

            verticesXY = app.maskToSimplifiedPolygonVertices(componentMask);
            verticesXY = app.cleanROIVertices(verticesXY);

            if size(verticesXY, 1) < 3
                error('DataViewer:InvalidROIComponent', ...
                    'Mask component could not be converted to a valid polygon.');
            end

            pgon = polyshape(verticesXY(:, 1), verticesXY(:, 2), 'Simplify', true);
            if isempty(pgon.Vertices)
                error('DataViewer:InvalidROIComponent', ...
                    'Mask component produced an invalid polyshape.');
            end

            ROI = templateROI;

            ROI.name = char(string(roiName));
            ROI.type = 'polygon';
            ROI.ID = roiID;
            ROI.modifiedOn = datetime('now');

            if ~isfield(ROI, 'DOC') || isempty(ROI.DOC)
                ROI.DOC = ROI.modifiedOn;
            end

            ROI.geometry = struct();
            ROI.geometry.polyshape = pgon;
            ROI.geometry.verticesXY_px = verticesXY;
            ROI.geometry.ROIType = 'polygon';
            ROI.geometry.ROIParameters = app.makePolygonROIParametersFromVertices(verticesXY);

            ROI.mask = logical(componentMask);
            ROI.stats = app.computeROIStatsFromMask(componentMask);

            visibleValue = true;
            if isfield(templateROI, 'runtime') && isfield(templateROI.runtime, 'visible')
                visibleValue = logical(templateROI.runtime.visible);
            end

            ROI.runtime = struct();
            ROI.runtime.visible = visibleValue;
            ROI.runtime.selected = false;
            ROI.runtime.ROIHandle = gobjects(1);
            ROI.runtime.editHandle = gobjects(1);
            ROI.runtime.traceLine = gobjects(1);
            ROI.runtime.tracePatch = gobjects(1);

            ROI.runtime.trace = struct();
            ROI.runtime.trace.XData = [];
            ROI.runtime.trace.Mean = [];
            ROI.runtime.trace.SEM = [];
            ROI.runtime.trace.Mode = '';
            ROI.runtime.trace.Status = 'not computed';

            if bCreateOverlay
                ROI.runtime.ROIHandle = app.createStaticROIOverlayFromROI(ROI);
            end

        end

        function nParts = replaceROIByMaskComponents(app, roiIdx, componentMasks, baseName, bUsePartNames)
            %REPLACEROIBYMASKCOMPONENTS Replace one ROI by one or more polygon components.
            %
            %   nParts = replaceROIByMaskComponents(app, roiIdx, componentMasks, baseName,
            %   bUsePartNames) replaces ROIList(roiIdx). If bUsePartNames is true, names
            %   are baseName_part<N>; otherwise the original baseName is preserved.

            nParts = 0;

            if isempty(roiIdx) || roiIdx < 1 || roiIdx > numel(app.ROIList)
                return
            end

            if isempty(componentMasks)
                return
            end

            templateROI = app.ROIList(roiIdx);
            app.deleteROIGraphicsByIndex(roiIdx);

            nComp = numel(componentMasks);
            newROIs = repmat(templateROI, 1, 0);

            existingIDs = [];
            if ~isempty(app.ROIList)
                existingIDs = [app.ROIList.ID];
                existingIDs = existingIDs(isfinite(existingIDs));
            end

            nextID = max([0 existingIDs]) + 1;

            for iComp = 1:nComp
                if bUsePartNames
                    requestedName = sprintf('%s_part%d', baseName, iComp);
                else
                    requestedName = baseName;
                end

                if iComp == 1 && isfield(templateROI, 'ID') && isfinite(templateROI.ID)
                    roiID = templateROI.ID;
                else
                    roiID = nextID;
                    nextID = nextID + 1;
                end

                % Avoid reusing the old ROI ID for newly appended parts.
                while any(existingIDs == roiID)
                    if iComp == 1 && roiID == templateROI.ID
                        break
                    end
                    roiID = nextID;
                    nextID = nextID + 1;
                end

                roiName = app.makeUniqueROIName(requestedName, roiID);

                try
                    newROI = app.makePolygonROIFromMaskComponent(templateROI, componentMasks{iComp}, roiName, roiID);
                catch ME
                    warning('DataViewer:SplitROICreateFailed', ...
                        'Skipping invalid split ROI component from "%s": %s', baseName, ME.message);
                    continue
                end

                if isempty(newROIs)
                    newROIs = newROI;
                else
                    newROIs(end+1) = newROI; %#ok<AGROW>
                end

                existingIDs(end+1) = roiID; %#ok<AGROW>
            end

            if isempty(newROIs)
                app.ROIList(roiIdx) = [];
                return
            end

            if roiIdx == 1
                beforeList = app.ROIList([]);
            else
                beforeList = app.ROIList(1:roiIdx-1);
            end

            if roiIdx == numel(app.ROIList)
                afterList = app.ROIList([]);
            else
                afterList = app.ROIList(roiIdx+1:end);
            end

            app.ROIList = [beforeList, newROIs, afterList];

            nParts = numel(newROIs);

        end

        function nAdded = addPolygonROIsFromMaskComponents(app, componentMasks, baseName, roiColor, bUsePartNames)
            %ADDPOLYGONROISFROMMASKCOMPONENTS Add one or more polygon ROIs from masks.
            %
            %   This is used when a newly drawn ROI is clipped by the logical mask into
            %   multiple connected components.

            nAdded = 0;

            if isempty(componentMasks)
                return
            end

            templateROI = struct();
            templateROI.name = char(string(baseName));
            templateROI.type = 'polygon';
            templateROI.DOC = datetime('now');
            templateROI.modifiedOn = templateROI.DOC;
            templateROI.color = roiColor;
            templateROI.notes = '';

            templateROI.geometry = struct();
            templateROI.mask = [];
            templateROI.stats = struct();

            templateROI.ID = NaN;
            templateROI.runtime = struct();
            templateROI.runtime.visible = true;
            templateROI.runtime.selected = false;
            templateROI.runtime.ROIHandle = gobjects(1);
            templateROI.runtime.editHandle = gobjects(1);
            templateROI.runtime.traceLine = gobjects(1);
            templateROI.runtime.tracePatch = gobjects(1);

            templateROI.runtime.trace = struct();
            templateROI.runtime.trace.XData = [];
            templateROI.runtime.trace.Mean = [];
            templateROI.runtime.trace.SEM = [];
            templateROI.runtime.trace.Mode = '';
            templateROI.runtime.trace.Status = 'not computed';

            nextID = app.getNextAvailableROIID();

            for iComp = 1:numel(componentMasks)
                if bUsePartNames
                    requestedName = sprintf('%s_part%d', baseName, iComp);
                else
                    requestedName = baseName;
                end

                roiID = nextID;
                nextID = nextID + 1;

                roiName = app.makeUniqueROIName(requestedName, roiID);

                try
                    ROI = app.makePolygonROIFromMaskComponent(templateROI, componentMasks{iComp}, roiName, roiID);
                catch ME
                    warning('DataViewer:SplitROICreateFailed', ...
                        'Skipping invalid split ROI component from "%s": %s', baseName, ME.message);
                    continue
                end

                if isempty(app.ROIList)
                    app.ROIList = ROI;
                else
                    app.ROIList(end+1) = ROI;
                end

                nAdded = nAdded + 1;
            end

        end

        function params = makePolygonROIParametersFromVertices(app, verticesXY) %#ok<INUSL>
            %MAKEPOLYGONROIPARAMETERSFROMVERTICES Build polygon ROIParameters struct.

            verticesXY = double(verticesXY);

            params = struct();
            params.ROIType = 'polygon';
            params.Position = verticesXY;
            params.Vertices = verticesXY;

        end

        function safeSetGraphicsContextMenu(~, hObj, hMenu)
            %SAFESETGRAPHICSCONTEXTMENU Assign a context menu only when figures match.
            %
            %   MATLAB errors if a UIContextMenu is parented to a different figure than
            %   the target graphics/component object. This helper skips invalid
            %   assignments and keeps drawing workflows from failing because of stale
            %   context-menu handles.

            if isempty(hObj) || isempty(hMenu)
                return
            end

            try
                if ~isvalid(hObj) || ~isvalid(hMenu)
                    return
                end
            catch
                return
            end

            try
                if ~isprop(hObj, 'ContextMenu')
                    return
                end

                objFig = ancestor(hObj, 'figure');
                menuFig = ancestor(hMenu, 'figure');

                if isempty(objFig) || isempty(menuFig) || ~isequal(objFig, menuFig)
                    return
                end

                hObj.ContextMenu = hMenu;
            catch
            end

        end

        function restoreGraphicsContextMenu(app, hObj, hMenu)
            %RESTOREGRAPHICSCONTEXTMENU Restore previous context menu if still compatible.
            %
            %   This is used after temporary logical-mask drawing. If the old menu is
            %   empty, the object context menu is cleared. If it is stale or parented to a
            %   different figure, the assignment is skipped safely.

            if isempty(hObj)
                return
            end

            try
                if ~isvalid(hObj) || ~isprop(hObj, 'ContextMenu')
                    return
                end
            catch
                return
            end

            if isempty(hMenu)
                try
                    hObj.ContextMenu = [];
                catch
                end
                return
            end

            app.safeSetGraphicsContextMenu(hObj, hMenu);

        end

        function applyROIBooleanOperation(app, operation)
            %APPLYROIBOOLEANOPERATION Create non-destructive Boolean ROI result.
            %
            %   operation:
            %       'intersect' - common pixels across all selected ROIs
            %       'merge'     - union of all selected ROIs
            %       'xor'       - pixels included by an odd number of selected ROIs
            %       'subtract'  - first selected ROI minus all other selected ROIs

            operation = lower(char(string(operation)));

            [roiIdxList, ok, msg] = app.getSelectedROIIndicesForBooleanOps();
            if ~ok
                app.setStatusMessage(msg);
                return
            end

            if ~app.selectedROIsHaveAnyOverlap(roiIdxList)
                app.setStatusMessage('Boolean operation cancelled: selected ROIs do not overlap.');
                return
            end

            [outMask, baseName, opLabel, ok, msg] = app.computeROIBooleanMask(roiIdxList, operation);
            if ~ok
                app.setStatusMessage(msg);
                return
            end

            outMask = app.clipROIMaskToActiveLogicalMask(outMask);

            if isempty(outMask) || ~any(outMask(:))
                app.setStatusMessage(sprintf('Boolean %s cancelled: result is empty.', opLabel));
                return
            end

            if app.roiMaskHasHoles(outMask)
                app.alertROIBooleanHolesUnsupported();
                app.setStatusMessage('Boolean operation cancelled: result contains holes.');
                return
            end

            componentMasks = app.maskToConnectedComponentMasks(outMask);

            if isempty(componentMasks)
                app.setStatusMessage(sprintf('Boolean %s cancelled: no valid output components.', opLabel));
                return
            end

            resultName = app.makeUniqueROIName(baseName, []);
            roiColor = app.getNextROIColorForBooleanResult(roiIdxList);

            nAdded = app.addPolygonROIsFromMaskComponents( ...
                componentMasks, ...
                resultName, ...
                roiColor, ...
                numel(componentMasks) > 1);

            if nAdded < 1
                app.setStatusMessage(sprintf('Boolean %s cancelled: output polygon conversion failed.', opLabel));
                return
            end

            app.clearROISelectionState();
            app.stopActiveROIEditForSelectionChange();

            app.refreshROITable();
            app.refreshROITraces();
            app.refreshEventPatches();
            app.stackROITraceGraphics();
            app.updateGUIEnabledState();

            if nAdded == 1
                app.setStatusMessage(sprintf('Created Boolean ROI "%s" using %s.', resultName, opLabel));
            else
                app.setStatusMessage(sprintf('Created Boolean ROI "%s" as %d parts using %s.', resultName, nAdded, opLabel));
            end

        end

        function [roiIdxList, ok, msg] = getSelectedROIIndicesForBooleanOps(app)
            %GETSELECTEDROIINDICESFORBOOLEANOPS Return selected ROI indices for Boolean Ops.

            roiIdxList = [];
            ok = false;
            msg = '';

            if isempty(app.ROIList)
                msg = 'Boolean operation cancelled: no ROIs available.';
                return
            end

            selectedIDs = [];

            if isfield(app.GroupROIEditState, 'roiIDList') && ~isempty(app.GroupROIEditState.roiIDList)
                selectedIDs = double(app.GroupROIEditState.roiIDList(:).');
            else
                selectedIDs = app.getSelectedROIIDList();
            end

            selectedIDs = unique(double(selectedIDs(:).'), 'stable');
            selectedIDs = selectedIDs(isfinite(selectedIDs));

            if numel(selectedIDs) < 2
                msg = 'Boolean operation cancelled: select at least two ROIs.';
                return
            end

            for iID = 1:numel(selectedIDs)
                roiIdx = app.getROIIndexByID(selectedIDs(iID));

                if ~isempty(roiIdx)
                    roiIdxList(end+1) = roiIdx; %#ok<AGROW>
                end
            end

            roiIdxList = unique(roiIdxList, 'stable');

            if numel(roiIdxList) < 2
                msg = 'Boolean operation cancelled: at least two valid selected ROIs are required.';
                return
            end

            sz = app.getDataSize();
            expectedSize = sz(1:2);

            for iROI = 1:numel(roiIdxList)
                roiIdx = roiIdxList(iROI);

                if ~isfield(app.ROIList(roiIdx), 'mask') || isempty(app.ROIList(roiIdx).mask) || ...
                        ~isequal(size(app.ROIList(roiIdx).mask), expectedSize)
                    msg = sprintf('Boolean operation cancelled: ROI "%s" has an invalid mask.', app.ROIList(roiIdx).name);
                    return
                end
            end

            ok = true;

        end

        function tf = selectedROIsHaveAnyOverlap(app, roiIdxList)
            %SELECTEDROISHAVEANYOVERLAP True when at least two selected ROI masks overlap.

            tf = false;

            if numel(roiIdxList) < 2
                return
            end

            sz = app.getDataSize();
            Ny = sz(1);
            Nx = sz(2);

            countMap = zeros(Ny, Nx, 'uint16');

            for iROI = 1:numel(roiIdxList)
                mask = logical(app.ROIList(roiIdxList(iROI)).mask);

                if ~isequal(size(mask), [Ny Nx])
                    continue
                end

                countMap(mask) = countMap(mask) + 1;

                if any(countMap(:) >= 2)
                    tf = true;
                    return
                end
            end

        end

        function [outMask, baseName, opLabel, ok, msg] = computeROIBooleanMask(app, roiIdxList, operation)
            %COMPUTEROIBOOLEANMASK Compute Boolean operation from selected ROI masks.

            outMask = [];
            baseName = '';
            opLabel = operation;
            ok = false;
            msg = '';

            if numel(roiIdxList) < 2
                msg = 'Boolean operation cancelled: select at least two ROIs.';
                return
            end

            masks = cell(1, numel(roiIdxList));

            for iROI = 1:numel(roiIdxList)
                masks{iROI} = logical(app.ROIList(roiIdxList(iROI)).mask);
            end

            firstName = char(string(app.ROIList(roiIdxList(1)).name));

            switch lower(operation)
                case {'intersect', 'intersection'}
                    outMask = masks{1};
                    for iROI = 2:numel(masks)
                        outMask = outMask & masks{iROI};
                    end

                    baseName = sprintf('%s_intersect', firstName);
                    opLabel = 'intersection';

                case {'merge', 'union'}
                    outMask = false(size(masks{1}));
                    for iROI = 1:numel(masks)
                        outMask = outMask | masks{iROI};
                    end

                    baseName = sprintf('%s_merge', firstName);
                    opLabel = 'merge';

                case 'xor'
                    outMask = false(size(masks{1}));
                    for iROI = 1:numel(masks)
                        outMask = xor(outMask, masks{iROI});
                    end

                    baseName = sprintf('%s_xor', firstName);
                    opLabel = 'XOR';

                case {'subtract', 'minus'}
                    baseMask = masks{1};
                    subtractMask = false(size(baseMask));

                    for iROI = 2:numel(masks)
                        subtractMask = subtractMask | masks{iROI};
                    end

                    outMask = baseMask & ~subtractMask;
                    baseName = sprintf('%s_subtract', firstName);
                    opLabel = sprintf('subtract from %s', firstName);

                otherwise
                    msg = sprintf('Unsupported Boolean operation: %s.', operation);
                    return
            end

            outMask = logical(outMask);
            ok = true;

        end

        function tf = roiMaskHasHoles(app, mask)
            %ROIMASKHASHOLES True when a mask contains internal holes.
            %
            %   Hole-containing Boolean results are rejected for now because the current
            %   ROI overlay/editing model does not display holes explicitly. This prevents
            %   a mismatch between what the user sees and what ROI.mask uses for analysis.

            tf = false;

            if isempty(mask) || ~any(mask(:))
                return
            end

            componentMasks = app.maskToConnectedComponentMasks(mask);

            for iComp = 1:numel(componentMasks)
                compMask = logical(componentMasks{iComp});

                try
                    filledMask = imfill(compMask, 'holes');
                catch
                    filledMask = compMask;
                end

                holeMask = filledMask & ~compMask;

                if any(holeMask(:))
                    tf = true;
                    return
                end
            end

        end

        function alertROIBooleanHolesUnsupported(app)
            %ALERTROIBOOLEANHOLESUNSUPPORTED Warn that hole-containing ROIs are rejected.

            msg = [ ...
                'The Boolean result contains internal holes. ' ...
                'Creating this ROI would make the visual overlay differ from the analysis mask. ' ...
                'Operation cancelled.'];

            try
                uialert(app.UIFigure, msg, ...
                    'Hole-containing ROI not supported', ...
                    'Icon', 'warning');
            catch
                warning('DataViewer:ROIBooleanHolesUnsupported', '%s', msg);
            end

        end

        function roiColor = getNextROIColorForBooleanResult(app, roiIdxList)
            %GETNEXTROICOLORFORBOOLEANRESULT Pick display color for Boolean output ROI.

            roiColor = [];

            if ~isempty(roiIdxList) && roiIdxList(1) >= 1 && roiIdxList(1) <= numel(app.ROIList) && ...
                    isfield(app.ROIList(roiIdxList(1)), 'color') && ...
                    numel(app.ROIList(roiIdxList(1)).color) == 3
                roiColor = double(app.ROIList(roiIdxList(1)).color(:).');
            end

            if isempty(roiColor) || any(~isfinite(roiColor))
                try
                    roiColor = app.ROIColorList(app.ROINextColorIndex, :);
                catch
                    roiColor = [1 0 0];
                end
            end

            roiColor = min(max(roiColor, 0), 1);

        end

        function updateROISelectionOrder(app, roiID, bSelected)
            %UPDATEROISELECTIONORDER Track chronological ROI Select-checkbox order.
            %
            %   This is used by order-dependent Boolean Ops. In particular:
            %       subtract => first selected ROI - all later selected ROIs.

            if isempty(roiID) || ~isfinite(roiID)
                return
            end

            roiID = double(roiID);
            app.ROISelectionOrder = double(app.ROISelectionOrder(:).');
            app.ROISelectionOrder = app.ROISelectionOrder(isfinite(app.ROISelectionOrder));

            % Remove stale duplicate entries first.
            app.ROISelectionOrder(app.ROISelectionOrder == roiID) = [];

            if bSelected
                app.ROISelectionOrder(end+1) = roiID;
            end

            % Keep only currently existing ROI IDs.
            if ~isempty(app.ROIList)
                validIDs = [app.ROIList.ID];
                app.ROISelectionOrder = app.ROISelectionOrder(ismember(app.ROISelectionOrder, validIDs));
            else
                app.ROISelectionOrder = [];
            end

        end

        function createEditableROIContextMenu(app, hEdit, roiID)
            %CREATEEDITABLEROICONTEXTMENU Modify built-in ROI edit context menu.
            %
            %   MATLAB ROI objects already create useful built-in menu items. Do not replace
            %   the full menu. Instead, rename the built-in Delete <shape> item to
            %   "Cancel ROI editing", then append DataViewer-specific actions.

            if ~app.isUsableGraphicsHandle(hEdit)
                return
            end

            roiIdx = app.getROIIndexByID(roiID);
            if isempty(roiIdx)
                return
            end

            ctx = app.getOrCreateROIObjectContextMenu(hEdit);

            if isempty(ctx) || ~isvalid(ctx)
                return
            end

            % Replace the built-in "Delete <shape>" action with cancellation. This preserves
            % other built-in ROI menu actions such as label/color/options when available.
            app.renameBuiltInROIDeleteMenuItem( ...
                ctx, ...
                'Cancel ROI editing', ...
                @(src, evt) app.cancelROIEditingByIDFromContextMenu(roiID));

            % Add/refresh app-specific edit actions.
            app.deleteMenuChildrenByTag(ctx, 'DataViewer_EditDeleteROI');
            app.deleteMenuChildrenByTag(ctx, 'DataViewer_EditConstrainedROI');

            deleteMenu = uimenu(ctx);
            deleteMenu.Text = 'Delete ROI';
            deleteMenu.Tag = 'DataViewer_EditDeleteROI';
            deleteMenu.Separator = 'on';
            deleteMenu.MenuSelectedFcn = @(src, evt) app.deleteROIByIDFromEditContextMenu(roiID);

            constrainedMenu = uimenu(ctx);
            constrainedMenu.Text = 'Constrained edit';
            constrainedMenu.Tag = 'DataViewer_EditConstrainedROI';
            constrainedMenu.MenuSelectedFcn = @(src, evt) app.startConstrainedSingleROIEditByID(roiID);

            app.ROIList(roiIdx).runtime.editContextMenu = ctx;

        end

        function cancelROIEditingByIDFromContextMenu(app, roiID)
            %CANCELROIEDITINGBYIDFROMCONTEXTMENU Cancel active single-ROI edit by ROI ID.

            roiIdx = app.getROIIndexByID(roiID);

            if isempty(roiIdx)
                return
            end

            app.cancelROIEditingByIndex(roiIdx, false);

        end

        function deleteROIByIDFromEditContextMenu(app, roiID)
            %DELETEROIBYIDFROMEDITCONTEXTMENU Delete active single ROI after confirmation.

            roiIdx = app.getROIIndexByID(roiID);

            if isempty(roiIdx)
                return
            end

            roiName = char(string(app.ROIList(roiIdx).name));
            promptText = sprintf('Delete ROI "%s"?', roiName);

            if ~app.confirmROIDeletion(promptText)
                app.setStatusMessage('ROI deletion cancelled.');
                return
            end

            app.deleteROIsByIDList(roiID, true);

        end

        function startConstrainedSingleROIEditByID(app, roiID)
            %STARTCONSTRAINEDSINGLEROIEDITBYID Edit one ROI with a transform bounding box.
            %
            %   This mode is intended for complex polygons. The user manipulates a single
            %   bounding rectangle that translates, scales, and rotates the ROI as a whole.
            %   The final commit reuses commitGroupROIEdit, so logical-mask clipping,
            %   invalid-ROI deletion, split ROI handling, and primitive approximation remain
            %   consistent with grouped ROI editing.

            roiIdx = app.getROIIndexByID(roiID);

            if isempty(roiIdx) || roiIdx < 1 || roiIdx > numel(app.ROIList)
                return
            end

            if ~app.hasData()
                return
            end

            % Remove active free-edit object if this command was launched from the
            % single-ROI edit context menu.
            if isfield(app.ROIList(roiIdx), 'runtime') && ...
                    isfield(app.ROIList(roiIdx).runtime, 'editHandle') && ...
                    app.isUsableGraphicsHandle(app.ROIList(roiIdx).runtime.editHandle)
                app.cleanupROIEditRuntimeByIndex(roiIdx);
            end

            app.deleteGroupEditRuntimeGraphics();

            [state, ok, msg] = app.buildGroupROIEditState(roiIdx);
            if ~ok
                app.setStatusMessage(msg);
                return
            end

            previousKeyFcn = app.UIFigure.WindowKeyPressFcn;
            state.previousKeyFcn = previousKeyFcn;
            state.previousMode = app.InteractionMode;
            state.editMode = 'singleConstrained';

            if isfield(app.ROIList(roiIdx).runtime, 'ROIHandle') && ...
                    app.isUsableGraphicsHandle(app.ROIList(roiIdx).runtime.ROIHandle)
                try
                    delete(app.ROIList(roiIdx).runtime.ROIHandle);
                catch
                end
            end

            app.ROIList(roiIdx).runtime.ROIHandle = gobjects(1);
            app.setROISelectedStateByIndex(roiIdx, true);

            bbox = state.originalBBox;

            hold(app.ImageAxes, 'on');

            hBox = drawrectangle(app.ImageAxes, ...
                'Position', bbox, ...
                'Color', [1 1 1], ...
                'FaceAlpha', 0, ...
                'LineWidth', 1.5, ...
                'Rotatable', true, ...
                'DrawingArea', 'unlimited', ...
                'InteractionsAllowed', 'all');

            app.setROIObjectPropertyIfAvailable(hBox, 'LineStyle', '--');
            app.updateGroupEditBoxSizeLabel(hBox);

            listeners = {};
            listeners{end+1} = addlistener(hBox, 'MovingROI', ...
                @(src, evt) app.onGroupEditBoxMoving(src, evt));
            listeners{end+1} = addlistener(hBox, 'ROIMoved', ...
                @(src, evt) app.onGroupEditBoxMoved(src, evt));
            listeners{end+1} = addlistener(hBox, 'ROIClicked', ...
                @(src, evt) app.onGroupEditBoxClicked(src, evt));

            try
                listeners{end+1} = addlistener(hBox, 'DeletingROI', ...
                    @(src, evt) app.onGroupEditBoxDeleting(src, evt));
            catch
            end

            state.listeners = listeners;
            app.GroupROIEditState = state;
            app.GroupEditBoxHandle = hBox;
            app.UIFigure.WindowKeyPressFcn = @(src, evt) app.onActiveGroupROIEditKeyPress(src, evt);

            app.createGroupROIEditContextMenu();
            app.createGroupEditPreviewGraphics();
            app.updateGroupEditPreview();

            app.setInteractionMode('editingROI');
            app.refreshROITable();
            app.setStatusMessage(sprintf( ...
                ['Constrained editing ROI "%s". Move/resize/rotate the box; ' ...
                'double-click or press Enter to confirm; press Escape or use the context menu to cancel.'], ...
                app.ROIList(roiIdx).name));

        end

        function createGroupROIEditContextMenu(app)
            %CREATEGROUPROIEDITCONTEXTMENU Create context menu for group edit box.
            %
            %   Menu order:
            %       Delete ROI(s)
            %       Cancel
            %       separator
            %       Boolean Ops
            %
            %   The menu is attached to the temporary group edit rectangle and destroyed
            %   with the group edit runtime.

            if ~app.isUsableGraphicsHandle(app.GroupEditBoxHandle)
                return
            end

            parentFig = ancestor(app.GroupEditBoxHandle, 'figure');
            if isempty(parentFig) || ~isvalid(parentFig)
                parentFig = app.UIFigure;
            end

            app.deleteGroupROIEditContextMenu();

            app.GroupROIEditContextMenu = uicontextmenu(parentFig);

            app.GroupROIDeleteMenu = uimenu(app.GroupROIEditContextMenu);
            if app.isSingleConstrainedGroupEdit()
                app.GroupROIDeleteMenu.Text = 'Delete ROI';
            else
                app.GroupROIDeleteMenu.Text = 'Delete ROIs';
            end
            app.GroupROIDeleteMenu.MenuSelectedFcn = @(src, evt) app.deleteGroupEditedROIsFromContextMenu();

            app.GroupROICancelMenu = uimenu(app.GroupROIEditContextMenu);
            app.GroupROICancelMenu.Text = 'Cancel';
            app.GroupROICancelMenu.MenuSelectedFcn = @(src, evt) app.cancelGroupROIEdit(false);

            app.GroupROIBooleanOpsMenu = uimenu(app.GroupROIEditContextMenu);
            app.GroupROIBooleanOpsMenu.Text = 'Boolean Ops';
            app.GroupROIBooleanOpsMenu.Separator = 'on';

            app.GroupROIIntersectMenu = uimenu(app.GroupROIBooleanOpsMenu);
            app.GroupROIIntersectMenu.Text = 'Intersect';
            app.GroupROIIntersectMenu.MenuSelectedFcn = @(src, evt) app.applyROIBooleanOperation('intersect');

            app.GroupROIMergeMenu = uimenu(app.GroupROIBooleanOpsMenu);
            app.GroupROIMergeMenu.Text = 'Merge';
            app.GroupROIMergeMenu.MenuSelectedFcn = @(src, evt) app.applyROIBooleanOperation('merge');

            app.GroupROIXORMenu = uimenu(app.GroupROIBooleanOpsMenu);
            app.GroupROIXORMenu.Text = 'XOR';
            app.GroupROIXORMenu.MenuSelectedFcn = @(src, evt) app.applyROIBooleanOperation('xor');

            app.GroupROISubtractMenu = uimenu(app.GroupROIBooleanOpsMenu);
            app.GroupROISubtractMenu.Text = 'Subtract selected ROIs from first selected ROI';
            app.GroupROISubtractMenu.MenuSelectedFcn = @(src, evt) app.applyROIBooleanOperation('subtract');

            % Boolean Ops require at least two selected ROIs. Keep the menu visible but
            % disabled in constrained single-ROI edit mode.
            if numel(app.getActiveGroupEditROIIDList()) < 2
                app.GroupROIBooleanOpsMenu.Enable = 'off';
            else
                app.GroupROIBooleanOpsMenu.Enable = 'on';
            end

            app.safeSetGraphicsContextMenu(app.GroupEditBoxHandle, app.GroupROIEditContextMenu);

        end

        function deleteGroupROIEditContextMenu(app)
            %DELETEGROUPROIEDITCONTEXTMENU Delete temporary group edit context menu.

            try
                if ~isempty(app.GroupROIEditContextMenu) && isvalid(app.GroupROIEditContextMenu)
                    delete(app.GroupROIEditContextMenu);
                end
            catch
            end

            app.GroupROIEditContextMenu = matlab.ui.container.ContextMenu.empty;
            app.GroupROIDeleteMenu = matlab.ui.container.Menu.empty;
            app.GroupROICancelMenu = matlab.ui.container.Menu.empty;
            app.GroupROIBooleanOpsMenu = matlab.ui.container.Menu.empty;
            app.GroupROIIntersectMenu = matlab.ui.container.Menu.empty;
            app.GroupROIMergeMenu = matlab.ui.container.Menu.empty;
            app.GroupROIXORMenu = matlab.ui.container.Menu.empty;
            app.GroupROISubtractMenu = matlab.ui.container.Menu.empty;

        end

        function deleteGroupEditedROIsFromContextMenu(app)
            %DELETEGROUPEDITEDROISFROMCONTEXTMENU Delete ROIs represented by group edit box.

            roiIDList = app.getActiveGroupEditROIIDList();
            roiIDList = unique(double(roiIDList(:).'), 'stable');
            roiIDList = roiIDList(isfinite(roiIDList));

            if isempty(roiIDList)
                app.setStatusMessage('No selected ROIs to delete.');
                return
            end

            if numel(roiIDList) == 1
                roiIdx = app.getROIIndexByID(roiIDList(1));
                if isempty(roiIdx)
                    app.setStatusMessage('Selected ROI was not found.');
                    return
                end
                promptText = sprintf('Delete ROI "%s"?', app.ROIList(roiIdx).name);
            else
                promptText = sprintf('Delete %d selected ROIs?', numel(roiIDList));
            end

            if ~app.confirmROIDeletion(promptText)
                app.setStatusMessage('ROI deletion cancelled.');
                return
            end

            app.deleteROIsByIDList(roiIDList, true);

        end

        function tf = isSingleConstrainedGroupEdit(app)
            %ISSINGLECONSTRAINEDGROUPEDIT True for constrained single-ROI box edit mode.

            tf = false;

            try
                tf = isfield(app.GroupROIEditState, 'editMode') && ...
                    strcmpi(app.GroupROIEditState.editMode, 'singleConstrained');
            catch
                tf = false;
            end

        end

        function ctx = getOrCreateROIObjectContextMenu(app, hROI)
            %GETORCREATEROIOBJECTCONTEXTMENU Return an ROI object's built-in context menu.
            %
            %   MATLAB ROI objects normally expose a built-in context menu. This helper
            %   retrieves that menu when available. If no menu exists, it creates one as a
            %   fallback so the app still has a place for cancellation/actions.

            ctx = [];

            if ~app.isUsableGraphicsHandle(hROI)
                return
            end

            try
                if isprop(hROI, 'ContextMenu') && ~isempty(hROI.ContextMenu) && ...
                        isvalid(hROI.ContextMenu)
                    ctx = hROI.ContextMenu;
                    return
                end
            catch
                ctx = [];
            end

            try
                parentFig = ancestor(hROI, 'figure');
                if isempty(parentFig) || ~isvalid(parentFig)
                    parentFig = app.UIFigure;
                end

                ctx = uicontextmenu(parentFig);
                app.safeSetGraphicsContextMenu(hROI, ctx);
            catch
                ctx = [];
            end

        end

        function renameBuiltInROIDeleteMenuItem(app, ctx, newText, callbackFcn)
            %RENAMEBUILTINROIDELETEMENUITEM Rename built-in Delete <shape> ROI menu item.
            %
            %   The goal is to preserve MATLAB's other built-in ROI context-menu items while
            %   replacing only the destructive Delete <shape> action with a DataViewer cancel
            %   action.

            if isempty(ctx) || ~isvalid(ctx)
                return
            end

            deleteMenu = app.findROIDeleteMenuItem(ctx);

            if isempty(deleteMenu) || ~isvalid(deleteMenu)
                deleteMenu = uimenu(ctx);
            end

            deleteMenu.Text = char(string(newText));
            deleteMenu.Tag = 'DataViewer_RenamedBuiltInDelete';
            deleteMenu.MenuSelectedFcn = callbackFcn;

        end

        function hMenu = findROIDeleteMenuItem(app, ctx) %#ok<INUSD>
            %FINDROIDELETEMENUITEM Find built-in ROI Delete <shape> menu item.

            hMenu = [];

            if isempty(ctx) || ~isvalid(ctx)
                return
            end

            try
                children = findall(ctx, 'Type', 'uimenu');
            catch
                children = matlab.ui.container.Menu.empty;
            end

            for iChild = 1:numel(children)
                try
                    menuText = strtrim(char(string(children(iChild).Text)));
                    menuTag = '';
                    if isprop(children(iChild), 'Tag')
                        menuTag = char(string(children(iChild).Tag));
                    end

                    if strcmp(menuTag, 'DataViewer_RenamedBuiltInDelete')
                        hMenu = children(iChild);
                        return
                    end

                    % Built-in ROI menus are usually "Delete Rectangle", "Delete Ellipse",
                    % or "Delete Polygon". Avoid matching app-specific delete actions.
                    if startsWith(lower(menuText), 'delete ') && ...
                            ~contains(lower(menuText), 'roi') && ...
                            ~contains(lower(menuText), 'selected')
                        hMenu = children(iChild);
                        return
                    end
                catch
                end
            end

        end

        function deleteMenuChildrenByTag(app, ctx, tagText) %#ok<INUSD>
            %DELETEMENUCHILDRENBYTAG Delete app-added menu children with a specific tag.

            if isempty(ctx) || ~isvalid(ctx)
                return
            end

            try
                children = findall(ctx, 'Type', 'uimenu', 'Tag', char(string(tagText)));
            catch
                children = matlab.ui.container.Menu.empty;
            end

            for iChild = 1:numel(children)
                try
                    if isvalid(children(iChild))
                        delete(children(iChild));
                    end
                catch
                end
            end

        end

        function modifyROICreationContextMenu(app, hROI, cancelFcn)
            %MODIFYROICREATIONCONTEXTMENU Preserve built-ins and replace Delete by Cancel.
            %
            %   During ROI creation, MATLAB's built-in context menu may contain useful
            %   editing options. Only rename the Delete <shape> action to "Cancel ROI
            %   creation" and redirect it to the supplied cancel function.

            if ~app.isUsableGraphicsHandle(hROI)
                return
            end

            ctx = app.getOrCreateROIObjectContextMenu(hROI);

            if isempty(ctx) || ~isvalid(ctx)
                return
            end

            app.renameBuiltInROIDeleteMenuItem(ctx, 'Cancel ROI creation', cancelFcn);

        end

        function saveROIsToFile(app)
            %SAVEROISTOFILE Save current ROIs to a UMIT .roi MAT-file.
            %
            %   The saved file is a MAT-file with extension ".roi". Runtime graphics handles
            %   and listeners are stripped. Geometry, masks, stats, colors, names, image size,
            %   and source metadata are preserved.

            if isempty(app.ROIList)
                app.setStatusMessage('No ROIs available to save.');
                return
            end

            if ~app.hasData()
                app.setStatusMessage('Load image data before saving ROIs.');
                return
            end

            startFolder = pwd;
            if ~isempty(app.CurrentFile)
                startFolder = fileparts(app.CurrentFile);
            end

            defaultName = app.makeDefaultROIFileName();

            [fileName, folderName] = uiputfile( ...
                {'*.roi', 'UMIT ROI file (*.roi)'}, ...
                'Save ROIs', ...
                fullfile(startFolder, defaultName));

            if isequal(fileName, 0)
                app.setStatusMessage('Save ROI cancelled.');
                return
            end

            [~, ~, ext] = fileparts(fileName);
            if isempty(ext)
                fileName = [fileName, '.roi'];
            end

            filePath = fullfile(folderName, fileName);

            ROIFile = app.buildROIFileStruct();

            try
                save(filePath, 'ROIFile', '-v7.3');
                app.setStatusMessage(sprintf('Saved %d ROIs to "%s".', numel(app.ROIList), filePath));
            catch ME
                app.setStatusMessage(sprintf('Save ROI failed: %s', ME.message));
                rethrow(ME)
            end

        end

        function insertMode = promptROILoadInsertMode(app, nIncoming)
            %PROMPTROILOADINSERTMODE Ask whether loaded ROIs replace or append current ROIs.
            %
            %   insertMode:
            %       'replace' - remove current ROIs and load the file ROIs
            %       'append'  - preserve current ROIs and append loaded ROIs
            %       'cancel'  - cancel loading

            insertMode = 'replace';

            if isempty(app.ROIList)
                return
            end

            promptText = sprintf([ ...
                'The current viewer already contains %d ROI(s).\n\n' ...
                'The selected file contains %d source ROI(s).\n\n' ...
                'Choose how to load the ROIs.'], ...
                numel(app.ROIList), nIncoming);

            try
                choice = uiconfirm(app.UIFigure, promptText, ...
                    'Load ROIs', ...
                    'Options', {'Replace current ROIs', 'Append to current ROIs', 'Cancel'}, ...
                    'DefaultOption', 'Append to current ROIs', ...
                    'CancelOption', 'Cancel', ...
                    'Icon', 'question');
            catch
                choice = questdlg(promptText, ...
                    'Load ROIs', ...
                    'Replace current ROIs', 'Append to current ROIs', 'Cancel', ...
                    'Append to current ROIs');

                if isempty(choice)
                    choice = 'Cancel';
                end
            end

            switch choice
                case 'Replace current ROIs'
                    insertMode = 'replace';

                case 'Append to current ROIs'
                    insertMode = 'append';

                otherwise
                    insertMode = 'cancel';
            end

        end

        function loadROIsFromFile(app)
            %LOADROISFROMFILE Load ROIs from a UMIT .roi file.
            %
            %   If current ROIs already exist, the user can replace them, append the
            %   loaded ROIs, or cancel. When appending, loaded ROI names and IDs are
            %   resolved against the current ROIList to avoid collisions.
            %
            %   Loaded ROIs are always forced to:
            %       runtime.visible  = true
            %       runtime.selected = false

            if ~app.hasData()
                app.setStatusMessage('Load image data before loading ROIs.');
                return
            end

            startFolder = pwd;
            if ~isempty(app.CurrentFile)
                startFolder = fileparts(app.CurrentFile);
            end

            [fileName, folderName] = uigetfile( ...
                {'*.roi', 'UMIT ROI file (*.roi)'}, ...
                'Load ROIs', ...
                startFolder);

            if isequal(fileName, 0)
                app.setStatusMessage('Load ROI cancelled.');
                return
            end

            filePath = fullfile(folderName, fileName);

            try
                ROIFile = app.readROIFileStruct(filePath);
                app.validateROIFileStruct(ROIFile);

                [loadAction, scaleFactor] = app.promptROILoadSizeMismatchAction(ROIFile);

                if strcmpi(loadAction, 'cancel')
                    app.setStatusMessage('Load ROI cancelled.');
                    return
                end

                insertMode = app.promptROILoadInsertMode(numel(ROIFile.ROIs));
                if strcmpi(insertMode, 'cancel')
                    app.setStatusMessage('Load ROI cancelled.');
                    return
                end

                [baseUsedIDs, baseUsedNames] = app.getExistingROIIdentifiersForLoad(insertMode);

                [newROIList, report] = app.rebuildLoadedROIListForCurrentImage( ...
                    ROIFile.ROIs, loadAction, scaleFactor, baseUsedIDs, baseUsedNames);

                if isempty(newROIList)
                    app.setStatusMessage('No valid ROIs were loaded from the selected file.');
                    return
                end

                switch lower(insertMode)
                    case 'replace'
                        app.clearAllROIsWithoutPrompt();
                        app.ROIList = newROIList;

                    case 'append'
                        app.stopActiveROIEditForSelectionChange();
                        app.ROIList = app.appendROIsToROIList(app.ROIList, newROIList);

                    otherwise
                        app.setStatusMessage('Load ROI cancelled.');
                        return
                end

                app.SelectedROIID = NaN;
                app.ROISelectionOrder = [];

                % Build display handles only after ROIList is assigned. This avoids
                % creating transient graphics for loaded ROIs that may later be skipped
                % or appended with adjusted struct fields.
                app.rebuildMissingROIOverlaysAfterLoad();

                if ~isempty(app.ROIColorList)
                    app.ROINextColorIndex = mod(numel(app.ROIList), size(app.ROIColorList, 1)) + 1;
                end

                app.refreshROITable();
                app.refreshROITraces();
                app.refreshEventPatches();
                app.stackROITraceGraphics();
                app.updateGUIEnabledState();

                switch lower(insertMode)
                    case 'replace'
                        statusText = sprintf('Loaded %d ROI(s) from "%s".', numel(newROIList), filePath);
                    case 'append'
                        statusText = sprintf('Appended %d ROI(s) from "%s".', numel(newROIList), filePath);
                end

                if strcmpi(loadAction, 'scale')
                    statusText = sprintf('%s Scaled by %.6g.', statusText, scaleFactor);
                end
                if report.nClipped > 0
                    statusText = sprintf('%s %d source ROI(s) clipped to logical mask.', statusText, report.nClipped);
                end
                if report.nSplit > 0
                    statusText = sprintf('%s %d source ROI(s) split into parts.', statusText, report.nSplit);
                end
                if report.nSkipped > 0
                    statusText = sprintf('%s %d invalid/outside source ROI(s) skipped.', statusText, report.nSkipped);
                end

                app.setStatusMessage(statusText);

            catch ME
                app.setStatusMessage(sprintf('Load ROI failed: %s', ME.message));
                rethrow(ME)
            end

        end

        function [baseUsedIDs, baseUsedNames] = getExistingROIIdentifiersForLoad(app, insertMode)
            %GETEXISTINGROIIDENTIFIERSFORLOAD Return IDs/names to avoid during ROI load.
            %
            %   For replace mode, only the loaded ROIs are checked against each other.
            %   For append mode, loaded ROIs are also checked against current ROI names/IDs.

            baseUsedIDs = [];
            baseUsedNames = {};

            if ~strcmpi(insertMode, 'append') || isempty(app.ROIList)
                return
            end

            for iROI = 1:numel(app.ROIList)
                if isfield(app.ROIList(iROI), 'ID') && isfinite(app.ROIList(iROI).ID)
                    baseUsedIDs(end+1) = double(app.ROIList(iROI).ID); %#ok<AGROW>
                end

                if isfield(app.ROIList(iROI), 'name') && ~isempty(app.ROIList(iROI).name)
                    baseUsedNames{end+1} = char(string(app.ROIList(iROI).name)); %#ok<AGROW>
                end
            end

            baseUsedIDs = unique(baseUsedIDs, 'stable');

        end

        function [newROIList, report] = rebuildLoadedROIListForCurrentImage(app, savedROIList, loadAction, scaleFactor, baseUsedIDs, baseUsedNames)
            %REBUILDLOADEDROILISTFORCURRENTIMAGE Rebuild loaded ROIs for current image.
            %
            %   Saved masks are intentionally not trusted. Geometry is the durable
            %   representation. Geometry is optionally scaled, rasterized against the
            %   current image size, clipped to the active logical mask, and converted back
            %   to ROI structs with fresh runtime state. Graphics handles are rebuilt
            %   after ROIList is assigned.
            %
            %   baseUsedIDs/baseUsedNames allow append mode to resolve collisions against
            %   current ROIs before adding loaded ROIs.

            if nargin < 5 || isempty(baseUsedIDs)
                baseUsedIDs = [];
            end

            if nargin < 6 || isempty(baseUsedNames)
                baseUsedNames = {};
            end

            templateROI = app.makeLoadedROITemplate();
            newROIList = repmat(templateROI, 1, 0);

            report = struct();
            report.nSkipped = 0;
            report.nClipped = 0;
            report.nSplit = 0;

            if isempty(savedROIList)
                return
            end

            usedIDs = double(baseUsedIDs(:).');
            usedIDs = usedIDs(isfinite(usedIDs));
            usedNames = cellstr(string(baseUsedNames(:)));

            for iROI = 1:numel(savedROIList)
                ROI = savedROIList(iROI);

                if strcmpi(loadAction, 'scale')
                    ROI = app.scaleLoadedROIForCurrentImage(ROI, scaleFactor);
                end

                [roiOutList, status] = app.rebuildOneLoadedROIForCurrentImage(ROI, usedIDs, usedNames);

                if isempty(roiOutList)
                    report.nSkipped = report.nSkipped + 1;
                    continue
                end

                if status.bClipped
                    report.nClipped = report.nClipped + 1;
                end

                if status.bSplit
                    report.nSplit = report.nSplit + 1;
                end

                for k = 1:numel(roiOutList)
                    newROIList(end+1) = roiOutList(k); %#ok<AGROW>
                    usedIDs(end+1) = roiOutList(k).ID; %#ok<AGROW>
                    usedNames{end+1} = roiOutList(k).name; %#ok<AGROW>
                end
            end

        end

        function roiListOut = appendROIsToROIList(app, currentROIList, roiListToAppend) %#ok<INUSL>
            %APPENDROISTOROILIST Append ROI structs after aligning field sets.
            %
            %   Older ROI structs may not have the exact same fields as newly loaded
            %   ROIs. MATLAB requires all elements of a struct array to have identical
            %   fields before concatenation, so field alignment is performed inline here.

            if isempty(currentROIList)
                roiListOut = roiListToAppend;
                return
            end

            if isempty(roiListToAppend)
                roiListOut = currentROIList;
                return
            end

            fieldsA = fieldnames(currentROIList);
            fieldsB = fieldnames(roiListToAppend);
            allFields = unique([fieldsA; fieldsB], 'stable');

            for iField = 1:numel(allFields)
                fieldName = allFields{iField};

                if ~isfield(currentROIList, fieldName)
                    for iROI = 1:numel(currentROIList)
                        currentROIList(iROI).(fieldName) = [];
                    end
                end

                if ~isfield(roiListToAppend, fieldName)
                    for iROI = 1:numel(roiListToAppend)
                        roiListToAppend(iROI).(fieldName) = [];
                    end
                end
            end

            currentROIList = orderfields(currentROIList, allFields);
            roiListToAppend = orderfields(roiListToAppend, allFields);

            roiListOut = [currentROIList, roiListToAppend];

        end

        function fileName = makeDefaultROIFileName(app)
            %MAKEDEFAULTROIFILENAME Build default .roi filename from current data source.

            if isempty(app.CurrentFile)
                baseName = 'ROIs';
            else
                [~, baseName] = fileparts(app.CurrentFile);
            end

            if ~isempty(app.CurrentEntry)
                entryText = matlab.lang.makeValidName(app.CurrentEntry);
                baseName = sprintf('%s_%s', baseName, entryText);
            end

            fileName = sprintf('%s_ROIs_%s.roi', baseName, datestr(now, 'yyyymmdd_HHMMSS'));

        end

        function ROIFile = buildROIFileStruct(app)
            %BUILDROIFILESTRUCT Build serializable ROI file structure.
            %
            %   ROIFile is saved as a MAT-file with .roi extension. Runtime graphics
            %   are stripped. ROI visibility/selection are intentionally not preserved:
            %   loaded ROIs are always visible and unselected.

            displaySize = app.getDataSize();
            sourceSize = double(app.getSourceDataSize());
            sourceSize = sourceSize(:).';
            sourceSizeDimNames = {'Y', 'X', 'T', 'E'};

            % In normal mode, omit a singleton event dimension from the user-facing
            % sourceSize. Store matching dimension labels to avoid ambiguity.
            if strcmpi(app.ViewMode, 'normal') && numel(sourceSize) >= 4 && sourceSize(4) == 1
                sourceSize = sourceSize(1:3);
                sourceSizeDimNames = sourceSizeDimNames(1:3);
            end

            ROIFile = struct();
            ROIFile.schemaVersion = 1;
            ROIFile.createdOn = datetime('now');
            ROIFile.createdBy = 'DataViewer';
            ROIFile.fileType = 'UMIT_ROI_FILE';

            ROIFile.sourceFile = app.CurrentFile;
            ROIFile.sourceEntry = app.CurrentEntry;
            ROIFile.sourceType = app.getSourceType();

            % Size of the currently displayed image coordinate system.
            ROIFile.imageSizeYX = double(displaySize(1:2));

            % Backend source size with explicit dimension labels.
            ROIFile.sourceSize = sourceSize;
            ROIFile.sourceSizeDimNames = sourceSizeDimNames;

            ROIFile.viewModeAtSave = app.ViewMode;
            ROIFile.currentFrameAtSave = app.CurrentFrame;

            % Event context at save time.
            ROIFile.eventModeInfo = struct();
            ROIFile.eventModeInfo.ConditionName = '';
            ROIFile.eventModeInfo.ConditionID = [];
            ROIFile.eventModeInfo.RepetitionID = [];
            ROIFile.eventModeInfo.RepetitionLabel = '';

            if strcmpi(app.ViewMode, 'event')
                ROIFile.eventModeInfo.ConditionName = app.CurrentCondition;
                ROIFile.eventModeInfo.ConditionID = app.getConditionID(app.CurrentCondition);
                ROIFile.eventModeInfo.RepetitionLabel = app.CurrentRepetition;

                repID = str2double(app.CurrentRepetition);
                if isfinite(repID)
                    ROIFile.eventModeInfo.RepetitionID = repID;
                else
                    ROIFile.eventModeInfo.RepetitionID = [];
                end
            end

            ROIFile.DataParamsView = struct();
            if isfield(app.DataParams, 'view')
                ROIFile.DataParamsView = app.DataParams.view;
            end

            ROIFile.HasUserLogicalMask = app.hasUserLogicalMask();

            % Snapshot of the displayed image used for the current ROI spatial stats/table
            % values. Metadata describing this snapshot is stored at the ROIFile top level.
            ROIFile.statsImage = single(app.getCurrentFrame());
            ROIFile.statsImageViewMode = app.ViewMode;
            ROIFile.statsImageFrameIdx = app.CurrentFrame;
            ROIFile.statsImageFrameTime = app.getCurrentFrameTime(app.CurrentFrame);

            ROIFile.ROIs = app.stripROIRuntimeFieldsForSave(app.ROIList);

        end

        function roiListOut = stripROIRuntimeFieldsForSave(app, roiListIn) %#ok<INUSL>
            %STRIPROIRUNTIMEFIELDSFORSAVE Remove unsavable runtime handles/listeners.
            %
            %   Runtime graphics are not serializable and are rebuilt after loading.
            %   ROI visibility/selection are not saved. On load, ROIs are forced to:
            %       visible  = true
            %       selected = false

            roiListOut = roiListIn;

            if isempty(roiListOut)
                return
            end

            % Remove runtime from the full struct array in one operation. Do not remove
            % fields element-by-element because MATLAB struct arrays require all elements
            % to have the same field set.
            if isfield(roiListOut, 'runtime')
                roiListOut = rmfield(roiListOut, 'runtime');
            end

            % Remove older saved/internal persistent state if present.
            if isfield(roiListOut, 'persistent')
                roiListOut = rmfield(roiListOut, 'persistent');
            end

        end

        function ROIFile = readROIFileStruct(app, filePath) %#ok<INUSL>
            %READROIFILESTRUCT Read ROIFile struct from .roi MAT-file.

            if ~isfile(filePath)
                error('DataViewer:ROIFileNotFound', ...
                    'ROI file not found: "%s".', filePath);
            end

            S = load(filePath, '-mat');

            if isfield(S, 'ROIFile')
                ROIFile = S.ROIFile;
                return
            end

            error('DataViewer:InvalidROIFile', ...
                'The selected file does not contain variable "ROIFile".');

        end

        function validateROIFileStruct(app, ROIFile) %#ok<INUSL>
            %VALIDATEROIFILESTRUCT Validate minimal .roi file structure.
            %
            %   This validator intentionally remains tolerant of older test ROI files.
            %   The load path only requires ROIFile.ROIs and geometry.verticesXY_px.

            if ~isstruct(ROIFile) || ~isscalar(ROIFile)
                error('DataViewer:InvalidROIFile', ...
                    'ROIFile must be a scalar struct.');
            end

            if ~isfield(ROIFile, 'ROIs') || isempty(ROIFile.ROIs) || ~isstruct(ROIFile.ROIs)
                error('DataViewer:InvalidROIFile', ...
                    'ROIFile.ROIs must be a non-empty struct array.');
            end

            for iROI = 1:numel(ROIFile.ROIs)
                ROI = ROIFile.ROIs(iROI);

                if ~isfield(ROI, 'geometry') || ~isstruct(ROI.geometry)
                    error('DataViewer:InvalidROIFile', ...
                        'ROI %d is missing geometry.', iROI);
                end

                if ~isfield(ROI.geometry, 'verticesXY_px') || isempty(ROI.geometry.verticesXY_px)
                    error('DataViewer:InvalidROIFile', ...
                        'ROI %d is missing geometry.verticesXY_px.', iROI);
                end

                verticesXY = double(ROI.geometry.verticesXY_px);

                if size(verticesXY, 2) ~= 2 || size(verticesXY, 1) < 3 || ...
                        any(~isfinite(verticesXY(:)))
                    error('DataViewer:InvalidROIFile', ...
                        'ROI %d has invalid geometry.verticesXY_px.', iROI);
                end
            end

        end

        function [loadAction, scaleFactor] = promptROILoadSizeMismatchAction(app, ROIFile)
            %PROMPTROILOADSIZEMISMATCHACTION Ask how to handle size mismatch on load.
            %
            %   loadAction:
            %       'load'   - load without scaling
            %       'scale'  - scale geometry uniformly to the current image size
            %       'cancel' - cancel loading
            %
            %   Symmetric scaling is available only when X and Y scale factors are equal
            %   within tolerance.

            loadAction = 'load';
            scaleFactor = 1;

            if ~isfield(ROIFile, 'imageSizeYX') || isempty(ROIFile.imageSizeYX)
                return
            end

            savedSizeYX = double(ROIFile.imageSizeYX(:).');
            if numel(savedSizeYX) ~= 2 || any(~isfinite(savedSizeYX)) || any(savedSizeYX <= 0)
                return
            end

            sz = app.getDataSize();
            currentSizeYX = double(sz(1:2));

            if isequal(round(savedSizeYX), round(currentSizeYX))
                return
            end

            savedY = savedSizeYX(1);
            savedX = savedSizeYX(2);
            currentY = currentSizeYX(1);
            currentX = currentSizeYX(2);

            scaleX = currentX ./ savedX;
            scaleY = currentY ./ savedY;

            symmetricTol = max(1e-6, 1e-3 * max(abs([scaleX, scaleY])));
            bCanSymmetricScale = abs(scaleX - scaleY) <= symmetricTol;

            if bCanSymmetricScale
                scaleFactor = mean([scaleX, scaleY]);

                promptText = sprintf([ ...
                    'The ROI file was created for a different image size.\n\n' ...
                    'Saved size:   Y = %.0f, X = %.0f\n' ...
                    'Current size: Y = %.0f, X = %.0f\n\n' ...
                    'Symmetric scale factor: %.6g\n\n' ...
                    'Choose how to load the ROIs.'], ...
                    savedY, savedX, currentY, currentX, scaleFactor);

                try
                    choice = uiconfirm(app.UIFigure, promptText, ...
                        'ROI image-size mismatch', ...
                        'Options', {'Load without scaling', 'Scale ROIs', 'Cancel'}, ...
                        'DefaultOption', 'Scale ROIs', ...
                        'CancelOption', 'Cancel', ...
                        'Icon', 'warning');
                catch
                    choice = questdlg(promptText, ...
                        'ROI image-size mismatch', ...
                        'Load without scaling', 'Scale ROIs', 'Cancel', 'Scale ROIs');

                    if isempty(choice)
                        choice = 'Cancel';
                    end
                end

                switch choice
                    case 'Scale ROIs'
                        loadAction = 'scale';

                    case 'Load without scaling'
                        loadAction = 'load';
                        scaleFactor = 1;

                    otherwise
                        loadAction = 'cancel';
                        scaleFactor = 1;
                end

                return
            end

            promptText = sprintf([ ...
                'The ROI file was created for a different image size, and symmetric scaling is not possible.\n\n' ...
                'Saved size:   Y = %.0f, X = %.0f\n' ...
                'Current size: Y = %.0f, X = %.0f\n\n' ...
                'Scale X = %.6g, Scale Y = %.6g'], ...
                savedY, savedX, currentY, currentX, scaleX, scaleY);

            try
                choice = uiconfirm(app.UIFigure, promptText, ...
                    'ROI image-size mismatch', ...
                    'Options', {'Load without scaling', 'Cancel'}, ...
                    'DefaultOption', 'Load without scaling', ...
                    'CancelOption', 'Cancel', ...
                    'Icon', 'warning');
            catch
                choice = questdlg(promptText, ...
                    'ROI image-size mismatch', ...
                    'Load without scaling', 'Cancel', 'Load without scaling');

                if isempty(choice)
                    choice = 'Cancel';
                end
            end

            if strcmp(choice, 'Cancel')
                loadAction = 'cancel';
            else
                loadAction = 'load';
            end

            scaleFactor = 1;

        end

        function clearAllROIsWithoutPrompt(app)
            %CLEARALLROISWITHOUTPROMPT Delete ROI graphics and clear ROIList.

            try
                app.stopActiveROIEditForSelectionChange();
            catch
            end

            if ~isempty(app.ROIList)
                for iROI = 1:numel(app.ROIList)
                    try
                        app.deleteROIGraphicsByIndex(iROI);
                    catch
                    end
                end
            end

            app.ROIList = struct([]);
            app.SelectedROIID = NaN;
            app.ROISelectionOrder = [];

        end

        function templateROI = makeLoadedROITemplate(app)
            %MAKELOADEDROITEMPLATE Return canonical ROI struct template for load rebuilds.

            templateROI = struct();
            templateROI.name = '';
            templateROI.type = 'polygon';
            templateROI.DOC = datetime('now');
            templateROI.modifiedOn = templateROI.DOC;
            templateROI.color = [1 0 0];
            templateROI.notes = '';
            templateROI.geometry = struct();
            templateROI.mask = [];
            templateROI.stats = struct();
            templateROI.ID = NaN;
            templateROI.runtime = struct();

        end

        function [roiOutList, status] = rebuildOneLoadedROIForCurrentImage(app, ROI, usedIDs, usedNames)
            %REBUILDONELOADEDROIFORCURRENTIMAGE Rebuild one loaded source ROI.

            templateROI = app.makeLoadedROITemplate();
            roiOutList = repmat(templateROI, 1, 0);

            status = struct();
            status.bClipped = false;
            status.bSplit = false;

            [loaded, ok] = app.normalizeLoadedROIRecord(ROI, usedIDs, usedNames);
            if ~ok
                return
            end

            verticesXY = loaded.verticesXY;

            rawMask = app.createMaskFromVertices(verticesXY);
            if isempty(rawMask) || ~any(rawMask(:))
                return
            end

            clippedMask = app.clipROIMaskToActiveLogicalMask(rawMask);
            if isempty(clippedMask) || ~any(clippedMask(:))
                return
            end

            status.bClipped = ~isequal(rawMask, clippedMask);

            if status.bClipped
                componentMasks = app.maskToConnectedComponentMasks(clippedMask);

                if isempty(componentMasks)
                    return
                end

                status.bSplit = numel(componentMasks) > 1;

                template = templateROI;
                template.name = loaded.name;
                template.type = 'polygon';
                template.DOC = loaded.DOC;
                template.modifiedOn = datetime('now');
                template.color = loaded.color;
                template.notes = loaded.notes;
                template.ID = loaded.ID;
                template.runtime = app.makeROIRuntimeStruct(true);

                pendingNames = cellstr(string(usedNames(:)));
                pendingIDs = double(usedIDs(:).');
                pendingIDs = pendingIDs(isfinite(pendingIDs));

                for iComp = 1:numel(componentMasks)
                    if numel(componentMasks) > 1
                        thisName = sprintf('%s_part%d', loaded.name, iComp);
                    else
                        thisName = loaded.name;
                    end

                    thisName = app.makeUniqueNameAgainstList(thisName, pendingNames);

                    if iComp == 1
                        thisID = loaded.ID;
                    else
                        thisID = app.getNextUniqueROIIDFromList(pendingIDs);
                    end

                    try
                        newROI = app.makePolygonROIFromMaskComponent( ...
                            template, componentMasks{iComp}, thisName, thisID, false);
                        newROI.runtime.visible = true;
                        newROI.runtime.selected = false;
                    catch
                        continue
                    end

                    roiOutList(end+1) = newROI; %#ok<AGROW>
                    pendingNames{end+1} = thisName; %#ok<AGROW>
                    pendingIDs(end+1) = thisID; %#ok<AGROW>
                end

                return
            end

            % No logical-mask clipping: preserve original ROI type/geometry.
            try
                pgon = polyshape(verticesXY(:, 1), verticesXY(:, 2), 'Simplify', true);
            catch
                pgon = polyshape();
            end

            if isempty(pgon.Vertices)
                return
            end

            roiType = loaded.type;
            roiParams = loaded.params;

            if strcmpi(roiType, 'polygon') || isempty(fieldnames(roiParams))
                roiType = 'polygon';
                roiParams = app.makePolygonROIParametersFromVertices(verticesXY);
            end

            newROI = templateROI;
            newROI.name = loaded.name;
            newROI.type = roiType;
            newROI.DOC = loaded.DOC;
            newROI.modifiedOn = datetime('now');
            newROI.color = loaded.color;
            newROI.notes = loaded.notes;
            newROI.ID = loaded.ID;

            newROI.geometry = struct();
            newROI.geometry.polyshape = pgon;
            newROI.geometry.verticesXY_px = verticesXY;
            newROI.geometry.ROIType = roiType;
            newROI.geometry.ROIParameters = roiParams;

            newROI.mask = logical(rawMask);
            newROI.stats = app.computeROIStatsFromMask(newROI.mask);
            newROI.runtime = app.makeROIRuntimeStruct(true);

            roiOutList = newROI;

        end

        function [loaded, ok] = normalizeLoadedROIRecord(app, ROI, usedIDs, usedNames)
            %NORMALIZELOADEDROIRECORD Normalize one saved ROI for current-session loading.
            %
            %   This centralizes all compatibility fallbacks for older .roi files.

            ok = false;

            loaded = struct();
            loaded.verticesXY = zeros(0, 2);
            loaded.type = 'polygon';
            loaded.params = struct();
            loaded.name = 'ROI';
            loaded.ID = app.getNextUniqueROIIDFromList(usedIDs);
            loaded.color = [1 0 0];
            loaded.DOC = datetime('now');
            loaded.notes = '';

            try
                if isfield(ROI, 'geometry') && isfield(ROI.geometry, 'verticesXY_px')
                    loaded.verticesXY = app.cleanROIVertices(double(ROI.geometry.verticesXY_px));
                end
            catch
                loaded.verticesXY = zeros(0, 2);
            end

            if size(loaded.verticesXY, 1) < 3
                return
            end

            try
                if isfield(ROI, 'type') && ~isempty(ROI.type)
                    loaded.type = lower(char(string(ROI.type)));
                elseif isfield(ROI, 'geometry') && isfield(ROI.geometry, 'ROIType') && ...
                        ~isempty(ROI.geometry.ROIType)
                    loaded.type = lower(char(string(ROI.geometry.ROIType)));
                end
            catch
                loaded.type = 'polygon';
            end

            if ~ismember(loaded.type, {'rectangle', 'ellipse', 'polygon'})
                loaded.type = 'polygon';
            end

            try
                if isfield(ROI, 'geometry') && isfield(ROI.geometry, 'ROIParameters') && ...
                        isstruct(ROI.geometry.ROIParameters)
                    loaded.params = ROI.geometry.ROIParameters;
                end
            catch
                loaded.params = struct();
            end

            try
                if isfield(ROI, 'name') && ~isempty(ROI.name)
                    loaded.name = strtrim(char(string(ROI.name)));
                end
            catch
                loaded.name = 'ROI';
            end

            if isempty(loaded.name)
                loaded.name = 'ROI';
            end
            loaded.name = app.makeUniqueNameAgainstList(loaded.name, usedNames);

            try
                if isfield(ROI, 'ID') && isfinite(ROI.ID) && ~any(double(usedIDs) == double(ROI.ID))
                    loaded.ID = double(ROI.ID);
                else
                    loaded.ID = app.getNextUniqueROIIDFromList(usedIDs);
                end
            catch
                loaded.ID = app.getNextUniqueROIIDFromList(usedIDs);
            end

            try
                if isfield(ROI, 'color') && numel(ROI.color) == 3
                    loaded.color = double(ROI.color(:).');
                elseif ~isempty(app.ROIColorList)
                    loaded.color = app.ROIColorList(app.ROINextColorIndex, :);
                end
            catch
                loaded.color = [1 0 0];
            end
            loaded.color = min(max(loaded.color, 0), 1);

            try
                if isfield(ROI, 'DOC') && ~isempty(ROI.DOC)
                    loaded.DOC = ROI.DOC;
                end
            catch
                loaded.DOC = datetime('now');
            end

            try
                if isfield(ROI, 'notes') && ~isempty(ROI.notes)
                    loaded.notes = char(string(ROI.notes));
                end
            catch
                loaded.notes = '';
            end

            ok = true;

        end

        function ROI = scaleLoadedROIForCurrentImage(app, ROI, scaleFactor) %#ok<INUSL>
            %SCALELOADEDROIFORCURRENTIMAGE Uniformly scale saved ROI geometry.
            %
            %   Pixel-boundary preserving transform:
            %       xNew = (xOld - 0.5) * scale + 0.5

            if isempty(scaleFactor) || ~isfinite(scaleFactor) || scaleFactor <= 0 || scaleFactor == 1
                return
            end

            if ~isfield(ROI, 'geometry') || ~isstruct(ROI.geometry)
                return
            end

            if isfield(ROI.geometry, 'verticesXY_px') && ~isempty(ROI.geometry.verticesXY_px)
                ROI.geometry.verticesXY_px = iScaleXY(ROI.geometry.verticesXY_px);
            end

            if isfield(ROI.geometry, 'polyshape')
                try
                    verticesXY = double(ROI.geometry.verticesXY_px);
                    ROI.geometry.polyshape = polyshape(verticesXY(:, 1), verticesXY(:, 2), 'Simplify', true);
                catch
                    ROI.geometry.polyshape = polyshape();
                end
            end

            if isfield(ROI.geometry, 'ROIParameters') && isstruct(ROI.geometry.ROIParameters)
                ROI.geometry.ROIParameters = iScaleParams(ROI.geometry.ROIParameters);
            end

            ROI.mask = [];
            ROI.stats = struct();

            function xyOut = iScaleXY(xyIn)
                xyOut = (double(xyIn) - 0.5) .* scaleFactor + 0.5;
            end

            function params = iScaleParams(params)
                if isempty(params) || ~isstruct(params) || ~isscalar(params)
                    return
                end

                if isfield(params, 'Position') && ~isempty(params.Position) && isnumeric(params.Position)
                    pos = double(params.Position);

                    if ismatrix(pos) && size(pos, 2) == 2
                        params.Position = iScaleXY(pos);
                    elseif isvector(pos) && numel(pos) == 4
                        pos = pos(:).';
                        xy = iScaleXY(pos(1:2));
                        params.Position = [xy, pos(3:4) .* scaleFactor];
                    end
                end

                if isfield(params, 'Vertices') && ~isempty(params.Vertices) && isnumeric(params.Vertices) && ...
                        size(params.Vertices, 2) == 2
                    params.Vertices = iScaleXY(params.Vertices);
                end

                if isfield(params, 'Center') && ~isempty(params.Center) && isnumeric(params.Center)
                    params.Center = iScaleXY(params.Center(:).');
                end

                if isfield(params, 'SemiAxes') && ~isempty(params.SemiAxes) && isnumeric(params.SemiAxes)
                    params.SemiAxes = double(params.SemiAxes) .* scaleFactor;
                end
            end

        end

        function uniqueName = makeUniqueNameAgainstList(app, requestedName, existingNames) %#ok<INUSL>
            %MAKEUNIQUENAMEAGAINSTLIST Return a unique name against supplied names.

            requestedName = strtrim(char(string(requestedName)));

            if isempty(requestedName)
                requestedName = 'ROI';
            end

            existingNames = string(existingNames(:));
            uniqueName = requestedName;

            if ~any(strcmp(existingNames, uniqueName))
                return
            end

            suffix = 2;
            while any(strcmp(existingNames, sprintf('%s_%d', requestedName, suffix)))
                suffix = suffix + 1;
            end

            uniqueName = sprintf('%s_%d', requestedName, suffix);

        end

        function roiID = getNextUniqueROIIDFromList(app, usedIDs) %#ok<INUSL>
            %GETNEXTUNIQUEROIIDFROMLIST Return next finite ROI ID outside usedIDs.

            usedIDs = double(usedIDs(:).');
            usedIDs = usedIDs(isfinite(usedIDs));

            roiID = max([0, usedIDs]) + 1;

            while any(usedIDs == roiID)
                roiID = roiID + 1;
            end

        end

        function rebuildMissingROIOverlaysAfterLoad(app)
            %REBUILDMISSINGROIOVERLAYSAFTERLOAD Recreate missing ROI display handles.
            %
            %   Load/import helpers build the ROI data model first. Runtime overlays are
            %   created only after ROIList has been assigned to avoid transient graphics.

            if isempty(app.ROIList)
                return
            end

            for iROI = 1:numel(app.ROIList)
                if ~isfield(app.ROIList(iROI), 'runtime') || ~isstruct(app.ROIList(iROI).runtime)
                    app.ROIList(iROI).runtime = app.makeROIRuntimeStruct(true);
                end

                app.ROIList(iROI).runtime.visible = true;
                app.ROIList(iROI).runtime.selected = false;

                if ~isfield(app.ROIList(iROI).runtime, 'ROIHandle') || ...
                        ~app.isUsableGraphicsHandle(app.ROIList(iROI).runtime.ROIHandle)
                    app.ROIList(iROI).runtime.ROIHandle = app.createStaticROIOverlayFromROI(app.ROIList(iROI));
                end
            end

        end

        function runtime = makeROIRuntimeStruct(app, visibleValue) %#ok<INUSL>
            %MAKEROIRUNTIMESTRUCT Return initialized ROI runtime struct.

            if nargin < 2 || isempty(visibleValue)
                visibleValue = true;
            end

            runtime = struct();
            runtime.visible = logical(visibleValue);
            runtime.selected = false;
            runtime.ROIHandle = gobjects(1);
            runtime.editHandle = gobjects(1);
            runtime.traceLine = gobjects(1);
            runtime.tracePatch = gobjects(1);

            runtime.trace = struct();
            runtime.trace.XData = [];
            runtime.trace.Mean = [];
            runtime.trace.SEM = [];
            runtime.trace.Mode = '';
            runtime.trace.Status = 'not computed';

        end

        function importROIsFromExternalSource(app)
            %IMPORTROISFROMEXTERNALSOURCE Import ROIs from a MAT-file or workspace variable.
            %
            %   Supported input formats follow the legacy ROImanager import rules:
            %       - 2D logical mask: each connected component becomes one ROI.
            %       - 3D logical stack: each Y-by-X frame is imported as one source ROI.
            %       - 2D label matrix: each positive integer label is imported as one source ROI.
            %
            %   Imported arrays may either match the current image Y/X size exactly or be
            %   symmetrically scalable to the current image size. For size mismatches, the
            %   user is prompted to either scale the imported masks or cancel. There is no
            %   "load without scaling" option for import because raw masks do not carry a
            %   durable geometry coordinate system like .roi files do.

            if ~app.hasData()
                app.setStatusMessage('Load image data before importing ROIs.');
                return
            end

            try
                [roiData, sourceInfo, ok] = app.readImportedROIVariable();

                if ~ok
                    app.setStatusMessage('Import ROI cancelled.');
                    return
                end

                [roiData, sourceInfo, ok] = app.prepareImportedROIDataForCurrentImage(roiData, sourceInfo);

                if ~ok
                    app.setStatusMessage('Import ROI cancelled.');
                    return
                end

                [newROIList, report] = app.buildImportedROIListFromArray(roiData, sourceInfo);

                if isempty(newROIList)
                    statusText = app.formatROIImportStatusText(0, report, sourceInfo);
                    app.setStatusMessage(statusText);
                    return
                end

                app.appendImportedROIsAndRefresh(newROIList);

                statusText = app.formatROIImportStatusText(numel(newROIList), report, sourceInfo);
                app.setStatusMessage(statusText);

            catch ME
                app.setStatusMessage(sprintf('Import ROI failed: %s', ME.message));
                rethrow(ME)
            end

        end

        function [roiData, sourceInfo, ok] = readImportedROIVariable(app)
            %READIMPORTEDROIVARIABLE Read one supported ROI variable from MAT-file/workspace.

            roiData = [];
            sourceInfo = struct();
            ok = false;

            importSource = app.promptROIImportSource();

            switch importSource
                case 'mat'
                    startFolder = pwd;
                    if ~isempty(app.CurrentFile)
                        startFolder = fileparts(app.CurrentFile);
                    end

                    [fileName, folderName] = uigetfile( ...
                        {'*.mat', 'MAT-file (*.mat)'}, ...
                        'Import ROIs from MAT-file', ...
                        startFolder);

                    if isequal(fileName, 0)
                        return
                    end

                    filePath = fullfile(folderName, fileName);

                    infoList = whos('-file', filePath);
                    candidateTable = app.buildROIImportCandidateTable(infoList);

                    if isempty(candidateTable) || height(candidateTable) == 0
                        app.setStatusMessage('No supported ROI variables were found in the selected MAT-file.');
                        return
                    end

                    varName = app.selectImportedROIVariableDialog(candidateTable, ...
                        sprintf('Select ROI variable from "%s"', fileName));

                    if isempty(varName)
                        return
                    end

                    S = load(filePath, varName);
                    roiData = S.(varName);

                    sourceInfo = struct();
                    sourceInfo.sourceType = 'mat';
                    sourceInfo.sourcePath = filePath;
                    sourceInfo.sourceLabel = filePath;
                    sourceInfo.variableName = varName;

                case 'workspace'
                    infoList = evalin('base', 'whos');
                    candidateTable = app.buildROIImportCandidateTable(infoList);

                    if isempty(candidateTable) || height(candidateTable) == 0
                        app.setStatusMessage('No supported ROI variables were found in the base workspace.');
                        return
                    end

                    varName = app.selectImportedROIVariableDialog(candidateTable, ...
                        'Select ROI variable from base workspace');

                    if isempty(varName)
                        return
                    end

                    roiData = evalin('base', varName);

                    sourceInfo = struct();
                    sourceInfo.sourceType = 'workspace';
                    sourceInfo.sourcePath = '';
                    sourceInfo.sourceLabel = 'base workspace';
                    sourceInfo.variableName = varName;

                otherwise
                    return
            end

            ok = true;

        end

        function importSource = promptROIImportSource(app)
            %PROMPTROIIMPORTSOURCE Ask whether to import from MAT-file or workspace.

            importSource = 'cancel';

            promptText = [ ...
                'Import ROIs from a non-.roi source.' newline newline ...
                'Supported formats:' newline ...
                '  - 2D logical mask' newline ...
                '  - 3D logical mask stack' newline ...
                '  - 2D positive-integer label matrix'];

            try
                choice = uiconfirm(app.UIFigure, promptText, ...
                    'Import ROIs', ...
                    'Options', {'MAT-file', 'Base workspace', 'Cancel'}, ...
                    'DefaultOption', 'MAT-file', ...
                    'CancelOption', 'Cancel', ...
                    'Icon', 'question');
            catch
                choice = questdlg(promptText, ...
                    'Import ROIs', ...
                    'MAT-file', 'Base workspace', 'Cancel', 'MAT-file');

                if isempty(choice)
                    choice = 'Cancel';
                end
            end

            switch choice
                case 'MAT-file'
                    importSource = 'mat';

                case 'Base workspace'
                    importSource = 'workspace';

                otherwise
                    importSource = 'cancel';
            end

        end

        function candidateTable = buildROIImportCandidateTable(app, infoList)
            %BUILDROIIMPORTCANDIDATETABLE Return supported ROI import variables.
            %
            %   Candidate filtering uses WHOS/WHOS('-file') metadata. Variables are listed
            %   when they are either:
            %       - exact Y-by-X matches for the current image, or
            %       - symmetrically scalable to the current image size.
            %
            %   EstimatedROIs is a string column. It reports the known source ROI count
            %   when this can be inferred from metadata. Otherwise it shows "unknown".

            candidateTable = table();

            if isempty(infoList)
                return
            end

            names = strings(0, 1);
            classes = strings(0, 1);
            sizeText = strings(0, 1);
            importTypes = strings(0, 1);
            sizeCompatibility = strings(0, 1);
            scaleFactors = nan(0, 1);
            estimatedROIs = strings(0, 1);

            for iVar = 1:numel(infoList)
                varSize = double(infoList(iVar).size);
                varClass = char(string(infoList(iVar).class));

                [isSupported, importType, estimate, compatibilityText, scaleFactor] = ...
                    iClassifyCandidate(varClass, varSize);

                if ~isSupported
                    continue
                end

                names(end+1, 1) = string(infoList(iVar).name); %#ok<AGROW>
                classes(end+1, 1) = string(varClass); %#ok<AGROW>
                sizeText(end+1, 1) = string(iFormatSize(varSize)); %#ok<AGROW>
                importTypes(end+1, 1) = string(importType); %#ok<AGROW>
                sizeCompatibility(end+1, 1) = string(compatibilityText); %#ok<AGROW>
                scaleFactors(end+1, 1) = scaleFactor; %#ok<AGROW>
                estimatedROIs(end+1, 1) = string(estimate); %#ok<AGROW>
            end

            if isempty(names)
                return
            end

            candidateTable = table( ...
                names, ...
                classes, ...
                sizeText, ...
                importTypes, ...
                sizeCompatibility, ...
                scaleFactors, ...
                estimatedROIs, ...
                'VariableNames', {'Name', 'Class', 'Size', 'ImportType', ...
                'SizeCompatibility', 'ScaleFactor', 'EstimatedROIs'});

            function [tf, importType, estimate, compatibilityText, scaleFactor] = iClassifyCandidate(varClass, varSize)
                tf = false;
                importType = '';
                estimate = "unknown";
                compatibilityText = '';
                scaleFactor = NaN;

                [sizeMode, scaleFactor, sizeMsg] = app.getROIImportSizeCompatibilityFromSize(varSize);

                if strcmpi(sizeMode, 'incompatible')
                    return
                end

                if strcmpi(sizeMode, 'exact')
                    compatibilityText = 'exact size';
                else
                    compatibilityText = sprintf('scale by %.6g', scaleFactor);
                end

                % Ignore trailing singleton dimensions in the candidate classification.
                classSize = varSize;
                while numel(classSize) > 2 && classSize(end) == 1
                    classSize(end) = [];
                end

                isLogical = strcmpi(varClass, 'logical');
                isNumeric = ismember(lower(varClass), { ...
                    'double', 'single', ...
                    'uint8', 'uint16', 'uint32', 'uint64', ...
                    'int8', 'int16', 'int32', 'int64'});

                if isLogical && numel(classSize) == 2
                    tf = true;
                    importType = '2D logical mask';
                    estimate = "unknown";
                    return
                end

                if isLogical && numel(classSize) == 3
                    tf = true;
                    importType = '3D logical mask stack';
                    estimate = string(classSize(3));
                    return
                end

                if isNumeric && numel(classSize) == 2
                    tf = true;
                    importType = '2D label matrix';
                    estimate = "unknown";
                    return
                end

                compatibilityText = sizeMsg;
            end

            function txt = iFormatSize(varSize)
                varSize = double(varSize(:).');
                parts = arrayfun(@(x) sprintf('%.0f', x), varSize, 'UniformOutput', false);
                txt = strjoin(parts, ' x ');
            end

        end

        function selectedName = selectImportedROIVariableDialog(app, candidateTable, dialogTitle)
            %SELECTIMPORTEDROIVARIABLEDIALOG Modal table for choosing an import variable.

            selectedName = '';

            if isempty(candidateTable) || height(candidateTable) == 0
                return
            end

            selectedRow = 1;

            dlg = uifigure( ...
                'Name', dialogTitle, ...
                'WindowStyle', 'modal', ...
                'Position', [100 100 640 360], ...
                'Visible', 'off', ...
                'CloseRequestFcn', @onCancel);

            grid = uigridlayout(dlg);
            grid.RowHeight = {28, '1x', 28, 36};
            grid.ColumnWidth = {'1x'};
            grid.Padding = [12 12 12 12];

            titleLabel = uilabel(grid);
            titleLabel.Text = 'Select the variable to import as ROIs:';
            titleLabel.FontWeight = 'bold';
            titleLabel.Layout.Row = 1;
            titleLabel.Layout.Column = 1;

            varTable = uitable(grid);
            varTable.Data = candidateTable;
            varTable.ColumnEditable = false(1, width(candidateTable));
            varTable.RowName = {};
            varTable.Layout.Row = 2;
            varTable.Layout.Column = 1;
            varTable.SelectionChangedFcn = @onSelectionChanged;

            try
                if isprop(varTable, 'SelectionType')
                    varTable.SelectionType = 'row';
                end
            catch
            end

            try
                if isprop(varTable, 'Multiselect')
                    varTable.Multiselect = 'off';
                end
            catch
            end

            try
                if isprop(varTable, 'SelectionType') && strcmpi(varTable.SelectionType, 'row')
                    varTable.Selection = 1;
                else
                    varTable.Selection = [1 1];
                end
            catch
            end

            statusLabel = uilabel(grid);
            statusLabel.Text = '';
            statusLabel.FontColor = [0.65 0 0];
            statusLabel.Layout.Row = 3;
            statusLabel.Layout.Column = 1;

            buttonGrid = uigridlayout(grid);
            buttonGrid.RowHeight = {'1x'};
            buttonGrid.ColumnWidth = {'1x', '1x'};
            buttonGrid.Padding = [0 0 0 0];
            buttonGrid.Layout.Row = 4;
            buttonGrid.Layout.Column = 1;

            okButton = uibutton(buttonGrid, 'push');
            okButton.Text = 'OK';
            okButton.Layout.Row = 1;
            okButton.Layout.Column = 1;
            okButton.ButtonPushedFcn = @onOK;

            cancelButton = uibutton(buttonGrid, 'push');
            cancelButton.Text = 'Cancel';
            cancelButton.Layout.Row = 1;
            cancelButton.Layout.Column = 2;
            cancelButton.ButtonPushedFcn = @onCancel;

            try
                placeAppInsideCaller(app, dlg, 'center');
            catch
            end
            dlg.Visible = 'on';

            uiwait(dlg);

            if isvalid(dlg)
                delete(dlg);
            end

            function onSelectionChanged(~, event)
                selectedRow = iGetSelectedRow(event);

                if isempty(selectedRow) || selectedRow < 1 || selectedRow > height(candidateTable)
                    selectedRow = [];
                    return
                end

                try
                    if isprop(varTable, 'SelectionType') && strcmpi(varTable.SelectionType, 'row')
                        varTable.Selection = selectedRow;
                    else
                        varTable.Selection = [selectedRow, 1];
                    end
                catch
                end

                statusLabel.Text = '';
            end

            function rowIdx = iGetSelectedRow(event)
                rowIdx = [];

                try
                    if isprop(event, 'Selection') && ~isempty(event.Selection)
                        selection = event.Selection;
                    elseif isprop(varTable, 'Selection') && ~isempty(varTable.Selection)
                        selection = varTable.Selection;
                    else
                        return
                    end

                    if isvector(selection)
                        rowIdx = selection(1);
                    else
                        rowIdx = selection(1, 1);
                    end

                    rowIdx = round(double(rowIdx));
                catch
                    rowIdx = [];
                end
            end

            function onOK(~, ~)
                if isempty(selectedRow) || selectedRow < 1 || selectedRow > height(candidateTable)
                    statusLabel.Text = 'Select one variable.';
                    return
                end

                selectedName = char(candidateTable.Name(selectedRow));
                uiresume(dlg);
            end

            function onCancel(~, ~)
                selectedName = '';
                uiresume(dlg);
            end

        end

        function [newROIList, report] = buildImportedROIListFromArray(app, roiData, sourceInfo)
            %BUILDIMPORTEDROILISTFROMARRAY Convert supported arrays into polygon ROIs.
            %
            %   Imported ROIs are built without graphics overlays. The caller appends the
            %   resulting structs to app.ROIList and then calls rebuildMissingROIOverlaysAfterLoad.

            templateROI = app.makeLoadedROITemplate();
            newROIList = repmat(templateROI, 1, 0);

            report = struct();
            report.nSourceMasks = 0;
            report.nCreated = 0;
            report.nSkippedEmpty = 0;
            report.nSkippedHoles = 0;
            report.nSkippedInvalid = 0;
            report.nClipped = 0;
            report.nSplit = 0;
            report.importType = '';

            if isempty(roiData)
                report.nSkippedInvalid = 1;
                return
            end

            if issparse(roiData)
                roiData = full(roiData);
            end

            sz = app.getDataSize();
            Ny = sz(1);
            Nx = sz(2);

            if size(roiData, 1) ~= Ny || size(roiData, 2) ~= Nx
                error('DataViewer:InvalidImportedROISize', ...
                    'Imported ROI variable must match current image size [Y X] = [%d %d].', Ny, Nx);
            end

            baseName = char(string(sourceInfo.variableName));
            if isempty(strtrim(baseName))
                baseName = 'ImportedROI';
            end

            existingIDs = [];
            existingNames = {};

            if ~isempty(app.ROIList)
                if isfield(app.ROIList, 'ID')
                    existingIDs = [app.ROIList.ID];
                    existingIDs = double(existingIDs(isfinite(existingIDs)));
                end

                if isfield(app.ROIList, 'name')
                    existingNames = cell(1, numel(app.ROIList));
                    for iROI = 1:numel(app.ROIList)
                        existingNames{iROI} = char(string(app.ROIList(iROI).name));
                    end
                end
            end

            pendingIDs = existingIDs(:).';
            pendingNames = cellstr(string(existingNames(:)));

            if islogical(roiData)
                if ndims(roiData) <= 2
                    report.importType = '2D logical mask';
                    rawComponents = app.maskToConnectedComponentMasks(logical(roiData));

                    if isempty(rawComponents)
                        report.nSkippedEmpty = report.nSkippedEmpty + 1;
                        return
                    end

                    for iComp = 1:numel(rawComponents)
                        if numel(rawComponents) == 1
                            sourceName = baseName;
                        else
                            sourceName = sprintf('%s_region%d', baseName, iComp);
                        end

                        iAddSourceMask(rawComponents{iComp}, sourceName);
                    end

                    return
                end

                if ndims(roiData) == 3
                    report.importType = '3D logical mask stack';
                    nFrames = size(roiData, 3);

                    for iFrame = 1:nFrames
                        sourceName = sprintf('%s_frame%d', baseName, iFrame);
                        iAddSourceMask(logical(roiData(:, :, iFrame)), sourceName);
                    end

                    return
                end

                report.nSkippedInvalid = report.nSkippedInvalid + 1;
                return
            end

            if isnumeric(roiData) && ismatrix(roiData)
                report.importType = '2D label matrix';

                labelMatrix = double(roiData);

                if any(~isfinite(labelMatrix(:))) || any(labelMatrix(:) < 0) || ...
                        any(abs(labelMatrix(:) - round(labelMatrix(:))) > 1e-9)
                    error('DataViewer:InvalidImportedLabelMatrix', ...
                        'Imported label matrix must contain finite non-negative integer labels.');
                end

                labels = unique(labelMatrix(:));
                labels = labels(labels > 0);

                if isempty(labels)
                    report.nSkippedEmpty = report.nSkippedEmpty + 1;
                    return
                end

                for iLabel = 1:numel(labels)
                    labelValue = labels(iLabel);
                    sourceName = sprintf('%s_label%d', baseName, labelValue);
                    iAddSourceMask(labelMatrix == labelValue, sourceName);
                end

                return
            end

            error('DataViewer:UnsupportedImportedROIFormat', ...
                'Unsupported ROI import variable. Use 2D/3D logical arrays or a 2D label matrix.');

            function iAddSourceMask(rawMask, sourceName)
                rawMask = logical(rawMask);
                report.nSourceMasks = report.nSourceMasks + 1;

                if isempty(rawMask) || ~any(rawMask(:))
                    report.nSkippedEmpty = report.nSkippedEmpty + 1;
                    return
                end

                clippedMask = app.clipROIMaskToActiveLogicalMask(rawMask);

                if isempty(clippedMask) || ~any(clippedMask(:))
                    report.nSkippedEmpty = report.nSkippedEmpty + 1;
                    return
                end

                if ~isequal(rawMask, clippedMask)
                    report.nClipped = report.nClipped + 1;
                end

                componentMasks = app.maskToConnectedComponentMasks(clippedMask);

                if isempty(componentMasks)
                    report.nSkippedEmpty = report.nSkippedEmpty + 1;
                    return
                end

                if numel(componentMasks) > 1
                    report.nSplit = report.nSplit + 1;
                end

                for iPart = 1:numel(componentMasks)
                    componentMask = logical(componentMasks{iPart});

                    if app.roiMaskHasHoles(componentMask)
                        report.nSkippedHoles = report.nSkippedHoles + 1;
                        continue
                    end

                    if numel(componentMasks) == 1
                        requestedName = sourceName;
                    else
                        requestedName = sprintf('%s_part%d', sourceName, iPart);
                    end

                    uniqueName = app.makeUniqueNameAgainstList(requestedName, pendingNames);
                    roiID = iNextID();
                    roiColor = iNextColor();

                    template = app.makeLoadedROITemplate();
                    template.name = uniqueName;
                    template.type = 'polygon';
                    template.DOC = datetime('now');
                    template.modifiedOn = template.DOC;
                    template.color = roiColor;
                    template.notes = sprintf('Imported from %s variable "%s".', ...
                        sourceInfo.sourceLabel, sourceInfo.variableName);
                    template.ID = roiID;
                    template.runtime = app.makeROIRuntimeStruct(true);

                    try
                        newROI = app.makePolygonROIFromMaskComponent( ...
                            template, componentMask, uniqueName, roiID, false);
                    catch
                        report.nSkippedInvalid = report.nSkippedInvalid + 1;
                        continue
                    end

                    newROIList(end+1) = newROI; %#ok<AGROW>
                    pendingIDs(end+1) = roiID; %#ok<AGROW>
                    pendingNames{end+1} = uniqueName; %#ok<AGROW>
                    report.nCreated = report.nCreated + 1;
                end
            end

            function roiID = iNextID()
                roiID = app.getNextUniqueROIIDFromList(pendingIDs);
            end

            function roiColor = iNextColor()
                if isempty(app.ROIColorList)
                    roiColor = [1 0 0];
                    return
                end

                nColors = size(app.ROIColorList, 1);
                nextOrdinal = numel(app.ROIList) + numel(newROIList) + 1;
                colorIdx = mod(nextOrdinal - 1, nColors) + 1;

                roiColor = app.ROIColorList(colorIdx, :);
                roiColor = min(max(double(roiColor(:).'), 0), 1);
            end

        end

        function appendImportedROIsAndRefresh(app, newROIList)
            %APPENDIMPORTEDROISANDREFRESH Append imported ROIs and refresh runtime views.

            if isempty(newROIList)
                return
            end

            app.stopActiveROIEditForSelectionChange();

            app.ROIList = app.appendROIsToROIList(app.ROIList, newROIList);
            app.SelectedROIID = NaN;
            app.ROISelectionOrder = [];

            app.rebuildMissingROIOverlaysAfterLoad();

            if ~isempty(app.ROIColorList)
                app.ROINextColorIndex = mod(numel(app.ROIList), size(app.ROIColorList, 1)) + 1;
            end

            app.refreshROITable();
            app.refreshROITraces();
            app.refreshEventPatches();
            app.stackROITraceGraphics();
            app.updateGUIEnabledState();

        end

        function statusText = formatROIImportStatusText(app, nAdded, report, sourceInfo) %#ok<INUSL>
            %FORMATROIIMPORTSTATUSTEXT Build compact import result message.

            statusText = sprintf('Imported %d ROI(s) from "%s" (%s).', ...
                nAdded, sourceInfo.variableName, report.importType);

            if isfield(sourceInfo, 'bScaled') && sourceInfo.bScaled
                statusText = sprintf( ...
                    '%s Scaled masks from [Y X] = [%.0f %.0f] to [%.0f %.0f] using scale %.6g.', ...
                    statusText, ...
                    sourceInfo.originalSizeYX(1), sourceInfo.originalSizeYX(2), ...
                    sourceInfo.targetSizeYX(1), sourceInfo.targetSizeYX(2), ...
                    sourceInfo.scaleFactor);
            end

            if report.nClipped > 0
                statusText = sprintf('%s %d source mask(s) clipped to the logical mask.', ...
                    statusText, report.nClipped);
            end

            if report.nSplit > 0
                statusText = sprintf('%s %d source mask(s) split into parts.', ...
                    statusText, report.nSplit);
            end

            if report.nSkippedEmpty > 0
                statusText = sprintf('%s %d empty/outside mask(s) skipped.', ...
                    statusText, report.nSkippedEmpty);
            end

            if report.nSkippedHoles > 0
                statusText = sprintf('%s %d hole-containing mask(s) skipped.', ...
                    statusText, report.nSkippedHoles);
            end

            if report.nSkippedInvalid > 0
                statusText = sprintf('%s %d invalid mask(s) skipped.', ...
                    statusText, report.nSkippedInvalid);
            end

        end

        function [roiData, sourceInfo, ok] = prepareImportedROIDataForCurrentImage(app, roiData, sourceInfo)
            %PREPAREIMPORTEDROIDATAFORCURRENTIMAGE Validate/import-size and scale if needed.
            %
            %   Import sources are raw masks/labels, so size mismatch handling is stricter
            %   than .roi loading:
            %       - exact size: import directly
            %       - symmetric size mismatch: prompt Scale ROIs / Cancel
            %       - non-symmetric mismatch: reject

            ok = false;

            if issparse(roiData)
                roiData = full(roiData);
            end

            rawSize = size(roiData);

            sz = app.getDataSize();
            targetSizeYX = double(sz(1:2));

            sourceInfo.originalSizeYX = double(rawSize(1:2));
            sourceInfo.targetSizeYX = targetSizeYX;
            sourceInfo.bScaled = false;
            sourceInfo.scaleFactor = 1;

            [sizeMode, scaleFactor, msg] = app.getROIImportSizeCompatibilityFromSize(rawSize);

            switch lower(sizeMode)
                case 'exact'
                    ok = true;
                    return

                case 'scale'
                    promptText = sprintf([ ...
                        'The imported ROI variable has a different image size.\n\n' ...
                        'Imported size: Y = %.0f, X = %.0f\n' ...
                        'Current size:  Y = %.0f, X = %.0f\n\n' ...
                        'Symmetric scale factor: %.6g\n\n' ...
                        'Scale imported ROIs to the current image size?'], ...
                        sourceInfo.originalSizeYX(1), sourceInfo.originalSizeYX(2), ...
                        targetSizeYX(1), targetSizeYX(2), ...
                        scaleFactor);

                    try
                        choice = uiconfirm(app.UIFigure, promptText, ...
                            'Import ROI size mismatch', ...
                            'Options', {'Scale ROIs', 'Cancel'}, ...
                            'DefaultOption', 'Scale ROIs', ...
                            'CancelOption', 'Cancel', ...
                            'Icon', 'warning');
                    catch
                        choice = questdlg(promptText, ...
                            'Import ROI size mismatch', ...
                            'Scale ROIs', 'Cancel', 'Scale ROIs');

                        if isempty(choice)
                            choice = 'Cancel';
                        end
                    end

                    if ~strcmp(choice, 'Scale ROIs')
                        return
                    end

                    roiData = app.scaleImportedROIArrayToCurrentImage(roiData, scaleFactor);

                    sourceInfo.bScaled = true;
                    sourceInfo.scaleFactor = scaleFactor;
                    ok = true;

                otherwise
                    try
                        uialert(app.UIFigure, msg, ...
                            'Import ROI size mismatch', ...
                            'Icon', 'warning');
                    catch
                        warning('DataViewer:ImportROISizeMismatch', '%s', msg);
                    end

                    ok = false;
            end

        end

        function [sizeMode, scaleFactor, msg] = getROIImportSizeCompatibilityFromSize(app, varSize)
            %GETROIIMPORTSIZECOMPATIBILITYFROMSIZE Return import size compatibility.
            %
            %   sizeMode:
            %       'exact'        - source Y/X matches current image Y/X
            %       'scale'        - source Y/X is symmetrically scalable to current Y/X
            %       'incompatible' - source Y/X is not compatible

            sizeMode = 'incompatible';
            scaleFactor = 1;
            msg = '';

            if numel(varSize) < 2
                msg = 'Imported ROI variable must have at least two dimensions.';
                return
            end

            sourceY = double(varSize(1));
            sourceX = double(varSize(2));

            if ~isfinite(sourceY) || ~isfinite(sourceX) || sourceY <= 0 || sourceX <= 0
                msg = 'Imported ROI variable has invalid Y/X dimensions.';
                return
            end

            sz = app.getDataSize();
            targetY = double(sz(1));
            targetX = double(sz(2));

            if sourceY == targetY && sourceX == targetX
                sizeMode = 'exact';
                scaleFactor = 1;
                return
            end

            scaleY = targetY ./ sourceY;
            scaleX = targetX ./ sourceX;

            symmetricTol = max(1e-6, 1e-3 * max(abs([scaleX, scaleY])));

            if isfinite(scaleX) && isfinite(scaleY) && scaleX > 0 && scaleY > 0 && ...
                    abs(scaleX - scaleY) <= symmetricTol
                sizeMode = 'scale';
                scaleFactor = mean([scaleX, scaleY]);
                return
            end

            msg = sprintf([ ...
                'Imported ROI variable size [Y X] = [%.0f %.0f] does not match current ' ...
                'image size [Y X] = [%.0f %.0f], and symmetric scaling is not possible.\n\n' ...
                'Scale Y = %.6g, Scale X = %.6g'], ...
                sourceY, sourceX, targetY, targetX, scaleY, scaleX);

        end

        function roiDataOut = scaleImportedROIArrayToCurrentImage(app, roiDataIn, scaleFactor)
            %SCALEIMPORTEDROIARRAYTOCURRENTIMAGE Scale imported ROI masks/labels.
            %
            %   Uses nearest-neighbor interpolation so logical masks stay logical and label
            %   matrices preserve integer labels. The target size is the current image Y/X
            %   size. The scaleFactor input is kept for reporting/sanity but the target size
            %   is authoritative.

            if isempty(scaleFactor) || ~isfinite(scaleFactor) || scaleFactor <= 0
                error('DataViewer:InvalidImportROIScaleFactor', ...
                    'Invalid ROI import scale factor.');
            end

            if issparse(roiDataIn)
                roiDataIn = full(roiDataIn);
            end

            sz = app.getDataSize();
            targetNy = sz(1);
            targetNx = sz(2);

            if islogical(roiDataIn)
                if ndims(roiDataIn) <= 2
                    roiDataOut = iResizeLogicalMask2D(roiDataIn);
                    return
                end

                if ndims(roiDataIn) == 3
                    nMasks = size(roiDataIn, 3);
                    roiDataOut = false(targetNy, targetNx, nMasks);

                    for iMask = 1:nMasks
                        roiDataOut(:, :, iMask) = iResizeLogicalMask2D(roiDataIn(:, :, iMask));
                    end

                    return
                end

                error('DataViewer:UnsupportedImportedROIFormat', ...
                    'Logical ROI imports must be 2D masks or 3D mask stacks.');
            end

            if isnumeric(roiDataIn) && ismatrix(roiDataIn)
                labelMatrix = double(roiDataIn);

                if any(~isfinite(labelMatrix(:))) || any(labelMatrix(:) < 0) || ...
                        any(abs(labelMatrix(:) - round(labelMatrix(:))) > 1e-9)
                    error('DataViewer:InvalidImportedLabelMatrix', ...
                        'Imported label matrix must contain finite non-negative integer labels before scaling.');
                end

                roiDataOut = iResizeNumericNearest2D(labelMatrix);
                roiDataOut = round(double(roiDataOut));
                return
            end

            error('DataViewer:UnsupportedImportedROIFormat', ...
                'Unsupported ROI import variable. Use 2D/3D logical arrays or a 2D label matrix.');

            function outMask = iResizeLogicalMask2D(inMask)
                outMask = iResizeNumericNearest2D(logical(inMask));
                outMask = logical(outMask);
            end

            function out = iResizeNumericNearest2D(in)
                try
                    out = imresize(in, [targetNy, targetNx], 'nearest');
                    return
                catch
                    % Fallback without imresize. This uses inverse pixel-center mapping
                    % consistent with the boundary-preserving scaling convention used by
                    % .roi geometry scaling.
                    oldNy = size(in, 1);
                    oldNx = size(in, 2);

                    scaleY = targetNy ./ oldNy;
                    scaleX = targetNx ./ oldNx;

                    ySrc = round(((1:targetNy) - 0.5) ./ scaleY + 0.5);
                    xSrc = round(((1:targetNx) - 0.5) ./ scaleX + 0.5);

                    ySrc = min(max(ySrc, 1), oldNy);
                    xSrc = min(max(xSrc, 1), oldNx);

                    out = in(ySrc, xSrc);
                end
            end

        end

        function exportROIDataFromDialog(app)
            %EXPORTROIDATAFROMDIALOG Export selected ROI data products to CSV files.
            %
            %   The user chooses one base filename. Each selected export product appends its
            %   own suffix to that base:
            %       <base>_ROI_spatialMeasurements.csv
            %       <base>_ROI_normalTimeVectors.csv
            %       <base>_ROI_eventTimeVectors.csv
            %
            %   Trace exports use long-table CSV format. The spatial mean trace column is
            %   named "trace" for both normal and event exports.

            if ~app.hasData()
                app.setStatusMessage('Load image data before exporting ROI data.');
                return
            end

            if isempty(app.ROIList)
                app.setStatusMessage('No ROIs available to export.');
                return
            end

            opts = app.promptROIDataExportOptions();

            if opts.bCancelled
                app.setStatusMessage('Export ROI data cancelled.');
                return
            end

            if ~(opts.bSpatialMeasurements || opts.bNormalTimeVectors || opts.bEventTimeVectors)
                app.setStatusMessage('Export ROI data cancelled: no export option selected.');
                return
            end

            [baseFilePath, bConfirmed] = app.selectROIDataExportBaseFile();

            if ~bConfirmed
                app.setStatusMessage('Export ROI data cancelled.');
                return
            end

            [outFolder, baseFileName, ~] = fileparts(baseFilePath);

            if isempty(outFolder)
                outFolder = pwd;
            end

            if isempty(baseFileName)
                app.setStatusMessage('Export ROI data cancelled: filename cannot be empty.');
                return
            end

            exportedFiles = {};
            nRows = [];

            try
                if opts.bSpatialMeasurements
                    outFile = fullfile(outFolder, sprintf('%s_ROI_spatialMeasurements.csv', baseFileName));
                    nRows(end+1) = app.exportROISpatialMeasurementsCSV(outFile); %#ok<AGROW>
                    exportedFiles{end+1} = outFile; %#ok<AGROW>
                end

                if opts.bNormalTimeVectors
                    outFile = fullfile(outFolder, sprintf('%s_ROI_normalTimeVectors.csv', baseFileName));
                    nRows(end+1) = app.exportROINormalTimeVectorsCSV(outFile); %#ok<AGROW>
                    exportedFiles{end+1} = outFile; %#ok<AGROW>
                end

                if opts.bEventTimeVectors
                    outFile = fullfile(outFolder, sprintf('%s_ROI_eventTimeVectors.csv', baseFileName));
                    nRows(end+1) = app.exportROIEventTimeVectorsCSV(outFile); %#ok<AGROW>
                    exportedFiles{end+1} = outFile; %#ok<AGROW>
                end

            catch ME
                app.setStatusMessage(sprintf('Export ROI data failed: %s', ME.message));
                rethrow(ME)
            end

            if isempty(exportedFiles)
                app.setStatusMessage('No ROI data files were exported.');
                return
            end

            app.setStatusMessage(sprintf( ...
                'Exported %d ROI data CSV file(s) using base "%s". Total rows: %d.', ...
                numel(exportedFiles), fullfile(outFolder, baseFileName), sum(nRows)));

        end

        function opts = promptROIDataExportOptions(app)
            %PROMPTROIDATAEXPORTOPTIONS Modal checkbox dialog for ROI CSV export options.

            opts = struct();
            opts.bCancelled = true;
            opts.bSpatialMeasurements = false;
            opts.bNormalTimeVectors = false;
            opts.bEventTimeVectors = false;

            bEventsAvailable = app.hasNormalizedEvents();

            dlg = uifigure( ...
                'Name', 'Export ROI data', ...
                'WindowStyle', 'modal', ...
                'Position', [100 100 430 245], ...
                'Visible', 'off', ...
                'CloseRequestFcn', @onCancel);

            grid = uigridlayout(dlg);
            grid.RowHeight = {28, 34, 34, 34, '1x', 36};
            grid.ColumnWidth = {'1x'};
            grid.Padding = [12 12 12 12];

            titleLabel = uilabel(grid);
            titleLabel.Text = 'Select ROI data to export as CSV files:';
            titleLabel.FontWeight = 'bold';
            titleLabel.Layout.Row = 1;
            titleLabel.Layout.Column = 1;

            spatialCB = uicheckbox(grid);
            spatialCB.Text = 'Spatial measurements from current frame';
            spatialCB.Value = true;
            spatialCB.Layout.Row = 2;
            spatialCB.Layout.Column = 1;

            normalCB = uicheckbox(grid);
            normalCB.Text = 'ROI time vectors - normal mode';
            normalCB.Value = true;
            normalCB.Layout.Row = 3;
            normalCB.Layout.Column = 1;

            eventCB = uicheckbox(grid);
            eventCB.Text = 'ROI time vectors - event mode';
            eventCB.Value = bEventsAvailable;
            eventCB.Enable = matlab.lang.OnOffSwitchState(bEventsAvailable);
            eventCB.Layout.Row = 4;
            eventCB.Layout.Column = 1;

            statusLabel = uilabel(grid);
            statusLabel.Text = '';
            statusLabel.FontColor = [0.65 0 0];
            statusLabel.Layout.Row = 5;
            statusLabel.Layout.Column = 1;

            if ~bEventsAvailable
                statusLabel.Text = 'Event-mode export disabled: no event metadata available.';
            end

            buttonGrid = uigridlayout(grid);
            buttonGrid.RowHeight = {'1x'};
            buttonGrid.ColumnWidth = {'1x', '1x'};
            buttonGrid.Padding = [0 0 0 0];
            buttonGrid.Layout.Row = 6;
            buttonGrid.Layout.Column = 1;

            okButton = uibutton(buttonGrid, 'push');
            okButton.Text = 'OK';
            okButton.Layout.Row = 1;
            okButton.Layout.Column = 1;
            okButton.ButtonPushedFcn = @onOK;

            cancelButton = uibutton(buttonGrid, 'push');
            cancelButton.Text = 'Cancel';
            cancelButton.Layout.Row = 1;
            cancelButton.Layout.Column = 2;
            cancelButton.ButtonPushedFcn = @onCancel;

            if exist('placeAppInsideCaller', 'file') == 2
                placeAppInsideCaller(app, dlg, 'center');
            end
            dlg.Visible = 'on';

            uiwait(dlg);

            if isvalid(dlg)
                delete(dlg);
            end

            function onOK(~, ~)
                if ~(spatialCB.Value || normalCB.Value || eventCB.Value)
                    statusLabel.Text = 'Select at least one export option.';
                    return
                end

                opts.bCancelled = false;
                opts.bSpatialMeasurements = logical(spatialCB.Value);
                opts.bNormalTimeVectors = logical(normalCB.Value);
                opts.bEventTimeVectors = logical(eventCB.Value) && bEventsAvailable;

                uiresume(dlg);
            end

            function onCancel(~, ~)
                opts.bCancelled = true;
                uiresume(dlg);
            end

        end

        function baseName = makeROIDataExportBaseName(app, timestamp)
            %MAKEROIDATAEXPORTBASENAME Build safe base name for ROI CSV exports.

            if nargin < 2 || isempty(timestamp)
                timestamp = datestr(now, 'yyyymmdd_HHMMSS');
            end

            if isempty(app.CurrentFile)
                baseName = 'ROIData';
            else
                [~, baseName] = fileparts(app.CurrentFile);
            end

            if ~isempty(app.CurrentEntry)
                entryName = matlab.lang.makeValidName(app.CurrentEntry);
                baseName = sprintf('%s_%s', baseName, entryName);
            end

            baseName = sprintf('%s_%s', matlab.lang.makeValidName(baseName), timestamp);

        end

        function [baseFilePath, bConfirmed] = selectROIDataExportBaseFile(app)
            %SELECTROIDATAEXPORTBASEFILE Ask user for ROI data export base filename.
            %
            %   The selected filename is used only as a base. Its extension is removed and
            %   each selected export product appends its own suffix before ".csv".

            baseFilePath = '';
            bConfirmed = false;

            startFolder = pwd;
            if ~isempty(app.CurrentFile)
                startFolder = fileparts(app.CurrentFile);
            end

            timestamp = datestr(now, 'yyyymmdd_HHMMSS');
            defaultBaseName = app.makeROIDataExportBaseName(timestamp);
            defaultFileName = sprintf('%s.csv', defaultBaseName);

            [fileName, folderName] = uiputfile( ...
                {'*.csv', 'CSV files (*.csv)'}, ...
                'Set ROI data export filename', ...
                fullfile(startFolder, defaultFileName));

            if isequal(fileName, 0)
                return
            end

            [~, baseName, ~] = fileparts(fileName);

            baseName = strtrim(char(string(baseName)));

            if isempty(baseName)
                return
            end

            baseFilePath = fullfile(folderName, baseName);
            bConfirmed = true;
        end

        function nRows = exportROISpatialMeasurementsCSV(app, outFile)
            %EXPORTROISPATIALMEASUREMENTSCSV Export current ROI table/stat data.

            % Ensure spatial stats reflect the currently displayed frame.
            app.updateROIStatsForCurrentFrame();

            T = app.buildROISpatialMeasurementsTable();

            writetable(T, outFile);

            nRows = height(T);

        end

        function T = buildROISpatialMeasurementsTable(app)
            %BUILDROISPATIALMEASUREMENTSTABLE Build long metadata table for current ROI stats.

            nROI = numel(app.ROIList);

            if nROI == 0
                T = table();
                return
            end

            ctx = app.getCurrentROIExportContext();

            sourceFile = strings(nROI, 1);
            sourceEntry = strings(nROI, 1);
            viewMode = strings(nROI, 1);
            conditionName = strings(nROI, 1);
            repetitionLabel = strings(nROI, 1);
            roiName = strings(nROI, 1);
            roiType = strings(nROI, 1);
            colorRGB = strings(nROI, 1);

            frameIndex = nan(nROI, 1);
            frameTime = nan(nROI, 1);
            conditionID = nan(nROI, 1);
            repetitionID = nan(nROI, 1);
            roiID = nan(nROI, 1);
            visible = false(nROI, 1);
            selected = false(nROI, 1);

            centroidX_px = nan(nROI, 1);
            centroidY_px = nan(nROI, 1);
            centroidX_mm = nan(nROI, 1);
            centroidY_mm = nan(nROI, 1);
            distanceFromOrigin_px = nan(nROI, 1);
            distanceFromOrigin_mm = nan(nROI, 1);
            areaPx2 = nan(nROI, 1);
            areaMM2 = nan(nROI, 1);
            nPixels = nan(nROI, 1);
            spatialMean = nan(nROI, 1);
            spatialStd = nan(nROI, 1);
            spatialMedian = nan(nROI, 1);
            spatialMin = nan(nROI, 1);
            spatialMax = nan(nROI, 1);

            for iROI = 1:nROI
                ROI = app.ROIList(iROI);

                sourceFile(iROI) = string(ctx.sourceFile);
                sourceEntry(iROI) = string(ctx.sourceEntry);
                viewMode(iROI) = string(ctx.viewMode);
                frameIndex(iROI) = ctx.frameIndex;
                frameTime(iROI) = ctx.frameTime;
                conditionName(iROI) = string(ctx.conditionName);
                conditionID(iROI) = ctx.conditionID;
                repetitionID(iROI) = ctx.repetitionID;
                repetitionLabel(iROI) = string(ctx.repetitionLabel);

                if isfield(ROI, 'ID') && ~isempty(ROI.ID) && isfinite(ROI.ID)
                    roiID(iROI) = ROI.ID;
                end

                if isfield(ROI, 'name') && ~isempty(ROI.name)
                    roiName(iROI) = string(ROI.name);
                end

                if isfield(ROI, 'type') && ~isempty(ROI.type)
                    roiType(iROI) = string(ROI.type);
                end

                if isfield(ROI, 'color') && numel(ROI.color) == 3
                    c = double(ROI.color(:).');
                    colorRGB(iROI) = sprintf('%.6g,%.6g,%.6g', c(1), c(2), c(3));
                end

                if isfield(ROI, 'runtime')
                    if isfield(ROI.runtime, 'visible') && ~isempty(ROI.runtime.visible)
                        visible(iROI) = logical(ROI.runtime.visible);
                    end
                    if isfield(ROI.runtime, 'selected') && ~isempty(ROI.runtime.selected)
                        selected(iROI) = logical(ROI.runtime.selected);
                    end
                end

                if ~isfield(ROI, 'stats') || ~isstruct(ROI.stats)
                    continue
                end

                stats = ROI.stats;

                [centroidX_px(iROI), centroidY_px(iROI)] = iGetPair(stats, 'centroidXY_px');
                [centroidX_mm(iROI), centroidY_mm(iROI)] = iGetPair(stats, 'centroidXY_mm');

                distanceFromOrigin_px(iROI) = iGetScalar(stats, 'distanceFromOrigin_px');
                distanceFromOrigin_mm(iROI) = iGetScalar(stats, 'distanceFromOrigin_mm');
                areaPx2(iROI) = iGetScalar(stats, 'areaPx2');
                areaMM2(iROI) = iGetScalar(stats, 'areaMM2');
                nPixels(iROI) = iGetScalar(stats, 'NPixels');
                spatialMean(iROI) = iGetScalar(stats, 'spatialMean');
                spatialStd(iROI) = iGetScalar(stats, 'spatialStd');
                spatialMedian(iROI) = iGetScalar(stats, 'spatialMedian');
                spatialMin(iROI) = iGetScalar(stats, 'spatialMin');
                spatialMax(iROI) = iGetScalar(stats, 'spatialMax');
            end

            T = table( ...
                sourceFile, sourceEntry, viewMode, frameIndex, frameTime, ...
                conditionName, conditionID, repetitionID, repetitionLabel, ...
                roiID, roiName, roiType, visible, selected, colorRGB, ...
                centroidX_px, centroidY_px, centroidX_mm, centroidY_mm, ...
                distanceFromOrigin_px, distanceFromOrigin_mm, areaPx2, areaMM2, ...
                nPixels, spatialMean, spatialStd, spatialMedian, spatialMin, spatialMax);

            function val = iGetScalar(S, fieldName)
                val = NaN;

                if ~isfield(S, fieldName) || isempty(S.(fieldName))
                    return
                end

                candidate = double(S.(fieldName));

                if isscalar(candidate) && isfinite(candidate)
                    val = candidate;
                end
            end

            function [x, y] = iGetPair(S, fieldName)
                x = NaN;
                y = NaN;

                if ~isfield(S, fieldName) || numel(S.(fieldName)) < 2
                    return
                end

                candidate = double(S.(fieldName)(:).');

                if numel(candidate) >= 2
                    x = candidate(1);
                    y = candidate(2);
                end
            end

        end

        function nRows = exportROINormalTimeVectorsCSV(app, outFile)
            %EXPORTROINORMALTIMEVECTORSCSV Export normal-mode ROI traces in long format.

            [roiMasks, roiIdxList] = app.collectROIMasksForTrace();

            if isempty(roiIdxList)
                T = table();
                writetable(T, outFile);
                nRows = 0;
                return
            end

            sourceSize = app.getSourceDataSize();
            nFrames = sourceSize(3);
            frameIdx = 1:nFrames;

            traceMatrix = app.computeROITraceMatrixForFrameVectorExport(roiMasks, frameIdx, []);
            xData = app.getNormalTraceTimeVectorForExport(nFrames);

            T = app.buildROINormalTraceExportTable(roiIdxList, traceMatrix, xData, frameIdx);

            writetable(T, outFile);

            nRows = height(T);

        end

        function T = buildROINormalTraceExportTable(app, roiIdxList, traceMatrix, xData, frameIdx)
            %BUILDROINORMALTRACEEXPORTTABLE Convert normal ROI trace matrix to long table.

            nROI = numel(roiIdxList);
            nFrames = numel(frameIdx);
            nTotal = nROI * nFrames;

            sourceFile = strings(nTotal, 1);
            sourceEntry = strings(nTotal, 1);
            viewMode = strings(nTotal, 1);
            roiName = strings(nTotal, 1);
            roiType = strings(nTotal, 1);

            frameIndex = nan(nTotal, 1);
            time_s = nan(nTotal, 1);
            roiID = nan(nTotal, 1);
            trace = nan(nTotal, 1);
            nPixels = nan(nTotal, 1);

            ctx = app.getCurrentROIExportContext();
            row = 0;

            for k = 1:nROI
                roiIdx = roiIdxList(k);
                ROI = app.ROIList(roiIdx);

                thisRoiID = NaN;
                if isfield(ROI, 'ID') && isfinite(ROI.ID)
                    thisRoiID = ROI.ID;
                end

                thisName = string(ROI.name);
                thisType = string(ROI.type);
                thisNPixels = sum(logical(ROI.mask(:)));

                for t = 1:nFrames
                    row = row + 1;

                    sourceFile(row) = string(ctx.sourceFile);
                    sourceEntry(row) = string(ctx.sourceEntry);
                    viewMode(row) = "normal";
                    frameIndex(row) = frameIdx(t);
                    time_s(row) = xData(t);
                    roiID(row) = thisRoiID;
                    roiName(row) = thisName;
                    roiType(row) = thisType;
                    trace(row) = traceMatrix(k, t);
                    nPixels(row) = thisNPixels;
                end
            end

            T = table( ...
                sourceFile, sourceEntry, viewMode, frameIndex, time_s, ...
                roiID, roiName, roiType, trace, nPixels);

        end

        function nRows = exportROIEventTimeVectorsCSV(app, outFile)
            %EXPORTROIEVENTTIMEVECTORSCSV Export all non-ignored event ROI traces.

            if ~app.hasNormalizedEvents()
                error('DataViewer:NoEventsForROIExport', ...
                    'Event-mode ROI export requires event metadata.');
            end

            [roiMasks, roiIdxList] = app.collectROIMasksForTrace();

            if isempty(roiIdxList)
                T = table();
                writetable(T, outFile);
                nRows = 0;
                return
            end

            switch lower(app.getSourceType())
                case 'dat'
                    T = app.buildDatROIEventTraceExportTable(roiMasks, roiIdxList);

                case 'umt'
                    T = app.buildUMTROIEventTraceExportTable(roiMasks, roiIdxList);

                otherwise
                    error('DataViewer:UnsupportedEventROIExport', ...
                        'Event ROI export is supported only for DAT and UMT sources.');
            end

            writetable(T, outFile);

            nRows = height(T);

        end

        function T = buildDatROIEventTraceExportTable(app, roiMasks, roiIdxList)
            %BUILDDATROIEVENTTRACEEXPORTTABLE Build DAT event-mode long trace table.
            %
            %   Exports every non-ignored event instance represented in normalized EventInfo.

            if isempty(app.EventSource) || ~app.hasNormalizedEvents()
                T = table();
                return
            end

            sourceSize = app.getSourceDataSize();
            datLength = sourceSize(3);

            eventInfo = app.EventInfo;
            nEvents = numel(eventInfo.EventIDs);

            Tall = table();

            for iEvent = 1:nEvents
                eventName = char(string(eventInfo.EventNamePerEvent{iEvent}));
                eventID = double(eventInfo.EventIDs(iEvent));
                repetitionID = double(eventInfo.RepetitionIndex(iEvent));

                try
                    [frMat, ~, repList] = app.EventSource.getFrameMatrix(datLength, eventName, repetitionID);
                catch
                    continue
                end

                if isempty(frMat)
                    continue
                end

                rowIdx = 1;

                if exist('repList', 'var') && ~isempty(repList)
                    idx = find(double(repList(:)) == repetitionID, 1, 'first');
                    if ~isempty(idx)
                        rowIdx = idx;
                    end
                end

                rowIdx = min(max(rowIdx, 1), size(frMat, 1));
                frameIdx = double(frMat(rowIdx, :));
                frameIdx = frameIdx(isfinite(frameIdx));
                frameIdx = frameIdx(frameIdx >= 1 & frameIdx <= datLength);

                if isempty(frameIdx)
                    continue
                end

                traceMatrix = app.computeROITraceMatrixForFrameVectorExport(roiMasks, frameIdx, []);
                xData = app.getEventTraceTimeVectorForExport(numel(frameIdx));

                thisTable = app.buildROIEventTraceRowsTable( ...
                    roiIdxList, traceMatrix, xData, frameIdx, ...
                    eventName, eventID, repetitionID, iEvent);

                Tall = [Tall; thisTable]; %#ok<AGROW>
            end

            T = Tall;

        end

        function T = buildUMTROIEventTraceExportTable(app, roiMasks, roiIdxList)
            %BUILDUMTROIEVENTTRACEEXPORTTABLE Build UMT event-mode long trace table.

            if ~app.hasNormalizedEvents()
                T = table();
                return
            end

            sourceSize = app.getSourceDataSize();
            nFrames = sourceSize(3);

            eventInfo = app.EventInfo;
            nEvents = numel(eventInfo.EventIDs);

            Tall = table();

            for eIdx = 1:nEvents
                eventName = char(string(eventInfo.EventNamePerEvent{eIdx}));
                eventID = double(eventInfo.EventIDs(eIdx));
                repetitionID = double(eventInfo.RepetitionIndex(eIdx));

                frameIdx = 1:nFrames;
                traceMatrix = app.computeROITraceMatrixForFrameVectorExport(roiMasks, frameIdx, eIdx);
                xData = app.getEventTraceTimeVectorForExport(nFrames);

                thisTable = app.buildROIEventTraceRowsTable( ...
                    roiIdxList, traceMatrix, xData, frameIdx, ...
                    eventName, eventID, repetitionID, eIdx);

                Tall = [Tall; thisTable]; %#ok<AGROW>
            end

            T = Tall;

        end

        function T = buildROIEventTraceRowsTable(app, roiIdxList, traceMatrix, xData, sourceFrameIdx, eventName, eventID, repetitionID, eventInstanceIndex)
            %BUILDROIEVENTTRACEROWSTABLE Convert one event ROI trace matrix to long table.

            nROI = numel(roiIdxList);
            nFrames = numel(sourceFrameIdx);
            nTotal = nROI * nFrames;

            sourceFile = strings(nTotal, 1);
            sourceEntry = strings(nTotal, 1);
            eventNameCol = strings(nTotal, 1);
            roiName = strings(nTotal, 1);
            roiType = strings(nTotal, 1);

            eventIDCol = nan(nTotal, 1);
            repetitionIDCol = nan(nTotal, 1);
            eventInstanceIndexCol = nan(nTotal, 1);
            eventFrameIndex = nan(nTotal, 1);
            sourceFrameIndex = nan(nTotal, 1);
            time_s = nan(nTotal, 1);
            roiID = nan(nTotal, 1);
            trace = nan(nTotal, 1);
            nPixels = nan(nTotal, 1);

            ctx = app.getCurrentROIExportContext();
            row = 0;

            for k = 1:nROI
                roiIdx = roiIdxList(k);
                ROI = app.ROIList(roiIdx);

                thisRoiID = NaN;
                if isfield(ROI, 'ID') && isfinite(ROI.ID)
                    thisRoiID = ROI.ID;
                end

                thisName = string(ROI.name);
                thisType = string(ROI.type);
                thisNPixels = sum(logical(ROI.mask(:)));

                for t = 1:nFrames
                    row = row + 1;

                    sourceFile(row) = string(ctx.sourceFile);
                    sourceEntry(row) = string(ctx.sourceEntry);
                    eventNameCol(row) = string(eventName);
                    eventIDCol(row) = eventID;
                    repetitionIDCol(row) = repetitionID;
                    eventInstanceIndexCol(row) = eventInstanceIndex;
                    eventFrameIndex(row) = t;
                    sourceFrameIndex(row) = sourceFrameIdx(t);
                    time_s(row) = xData(t);
                    roiID(row) = thisRoiID;
                    roiName(row) = thisName;
                    roiType(row) = thisType;
                    trace(row) = traceMatrix(k, t);
                    nPixels(row) = thisNPixels;
                end
            end

            T = table( ...
                sourceFile, sourceEntry, eventNameCol, eventIDCol, repetitionIDCol, ...
                eventInstanceIndexCol, eventFrameIndex, sourceFrameIndex, time_s, ...
                roiID, roiName, roiType, trace, nPixels, ...
                'VariableNames', {'sourceFile', 'sourceEntry', 'eventName', 'eventID', ...
                'repetitionID', 'eventInstanceIndex', 'eventFrameIndex', ...
                'sourceFrameIndex', 'time_s', 'roiID', 'roiName', ...
                'roiType', 'trace', 'nPixels'});

        end

        function traceMatrix = computeROITraceMatrixForFrameVectorExport(app, roiMasks, frameIdx, eventIndex)
            %COMPUTEROITRACEMATRIXFORFRAMEVECTOREXPORT Compute ROI means for explicit frames.
            %
            %   eventIndex:
            %       []      - call DataSource.getFrame(frameIdx)
            %       scalar  - call DataSource.getFrame(frameIdx, eventIndex), used for UMT E.

            frameIdx = double(frameIdx(:).');
            frameIdx = frameIdx(isfinite(frameIdx));

            nROI = numel(roiMasks);
            nFrames = numel(frameIdx);

            traceMatrix = nan(nROI, nFrames);

            if isempty(frameIdx)
                return
            end

            if nargin < 4
                eventIndex = [];
            end

            % DAT backend optimized path. This should be used only when no explicit UMT
            % event index is requested.
            if isempty(eventIndex) && strcmpi(app.getSourceType(), 'dat') && ...
                    ismethod(app.DataSource, 'getROIMeanTraceMatrix')
                try
                    [traceMatrix, ~] = app.DataSource.getROIMeanTraceMatrix(roiMasks, frameIdx);
                    traceMatrix = double(traceMatrix);
                    return
                catch
                    % Fall back to frame-by-frame correctness path.
                end
            end

            for t = 1:nFrames
                if isempty(eventIndex)
                    frame = app.DataSource.getFrame(frameIdx(t));
                else
                    frame = app.DataSource.getFrame(frameIdx(t), eventIndex);
                end

                traceMatrix(:, t) = app.computeROIMeansFromFrame(roiMasks, frame);
            end

        end

        function xData = getNormalTraceTimeVectorForExport(app, nSamples)
            %GETNORMALTRACETIMEVECTORFOREXPORT Return normal-mode trace times.

            nSamples = max(1, round(nSamples));
            frameRateHz = app.getFrameRateHz();

            if isempty(frameRateHz)
                xData = (1:nSamples).';
            else
                xData = ((0:nSamples-1) ./ frameRateHz).';
            end

        end

        function xData = getEventTraceTimeVectorForExport(app, nSamples)
            %GETEVENTTRACETIMEVECTORFOREXPORT Return event-aligned trace times.

            nSamples = max(1, round(nSamples));
            frameRateHz = app.getFrameRateHz();

            if isempty(frameRateHz)
                xData = (1:nSamples).';
                return
            end

            baselinePeriod = app.getActiveBaselinePeriod();

            if ~isempty(baselinePeriod)
                eventOnsetFrame = round(baselinePeriod * frameRateHz) + 1;
                xData = (((1:nSamples) - eventOnsetFrame) ./ frameRateHz).';
            else
                xData = ((0:nSamples-1) ./ frameRateHz).';
            end

        end

        function ctx = getCurrentROIExportContext(app)
            %GETCURRENTROIEXPORTCONTEXT Return source/view metadata for ROI exports.

            ctx = struct();
            ctx.sourceFile = app.CurrentFile;
            ctx.sourceEntry = app.CurrentEntry;
            ctx.viewMode = app.ViewMode;
            ctx.frameIndex = app.CurrentFrame;
            ctx.frameTime = app.getCurrentFrameTime(app.CurrentFrame);
            ctx.conditionName = '';
            ctx.conditionID = NaN;
            ctx.repetitionID = NaN;
            ctx.repetitionLabel = '';

            if strcmpi(app.ViewMode, 'event')
                ctx.conditionName = app.CurrentCondition;
                condID = app.getConditionID(app.CurrentCondition);

                if ~isempty(condID)
                    ctx.conditionID = condID;
                end

                ctx.repetitionLabel = app.CurrentRepetition;

                repID = str2double(app.CurrentRepetition);
                if isfinite(repID)
                    ctx.repetitionID = repID;
                end
            end

        end

        function openThresholdROICreatorDialog(app)
            %OPENTHRESHOLDROICREATORDIALOG Create ROIs from thresholded current frame.
            %
            %   The dialog thresholds the currently displayed image frame, intersects the
            %   candidate mask with the active logical mask, labels connected components,
            %   and appends the resulting polygon ROIs through the same import pipeline used
            %   by Import ROI.
            %
            %   The dialog uses nested grid layouts directly. No controls are placed inside
            %   an intermediate panel, which keeps layout behavior more predictable across
            %   MATLAB releases and display scaling settings.

            if ~app.hasData()
                app.setStatusMessage('Load image data before creating ROIs by threshold.');
                return
            end

            rawImage = single(app.getCurrentFrame());
            activeMask = app.getActiveLogicalMask();

            if isempty(activeMask) || ~isequal(size(activeMask), size(rawImage))
                activeMask = true(size(rawImage));
            else
                activeMask = logical(activeMask);
            end

            [imageOriginal, ok, msg] = iSanitizeImage(rawImage, activeMask);

            if ~ok
                app.setStatusMessage(msg);
                try
                    uialert(app.UIFigure, msg, 'Threshold ROI creation failed', 'Icon', 'warning');
                catch
                    warning('DataViewer:ThresholdROICreationFailed', '%s', msg);
                end
                return
            end

            [Ny, Nx] = size(imageOriginal);
            validMask = activeMask & isfinite(imageOriginal);

            if ~any(validMask(:))
                app.setStatusMessage('Threshold ROI creation failed: no finite pixels inside the logical mask.');
                return
            end

            originalLimits = double([min(imageOriginal(validMask), [], 'all'), ...
                max(imageOriginal(validMask), [], 'all')]);

            if originalLimits(1) == originalLimits(2)
                originalLimits = originalLimits + [-1 1];
            end

            dlg = uifigure( ...
                'Name', 'Create ROIs by threshold', ...
                'WindowStyle', 'modal', ...
                'Position', [100 100 1040 640], ...
                'Visible', 'off', ...
                'CloseRequestFcn', @onCancel);

            mainGrid = uigridlayout(dlg);
            mainGrid.RowHeight = {28, '1x', 40};
            mainGrid.ColumnWidth = {'1x', '1x', 310};
            mainGrid.Padding = [12 12 12 12];
            mainGrid.RowSpacing = 8;
            mainGrid.ColumnSpacing = 10;

            titleLabel = uilabel(mainGrid);
            titleLabel.Text = 'Create ROIs by threshold from the current displayed frame';
            titleLabel.FontWeight = 'bold';
            titleLabel.Layout.Row = 1;
            titleLabel.Layout.Column = [1 3];

            imgAx = uiaxes(mainGrid);
            imgAx.Layout.Row = 2;
            imgAx.Layout.Column = 1;
            imgAx.Box = 'on';
            imgAx.YDir = 'reverse';
            imgAx.DataAspectRatio = [1 1 1];
            title(imgAx, 'Image + threshold overlay');

            maskAx = uiaxes(mainGrid);
            maskAx.Layout.Row = 2;
            maskAx.Layout.Column = 2;
            maskAx.Box = 'on';
            maskAx.YDir = 'reverse';
            maskAx.DataAspectRatio = [1 1 1];
            title(maskAx, 'Labeled candidate ROIs');

            % Controls are placed directly in a nested gridlayout instead of inside a
            % panel. This avoids panel-title spacing and improves resizing behavior.
            controlGrid = uigridlayout(mainGrid);
            controlGrid.Layout.Row = 2;
            controlGrid.Layout.Column = 3;
            controlGrid.RowHeight = {24, 28, 28, 28, 24, 58, 54, 28, 28, 28, 36, '1x'};
            controlGrid.ColumnWidth = {110, '1x'};
            controlGrid.Padding = [4 0 4 0];
            controlGrid.RowSpacing = 7;
            controlGrid.ColumnSpacing = 8;

            controlsTitle = uilabel(controlGrid);
            controlsTitle.Text = 'Controls';
            controlsTitle.FontWeight = 'bold';
            controlsTitle.Layout.Row = 1;
            controlsTitle.Layout.Column = [1 2];

            nameLabel = uilabel(controlGrid);
            nameLabel.Text = 'Base name';
            nameLabel.Layout.Row = 2;
            nameLabel.Layout.Column = 1;

            nameField = uieditfield(controlGrid, 'text');
            nameField.Value = 'ThresholdROI';
            nameField.Layout.Row = 2;
            nameField.Layout.Column = 2;

            normalizeCheck = uicheckbox(controlGrid);
            normalizeCheck.Text = 'Normalize image values';
            normalizeCheck.Value = false;
            normalizeCheck.Layout.Row = 3;
            normalizeCheck.Layout.Column = [1 2];
            normalizeCheck.ValueChangedFcn = @onNormalizeChanged;

            autoButton = uibutton(controlGrid, 'push');
            autoButton.Text = 'Auto set';
            autoButton.Tooltip = {'Set threshold to mean ± 2 SD inside the active logical mask.'};
            autoButton.Layout.Row = 4;
            autoButton.Layout.Column = 1;
            autoButton.ButtonPushedFcn = @onAutoSet;

            thresholdField = uieditfield(controlGrid, 'numeric');
            thresholdField.ValueDisplayFormat = '%.6g';
            thresholdField.Layout.Row = 4;
            thresholdField.Layout.Column = 2;
            thresholdField.ValueChangedFcn = @onThresholdFieldChanged;

            thresholdLabel = uilabel(controlGrid);
            thresholdLabel.Text = 'Threshold level';
            thresholdLabel.FontWeight = 'bold';
            thresholdLabel.Layout.Row = 5;
            thresholdLabel.Layout.Column = [1 2];

            thresholdSlider = uislider(controlGrid);
            thresholdSlider.Layout.Row = 6;
            thresholdSlider.Layout.Column = [1 2];
            thresholdSlider.ValueChangingFcn = @onThresholdSliderChanging;
            thresholdSlider.ValueChangedFcn = @onThresholdSliderChanged;

            directionGroup = uibuttongroup(controlGrid);
            directionGroup.Title = 'Selection direction';
            directionGroup.Layout.Row = 7;
            directionGroup.Layout.Column = [1 2];
            directionGroup.SelectionChangedFcn = @onControlChanged;

            aboveButton = uiradiobutton(directionGroup);
            aboveButton.Text = 'above threshold';
            aboveButton.Position = [12 6 125 22];
            aboveButton.Value = true;

            belowButton = uiradiobutton(directionGroup);
            belowButton.Text = 'below threshold';
            belowButton.Position = [148 6 125 22];

            minPixelLabel = uilabel(controlGrid);
            minPixelLabel.Text = 'Min pixels';
            minPixelLabel.Layout.Row = 8;
            minPixelLabel.Layout.Column = 1;

            minPixelSpinner = uispinner(controlGrid);
            minPixelSpinner.Limits = [1, max(1, numel(imageOriginal))];
            minPixelSpinner.RoundFractionalValues = 'on';
            minPixelSpinner.Value = 1;
            minPixelSpinner.Tooltip = {'Remove connected candidate regions with fewer pixels than this value.'};
            minPixelSpinner.Layout.Row = 8;
            minPixelSpinner.Layout.Column = 2;
            minPixelSpinner.ValueChangedFcn = @onControlChanged;

            removeHolesCheck = uicheckbox(controlGrid);
            removeHolesCheck.Text = 'Remove holes';
            removeHolesCheck.Value = true;
            removeHolesCheck.Tooltip = {'Recommended. Hole-containing ROIs are rejected by the current DataViewer ROI model.'};
            removeHolesCheck.Layout.Row = 9;
            removeHolesCheck.Layout.Column = [1 2];
            removeHolesCheck.ValueChangedFcn = @onControlChanged;

            smoothCheck = uicheckbox(controlGrid);
            smoothCheck.Text = 'Smooth borders / remove narrow bridges';
            smoothCheck.Value = false;
            smoothCheck.Tooltip = {'Applies bwmorph(mask, ''open'', Inf) to each candidate region.'};
            smoothCheck.Layout.Row = 10;
            smoothCheck.Layout.Column = [1 2];
            smoothCheck.ValueChangedFcn = @onControlChanged;

            countLabel = uilabel(controlGrid);
            countLabel.Text = 'Candidate ROIs: 0';
            countLabel.FontWeight = 'bold';
            countLabel.Layout.Row = 11;
            countLabel.Layout.Column = [1 2];

            infoLabel = uilabel(controlGrid);
            infoLabel.Text = '';
            infoLabel.WordWrap = 'on';
            infoLabel.FontColor = [0.4 0.4 0.4];
            infoLabel.Layout.Row = 12;
            infoLabel.Layout.Column = [1 2];

            buttonGrid = uigridlayout(mainGrid);
            buttonGrid.RowHeight = {'1x'};
            buttonGrid.ColumnWidth = {'1x', 125, 125};
            buttonGrid.Padding = [0 0 0 0];
            buttonGrid.ColumnSpacing = 8;
            buttonGrid.Layout.Row = 3;
            buttonGrid.Layout.Column = [1 3];

            statusLabel = uilabel(buttonGrid);
            statusLabel.Text = '';
            statusLabel.FontColor = [0.65 0 0];
            statusLabel.Layout.Row = 1;
            statusLabel.Layout.Column = 1;

            createButton = uibutton(buttonGrid, 'push');
            createButton.Text = 'Create ROIs';
            createButton.FontWeight = 'bold';
            createButton.Layout.Row = 1;
            createButton.Layout.Column = 2;
            createButton.ButtonPushedFcn = @onCreateROIs;

            cancelButton = uibutton(buttonGrid, 'push');
            cancelButton.Text = 'Cancel';
            cancelButton.Layout.Row = 1;
            cancelButton.Layout.Column = 3;
            cancelButton.ButtonPushedFcn = @onCancel;

            hold(imgAx, 'on');
            imageHandle = imagesc(imgAx, imageOriginal);
            overlayHandle = image(imgAx, ...
                'CData', zeros(Ny, Nx, 3), ...
                'AlphaData', zeros(Ny, Nx), ...
                'XData', [1 Nx], ...
                'YData', [1 Ny], ...
                'HitTest', 'off', ...
                'PickableParts', 'none');
            colormap(imgAx, 'gray');
            axis(imgAx, 'image');
            imgAx.YDir = 'reverse';

            maskHandle = image(maskAx, ...
                'CData', zeros(Ny, Nx, 3), ...
                'XData', [1 Nx], ...
                'YData', [1 Ny], ...
                'HitTest', 'off', ...
                'PickableParts', 'none');
            axis(maskAx, 'image');
            maskAx.YDir = 'reverse';

            labelTextHandles = gobjects(0);
            currentLabelMatrix = zeros(Ny, Nx, 'uint16');

            iSetThresholdLimitsAndValue(mean(originalLimits));
            onAutoSet();

            try
                if exist('placeAppInsideCaller', 'file') == 2
                    placeAppInsideCaller(app, dlg, 'center');
                end
            catch
            end

            dlg.Visible = 'on';
            uiwait(dlg);

            if isvalid(dlg)
                delete(dlg);
            end

            % ---------------------------------------------------------------------
            % Nested callbacks
            % ---------------------------------------------------------------------
            function onNormalizeChanged(~, ~)
                oldValue = double(thresholdField.Value);
                oldLimits = double(thresholdSlider.Limits);

                if diff(oldLimits) <= 0
                    frac = 0.5;
                else
                    frac = (oldValue - oldLimits(1)) ./ diff(oldLimits);
                    frac = min(max(frac, 0), 1);
                end

                imgNow = iGetDisplayImage();
                limNow = iGetDisplayLimits(imgNow);

                newValue = limNow(1) + frac .* diff(limNow);
                iSetThresholdLimitsAndValue(newValue);
                iUpdatePreview(true);
            end

            function onAutoSet(~, ~)
                imgNow = iGetDisplayImage();
                maskedValues = double(imgNow(validMask));
                maskedValues = maskedValues(isfinite(maskedValues));

                if isempty(maskedValues)
                    statusLabel.Text = 'No finite pixels inside the active logical mask.';
                    return
                end

                mu = mean(maskedValues, 'omitnan');
                sigma = std(maskedValues, 0, 'omitnan');

                if aboveButton.Value
                    thr = mu + 2 .* sigma;
                else
                    thr = mu - 2 .* sigma;
                end

                limNow = iGetDisplayLimits(imgNow);
                thr = min(max(thr, limNow(1)), limNow(2));

                iSetThresholdLimitsAndValue(thr);
                iUpdatePreview(true);
            end

            function onThresholdFieldChanged(~, ~)
                value = double(thresholdField.Value);
                lim = double(thresholdSlider.Limits);
                value = min(max(value, lim(1)), lim(2));

                thresholdField.Value = value;
                thresholdSlider.Value = value;

                iUpdatePreview(true);
            end

            function onThresholdSliderChanging(~, event)
                value = double(event.Value);
                thresholdField.Value = value;
                iUpdatePreview(false);
            end

            function onThresholdSliderChanged(~, ~)
                thresholdField.Value = double(thresholdSlider.Value);
                iUpdatePreview(true);
            end

            function onControlChanged(~, ~)
                iUpdatePreview(true);
            end

            function onCreateROIs(~, ~)
                labelMatrix = iBuildLabelMatrix(true);

                if isempty(labelMatrix) || ~any(labelMatrix(:))
                    statusLabel.Text = 'No ROIs found. Adjust threshold or minimum-pixel settings.';
                    try
                        uialert(dlg, statusLabel.Text, 'No ROIs found', 'Icon', 'warning');
                    catch
                    end
                    return
                end

                baseName = strtrim(char(string(nameField.Value)));
                if isempty(baseName)
                    baseName = 'ThresholdROI';
                end

                sourceInfo = struct();
                sourceInfo.variableName = baseName;
                sourceInfo.sourceLabel = 'threshold ROI creator';
                sourceInfo.sourceType = 'threshold';
                sourceInfo.bScaled = false;
                sourceInfo.originalSizeYX = double([Ny Nx]);
                sourceInfo.targetSizeYX = double([Ny Nx]);
                sourceInfo.scaleFactor = 1;

                [newROIList, report] = app.buildImportedROIListFromArray(double(labelMatrix), sourceInfo);

                if isempty(newROIList)
                    statusText = app.formatROIImportStatusText(0, report, sourceInfo);
                    statusLabel.Text = statusText;
                    app.setStatusMessage(statusText);
                    try
                        uialert(dlg, statusText, 'No valid ROIs created', 'Icon', 'warning');
                    catch
                    end
                    return
                end

                app.appendImportedROIsAndRefresh(newROIList);

                statusText = app.formatROIImportStatusText(numel(newROIList), report, sourceInfo);
                app.setStatusMessage(statusText);

                uiresume(dlg);
            end

            function onCancel(~, ~)
                uiresume(dlg);
            end

            % ---------------------------------------------------------------------
            % Nested utilities
            % ---------------------------------------------------------------------
            function iSetThresholdLimitsAndValue(value)
                imgNow = iGetDisplayImage();
                limNow = iGetDisplayLimits(imgNow);

                thresholdSlider.Limits = limNow;
                thresholdField.Limits = limNow;

                value = min(max(double(value), limNow(1)), limNow(2));

                thresholdSlider.Value = value;
                thresholdField.Value = value;

                imageHandle.CData = imgNow;
                imgAx.CLim = limNow;
            end

            function limNow = iGetDisplayLimits(imgNow)
                vals = double(imgNow(validMask));
                vals = vals(isfinite(vals));

                if isempty(vals)
                    limNow = [0 1];
                    return
                end

                limNow = [min(vals, [], 'all'), max(vals, [], 'all')];

                if limNow(1) == limNow(2)
                    limNow = limNow + [-1 1];
                end
            end

            function imgNow = iGetDisplayImage()
                if ~normalizeCheck.Value
                    imgNow = imageOriginal;
                    return
                end

                denom = originalLimits(2) - originalLimits(1);

                if denom <= 0 || ~isfinite(denom)
                    imgNow = zeros(size(imageOriginal), 'single');
                else
                    imgNow = single((double(imageOriginal) - originalLimits(1)) ./ denom);
                end
            end

            function iUpdatePreview(bApplyMorphology)
                if ~isvalid(dlg)
                    return
                end

                currentLabelMatrix = iBuildLabelMatrix(bApplyMorphology);

                rgbMask = iLabelMatrixToRGB(currentLabelMatrix);
                overlayHandle.CData = rgbMask;
                overlayHandle.AlphaData = 0.45 .* double(currentLabelMatrix > 0);
                maskHandle.CData = rgbMask;

                iUpdateLabelText(currentLabelMatrix);

                nRegions = double(max(currentLabelMatrix(:)));
                nHoleRegions = iCountHoleRegions(currentLabelMatrix);

                if nHoleRegions > 0
                    countLabel.Text = sprintf('Candidate ROIs: %d (%d with holes)', nRegions, nHoleRegions);
                    infoLabel.Text = ['Hole-containing regions will be skipped at creation time. ' ...
                        'Enable "Remove holes" to preserve those regions.'];
                else
                    countLabel.Text = sprintf('Candidate ROIs: %d', nRegions);
                    infoLabel.Text = sprintf('Logical mask active: %s', string(app.hasUserLogicalMask()));
                end

                title(maskAx, sprintf('Labeled candidate ROIs: %d', nRegions));
                drawnow limitrate
            end

            function labelMatrix = iBuildLabelMatrix(bApplyMorphology)
                imgNow = iGetDisplayImage();
                thr = double(thresholdField.Value);

                if aboveButton.Value
                    BW = imgNow > thr;
                else
                    BW = imgNow < thr;
                end

                BW = logical(BW) & activeMask;

                minPixels = max(1, round(double(minPixelSpinner.Value)));

                try
                    BW = bwareaopen(BW, minPixels, 8);
                catch
                    BW = iRemoveSmallComponents(BW, minPixels);
                end

                labelMatrix = bwlabel(BW, 8);

                if ~bApplyMorphology || ...
                        ~(removeHolesCheck.Value || smoothCheck.Value) || ...
                        ~any(labelMatrix(:))
                    return
                end

                mergedMask = false(size(labelMatrix));
                nLabels = max(labelMatrix(:));

                for iLabel = 1:nLabels
                    componentMask = labelMatrix == iLabel;

                    if removeHolesCheck.Value
                        try
                            componentMask = imfill(componentMask, 'holes');
                        catch
                        end
                    end

                    if smoothCheck.Value
                        try
                            componentMask = bwmorph(componentMask, 'open', Inf);
                        catch
                        end
                    end

                    mergedMask = mergedMask | componentMask;
                end

                try
                    mergedMask = bwareaopen(mergedMask, minPixels, 8);
                catch
                    mergedMask = iRemoveSmallComponents(mergedMask, minPixels);
                end

                labelMatrix = bwlabel(mergedMask, 8);
            end

            function BWout = iRemoveSmallComponents(BWin, minPixels)
                BWout = false(size(BWin));
                L = bwlabel(logical(BWin), 8);
                n = max(L(:));

                for iComp = 1:n
                    comp = L == iComp;
                    if nnz(comp) >= minPixels
                        BWout = BWout | comp;
                    end
                end
            end

            function rgb = iLabelMatrixToRGB(labelMatrix)
                labelMatrix = double(labelMatrix);
                nLabels = max(labelMatrix(:));

                rgb = zeros([size(labelMatrix), 3], 'single');

                if nLabels < 1
                    return
                end

                try
                    colors = distinguishable_colors(nLabels, {'w', 'k'});
                catch
                    colors = lines(nLabels);
                end

                colors = min(max(double(colors), 0), 1);

                for iLabel = 1:nLabels
                    idx = labelMatrix == iLabel;

                    for iChannel = 1:3
                        tmp = rgb(:, :, iChannel);
                        tmp(idx) = single(colors(iLabel, iChannel));
                        rgb(:, :, iChannel) = tmp;
                    end
                end
            end

            function iUpdateLabelText(labelMatrix)
                if ~isempty(labelTextHandles)
                    try
                        delete(labelTextHandles(isvalid(labelTextHandles)));
                    catch
                    end
                end

                labelTextHandles = gobjects(0);

                nLabels = max(labelMatrix(:));
                if nLabels < 1
                    return
                end

                try
                    props = regionprops(labelMatrix, 'Centroid');
                catch
                    props = struct('Centroid', {});
                end

                hold(maskAx, 'on');

                for iLabel = 1:min(nLabels, numel(props))
                    centroid = props(iLabel).Centroid;

                    if isempty(centroid) || numel(centroid) ~= 2 || any(~isfinite(centroid))
                        continue
                    end

                    labelTextHandles(end+1) = text(maskAx, ... %#ok<AGROW>
                        centroid(1), centroid(2), sprintf('%d', iLabel), ...
                        'Color', 'w', ...
                        'FontWeight', 'bold', ...
                        'HorizontalAlignment', 'center', ...
                        'HitTest', 'off', ...
                        'PickableParts', 'none');
                end

                hold(maskAx, 'off');
            end

            function nHoleRegions = iCountHoleRegions(labelMatrix)
                nHoleRegions = 0;
                nLabels = max(labelMatrix(:));

                for iLabel = 1:nLabels
                    componentMask = labelMatrix == iLabel;

                    try
                        if app.roiMaskHasHoles(componentMask)
                            nHoleRegions = nHoleRegions + 1;
                        end
                    catch
                    end
                end
            end

            function [imgOut, tf, message] = iSanitizeImage(imgIn, maskIn)
                tf = false;
                message = '';

                if ~ismatrix(imgIn) || ~isnumeric(imgIn)
                    imgOut = [];
                    message = 'Threshold ROI creation requires a 2D numeric image.';
                    return
                end

                imgOut = single(imgIn);

                finiteInsideMask = isfinite(imgOut) & logical(maskIn);

                if ~any(finiteInsideMask(:))
                    finiteAnywhere = isfinite(imgOut);
                    if ~any(finiteAnywhere(:))
                        message = 'Threshold ROI creation failed: current image contains no finite pixels.';
                        return
                    end

                    finiteInsideMask = finiteAnywhere;
                end

                finiteValues = imgOut(finiteInsideMask);
                minFinite = min(finiteValues, [], 'all');
                maxFinite = max(finiteValues, [], 'all');

                imgOut(isnan(imgOut)) = 0;
                imgOut(isinf(imgOut) & imgOut < 0) = minFinite;
                imgOut(isinf(imgOut) & imgOut > 0) = maxFinite;

                tf = true;
            end

        end

        function openAllenAtlasROICreator(app)
            %OPENALLENATLASROICREATOR Open child app for Allen atlas region selection.

            if ~app.hasData()
                app.setStatusMessage('Load image data before creating Allen atlas ROIs.');
                return
            end

            try
                app.AllenAtlasROICreatorApp = AllenAtlasROICreator(app);
                % waitfor
            catch ME
                app.setStatusMessage(sprintf('Could not open Allen Atlas ROI Creator: %s', ME.message));
                rethrow(ME)
            end

        end

        function [state, ok, msg] = buildAllenAtlasROIPlacementState(app, atlasResult, previousKeyFcn)
            %BUILDALLENATLASROIPLACEMENTSTATE Build transformed initial atlas placement state.
            %
            %   The state stores two geometry resolutions:
            %       - regionVertices0 / boundaryVertices0:
            %           full-resolution vertices used for final ROI creation.
            %
            %       - regionPreviewVertices0 / boundaryPreviewVertices0:
            %           reduced vertices used only for live fitting preview.
            %
            %   This keeps final ROI geometry accurate while making the bounding-box
            %   fitting interaction closer in speed to grouped ROI editing.

            ok = false;
            msg = '';
            state = struct();

            if ~isstruct(atlasResult) || ~isfield(atlasResult, 'selectedRegions')
                msg = 'Invalid Allen atlas selection result.';
                return
            end

            selectedRegions = atlasResult.selectedRegions;

            if isempty(selectedRegions)
                msg = 'No Allen atlas regions were selected.';
                return
            end

            requiredFields = {'BregmaCoordsXY', 'LambdaCoordsXY', 'pixels_per_mm', 'atlasImageSizeYX'};
            for iField = 1:numel(requiredFields)
                if ~isfield(atlasResult, requiredFields{iField})
                    msg = sprintf('Allen atlas result is missing "%s".', requiredFields{iField});
                    return
                end
            end

            [scaleXY, originXY, bHasCalibration, calibrationMsg] = ...
                app.getAllenAtlasInitialPlacementScaleAndOrigin(atlasResult);

            regionVertices0 = {};
            regionPreviewVertices0 = {};
            regionTag = {};
            regionType = {};
            regionName = {};
            regionColor = {};

            allVertices = zeros(0, 2);

            maxPreviewRegionVertices = 600;
            maxPreviewBoundaryVertices = 1500;

            for iRegion = 1:numel(selectedRegions)
                if ~isfield(selectedRegions(iRegion), 'verticesXY')
                    continue
                end

                verticesAtlas = double(selectedRegions(iRegion).verticesXY);

                if isempty(verticesAtlas) || size(verticesAtlas, 2) ~= 2
                    continue
                end

                vertices0 = app.transformAllenAtlasXYToViewerInitial( ...
                    verticesAtlas, ...
                    atlasResult.BregmaCoordsXY, ...
                    scaleXY, ...
                    originXY);

                finiteRows = all(isfinite(vertices0), 2);
                if nnz(finiteRows) < 3
                    continue
                end

                regionVertices0{end+1} = vertices0; %#ok<AGROW>
                regionPreviewVertices0{end+1} = iReducePreviewPolyline(vertices0, maxPreviewRegionVertices); %#ok<AGROW>

                if isfield(selectedRegions(iRegion), 'tag') && ~isempty(selectedRegions(iRegion).tag)
                    regionTag{end+1} = char(string(selectedRegions(iRegion).tag)); %#ok<AGROW>
                else
                    regionTag{end+1} = sprintf('AtlasROI_%d', numel(regionTag) + 1); %#ok<AGROW>
                end

                if isfield(selectedRegions(iRegion), 'type') && ~isempty(selectedRegions(iRegion).type)
                    regionType{end+1} = char(string(selectedRegions(iRegion).type)); %#ok<AGROW>
                else
                    regionType{end+1} = 'area'; %#ok<AGROW>
                end

                if isfield(selectedRegions(iRegion), 'name') && ~isempty(selectedRegions(iRegion).name)
                    regionName{end+1} = char(string(selectedRegions(iRegion).name)); %#ok<AGROW>
                else
                    regionName{end+1} = regionTag{end}; %#ok<AGROW>
                end

                if isfield(selectedRegions(iRegion), 'color') && numel(selectedRegions(iRegion).color) == 3
                    c = double(selectedRegions(iRegion).color(:).');
                else
                    c = [0.5 0.5 0.5];
                end

                regionColor{end+1} = min(max(c, 0), 1); %#ok<AGROW>
                allVertices = [allVertices; vertices0(finiteRows, :)]; %#ok<AGROW>
            end

            if isempty(regionVertices0) || isempty(allVertices)
                msg = 'Selected Allen atlas regions do not contain valid vertices.';
                return
            end

            boundaryVertices0 = {};
            boundaryPreviewVertices0 = {};

            if isfield(atlasResult, 'boundaryVerticesXY') && ~isempty(atlasResult.boundaryVerticesXY)
                for iB = 1:numel(atlasResult.boundaryVerticesXY)
                    xyAtlas = double(atlasResult.boundaryVerticesXY{iB});

                    if isempty(xyAtlas) || size(xyAtlas, 2) ~= 2
                        continue
                    end

                    xy0 = app.transformAllenAtlasXYToViewerInitial( ...
                        xyAtlas, ...
                        atlasResult.BregmaCoordsXY, ...
                        scaleXY, ...
                        originXY);

                    finiteRows = all(isfinite(xy0), 2);
                    if nnz(finiteRows) < 2
                        continue
                    end

                    boundaryVertices0{end+1} = xy0; %#ok<AGROW>
                    boundaryPreviewVertices0{end+1} = iReducePreviewPolyline(xy0, maxPreviewBoundaryVertices); %#ok<AGROW>
                end
            end

            bregma0 = app.transformAllenAtlasXYToViewerInitial( ...
                double(atlasResult.BregmaCoordsXY(:).'), ...
                atlasResult.BregmaCoordsXY, ...
                scaleXY, ...
                originXY);

            lambda0 = app.transformAllenAtlasXYToViewerInitial( ...
                double(atlasResult.LambdaCoordsXY(:).'), ...
                atlasResult.BregmaCoordsXY, ...
                scaleXY, ...
                originXY);

            minXY = min(allVertices, [], 1);
            maxXY = max(allVertices, [], 1);
            bbox = [minXY, maxXY - minXY];

            if any(~isfinite(bbox)) || bbox(3) <= 0 || bbox(4) <= 0
                msg = 'Allen atlas placement failed: invalid bounding box.';
                return
            end

            state.atlasResult = atlasResult;
            state.previousKeyFcn = previousKeyFcn;
            state.previousMode = app.InteractionMode;

            state.bHasCalibration = bHasCalibration;
            state.calibrationMsg = calibrationMsg;

            state.scaleXY = scaleXY;
            state.originXY = originXY;

            % Full-resolution geometry. Use this for final ROI creation.
            state.regionVertices0 = regionVertices0;
            state.boundaryVertices0 = boundaryVertices0;

            % Preview-resolution geometry. Use this for live fitting graphics.
            state.regionPreviewVertices0 = regionPreviewVertices0;
            state.boundaryPreviewVertices0 = boundaryPreviewVertices0;

            state.regionTag = regionTag;
            state.regionType = regionType;
            state.regionName = regionName;
            state.regionColor = regionColor;

            state.bregma0 = bregma0;
            state.lambda0 = lambda0;

            state.originalBBox = bbox;
            state.originalCenter = [bbox(1) + bbox(3) / 2, bbox(2) + bbox(4) / 2];
            state.originalSize = bbox(3:4);

            state.listeners = {};

            ok = true;

            if ~bHasCalibration
                msg = calibrationMsg;
            end

            % ---------------------------------------------------------------------
            % Local helper
            % ---------------------------------------------------------------------
            function xyOut = iReducePreviewPolyline(xyIn, maxVertices)
                %IREDUCEPREVIEWPOLYLINE Downsample visual-only atlas preview vertices.
                %
                %   NaN-separated polygon parts are reduced independently. This is used
                %   only for preview graphics, never for final ROI creation.

                xyIn = double(xyIn);

                if nargin < 2 || isempty(maxVertices)
                    maxVertices = 600;
                end

                if isempty(xyIn) || size(xyIn, 2) ~= 2 || size(xyIn, 1) <= maxVertices
                    xyOut = xyIn;
                    return
                end

                nanRows = any(~isfinite(xyIn), 2);

                if ~any(nanRows)
                    xyOut = iReduceOneSegment(xyIn, maxVertices);
                    return
                end

                xyOut = zeros(0, 2);
                breakIdx = [0; find(nanRows); size(xyIn, 1) + 1];

                for iSeg = 1:numel(breakIdx)-1
                    i1 = breakIdx(iSeg) + 1;
                    i2 = breakIdx(iSeg+1) - 1;

                    if i2 < i1
                        continue
                    end

                    seg = xyIn(i1:i2, :);
                    seg = seg(all(isfinite(seg), 2), :);

                    if isempty(seg)
                        continue
                    end

                    seg = iReduceOneSegment(seg, maxVertices);

                    if ~isempty(xyOut)
                        xyOut(end+1, :) = [NaN NaN]; %#ok<AGROW>
                    end

                    xyOut = [xyOut; seg]; %#ok<AGROW>
                end
            end

            function segOut = iReduceOneSegment(segIn, maxN)
                %IREDUCEONESEGMENT Uniformly downsample one finite vertex segment.

                if isempty(segIn) || size(segIn, 1) <= maxN
                    segOut = segIn;
                    return
                end

                step = max(1, ceil(size(segIn, 1) / maxN));
                keepIdx = 1:step:size(segIn, 1);

                if keepIdx(end) ~= size(segIn, 1)
                    keepIdx(end+1) = size(segIn, 1); %#ok<AGROW>
                end

                segOut = segIn(keepIdx, :);
            end

        end

        function [scaleXY, originXY, bHasCalibration, msg] = getAllenAtlasInitialPlacementScaleAndOrigin(app, atlasResult)
            %GETALLENATLASINITIALPLACEMENTSCALEANDORIGIN Return initial atlas->image transform.
            %
            %   If DataParams.view.pixelSize_px_per_mm and origin_xy_px are available,
            %   Bregma is aligned to the DataViewer origin. Otherwise, the atlas is placed
            %   near image center and can be manually fit with the bounding box.

            bHasCalibration = false;
            msg = '';

            sz = app.getDataSize();
            Ny = double(sz(1));
            Nx = double(sz(2));

            atlasSizeYX = double(atlasResult.atlasImageSizeYX(:).');
            if numel(atlasSizeYX) ~= 2 || any(~isfinite(atlasSizeYX)) || any(atlasSizeYX <= 0)
                atlasSizeYX = [Ny Nx];
            end

            originXY = [Nx, Ny] ./ 2;
            atlasScaleFallback = 0.85 .* min([Nx ./ atlasSizeYX(2), Ny ./ atlasSizeYX(1)]);
            scaleXY = [atlasScaleFallback, atlasScaleFallback];

            if isempty(app.DataParams) || ~isfield(app.DataParams, 'view')
                msg = 'View calibration not found. Atlas was placed near image center for manual fitting.';
                return
            end

            viewInfo = app.DataParams.view;

            if ~isfield(viewInfo, 'pixelSize_px_per_mm') || isempty(viewInfo.pixelSize_px_per_mm) || ...
                    ~isfield(viewInfo, 'origin_xy_px') || isempty(viewInfo.origin_xy_px)
                msg = 'Pixel size or origin is missing. Atlas was placed near image center for manual fitting.';
                return
            end

            pxPerMm = double(viewInfo.pixelSize_px_per_mm(:).');

            if isscalar(pxPerMm)
                pxPerMm = [pxPerMm, pxPerMm];
            end

            originCandidate = double(viewInfo.origin_xy_px(:).');

            if numel(pxPerMm) ~= 2 || any(~isfinite(pxPerMm)) || any(pxPerMm <= 0) || ...
                    numel(originCandidate) ~= 2 || any(~isfinite(originCandidate))
                msg = 'Invalid pixel size or origin. Atlas was placed near image center for manual fitting.';
                return
            end

            atlasPxPerMm = double(atlasResult.pixels_per_mm);

            if ~isscalar(atlasPxPerMm) || ~isfinite(atlasPxPerMm) || atlasPxPerMm <= 0
                msg = 'Invalid atlas pixels_per_mm. Atlas was placed near image center for manual fitting.';
                return
            end

            scaleXY = pxPerMm ./ atlasPxPerMm;
            originXY = originCandidate;

            if isequal(round(originXY), [1 1])
                bHasCalibration = true;
                msg = 'Origin is [1 1]. Atlas was scaled anatomically but may require manual positioning.';
                return
            end

            bHasCalibration = true;

        end

        function xyViewer = transformAllenAtlasXYToViewerInitial(app, xyAtlas, bregmaXY, scaleXY, originXY) %#ok<INUSL>
            %TRANSFORMALLENATLASXYTOVIEWERINITIAL Apply initial atlas->viewer transform.

            xyAtlas = double(xyAtlas);
            bregmaXY = double(bregmaXY(:).');
            scaleXY = double(scaleXY(:).');
            originXY = double(originXY(:).');

            xyViewer = (xyAtlas - bregmaXY) .* scaleXY + originXY;

        end

        function onAllenAtlasROIPlacementBoxMoving(app, src, evt) %#ok<INUSD>
            %ONALLENATLASROIPLACEMENTBOXMOVING Update atlas preview during box movement.

            app.updateAllenAtlasROIPlacementBoxSizeLabel(src);
            app.updateAllenAtlasROIPlacementPreview();

        end

        function onAllenAtlasROIPlacementBoxMoved(app, src, evt) %#ok<INUSD>
            %ONALLENATLASROIPLACEMENTBOXMOVED Update atlas preview after box movement.

            app.updateAllenAtlasROIPlacementBoxSizeLabel(src);
            app.updateAllenAtlasROIPlacementPreview();

        end

        function updateAllenAtlasROIPlacementBoxSizeLabel(app, hBox)
            %UPDATEALLENATLASROIPLACEMENTBOXSIZELABEL Show current atlas fit-box size.
            %
            %   Uses the same unit convention as grouped ROI editing.

            if ~app.isUsableGraphicsHandle(hBox)
                return
            end

            try
                pos = double(hBox.Position(:).');
            catch
                return
            end

            if numel(pos) < 4 || any(~isfinite(pos(3:4)))
                return
            end

            [scaleFactor, unitText] = app.getROIAxisUnitScale();
            sizeXY = abs(pos(3:4)) .* scaleFactor;

            labelText = sprintf('Atlas | %.3g x %.3g %s', ...
                sizeXY(1), sizeXY(2), unitText);

            app.setROIObjectPropertyIfAvailable(hBox, 'Label', labelText);
            app.setROIObjectPropertyIfAvailable(hBox, 'LabelVisible', 'on');

        end

        function onAllenAtlasROIPlacementBoxClicked(app, src, evt)
            %ONALLENATLASROIPLACEMENTBOXCLICKED Commit on double-click.

            try
                if isprop(evt, 'SelectionType') && strcmpi(evt.SelectionType, 'double')
                    app.commitAllenAtlasROIPlacement();
                end
            catch
            end

        end

        function onAllenAtlasROIPlacementBoxDeleting(app, src, evt) %#ok<INUSD>
            %ONALLENATLASROIPLACEMENTBOXDELETING Cancel if user deletes the placement box.

            if strcmpi(app.InteractionMode, 'editingROI')
                app.cancelAllenAtlasROIPlacement(false);
            end

        end

        function onAllenAtlasROIPlacementKeyPress(app, src, evt) %#ok<INUSL>
            %ONALLENATLASROIPLACEMENTKEYPRESS Keyboard commit/cancel for atlas placement.

            switch lower(evt.Key)
                case {'return', 'enter'}
                    app.commitAllenAtlasROIPlacement();

                case 'escape'
                    app.cancelAllenAtlasROIPlacement(true);
            end

        end

        function updateAllenAtlasROIPlacementPreview(app)
            %UPDATEALLENATLASROIPLACEMENTPREVIEW Transform temporary atlas preview graphics.

            if isempty(fieldnames(app.AtlasROIPlacementState)) || ...
                    ~app.isUsableGraphicsHandle(app.AtlasROIPlacementBoxHandle)
                return
            end

            state = app.AtlasROIPlacementState;

            % Selected atlas regions: preview-resolution vertices.
            for iRegion = 1:numel(state.regionPreviewVertices0)
                if iRegion > numel(app.AtlasROIPlacementRegionHandles) || ...
                        ~app.isUsableGraphicsHandle(app.AtlasROIPlacementRegionHandles(iRegion))
                    continue
                end

                xy = app.transformAllenAtlasPlacementVertices(state.regionPreviewVertices0{iRegion});

                app.AtlasROIPlacementRegionHandles(iRegion).XData = xy(:, 1);
                app.AtlasROIPlacementRegionHandles(iRegion).YData = xy(:, 2);
            end

            % Atlas border/bound_mask: preview-resolution vertices.
            if isfield(state, 'boundaryPreviewVertices0')
                boundaryVertices = state.boundaryPreviewVertices0;
            else
                boundaryVertices = state.boundaryVertices0;
            end

            for iB = 1:numel(boundaryVertices)
                if iB > numel(app.AtlasROIPlacementBoundaryHandles) || ...
                        ~app.isUsableGraphicsHandle(app.AtlasROIPlacementBoundaryHandles(iB))
                    continue
                end

                xy = app.transformAllenAtlasPlacementVertices(boundaryVertices{iB});

                app.AtlasROIPlacementBoundaryHandles(iB).XData = xy(:, 1);
                app.AtlasROIPlacementBoundaryHandles(iB).YData = xy(:, 2);
            end

            if app.isUsableGraphicsHandle(app.AtlasROIPlacementBregmaHandle)
                xy = app.transformAllenAtlasPlacementVertices(state.bregma0);
                app.AtlasROIPlacementBregmaHandle.XData = xy(1);
                app.AtlasROIPlacementBregmaHandle.YData = xy(2);
            end

            if app.isUsableGraphicsHandle(app.AtlasROIPlacementLambdaHandle)
                xy = app.transformAllenAtlasPlacementVertices(state.lambda0);
                app.AtlasROIPlacementLambdaHandle.XData = xy(1);
                app.AtlasROIPlacementLambdaHandle.YData = xy(2);
            end

            drawnow limitrate

        end

        function xyOut = reduceAllenAtlasPreviewPolyline(app, xyIn, maxVertices) %#ok<INUSL>
            %REDUCEALLENATLASPREVIEWPOLYLINE Downsample visual-only atlas preview lines.
            %
            %   This does not alter final ROI creation. Full-resolution vertices are kept in
            %   state.regionVertices0 and used by commitAllenAtlasROIPlacement.

            xyIn = double(xyIn);

            if nargin < 3 || isempty(maxVertices)
                maxVertices = 600;
            end

            if isempty(xyIn) || size(xyIn, 2) ~= 2 || size(xyIn, 1) <= maxVertices
                xyOut = xyIn;
                return
            end

            nanRows = any(~isfinite(xyIn), 2);

            if ~any(nanRows)
                xyOut = iReduceOneSegment(xyIn, maxVertices);
                return
            end

            xyOut = zeros(0, 2);
            breakIdx = [0; find(nanRows); size(xyIn, 1) + 1];

            for iSeg = 1:numel(breakIdx)-1
                i1 = breakIdx(iSeg) + 1;
                i2 = breakIdx(iSeg+1) - 1;

                if i2 < i1
                    continue
                end

                seg = xyIn(i1:i2, :);
                seg = iReduceOneSegment(seg, maxVertices);

                if ~isempty(xyOut)
                    xyOut(end+1, :) = [NaN NaN]; %#ok<AGROW>
                end

                xyOut = [xyOut; seg]; %#ok<AGROW>
            end

            function segOut = iReduceOneSegment(segIn, maxN)
                if size(segIn, 1) <= maxN
                    segOut = segIn;
                    return
                end

                step = max(1, ceil(size(segIn, 1) ./ maxN));
                keepIdx = 1:step:size(segIn, 1);

                if keepIdx(end) ~= size(segIn, 1)
                    keepIdx(end+1) = size(segIn, 1); %#ok<AGROW>
                end

                segOut = segIn(keepIdx, :);
            end

        end

        function verticesOut = transformAllenAtlasPlacementVertices(app, verticesIn)
            %TRANSFORMALLENATLASPLACEMENTVERTICES Apply current atlas placement transform.
            %
            %   Mirrors grouped ROI edit transform so atlas placement behaves like grouped
            %   ROI resizing.

            if isempty(verticesIn)
                verticesOut = verticesIn;
                return
            end

            transformInfo = app.getCurrentAllenAtlasPlacementTransform();
            verticesOut = app.transformVerticesFromGroupEdit(verticesIn, transformInfo);

        end

        function transformInfo = getCurrentAllenAtlasPlacementTransform(app)
            %GETCURRENTALLENATLASPLACEMENTTRANSFORM Return current atlas box transform.
            %
            %   Matches getCurrentGroupEditTransform.

            state = app.AtlasROIPlacementState;
            hBox = app.AtlasROIPlacementBoxHandle;

            pos = double(hBox.Position(:).');

            newSize = max(abs(pos(3:4)), eps);
            newCenter = [pos(1) + pos(3)/2, pos(2) + pos(4)/2];

            angleDeg = app.getNumericROIObjectProperty(hBox, 'RotationAngle', 0);

            transformInfo = struct();
            transformInfo.originalCenter = state.originalCenter;
            transformInfo.newCenter = newCenter;
            transformInfo.scaleXY = newSize ./ max(state.originalSize, eps);
            transformInfo.angleDeg = angleDeg;

        end

        function stackAllenAtlasROIPlacementGraphics(app)
            %STACKALLENATLASROIPLACEMENTGRAPHICS Keep placement controls/markers visible.

            hList = [ ...
                app.AtlasROIPlacementRegionHandles(:); ...
                app.AtlasROIPlacementBoundaryHandles(:); ...
                app.AtlasROIPlacementBregmaHandle(:); ...
                app.AtlasROIPlacementLambdaHandle(:); ...
                app.AtlasROIPlacementBoxHandle(:)];

            for iH = 1:numel(hList)
                try
                    if app.isUsableGraphicsHandle(hList(iH))
                        uistack(hList(iH), 'top');
                    end
                catch
                end
            end

        end

        function createAllenAtlasROIPlacementContextMenu(app)
            %CREATEALLENATLASROIPLACEMENTCONTEXTMENU Create placement box context menu.

            if ~app.isUsableGraphicsHandle(app.AtlasROIPlacementBoxHandle)
                return
            end

            try
                if ~isempty(app.AtlasROIPlacementContextMenu) && isvalid(app.AtlasROIPlacementContextMenu)
                    delete(app.AtlasROIPlacementContextMenu);
                end
            catch
            end

            app.AtlasROIPlacementContextMenu = uicontextmenu(app.UIFigure);

            app.AtlasROIPlacementConfirmMenu = uimenu(app.AtlasROIPlacementContextMenu);
            app.AtlasROIPlacementConfirmMenu.Text = 'Create atlas ROIs';
            app.AtlasROIPlacementConfirmMenu.MenuSelectedFcn = @(src, evt) app.commitAllenAtlasROIPlacement();

            app.AtlasROIPlacementCancelMenu = uimenu(app.AtlasROIPlacementContextMenu);
            app.AtlasROIPlacementCancelMenu.Text = 'Cancel';
            app.AtlasROIPlacementCancelMenu.Separator = 'on';
            app.AtlasROIPlacementCancelMenu.MenuSelectedFcn = @(src, evt) app.cancelAllenAtlasROIPlacement(true);

            app.AtlasROIPlacementBoxHandle.ContextMenu = app.AtlasROIPlacementContextMenu;

        end

        function commitAllenAtlasROIPlacement(app)
            %COMMITALLENATLASROIPLACEMENT Create final polygon ROIs from atlas placement.

            if isempty(fieldnames(app.AtlasROIPlacementState))
                return
            end
            polyshapeWarningCleanup = app.suppressExpectedPolyshapeWarnings(); %#ok<NASGU>
            state = app.AtlasROIPlacementState;

            [usedIDs, usedNames] = app.getExistingROIIdentifiersForLoad('append');

            template = app.makeLoadedROITemplate();
            template.type = 'polygon';
            template.DOC = datetime('now');
            template.modifiedOn = template.DOC;
            template.notes = 'Created from Allen Brain Atlas';

            newROIList = repmat(template, 1, 0);

            nSourceSkipped = 0;
            nHoleSkipped = 0;
            nSplit = 0;

            for iRegion = 1:numel(state.regionVertices0)
                verticesXY = app.transformAllenAtlasPlacementVertices(state.regionVertices0{iRegion});
                regionMask = app.createMaskFromAtlasVertices(verticesXY);

                if isempty(regionMask) || ~any(regionMask(:))
                    nSourceSkipped = nSourceSkipped + 1;
                    continue
                end

                regionMask = app.clipROIMaskToActiveLogicalMask(regionMask);

                if isempty(regionMask) || ~any(regionMask(:))
                    nSourceSkipped = nSourceSkipped + 1;
                    continue
                end

                if app.roiMaskHasHoles(regionMask)
                    nHoleSkipped = nHoleSkipped + 1;
                    continue
                end

                componentMasks = app.maskToConnectedComponentMasks(regionMask);

                if isempty(componentMasks)
                    nSourceSkipped = nSourceSkipped + 1;
                    continue
                end

                if numel(componentMasks) > 1
                    nSplit = nSplit + 1;
                end

                baseName = char(string(state.regionTag{iRegion}));

                template.color = state.regionColor{iRegion};
                template.notes = sprintf('Created from Allen Brain Atlas | type=%s | name=%s', ...
                    state.regionType{iRegion}, state.regionName{iRegion});

                for iComp = 1:numel(componentMasks)
                    if numel(componentMasks) > 1
                        requestedName = sprintf('%s_part%d', baseName, iComp);
                    else
                        requestedName = baseName;
                    end

                    roiName = app.makeUniqueNameAgainstList(requestedName, usedNames);
                    roiID = app.getNextUniqueROIIDFromList(usedIDs);
                    try
                        newROI = app.makePolygonROIFromMaskComponent( ...
                            template, ...
                            componentMasks{iComp}, ...
                            roiName, ...
                            roiID, ...
                            false);
                    catch
                        nSourceSkipped = nSourceSkipped + 1;
                        continue
                    end

                    newROI.runtime = app.makeROIRuntimeStruct(true);
                    newROI.runtime.selected = false;
                    newROI.runtime.ROIHandle = gobjects(1);

                    newROIList(end+1) = newROI; %#ok<AGROW>
                    usedIDs(end+1) = roiID; %#ok<AGROW>
                    usedNames{end+1} = roiName; %#ok<AGROW>
                end
            end

            if isempty(newROIList)
                app.setStatusMessage('Allen atlas ROI creation failed: no valid ROIs after fitting/clipping.');
                return
            end

            app.cancelAllenAtlasROIPlacement(false);

            app.ROIList = app.appendROIsToROIList(app.ROIList, newROIList);
            app.rebuildMissingROIOverlaysAfterLoad();

            app.SelectedROIID = NaN;
            app.ROISelectionOrder = [];

            app.refreshROITable();
            app.refreshROITraces();
            app.refreshEventPatches();
            app.stackROITraceGraphics();
            app.updateGUIEnabledState();

            msg = sprintf('Created %d Allen atlas ROI(s).', numel(newROIList));

            if nSplit > 0
                msg = sprintf('%s %d source region(s) split into parts.', msg, nSplit);
            end

            if nHoleSkipped > 0
                msg = sprintf('%s %d hole-containing source region(s) skipped.', msg, nHoleSkipped);
            end

            if nSourceSkipped > 0
                msg = sprintf('%s %d invalid/outside source region(s) skipped.', msg, nSourceSkipped);
            end

            app.setStatusMessage(msg);

        end

        function mask = createMaskFromAtlasVertices(app, verticesXY)
            %CREATEMASKFROMATLASVERTICES Rasterize atlas polygon vertices into image mask.
            %
            %   Handles NaN-separated polyshape boundaries by rasterizing each finite
            %   contour and taking the union.

            mask = [];

            if isempty(verticesXY) || size(verticesXY, 2) ~= 2
                return
            end

            verticesXY = double(verticesXY);
            mask = false(app.getDataSize2D());

            nanRows = any(~isfinite(verticesXY), 2);
            breakIdx = [0; find(nanRows); size(verticesXY, 1) + 1];

            for iSeg = 1:numel(breakIdx)-1
                idx1 = breakIdx(iSeg) + 1;
                idx2 = breakIdx(iSeg+1) - 1;

                if idx2 < idx1
                    continue
                end

                contourXY = verticesXY(idx1:idx2, :);
                contourXY = contourXY(all(isfinite(contourXY), 2), :);

                if size(contourXY, 1) < 3
                    continue
                end

                try
                    mask = mask | app.createMaskFromVertices(contourXY);
                catch
                end
            end

        end

        function sz2 = getDataSize2D(app)
            %GETDATASIZE2D Return current image [Y X] size.

            sz = app.getDataSize();
            sz2 = sz(1:2);

        end

        function cancelAllenAtlasROIPlacement(app, bShowStatus)
            %CANCELALLENATLASROIPLACEMENT Delete temporary atlas placement graphics.

            if nargin < 2
                bShowStatus = true;
            end

            previousKeyFcn = [];

            try
                if isfield(app.AtlasROIPlacementState, 'previousKeyFcn')
                    previousKeyFcn = app.AtlasROIPlacementState.previousKeyFcn;
                end
            catch
            end

            % Delete listeners first.
            try
                if isfield(app.AtlasROIPlacementState, 'listeners')
                    listeners = app.AtlasROIPlacementState.listeners;

                    for iL = 1:numel(listeners)
                        try
                            if isvalid(listeners{iL})
                                delete(listeners{iL});
                            end
                        catch
                        end
                    end
                end
            catch
            end

            hList = [ ...
                app.AtlasROIPlacementRegionHandles(:); ...
                app.AtlasROIPlacementBoundaryHandles(:); ...
                app.AtlasROIPlacementBregmaHandle(:); ...
                app.AtlasROIPlacementLambdaHandle(:); ...
                app.AtlasROIPlacementBoxHandle(:)];

            for iH = 1:numel(hList)
                try
                    if app.isUsableGraphicsHandle(hList(iH))
                        delete(hList(iH));
                    end
                catch
                end
            end

            try
                if ~isempty(app.AtlasROIPlacementContextMenu) && isvalid(app.AtlasROIPlacementContextMenu)
                    delete(app.AtlasROIPlacementContextMenu);
                end
            catch
            end

            app.AtlasROIPlacementState = struct();
            app.AtlasROIPlacementRegionHandles = gobjects(0);
            app.AtlasROIPlacementBoundaryHandles = gobjects(0);
            app.AtlasROIPlacementBregmaHandle = gobjects(1);
            app.AtlasROIPlacementLambdaHandle = gobjects(1);
            app.AtlasROIPlacementBoxHandle = gobjects(1);
            app.AtlasROIPlacementContextMenu = matlab.ui.container.ContextMenu.empty;
            app.AtlasROIPlacementConfirmMenu = matlab.ui.container.Menu.empty;
            app.AtlasROIPlacementCancelMenu = matlab.ui.container.Menu.empty;

            try
                app.UIFigure.WindowKeyPressFcn = previousKeyFcn;
            catch
            end

            app.setInteractionMode('idle');
            app.updateGUIEnabledState();

            if bShowStatus
                app.setStatusMessage('Allen atlas ROI placement cancelled.');
            end

        end

        function [rawFolder,saveFolder] = selectRawAndSaveFolderForDataImport(app,fileExt)

            rawFolder = [];
            saveFolder = [];

            if ~isempty(app.CurrentFile)
                initialPath = fileparts(app.CurrentFile);
            else
                initialPath = pwd;
            end
            rawFolder = uigetdir(initialPath, ['Select the raw folder containing the ' fileExt ' file(s)']);
            if ~rawFolder
                rawFolder = [];
                return
            end

            % SaveFolder
            saveFolder = uigetdir(rawFolder, 'Select the folder to save the .dat file(s)');
            if ~saveFolder
                saveFolder = [];
                rawFolder = [];
                return
            end
        end

        function startSetImageReference(app)
            %STARTSETIMAGEREFERENCE Create a manually transformed reference image.
            %
            %   startSetImageReference(app) captures the current displayed frame,
            %   lets the user translate and rotate it through a full-frame bounding
            %   box, applies the temporary transform to the captured frame, previews
            %   the transformed frame live, and saves only the resulting 2D image as:
            %
            %       DataParams.registration.referenceImage
            %
            %   The temporary transform is intentionally not stored. This method does
            %   not modify DataParams.registration.tform and does not mark the data as
            %   registered.

            if ~app.hasData()
                app.setStatusMessage('Load image data before setting an image reference.');
                return
            end

            if ~strcmp(app.InteractionMode, 'idle')
                app.setStatusMessage('Finish the current interaction before setting an image reference.');
                return
            end

            if app.hasExistingImageReference()
                bOverwrite = app.confirmOverwriteImageReference();

                if ~bOverwrite
                    app.setStatusMessage('Set Reference cancelled. Existing image reference was kept.');
                    return
                end
            end

            sourceFrame = app.getCurrentFrame();
            sourceFrame = squeeze(sourceFrame);

            if ndims(sourceFrame) ~= 2 || isempty(sourceFrame)
                app.setStatusMessage('Current frame is not a valid 2D image.');
                return
            end

            sourceFrame = single(sourceFrame);
            originalDisplayedFrame = sourceFrame;

            [Ny, Nx] = size(sourceFrame);
            imageSizeYX = double([Ny Nx]);

            % Use pixel-edge coordinates so the manipulator spans the displayed image.
            originalBoxPosition = [0.5, 0.5, Nx, Ny];
            originalCenterXY = [ ...
                originalBoxPosition(1) + originalBoxPosition(3) / 2, ...
                originalBoxPosition(2) + originalBoxPosition(4) / 2];

            originalSizeXY = originalBoxPosition(3:4);

            previousMode = app.InteractionMode;
            previousKeyFcn = app.UIFigure.WindowKeyPressFcn;

            hBox = [];
            listeners = {};
            bDone = false;
            bSaved = false;
            bCleaningUp = false;
            bUpdatingBox = false;

            try
                app.setInteractionMode('settingViewCalibration');
                app.setStatusMessage(['Move/rotate the reference box. ' ...
                    'Double-click or press Enter to preview/save; press Escape to cancel.']);

                figure(app.UIFigure);
                hold(app.ImageAxes, 'on');

                hBox = drawrectangle(app.ImageAxes, ...
                    'Position', originalBoxPosition, ...
                    'Color', [0 1 1], ...
                    'FaceAlpha', 0.08, ...
                    'LineWidth', 1.5, ...
                    'DrawingArea', 'unlimited', ...
                    'InteractionsAllowed', 'all', ...
                    'Rotatable', true);

                updateBoxLabel(hBox);
                updateReferencePreview(hBox);

                listeners{end+1} = addlistener(hBox, 'MovingROI', ...
                    @(src, evt) onBoxMoving(src, evt)); %#ok<AGROW>
                listeners{end+1} = addlistener(hBox, 'ROIMoved', ...
                    @(src, evt) onBoxMoved(src, evt)); %#ok<AGROW>
                listeners{end+1} = addlistener(hBox, 'ROIClicked', ...
                    @(src, evt) onBoxClicked(src, evt)); %#ok<AGROW>

                try
                    listeners{end+1} = addlistener(hBox, 'DeletingROI', ...
                        @(src, evt) onBoxDeleting(src, evt)); %#ok<AGROW>
                catch
                    % Some MATLAB releases do not expose DeletingROI for all ROI objects.
                end

                app.UIFigure.WindowKeyPressFcn = @(src, evt) onKeyPress(src, evt);

                while ~bDone && isvalid(app.UIFigure) && app.isUsableGraphicsHandle(hBox)
                    drawnow
                    pause(0.02)
                end

                if bSaved
                    app.setStatusMessage('Image reference saved in DataParams.mat.');
                else
                    app.setStatusMessage('Set Reference cancelled.');
                end

            catch ME
                cleanupRuntime();
                rethrow(ME)
            end

            cleanupRuntime();

            % ---------------------------------------------------------------------
            % Nested callbacks and local utilities
            % ---------------------------------------------------------------------
            function onBoxMoving(src, ~)
                lockBoxSize(src);
                updateBoxLabel(src);
                updateReferencePreview(src);
            end

            function onBoxMoved(src, ~)
                lockBoxSize(src);
                updateBoxLabel(src);
                updateReferencePreview(src);
            end

            function onBoxClicked(~, evt)
                try
                    if isprop(evt, 'SelectionType') && strcmpi(evt.SelectionType, 'double')
                        attemptConfirmReference();
                    end
                catch ME
                    app.setStatusMessage(sprintf('Reference confirmation failed: %s', ME.message));
                end
            end

            function onBoxDeleting(~, ~)
                if bCleaningUp
                    return
                end

                bSaved = false;
                bDone = true;
            end

            function onKeyPress(~, evt)
                switch lower(evt.Key)
                    case {'return', 'enter'}
                        attemptConfirmReference();

                    case 'escape'
                        bSaved = false;
                        bDone = true;
                end
            end

            function lockBoxSize(src)
                %LOCKBOXSIZE Keep the manipulator full-frame sized.
                %
                %   This utility supports translation and rotation only. Resize
                %   attempts are ignored because reference creation should not scale
                %   the image.

                if bUpdatingBox || ~app.isUsableGraphicsHandle(src)
                    return
                end

                bUpdatingBox = true;
                cleanupSizeLock = onCleanup(@() setUpdatingBoxFalse());

                try
                    pos = double(src.Position(:).');

                    if numel(pos) < 4 || any(~isfinite(pos(1:4)))
                        return
                    end

                    centerXY = [pos(1) + pos(3) / 2, pos(2) + pos(4) / 2];

                    src.Position = [ ...
                        centerXY(1) - originalSizeXY(1) / 2, ...
                        centerXY(2) - originalSizeXY(2) / 2, ...
                        originalSizeXY(1), ...
                        originalSizeXY(2)];
                catch
                end

                function setUpdatingBoxFalse()
                    bUpdatingBox = false;
                end
            end

            function updateBoxLabel(src)
                if ~app.isUsableGraphicsHandle(src)
                    return
                end

                try
                    pos = double(src.Position(:).');
                    centerXY = [pos(1) + pos(3) / 2, pos(2) + pos(4) / 2];
                    translationXY = centerXY - originalCenterXY;
                    rotationDeg = getBoxRotationDeg(src);

                    labelText = sprintf( ...
                        'Reference | move [%.1f %.1f] px | rotate %.1f deg', ...
                        translationXY(1), translationXY(2), rotationDeg);

                    if isprop(src, 'Label')
                        src.Label = labelText;
                    end

                    if isprop(src, 'LabelVisible')
                        src.LabelVisible = 'on';
                    end
                catch
                end
            end

            function rotationDeg = getBoxRotationDeg(src)
                rotationDeg = 0;

                if ~app.isUsableGraphicsHandle(src)
                    return
                end

                try
                    if isprop(src, 'RotationAngle') && ~isempty(src.RotationAngle)
                        rotationDeg = double(src.RotationAngle);
                    end
                catch
                    rotationDeg = 0;
                end

                if ~isscalar(rotationDeg) || ~isfinite(rotationDeg)
                    rotationDeg = 0;
                end
            end

            function updateReferencePreview(src)
                %UPDATEREFERENCEPREVIEW Show the transformed frame during interaction.

                if ~app.isUsableGraphicsHandle(src) || ...
                        isempty(app.ImageHandle) || ~isvalid(app.ImageHandle)
                    return
                end

                try
                    previewFrame = makeReferenceImageFromBox(src);
                    previewFrame(~isfinite(previewFrame)) = 0;

                    app.ImageHandle.CData = previewFrame;
                    drawnow limitrate
                catch ME
                    app.setStatusMessage(sprintf('Reference preview failed: %s', ME.message));
                end
            end

            function referenceImage = makeReferenceImageFromBox(src)
                %MAKEREFERENCEIMAGEFROMBOX Apply the temporary transform to sourceFrame.

                pos = double(src.Position(:).');

                if numel(pos) < 4 || any(~isfinite(pos(1:4)))
                    error('DataViewer:InvalidReferenceBox', ...
                        'Reference box position is invalid.');
                end

                currentCenterXY = [pos(1) + pos(3) / 2, pos(2) + pos(4) / 2];
                rotationDeg = getBoxRotationDeg(src);

                % Match the image-coordinate rotation convention already used by
                % the app for ROI group transforms.
                angleRad = -deg2rad(rotationDeg);
                R = [cos(angleRad), -sin(angleRad); ...
                    sin(angleRad),  cos(angleRad)];

                Rrow = R.';
                translationXY = currentCenterXY - originalCenterXY * Rrow;

                T = [ ...
                    Rrow(1, 1),        Rrow(1, 2),        0; ...
                    Rrow(2, 1),        Rrow(2, 2),        0; ...
                    translationXY(1),  translationXY(2),  1];

                manualReferenceTransform = affine2d(T);
                outputRef = imref2d([Ny Nx]);

                referenceImage = imwarp( ...
                    sourceFrame, ...
                    manualReferenceTransform, ...
                    'OutputView', outputRef, ...
                    'FillValues', NaN);

                referenceImage = single(referenceImage);
            end

            function attemptConfirmReference()
                %ATTEMPTCONFIRMREFERENCE Preview transformed image before saving.

                if ~app.isUsableGraphicsHandle(hBox)
                    return
                end

                referenceImage = makeReferenceImageFromBox(hBox);
                referenceImage(~isfinite(referenceImage)) = 0;
                referenceImage = single(referenceImage);

                choice = app.confirmReferenceImageDialog(referenceImage);

                switch choice
                    case 'save'
                        saveReferenceImage(referenceImage);
                        bSaved = true;
                        bDone = true;

                    case 'resume'
                        updateReferencePreview(hBox);
                        app.setStatusMessage('Reference editing resumed.');

                    case 'cancel'
                        bSaved = false;
                        bDone = true;
                end
            end

            function saveReferenceImage(referenceImage)
                %SAVEREFERENCEIMAGE Save transformed reference image to DataParams.
                %
                %   The affine transform used to generate referenceImage is temporary
                %   and local to reference image generation. It is intentionally not
                %   saved to DataParams.

                referenceImage(~isfinite(referenceImage)) = 0;
                referenceImage = single(referenceImage);

                folderPath = app.getCurrentDataFolder();

                % updateDataParam requires intermediate fields to already exist, so
                % ensure the folder-level DataParams schema exists before writing.
                app.DataParams = app.ensureDataParamsFileForCurrentFolder();
                app.ensureDataParamsViewFields();
                app.ensureDataParamsMaskFields();

                sourceFile = char(string(app.CurrentFile));

                if isempty(app.CurrentEntry)
                    entryText = '';
                else
                    entryText = sprintf(' | entry: %s', app.CurrentEntry);
                end

                if strcmpi(app.ViewMode, 'event')
                    eventText = sprintf(' | condition: %s | repetition: %s', ...
                        app.CurrentCondition, app.CurrentRepetition);
                else
                    eventText = '';
                end

                referenceDescription = sprintf([ ...
                    'Reference image manually transformed from original data using ' ...
                    'DataViewer Set Reference utility. Source: %s%s | frame: %d | mode: %s%s.'], ...
                    sourceFile, entryText, app.CurrentFrame, app.ViewMode, eventText);

                % Grouped update. Validate only after the final related field is set.
                DataParams = updateDataParam(folderPath, ...
                    'view.imageSizeYX', imageSizeYX, ...
                    'validateAfterSet', false);

                DataParams = updateDataParam(folderPath, ...
                    'registration.referenceImage', referenceImage, ...
                    'validateAfterSet', false);

                DataParams = updateDataParam(folderPath, ...
                    'registration.referenceFile', sourceFile, ...
                    'validateAfterSet', false);

                DataParams = updateDataParam(folderPath, ...
                    'registration.referenceDescription', referenceDescription, ...
                    'validateAfterSet', true);

                app.DataParams = DataParams;
                app.ensureDataParamsViewFields();
                app.ensureDataParamsMaskFields();
                app.updateImageStatusLabel();
                app.refreshSetReferenceButtonContextMenuState();
            end

            function cleanupRuntime()
                bCleaningUp = true;

                for iListener = 1:numel(listeners)
                    try
                        if ~isempty(listeners{iListener}) && isvalid(listeners{iListener})
                            delete(listeners{iListener});
                        end
                    catch
                    end
                end

                listeners = {};

                try
                    if app.isUsableGraphicsHandle(hBox)
                        delete(hBox);
                    end
                catch
                end

                hBox = [];

                % Restore the viewer to the real source frame. The saved reference is
                % stored in DataParams, but the DataViewer display should continue to
                % represent the loaded data source.
                try
                    if ~isempty(app.ImageHandle) && isvalid(app.ImageHandle)
                        app.ImageHandle.CData = originalDisplayedFrame;
                    end
                catch
                end

                try
                    if isvalid(app.UIFigure)
                        app.UIFigure.WindowKeyPressFcn = previousKeyFcn;
                    end
                catch
                end

                try
                    app.setInteractionMode(previousMode);
                catch
                    app.InteractionMode = previousMode;
                    app.updateGUIEnabledState();
                end
            end

        end

        function choice = confirmReferenceImageDialog(app, referenceImage)
            %CONFIRMREFERENCEIMAGEDIALOG Preview transformed reference before saving.
            %
            %   choice = confirmReferenceImageDialog(app, referenceImage) opens a
            %   modal preview dialog and returns one of:
            %
            %       'save'    - Save reference image to DataParams.mat
            %       'resume'  - Return to live reference editing
            %       'cancel'  - Abort Set Reference utility
            %
            %   The preview uses the current image colormap and CLim for visual
            %   consistency. Display settings are not saved with the reference image.

            choice = 'resume';

            referenceImage = squeeze(referenceImage);

            if ndims(referenceImage) ~= 2 || isempty(referenceImage)
                error('DataViewer:InvalidReferencePreview', ...
                    'Reference image preview must be a non-empty 2D image.');
            end

            referenceImage = single(referenceImage);
            referenceImage(~isfinite(referenceImage)) = 0;

            dlg = uifigure( ...
                'Name', 'Confirm Image Reference', ...
                'WindowStyle', 'modal', ...
                'Position', [100 100 520 520], ...
                'Visible', 'off', ...
                'CloseRequestFcn', @onResumeEditing);

            grid = uigridlayout(dlg);
            grid.RowHeight = {28, '1x', 42, 38};
            grid.ColumnWidth = {'1x'};
            grid.Padding = [12 12 12 12];

            titleLabel = uilabel(grid);
            titleLabel.Text = 'Preview of image reference to save';
            titleLabel.FontWeight = 'bold';
            titleLabel.Layout.Row = 1;
            titleLabel.Layout.Column = 1;

            previewAxes = uiaxes(grid);
            previewAxes.Layout.Row = 2;
            previewAxes.Layout.Column = 1;

            imagesc(previewAxes, referenceImage);
            axis(previewAxes, 'image');
            axis(previewAxes, 'off');

            try
                colormap(previewAxes, app.ImageAxes.Colormap);
            catch
            end

            try
                previewAxes.CLim = app.ImageAxes.CLim;
            catch
            end

            infoLabel = uilabel(grid);
            infoLabel.Text = 'Empty pixels from rotation/translation will be saved as 0.';
            infoLabel.FontAngle = 'italic';
            infoLabel.Layout.Row = 3;
            infoLabel.Layout.Column = 1;

            buttonGrid = uigridlayout(grid);
            buttonGrid.RowHeight = {'1x'};
            buttonGrid.ColumnWidth = {'1x', '1x', '1x'};
            buttonGrid.Padding = [0 0 0 0];
            buttonGrid.ColumnSpacing = 8;
            buttonGrid.Layout.Row = 4;
            buttonGrid.Layout.Column = 1;

            saveButton = uibutton(buttonGrid, 'push');
            saveButton.Text = 'Save Reference';
            saveButton.Layout.Row = 1;
            saveButton.Layout.Column = 1;
            saveButton.ButtonPushedFcn = @onSaveReference;

            resumeButton = uibutton(buttonGrid, 'push');
            resumeButton.Text = 'Resume Editing';
            resumeButton.Layout.Row = 1;
            resumeButton.Layout.Column = 2;
            resumeButton.ButtonPushedFcn = @onResumeEditing;

            cancelButton = uibutton(buttonGrid, 'push');
            cancelButton.Text = 'Cancel';
            cancelButton.Layout.Row = 1;
            cancelButton.Layout.Column = 3;
            cancelButton.ButtonPushedFcn = @onCancel;

            if exist('placeAppInsideCaller', 'file') == 2
                placeAppInsideCaller(app, dlg, 'center');
            end

            dlg.Visible = 'on';
            uiwait(dlg);

            if isvalid(dlg)
                delete(dlg);
            end

            function onSaveReference(~, ~)
                choice = 'save';
                uiresume(dlg);
            end

            function onResumeEditing(~, ~)
                choice = 'resume';
                uiresume(dlg);
            end

            function onCancel(~, ~)
                choice = 'cancel';
                uiresume(dlg);
            end

        end

        function tf = hasExistingImageReference(app)
            %HASEXISTINGIMAGEREFERENCE Return true when DataParams has a reference image.
            %
            %   This method checks the folder-level DataParams.mat directly. It does
            %   not create, validate, or modify DataParams.mat.

            [tf, ~, ~] = app.readExistingImageReference();

        end

        function tf = confirmOverwriteImageReference(app)
            %CONFIRMOVERWRITEIMAGEREFERENCE Ask before replacing an existing reference image.
            %
            %   tf = confirmOverwriteImageReference(app) returns true only when the
            %   user explicitly chooses to overwrite the existing image reference.

            tf = false;

            msg = ['An image reference already exists in DataParams.mat. ' ...
                'Do you want to overwrite it?'];

            try
                selection = uiconfirm(app.UIFigure, ...
                    msg, ...
                    'Overwrite image reference?', ...
                    'Options', {'Overwrite', 'Cancel'}, ...
                    'DefaultOption', 'Cancel', ...
                    'CancelOption', 'Cancel', ...
                    'Icon', 'warning');

                tf = strcmp(selection, 'Overwrite');

            catch
                % Fallback for contexts where uiconfirm is unavailable.
                answer = questdlg( ...
                    msg, ...
                    'Overwrite image reference?', ...
                    'Overwrite', ...
                    'Cancel', ...
                    'Cancel');

                tf = strcmp(answer, 'Overwrite');
            end

        end

        function createSetReferenceButtonContextMenu(app)
            %CREATESETREFERENCEBUTTONCONTEXTMENU Create Set Reference button menu.
            %
            %   The context menu provides read-only utilities related to the saved
            %   image reference. Creating or refreshing the menu does not create or
            %   modify DataParams.mat.

            if isempty(app.ImageRefButton) || ~isvalid(app.ImageRefButton)
                return
            end

            if ~isempty(app.SetReferenceButtonContextMenu) && ...
                    isvalid(app.SetReferenceButtonContextMenu)
                try
                    app.ImageRefButton.ContextMenu = app.SetReferenceButtonContextMenu;
                catch
                end

                app.refreshSetReferenceButtonContextMenuState();
                return
            end

            app.SetReferenceButtonContextMenu = uicontextmenu(app.UIFigure);

            app.PreviewImageReferenceMenu = uimenu(app.SetReferenceButtonContextMenu);
            app.PreviewImageReferenceMenu.Text = 'Preview Existing Reference';
            app.PreviewImageReferenceMenu.MenuSelectedFcn = ...
                @(src, evt) app.previewExistingImageReference();

            try
                app.ImageRefButton.ContextMenu = app.SetReferenceButtonContextMenu;
            catch
            end

            app.refreshSetReferenceButtonContextMenuState();

        end

        function refreshSetReferenceButtonContextMenuState(app)
            %REFRESHSETREFERENCEBUTTONCONTEXTMENUSTATE Enable reference preview menu.

            if isempty(app.PreviewImageReferenceMenu) || ~isvalid(app.PreviewImageReferenceMenu)
                return
            end

            if app.hasExistingImageReference()
                app.PreviewImageReferenceMenu.Enable = 'on';
            else
                app.PreviewImageReferenceMenu.Enable = 'off';
            end

        end

        function [tf, referenceImage, DataParams] = readExistingImageReference(app)
            %READEXISTINGIMAGEREFERENCE Load saved image reference from DataParams.mat.
            %
            %   [tf, referenceImage, DataParams] = readExistingImageReference(app)
            %   returns true only when DataParams.registration.referenceImage exists
            %   and is non-empty.
            %
            %   This method is read-only. It does not create DataParams.mat and does
            %   not modify app.DataParams.

            tf = false;
            referenceImage = [];
            DataParams = struct();

            if isempty(app.CurrentFile)
                return
            end

            folderPath = fileparts(app.CurrentFile);
            dataParamsPath = fullfile(folderPath, 'DataParams.mat');

            if ~isfile(dataParamsPath)
                return
            end

            try
                S = load(dataParamsPath, 'DataParams');

                if ~isfield(S, 'DataParams') || ...
                        ~isstruct(S.DataParams) || ...
                        ~isscalar(S.DataParams)
                    return
                end

                DataParams = S.DataParams;

                if ~isfield(DataParams, 'registration') || ...
                        ~isstruct(DataParams.registration) || ...
                        ~isscalar(DataParams.registration) || ...
                        ~isfield(DataParams.registration, 'referenceImage') || ...
                        isempty(DataParams.registration.referenceImage)
                    return
                end

                referenceImage = DataParams.registration.referenceImage;

                if ~(isnumeric(referenceImage) || islogical(referenceImage)) || ...
                        ndims(referenceImage) ~= 2
                    referenceImage = [];
                    return
                end

                referenceImage = single(referenceImage);
                referenceImage(~isfinite(referenceImage)) = 0;

                tf = true;

            catch
                tf = false;
                referenceImage = [];
                DataParams = struct();
            end

        end

        function previewExistingImageReference(app)
            %PREVIEWEXISTINGIMAGEREFERENCE Show saved reference image and metadata.

            [bHasReference, referenceImage, DataParams] = app.readExistingImageReference();

            if ~bHasReference
                app.setStatusMessage('No image reference found in DataParams.mat.');

                try
                    uialert(app.UIFigure, ...
                        'No image reference was found in DataParams.mat.', ...
                        'No image reference', ...
                        'Icon', 'info');
                catch
                end
                return
            end

            metadataRows = app.buildImageReferenceMetadataRows(referenceImage, DataParams);
            app.showExistingImageReferencePreviewDialog(referenceImage, metadataRows);

        end

        function metadataRows = buildImageReferenceMetadataRows(app, referenceImage, DataParams)
            %BUILDIMAGEREFERENCEMETADATAROWS Build display metadata for reference preview.

            metadataRows = cell(0, 2);

            if nargin < 3 || isempty(DataParams) || ~isstruct(DataParams) || ...
                    ~isfield(DataParams, 'registration')
                return
            end

            reg = DataParams.registration;

            metadataRows(end+1, :) = {'Current folder', app.getCurrentDataFolder()}; %#ok<AGROW>
            metadataRows(end+1, :) = {'Image size [Y X]', sprintf('%d x %d', size(referenceImage, 1), size(referenceImage, 2))}; %#ok<AGROW>
            metadataRows(end+1, :) = {'Image class', class(referenceImage)}; %#ok<AGROW>

            finitePixels = referenceImage(isfinite(referenceImage));
            if isempty(finitePixels)
                metadataRows(end+1, :) = {'Intensity range', 'No finite pixels'}; %#ok<AGROW>
            else
                metadataRows(end+1, :) = {'Intensity range', ...
                    sprintf('%.6g to %.6g', min(finitePixels(:)), max(finitePixels(:)))}; %#ok<AGROW>
            end

            fieldList = { ...
                'referenceFile',        'Reference file'; ...
                'referenceDescription', 'Description'; ...
                'createdOn',            'Created on'; ...
                'source',               'Source'; ...
                'method',               'Method'; ...
                'notes',                'Notes'};

            for iField = 1:size(fieldList, 1)
                fieldName = fieldList{iField, 1};
                displayName = fieldList{iField, 2};

                if ~isfield(reg, fieldName)
                    continue
                end

                valueText = app.formatMetadataValueForDisplay(reg.(fieldName));
                if isempty(valueText)
                    continue
                end

                metadataRows(end+1, :) = {displayName, valueText}; %#ok<AGROW>
            end

        end

        function valueText = formatMetadataValueForDisplay(app, value) %#ok<INUSD>
            %FORMATMETADATAVALUEFORDISPLAY Convert metadata value to compact text.

            valueText = '';

            if isempty(value)
                return
            end

            if ischar(value)
                valueText = value;
                return
            end

            if isstring(value)
                valueText = char(strjoin(value(:), ', '));
                return
            end

            if isdatetime(value)
                valueText = char(string(value));
                return
            end

            if isduration(value)
                valueText = char(string(value));
                return
            end

            if islogical(value) && isscalar(value)
                valueText = mat2str(value);
                return
            end

            if isnumeric(value)
                if isscalar(value)
                    valueText = sprintf('%.6g', value);
                else
                    sz = size(value);
                    valueText = sprintf('%s numeric [%s]', ...
                        class(value), strjoin(cellstr(string(sz)), ' x '));
                end
                return
            end

            try
                valueText = char(string(value));
            catch
                valueText = sprintf('%s value', class(value));
            end

        end

        function showExistingImageReferencePreviewDialog(app, referenceImage, metadataRows)
            %SHOWEXISTINGIMAGEREFERENCEPREVIEWDIALOG Show saved reference image.
            %
            %   The dialog is read-only and does not modify DataParams.mat.

            referenceImage = squeeze(referenceImage);

            if ndims(referenceImage) ~= 2 || isempty(referenceImage)
                error('DataViewer:InvalidReferencePreview', ...
                    'Reference image preview must be a non-empty 2D image.');
            end

            referenceImage = single(referenceImage);
            referenceImage(~isfinite(referenceImage)) = 0;

            if nargin < 3 || isempty(metadataRows)
                metadataRows = {'Metadata', 'No metadata available'};
            end

            dlg = uifigure( ...
                'Name', 'Existing Image Reference', ...
                'WindowStyle', 'modal', ...
                'Position', [100 100 640 640], ...
                'Visible', 'off', ...
                'CloseRequestFcn', @onClose);

            grid = uigridlayout(dlg);
            grid.RowHeight = {28, '1x', 150, 38};
            grid.ColumnWidth = {'1x'};
            grid.Padding = [12 12 12 12];

            titleLabel = uilabel(grid);
            titleLabel.Text = 'Saved image reference';
            titleLabel.FontWeight = 'bold';
            titleLabel.Layout.Row = 1;
            titleLabel.Layout.Column = 1;

            previewAxes = uiaxes(grid);
            previewAxes.Layout.Row = 2;
            previewAxes.Layout.Column = 1;

            imagesc(previewAxes, referenceImage);
            axis(previewAxes, 'image');
            axis(previewAxes, 'off');

            try
                colormap(previewAxes, app.ImageAxes.Colormap);
            catch
            end

            try
                previewAxes.CLim = app.ImageAxes.CLim;
            catch
            end

            metaTable = uitable(grid);
            metaTable.Data = metadataRows;
            metaTable.ColumnName = {'Field', 'Value'};
            metaTable.ColumnEditable = [false false];
            metaTable.RowName = {};
            metaTable.Layout.Row = 3;
            metaTable.Layout.Column = 1;

            try
                metaTable.ColumnWidth = {150, 'auto'};
            catch
            end

            buttonGrid = uigridlayout(grid);
            buttonGrid.RowHeight = {'1x'};
            buttonGrid.ColumnWidth = {'1x', 120};
            buttonGrid.Padding = [0 0 0 0];
            buttonGrid.Layout.Row = 4;
            buttonGrid.Layout.Column = 1;

            closeButton = uibutton(buttonGrid, 'push');
            closeButton.Text = 'Close';
            closeButton.Layout.Row = 1;
            closeButton.Layout.Column = 2;
            closeButton.ButtonPushedFcn = @onClose;

            if exist('placeAppInsideCaller', 'file') == 2
                placeAppInsideCaller(app, dlg, 'center');
            end

            dlg.Visible = 'on';
            uiwait(dlg);

            if isvalid(dlg)
                delete(dlg);
            end

            function onClose(~, ~)
                uiresume(dlg);
            end

        end

        function showDataHistoryDialog(app, historyTable, sourceText, historyPath)
            %SHOWDATAHISTORYDIALOG Display data-history field/value table.

            if nargin < 3 || isempty(sourceText)
                sourceText = 'dataHistory.mat';
            end

            if nargin < 4 || isempty(historyPath)
                historyPath = '';
            end

            dlg = uifigure( ...
                'Name', 'Data History', ...
                'WindowStyle', 'modal', ...
                'Position', [100 100 760 560], ...
                'Visible', 'off', ...
                'CloseRequestFcn', @onClose);

            grid = uigridlayout(dlg);
            grid.RowHeight = {28, 24, 24, '1x', 38};
            grid.ColumnWidth = {'1x'};
            grid.Padding = [12 12 12 12];

            titleLabel = uilabel(grid);
            titleLabel.Text = 'Data history for current file';
            titleLabel.FontWeight = 'bold';
            titleLabel.Layout.Row = 1;
            titleLabel.Layout.Column = 1;

            fileLabel = uilabel(grid);
            fileLabel.Text = sprintf('Current file: %s', app.CurrentFile);
            fileLabel.Tooltip = app.CurrentFile;
            fileLabel.Layout.Row = 2;
            fileLabel.Layout.Column = 1;

            sourceLabel = uilabel(grid);
            if isempty(historyPath)
                sourceLabel.Text = sprintf('Matched entry: %s', sourceText);
            else
                sourceLabel.Text = sprintf('History file: %s | Matched entry: %s', ...
                    historyPath, sourceText);
                sourceLabel.Tooltip = historyPath;
            end
            sourceLabel.Layout.Row = 3;
            sourceLabel.Layout.Column = 1;

            historyUITable = uitable(grid);
            historyUITable.Data = historyTable;
            historyUITable.ColumnName = historyTable.Properties.VariableNames;
            historyUITable.ColumnEditable = false(1, width(historyTable));
            historyUITable.RowName = {};
            historyUITable.Layout.Row = 4;
            historyUITable.Layout.Column = 1;

            try
                historyUITable.ColumnWidth = {220, 'auto'};
            catch
            end

            buttonGrid = uigridlayout(grid);
            buttonGrid.RowHeight = {'1x'};
            buttonGrid.ColumnWidth = {'1x', 120};
            buttonGrid.Padding = [0 0 0 0];
            buttonGrid.Layout.Row = 5;
            buttonGrid.Layout.Column = 1;

            closeButton = uibutton(buttonGrid, 'push');
            closeButton.Text = 'Close';
            closeButton.Layout.Row = 1;
            closeButton.Layout.Column = 2;
            closeButton.ButtonPushedFcn = @onClose;

            if exist('placeAppInsideCaller', 'file') == 2
                placeAppInsideCaller(app, dlg, 'center');
            end

            dlg.Visible = 'on';
            uiwait(dlg);

            if isvalid(dlg)
                delete(dlg);
            end

            function onClose(~, ~)
                uiresume(dlg);
            end

        end

        function openCurrentFileDataHistoryDialog(app)
            %OPENCURRENTFILEDATAHISTORYDIALOG Display pipeline history for current file.
            %
            %   openCurrentFileDataHistoryDialog(app) loads dataHistory.mat from the
            %   current data folder, finds the dataHistory element associated with the
            %   currently opened file, and displays the function execution history in
            %   chronological order.
            %
            %   Expected schema:
            %       dataHistory(i).file
            %       dataHistory(i).timestamp
            %       dataHistory(i).pipelineVersion
            %       dataHistory(i).info(j).name
            %       dataHistory(i).info(j).version
            %       dataHistory(i).info(j).inputs(k).name
            %       dataHistory(i).info(j).inputs(k).source
            %       dataHistory(i).info(j).parameters(k).name
            %       dataHistory(i).info(j).parameters(k).value
            %
            %   The dialog is read-only and modal. Each table row represents one
            %   executed function from dataHistory(i).info.

            if ~app.hasData() || isempty(app.CurrentFile)
                app.setStatusMessage('Load image data before viewing data history.');
                return
            end

            folderPath = app.getCurrentDataFolder();
            dataHistoryPath = fullfile(folderPath, 'dataHistory.mat');

            if ~isfile(dataHistoryPath)
                app.setStatusMessage('dataHistory.mat not found for current data folder.');

                uialert(app.UIFigure, ...
                    sprintf('dataHistory.mat was not found in:\n\n%s', folderPath), ...
                    'Data History Not Found', ...
                    'Icon', 'warning');
                return
            end

            S = load(dataHistoryPath, 'dataHistory');

            if ~isfield(S, 'dataHistory') || ~isstruct(S.dataHistory) || isempty(S.dataHistory)
                app.setStatusMessage('dataHistory.mat does not contain a valid dataHistory struct.');

                uialert(app.UIFigure, ...
                    'dataHistory.mat does not contain a valid non-empty dataHistory struct.', ...
                    'Invalid Data History', ...
                    'Icon', 'warning');
                return
            end

            dataHistory = S.dataHistory(:);

            [~, currentName, currentExt] = fileparts(app.CurrentFile);
            currentFileName = [currentName currentExt];
            currentFullPath = app.CurrentFile;

            matchIdx = [];

            for iEntry = 1:numel(dataHistory)
                if ~isfield(dataHistory(iEntry), 'file') || isempty(dataHistory(iEntry).file)
                    continue
                end

                historyFile = char(string(dataHistory(iEntry).file));
                [~, historyName, historyExt] = fileparts(historyFile);
                historyFileName = [historyName historyExt];

                if strcmpi(historyFile, currentFullPath) || ...
                        strcmpi(historyFile, currentFileName) || ...
                        strcmpi(historyFileName, currentFileName)
                    matchIdx = iEntry;
                    break
                end
            end

            if isempty(matchIdx)
                app.setStatusMessage('No dataHistory entry found for the current file.');

                uialert(app.UIFigure, ...
                    sprintf('No dataHistory entry matched the currently opened file:\n\n%s', ...
                    app.CurrentFile), ...
                    'Data History Entry Not Found', ...
                    'Icon', 'warning');
                return
            end

            entry = dataHistory(matchIdx);

            if ~isfield(entry, 'info') || ~isstruct(entry.info) || isempty(entry.info)
                app.setStatusMessage('Matched dataHistory entry has no valid info struct.');

                uialert(app.UIFigure, ...
                    'The matched dataHistory entry does not contain a valid info struct.', ...
                    'Invalid Data History Entry', ...
                    'Icon', 'warning');
                return
            end

            infoList = entry.info(:);
            nSteps = numel(infoList);

            tableData = cell(nSteps, 5);

            for iStep = 1:nSteps
                thisInfo = infoList(iStep);

                functionName = iGetTextField(thisInfo, 'name', '');
                functionVersion = iGetTextField(thisInfo, 'version', '');

                inputText = '';
                if isfield(thisInfo, 'inputs') && ~isempty(thisInfo.inputs)
                    inputText = iFormatInputs(thisInfo.inputs);
                end

                parameterText = '';
                if isfield(thisInfo, 'parameters') && ~isempty(thisInfo.parameters)
                    parameterText = iFormatParameters(thisInfo.parameters);
                end

                tableData(iStep, :) = { ...
                    iStep, ...
                    functionName, ...
                    functionVersion, ...
                    inputText, ...
                    parameterText};
            end

            displayTable = cell2table(tableData, ...
                'VariableNames', {'Step', 'Function', 'Version', 'Inputs', 'Parameters'});

            fileText = iGetTextField(entry, 'file', currentFileName);
            timestampText = iFormatValue(iGetFieldValue(entry, 'timestamp', ''));
            pipelineVersionText = iGetTextField(entry, 'pipelineVersion', '');

            dlg = uifigure( ...
                'Name', 'Data History', ...
                'WindowStyle', 'modal', ...
                'Position', [100 100 900 600], ...
                'Visible', 'off', ...
                'CloseRequestFcn', @onClose);

            grid = uigridlayout(dlg);
            grid.RowHeight = {28, 26, 26, 26, '1x', 38};
            grid.ColumnWidth = {'1x'};
            grid.Padding = [12 12 12 12];

            titleLabel = uilabel(grid);
            titleLabel.Text = 'Data history for current file';
            titleLabel.FontWeight = 'bold';
            titleLabel.Layout.Row = 1;
            titleLabel.Layout.Column = 1;

            fileLabel = uilabel(grid);
            fileLabel.Text = sprintf('File: %s', fileText);
            fileLabel.Tooltip = fileText;
            fileLabel.Layout.Row = 2;
            fileLabel.Layout.Column = 1;

            timestampLabel = uilabel(grid);
            timestampLabel.Text = sprintf('Created on: %s', timestampText);
            timestampLabel.Layout.Row = 3;
            timestampLabel.Layout.Column = 1;

            versionLabel = uilabel(grid);
            versionLabel.Text = sprintf('Pipeline version: %s', pipelineVersionText);
            versionLabel.Layout.Row = 4;
            versionLabel.Layout.Column = 1;

            historyUITable = uitable(grid);
            historyUITable.Data = displayTable;
            historyUITable.ColumnName = displayTable.Properties.VariableNames;
            historyUITable.ColumnEditable = false(1, width(displayTable));
            historyUITable.RowName = {};
            historyUITable.Layout.Row = 5;
            historyUITable.Layout.Column = 1;

            try
                historyUITable.ColumnWidth = {55, 180, 90, 260, 300};
            catch
            end

            buttonGrid = uigridlayout(grid);
            buttonGrid.RowHeight = {'1x'};
            buttonGrid.ColumnWidth = {'1x', 120};
            buttonGrid.Padding = [0 0 0 0];
            buttonGrid.Layout.Row = 6;
            buttonGrid.Layout.Column = 1;

            closeButton = uibutton(buttonGrid, 'push');
            closeButton.Text = 'Close';
            closeButton.Layout.Row = 1;
            closeButton.Layout.Column = 2;
            closeButton.ButtonPushedFcn = @onClose;

            if exist('placeAppInsideCaller', 'file') == 2
                placeAppInsideCaller(app, dlg, 'center');
            end

            dlg.Visible = 'on';
            uiwait(dlg);

            if isvalid(dlg)
                delete(dlg);
            end

            app.setStatusMessage('Displayed data history for current file.');

            function onClose(~, ~)
                uiresume(dlg);
            end

            function value = iGetFieldValue(s, fieldName, defaultValue)
                value = defaultValue;

                if isstruct(s) && isfield(s, fieldName) && ~isempty(s.(fieldName))
                    value = s.(fieldName);
                end
            end

            function txt = iGetTextField(s, fieldName, defaultValue)
                txt = iFormatValue(iGetFieldValue(s, fieldName, defaultValue));
            end

            function txt = iFormatInputs(inputs)
                inputs = iAsStructArray(inputs);
                parts = strings(0, 1);

                for iInput = 1:numel(inputs)
                    inputName = iGetTextField(inputs(iInput), 'name', '');
                    inputSource = iGetTextField(inputs(iInput), 'source', '');

                    if isempty(inputName) && isempty(inputSource)
                        continue
                    elseif isempty(inputName)
                        parts(end+1, 1) = string(inputSource); %#ok<AGROW>
                    elseif isempty(inputSource)
                        parts(end+1, 1) = string(inputName); %#ok<AGROW>
                    else
                        parts(end+1, 1) = sprintf('%s: %s', inputName, inputSource); %#ok<AGROW>
                    end
                end

                txt = char(strjoin(parts, '; '));
            end

            function txt = iFormatParameters(parameters)
                parameters = iAsStructArray(parameters);
                parts = strings(0, 1);

                for iParam = 1:numel(parameters)
                    paramName = iGetTextField(parameters(iParam), 'name', '');
                    paramValue = iFormatValue(iGetFieldValue(parameters(iParam), 'value', ''));

                    if isempty(paramName) && isempty(paramValue)
                        continue
                    elseif isempty(paramName)
                        parts(end+1, 1) = string(paramValue); %#ok<AGROW>
                    elseif isempty(paramValue)
                        parts(end+1, 1) = string(paramName); %#ok<AGROW>
                    else
                        parts(end+1, 1) = sprintf('%s = %s', paramName, paramValue); %#ok<AGROW>
                    end
                end

                txt = char(strjoin(parts, '; '));
            end

            function structArray = iAsStructArray(value)
                if istable(value)
                    structArray = table2struct(value);
                elseif isstruct(value)
                    structArray = value(:);
                else
                    structArray = struct([]);
                end
            end

            function txt = iFormatValue(value)
                if isempty(value)
                    txt = '';
                    return
                end

                if ischar(value)
                    txt = value;
                    return
                end

                if isstring(value)
                    txt = char(strjoin(value(:), ', '));
                    return
                end

                if isdatetime(value) || isduration(value)
                    txt = char(strjoin(string(value(:)), ', '));
                    return
                end

                if islogical(value)
                    txt = mat2str(value);
                    return
                end

                if isnumeric(value)
                    if isscalar(value)
                        txt = sprintf('%.8g', value);
                    else
                        txt = mat2str(value);
                    end
                    return
                end

                if iscell(value)
                    if iscellstr(value) || all(cellfun(@(x) ischar(x) || ...
                            (isstring(x) && isscalar(x)), value(:)))
                        txt = char(strjoin(string(value(:)), ', '));
                    else
                        txt = sprintf('cell [%s]', ...
                            strjoin(cellstr(string(size(value))), ' x '));
                    end
                    return
                end

                if isstruct(value)
                    txt = sprintf('struct [%s]', ...
                        strjoin(cellstr(string(size(value))), ' x '));
                    return
                end

                try
                    txt = char(string(value));
                catch
                    txt = sprintf('%s value', class(value));
                end
            end

        end

        function ctx = getImageReferenceManagerContext(app)
            %GETIMAGEREFERENCEMANAGERCONTEXT Build current-frame context for ImageReferenceManager.
            %
            % Output:
            %   ctx - Scalar struct with the current displayed frame and source
            %         metadata used to create ImageReference files.
            %
            % Notes:
            %   This uses the viewer's current display state. In event mode,
            %   app.getCurrentFrame() already returns the currently displayed
            %   event/average frame.

            frame = app.getCurrentFrame();

            if ~(isnumeric(frame) || islogical(frame)) || ndims(frame) ~= 2
                error('DataViewer:InvalidImageReferenceFrame', ...
                    'Current frame must be a 2D numeric or logical image.');
            end

            ctx = struct();
            ctx.currentFrameImage = single(frame);
            ctx.sourceFile = app.CurrentFile;
            ctx.sourceFolder = app.getCurrentDataFolder();
            ctx.sourceFrame = app.CurrentFrame;
            ctx.sourceType = app.getSourceType();
        end

        function ctx = getImageAlignmentToolContext(app)
            %GETIMAGEALIGNMENTTOOLCONTEXT Build ImageAlignmentTool launch context.

            frame = app.getCurrentFrame();

            if ~(isnumeric(frame) || islogical(frame)) || ndims(frame) ~= 2
                error('DataViewer:InvalidAlignmentMovingFrame', ...
                    'Current DataViewer frame must be a 2D numeric or logical image.');
            end

            ctx = struct();
            ctx.movingImage = single(frame);
            ctx.sourceFile = app.CurrentFile;
            ctx.targetFolder = app.getCurrentDataFolder();

            try
                ctx.sourceFrame = app.CurrentFrame;
            catch
                ctx.sourceFrame = NaN;
            end

            ctx.sourceType = app.getSourceType();

            if isempty(ctx.sourceFile)
                sourceText = 'current DataViewer source';
            else
                sourceText = ctx.sourceFile;
            end

            if isfinite(ctx.sourceFrame)
                ctx.description = sprintf('Current DataViewer frame %g from %s', ...
                    ctx.sourceFrame, sourceText);
            else
                ctx.description = sprintf('Current DataViewer frame from %s', sourceText);
            end
        end

        function saveFolder = getPipelineSaveFolderForCurrentData(app)
            %GETPIPELINESAVEFOLDERFORCURRENTDATA Backward-compatible wrapper.

            saveFolder = app.getCurrentSaveFolderForData();
        end

        function saveFolder = getCurrentSaveFolderForData(app)
            %GETCURRENTSAVEFOLDERFORDATA Return SaveFolder for current loaded data.

            saveFolder = '';

            if isempty(app.CurrentFile)
                return
            end

            try
                saveFolder = fileparts(app.CurrentFile);
            catch
                saveFolder = '';
            end

            if ~isempty(saveFolder) && ~isfolder(saveFolder)
                saveFolder = '';
            end
        end

        function printExceptionToConsole(app, ME, contextLabel) %#ok<INUSL>
            %PRINTEXCEPTIONTOCONSOLE Print full MATLAB exception report.
            %
            %   App Designer alerts usually hide the call stack. This helper prints
            %   the full error report to the MATLAB Command Window for debugging.

            if nargin < 3 || isempty(contextLabel)
                contextLabel = 'DataViewer error';
            end

            fprintf(2, '\n');
            fprintf(2, '============================================================\n');
            fprintf(2, '%s\n', char(string(contextLabel)));
            fprintf(2, '============================================================\n');
            fprintf(2, '%s\n', getReport(ME, 'extended', 'hyperlinks', 'on'));
            fprintf(2, '============================================================\n\n');
        end

        function [rawFolder, didResolve] = resolvePipelineRawFolderForLaunch(app, saveFolder)
            %RESOLVEPIPELINERAWFOLDERFORLAUNCH Backward-compatible wrapper.
            %
            % The RawFolder is now a DataViewer data-context property, not a pipeline-only
            % launch prompt. This wrapper remains only to avoid breaking older callbacks.

            rawFolder = app.resolveDataRawFolderForCurrentData(saveFolder);
            didResolve = true;
        end

        function rawFolder = resolveDataRawFolderForCurrentData(app, saveFolder)
            %RESOLVEDATARAWFOLDERFORCURRENTDATA Return RawFolder for current data context.
            %
            % DataViewer owns RawFolder as a shared SaveFolder-level data-context value.
            % The persistent source of truth is DataParams.folders.RawFolder. Runtime
            % app.DataRawFolder is only a cache bound to one SaveFolder.
            %
            % No modal prompt is shown here. Invalid stored paths resolve to "Missing".

            if nargin < 2 || isempty(saveFolder)
                saveFolder = app.getCurrentSaveFolderForData();
            end

            saveFolder = char(string(saveFolder));
            rawFolder = 'Missing';

            if isempty(saveFolder)
                return
            end

            % Reuse explicit app-level RawFolder only when it belongs to this SaveFolder.
            try
                if ~isempty(app.DataRawFolder) && ...
                        strcmpi(char(string(app.DataRawFolderSaveFolder)), saveFolder)

                    rawCandidate = char(string(app.DataRawFolder));

                    if strcmpi(rawCandidate, 'Missing')
                        rawFolder = 'Missing';
                        return
                    end

                    if isfolder(rawCandidate)
                        rawFolder = app.normalizeFolderPath(rawCandidate);
                        return
                    end
                end
            catch
            end

            % Preferred persistent source: DataParams.folders.RawFolder.
            try
                app.ensureDataParamsFolderFields();
                [storedRawFolder, storedStatus] = app.getStoredRawFolderFromDataParams();

                if strcmpi(storedRawFolder, 'Missing') || strcmpi(storedStatus, 'missing')
                    rawFolder = 'Missing';
                    app.DataRawFolder = 'Missing';
                    app.DataRawFolderSaveFolder = saveFolder;
                    return
                end

                if isfolder(storedRawFolder)
                    rawFolder = app.normalizeFolderPath(storedRawFolder);
                    app.DataRawFolder = rawFolder;
                    app.DataRawFolderSaveFolder = saveFolder;
                    return
                end
            catch
            end

            % Legacy fallback for older in-memory DataParams structs.
            try
                if isfield(app.DataParams, 'RawFolder') && ~isempty(app.DataParams.RawFolder)
                    legacyRawFolder = char(string(app.DataParams.RawFolder));
                elseif isfield(app.DataParams, 'rawFolder') && ~isempty(app.DataParams.rawFolder)
                    legacyRawFolder = char(string(app.DataParams.rawFolder));
                else
                    legacyRawFolder = '';
                end

                if strcmpi(legacyRawFolder, 'Missing')
                    rawFolder = 'Missing';
                    app.DataRawFolder = 'Missing';
                    app.DataRawFolderSaveFolder = saveFolder;
                    return
                end

                if isfolder(legacyRawFolder)
                    rawFolder = app.normalizeFolderPath(legacyRawFolder);
                    app.DataRawFolder = rawFolder;
                    app.DataRawFolderSaveFolder = saveFolder;
                    return
                end
            catch
            end

            rawFolder = 'Missing';
        end

        function [tf, rawFolder, reasonText] = canOpenEventsManager(app)
            %CANOPENEVENTSMANAGER True when Events Manager can be launched.

            tf = false;
            rawFolder = 'Missing';
            reasonText = '';

            if ~app.hasData()
                reasonText = 'Load .dat data before opening Events Manager.';
                return
            end

            if ~strcmpi(app.getSourceType(), 'dat')
                reasonText = 'Events Manager is only available for .dat data.';
                return
            end

            saveFolder = app.getCurrentSaveFolderForData();
            if isempty(saveFolder)
                reasonText = 'The current data file does not have a valid SaveFolder.';
                return
            end

            [rawFolder, bRawFolderSet] = app.getExplicitRawFolderForCurrentData(saveFolder);

            if ~bRawFolderSet
                reasonText = 'Set RawFolder before opening Events Manager.';
                return
            end

            if strcmpi(rawFolder, 'Missing')
                tf = true;
                return
            end

            if isempty(rawFolder) || ~isfolder(rawFolder)
                reasonText = 'The configured RawFolder is not a valid folder.';
                return
            end

            if ~app.rawFolderContainsAnalogInputBins(rawFolder)
                reasonText = 'Events Manager requires ai_(number).bin files in RawFolder, or RawFolder set to Missing.';
                return
            end

            tf = true;

        end

        function [tf, rawFolder, reasonText] = canOpenOiSDualCamCoreg(app)
            %CANOPENOISDUALCAMCOREG True when dual-camera coregistration can launch.

            tf = false;
            rawFolder = 'Missing';
            reasonText = '';

            if ~app.hasData()
                reasonText = 'Load dual-camera .dat data before opening coregistration.';
                return
            end

            if ~strcmpi(app.getSourceType(), 'dat')
                reasonText = 'Dual-camera coregistration is only available for .dat data.';
                return
            end

            saveFolder = app.getCurrentSaveFolderForData();
            if isempty(saveFolder) || ~isfolder(saveFolder)
                reasonText = 'The current data file does not have a valid SaveFolder.';
                return
            end

            rawFolder = app.resolveDataRawFolderForCurrentData(saveFolder);

            if isempty(rawFolder) || strcmpi(rawFolder, 'Missing') || ~isfolder(rawFolder)
                reasonText = 'Dual-camera coregistration requires a valid RawFolder.';
                return
            end

            if ~app.currentDatSourceIsMultiCam()
                reasonText = 'The current .dat metadata do not report MultiCam = true.';
                return
            end

            tf = true;

        end

        function [rawFolder, bIsExplicit] = getExplicitRawFolderForCurrentData(app, saveFolder)
            %GETEXPLICITRAWFOLDERFORCURRENTDATA Return explicitly configured RawFolder.
            %
            %   Unlike resolveDataRawFolderForCurrentData, this method distinguishes
            %   between an explicit "Missing" RawFolder and the default unresolved state.

            if nargin < 2 || isempty(saveFolder)
                saveFolder = app.getCurrentSaveFolderForData();
            end

            saveFolder = char(string(saveFolder));
            rawFolder = 'Missing';
            bIsExplicit = false;

            if isempty(saveFolder)
                return
            end

            % Preferred app-level SaveFolder-bound RawFolder.
            try
                if ~isempty(app.DataRawFolder) && ...
                        strcmpi(char(string(app.DataRawFolderSaveFolder)), saveFolder)

                    rawCandidate = char(string(app.DataRawFolder));

                    if strcmpi(rawCandidate, 'Missing')
                        rawFolder = 'Missing';
                        bIsExplicit = true;
                        return
                    end

                    if isfolder(rawCandidate)
                        rawFolder = rawCandidate;
                        bIsExplicit = true;
                        return
                    end

                    rawFolder = rawCandidate;
                    bIsExplicit = true;
                    return
                end
            catch
            end

            % Transitional fallback: DataParams may carry a persisted RawFolder value.
            try
                if isstruct(app.DataParams)
                    if isfield(app.DataParams, 'RawFolder') && ~isempty(app.DataParams.RawFolder)
                        rawCandidate = char(string(app.DataParams.RawFolder));
                        bIsExplicit = true;
                    elseif isfield(app.DataParams, 'rawFolder') && ~isempty(app.DataParams.rawFolder)
                        rawCandidate = char(string(app.DataParams.rawFolder));
                        bIsExplicit = true;
                    else
                        rawCandidate = '';
                    end

                    if bIsExplicit
                        if strcmpi(rawCandidate, 'Missing')
                            rawFolder = 'Missing';
                            app.DataRawFolder = 'Missing';
                            app.DataRawFolderSaveFolder = saveFolder;
                            return
                        end

                        if isfolder(rawCandidate)
                            rawFolder = rawCandidate;
                            app.DataRawFolder = rawFolder;
                            app.DataRawFolderSaveFolder = saveFolder;
                            return
                        end

                        rawFolder = rawCandidate;
                        return
                    end
                end
            catch
            end

        end

        function ensureDataParamsFolderFields(app)
            %ENSUREDATAPARAMSFOLDERFIELDS Ensure DataParams.folders RawFolder schema.

            if isempty(app.DataParams) || ~isstruct(app.DataParams)
                app.DataParams = struct();
            end

            if ~isfield(app.DataParams, 'folders') || ...
                    ~isstruct(app.DataParams.folders) || ...
                    ~isscalar(app.DataParams.folders)
                app.DataParams.folders = struct();
            end

            if ~isfield(app.DataParams.folders, 'RawFolder') || ...
                    isempty(app.DataParams.folders.RawFolder)

                legacyRawFolder = '';
                try
                    if isfield(app.DataParams, 'RawFolder') && ~isempty(app.DataParams.RawFolder)
                        legacyRawFolder = app.DataParams.RawFolder;
                    elseif isfield(app.DataParams, 'rawFolder') && ~isempty(app.DataParams.rawFolder)
                        legacyRawFolder = app.DataParams.rawFolder;
                    end
                catch
                    legacyRawFolder = '';
                end

                if isempty(legacyRawFolder)
                    app.DataParams.folders.RawFolder = 'Missing';
                else
                    app.DataParams.folders.RawFolder = char(string(legacyRawFolder));
                end
            end

            app.DataParams.folders.RawFolder = char(string(app.DataParams.folders.RawFolder));

            if ~isfield(app.DataParams.folders, 'RawFolderStatus') || ...
                    isempty(app.DataParams.folders.RawFolderStatus)
                [~, statusText] = app.normalizeRawFolderValue(app.DataParams.folders.RawFolder);
                app.DataParams.folders.RawFolderStatus = statusText;
            else
                app.DataParams.folders.RawFolderStatus = lower(strtrim(char(string(app.DataParams.folders.RawFolderStatus))));
            end

            if strcmpi(app.DataParams.folders.RawFolder, 'Missing')
                app.DataParams.folders.RawFolder = 'Missing';
                app.DataParams.folders.RawFolderStatus = 'missing';
            elseif strcmpi(app.DataParams.folders.RawFolderStatus, 'missing')
                [~, statusText] = app.normalizeRawFolderValue(app.DataParams.folders.RawFolder);
                app.DataParams.folders.RawFolderStatus = statusText;
            end

            if ~isfield(app.DataParams.folders, 'RawFolderSetOn')
                app.DataParams.folders.RawFolderSetOn = [];
            end

            if ~isfield(app.DataParams.folders, 'RawFolderSetBy') || ...
                    isempty(app.DataParams.folders.RawFolderSetBy)
                app.DataParams.folders.RawFolderSetBy = '';
            else
                app.DataParams.folders.RawFolderSetBy = char(string(app.DataParams.folders.RawFolderSetBy));
            end

        end

        function [rawFolder, statusText] = getStoredRawFolderFromDataParams(app)
            %GETSTOREDRAWFOLDERFROMDATAPARAMS Return DataParams.folders RawFolder.

            app.ensureDataParamsFolderFields();

            rawFolder = char(string(app.DataParams.folders.RawFolder));
            statusText = lower(strtrim(char(string(app.DataParams.folders.RawFolderStatus))));

            if isempty(rawFolder)
                rawFolder = 'Missing';
                statusText = 'missing';
            end

            if isempty(statusText)
                [~, statusText] = app.normalizeRawFolderValue(rawFolder);
            end

        end

        function [rawFolder, statusText] = normalizeRawFolderValue(app, rawFolder) %#ok<INUSL>
            %NORMALIZERAWFOLDERVALUE Normalize RawFolder text and infer status.

            if nargin < 2 || isempty(rawFolder)
                rawFolder = 'Missing';
            end

            rawFolder = char(string(rawFolder));
            rawFolder = strtrim(rawFolder);

            if isempty(rawFolder) || strcmpi(rawFolder, 'Missing')
                rawFolder = 'Missing';
                statusText = 'missing';
                return
            end

            if isfolder(rawFolder)
                try
                    [ok, folderInfo] = fileattrib(rawFolder);
                    if ok && isfield(folderInfo, 'Name') && ~isempty(folderInfo.Name)
                        rawFolder = folderInfo.Name;
                    end
                catch
                end

                statusText = 'valid';
                return
            end

            statusText = 'invalid';

        end

        function folderPath = normalizeFolderPath(app, folderPath) %#ok<INUSL>
            %NORMALIZEFOLDERPATH Return canonical folder path when possible.

            folderPath = char(string(folderPath));

            if isempty(folderPath) || ~isfolder(folderPath)
                return
            end

            try
                [ok, folderInfo] = fileattrib(folderPath);
                if ok && isfield(folderInfo, 'Name') && ~isempty(folderInfo.Name)
                    folderPath = folderInfo.Name;
                end
            catch
            end
        end

        function setRawFolderForCurrentSaveFolder(app, rawFolder, setBy, bShowStatus, statusOverride)
            %SETRAWFOLDERFORCURRENTSAVEFOLDER Persist and apply current RawFolder.

            if nargin < 3 || isempty(setBy)
                setBy = 'DataViewer';
            end

            if nargin < 4 || isempty(bShowStatus)
                bShowStatus = true;
            end

            if nargin < 5
                statusOverride = '';
            end

            saveFolder = app.getCurrentSaveFolderForData();

            if isempty(saveFolder) || ~isfolder(saveFolder)
                app.setStatusMessage('RawFolder was not updated: current SaveFolder is invalid.');
                return
            end

            if isempty(statusOverride)
                [rawFolderStored, statusText] = app.normalizeRawFolderValue(rawFolder);
            else
                [rawFolderStored, statusText] = app.normalizeRawFolderValue(rawFolder);
                statusText = lower(strtrim(char(string(statusOverride))));
                if strcmp(statusText, 'missing')
                    rawFolderStored = 'Missing';
                end
            end

            app.persistRawFolderToDataParams(saveFolder, rawFolderStored, setBy, statusText);

            if strcmp(statusText, 'valid')
                app.DataRawFolder = rawFolderStored;
            else
                app.DataRawFolder = 'Missing';
            end

            app.DataRawFolderSaveFolder = saveFolder;
            app.updatePipelineRawFolderContext(saveFolder, app.DataRawFolder);

            app.updateDataFolderContextLabel();
            app.updatePipelineTabState();
            app.updateGUIEnabledState();

            if bShowStatus
                switch statusText
                    case 'valid'
                        if strcmpi(rawFolderStored, saveFolder)
                            app.setStatusMessage('RawFolder set to current SaveFolder.');
                        else
                            app.setStatusMessage(sprintf('RawFolder set: %s', rawFolderStored));
                        end
                    case 'missing'
                        app.setStatusMessage('RawFolder set to Missing.');
                    otherwise
                        app.setStatusMessage('Stored RawFolder is invalid. Runtime RawFolder set to Missing.');
                end
            end

        end

        function persistRawFolderToDataParams(app, saveFolder, rawFolder, setBy, statusText)
            %PERSISTRAWFOLDERTODATAPARAMS Save RawFolder context to DataParams.mat.

            if nargin < 4 || isempty(setBy)
                setBy = 'DataViewer';
            end

            saveFolder = char(string(saveFolder));

            if isempty(saveFolder) || ~isfolder(saveFolder)
                error('DataViewer:InvalidSaveFolder', ...
                    'Cannot persist RawFolder because SaveFolder is invalid.');
            end

            if nargin < 5 || isempty(statusText)
                [rawFolder, statusText] = app.normalizeRawFolderValue(rawFolder);
            else
                [rawFolder, inferredStatus] = app.normalizeRawFolderValue(rawFolder);
                statusText = lower(strtrim(char(string(statusText))));

                if isempty(statusText)
                    statusText = inferredStatus;
                end
            end

            if ~ismember(statusText, {'valid', 'missing', 'invalid'})
                error('DataViewer:InvalidRawFolderStatus', ...
                    'RawFolderStatus must be valid, missing, or invalid.');
            end

            if strcmp(statusText, 'missing')
                rawFolder = 'Missing';
            end

            if ~isfile(fullfile(saveFolder, 'DataParams.mat'))
                DataParams = createDataParams(saveFolder);
            else
                DataParams = loadDataParams(saveFolder);
            end

            app.DataParams = DataParams;
            app.ensureDataParamsViewFields();
            app.ensureDataParamsMaskFields();
            app.ensureDataParamsFolderFields();

            app.DataParams.folders.RawFolder = char(string(rawFolder));
            app.DataParams.folders.RawFolderStatus = statusText;
            app.DataParams.folders.RawFolderSetOn = datetime('now');
            app.DataParams.folders.RawFolderSetBy = char(string(setBy));

            saveDataParams(saveFolder, app.DataParams);

            % Keep runtime state aligned with the saved, normalized file.
            currentSaveFolder = app.getCurrentSaveFolderForData();
            if ~isempty(currentSaveFolder) && strcmpi(currentSaveFolder, saveFolder)
                app.DataParams = loadDataParams(saveFolder);
                app.ensureDataParamsViewFields();
                app.ensureDataParamsMaskFields();
                app.ensureDataParamsFolderFields();
            end

        end

        function synchronizeRawFolderFromDataParams(app)
            %SYNCHRONIZERAWFOLDERFROMDATAPARAMS Adopt persisted RawFolder at runtime.

            if ~app.hasData() || isempty(app.CurrentFile)
                return
            end

            saveFolder = app.getCurrentSaveFolderForData();
            if isempty(saveFolder) || ~isfolder(saveFolder)
                return
            end

            app.ensureDataParamsFolderFields();
            [storedRawFolder, storedStatus] = app.getStoredRawFolderFromDataParams();

            if strcmpi(storedRawFolder, 'Missing') || strcmpi(storedStatus, 'missing')
                app.DataRawFolder = 'Missing';
                app.DataRawFolderSaveFolder = saveFolder;
                app.updatePipelineRawFolderContext(saveFolder, 'Missing');
                return
            end

            if isfolder(storedRawFolder)
                normalizedRawFolder = app.normalizeFolderPath(storedRawFolder);
                app.DataRawFolder = normalizedRawFolder;
                app.DataRawFolderSaveFolder = saveFolder;
                app.updatePipelineRawFolderContext(saveFolder, normalizedRawFolder);

                if ~strcmp(storedRawFolder, normalizedRawFolder) || ~strcmpi(storedStatus, 'valid')
                    app.persistRawFolderToDataParams(saveFolder, normalizedRawFolder, 'DataViewer', 'valid');
                end
                return
            end

            % Stored path is no longer reachable. Do not use it at runtime, but retain
            % it in DataParams with invalid status so the old location remains visible.
            app.DataRawFolder = 'Missing';
            app.DataRawFolderSaveFolder = saveFolder;
            app.updatePipelineRawFolderContext(saveFolder, 'Missing');

            try
                app.persistRawFolderToDataParams(saveFolder, storedRawFolder, 'DataViewer', 'invalid');
            catch
            end

            if ~strcmpi(app.RawFolderInvalidWarningSaveFolder, saveFolder)
                app.RawFolderInvalidWarningSaveFolder = saveFolder;
                app.promptForInvalidStoredRawFolder(storedRawFolder, saveFolder);
            end

        end

        function promptForInvalidStoredRawFolder(app, storedRawFolder, saveFolder)
            %PROMPTFORINVALIDSTOREDRAWFOLDER Warn and optionally repair invalid RawFolder.

            storedRawFolder = char(string(storedRawFolder));
            saveFolder = char(string(saveFolder));

            msg = sprintf(['The stored RawFolder path no longer exists:\n\n%s\n\n' ...
                'This can happen when raw data are on an external drive and the drive letter changes.\n\n' ...
                'Select the current RawFolder, set it to Missing, or ignore for now.'], ...
                storedRawFolder);

            try
                choice = uiconfirm(app.UIFigure, ...
                    msg, ...
                    'RawFolder not found', ...
                    'Options', {'Select Folder', 'Set Missing', 'Ignore'}, ...
                    'DefaultOption', 1, ...
                    'CancelOption', 3, ...
                    'Icon', 'warning');
            catch
                app.setStatusMessage('Stored RawFolder is invalid. Runtime RawFolder set to Missing.');
                return
            end

            switch choice
                case 'Select Folder'
                    startFolder = saveFolder;

                    selectedFolder = uigetdir(startFolder, 'Select current RawFolder');
                    if isequal(selectedFolder, 0)
                        app.setStatusMessage('RawFolder repair cancelled. Runtime RawFolder set to Missing.');
                        return
                    end

                    if isempty(selectedFolder) || ~isfolder(selectedFolder)
                        app.setStatusMessage('Selected RawFolder is invalid. Runtime RawFolder set to Missing.');
                        return
                    end

                    app.setRawFolderForCurrentSaveFolder(selectedFolder, 'DataViewer', true);

                case 'Set Missing'
                    app.setRawFolderForCurrentSaveFolder('Missing', 'DataViewer', true);

                otherwise
                    app.setStatusMessage('Stored RawFolder is invalid. Runtime RawFolder set to Missing.');
            end

        end

        function updatePipelineRawFolderContext(app, saveFolder, rawFolder)
            %UPDATEPIPELINERAWFOLDERCONTEXT Synchronize PipelineManager RawFolder context.

            if isempty(app.PipelineManagerObj) || ~isa(app.PipelineManagerObj, 'PipelineManager')
                return
            end

            saveFolder = char(string(saveFolder));
            rawFolder = char(string(rawFolder));

            try
                if isprop(app.PipelineManagerObj, 'SaveFolderList') && ...
                        isprop(app.PipelineManagerObj, 'RawFolderList') && ...
                        ~isempty(app.PipelineManagerObj.SaveFolderList) && ...
                        ~isempty(app.PipelineManagerObj.RawFolderList)

                    for iSF = 1:numel(app.PipelineManagerObj.SaveFolderList)
                        try
                            if strcmpi(char(string(app.PipelineManagerObj.SaveFolderList{iSF})), saveFolder)
                                app.PipelineManagerObj.RawFolderList{iSF} = rawFolder;
                            end
                        catch
                        end
                    end
                end
            catch
            end

        end

        function tf = rawFolderContainsAnalogInputBins(app, rawFolder) %#ok<INUSL>
            %RAWFOLDERCONTAINSANALOGINPUTBINS True for ai_(number).bin files.

            tf = false;

            if isempty(rawFolder) || ~isfolder(rawFolder)
                return
            end

            fileList = dir(fullfile(rawFolder, 'ai_*.bin'));

            if isempty(fileList)
                return
            end

            for iFile = 1:numel(fileList)
                if fileList(iFile).isdir
                    continue
                end

                fileName = char(string(fileList(iFile).name));

                if ~isempty(regexpi(fileName, '^ai_(\d+)\.bin$', 'once'))
                    tf = true;
                    return
                end
            end

        end

        function tf = currentDatSourceIsMultiCam(app)
            %CURRENTDATSOURCEISMULTICAM True when active DatImageSource reports MultiCam.

            tf = false;

            if ~app.hasData() || ~strcmpi(app.getSourceType(), 'dat')
                return
            end

            Info = app.getSourceInfo();

            if isempty(Info) || ~isstruct(Info)
                return
            end

            if iHasTrueMultiCamField(Info)
                tf = true;
                return
            end

            if isfield(Info, 'AcqInfoStream') && iHasTrueMultiCamField(Info.AcqInfoStream)
                tf = true;
                return
            end

            if isfield(Info, 'AcqInfos') && iHasTrueMultiCamField(Info.AcqInfos)
                tf = true;
                return
            end

            if isfield(Info, 'img_info') && iHasTrueMultiCamField(Info.img_info)
                tf = true;
                return
            end

            function tfLocal = iHasTrueMultiCamField(S)
                tfLocal = false;

                if ~isstruct(S) || ~isfield(S, 'MultiCam') || isempty(S.MultiCam)
                    return
                end

                val = S.MultiCam;

                try
                    if islogical(val) || isnumeric(val)
                        tfLocal = isscalar(val) && logical(val);
                        return
                    end

                    if ischar(val) || isstring(val)
                        txt = lower(strtrim(char(string(val))));
                        tfLocal = any(strcmp(txt, {'true', '1', 'yes', 'on'}));
                    end
                catch
                    tfLocal = false;
                end
            end

        end

        function tf = ensurePipelineManagerForCurrentData(app)
            %ENSUREPIPELINEMANAGERFORCURRENTDATA Create/reuse persistent PM context.
            %
            % DataViewer owns the PipelineManager object. The PM Tool receives this handle
            % and acts as a temporary editor/executor window.
            %
            % Important:
            %   Do NOT create/update the DataViewer FileSource here. The current input
            %   file is PM Tool UI state and is updated by PipelineManagerTool on every
            %   DataViewer-mode opening. This prevents a persistent pipeline from
            %   accumulating one FileSource node per execution turn.

            tf = false;

            if ~app.hasData() || isempty(app.CurrentFile) || ~isfile(app.CurrentFile)
                return
            end

            saveFolder = app.getCurrentSaveFolderForData();

            if isempty(saveFolder) || ~isfolder(saveFolder)
                return
            end

            rawFolder = app.resolveDataRawFolderForCurrentData(saveFolder);

            bNeedNew = isempty(app.PipelineManagerObj);

            if ~bNeedNew
                try
                    bNeedNew = ~isa(app.PipelineManagerObj, 'PipelineManager') || ...
                        isempty(app.PipelineManagerObj.SaveFolderList) || ...
                        ~strcmpi(app.PipelineManagerObj.SaveFolderList{1}, saveFolder);
                catch
                    bNeedNew = true;
                end
            end

            if bNeedNew
                app.PipelineManagerObj = PipelineManager(saveFolder, rawFolder);
                app.PipelineDataViewerFileSourceNodeID = NaN;

                app.PipelineManagerObj.setLeafOutputPolicy('viewerTemp');

                try
                    app.PipelineManagerObj.b_inputFromDataViewer = true;
                catch
                end
            else
                % Keep existing graph, but enforce DataViewer execution policy.
                try
                    app.PipelineManagerObj.setLeafOutputPolicy('viewerTemp');
                catch
                end

                % Keep folder context synchronized if the RawFolder was changed.
                try
                    if isprop(app.PipelineManagerObj, 'RawFolderList') && ...
                            ~isempty(app.PipelineManagerObj.RawFolderList)
                        app.PipelineManagerObj.RawFolderList{1} = rawFolder;
                    end
                catch
                end

                try
                    app.PipelineManagerObj.b_inputFromDataViewer = true;
                catch
                end
            end

            tf = true;
        end

        function semanticTypes = getDataViewerCurrentFileSemanticTypes(app, filePath) %#ok<INUSD>
            %GETDATAVIEWERCURRENTFILESEMANTICTYPES Return permissive source types.

            [~, ~, ext] = fileparts(char(string(filePath)));

            switch lower(ext)
                case '.dat'
                    semanticTypes = {'ImageTimeSeries', 'Image', 'UnknownDataType'};

                case '.umt'
                    semanticTypes = {'ProcessedData', 'ImageTimeSeries', 'Image', 'UnknownDataType'};

                otherwise
                    semanticTypes = {'UnknownDataType'};
            end
        end

        function resetPipelineContextForNewSaveFolder(app, previousFile, newFile)
            %RESETPIPELINECONTEXTFORNEWSAVEFOLDER Reset PM state when SaveFolder changes.

            previousFolder = '';
            newFolder = '';

            try
                if ~isempty(previousFile)
                    previousFolder = fileparts(char(string(previousFile)));
                end
            catch
                previousFolder = '';
            end

            try
                if ~isempty(newFile)
                    newFolder = fileparts(char(string(newFile)));
                end
            catch
                newFolder = '';
            end

            if isempty(newFolder)
                return
            end

            if ~isempty(previousFolder) && strcmpi(previousFolder, newFolder)
                return
            end

            % A different SaveFolder means a different dataset context. Do not reuse the
            % previous graph, DataViewer FileSource ID, temporary-file registry, or run summaries.
            app.PipelineManagerObj = [];
            app.PipelineManagerToolApp = [];
            app.PipelineDataViewerFileSourceNodeID = NaN;
            app.LastPipelineResult = [];
            app.PipelineRunHistory = cell(0,1);

            try
                app.PipelineTemporaryFiles = app.PipelineTemporaryFiles([], :);
            catch
                app.PipelineTemporaryFiles = table( ...
                    strings(0,1), ...
                    datetime.empty(0,1), ...
                    strings(0,1), ...
                    'VariableNames', {'FilePath','CreatedOn','Source'});
            end

            % RawFolder is SaveFolder-bound. Clear runtime cache unless it explicitly
            % belongs to the newly opened folder. It will be restored from DataParams.mat
            % during refreshViewerAfterLoad.
            try
                if ~strcmpi(char(string(app.DataRawFolderSaveFolder)), newFolder)
                    app.DataRawFolder = '';
                    app.DataRawFolderSaveFolder = '';
                end
            catch
                app.DataRawFolder = '';
                app.DataRawFolderSaveFolder = '';
            end

            app.RawFolderInvalidWarningSaveFolder = '';
        end

        function updateDataFolderContextLabel(app)
            %UPDATEDATAFOLDERCONTEXTLABEL Update compact SaveFolder | RawFolder label.

            if isempty(app.DataFolderContextLabel) || ~isvalid(app.DataFolderContextLabel)
                return
            end

            if ~app.hasData() || isempty(app.CurrentFile)
                app.DataFolderContextLabel.Text = 'SaveFolder: <none> | RawFolder: <none>';
                return
            end

            saveFolder = app.getCurrentSaveFolderForData();
            rawFolder = app.resolveDataRawFolderForCurrentData(saveFolder);

            app.DataFolderContextLabel.Text = sprintf('SaveFolder: %s | RawFolder: %s', ...
                app.shortenPathForLabel(saveFolder, 70), ...
                app.shortenPathForLabel(rawFolder, 70));
        end

        function txt = shortenPathForLabel(app, pathValue, maxChars) %#ok<INUSD>
            %SHORTENPATHFORLABEL Compact long paths for labels.

            if nargin < 3 || isempty(maxChars)
                maxChars = 70;
            end

            txt = char(string(pathValue));

            if isempty(txt)
                txt = '<none>';
                return
            end

            if numel(txt) <= maxChars
                return
            end

            keepTail = max(10, maxChars - 10);
            txt = ['...', txt(end-keepTail+1:end)];
        end

        function appendPipelineRunHistory(app, result)
            %APPENDPIPELINERUNHISTORY Store one execution result.

            if ~isstruct(result)
                return
            end

            try
                app.PipelineRunHistory{end+1,1} = result;
            catch
                app.PipelineRunHistory = {result};
            end
        end

        function updatePipelineTabState(app)
            %UPDATEPIPELINETABSTATE Refresh read-only Pipeline tab labels.

            if ~isempty(app.PipelineStatusLabel) && isvalid(app.PipelineStatusLabel)
                app.PipelineStatusLabel.Text = app.getPipelineStatusLabelText();
            end

            if ~isempty(app.LastrunoutputsLabel) && isvalid(app.LastrunoutputsLabel)
                app.LastrunoutputsLabel.Text = app.getPipelineOutputSummaryLabelText();
            end

            if ~isempty(app.CurrentpipelinefileLabel) && isvalid(app.CurrentpipelinefileLabel)
                app.CurrentpipelinefileLabel.Text = app.getCurrentPipelineFileLabelText();
            end

            if ~isempty(app.LatestsummaryLabel) && isvalid(app.LatestsummaryLabel)
                app.LatestsummaryLabel.Text = app.getLatestPipelineSummaryLabelText();
            end
        end

        function txt = getPipelineStatusLabelText(app)
            %GETPIPELINESTATUSLABELTEXT Build Pipeline tab status text.

            if isempty(app.PipelineManagerObj)
                pmText = 'No pipeline object';
            else
                pmText = 'Pipeline object ready';
            end

            if isempty(app.LastPipelineResult) || ~isstruct(app.LastPipelineResult) || ...
                    ~isfield(app.LastPipelineResult, 'status')
                txt = sprintf('Pipeline status: No run yet | %s', pmText);
                return
            end

            status = lower(char(string(app.LastPipelineResult.status)));

            switch status
                case {'completed', 'complete', 'passed', 'success'}
                    statusText = 'Passed';

                case {'completed_with_errors', 'partial', 'partially executed'}
                    statusText = 'Completed with errors';

                case {'failed', 'error'}
                    statusText = 'Failed';

                case {'cancelled', 'canceled'}
                    statusText = 'Cancelled';

                otherwise
                    statusText = char(string(app.LastPipelineResult.status));
            end

            txt = sprintf('Pipeline status: %s | %s', statusText, pmText);
        end

        function txt = getPipelineOutputSummaryLabelText(app)
            %GETPIPELINEOUTPUTSUMMARYLABELTEXT Build last-run output count text.

            result = app.LastPipelineResult;

            if isempty(result) || ~isstruct(result)
                txt = 'Last run outputs: none';
                return
            end

            [nCreated, nTemporary, nPermanent, nViewerCompatible] = ...
                app.getPipelineResultOutputCounts(result);

            txt = sprintf(['Last run outputs: Created %d | Temporary %d | ' ...
                'Permanent %d | Viewer-compatible %d'], ...
                nCreated, nTemporary, nPermanent, nViewerCompatible);
        end

        function [nCreated, nTemporary, nPermanent, nViewerCompatible] = getPipelineResultOutputCounts(app, result)
            %GETPIPELINERESULTOUTPUTCOUNTS Count output classes from a PM result.

            nCreated = 0;
            nTemporary = 0;
            nPermanent = 0;
            nViewerCompatible = 0;

            if ~isstruct(result)
                return
            end

            % Created files.
            try
                if isfield(result, 'createdFiles') && istable(result.createdFiles) && ...
                        ~isempty(result.createdFiles)
                    T = result.createdFiles;
                    nCreated = height(T);

                    if ismember('IsTemporary', T.Properties.VariableNames)
                        nTemporary = nnz(logical(T.IsTemporary));
                    elseif ismember('isTemporary', T.Properties.VariableNames)
                        nTemporary = nnz(logical(T.isTemporary));
                    else
                        files = app.getFilePathColumnFromTable(T);
                        for iFile = 1:numel(files)
                            [~, f, e] = fileparts(char(files(iFile)));
                            nTemporary = nTemporary + startsWith(string(f) + string(e), ...
                                "PMTMP_", 'IgnoreCase', true);
                        end
                    end

                    nPermanent = max(0, nCreated - nTemporary);
                end
            catch
                nCreated = 0;
                nTemporary = 0;
                nPermanent = 0;
            end

            % Viewer-compatible files.
            try
                compatibleFiles = app.getViewerCompatiblePipelineOutputs(result);
                nViewerCompatible = numel(compatibleFiles);
            catch
                nViewerCompatible = 0;
            end
        end

        function files = getFilePathColumnFromTable(app, T) %#ok<INUSD>
            %GETFILEPATHCOLUMNFROMTABLE Return FilePath/filePath column as string vector.

            files = strings(0,1);

            if ~istable(T) || isempty(T)
                return
            end

            if ismember('FilePath', T.Properties.VariableNames)
                files = string(T.FilePath(:));
            elseif ismember('filePath', T.Properties.VariableNames)
                files = string(T.filePath(:));
            end

            files = files(strlength(strtrim(files)) > 0);
        end

        function txt = getCurrentPipelineFileLabelText(app)
            %GETCURRENTPIPELINEFILELABELTEXT Build current-file Pipeline tab text.

            if isempty(app.CurrentFile)
                txt = 'Current pipeline file: none';
                return
            end

            [~, name, ext] = fileparts(app.CurrentFile);
            fileName = [name ext];

            if app.isCurrentFileRegisteredPipelineTemp()
                statusText = 'Temporary pipeline output';
            else
                statusText = 'Permanent / source file';
            end

            txt = sprintf('Current pipeline file: %s | %s', fileName, statusText);
        end

        function tf = isCurrentFileRegisteredPipelineTemp(app)
            %ISCURRENTFILEREGISTEREDPIPELINETEMP True if CurrentFile is registered temp.

            tf = false;

            if isempty(app.CurrentFile) || isempty(app.PipelineTemporaryFiles) || ...
                    height(app.PipelineTemporaryFiles) == 0
                return
            end

            try
                files = string(app.PipelineTemporaryFiles.FilePath(:));
                tf = any(strcmpi(files, string(app.CurrentFile)));
            catch
                tf = false;
            end
        end

        function txt = getLatestPipelineSummaryLabelText(app)
            %GETLATESTPIPELINESUMMARYLABELTEXT Build compact latest-run text.

            if isempty(app.LastPipelineResult) || ~isstruct(app.LastPipelineResult)
                txt = 'Latest summary: none';
                return
            end

            result = app.LastPipelineResult;

            timeText = '';
            durationText = '';

            try
                if isfield(result, 'finishedOn') && ~isempty(result.finishedOn)
                    timeText = char(string(result.finishedOn));
                elseif isfield(result, 'FinishedOn') && ~isempty(result.FinishedOn)
                    timeText = char(string(result.FinishedOn));
                elseif isfield(result, 'globalPipeLog') && istable(result.globalPipeLog) && ...
                        ismember('FinishedOn', result.globalPipeLog.Properties.VariableNames) && ...
                        ~isempty(result.globalPipeLog.FinishedOn)
                    timeText = char(string(result.globalPipeLog.FinishedOn(end)));
                end
            catch
                timeText = '';
            end

            try
                if isfield(result, 'duration_s') && ~isempty(result.duration_s)
                    durationText = sprintf('%.2f s', double(result.duration_s));
                elseif isfield(result, 'Duration_s') && ~isempty(result.Duration_s)
                    durationText = sprintf('%.2f s', double(result.Duration_s));
                elseif isfield(result, 'globalPipeLog') && istable(result.globalPipeLog) && ...
                        ismember('Duration_s', result.globalPipeLog.Properties.VariableNames) && ...
                        ~isempty(result.globalPipeLog.Duration_s)
                    durationText = sprintf('%.2f s', double(result.globalPipeLog.Duration_s(end)));
                end
            catch
                durationText = '';
            end

            nRuns = 0;
            try
                nRuns = numel(app.PipelineRunHistory);
            catch
                nRuns = 0;
            end

            if isempty(timeText) && isempty(durationText)
                txt = sprintf('Latest summary: Run history %d', nRuns);
            elseif isempty(durationText)
                txt = sprintf('Latest summary: %s | Run history %d', timeText, nRuns);
            elseif isempty(timeText)
                txt = sprintf('Latest summary: Duration %s | Run history %d', durationText, nRuns);
            else
                txt = sprintf('Latest summary: %s | Duration %s | Run history %d', ...
                    timeText, durationText, nRuns);
            end
        end

        function onPipelineManagerFinished(app, result)
            %ONPIPELINEMANAGERFINISHED Handle PipelineManagerTool execution result.

            app.LastPipelineResult = result;
            app.appendPipelineRunHistory(result);
            app.registerPipelineTemporaryFiles(result);
            app.updatePipelineTabState();

            if ~isstruct(result) || ~isfield(result, 'status')
                app.setStatusMessage('Pipeline finished, but returned an invalid result.');
                return
            end

            compatibleFiles = app.getViewerCompatiblePipelineOutputs(result);
            compatibleFiles = string(compatibleFiles(:));
            compatibleFiles = compatibleFiles(strlength(strtrim(compatibleFiles)) > 0);
            compatibleFiles = unique(compatibleFiles, 'stable');

            if ~strcmpi(string(result.status), "completed")
                app.setStatusMessage(sprintf('Pipeline finished with status: %s.', ...
                    char(string(result.status))));
                return
            end

            if isempty(compatibleFiles)
                uialert(app.UIFigure, ...
                    'The pipeline finished, but no viewer-compatible output was found.', ...
                    'Pipeline Finished', ...
                    'Icon', 'warning');
                return
            end

            if numel(compatibleFiles) == 1
                app.openPipelineOutputFile(compatibleFiles(1));
                return
            end

            selectedFile = app.promptPipelineOutputSelection(compatibleFiles);

            if strlength(selectedFile) > 0
                app.openPipelineOutputFile(selectedFile);
            end
        end

        function onPipelineManagerClosed(app, toolState)
            %ONPIPELINEMANAGERCLOSED Handle PM Tool close callback.
            %
            % The PipelineManager object intentionally survives. Only the transient PM Tool
            % app handle is cleared. The DataViewer FileSource node id is stored so the
            % next PM Tool opening updates the same node to the current DataViewer file.

            app.PipelineManagerToolApp = [];

            if nargin >= 2 && isstruct(toolState)
                try
                    if isfield(toolState, 'dataViewerFileSourceNodeID') && ...
                            ~isempty(toolState.dataViewerFileSourceNodeID) && ...
                            isfinite(double(toolState.dataViewerFileSourceNodeID))
                        app.PipelineDataViewerFileSourceNodeID = double(toolState.dataViewerFileSourceNodeID);
                    elseif isfield(toolState, 'DataViewerFileSourceNodeID') && ...
                            ~isempty(toolState.DataViewerFileSourceNodeID) && ...
                            isfinite(double(toolState.DataViewerFileSourceNodeID))
                        app.PipelineDataViewerFileSourceNodeID = double(toolState.DataViewerFileSourceNodeID);
                    end
                catch
                end

                try
                    if isfield(toolState, 'lastExecutionResult') && ~isempty(toolState.lastExecutionResult)
                        app.LastPipelineResult = toolState.lastExecutionResult;
                        app.appendPipelineRunHistory(toolState.lastExecutionResult);
                    elseif isfield(toolState, 'LastExecutionResult') && ~isempty(toolState.LastExecutionResult)
                        app.LastPipelineResult = toolState.LastExecutionResult;
                        app.appendPipelineRunHistory(toolState.LastExecutionResult);
                    end
                catch
                end
            end

            try
                app.UIFigure.Visible = 'on';
                figure(app.UIFigure);
            catch
            end

            try
                app.setInteractionMode('idle');
            catch
            end

            app.updatePipelineTabState();
            app.updateDataFolderContextLabel();
            app.setStatusMessage('Pipeline Manager Tool closed. Pipeline object preserved.');
        end

        function openPipelineOutputFile(app, filePath)
            %OPENPIPELINEOUTPUTFILE Open one compatible PM output in DataViewer.
            %
            % This method intentionally does not clean temporary files so pipeline outputs
            % can be chained: source -> PMTMP_1 -> PMTMP_2 -> ...

            filePath = char(string(filePath));

            if isempty(strtrim(filePath)) || ~isfile(filePath)
                error('DataViewer:PipelineOutputNotFound', ...
                    'Pipeline output file not found: "%s".', filePath);
            end

            bLoaded = app.loadDataSource(filePath);

            if bLoaded
                app.refreshViewerAfterLoad();
                app.updatePipelineTabState();
                app.updateDataFolderContextLabel();
            end
        end

        function compatibleFiles = getViewerCompatiblePipelineOutputs(app, result) %#ok<INUSD>
            %GETVIEWERCOMPATIBLEPIPELINEOUTPUTS Return files DataViewer can open.
            %
            % Priority:
            %   1) result.viewerCompatibleFinalOutputs
            %   2) result.outputManifest rows marked viewer-compatible
            %   3) result.createdFiles fallback filtered by extension/existence

            compatibleFiles = strings(0,1);

            if ~isstruct(result)
                return
            end

            % ---------------------------------------------------------------------
            % Preferred backend report: viewer-compatible final outputs.
            % ---------------------------------------------------------------------
            try
                if isfield(result, 'viewerCompatibleFinalOutputs') && ...
                        istable(result.viewerCompatibleFinalOutputs) && ...
                        ~isempty(result.viewerCompatibleFinalOutputs)

                    T = result.viewerCompatibleFinalOutputs;

                    if ismember('FilePath', T.Properties.VariableNames)
                        compatibleFiles = string(T.FilePath(:));
                    elseif ismember('filePath', T.Properties.VariableNames)
                        compatibleFiles = string(T.filePath(:));
                    end
                end
            catch
                compatibleFiles = strings(0,1);
            end

            % ---------------------------------------------------------------------
            % Fallback: output manifest.
            % ---------------------------------------------------------------------
            if isempty(compatibleFiles)
                try
                    if isfield(result, 'outputManifest') && istable(result.outputManifest) && ...
                            ~isempty(result.outputManifest)

                        T = result.outputManifest;

                        if all(ismember({'FilePath','FileExists','IsViewerCompatible'}, T.Properties.VariableNames))
                            keep = logical(T.FileExists) & logical(T.IsViewerCompatible);
                            compatibleFiles = string(T.FilePath(keep));
                        elseif all(ismember({'filePath','fileExists','isViewerCompatible'}, T.Properties.VariableNames))
                            keep = logical(T.fileExists) & logical(T.isViewerCompatible);
                            compatibleFiles = string(T.filePath(keep));
                        end
                    end
                catch
                    compatibleFiles = strings(0,1);
                end
            end

            % ---------------------------------------------------------------------
            % Fallback: created files. Only extension/existence filtered here.
            % ---------------------------------------------------------------------
            if isempty(compatibleFiles)
                try
                    if isfield(result, 'createdFiles') && istable(result.createdFiles) && ...
                            ~isempty(result.createdFiles)

                        T = result.createdFiles;

                        if ismember('FilePath', T.Properties.VariableNames)
                            files = string(T.FilePath(:));
                        elseif ismember('filePath', T.Properties.VariableNames)
                            files = string(T.filePath(:));
                        else
                            files = strings(0,1);
                        end

                        files = files(strlength(strtrim(files)) > 0);

                        keep = false(size(files));
                        for iFile = 1:numel(files)
                            [~, ~, ext] = fileparts(char(files(iFile)));
                            keep(iFile) = isfile(files(iFile)) && ...
                                any(strcmpi(ext, {'.dat', '.umt'}));
                        end

                        compatibleFiles = files(keep);
                    end
                catch
                    compatibleFiles = strings(0,1);
                end
            end

            compatibleFiles = compatibleFiles(:);
            compatibleFiles = compatibleFiles(strlength(strtrim(compatibleFiles)) > 0);
            compatibleFiles = compatibleFiles(arrayfun(@(x) isfile(char(x)), compatibleFiles));
            compatibleFiles = unique(compatibleFiles, 'stable');
        end

        function tf = isViewerCompatiblePipelineFile(app, filePath) %#ok<INUSD>
            %ISVIEWERCOMPATIBLEPIPELINEFILE True for files DataViewer can open.

            tf = false;
            filePath = char(string(filePath));

            if ~isfile(filePath)
                return
            end

            [~, ~, ext] = fileparts(filePath);

            switch lower(ext)
                case '.dat'
                    tf = true;

                case '.umt'
                    try
                        entrySummary = UMTImageSource.inspectEntries(filePath);
                        tf = istable(entrySummary) && height(entrySummary) >= 1;
                    catch
                        tf = false;
                    end
            end
        end

        function selectedFile = promptPipelineOutputSelection(app, outputTable)
            %PROMPTPIPELINEOUTPUTSELECTION Ask user to choose one PM output.

            selectedFile = "";

            if isempty(outputTable) || height(outputTable) == 0
                return
            end

            selectedRow = 1;

            dlg = uifigure( ...
                'Name', 'Select pipeline output', ...
                'WindowStyle', 'modal', ...
                'Position', [100 100 760 360], ...
                'Visible', 'off', ...
                'CloseRequestFcn', @onCancel);

            grid = uigridlayout(dlg);
            grid.RowHeight = {25, '1x', 35};
            grid.ColumnWidth = {'1x'};
            grid.Padding = [12 12 12 12];

            titleLabel = uilabel(grid);
            titleLabel.Text = 'Select the pipeline output to open in DataViewer:';
            titleLabel.FontWeight = 'bold';
            titleLabel.Layout.Row = 1;
            titleLabel.Layout.Column = 1;

            displayTable = outputTable;
            try
                displayTable.FilePath = string(displayTable.FilePath);
            catch
            end

            tbl = uitable(grid);
            tbl.Data = displayTable;
            tbl.ColumnEditable = false(1, width(displayTable));
            tbl.RowName = {};
            tbl.Layout.Row = 2;
            tbl.Layout.Column = 1;
            tbl.SelectionChangedFcn = @onSelectionChanged;

            try
                if isprop(tbl, 'SelectionType')
                    tbl.SelectionType = 'row';
                end
                if isprop(tbl, 'Multiselect')
                    tbl.Multiselect = 'off';
                end
                tbl.Selection = 1;
            catch
            end

            buttonGrid = uigridlayout(grid);
            buttonGrid.RowHeight = {'1x'};
            buttonGrid.ColumnWidth = {'1x', '1x'};
            buttonGrid.Padding = [0 0 0 0];
            buttonGrid.Layout.Row = 3;
            buttonGrid.Layout.Column = 1;

            okButton = uibutton(buttonGrid, 'push');
            okButton.Text = 'Open';
            okButton.Layout.Row = 1;
            okButton.Layout.Column = 1;
            okButton.ButtonPushedFcn = @onOK;

            cancelButton = uibutton(buttonGrid, 'push');
            cancelButton.Text = 'Cancel';
            cancelButton.Layout.Row = 1;
            cancelButton.Layout.Column = 2;
            cancelButton.ButtonPushedFcn = @onCancel;

            if exist('placeAppInsideCaller', 'file') == 2
                placeAppInsideCaller(app, dlg, 'center');
            end
            dlg.Visible = 'on';

            uiwait(dlg);

            if isvalid(dlg)
                delete(dlg);
            end

            function onSelectionChanged(~, event)
                try
                    if isprop(event, 'Selection') && ~isempty(event.Selection)
                        selection = event.Selection;
                    else
                        selection = tbl.Selection;
                    end

                    if isvector(selection)
                        selectedRow = round(double(selection(1)));
                    else
                        selectedRow = round(double(selection(1,1)));
                    end
                catch
                    selectedRow = 1;
                end
            end

            function onOK(~, ~)
                if isempty(selectedRow) || selectedRow < 1 || selectedRow > height(outputTable)
                    return
                end
                selectedFile = string(outputTable.FilePath(selectedRow));
                uiresume(dlg);
            end

            function onCancel(~, ~)
                selectedFile = "";
                uiresume(dlg);
            end
        end

        function registerPipelineTemporaryFiles(app, result)
            %REGISTERPIPELINETEMPORARYFILES Add PM temporary outputs to registry.

            if ~isstruct(result) || ~isfield(result, 'createdFiles') || ...
                    ~istable(result.createdFiles) || isempty(result.createdFiles)
                app.updatePipelineTabState();
                return
            end

            T = result.createdFiles;
            filePaths = app.getFilePathColumnFromTable(T);

            if isempty(filePaths)
                app.updatePipelineTabState();
                return
            end

            isTemp = false(numel(filePaths), 1);

            if ismember('IsTemporary', T.Properties.VariableNames)
                try
                    isTemp = logical(T.IsTemporary(:));
                catch
                    isTemp = false(numel(filePaths), 1);
                end
            elseif ismember('isTemporary', T.Properties.VariableNames)
                try
                    isTemp = logical(T.isTemporary(:));
                catch
                    isTemp = false(numel(filePaths), 1);
                end
            else
                for iFile = 1:numel(filePaths)
                    [~, f, e] = fileparts(filePaths(iFile));
                    isTemp(iFile) = startsWith(string(f) + string(e), "PMTMP_", ...
                        'IgnoreCase', true);
                end
            end

            filePaths = filePaths(isTemp & arrayfun(@(x) isfile(char(x)), filePaths));
            filePaths = unique(filePaths, 'stable');

            if isempty(filePaths)
                app.updatePipelineTabState();
                return
            end

            newRows = table( ...
                filePaths(:), ...
                repmat(datetime('now'), numel(filePaths), 1), ...
                repmat("PipelineManagerTool", numel(filePaths), 1), ...
                'VariableNames', {'FilePath','CreatedOn','Source'});

            app.PipelineTemporaryFiles = unique([app.PipelineTemporaryFiles; newRows], 'rows');
            app.updatePipelineTabState();
        end

        function canContinue = cleanupPipelineTemporaryFiles(app, varargin)
            %CLEANUPPIPELINETEMPORARYFILES Delete registered PM temporary files.

            p = inputParser;
            addParameter(p, 'PreserveFile', '', @(x)ischar(x)||isstring(x));
            addParameter(p, 'AskForCurrentTemporaryFile', true, @(x)islogical(x)&&isscalar(x));
            parse(p, varargin{:});

            preserveFile = string(p.Results.PreserveFile);
            askForCurrentTemporaryFile = logical(p.Results.AskForCurrentTemporaryFile);
            canContinue = true;

            if isempty(app.PipelineTemporaryFiles) || height(app.PipelineTemporaryFiles) == 0
                app.updatePipelineTabState();
                return
            end

            files = string(app.PipelineTemporaryFiles.FilePath(:));
            files = unique(files(strlength(files) > 0), 'stable');

            if isempty(files)
                app.PipelineTemporaryFiles = app.PipelineTemporaryFiles([], :);
                app.updatePipelineTabState();
                return
            end

            currentFile = "";
            try
                currentFile = string(app.CurrentFile);
            catch
            end

            currentIsRegisteredTemp = strlength(currentFile) > 0 && ...
                any(strcmpi(files, currentFile)) && isfile(currentFile);

            if currentIsRegisteredTemp

                if askForCurrentTemporaryFile && ~strcmpi(currentFile, preserveFile)

                    action = app.promptCurrentTemporaryFileCleanup(currentFile);

                    switch action
                        case 'saveas'
                            savedPath = app.saveTemporaryPipelineFileAs(currentFile);

                            if strlength(savedPath) == 0
                                canContinue = false;
                                app.updatePipelineTabState();
                                return
                            end

                            if ~strcmpi(currentFile, preserveFile) && isfile(currentFile)
                                try
                                    delete(currentFile);
                                catch ME
                                    warning('DataViewer:PipelineTempCleanupFailed', ...
                                        'Could not delete temporary pipeline file "%s": %s', ...
                                        currentFile, ME.message);
                                end
                            end

                        case 'delete'
                            % Continue into normal deletion loop.

                        otherwise
                            canContinue = false;
                            app.updatePipelineTabState();
                            return
                    end

                elseif ~askForCurrentTemporaryFile
                    % Important for pipeline chaining:
                    % do not silently delete the current pipeline temp file.
                    if strlength(preserveFile) == 0
                        preserveFile = currentFile;
                    end
                end
            end

            for iFile = 1:numel(files)
                filePath = files(iFile);

                if strlength(preserveFile) > 0 && strcmpi(filePath, preserveFile)
                    continue
                end

                if isfile(filePath)
                    try
                        delete(filePath);
                    catch ME
                        warning('DataViewer:PipelineTempCleanupFailed', ...
                            'Could not delete temporary pipeline file "%s": %s', ...
                            filePath, ME.message);
                    end
                end
            end

            keepMask = false(height(app.PipelineTemporaryFiles), 1);

            for iRow = 1:height(app.PipelineTemporaryFiles)
                f = string(app.PipelineTemporaryFiles.FilePath(iRow));

                if strlength(preserveFile) > 0 && strcmpi(f, preserveFile) && isfile(f)
                    keepMask(iRow) = true;
                elseif isfile(f)
                    keepMask(iRow) = true;
                end
            end

            app.PipelineTemporaryFiles = app.PipelineTemporaryFiles(keepMask, :);
            app.updatePipelineTabState();
        end

        function action = promptCurrentTemporaryFileCleanup(app, tempFile)
            %PROMPTCURRENTTEMPORARYFILECLEANUP Ask what to do with current PM temp file.
            %
            % Output:
            %   'saveas'  - copy current temp file to a permanent file, then continue cleanup
            %   'delete'  - delete current temp file and continue cleanup
            %   'cancel'  - abort the caller action

            tempFile = string(tempFile);
            [~, tempName, tempExt] = fileparts(char(tempFile));
            tempDisplay = [tempName tempExt];

            nOtherTmp = 0;
            try
                if ~isempty(app.PipelineTemporaryFiles) && height(app.PipelineTemporaryFiles) > 0
                    tmpFiles = string(app.PipelineTemporaryFiles.FilePath(:));
                    tmpFiles = tmpFiles(strlength(strtrim(tmpFiles)) > 0);
                    tmpFiles = unique(tmpFiles, 'stable');
                    tmpFiles = tmpFiles(arrayfun(@(x) isfile(char(x)), tmpFiles));
                    nOtherTmp = nnz(~strcmpi(tmpFiles, tempFile));
                end
            catch
                nOtherTmp = 0;
            end

            if nOtherTmp == 0
                msg = sprintf([ ...
                    'The current file is a temporary PipelineManager output:\n\n' ...
                    '%s\n\n' ...
                    'Save it as a permanent file before continuing?'], ...
                    tempDisplay);
            elseif nOtherTmp == 1
                msg = sprintf([ ...
                    'The current file is a temporary PipelineManager output:\n\n' ...
                    '%s\n\n' ...
                    'Save it as a permanent file before continuing?\n\n' ...
                    'Note: 1 other temporary PipelineManager file from this chain will be deleted.'], ...
                    tempDisplay);
            else
                msg = sprintf([ ...
                    'The current file is a temporary PipelineManager output:\n\n' ...
                    '%s\n\n' ...
                    'Save it as a permanent file before continuing?\n\n' ...
                    'Note: %d other temporary PipelineManager files from this chain will be deleted.'], ...
                    tempDisplay, nOtherTmp);
            end

            choice = uiconfirm(app.UIFigure, ...
                msg, ...
                'Temporary Pipeline Output', ...
                'Options', {'Save As...', 'Delete Temporary', 'Cancel'}, ...
                'DefaultOption', 1, ...
                'CancelOption', 3, ...
                'Icon', 'warning');

            switch choice
                case 'Save As...'
                    action = 'saveas';

                case 'Delete Temporary'
                    action = 'delete';

                otherwise
                    action = 'cancel';
            end
        end

        function savedPath = saveTemporaryPipelineFileAs(app, sourceFile)
            %SAVETEMPORARYPIPELINEFILEAS Copy current temp output to permanent file.

            savedPath = "";
            sourceFile = char(string(sourceFile));

            if ~isfile(sourceFile)
                return
            end

            [srcFolder, srcName, srcExt] = fileparts(sourceFile);

            [fileName, folderName] = uiputfile( ...
                {'*.dat;*.umt', 'Viewer data (*.dat, *.umt)'}, ...
                'Save temporary pipeline output as', ...
                fullfile(srcFolder, [srcName srcExt]));

            if isequal(fileName, 0)
                return
            end

            targetFile = fullfile(folderName, fileName);

            if strcmpi(sourceFile, targetFile)
                uialert(app.UIFigure, ...
                    'Choose a different filename or folder for the permanent copy.', ...
                    'Invalid Save As Target', ...
                    'Icon', 'warning');
                return
            end

            copyfile(sourceFile, targetFile);
            savedPath = string(targetFile);
        end

    end

    methods (Access = public)

        function startAllenAtlasROIPlacement(app, atlasResult)
            %STARTALLENATLASROIPLACEMENT Start temporary atlas ROI fit/placement.
            %
            %   Called by AllenAtlasROICreator_exported after OK. The selected atlas
            %   polygons are projected into DataViewer image coordinates, drawn as a
            %   temporary overlay, and controlled by one unbounded drawrectangle.
            %
            %   Double-click or Enter commits polygon ROIs.
            %   Escape or context menu cancels.

            if ~app.hasData()
                app.setStatusMessage('Load image data before placing Allen atlas ROIs.');
                return
            end

            if ~isstruct(atlasResult) || ~isfield(atlasResult, 'selectedRegions') || ...
                    isempty(atlasResult.selectedRegions)
                app.setStatusMessage('No Allen atlas regions were selected.');
                return
            end

            app.cancelAllenAtlasROIPlacement(false);

            previousKeyFcn = app.UIFigure.WindowKeyPressFcn;

            [state, ok, msg] = app.buildAllenAtlasROIPlacementState(atlasResult, previousKeyFcn);

            if ~ok
                app.setStatusMessage(msg);
                return
            end

            app.AtlasROIPlacementState = state;

            hold(app.ImageAxes, 'on');

            % Preview selected atlas regions.
            nRegions = numel(state.regionVertices0);
            app.AtlasROIPlacementRegionHandles = gobjects(0);

            for iRegion = 1:nRegions
                verticesXY = state.regionVertices0{iRegion};

                if isempty(verticesXY)
                    continue
                end

                h = patch(app.ImageAxes, ...
                    'XData', verticesXY(:, 1), ...
                    'YData', verticesXY(:, 2), ...
                    'FaceColor', state.regionColor{iRegion}, ...
                    'FaceAlpha', 0.35, ...
                    'EdgeColor', [0 0 0], ...
                    'LineWidth', 0.75, ...
                    'HitTest', 'off', ...
                    'PickableParts', 'none');

                app.AtlasROIPlacementRegionHandles(end+1) = h; %#ok<AGROW>
            end

            % Preview full atlas area boundary mask in black.
            app.AtlasROIPlacementBoundaryHandles = gobjects(0);

            for iB = 1:numel(state.boundaryVertices0)
                xy = state.boundaryVertices0{iB};

                if isempty(xy)
                    continue
                end

                h = line(app.ImageAxes, ...
                    xy(:, 1), xy(:, 2), ...
                    'Color', [0 0 0], ...
                    'LineWidth', 0.75, ...
                    'HitTest', 'off', ...
                    'PickableParts', 'none');

                app.AtlasROIPlacementBoundaryHandles(end+1) = h; %#ok<AGROW>
            end

            % Bregma/Lambda markers are shown only in DataViewer placement.
            app.AtlasROIPlacementBregmaHandle = line(app.ImageAxes, ...
                state.bregma0(1), state.bregma0(2), ...
                'Marker', 'o', ...
                'MarkerFaceColor', [1 0 0], ...
                'MarkerEdgeColor', [1 0 0], ...
                'LineStyle', 'none', ...
                'MarkerSize', 7, ...
                'HitTest', 'off', ...
                'PickableParts', 'none');

            app.AtlasROIPlacementLambdaHandle = line(app.ImageAxes, ...
                state.lambda0(1), state.lambda0(2), ...
                'Marker', 'o', ...
                'MarkerFaceColor', [0 1 0], ...
                'MarkerEdgeColor', [0 1 0], ...
                'LineStyle', 'none', ...
                'MarkerSize', 7, ...
                'HitTest', 'off', ...
                'PickableParts', 'none');

            bbox = state.originalBBox;

            app.AtlasROIPlacementBoxHandle = drawrectangle(app.ImageAxes, ...
                'Position', bbox, ...
                'Color', [1 1 1], ...
                'FaceAlpha', 0, ...
                'LineWidth', 1.5, ...
                'Rotatable', true, ...
                'DrawingArea', 'unlimited', ...
                'InteractionsAllowed', 'all');

            app.setROIObjectPropertyIfAvailable(app.AtlasROIPlacementBoxHandle, 'LineStyle', '--');
            app.updateAllenAtlasROIPlacementBoxSizeLabel(app.AtlasROIPlacementBoxHandle);

            listeners = {};
            listeners{end+1} = addlistener(app.AtlasROIPlacementBoxHandle, 'MovingROI', ...
                @(src, evt) app.onAllenAtlasROIPlacementBoxMoving(src, evt));
            listeners{end+1} = addlistener(app.AtlasROIPlacementBoxHandle, 'ROIMoved', ...
                @(src, evt) app.onAllenAtlasROIPlacementBoxMoved(src, evt));
            listeners{end+1} = addlistener(app.AtlasROIPlacementBoxHandle, 'ROIClicked', ...
                @(src, evt) app.onAllenAtlasROIPlacementBoxClicked(src, evt));

            try
                listeners{end+1} = addlistener(app.AtlasROIPlacementBoxHandle, 'DeletingROI', ...
                    @(src, evt) app.onAllenAtlasROIPlacementBoxDeleting(src, evt));
            catch
            end

            app.AtlasROIPlacementState.listeners = listeners;

            app.createAllenAtlasROIPlacementContextMenu();

            app.UIFigure.WindowKeyPressFcn = @(src, evt) app.onAllenAtlasROIPlacementKeyPress(src, evt);
            app.setInteractionMode('editingROI');

            app.stackAllenAtlasROIPlacementGraphics();
            app.setStatusMessage( ...
                'Allen atlas ROI placement. Move/resize/rotate the box; double-click or press Enter to create ROIs; Escape to cancel.');

        end

        function cleanupObj = suppressExpectedPolyshapeWarnings(app)
            %SUPPRESSEXPECTEDPOLYSHAPEWARNINGS Temporarily suppress noisy polyshape warnings.
            %
            %   Atlas boundaries can contain duplicate vertices, tiny fragments, or
            %   self-intersections. MATLAB's polyshape repair warnings are expected during
            %   mask-to-polygon conversion and do not need to be shown repeatedly to users.

            warnState = warning();

            warning('off', 'MATLAB:polyshape:repairedBySimplify');
            warning('off', 'MATLAB:polyshape:boolOperationFailed');
            warning('off', 'MATLAB:polyshape:boundary3Points');
            warning('off', 'MATLAB:polyshape:emptyShape');

            % Fallback for release/version-dependent warning identifiers. Scope is still
            % local because the warning state is restored by onCleanup.
            warning('off', 'backtrace');

            cleanupObj = onCleanup(@() warning(warnState));

        end

        function n = getLoadedROICountForExternalTool(app)
            %GETLOADEDROICOUNTFOREXTERNALTOOL Return number of loaded in-memory ROIs.
            %
            %   This method is intentionally public so standalone utility apps can
            %   warn the user without directly accessing DataViewer private ROIList.

            if isempty(app.ROIList)
                n = 0;
            else
                n = numel(app.ROIList);
            end
        end

        function refreshAfterImageAlignment(app, targetFolder)
            %REFRESHAFTERIMAGEALIGNMENT Reload viewer state after folder registration.
            %
            %   refreshAfterImageAlignment(app, targetFolder) is called by
            %   ImageAlignmentTool after applyImageAlignmentToFolder succeeds.
            %
            %   Expected behavior:
            %       1) clear stale in-memory ROIs
            %       2) reload the current .dat file, or first .dat in targetFolder
            %       3) reload DataParams.mat through loadDataSource()
            %       4) refresh the viewer and GUI state
            %
            %   Notes:
            %       - loadDataSource() already reloads folder-global DataParams.mat.
            %       - refreshViewerAfterLoad() already refreshes image, frame controls,
            %         temporal plot, cache, event patches, ROI traces, and GUI state.
            %       - ROIs are cleared rather than transformed in memory because the
            %         backend transforms saved .roi files. Keeping old in-memory ROI
            %         geometry would be misleading.

            if nargin < 2 || isempty(targetFolder)
                targetFolder = app.getCurrentDataFolder();
            end

            targetFolder = char(string(targetFolder));

            if ~isfolder(targetFolder)
                error('DataViewer:InvalidAlignmentTargetFolder', ...
                    'Image alignment target folder does not exist: %s', targetFolder);
            end

            if ~app.hasData()
                app.setStatusMessage('Image registration completed, but no data are loaded.');
                return
            end

            currentFile = app.CurrentFile;
            reloadFile = app.resolveImageAlignmentReloadFile(targetFolder, currentFile);

            hadROIs = ~isempty(app.ROIList);
            nROIs = numel(app.ROIList);

            if hadROIs
                app.clearLoadedROIsAfterImageAlignment();
            end

            previousMode = app.ViewMode;

            try
                if strcmpi(previousMode, 'event')
                    app.ViewMode = 'normal';
                    app.syncSwitchToViewMode();
                end
            catch
                app.ViewMode = 'normal';
            end

            bLoaded = app.loadDataSource(reloadFile);

            if bLoaded
                app.refreshViewerAfterLoad();

                if hadROIs
                    app.setStatusMessage(sprintf( ...
                        ['Image registration applied. DataParams and data were reloaded. ' ...
                        '%d stale in-memory ROI(s) were cleared; reload transformed .roi files if needed.'], ...
                        nROIs));
                else
                    app.setStatusMessage('Image registration applied. DataParams and data were reloaded.');
                end
            else
                app.setStatusMessage('Image registration applied, but DataViewer reload was cancelled.');
            end
        end

        function reloadFile = resolveImageAlignmentReloadFile(app, targetFolder, currentFile)
            %RESOLVEIMAGEALIGNMENTRELOADFILE Choose file to reload after registration.
            %
            %   Prefer the current .dat file if it is still inside targetFolder.
            %   Otherwise reload the first .dat file in the registered folder.

            reloadFile = '';

            if nargin < 3
                currentFile = '';
            end

            currentFile = char(string(currentFile));

            if ~isempty(currentFile) && isfile(currentFile)
                [currentFolder, ~, ext] = fileparts(currentFile);

                if strcmpi(ext, '.dat') && strcmpi(currentFolder, targetFolder)
                    reloadFile = currentFile;
                    return
                end
            end

            datFiles = dir(fullfile(targetFolder, '*.dat'));

            if isempty(datFiles)
                error('DataViewer:NoDatFileAfterImageAlignment', ...
                    ['Image alignment completed, but no .dat file was found in the ' ...
                    'target folder for DataViewer reload: %s'], targetFolder);
            end

            reloadFile = fullfile(datFiles(1).folder, datFiles(1).name);
        end

        function clearLoadedROIsAfterImageAlignment(app)
            %CLEARLOADEDROISAFTERIMAGEALIGNMENT Clear stale in-memory ROI state.
            %
            %   Saved .roi files are transformed by the backend. In-memory ROIs in
            %   DataViewer are not transformed, so clear them after successful folder
            %   registration to avoid showing stale geometry.

            if isempty(app.ROIList)
                app.refreshROITable();
                app.updateGUIEnabledState();
                return
            end

            try
                app.deleteGroupEditRuntimeGraphics();
            catch
            end

            for iROI = 1:numel(app.ROIList)
                try
                    app.deleteROIGraphicsByIndex(iROI);
                catch
                end
            end

            app.ROIList = struct([]);
            app.SelectedROIID = NaN;
            app.ROISelectionOrder = [];

            try
                app.GroupROIEditState = struct();
                app.GroupEditBoxHandle = gobjects(1);
                app.GroupEditPreviewHandles = gobjects(0);
            catch
            end

            try
                app.refreshROITable();
            catch
            end

            try
                app.refreshDeleteROIContextMenuState();
            catch
            end

            try
                app.refreshTemporalPlot();
            catch
            end

            app.updateGUIEnabledState();
        end

        function setRawFolderFromChildTool(app, rawFolder, saveFolder, setBy)
            %SETRAWFOLDERFROMCHILDTOOL Update RawFolder from a child utility app.
            %
            %   Child tools such as DataViewer_EventsManager should call this when the
            %   user sets or changes RawFolder inside the child tool.

            if nargin < 3 || isempty(saveFolder)
                saveFolder = app.getCurrentSaveFolderForData();
            end

            if nargin < 4 || isempty(setBy)
                setBy = 'DataViewer_childTool';
            end

            saveFolder = char(string(saveFolder));

            if isempty(saveFolder) || ~isfolder(saveFolder)
                return
            end

            currentSaveFolder = app.getCurrentSaveFolderForData();

            if ~isempty(currentSaveFolder) && strcmpi(currentSaveFolder, saveFolder)
                app.setRawFolderForCurrentSaveFolder(rawFolder, setBy, true);
                return
            end

            % Fallback for non-current SaveFolder updates. This persists the value but
            % does not alter the current viewer runtime context.
            app.persistRawFolderToDataParams(saveFolder, rawFolder, setBy);
        end

    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, filePath)
            app.StatusLabel.Tag = 'AlwaysEnabled';

            app.initializeViewerGraphics();

            app.populateColormapDropdown();
            app.applySelectedColormap();
            app.updateMovieSpeedLabel();

            app.updateImageStatusLabel();
            app.updateDataFolderContextLabel();
            app.updatePipelineTabState();

            % Configure the ROI table before the first GUI state evaluation so
            % updateGUIEnabledState sees a valid table setup.
            app.setupROITable();
            app.createLogicalMaskButtonContextMenu();
            app.createSetReferenceButtonContextMenu();
            app.setInteractionMode('idle');

            % Load input data
            if exist('filePath','var')
                [Path,fileName,ext] = fileparts(char(filePath));
                if isempty(Path)
                    Path = pwd;
                end
                % Rebuild full file path
                filePath = fullfile(Path,[fileName,ext]);
                bLoaded = app.loadDataSource(filePath);

                if bLoaded
                    app.refreshViewerAfterLoad();
                end
            end

        end

        % Menu selected function: OpenMenu
        function OpenMenuSelected(app, event)
            app.openDataFromDialog();
        end

        % Value changed function: Slider
        function SliderValueChanged(app, event)
            %SLIDERVALUECHANGED Commit frame change after slider release.

            if ~app.hasData()
                return
            end

            app.setCurrentFrame(app.Slider.Value);


        end

        % Button pushed function: PreviousFrameButton
        function PreviousFrameButtonPushed(app, event)
            app.goToPreviousFrame();
        end

        % Button pushed function: NextFrameButton
        function NextFrameButtonPushed(app, event)
            app.goToNextFrame();
        end

        % Value changing function: Slider
        function SliderValueChanging(app, event)
            %SLIDERVALUECHANGING Preview target frame while dragging the slider.
            %
            %   This callback is intentionally lightweight. It updates only the title
            %   and timebar, without loading the image frame or updating the trace.

            if ~app.hasData()
                return
            end

            sz = app.getDataSize();
            previewFrame = min(max(round(event.Value), 1), sz(3));

            [~, name, ext] = fileparts(app.CurrentFile);
            title(app.ImageAxes, sprintf('%s%s | moving to frame %d', ...
                name, ext, previewFrame));

            app.refreshTimeBar(previewFrame);

        end

        % Button pushed function: PlayMovieButton
        function PlayMovieButtonPushed(app, event)
            %PLAYMOVIEBUTTONPUSHED Toggle movie playback.

            app.toggleMoviePlayback();


        end

        % Value changed function: MovieSpeedDropDown
        function MovieSpeedDropDownValueChanged(app, event)
            %MOVIESPEEDDROPDOWNVALUECHANGED Update target playback speed.

            app.updateMovieSpeedLabel();

            if strcmp(app.InteractionMode, 'playingMovie')
                % Restart timing from the current frame so speed changes are applied
                % immediately without using the old timing reference.
                if ~isempty(app.MovieTimer) && isvalid(app.MovieTimer) && ...
                        strcmpi(app.MovieTimer.Running, 'on')
                    stop(app.MovieTimer);
                end

                app.setInteractionMode('idle');
                app.toggleMoviePlayback();
            end
        end

        % Close request function: UIFigure
        function UIFigureCloseRequest(app, event)
            %UIFIGURECLOSEREQUEST Safely close DataViewer.

            canClose = true;

            try
                canClose = app.cleanupPipelineTemporaryFiles( ...
                    'AskForCurrentTemporaryFile', true);
            catch ME
                warning('DataViewer:PipelineTempCleanupOnCloseFailed', ...
                    'Pipeline temporary-file cleanup failed before close.\n%s', ME.message);
                canClose = true;
            end

            if ~canClose
                return
            end

            app.cleanupAppResources();
            delete(app);

        end

        % Value changed function: ColormapDropDown
        function ColormapDropDownValueChanged(app, event)
            app.applySelectedColormap();
        end

        % Value changed function: InvertCheckBox
        function InvertCheckBoxValueChanged(app, event)
            app.applySelectedColormap();

        end

        % Value changed function: ClipSliderRange
        function ClipSliderRangeValueChanged(app, event)
            %CLIPSLIDERRANGEVALUECHANGING Preview CLim while dragging range slider.
            %CLIPSLIDERRANGEVALUECHANGED Commit CLim after range slider release.

            clim = double(app.ClipSliderRange.Value(:).');

            if numel(clim) ~= 2 || any(~isfinite(clim)) || clim(1) >= clim(2)
                app.ClipSliderRange.Value = app.ImageAxes.CLim;
                return
            end

            try
                app.setDisplayCLim(clim, false);
            catch ME
                app.ClipSliderRange.Value = app.ImageAxes.CLim;
                app.setStatusMessage(ME.message);
            end


        end

        % Value changing function: ClipSliderRange
        function ClipSliderRangeValueChanging(app, event)
            clim = double(event.Value(:).');

            if numel(clim) ~= 2 || any(~isfinite(clim)) || clim(1) >= clim(2)
                return
            end

            try
                app.setDisplayCLim(clim, false,true);
            catch ME
                app.setStatusMessage(ME.message);
            end


        end

        % Button pushed function: AutoButton
        function AutoButtonPushed(app, event)
            %AUTOBUTTONPUSHED Set CLim from center 75 percent of current frame.

            app.setAutoCLimFromCurrentFrame();
        end

        % Button pushed function: SetClipButton
        function SetClipButtonPushed(app, event)
            %SETCLIPBUTTONPUSHED Open modal editor for ClipSliderRange limits.
            %
            %   This callback edits the RangeSlider Limits property. Valid values are
            %   previewed live while the user types. OK confirms the preview. Cancel
            %   restores the limits/value present when the dialog opened. Restore uses
            %   the persistent data-loading limits stored in app.OriginalClipSliderLimits.

            previousMode = app.InteractionMode;
            app.setInteractionMode('settingClip');

            cleanupMode = onCleanup(@() app.setInteractionMode(previousMode));

            dialogStartLimits = double(app.ClipSliderRange.Limits);
            dialogStartValue  = double(app.ClipSliderRange.Value);

            restoreLimits = double(app.OriginalClipSliderLimits(:).');

            if numel(restoreLimits) ~= 2 || any(~isfinite(restoreLimits)) || restoreLimits(1) >= restoreLimits(2)
                restoreLimits = dialogStartLimits;
            end

            dlg = uifigure( ...
                'Name', 'Set clip slider limits', ...
                'WindowStyle', 'modal', ...
                'Position', [100 100 330 175], ...
                'Visible','off',...
                'CloseRequestFcn', @onCancel);
            placeAppInsideCaller(app,dlg,'center');
            dlg.Visible = 'on';
            grid = uigridlayout(dlg);
            grid.RowHeight = {30, 30, 25, 35};
            grid.ColumnWidth = {75, '1x'};
            grid.Padding = [12 12 12 12];

            minLabel = uilabel(grid);
            minLabel.Text = 'Min limit';
            minLabel.Layout.Row = 1;
            minLabel.Layout.Column = 1;

            minField = uieditfield(grid, 'numeric');
            minField.Value = dialogStartLimits(1);
            minField.Layout.Row = 1;
            minField.Layout.Column = 2;

            maxLabel = uilabel(grid);
            maxLabel.Text = 'Max limit';
            maxLabel.Layout.Row = 2;
            maxLabel.Layout.Column = 1;

            maxField = uieditfield(grid, 'numeric');
            maxField.Value = dialogStartLimits(2);
            maxField.Layout.Row = 2;
            maxField.Layout.Column = 2;

            statusLabel = uilabel(grid);
            statusLabel.Text = '';
            statusLabel.FontColor = [0.65 0 0];
            statusLabel.Layout.Row = 3;
            statusLabel.Layout.Column = [1 2];

            buttonGrid = uigridlayout(grid);
            buttonGrid.RowHeight = {'1x'};
            buttonGrid.ColumnWidth = {'1x', '1x', '1x'};
            buttonGrid.Padding = [0 0 0 0];
            buttonGrid.Layout.Row = 4;
            buttonGrid.Layout.Column = [1 2];

            okButton = uibutton(buttonGrid, 'push');
            okButton.Text = 'OK';
            okButton.Layout.Row = 1;
            okButton.Layout.Column = 1;
            okButton.ButtonPushedFcn = @onOK;

            restoreButton = uibutton(buttonGrid, 'push');
            restoreButton.Text = 'Restore';
            restoreButton.Layout.Row = 1;
            restoreButton.Layout.Column = 2;
            restoreButton.ButtonPushedFcn = @onRestore;

            cancelButton = uibutton(buttonGrid, 'push');
            cancelButton.Text = 'Cancel';
            cancelButton.Layout.Row = 1;
            cancelButton.Layout.Column = 3;
            cancelButton.ButtonPushedFcn = @onCancel;

            if isprop(minField, 'ValueChangingFcn')
                minField.ValueChangingFcn = @(src, evt) previewLimits(evt.Value, maxField.Value);
            end

            if isprop(maxField, 'ValueChangingFcn')
                maxField.ValueChangingFcn = @(src, evt) previewLimits(minField.Value, evt.Value);
            end

            minField.ValueChangedFcn = @(src, evt) previewLimits(minField.Value, maxField.Value);
            maxField.ValueChangedFcn = @(src, evt) previewLimits(minField.Value, maxField.Value);

            uiwait(dlg);

            if isvalid(dlg)
                delete(dlg);
            end

            function tf = validateLimits(minVal, maxVal)
                minVal = double(minVal);
                maxVal = double(maxVal);

                tf = isfinite(minVal) && isfinite(maxVal) && minVal < maxVal;

                if tf
                    statusLabel.Text = '';
                else
                    statusLabel.Text = 'Invalid limits: Min must be < Max.';
                end
            end

            function previewLimits(minVal, maxVal)
                if ~validateLimits(minVal, maxVal)
                    return
                end

                newLimits = double([minVal, maxVal]);
                applySliderLimits(newLimits, dialogStartValue);
            end

            function applySliderLimits(newLimits, baseValue)
                newValue = double(baseValue);
                newValue(1) = max(newValue(1), newLimits(1));
                newValue(2) = min(newValue(2), newLimits(2));

                if newValue(1) >= newValue(2)
                    newValue = newLimits;
                end

                currentLimits = double(app.ClipSliderRange.Limits);
                currentValue  = double(app.ClipSliderRange.Value);

                tempLimits = [ ...
                    min([currentLimits(1), newLimits(1), currentValue(1), newValue(1)]), ...
                    max([currentLimits(2), newLimits(2), currentValue(2), newValue(2)])];

                app.ClipSliderRange.Limits = tempLimits;
                app.ClipSliderRange.Value  = newValue;
                app.ClipSliderRange.Limits = newLimits;
                app.ClipSliderRange.Value  = newValue;

                app.ImageAxes.CLim = newValue;
                app.PlotAxes.YLimMode = 'manual';
                app.PlotAxes.YLim = newValue;

                drawnow limitrate
            end

            function onOK(~, ~)
                newLimits = double([minField.Value, maxField.Value]);

                if ~validateLimits(newLimits(1), newLimits(2))
                    return
                end

                applySliderLimits(newLimits, dialogStartValue);

                app.setStatusMessage(sprintf( ...
                    'Clip slider limits set to [%.4g, %.4g].', ...
                    newLimits(1), newLimits(2)));

                uiresume(dlg);
            end

            function onRestore(~, ~)
                minField.Value = restoreLimits(1);
                maxField.Value = restoreLimits(2);
                statusLabel.Text = '';

                applySliderLimits(restoreLimits, dialogStartValue);
            end

            function onCancel(~, ~)
                applySliderLimits(dialogStartLimits, dialogStartValue);
                uiresume(dlg);
            end
        end

        % Value changed function: HidecrosshairCheckBox
        function HidecrosshairCheckBoxValueChanged(app, event)
            %HIDECROSSHAIRCHECKBOXVALUECHANGED Toggle active crosshair visibility.
            %
            %   The checkbox text is "Hide crosshair", so checked means hidden.

            app.setCrosshairVisibility(~app.HidecrosshairCheckBox.Value);

        end

        % Value changed function: Switch
        function SwitchValueChanged(app, event)
            %SWITCHVALUECHANGED Switch between normal and event modes when allowed.

            if ~app.hasData()
                return
            end

            if strcmpi(app.getSourceType(), 'umt')
                % UMT shape fixes the mode. Revert the switch if user interaction occurs.
                app.syncSwitchToViewMode();
                return
            end

            if ~app.hasNormalizedEvents()
                app.ViewMode = 'normal';
                app.syncSwitchToViewMode();
                return
            end

            switch app.Switch.Value
                case 'Event mode'
                    app.ViewMode = 'event';
                otherwise
                    app.ViewMode = 'normal';
            end

            app.CurrentFrame = 1;
            app.rebuildEventFrameMatrix();
            app.refreshViewerForModeChange();



        end

        % Value changed function: ConditionDropDown
        function ConditionDropDownValueChanged(app, event)
            %CONDITIONDROPDOWNVALUECHANGED Update event selection condition.

            if ~app.hasNormalizedEvents()
                return
            end

            app.CurrentCondition = char(string(app.ConditionDropDown.Value));
            app.CurrentRepetition = 'AVERAGE';
            app.populateRepetitionDropDownForCurrentCondition();
            app.CurrentFrame = 1;
            app.rebuildEventFrameMatrix();
            app.refreshViewerForModeChange();



        end

        % Value changed function: RepetitionDropDown
        function RepetitionDropDownValueChanged(app, event)
            %REPETITIONDROPDOWNVALUECHANGED Update event selection repetition.

            if ~app.hasNormalizedEvents()
                return
            end

            app.CurrentRepetition = char(string(app.RepetitionDropDown.Value));
            app.CurrentFrame = 1;
            app.rebuildEventFrameMatrix();
            app.refreshViewerForModeChange();


        end

        % Button pushed function: DeleteConditionButton
        function DeleteConditionButtonPushed(app, event)

            %DELETECONDITIONBUTTONPUSHED Ignore the selected event condition for .dat data.

            if ~strcmpi(app.getSourceType(), 'dat') || isempty(app.EventSource)
                app.setStatusMessage('Condition deletion is only available for .dat event metadata.');
                return
            end

            if ~app.hasNormalizedEvents() || isempty(app.CurrentCondition)
                app.setStatusMessage('No event condition is available to delete.');
                return
            end

            conditionName = char(string(app.ConditionDropDown.Value));

            if isempty(conditionName) || strcmpi(conditionName, 'No events')
                app.setStatusMessage('No valid condition selected.');
                return
            end

            try
                app.EventSource.removeCondition(conditionName);
                app.setStatusMessage(sprintf( ...
                    'Condition "%s" ignored for display.', conditionName));

                app.refreshDatEventsAfterEdit();

            catch ME
                uialert(app.UIFigure, ME.message, 'Delete condition failed');
            end


        end

        % Button pushed function: DeleteRepetitionButton
        function DeleteRepetitionButtonPushed(app, event)

            %DELETEREPETITIONBUTTONPUSHED Ignore the selected repetition for .dat data.

            if ~strcmpi(app.getSourceType(), 'dat') || isempty(app.EventSource)
                app.setStatusMessage('Repetition deletion is only available for .dat event metadata.');
                return
            end

            if ~app.hasNormalizedEvents() || isempty(app.CurrentCondition)
                app.setStatusMessage('No event condition is available.');
                return
            end

            conditionName = char(string(app.ConditionDropDown.Value));
            repetitionValue = char(string(app.RepetitionDropDown.Value));

            if isempty(conditionName) || strcmpi(conditionName, 'No events')
                app.setStatusMessage('No valid condition selected.');
                return
            end

            if strcmpi(repetitionValue, 'AVERAGE')
                app.setStatusMessage('Select an individual repetition to delete.');
                return
            end

            repetitionIndex = str2double(repetitionValue);

            if ~isfinite(repetitionIndex) || repetitionIndex < 1 || repetitionIndex ~= round(repetitionIndex)
                app.setStatusMessage('Invalid repetition selected.');
                return
            end

            try
                app.EventSource.removeRepetition(conditionName, repetitionIndex);
                app.setStatusMessage(sprintf( ...
                    'Repetition %d from condition "%s" ignored for display.', ...
                    repetitionIndex, conditionName));

                app.refreshDatEventsAfterEdit();

            catch ME
                uialert(app.UIFigure, ME.message, 'Delete repetition failed');
            end


        end

        % Button pushed function: RestoreButton
        function RestoreButtonPushed(app, event)

            %RESTOREBUTTONPUSHED Restore all ignored .dat events.

            if ~strcmpi(app.getSourceType(), 'dat') || isempty(app.EventSource)
                app.setStatusMessage('Event restore is only available for .dat event metadata.');
                return
            end

            try
                app.EventSource.clearIgnoredEvents();

                if app.hasNormalizedEvents() && ~isempty(app.EventInfo.EventNames)
                    app.CurrentCondition = app.EventInfo.EventNames{1};
                end
                app.CurrentRepetition = 'AVERAGE';

                app.setStatusMessage('All ignored events restored.');
                app.refreshDatEventsAfterEdit();

            catch ME
                uialert(app.UIFigure, ME.message, 'Restore events failed');
            end


        end

        % Button pushed function: ViewCalibButton
        function ViewCalibButtonPushed(app, event)
            app.openViewCalibrationDialog();
        end

        % Button pushed function: DrawRectangleButton
        function DrawRectangleButtonPushed(app, event)
            app.createNewROI('rectangle');
        end

        % Button pushed function: DrawEllipseButton
        function DrawEllipseButtonPushed(app, event)
            app.createNewROI('ellipse');
        end

        % Button pushed function: DrawPolygonButton
        function DrawPolygonButtonPushed(app, event)
            app.createNewROI('polygon');
        end

        % Button pushed function: DeleteROIButton
        function DeleteROIButtonPushed(app, event)
            %DELETEROIBUTTONPUSHED Delete currently selected ROI(s).

            app.deleteSelectedROIs();
        end

        % Button pushed function: LogicalMaskButton
        function LogicalMaskButtonPushed(app, event)
            %CREATELOGICALMASKBUTTONPUSHED Start user-drawn logical-mask creation.

            app.startLogicalMaskDrawing();

        end

        % Button pushed function: SaveROIButton
        function SaveROIButtonPushed(app, event)
            %SAVEROIBUTTONPUSHED Save current ROIList to a .roi file.

            app.saveROIsToFile();

        end

        % Button pushed function: LoadROIButton
        function LoadROIButtonPushed(app, event)
            %LOADROIBUTTONPUSHED Load ROIList from a .roi file.

            app.loadROIsFromFile();

        end

        % Button pushed function: ImportROIButton
        function ImportROIButtonPushed(app, event)
            %IMPORTROIBUTTONPUSHED Import ROIs from a non-.roi source.

            app.importROIsFromExternalSource();

        end

        % Button pushed function: ExportROIDataButton
        function ExportROIDataButtonPushed(app, event)
            %EXPORTROIDATABUTTONPUSHED Export ROI measurements/traces to CSV files.

            app.exportROIDataFromDialog();


        end

        % Button pushed function: ROIbyTresholdButton
        function ROIbyTresholdButtonPushed(app, event)
            %ROIBYTHRESHOLDBUTTONPUSHED Open threshold-based ROI creation dialog.

            app.openThresholdROICreatorDialog();

        end

        % Button pushed function: ROIAllenButton
        function ROIAllenButtonPushed(app, event)
            %ALLENATLASROIBUTTONPUSHED Launch the Allen Atlas ROI selector child app.

            app.openAllenAtlasROICreator();

        end

        % Menu selected function: frombinMenu
        function frombinMenuSelected(app, event)
            %FROMBINMENUSELECTED Import raw .bin imaging data through PipelineManager.
            %
            %   This callback runs a one-time PipelineManager import workflow:
            %       1) Select raw and save folders.
            %       2) Run run_ImagesClassification.
            %       3) Optionally run getEvents for analog-input event import.
            %       4) Check PipelineManager.globalPipeLog for final status.
            %       5) Optionally open imported data from the selected save folder.
            %
            %   Pipeline success is determined only from globalPipeLog.Status.

            [rawFolder, saveFolder] = app.selectRawAndSaveFolderForDataImport('img_0000x.bin');

            if isempty(rawFolder)
                app.setStatusMessage('Data import cancelled.');
                return
            end

            previousMode = app.InteractionMode;
            app.setInteractionMode('loading');

            cleanupObj = onCleanup(@() app.setInteractionMode(previousMode)); %#ok<NASGU>

            try
                % -----------------------------------------------------------------
                % Build one-time import pipeline.
                % -----------------------------------------------------------------
                ppln = PipelineManager(saveFolder, rawFolder);

                importStep = ppln.addStep('run_ImagesClassification');
                if isempty(importStep)
                    app.setStatusMessage('Data import cancelled: run_ImagesClassification was not added.');
                    return
                end

                paramsOK = ppln.setParameters(importStep);
                if ~paramsOK
                    app.setStatusMessage('Data import cancelled: import parameters were not confirmed.');
                    return
                end

                % -----------------------------------------------------------------
                % Optional event import.
                % -----------------------------------------------------------------
                answer = uiconfirm(app.UIFigure, ...
                    'Import events from analog inputs (ai_0000x.bin)?', ...
                    'Import events', ...
                    'Options', {'Yes','No','Cancel'}, ...
                    'DefaultOption', 'No', ...
                    'CancelOption', 'Cancel');

                if strcmpi(answer, 'Cancel')
                    app.setStatusMessage('Data import cancelled.');
                    return
                end

                if strcmpi(answer, 'Yes')
                    eventStep = ppln.addStep('getEvents');
                    if isempty(eventStep)
                        app.setStatusMessage('Data import cancelled: getEvents was not added.');
                        return
                    end

                    paramsOK = ppln.setParameters(eventStep);
                    if ~paramsOK
                        app.setStatusMessage('Data import cancelled: event parameters were not confirmed.');
                        return
                    end
                end

                % -----------------------------------------------------------------
                % Execute pipeline.
                % -----------------------------------------------------------------
                app.setStatusMessage('Running data import pipeline...');
                drawnow limitrate

                ppln.executePipeline();

                % -----------------------------------------------------------------
                % Validate execution result from globalPipeLog only.
                % -----------------------------------------------------------------
                if isempty(ppln.globalPipeLog) || ~istable(ppln.globalPipeLog) || ...
                        ~ismember('Status', ppln.globalPipeLog.Properties.VariableNames)

                    errMsg = 'Pipeline execution did not return a valid globalPipeLog with a Status column.';

                    app.setStatusMessage(['Data import failed: ' errMsg]);

                    uialert(app.UIFigure, ...
                        errMsg, ...
                        'Data import failed', ...
                        'Icon', 'error');
                    return
                end

                statusVals = strip(string(ppln.globalPipeLog.Status));
                bImportSucceeded = all(strcmpi(statusVals, 'Completed'));

                if ~bImportSucceeded

                    if ismember('Messages_short', ppln.globalPipeLog.Properties.VariableNames)
                        msgVals = strip(string(ppln.globalPipeLog.Messages_short));
                        msgVals = msgVals(strlength(msgVals) > 0);
                        msgVals = msgVals(~strcmpi(msgVals, 'No Errors'));
                        msgVals = unique(msgVals, 'stable');
                    else
                        msgVals = strings(0, 1);
                    end

                    statusSummary = unique(statusVals, 'stable');
                    statusSummary = statusSummary(strlength(statusSummary) > 0);

                    if isempty(msgVals)
                        errMsg = sprintf('Pipeline status: %s', ...
                            char(strjoin(statusSummary, ', ')));
                    else
                        errMsg = char(strjoin(msgVals, newline));
                    end

                    app.setStatusMessage('Data import failed.');

                    uialert(app.UIFigure, ...
                        errMsg, ...
                        'Data import failed', ...
                        'Icon', 'error');
                    return
                end

                app.setStatusMessage('Data import completed.');
                app.openDataFromDialog(saveFolder);

            catch ME
                app.setStatusMessage(sprintf('Data import failed: %s', ME.message));

                try
                    uialert(app.UIFigure, ...
                        getReport(ME, 'basic', 'hyperlinks', 'off'), ...
                        'Data import failed', ...
                        'Icon', 'error');
                catch
                    warning('DataViewer:DataImportFailed', '%s', ...
                        getReport(ME, 'extended', 'hyperlinks', 'off'));
                end
            end

        end

        % Menu selected function: fromtifMenu
        function fromtifMenuSelected(app, event)
            %FROMTIFMENUSELECTED Import TIFF imaging data through PipelineManager.
            %
            %   This callback runs a one-time PipelineManager import workflow:
            %       1) Select raw and save folders.
            %       2) Run importFromTif.
            %       3) Optionally run getEvents for analog-input event import.
            %       4) Check PipelineManager.globalPipeLog for final status.
            %       5) Open imported data from the selected save folder.
            %
            %   Pipeline success is determined only from globalPipeLog.Status.

            [rawFolder, saveFolder] = app.selectRawAndSaveFolderForDataImport('info.json and .tif/.tiff');

            if isempty(rawFolder)
                app.setStatusMessage('TIFF import cancelled.');
                return
            end

            previousMode = app.InteractionMode;
            app.setInteractionMode('loading');

            cleanupObj = onCleanup(@() app.setInteractionMode(previousMode)); %#ok<NASGU>

            try
                % -----------------------------------------------------------------
                % Build one-time import pipeline.
                % -----------------------------------------------------------------
                ppln = PipelineManager(saveFolder, rawFolder);

                importStep = ppln.addStep('importFromTif');
                if isempty(importStep)
                    app.setStatusMessage('TIFF import cancelled: importFromTif was not added.');
                    return
                end

                paramsOK = ppln.setParameters(importStep);
                if ~paramsOK
                    app.setStatusMessage('TIFF import cancelled: import parameters were not confirmed.');
                    return
                end

                % -----------------------------------------------------------------
                % Execute pipeline.
                % -----------------------------------------------------------------
                app.setStatusMessage('Running TIFF import pipeline...');
                drawnow limitrate

                ppln.executePipeline();

                % -----------------------------------------------------------------
                % Validate execution result from globalPipeLog only.
                % -----------------------------------------------------------------
                if isempty(ppln.globalPipeLog) || ~istable(ppln.globalPipeLog) || ...
                        ~ismember('Status', ppln.globalPipeLog.Properties.VariableNames)

                    errMsg = 'Pipeline execution did not return a valid globalPipeLog with a Status column.';

                    app.setStatusMessage(['TIFF import failed: ' errMsg]);

                    uialert(app.UIFigure, ...
                        errMsg, ...
                        'TIFF import failed', ...
                        'Icon', 'error');
                    return
                end

                statusVals = strip(string(ppln.globalPipeLog.Status));
                bImportSucceeded = all(strcmpi(statusVals, 'Completed'));

                if ~bImportSucceeded

                    if ismember('Messages_short', ppln.globalPipeLog.Properties.VariableNames)
                        msgVals = strip(string(ppln.globalPipeLog.Messages_short));
                        msgVals = msgVals(strlength(msgVals) > 0);
                        msgVals = msgVals(~strcmpi(msgVals, 'No Errors'));
                        msgVals = unique(msgVals, 'stable');
                    else
                        msgVals = strings(0, 1);
                    end

                    statusSummary = unique(statusVals, 'stable');
                    statusSummary = statusSummary(strlength(statusSummary) > 0);

                    if isempty(msgVals)
                        errMsg = sprintf('Pipeline status: %s', ...
                            char(strjoin(statusSummary, ', ')));
                    else
                        errMsg = char(strjoin(msgVals, newline));
                    end

                    app.setStatusMessage('TIFF import failed.');

                    uialert(app.UIFigure, ...
                        errMsg, ...
                        'TIFF import failed', ...
                        'Icon', 'error');
                    return
                end

                app.setStatusMessage('TIFF import completed.');
                app.openDataFromDialog(saveFolder);

            catch ME
                app.setStatusMessage(sprintf('TIFF import failed: %s', ME.message));

                try
                    uialert(app.UIFigure, ...
                        getReport(ME, 'basic', 'hyperlinks', 'off'), ...
                        'TIFF import failed', ...
                        'Icon', 'error');
                catch
                    warning('DataViewer:TiffImportFailed', '%s', ...
                        getReport(ME, 'extended', 'hyperlinks', 'off'));
                end
            end
        end

        % Menu selected function: totifMenu
        function totifMenuSelected(app, event)

            %TOTIFMENUSELECTED Export the currently loaded image data to TIFF.
            %
            %   This callback exports the current .dat or .umt file using
            %   run_ConvertToTiff. Continuous YXT image data are exported as one
            %   TIFF stack. Event-split YXTE UMT image data are exported as one TIFF
            %   stack per event instance plus an event-info text manifest.
            %
            %   Export success is determined by successful completion of
            %   run_ConvertToTiff.

            if ~app.hasData() || isempty(app.CurrentFile)
                app.setStatusMessage('Load image data before exporting to TIFF.');

                try
                    uialert(app.UIFigure, ...
                        'Load image data before exporting to TIFF.', ...
                        'Export to TIFF', ...
                        'Icon', 'warning');
                catch
                    warning('DataViewer:NoDataForTiffExport', ...
                        'Load image data before exporting to TIFF.');
                end

                return
            end

            if ~isfile(app.CurrentFile)
                app.setStatusMessage('TIFF export failed: current data file was not found.');

                try
                    uialert(app.UIFigure, ...
                        sprintf('Current data file was not found:\n%s', app.CurrentFile), ...
                        'TIFF export failed', ...
                        'Icon', 'error');
                catch
                    warning('DataViewer:CurrentFileMissing', ...
                        'Current data file was not found: %s', app.CurrentFile);
                end

                return
            end

            [~, ~, ext] = fileparts(app.CurrentFile);
            ext = lower(ext);

            if ~ismember(ext, {'.dat', '.umt'})
                app.setStatusMessage('TIFF export failed: unsupported current file type.');

                try
                    uialert(app.UIFigure, ...
                        sprintf('Only .dat and .umt files can be exported to TIFF.\nCurrent file:\n%s', app.CurrentFile), ...
                        'TIFF export failed', ...
                        'Icon', 'error');
                catch
                    warning('DataViewer:UnsupportedTiffExportInput', ...
                        'Only .dat and .umt files can be exported to TIFF.');
                end

                return
            end

            if strcmpi(ext, '.umt') && ~isempty(app.CurrentEntry)
                % run_ConvertToTiff currently resolves .umt files by exporting the
                % first image entry. Warn because this may differ from the selected
                % displayed entry.
                answer = uiconfirm(app.UIFigure, ...
                    sprintf(['The loaded UMT entry is "%s".\n\n' ...
                    'run_ConvertToTiff currently exports the first image entry ' ...
                    'stored in the UMT file, which may not be the displayed entry.\n\n' ...
                    'Continue export?'], app.CurrentEntry), ...
                    'Confirm UMT TIFF export', ...
                    'Options', {'Continue','Cancel'}, ...
                    'DefaultOption', 'Continue', ...
                    'CancelOption', 'Cancel');

                if strcmpi(answer, 'Cancel')
                    app.setStatusMessage('TIFF export cancelled.');
                    return
                end
            end

            currentFolder = fileparts(app.CurrentFile);
            exportFolder = uigetdir(currentFolder, 'Select folder to save TIFF export');

            if isequal(exportFolder, 0)
                app.setStatusMessage('TIFF export cancelled.');
                return
            end

            previousMode = app.InteractionMode;
            app.setInteractionMode('loading');

            cleanupObj = onCleanup(@() app.setInteractionMode(previousMode)); %#ok<NASGU>

            try
                app.setStatusMessage('Exporting data to TIFF...');
                drawnow limitrate

                outFile = run_ConvertToTiff(app.CurrentFile, exportFolder);

                if isempty(outFile)
                    outText = 'No output files were reported.';
                else
                    outText = strjoin(string(outFile(:)), ', ');
                end

                app.setStatusMessage(sprintf('TIFF export completed: %s', char(outText)));

                try
                    uialert(app.UIFigure, ...
                        sprintf('TIFF export completed.\n\nOutput folder:\n%s\n\nFiles:\n%s', ...
                        exportFolder, char(outText)), ...
                        'TIFF export completed', ...
                        'Icon', 'success');
                catch
                    fprintf('[DataViewer] TIFF export completed.\n');
                    fprintf('Output folder: %s\n', exportFolder);
                    fprintf('Files: %s\n', char(outText));
                end

            catch ME
                app.setStatusMessage(sprintf('TIFF export failed: %s', ME.message));

                try
                    uialert(app.UIFigure, ...
                        getReport(ME, 'basic', 'hyperlinks', 'off'), ...
                        'TIFF export failed', ...
                        'Icon', 'error');
                catch
                    warning('DataViewer:TiffExportFailed', '%s', ...
                        getReport(ME, 'extended', 'hyperlinks', 'off'));
                end
            end


        end

        % Button pushed function: ImageRefButton
        function ImageRefButtonPushed(app, event)
            %IMAGEREFERENCEBUTTONPUSHED Open Image Reference Manager.

            try
                if app.hasData()
                    ctx = app.getImageReferenceManagerContext();

                    ImageReferenceManager(app, 'manage', ...
                        'currentFrameImage', ctx.currentFrameImage, ...
                        'sourceFile', ctx.sourceFile, ...
                        'sourceFolder', ctx.sourceFolder, ...
                        'sourceFrame', ctx.sourceFrame, ...
                        'sourceType', ctx.sourceType);
                else
                    % Manager can still preview/archive existing references. Creation
                    % from current frame will be disabled inside the manager.
                    ImageReferenceManager(app, 'manage');
                end

            catch ME
                app.setStatusMessage(sprintf( ...
                    'Failed to open Image Reference Manager: %s', ME.message));

                uialert(app.UIFigure, ...
                    getReport(ME, 'basic', 'hyperlinks', 'off'), ...
                    'Image Reference Manager failed', ...
                    'Icon', 'error');
            end
        end

        % Button pushed function: DataHistoryButton
        function DataHistoryButtonPushed(app, event)


            %VIEWDATAHISTORYBUTTONPUSHED Display dataHistory entry for current file.

            try
                app.openCurrentFileDataHistoryDialog();
            catch ME
                app.setStatusMessage(sprintf('View Data History failed: %s', ME.message));

                try
                    uialert(app.UIFigure, ...
                        getReport(ME, 'basic', 'hyperlinks', 'off'), ...
                        'View Data History failed', ...
                        'Icon', 'error');
                catch
                    warning('DataViewer:ViewDataHistoryFailed', '%s', ...
                        getReport(ME, 'extended', 'hyperlinks', 'off'));
                end
            end


        end

        % Button pushed function: EventsManagerButton
        function EventsManagerButtonPushed(app, event)
            %EVENTSMANAGERBUTTONPUSHED Open event manager for .dat data.
            %
            %   RawFolder is passed when available, but DataViewer-level launch is
            %   intentionally permissive because the Events Manager also supports manual
            %   event loading/editing workflows.

            if ~app.hasData()
                app.setStatusMessage('Load .dat data before opening Events Manager.');
                return
            end

            if ~strcmpi(app.getSourceType(), 'dat')
                uialert(app.UIFigure, ...
                    'Events Manager is only available for .dat data.', ...
                    'Unsupported data type', ...
                    'Icon', 'warning');
                return
            end

            saveFolder = fileparts(app.CurrentFile);
            rawFolder = app.resolveDataRawFolderForCurrentData(saveFolder);

            previousMode = app.ViewMode;

            try
                managerApp = DataViewer_EventsManager(app, saveFolder, rawFolder);

                % Wait until the modal Events Manager is actually closed before
                % reloading events.mat and refreshing the main DataViewer state.
                if ~isempty(managerApp) && isvalid(managerApp) && ...
                        isprop(managerApp, 'UIFigure') && isvalid(managerApp.UIFigure)
                    managerFig = managerApp.UIFigure;
                    waitfor(managerFig);
                end

                if ~app.hasData() || ~strcmpi(app.getSourceType(), 'dat')
                    return
                end

                % Reload events.mat and rebuild all event metadata/control state.
                app.initializeEventsForLoadedData('.dat');

                % Recover to normal source-frame mode after event metadata edits.
                app.ViewMode = 'normal';
                app.syncSwitchToViewMode();

                sourceSize = app.getSourceDataSize();
                if numel(sourceSize) >= 3
                    app.CurrentFrame = min(max(1, app.CurrentFrame), sourceSize(3));
                else
                    app.CurrentFrame = 1;
                end

                app.refreshFrameControls();
                app.refreshImageFrame();
                title(app.ImageAxes, app.getImageTitle());
                app.refreshTemporalProfile();
                app.refreshEventPatches();
                app.refreshROITraces();
                app.updateROIStatsForCurrentFrame();
                app.updateGUIEnabledState();

                if app.hasNormalizedEvents()
                    app.setStatusMessage('Events metadata reloaded. Event splitting is available from Normal Mode.');
                else
                    app.setStatusMessage('Events Manager closed. No valid events.mat metadata is available.');
                end

            catch ME
                app.setStatusMessage(sprintf('Events Manager failed: %s', ME.message));

                try
                    uialert(app.UIFigure, ...
                        getReport(ME, 'basic', 'hyperlinks', 'off'), ...
                        'Events Manager failed', ...
                        'Icon', 'error');
                catch
                    warning('DataViewer:ManageEventsFailed', '%s', ...
                        getReport(ME, 'extended', 'hyperlinks', 'off'));
                end

                if strcmpi(previousMode, 'event') && app.hasData()
                    try
                        app.ViewMode = 'normal';
                        app.syncSwitchToViewMode();
                        app.refreshFrameControls();
                        app.refreshImageFrame();
                        app.updateGUIEnabledState();
                    catch
                    end
                end
            end

        end

        % Button pushed function: OiSDUalCamCoregButton
        function OiSDUalCamCoregButtonPushed(app, event)
            %OISDUALCAMCOREGBUTTONPUSHED Open dual-camera calibration utility.

            [canOpen, ~, reasonText] = app.canOpenOiSDualCamCoreg();

            if ~canOpen
                app.setStatusMessage(reasonText);

                try
                    uialert(app.UIFigure, ...
                        reasonText, ...
                        'Dual-camera coregistration unavailable', ...
                        'Icon', 'warning');
                catch
                end
                return
            end

            dataFolder = fileparts(app.CurrentFile);

            try
                regApp = DataViewer_Coreg2Cams(app, dataFolder, 'CalibrationUtility');

                if ~isempty(regApp) && isvalid(regApp) && isvalid(regApp.UIFigure)
                    waitfor(regApp.UIFigure);
                end

                app.setStatusMessage('Dual-camera calibration utility closed.');

            catch ME
                app.setStatusMessage(sprintf('Dual-camera coregistration failed: %s', ME.message));

                try
                    uialert(app.UIFigure, ...
                        getReport(ME, 'basic', 'hyperlinks', 'off'), ...
                        'Dual-camera coregistration failed', ...
                        'Icon', 'error');
                catch
                    warning('DataViewer:DualCamCoregFailed', '%s', ...
                        getReport(ME, 'extended', 'hyperlinks', 'off'));
                end
            end

        end

        % Button pushed function: ImageAlignButton
        function ImageAlignButtonPushed(app, event)
            %IMAGEALIGNMENTTOOLBUTTONPUSHED Open ImageAlignmentTool from DataViewer.
            %
            %   The current displayed DataViewer frame is used as the moving source.
            %   The moving-source selector is disabled inside ImageAlignmentTool.

            if ~app.hasData()
                app.setStatusMessage('Load image data before opening Image Alignment Tool.');
                return
            end

            sourceType = lower(string(app.getSourceType()));

            if ~ismember(sourceType, ["dat", "umt"])
                uialert(app.UIFigure, ...
                    'Image Alignment Tool requires image data loaded from .dat or image-kind .umt.', ...
                    'Unsupported source type', ...
                    'Icon', 'warning');
                return
            end

            try
                ctx = app.getImageAlignmentToolContext();

                ImageAlignmentTool(app, ctx.targetFolder, 'dataviewer', ...
                    'movingImage', ctx.movingImage, ...
                    'movingSourceFile', ctx.sourceFile, ...
                    'movingSourceDescription', ctx.description, ...
                    'targetFolder', ctx.targetFolder, ...
                    'lockMovingSource', true);

                app.setStatusMessage('Image Alignment Tool opened.');

            catch ME
                app.printExceptionToConsole(ME, 'ImageAlignmentToolButtonPushed');

                app.setStatusMessage(sprintf( ...
                    'Failed to open Image Alignment Tool: %s', ME.message));

                uialert(app.UIFigure, ...
                    getReport(ME, 'basic', 'hyperlinks', 'off'), ...
                    'Image Alignment Tool failed', ...
                    'Icon', 'error');
            end

        end

        % Button pushed function: PipeLauncherButton
        function PipeLauncherButtonPushed(app, event)
            %PIPELINEMANAGERTOOLBUTTONPUSHED Launch PM Tool using persistent PM object.

            if ~app.hasData() || isempty(app.CurrentFile) || ~isfile(app.CurrentFile)
                uialert(app.UIFigure, ...
                    'Open a .dat or .umt file before launching Pipeline Manager.', ...
                    'No Data Loaded', ...
                    'Icon', 'warning');
                return
            end

            saveFolder = app.getCurrentSaveFolderForData();

            if isempty(saveFolder) || ~isfolder(saveFolder)
                uialert(app.UIFigure, ...
                    'The current data file does not have a valid SaveFolder.', ...
                    'Invalid Data Context', ...
                    'Icon', 'error');
                return
            end

            if ~app.ensurePipelineManagerForCurrentData()
                uialert(app.UIFigure, ...
                    'Could not create or reuse the PipelineManager object for the current data.', ...
                    'Pipeline Manager Error', ...
                    'Icon', 'error');
                return
            end

            rawFolder = app.resolveDataRawFolderForCurrentData(saveFolder);

            try
                if ~isempty(app.PipelineManagerToolApp) && isvalid(app.PipelineManagerToolApp)
                    try
                        figure(app.PipelineManagerToolApp.UIFigure);
                    catch
                    end
                    return
                end
            catch
                app.PipelineManagerToolApp = [];
            end

            try
                app.setInteractionMode('runningPipeline');

                try
                    app.UIFigure.Visible = 'off';
                catch
                end

                app.PipelineManagerToolApp = PipelineManagerTool( ...
                    'Mode', 'DataViewer', ...
                    'SaveFolder', saveFolder, ...
                    'RawFolder', rawFolder, ...
                    'CurrentFilePath', app.CurrentFile, ...
                    'PipelineManager', app.PipelineManagerObj, ...
                    'LastExecutionResult', app.LastPipelineResult, ...
                    'DataViewerFileSourceNodeID', app.PipelineDataViewerFileSourceNodeID, ...
                    'ExecutionFinishedFcn', @(result) app.onPipelineManagerFinished(result), ...
                    'ToolClosedFcn', @(toolState) app.onPipelineManagerClosed(toolState));

                app.updatePipelineTabState();
                app.updateDataFolderContextLabel();
                app.setStatusMessage('Pipeline Manager Tool opened.');

            catch ME
                app.PipelineManagerToolApp = [];

                try
                    app.UIFigure.Visible = 'on';
                    figure(app.UIFigure);
                catch
                end

                app.setInteractionMode('idle');
                app.updatePipelineTabState();
                app.setStatusMessage(sprintf('Failed to open Pipeline Manager Tool: %s', ME.message));
                uialert(app.UIFigure, ME.message, 'Pipeline Manager Tool Error', 'Icon', 'error');
            end


        end

        % Menu selected function: SetRawFolderMenu
        function SetRawFolderMenuSelected(app, event)
            %SETRAWFOLDERMENUSELECTED Set the RawFolder for the current data context.
            %
            % RawFolder is owned by DataViewer and persisted to DataParams.folders for
            % the current SaveFolder. PipelineManagerTool receives this value read-only
            % on launch.

            if ~app.hasData() || isempty(app.CurrentFile)
                uialert(app.UIFigure, ...
                    'Open a data file before setting the RawFolder.', ...
                    'No Data Loaded', ...
                    'Icon', 'warning');
                return
            end

            saveFolder = app.getCurrentSaveFolderForData();

            if isempty(saveFolder) || ~isfolder(saveFolder)
                uialert(app.UIFigure, ...
                    'The current data file does not have a valid SaveFolder.', ...
                    'Invalid Data Context', ...
                    'Icon', 'error');
                return
            end

            currentRawFolder = app.resolveDataRawFolderForCurrentData(saveFolder);

            msg = sprintf(['Set the RawFolder associated with the current SaveFolder.\n\n' ...
                'SaveFolder:\n%s\n\n' ...
                'Current RawFolder:\n%s\n\n' ...
                'Choose "Select Folder" to assign a RawFolder, "Use SaveFolder" when raw files are stored with the processed data, or "Set Missing" for datasets/pipelines that do not require raw data.'], ...
                saveFolder, currentRawFolder);

            choice = uiconfirm(app.UIFigure, ...
                msg, ...
                'Set RawFolder', ...
                'Options', {'Select Folder', 'Use SaveFolder', 'Set Missing', 'Cancel'}, ...
                'DefaultOption', 1, ...
                'CancelOption', 4);

            if strcmp(choice, 'Cancel')
                app.setStatusMessage('Set RawFolder cancelled.');
                return
            end

            if strcmp(choice, 'Set Missing')

                rawFolder = 'Missing';

            elseif strcmp(choice, 'Use SaveFolder')

                rawFolder = saveFolder;

            else

                startFolder = saveFolder;

                if ~isempty(currentRawFolder) && ...
                        ~strcmpi(char(string(currentRawFolder)), 'Missing') && ...
                        isfolder(currentRawFolder)
                    startFolder = currentRawFolder;
                end

                selectedFolder = uigetdir(startFolder, 'Select RawFolder');

                if isequal(selectedFolder, 0)
                    app.setStatusMessage('Set RawFolder cancelled.');
                    return
                end

                if isempty(selectedFolder) || ~isfolder(selectedFolder)
                    uialert(app.UIFigure, ...
                        'The selected RawFolder is not a valid folder.', ...
                        'Invalid RawFolder', ...
                        'Icon', 'error');
                    return
                end

                rawFolder = char(string(selectedFolder));
            end

            app.setRawFolderForCurrentSaveFolder(rawFolder, 'DataViewer', true);

        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [92 92 1400 900];
            app.UIFigure.Name = 'DataViewer';
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, @UIFigureCloseRequest, true);

            % Create FileMenu
            app.FileMenu = uimenu(app.UIFigure);
            app.FileMenu.Text = 'File';

            % Create OpenMenu
            app.OpenMenu = uimenu(app.FileMenu);
            app.OpenMenu.MenuSelectedFcn = createCallbackFcn(app, @OpenMenuSelected, true);
            app.OpenMenu.Accelerator = 'o';
            app.OpenMenu.Text = 'Open';

            % Create PreviewRawMenu
            app.PreviewRawMenu = uimenu(app.FileMenu);
            app.PreviewRawMenu.Enable = 'off';
            app.PreviewRawMenu.Text = 'Preview Raw';

            % Create SaveMenu
            app.SaveMenu = uimenu(app.FileMenu);
            app.SaveMenu.Accelerator = 's';
            app.SaveMenu.Text = 'Save';

            % Create ImportMenu
            app.ImportMenu = uimenu(app.FileMenu);
            app.ImportMenu.Text = 'Import';

            % Create frombinMenu
            app.frombinMenu = uimenu(app.ImportMenu);
            app.frombinMenu.MenuSelectedFcn = createCallbackFcn(app, @frombinMenuSelected, true);
            app.frombinMenu.Text = 'from .bin';

            % Create fromtifMenu
            app.fromtifMenu = uimenu(app.ImportMenu);
            app.fromtifMenu.MenuSelectedFcn = createCallbackFcn(app, @fromtifMenuSelected, true);
            app.fromtifMenu.Text = 'from .tif';

            % Create ExportMenu
            app.ExportMenu = uimenu(app.FileMenu);
            app.ExportMenu.Text = 'Export';

            % Create totifMenu
            app.totifMenu = uimenu(app.ExportMenu);
            app.totifMenu.MenuSelectedFcn = createCallbackFcn(app, @totifMenuSelected, true);
            app.totifMenu.Text = 'to tif';

            % Create EditmetadataMenu
            app.EditmetadataMenu = uimenu(app.FileMenu);
            app.EditmetadataMenu.Enable = 'off';
            app.EditmetadataMenu.Separator = 'on';
            app.EditmetadataMenu.Text = 'Edit meta data';

            % Create SetRawFolderMenu
            app.SetRawFolderMenu = uimenu(app.FileMenu);
            app.SetRawFolderMenu.MenuSelectedFcn = createCallbackFcn(app, @SetRawFolderMenuSelected, true);
            app.SetRawFolderMenu.Text = 'Set Raw Folder...';

            % Create HelpMenu_2
            app.HelpMenu_2 = uimenu(app.UIFigure);
            app.HelpMenu_2.Text = 'Help';

            % Create PreferencesMenu
            app.PreferencesMenu = uimenu(app.HelpMenu_2);
            app.PreferencesMenu.Enable = 'off';
            app.PreferencesMenu.Text = 'Preferences';

            % Create OnlinedocumentationMenu
            app.OnlinedocumentationMenu = uimenu(app.HelpMenu_2);
            app.OnlinedocumentationMenu.Enable = 'off';
            app.OnlinedocumentationMenu.Text = 'Online documentation';

            % Create GridLayoutMain
            app.GridLayoutMain = uigridlayout(app.UIFigure);
            app.GridLayoutMain.ColumnWidth = {'1x'};
            app.GridLayoutMain.RowHeight = {140, '1x'};

            % Create TabGroup
            app.TabGroup = uitabgroup(app.GridLayoutMain);
            app.TabGroup.Layout.Row = 1;
            app.TabGroup.Layout.Column = 1;

            % Create UtilitiesTab
            app.UtilitiesTab = uitab(app.TabGroup);
            app.UtilitiesTab.Title = 'Utilities';

            % Create GridUtilitiesTab
            app.GridUtilitiesTab = uigridlayout(app.UtilitiesTab);
            app.GridUtilitiesTab.ColumnWidth = {60, 60, 60, 60, 60, 60, 60, 60};
            app.GridUtilitiesTab.RowHeight = {60, 30};
            app.GridUtilitiesTab.ColumnSpacing = 30;
            app.GridUtilitiesTab.RowSpacing = 5;

            % Create ViewCalibButton
            app.ViewCalibButton = uibutton(app.GridUtilitiesTab, 'push');
            app.ViewCalibButton.ButtonPushedFcn = createCallbackFcn(app, @ViewCalibButtonPushed, true);
            app.ViewCalibButton.Icon = 'ViewCalibration.png';
            app.ViewCalibButton.Tooltip = {'View or edit the spatial calibration used to convert pixels to physical units'};
            app.ViewCalibButton.Layout.Row = 1;
            app.ViewCalibButton.Layout.Column = 1;
            app.ViewCalibButton.Text = '';

            % Create LogicalMaskButton
            app.LogicalMaskButton = uibutton(app.GridUtilitiesTab, 'push');
            app.LogicalMaskButton.ButtonPushedFcn = createCallbackFcn(app, @LogicalMaskButtonPushed, true);
            app.LogicalMaskButton.Icon = 'CreateLogicalMask.png';
            app.LogicalMaskButton.IconAlignment = 'center';
            app.LogicalMaskButton.Tooltip = {'Create Logical mask'};
            app.LogicalMaskButton.Layout.Row = 1;
            app.LogicalMaskButton.Layout.Column = 2;
            app.LogicalMaskButton.Text = '';

            % Create ImageRefButton
            app.ImageRefButton = uibutton(app.GridUtilitiesTab, 'push');
            app.ImageRefButton.ButtonPushedFcn = createCallbackFcn(app, @ImageRefButtonPushed, true);
            app.ImageRefButton.Icon = 'ImageReferenceTool.png';
            app.ImageRefButton.IconAlignment = 'center';
            app.ImageRefButton.Tooltip = {'Create, load, or update the reference image used to align imaging data.'};
            app.ImageRefButton.Layout.Row = 1;
            app.ImageRefButton.Layout.Column = 3;
            app.ImageRefButton.Text = '';

            % Create DataHistoryButton
            app.DataHistoryButton = uibutton(app.GridUtilitiesTab, 'push');
            app.DataHistoryButton.ButtonPushedFcn = createCallbackFcn(app, @DataHistoryButtonPushed, true);
            app.DataHistoryButton.Icon = 'ViewDataHistory.png';
            app.DataHistoryButton.IconAlignment = 'center';
            app.DataHistoryButton.Tooltip = {'View the processing history and file lineage for the current dataset.'};
            app.DataHistoryButton.Layout.Row = 1;
            app.DataHistoryButton.Layout.Column = 5;
            app.DataHistoryButton.Text = '';

            % Create EventsManagerButton
            app.EventsManagerButton = uibutton(app.GridUtilitiesTab, 'push');
            app.EventsManagerButton.ButtonPushedFcn = createCallbackFcn(app, @EventsManagerButtonPushed, true);
            app.EventsManagerButton.Icon = 'EventsManagerTool.png';
            app.EventsManagerButton.IconAlignment = 'center';
            app.EventsManagerButton.Tooltip = {'View / edit experimental event timestamps for the current dataset.'};
            app.EventsManagerButton.Layout.Row = 1;
            app.EventsManagerButton.Layout.Column = 6;
            app.EventsManagerButton.Text = '';

            % Create ViewCalibrationLabel
            app.ViewCalibrationLabel = uilabel(app.GridUtilitiesTab);
            app.ViewCalibrationLabel.HorizontalAlignment = 'center';
            app.ViewCalibrationLabel.FontSize = 10;
            app.ViewCalibrationLabel.Layout.Row = 2;
            app.ViewCalibrationLabel.Layout.Column = 1;
            app.ViewCalibrationLabel.Text = {'View '; 'Calibration'};

            % Create LogicalMaskLabel
            app.LogicalMaskLabel = uilabel(app.GridUtilitiesTab);
            app.LogicalMaskLabel.HorizontalAlignment = 'center';
            app.LogicalMaskLabel.Layout.Row = 2;
            app.LogicalMaskLabel.Layout.Column = 2;
            app.LogicalMaskLabel.Text = {'Logical'; 'Mask'};

            % Create ImageReferenceLabel
            app.ImageReferenceLabel = uilabel(app.GridUtilitiesTab);
            app.ImageReferenceLabel.HorizontalAlignment = 'center';
            app.ImageReferenceLabel.Layout.Row = 2;
            app.ImageReferenceLabel.Layout.Column = 3;
            app.ImageReferenceLabel.Text = {'Image '; 'Reference'};

            % Create DataHistoryLabel
            app.DataHistoryLabel = uilabel(app.GridUtilitiesTab);
            app.DataHistoryLabel.HorizontalAlignment = 'center';
            app.DataHistoryLabel.Layout.Row = 2;
            app.DataHistoryLabel.Layout.Column = 5;
            app.DataHistoryLabel.Text = {'Data'; ' History'};

            % Create EventsManagerLabel
            app.EventsManagerLabel = uilabel(app.GridUtilitiesTab);
            app.EventsManagerLabel.HorizontalAlignment = 'center';
            app.EventsManagerLabel.Layout.Row = 2;
            app.EventsManagerLabel.Layout.Column = 6;
            app.EventsManagerLabel.Text = {'Events '; 'Manager'};

            % Create OiSDUalCamCoregButton
            app.OiSDUalCamCoregButton = uibutton(app.GridUtilitiesTab, 'push');
            app.OiSDUalCamCoregButton.ButtonPushedFcn = createCallbackFcn(app, @OiSDUalCamCoregButtonPushed, true);
            app.OiSDUalCamCoregButton.Icon = 'OiSDualCamCoregistration.png';
            app.OiSDUalCamCoregButton.IconAlignment = 'center';
            app.OiSDUalCamCoregButton.Layout.Row = 1;
            app.OiSDUalCamCoregButton.Layout.Column = 7;
            app.OiSDUalCamCoregButton.Text = '';

            % Create OiS2CamAlignmentLabel
            app.OiS2CamAlignmentLabel = uilabel(app.GridUtilitiesTab);
            app.OiS2CamAlignmentLabel.HorizontalAlignment = 'center';
            app.OiS2CamAlignmentLabel.Layout.Row = 2;
            app.OiS2CamAlignmentLabel.Layout.Column = 7;
            app.OiS2CamAlignmentLabel.Text = {'OiS 2-Cam'; 'Alignment'};

            % Create ImageAlignButton
            app.ImageAlignButton = uibutton(app.GridUtilitiesTab, 'push');
            app.ImageAlignButton.ButtonPushedFcn = createCallbackFcn(app, @ImageAlignButtonPushed, true);
            app.ImageAlignButton.Icon = 'ImageAlignmentTool.png';
            app.ImageAlignButton.IconAlignment = 'center';
            app.ImageAlignButton.Tooltip = {'''Open the image alignment tool to register the current widefield image to a reference image.'};
            app.ImageAlignButton.Layout.Row = 1;
            app.ImageAlignButton.Layout.Column = 4;
            app.ImageAlignButton.Text = '';

            % Create AlignmentToolLabel
            app.AlignmentToolLabel = uilabel(app.GridUtilitiesTab);
            app.AlignmentToolLabel.HorizontalAlignment = 'center';
            app.AlignmentToolLabel.Layout.Row = 2;
            app.AlignmentToolLabel.Layout.Column = 4;
            app.AlignmentToolLabel.Text = {'Alignment'; 'Tool'};

            % Create PipelineTab
            app.PipelineTab = uitab(app.TabGroup);
            app.PipelineTab.Title = 'Pipeline';

            % Create GridPipelineTab
            app.GridPipelineTab = uigridlayout(app.PipelineTab);
            app.GridPipelineTab.ColumnWidth = {60, '1x'};
            app.GridPipelineTab.RowHeight = {20, 20, 20, 30};
            app.GridPipelineTab.ColumnSpacing = 30;
            app.GridPipelineTab.RowSpacing = 5;

            % Create PipeLauncherButton
            app.PipeLauncherButton = uibutton(app.GridPipelineTab, 'push');
            app.PipeLauncherButton.ButtonPushedFcn = createCallbackFcn(app, @PipeLauncherButtonPushed, true);
            app.PipeLauncherButton.Icon = 'PipelineManagerTool.png';
            app.PipeLauncherButton.IconAlignment = 'center';
            app.PipeLauncherButton.Layout.Row = [1 3];
            app.PipeLauncherButton.Layout.Column = 1;
            app.PipeLauncherButton.Text = '';

            % Create PipelineManagerLabel
            app.PipelineManagerLabel = uilabel(app.GridPipelineTab);
            app.PipelineManagerLabel.HorizontalAlignment = 'center';
            app.PipelineManagerLabel.Layout.Row = 4;
            app.PipelineManagerLabel.Layout.Column = 1;
            app.PipelineManagerLabel.Text = {'Pipeline '; 'Manager'};

            % Create PipelineStatusLabel
            app.PipelineStatusLabel = uilabel(app.GridPipelineTab);
            app.PipelineStatusLabel.Layout.Row = 1;
            app.PipelineStatusLabel.Layout.Column = 2;
            app.PipelineStatusLabel.Text = 'Pipeline Status';

            % Create LastrunoutputsLabel
            app.LastrunoutputsLabel = uilabel(app.GridPipelineTab);
            app.LastrunoutputsLabel.Layout.Row = 2;
            app.LastrunoutputsLabel.Layout.Column = 2;
            app.LastrunoutputsLabel.Text = 'Last run outputs';

            % Create CurrentpipelinefileLabel
            app.CurrentpipelinefileLabel = uilabel(app.GridPipelineTab);
            app.CurrentpipelinefileLabel.Layout.Row = 3;
            app.CurrentpipelinefileLabel.Layout.Column = 2;
            app.CurrentpipelinefileLabel.Text = 'Current pipeline file';

            % Create LatestsummaryLabel
            app.LatestsummaryLabel = uilabel(app.GridPipelineTab);
            app.LatestsummaryLabel.Layout.Row = 4;
            app.LatestsummaryLabel.Layout.Column = 2;
            app.LatestsummaryLabel.Text = 'Latest summary';

            % Create ROIsTab
            app.ROIsTab = uitab(app.TabGroup);
            app.ROIsTab.Title = 'ROIs';

            % Create GridROIsTab
            app.GridROIsTab = uigridlayout(app.ROIsTab);
            app.GridROIsTab.ColumnWidth = {60, 60, 60, 60, 60, 60, 5, 60, 60, 60, 60};
            app.GridROIsTab.RowHeight = {60, 30};
            app.GridROIsTab.ColumnSpacing = 30;
            app.GridROIsTab.RowSpacing = 5;

            % Create DrawRectangleButton
            app.DrawRectangleButton = uibutton(app.GridROIsTab, 'push');
            app.DrawRectangleButton.ButtonPushedFcn = createCallbackFcn(app, @DrawRectangleButtonPushed, true);
            app.DrawRectangleButton.Icon = 'ROI_drawRectangle.png';
            app.DrawRectangleButton.IconAlignment = 'center';
            app.DrawRectangleButton.Tooltip = {'Draw Rectangle'};
            app.DrawRectangleButton.Layout.Row = 1;
            app.DrawRectangleButton.Layout.Column = 1;
            app.DrawRectangleButton.Text = '';

            % Create DrawEllipseButton
            app.DrawEllipseButton = uibutton(app.GridROIsTab, 'push');
            app.DrawEllipseButton.ButtonPushedFcn = createCallbackFcn(app, @DrawEllipseButtonPushed, true);
            app.DrawEllipseButton.Icon = 'ROI_drawEllipse.png';
            app.DrawEllipseButton.Tooltip = {'Draw Ellipse or Circle'};
            app.DrawEllipseButton.Layout.Row = 1;
            app.DrawEllipseButton.Layout.Column = 2;
            app.DrawEllipseButton.Text = '';

            % Create DrawPolygonButton
            app.DrawPolygonButton = uibutton(app.GridROIsTab, 'push');
            app.DrawPolygonButton.ButtonPushedFcn = createCallbackFcn(app, @DrawPolygonButtonPushed, true);
            app.DrawPolygonButton.Icon = 'ROI_drawPolygon.png';
            app.DrawPolygonButton.Tooltip = {'Draw Polygon'};
            app.DrawPolygonButton.Layout.Row = 1;
            app.DrawPolygonButton.Layout.Column = 3;
            app.DrawPolygonButton.Text = '';

            % Create DeleteROIButton
            app.DeleteROIButton = uibutton(app.GridROIsTab, 'push');
            app.DeleteROIButton.ButtonPushedFcn = createCallbackFcn(app, @DeleteROIButtonPushed, true);
            app.DeleteROIButton.Icon = 'ROI_delete.png';
            app.DeleteROIButton.IconAlignment = 'center';
            app.DeleteROIButton.Tooltip = {'Delete selected ROI(s)'};
            app.DeleteROIButton.Layout.Row = 1;
            app.DeleteROIButton.Layout.Column = 6;
            app.DeleteROIButton.Text = '';

            % Create SaveROIButton
            app.SaveROIButton = uibutton(app.GridROIsTab, 'push');
            app.SaveROIButton.ButtonPushedFcn = createCallbackFcn(app, @SaveROIButtonPushed, true);
            app.SaveROIButton.Icon = 'ROI_save.png';
            app.SaveROIButton.IconAlignment = 'center';
            app.SaveROIButton.Tooltip = {'Save the current ROI set to file.'};
            app.SaveROIButton.Layout.Row = 1;
            app.SaveROIButton.Layout.Column = 8;
            app.SaveROIButton.Text = '';

            % Create RectangleLabel
            app.RectangleLabel = uilabel(app.GridROIsTab);
            app.RectangleLabel.HorizontalAlignment = 'center';
            app.RectangleLabel.Layout.Row = 2;
            app.RectangleLabel.Layout.Column = 1;
            app.RectangleLabel.Text = 'Rectangle';

            % Create EllipeLabel
            app.EllipeLabel = uilabel(app.GridROIsTab);
            app.EllipeLabel.HorizontalAlignment = 'center';
            app.EllipeLabel.Layout.Row = 2;
            app.EllipeLabel.Layout.Column = 2;
            app.EllipeLabel.Text = 'Ellipe';

            % Create PolygonLabel
            app.PolygonLabel = uilabel(app.GridROIsTab);
            app.PolygonLabel.HorizontalAlignment = 'center';
            app.PolygonLabel.Tooltip = {'Draw a polygon'};
            app.PolygonLabel.Layout.Row = 2;
            app.PolygonLabel.Layout.Column = 3;
            app.PolygonLabel.Text = 'Polygon';

            % Create DeleteLabel
            app.DeleteLabel = uilabel(app.GridROIsTab);
            app.DeleteLabel.HorizontalAlignment = 'center';
            app.DeleteLabel.Layout.Row = 2;
            app.DeleteLabel.Layout.Column = 6;
            app.DeleteLabel.Text = 'Delete ';

            % Create SaveLabel
            app.SaveLabel = uilabel(app.GridROIsTab);
            app.SaveLabel.HorizontalAlignment = 'center';
            app.SaveLabel.Layout.Row = 2;
            app.SaveLabel.Layout.Column = 8;
            app.SaveLabel.Text = 'Save ';

            % Create LoadROIButton
            app.LoadROIButton = uibutton(app.GridROIsTab, 'push');
            app.LoadROIButton.ButtonPushedFcn = createCallbackFcn(app, @LoadROIButtonPushed, true);
            app.LoadROIButton.Icon = 'ROI_load.png';
            app.LoadROIButton.IconAlignment = 'center';
            app.LoadROIButton.Tooltip = {'Load ROI definitions from file.'};
            app.LoadROIButton.Layout.Row = 1;
            app.LoadROIButton.Layout.Column = 9;
            app.LoadROIButton.Text = '';

            % Create LoadLabel
            app.LoadLabel = uilabel(app.GridROIsTab);
            app.LoadLabel.HorizontalAlignment = 'center';
            app.LoadLabel.Layout.Row = 2;
            app.LoadLabel.Layout.Column = 9;
            app.LoadLabel.Text = 'Load ';

            % Create ImportROIButton
            app.ImportROIButton = uibutton(app.GridROIsTab, 'push');
            app.ImportROIButton.ButtonPushedFcn = createCallbackFcn(app, @ImportROIButtonPushed, true);
            app.ImportROIButton.Icon = 'ROI_import.png';
            app.ImportROIButton.Tooltip = {'Import ROs from an external source.'};
            app.ImportROIButton.Layout.Row = 1;
            app.ImportROIButton.Layout.Column = 10;
            app.ImportROIButton.Text = '';

            % Create ExportROIDataButton
            app.ExportROIDataButton = uibutton(app.GridROIsTab, 'push');
            app.ExportROIDataButton.ButtonPushedFcn = createCallbackFcn(app, @ExportROIDataButtonPushed, true);
            app.ExportROIDataButton.Icon = 'ROI_export.png';
            app.ExportROIDataButton.Tooltip = {'Export extracted ROI time series or summary measurements to CSV file(s)'};
            app.ExportROIDataButton.Layout.Row = 1;
            app.ExportROIDataButton.Layout.Column = 11;
            app.ExportROIDataButton.Text = '';

            % Create ROIbyTresholdButton
            app.ROIbyTresholdButton = uibutton(app.GridROIsTab, 'push');
            app.ROIbyTresholdButton.ButtonPushedFcn = createCallbackFcn(app, @ROIbyTresholdButtonPushed, true);
            app.ROIbyTresholdButton.Icon = 'ROI_byThreshold.png';
            app.ROIbyTresholdButton.IconAlignment = 'center';
            app.ROIbyTresholdButton.Tooltip = {'Create ROI by Threshold'};
            app.ROIbyTresholdButton.Layout.Row = 1;
            app.ROIbyTresholdButton.Layout.Column = 4;
            app.ROIbyTresholdButton.Text = '';

            % Create ROIAllenButton
            app.ROIAllenButton = uibutton(app.GridROIsTab, 'push');
            app.ROIAllenButton.ButtonPushedFcn = createCallbackFcn(app, @ROIAllenButtonPushed, true);
            app.ROIAllenButton.Icon = 'ROI_AllenAtlas.png';
            app.ROIAllenButton.Tooltip = {'Create ROIs using the Mouse Allen Brain Atlas.'};
            app.ROIAllenButton.Layout.Row = 1;
            app.ROIAllenButton.Layout.Column = 5;
            app.ROIAllenButton.Text = '';

            % Create ImportLabel
            app.ImportLabel = uilabel(app.GridROIsTab);
            app.ImportLabel.HorizontalAlignment = 'center';
            app.ImportLabel.Layout.Row = 2;
            app.ImportLabel.Layout.Column = 10;
            app.ImportLabel.Text = 'Import ';

            % Create ExportROIdataLabel
            app.ExportROIdataLabel = uilabel(app.GridROIsTab);
            app.ExportROIdataLabel.HorizontalAlignment = 'center';
            app.ExportROIdataLabel.Layout.Row = 2;
            app.ExportROIdataLabel.Layout.Column = 11;
            app.ExportROIdataLabel.Text = {'Export ROI'; 'data'};

            % Create ThresholdLabel
            app.ThresholdLabel = uilabel(app.GridROIsTab);
            app.ThresholdLabel.HorizontalAlignment = 'center';
            app.ThresholdLabel.Layout.Row = 2;
            app.ThresholdLabel.Layout.Column = 4;
            app.ThresholdLabel.Text = 'Threshold';

            % Create AllenAtlasLabel
            app.AllenAtlasLabel = uilabel(app.GridROIsTab);
            app.AllenAtlasLabel.HorizontalAlignment = 'center';
            app.AllenAtlasLabel.Layout.Row = 2;
            app.AllenAtlasLabel.Layout.Column = 5;
            app.AllenAtlasLabel.Text = 'Allen Atlas';

            % Create Sep
            app.Sep = uipanel(app.GridROIsTab);
            app.Sep.BorderType = 'none';
            app.Sep.BackgroundColor = [0.8 0.8 0.8];
            app.Sep.Layout.Row = [1 2];
            app.Sep.Layout.Column = 7;

            % Create GridLayoutSub1
            app.GridLayoutSub1 = uigridlayout(app.GridLayoutMain);
            app.GridLayoutSub1.ColumnWidth = {'1.5x', '1x'};
            app.GridLayoutSub1.RowHeight = {'1x'};
            app.GridLayoutSub1.ColumnSpacing = 0;
            app.GridLayoutSub1.RowSpacing = 0;
            app.GridLayoutSub1.Padding = [0 0 0 0];
            app.GridLayoutSub1.Layout.Row = 2;
            app.GridLayoutSub1.Layout.Column = 1;

            % Create GridLayoutLeft
            app.GridLayoutLeft = uigridlayout(app.GridLayoutSub1);
            app.GridLayoutLeft.ColumnWidth = {'1x'};
            app.GridLayoutLeft.RowHeight = {20, 20, 20, 80, '1x', 50};
            app.GridLayoutLeft.RowSpacing = 2;
            app.GridLayoutLeft.Layout.Row = 1;
            app.GridLayoutLeft.Layout.Column = 1;

            % Create ImageAxes
            app.ImageAxes = uiaxes(app.GridLayoutLeft);
            title(app.ImageAxes, 'Title')
            xlabel(app.ImageAxes, 'X')
            ylabel(app.ImageAxes, 'Y')
            zlabel(app.ImageAxes, 'Z')
            app.ImageAxes.Layout.Row = 5;
            app.ImageAxes.Layout.Column = 1;

            % Create BottomGrid
            app.BottomGrid = uigridlayout(app.GridLayoutLeft);
            app.BottomGrid.ColumnWidth = {60, '1x', 60, 50, '1x', 50, 60, 100};
            app.BottomGrid.RowHeight = {30, 20};
            app.BottomGrid.RowSpacing = 3;
            app.BottomGrid.Padding = [10 0 10 0];
            app.BottomGrid.Layout.Row = 6;
            app.BottomGrid.Layout.Column = 1;

            % Create ColormapDropDownLabel
            app.ColormapDropDownLabel = uilabel(app.BottomGrid);
            app.ColormapDropDownLabel.Layout.Row = 1;
            app.ColormapDropDownLabel.Layout.Column = 1;
            app.ColormapDropDownLabel.Text = 'Colormap';

            % Create ColormapDropDown
            app.ColormapDropDown = uidropdown(app.BottomGrid);
            app.ColormapDropDown.Items = {'parula', 'turbo', 'jet', 'hsv', 'hot', 'gray'};
            app.ColormapDropDown.ValueChangedFcn = createCallbackFcn(app, @ColormapDropDownValueChanged, true);
            app.ColormapDropDown.Layout.Row = 1;
            app.ColormapDropDown.Layout.Column = 2;
            app.ColormapDropDown.Value = 'parula';

            % Create InvertCheckBox
            app.InvertCheckBox = uicheckbox(app.BottomGrid);
            app.InvertCheckBox.ValueChangedFcn = createCallbackFcn(app, @InvertCheckBoxValueChanged, true);
            app.InvertCheckBox.Tooltip = {'Inverts the colormap'};
            app.InvertCheckBox.Text = 'Invert';
            app.InvertCheckBox.Layout.Row = 1;
            app.InvertCheckBox.Layout.Column = 3;

            % Create AutoButton
            app.AutoButton = uibutton(app.BottomGrid, 'push');
            app.AutoButton.ButtonPushedFcn = createCallbackFcn(app, @AutoButtonPushed, true);
            app.AutoButton.Layout.Row = 1;
            app.AutoButton.Layout.Column = 6;
            app.AutoButton.Text = 'Auto';

            % Create Panel
            app.Panel = uipanel(app.BottomGrid);
            app.Panel.BorderType = 'none';
            app.Panel.Layout.Row = 1;
            app.Panel.Layout.Column = 5;

            % Create ClipSliderRange
            app.ClipSliderRange = uislider(app.Panel, 'range');
            app.ClipSliderRange.MajorTicks = [];
            app.ClipSliderRange.ValueChangedFcn = createCallbackFcn(app, @ClipSliderRangeValueChanged, true);
            app.ClipSliderRange.ValueChangingFcn = createCallbackFcn(app, @ClipSliderRangeValueChanging, true);
            app.ClipSliderRange.MinorTicks = [0 4 8 12 16 20 24 28 32 36 40 44 48 52 56 60 64 68 72 76 80 84 88 92 96 100];
            app.ClipSliderRange.Position = [16 13 145 3];

            % Create HidecrosshairCheckBox
            app.HidecrosshairCheckBox = uicheckbox(app.BottomGrid);
            app.HidecrosshairCheckBox.ValueChangedFcn = createCallbackFcn(app, @HidecrosshairCheckBoxValueChanged, true);
            app.HidecrosshairCheckBox.Text = 'Hide crosshair';
            app.HidecrosshairCheckBox.Layout.Row = 1;
            app.HidecrosshairCheckBox.Layout.Column = 8;

            % Create ClipSliderLabel
            app.ClipSliderLabel = uilabel(app.BottomGrid);
            app.ClipSliderLabel.HorizontalAlignment = 'center';
            app.ClipSliderLabel.Layout.Row = 1;
            app.ClipSliderLabel.Layout.Column = 4;
            app.ClipSliderLabel.Text = 'Clip';

            % Create SetClipButton
            app.SetClipButton = uibutton(app.BottomGrid, 'push');
            app.SetClipButton.ButtonPushedFcn = createCallbackFcn(app, @SetClipButtonPushed, true);
            app.SetClipButton.Tooltip = {'Set new range for Clip slider'};
            app.SetClipButton.Layout.Row = 1;
            app.SetClipButton.Layout.Column = 7;
            app.SetClipButton.Text = 'Set Limits';

            % Create ImageStatusLabel
            app.ImageStatusLabel = uilabel(app.BottomGrid);
            app.ImageStatusLabel.FontWeight = 'bold';
            app.ImageStatusLabel.Layout.Row = 2;
            app.ImageStatusLabel.Layout.Column = [1 4];
            app.ImageStatusLabel.Text = 'Image status';

            % Create TopGrid
            app.TopGrid = uigridlayout(app.GridLayoutLeft);
            app.TopGrid.ColumnWidth = {50, 35, 50, '1x', '1x', 35};
            app.TopGrid.RowHeight = {35, 30};
            app.TopGrid.Padding = [0 0 0 0];
            app.TopGrid.Layout.Row = 4;
            app.TopGrid.Layout.Column = 1;

            % Create PlayMovieButton
            app.PlayMovieButton = uibutton(app.TopGrid, 'push');
            app.PlayMovieButton.ButtonPushedFcn = createCallbackFcn(app, @PlayMovieButtonPushed, true);
            app.PlayMovieButton.Icon = 'icon_rightTriangle_Black.png';
            app.PlayMovieButton.Layout.Row = 1;
            app.PlayMovieButton.Layout.Column = 1;
            app.PlayMovieButton.Text = '';

            % Create PreviousFrameButton
            app.PreviousFrameButton = uibutton(app.TopGrid, 'push');
            app.PreviousFrameButton.ButtonPushedFcn = createCallbackFcn(app, @PreviousFrameButtonPushed, true);
            app.PreviousFrameButton.Layout.Row = 1;
            app.PreviousFrameButton.Layout.Column = 2;
            app.PreviousFrameButton.Text = '<';

            % Create MovieSpeedDropDown
            app.MovieSpeedDropDown = uidropdown(app.TopGrid);
            app.MovieSpeedDropDown.Items = {'0.25x', '0.5x', '1x', '2x', '4x'};
            app.MovieSpeedDropDown.ValueChangedFcn = createCallbackFcn(app, @MovieSpeedDropDownValueChanged, true);
            app.MovieSpeedDropDown.Layout.Row = 2;
            app.MovieSpeedDropDown.Layout.Column = 1;
            app.MovieSpeedDropDown.Value = '1x';

            % Create MovieSpeedLabel
            app.MovieSpeedLabel = uilabel(app.TopGrid);
            app.MovieSpeedLabel.Layout.Row = 2;
            app.MovieSpeedLabel.Layout.Column = [2 3];
            app.MovieSpeedLabel.Text = 'FPS';

            % Create Slider
            app.Slider = uislider(app.TopGrid);
            app.Slider.ValueChangedFcn = createCallbackFcn(app, @SliderValueChanged, true);
            app.Slider.ValueChangingFcn = createCallbackFcn(app, @SliderValueChanging, true);
            app.Slider.FontSize = 10;
            app.Slider.Layout.Row = 1;
            app.Slider.Layout.Column = [3 5];

            % Create NextFrameButton
            app.NextFrameButton = uibutton(app.TopGrid, 'push');
            app.NextFrameButton.ButtonPushedFcn = createCallbackFcn(app, @NextFrameButtonPushed, true);
            app.NextFrameButton.Layout.Row = 1;
            app.NextFrameButton.Layout.Column = 6;
            app.NextFrameButton.Text = '>';

            % Create ImageViewLabel
            app.ImageViewLabel = uilabel(app.GridLayoutLeft);
            app.ImageViewLabel.FontWeight = 'bold';
            app.ImageViewLabel.Layout.Row = 3;
            app.ImageViewLabel.Layout.Column = 1;
            app.ImageViewLabel.Text = 'Image View';

            % Create StatusLabel
            app.StatusLabel = uilabel(app.GridLayoutLeft);
            app.StatusLabel.Layout.Row = 2;
            app.StatusLabel.Layout.Column = 1;
            app.StatusLabel.Text = 'Status';

            % Create DataFolderContextLabel
            app.DataFolderContextLabel = uilabel(app.GridLayoutLeft);
            app.DataFolderContextLabel.Layout.Row = 1;
            app.DataFolderContextLabel.Layout.Column = 1;
            app.DataFolderContextLabel.Text = '';

            % Create GridLayoutRight
            app.GridLayoutRight = uigridlayout(app.GridLayoutSub1);
            app.GridLayoutRight.ColumnWidth = {'1x'};
            app.GridLayoutRight.RowHeight = {20, '1x', 20, '1x', 20, 40, 100};
            app.GridLayoutRight.RowSpacing = 5;
            app.GridLayoutRight.Layout.Row = 1;
            app.GridLayoutRight.Layout.Column = 2;

            % Create PlotAxes
            app.PlotAxes = uiaxes(app.GridLayoutRight);
            xlabel(app.PlotAxes, 'Time (s)')
            ylabel(app.PlotAxes, 'Amplitude')
            app.PlotAxes.Layout.Row = 2;
            app.PlotAxes.Layout.Column = 1;

            % Create TemporalProfileLabel
            app.TemporalProfileLabel = uilabel(app.GridLayoutRight);
            app.TemporalProfileLabel.FontWeight = 'bold';
            app.TemporalProfileLabel.Layout.Row = 1;
            app.TemporalProfileLabel.Layout.Column = 1;
            app.TemporalProfileLabel.Text = 'Temporal Profile';

            % Create ROItableLabel
            app.ROItableLabel = uilabel(app.GridLayoutRight);
            app.ROItableLabel.FontWeight = 'bold';
            app.ROItableLabel.Layout.Row = 3;
            app.ROItableLabel.Layout.Column = 1;
            app.ROItableLabel.Text = 'ROI table';

            % Create UITable
            app.UITable = uitable(app.GridLayoutRight);
            app.UITable.ColumnName = {'Column 1'; 'Column 2'; 'Column 3'; 'Column 4'};
            app.UITable.RowName = {};
            app.UITable.Layout.Row = 4;
            app.UITable.Layout.Column = 1;

            % Create EventsGrid
            app.EventsGrid = uigridlayout(app.GridLayoutRight);
            app.EventsGrid.ColumnWidth = {80, '1x', 80, 60};
            app.EventsGrid.Layout.Row = 7;
            app.EventsGrid.Layout.Column = 1;

            % Create ConditionDropDownLabel
            app.ConditionDropDownLabel = uilabel(app.EventsGrid);
            app.ConditionDropDownLabel.HorizontalAlignment = 'center';
            app.ConditionDropDownLabel.Layout.Row = 1;
            app.ConditionDropDownLabel.Layout.Column = 1;
            app.ConditionDropDownLabel.Text = 'Condition';

            % Create ConditionDropDown
            app.ConditionDropDown = uidropdown(app.EventsGrid);
            app.ConditionDropDown.ValueChangedFcn = createCallbackFcn(app, @ConditionDropDownValueChanged, true);
            app.ConditionDropDown.Layout.Row = 1;
            app.ConditionDropDown.Layout.Column = 2;

            % Create RepetitionDropDownLabel
            app.RepetitionDropDownLabel = uilabel(app.EventsGrid);
            app.RepetitionDropDownLabel.HorizontalAlignment = 'center';
            app.RepetitionDropDownLabel.Layout.Row = 2;
            app.RepetitionDropDownLabel.Layout.Column = 1;
            app.RepetitionDropDownLabel.Text = 'Repetition';

            % Create RepetitionDropDown
            app.RepetitionDropDown = uidropdown(app.EventsGrid);
            app.RepetitionDropDown.ValueChangedFcn = createCallbackFcn(app, @RepetitionDropDownValueChanged, true);
            app.RepetitionDropDown.Layout.Row = 2;
            app.RepetitionDropDown.Layout.Column = 2;

            % Create DeleteConditionButton
            app.DeleteConditionButton = uibutton(app.EventsGrid, 'push');
            app.DeleteConditionButton.ButtonPushedFcn = createCallbackFcn(app, @DeleteConditionButtonPushed, true);
            app.DeleteConditionButton.Layout.Row = 1;
            app.DeleteConditionButton.Layout.Column = 3;
            app.DeleteConditionButton.Text = 'Delete';

            % Create DeleteRepetitionButton
            app.DeleteRepetitionButton = uibutton(app.EventsGrid, 'push');
            app.DeleteRepetitionButton.ButtonPushedFcn = createCallbackFcn(app, @DeleteRepetitionButtonPushed, true);
            app.DeleteRepetitionButton.Layout.Row = 2;
            app.DeleteRepetitionButton.Layout.Column = 3;
            app.DeleteRepetitionButton.Text = 'Delete';

            % Create RestoreButton
            app.RestoreButton = uibutton(app.EventsGrid, 'push');
            app.RestoreButton.ButtonPushedFcn = createCallbackFcn(app, @RestoreButtonPushed, true);
            app.RestoreButton.Layout.Row = [1 2];
            app.RestoreButton.Layout.Column = 4;
            app.RestoreButton.Text = 'Restore';

            % Create EventsLabel
            app.EventsLabel = uilabel(app.GridLayoutRight);
            app.EventsLabel.FontWeight = 'bold';
            app.EventsLabel.Layout.Row = 5;
            app.EventsLabel.Layout.Column = 1;
            app.EventsLabel.Text = 'Events ';

            % Create EventSwitchGrid
            app.EventSwitchGrid = uigridlayout(app.GridLayoutRight);
            app.EventSwitchGrid.ColumnWidth = {200, '1x', '1x'};
            app.EventSwitchGrid.RowHeight = {'1x'};
            app.EventSwitchGrid.Layout.Row = 6;
            app.EventSwitchGrid.Layout.Column = 1;

            % Create Switch
            app.Switch = uiswitch(app.EventSwitchGrid, 'slider');
            app.Switch.Items = {'Normal mode', 'Event mode'};
            app.Switch.ValueChangedFcn = createCallbackFcn(app, @SwitchValueChanged, true);
            app.Switch.FontWeight = 'bold';
            app.Switch.Layout.Row = 1;
            app.Switch.Layout.Column = 1;
            app.Switch.Value = 'Normal mode';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = DataViewer_exported(varargin)

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @(app)startupFcn(app, varargin{:}))

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end