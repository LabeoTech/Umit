classdef DataViewer_EventsManager_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                      matlab.ui.Figure
        GridLayout                    matlab.ui.container.GridLayout
        DetectionSummaryPanel         matlab.ui.container.Panel
        GridDetectionSummary          matlab.ui.container.GridLayout
        EditConditionNamesLabel       matlab.ui.control.Label
        RestoreConditionNamesButton   matlab.ui.control.Button
        EditConditionNamesButton      matlab.ui.control.Button
        ConditionSummaryTable         matlab.ui.control.Table
        DetectionSummaryTextArea      matlab.ui.control.TextArea
        TrialSplittingPanel           matlab.ui.container.Panel
        GridTrialSplitting            matlab.ui.container.GridLayout
        SetBaselinePeriodButton       matlab.ui.control.Button
        AutoBaselineButton            matlab.ui.control.Button
        BaselinesEditField            matlab.ui.control.NumericEditField
        BaselinesEditFieldLabel       matlab.ui.control.Label
        GridFooter                    matlab.ui.container.GridLayout
        StatusLabel                   matlab.ui.control.Label
        CloseButton                   matlab.ui.control.Button
        SaveeventsmatButton           matlab.ui.control.Button
        DetectButton                  matlab.ui.control.Button
        ConditionFilePanel            matlab.ui.container.Panel
        GridConditionFile             matlab.ui.container.GridLayout
        ChoosefileButton              matlab.ui.control.Button
        UpdateConditionsButton        matlab.ui.control.Button
        CSVColumnsTable               matlab.ui.control.Table
        CSVColumnsLabel               matlab.ui.control.Label
        ConditionsourceDropDown       matlab.ui.control.DropDown
        ConditionsourceDropDownLabel  matlab.ui.control.Label
        PreviewPanel                  matlab.ui.container.Panel
        GridPreviewPanel              matlab.ui.container.GridLayout
        PlotAllButton                 matlab.ui.control.Button
        UIAxes                        matlab.ui.control.UIAxes
        TriggerChannelPanel           matlab.ui.container.Panel
        TriggerChannelGrid            matlab.ui.container.GridLayout
        EmptyLabel                    matlab.ui.control.Label
        DetectionParametersPanel      matlab.ui.container.Panel
        GridDetectionParameters       matlab.ui.container.GridLayout
        ThresholdEditField            matlab.ui.control.EditField
        ThresholdEditFieldLabel       matlab.ui.control.Label
        TriggertypeDropDown           matlab.ui.control.DropDown
        TriggertypeDropDownLabel      matlab.ui.control.Label
        PolarityDropDown              matlab.ui.control.DropDown
        PolarityLabel                 matlab.ui.control.Label
        MininterstimSpinner           matlab.ui.control.Spinner
        MininterstimLabel             matlab.ui.control.Label
        LowpassHzEditField            matlab.ui.control.NumericEditField
        LowpassHzLabel                matlab.ui.control.Label
        AnalogInputsPanel             matlab.ui.container.Panel
        GridAnalogInputs              matlab.ui.container.GridLayout
        SourceStatusLabel             matlab.ui.control.Label
        LoadAnalogInputsButton        matlab.ui.control.Button
    end


    properties (Access = private)
        ParentApp = []
        SaveFolder char = ''
        RawFolder char = ''
        evObj = []

        TriggerChannelCheckBoxes matlab.ui.control.CheckBox = matlab.ui.control.CheckBox.empty
        PreviewAxes matlab.ui.control.UIAxes = matlab.ui.control.UIAxes.empty

        ConditionFilePath char = ''
        HasUnsavedChanges logical = false
        OriginalEventNameList cell = {}
        IsEditingConditionNames logical = false
        ConditionNameColumnStyle matlab.ui.style.Style = matlab.ui.style.Style.empty

    end

    methods (Access = private)

        function configureInitialGUIState(app)
            %CONFIGUREINITIALGUISTATE Set static defaults and disabled states.

            app.StatusLabel.Text = 'Ready.';
            app.SourceStatusLabel.Text = 'Analog inputs not loaded.';

            app.ThresholdEditField.Value = 'auto';

            app.TriggertypeDropDown.Items = {'EdgeSet', 'EdgeToggle'};
            app.TriggertypeDropDown.Value = 'EdgeSet';

            app.PolarityDropDown.Items = {'positive', 'negative'};
            app.PolarityDropDown.Value = 'positive';

            app.MininterstimSpinner.Limits = [eps Inf];
            app.MininterstimSpinner.Value = 2;

            app.LowpassHzEditField.Limits = [0 Inf];
            app.LowpassHzEditField.Value = 0;

            app.BaselinesEditField.Limits = [eps Inf];
            app.BaselinesEditField.Value = 1;

            app.ConditionsourceDropDown.Items = {'none', 'csv', 'vpixx'};
            app.ConditionsourceDropDown.Value = 'none';

            app.CSVColumnsTable.Data = table(false(0, 1), strings(0, 1), ...
                'VariableNames', {'Use', 'Column'});
            app.CSVColumnsTable.ColumnName = {'Use', 'Column'};
            app.CSVColumnsTable.ColumnEditable = [true false];
            app.CSVColumnsTable.ColumnWidth = {40,'1x'};
            app.CSVColumnsTable.RowName = {};

            app.ConditionSummaryTable.Data = table(strings(0, 1), zeros(0, 1), zeros(0, 1), ...
                'VariableNames', {'ConditionName', 'NRepeats', 'EventDuration_s'});
            app.ConditionSummaryTable.ColumnName = {'ConditionName', 'NRepeats', 'EventLength_s'};
            app.ConditionSummaryTable.ColumnEditable = [false false false];
            app.ConditionSummaryTable.RowName = {};

            app.DetectionSummaryTextArea.Editable = 'off';
            app.DetectionSummaryTextArea.Value = {'No events detected.'};

            app.SaveeventsmatButton.Enable = 'off';
            app.DetectButton.Enable = 'off';
            app.ChoosefileButton.Enable = 'off';
            app.UpdateConditionsButton.Enable = 'off';
            app.CSVColumnsTable.Enable = 'off';
            app.CSVColumnsLabel.Enable = 'off';
            app.SetBaselinePeriodButton.Enable = 'off';
            app.AutoBaselineButton.Enable = 'off';
            cla(app.UIAxes);
            title(app.UIAxes, 'Analog trigger preview');
            xlabel(app.UIAxes, 'Time (s)');
            ylabel(app.UIAxes, 'Amplitude (V)');
            grid(app.UIAxes, 'on');

            app.PreviewAxes = app.UIAxes;

            app.EditConditionNamesButton.Enable = 'off';
            app.RestoreConditionNamesButton.Text = 'Restore';
            app.RestoreConditionNamesButton.Enable = 'off';

        end
        function createOrRefreshEventsObject(app)
            %CREATEORREFRESHEVENTSOBJECT Create local EventsManager instance.
            %
            %   If events.mat exists, EventsManager loads it during construction.
            %   In that case, synchronize the GUI from the loaded object so the
            %   startup state reflects existing event metadata.

            constructorSaveFolder = app.SaveFolder;
            if isempty(constructorSaveFolder) || ~isfolder(constructorSaveFolder)
                constructorSaveFolder = pwd;
            end

            if isempty(app.RawFolder) || ~isfolder(app.RawFolder)
                app.evObj = EventsManager(constructorSaveFolder);
            else
                app.evObj = EventsManager(constructorSaveFolder, app.RawFolder, 'none');

                % EventsManager may return early after loading events.mat, before
                % applying the constructor RawFolder input.
                app.evObj.RawFolder = app.RawFolder;
            end

            if isempty(app.RawFolder) && ~isempty(app.evObj.RawFolder)
                app.RawFolder = app.evObj.RawFolder;
            end

            app.syncGUIFromEventsObject();
            app.storeOriginalEventNameList();
            app.refreshSourceStatusLabel();
            app.refreshDetectionSummary();
            app.refreshConditionSummaryTable();
            app.refreshConditionNameEditButtonState();

            if ~isempty(app.evObj) && ~isempty(app.evObj.AnalogIN)
                app.configureLimitsAfterAnalogLoad();
                app.populateTriggerChannelCheckboxes();
                app.DetectButton.Enable = 'on';
                app.refreshPreviewAxes();
                app.SaveeventsmatButton.Enable = 'off';
                app.AutoBaselineButton.Enable = 'on';
                app.SetBaselinePeriodButton.Enable = 'off';
            end
        end

        function configureLimitsAfterAnalogLoad(app)
            %CONFIGURELIMITSAFTERANALOGLOAD Fine-tune limits using loaded AnalogIN.

            if isempty(app.evObj) || isempty(app.evObj.AnalogIN)
                return
            end

            nyquistHz = double(app.evObj.sr) / 2;

            app.LowpassHzEditField.Limits = [0 nyquistHz];
            app.LowpassHzEditField.Value = min(max(app.LowpassHzEditField.Value, 0), nyquistHz);

            app.MininterstimSpinner.Limits = [eps Inf];
            app.MininterstimSpinner.Value = max(app.MininterstimSpinner.Value, eps);

            app.BaselinesEditField.Limits = [eps Inf];
            app.BaselinesEditField.Value = max(app.BaselinesEditField.Value, eps);
        end

        function populateTriggerChannelCheckboxes(app)
            %POPULATETRIGGERCHANNELCHECKBOXES Replace placeholder with channel checkboxes.

            delete(app.TriggerChannelGrid.Children);
            app.TriggerChannelCheckBoxes = matlab.ui.control.CheckBox.empty;

            if isempty(app.evObj) || isempty(app.evObj.AIChanList)
                app.TriggerChannelGrid.RowHeight = {'1x'};
                app.EmptyLabel = uilabel(app.TriggerChannelGrid);
                app.EmptyLabel.Text = 'Empty';
                app.EmptyLabel.HorizontalAlignment = 'center';
                app.EmptyLabel.FontAngle = 'italic';
                app.EmptyLabel.Layout.Row = 1;
                app.EmptyLabel.Layout.Column = 1;
                return
            end

            chanList = setdiff(app.evObj.AIChanList, {'CameraTrig', 'CameraTrig2'}, 'stable');

            if isempty(chanList)
                chanList = app.evObj.AIChanList;
            end

            app.TriggerChannelGrid.RowHeight = repmat({24}, 1, numel(chanList));
            app.TriggerChannelGrid.ColumnWidth = {'1x'};

            selectedNames = app.evObj.trigChanName;
            if isempty(selectedNames) || isempty([selectedNames{:}])
                selectedNames = chanList(1);
            end

            for iChan = 1:numel(chanList)
                cb = uicheckbox(app.TriggerChannelGrid);
                cb.Text = chanList{iChan};
                cb.Value = any(strcmpi(chanList{iChan}, selectedNames));
                cb.ValueChangedFcn = @(src, evt) app.onTriggerChannelCheckboxChanged(src);
                cb.Layout.Row = iChan;
                cb.Layout.Column = 1;

                app.TriggerChannelCheckBoxes(end+1) = cb;
            end

            if ~any([app.TriggerChannelCheckBoxes.Value])
                app.TriggerChannelCheckBoxes(1).Value = true;
            end
        end

        function onTriggerChannelCheckboxChanged(app, src)
            %ONTRIGGERCHANNELCHECKBOXCHANGED Prevent empty channel selection.

            if isempty(app.TriggerChannelCheckBoxes)
                return
            end

            if ~any([app.TriggerChannelCheckBoxes.Value])
                src.Value = true;
                return
            end

            app.refreshPreviewAxes();
        end

        function selectedChannels = getSelectedTriggerChannels(app)
            %GETSELECTEDTRIGGERCHANNELS Return checked trigger channel names.

            selectedChannels = {};

            if isempty(app.TriggerChannelCheckBoxes)
                return
            end

            checked = [app.TriggerChannelCheckBoxes.Value];

            if any(checked)
                selectedChannels = {app.TriggerChannelCheckBoxes(checked).Text};
            end
        end

        function trigThr = parseThresholdFromGUI(app)
            %PARSETHRESHOLDFROMGUI Return threshold as 'auto' or finite numeric scalar.

            rawText = strtrim(app.ThresholdEditField.Value);

            if isempty(rawText) || strcmpi(rawText, 'auto')
                trigThr = 'auto';
                app.ThresholdEditField.Value = 'auto';
                return
            end

            trigThr = str2double(rawText);

            if ~isscalar(trigThr) || ~isfinite(trigThr)
                error('DataViewerEventsManager:InvalidThreshold', ...
                    'Threshold must be "auto" or a finite numeric value.');
            end
        end

        function applyConditionFileToDetectedEvents(app)
            %APPLYCONDITIONFILETODETECTEDEVENTS Apply csv/vpixx condition metadata.

            sourceType = lower(app.ConditionsourceDropDown.Value);

            if strcmpi(sourceType, 'none')
                return
            end

            if isempty(app.ConditionFilePath) || ~isfile(app.ConditionFilePath)
                error('DataViewerEventsManager:MissingConditionFile', ...
                    'Select a valid condition file first.');
            end

            app.evObj.EventFileParseMethod = sourceType;

            csvCols = {''};

            if strcmpi(sourceType, 'csv')
                T = app.CSVColumnsTable.Data;

                if istable(T) && all(ismember({'Use', 'Column'}, T.Properties.VariableNames))
                    useMask = logical(T.Use);
                    csvCols = cellstr(string(T.Column(useMask)));
                else
                    csvCols = {};
                end

                if isempty(csvCols)
                    error('DataViewerEventsManager:NoCSVColumns', ...
                        'Select at least one CSV column.');
                end
            end

            bOK = app.evObj.readConditionFile(app.ConditionFilePath, 'CSVcols', csvCols);

            if ~bOK
                error('DataViewerEventsManager:ConditionFileFailed', ...
                    'EventsManager failed to read the condition file.');
            end

            app.UpdateConditionsButton.Enable = 'on';
        end

        function refreshPreviewAxes(app)
            %REFRESHPREVIEWAXES Plot selected trigger channels in embedded axes.

            if isempty(app.evObj) || isempty(app.evObj.AnalogIN)
                cla(app.UIAxes);
                title(app.UIAxes, 'Analog trigger preview');
                return
            end

            selectedChannels = app.getSelectedTriggerChannels();

            if isempty(selectedChannels) && ~isempty(app.evObj.trigChanName) && ...
                    ~isempty([app.evObj.trigChanName{:}])
                selectedChannels = app.evObj.trigChanName;
            end

            if isempty(selectedChannels)
                selectedChannels = app.evObj.AIChanList(1);
            end

            app.rebuildPreviewAxes(numel(selectedChannels));
            app.evObj.plot(selectedChannels, app.PreviewAxes);
        end

        function rebuildPreviewAxes(app, nAxes)
            %REBUILDPREVIEWAXES Create one preview axes per selected trigger channel.

            nAxes = max(1, round(nAxes));

            delete(findall(app.GridPreviewPanel, 'Type', 'axes'));

            app.GridPreviewPanel.RowHeight = [{30}, repmat({'1x'}, 1, nAxes)];
            app.GridPreviewPanel.ColumnWidth = {80, '1x'};

            app.PreviewAxes = matlab.ui.control.UIAxes.empty;

            for iAxes = 1:nAxes
                ax = uiaxes(app.GridPreviewPanel);
                ax.Layout.Row = iAxes + 1;
                ax.Layout.Column = [1 2];
                ax.Box = 'on';
                grid(ax, 'on');

                app.PreviewAxes(end+1) = ax;
            end
        end



        function refreshSourceStatusLabel(app)
            %REFRESHSOURCESTATUSLABEL Update compact Analog Inputs status.

            if isempty(app.evObj) || isempty(app.evObj.AnalogIN)
                nFiles = numel(dir(fullfile(app.RawFolder, 'ai_*.bin')));
                app.SourceStatusLabel.Text = sprintf('Raw: %s | AI files: %d | Analog inputs not loaded', ...
                    app.safeFolderText(app.RawFolder), nFiles);
                return
            end

            nFiles = numel(dir(fullfile(app.RawFolder, 'ai_*.bin')));
            nChan = size(app.evObj.AnalogIN, 2);
            recLength = size(app.evObj.AnalogIN, 1) / double(app.evObj.sr);

            app.SourceStatusLabel.Text = sprintf( ...
                'Raw: %s | AI files: %d | channels: %d | AI sr: %.3g Hz | recording: %.3f s', ...
                app.safeFolderText(app.RawFolder), nFiles, nChan, double(app.evObj.sr), recLength);

            app.SourceStatusLabel.Tooltip = app.SourceStatusLabel.Text;
        end
        function refreshDetectionSummary(app)
            %REFRESHDETECTIONSUMMARY Update compact detection diagnostic text.

            eventsPath = fullfile(app.SaveFolder, 'events.mat');
            if ~isempty(app.SaveFolder) && isfile(eventsPath)
                eventsText = 'yes';
            else
                eventsText = 'no';
            end

            if isempty(app.evObj)
                app.DetectionSummaryTextArea.Value = { ...
                    sprintf('events.mat exists: %s', eventsText); ...
                    'Detected events: 0'; ...
                    'Conditions: 0'};
                return
            end

            nEvents = 0;
            nConditions = 0;

            if ~isempty(app.evObj.state)
                nEvents = nnz(app.evObj.state);
            end

            if ~isempty(app.evObj.eventNameList)
                nConditions = numel(app.evObj.eventNameList);
            end

            if isnumeric(app.evObj.trigThr)
                thresholdText = sprintf('%.6g V', double(app.evObj.trigThr));
            else
                thresholdText = char(string(app.evObj.trigThr));
            end

            if isempty(app.evObj.baselinePeriod)
                baselineText = 'not set';
            else
                baselineText = sprintf('%.6g s', double(app.evObj.baselinePeriod));
            end

            if isempty(app.ConditionFilePath)
                condFileText = 'none';
            else
                [~, f, e] = fileparts(app.ConditionFilePath);
                condFileText = [f e];
            end

            app.DetectionSummaryTextArea.Value = { ...
                sprintf('events.mat exists: %s', eventsText); ...
                sprintf('Detected events: %d', nEvents); ...
                sprintf('Conditions: %d', nConditions); ...
                sprintf('Threshold: %s', thresholdText); ...
                sprintf('Trigger type: %s', app.evObj.trigType); ...
                sprintf('Polarity: %s', app.evObj.trigPolarity); ...
                sprintf('Baseline period: %s', baselineText); ...
                sprintf('Condition source: %s', app.ConditionsourceDropDown.Value); ...
                sprintf('Condition file: %s', condFileText)};
        end
        function refreshConditionSummaryTable(app)
            %REFRESHCONDITIONSUMMARYTABLE Summarize detected conditions.

            T = table(strings(0, 1), zeros(0, 1), zeros(0, 1), ...
                'VariableNames', {'ConditionName', 'NRepeats', 'EventLength_s'});

            if isempty(app.evObj) || isempty(app.evObj.eventID) || isempty(app.evObj.state)
                app.ConditionSummaryTable.Data = T;
                return
            end

            condIDs = unique(double(app.evObj.eventID(app.evObj.state == 1)), 'stable');
            nCond = numel(condIDs);

            conditionName = strings(nCond, 1);
            nRepeats = zeros(nCond, 1);
            EventLength_s = nan(nCond, 1);

            for iCond = 1:nCond
                condID = condIDs(iCond);

                if condID >= 1 && condID <= numel(app.evObj.eventNameList)
                    conditionName(iCond) = string(app.evObj.eventNameList{condID});
                else
                    conditionName(iCond) = sprintf('Event %d', condID);
                end

                onTimes = double(app.evObj.timestamps(app.evObj.eventID == condID & app.evObj.state == 1));
                offTimes = double(app.evObj.timestamps(app.evObj.eventID == condID & app.evObj.state == 0));

                nPairs = min(numel(onTimes), numel(offTimes));
                nRepeats(iCond) = nPairs;

                if nPairs > 0
                    dur = offTimes(1:nPairs) - onTimes(1:nPairs);
                    dur = dur(isfinite(dur) & dur > 0);
                    if ~isempty(dur)
                        EventLength_s(iCond) = mean(dur, 'omitnan');
                    end
                end
            end

            app.ConditionSummaryTable.Data = table(conditionName, nRepeats, EventLength_s, ...
                'VariableNames', {'ConditionName', 'NRepeats', 'EventLength_s'});
        end



        function closeEventsManagerApp(app)
            %CLOSEEVENTSMANAGERAPP Close app, optionally handling unsaved changes.

            if app.HasUnsavedChanges
                answer = uiconfirm(app.UIFigure, ...
                    'Detected events have not been saved. Save before closing?', ...
                    'Unsaved event changes', ...
                    'Options', {'Save', 'Discard', 'Cancel'}, ...
                    'DefaultOption', 'Save', ...
                    'CancelOption', 'Cancel', ...
                    'Icon', 'warning');

                switch answer
                    case 'Save'
                        app.SaveeventsmatButtonPushed([]);
                        if app.HasUnsavedChanges
                            return
                        end

                    case 'Cancel'
                        return

                    case 'Discard'
                        % Close without saving.
                end
            end

            delete(app);
        end

        function refreshParentDataViewerEvents(app)
            %REFRESHPARENTDATAVIEWEREVENTS Reload events in parent DataViewer after save.

            if isempty(app.ParentApp) || ~isvalid(app.ParentApp)
                return
            end

            try
                if ismethod(app.ParentApp, 'initializeEventsForLoadedData')
                    app.ParentApp.initializeEventsForLoadedData('.dat');
                end

                if ismethod(app.ParentApp, 'refreshViewerAfterLoad')
                    app.ParentApp.refreshViewerAfterLoad();
                end

                if ismethod(app.ParentApp, 'setStatusMessage')
                    app.ParentApp.setStatusMessage('events.mat saved and event metadata reloaded.');
                end
            catch ME
                warning('DataViewerEventsManager:ParentRefreshFailed', ...
                    'Failed to refresh parent DataViewer events.\n%s', ME.message);
            end
        end


        function txt = safeFolderText(app, folderPath) %#ok<INUSD>
            %SAFEFOLDERTEXT Return compact folder text.

            if isempty(folderPath)
                txt = '<not set>';
                return
            end

            txt = char(string(folderPath));
        end

        function setLocalStatus(app, msg)
            %SETLOCALSTATUS Update footer status text.

            app.StatusLabel.Text = char(string(msg));
            drawnow limitrate
        end

        function applyDetectionParametersToEventsObject(app)
            %APPLYDETECTIONPARAMETERSTOEVENTSOBJECT Push GUI values into EventsManager.

            selectedChannels = app.getSelectedTriggerChannels();
            if isempty(selectedChannels)
                error('DataViewerEventsManager:NoTriggerChannel', ...
                    'Select at least one trigger channel.');
            end

            app.evObj.trigChanName = selectedChannels;
            app.evObj.trigThr = app.parseThresholdFromGUI();
            app.evObj.trigType = app.TriggertypeDropDown.Value;
            app.evObj.trigPolarity = lower(app.PolarityDropDown.Value);
            app.evObj.minInterStim = app.MininterstimSpinner.Value;

            if isempty(app.evObj.AnalogIN)
                error('DataViewerEventsManager:NoAnalogInput', ...
                    'Load analog inputs before detecting events.');
            end

            nyquistHz = double(app.evObj.sr) / 2;
            filterHz = app.LowpassHzEditField.Value;

            if filterHz < 0 || filterHz >= nyquistHz
                error('DataViewerEventsManager:InvalidLowPass', ...
                    'Low-pass cutoff must be >= 0 and smaller than Nyquist (%0.3f Hz).', ...
                    nyquistHz);
            end
        end

        function syncGUIFromEventsObject(app)
            %SYNCGUIFROMEVENTSOBJECT Update controls from current EventsManager state.

            if isempty(app.evObj)
                return
            end

            % Detection parameters.
            if isnumeric(app.evObj.trigThr) && isscalar(app.evObj.trigThr)
                app.ThresholdEditField.Value = sprintf('%.6g', double(app.evObj.trigThr));
            else
                app.ThresholdEditField.Value = char(string(app.evObj.trigThr));
            end

            if any(strcmpi(app.TriggertypeDropDown.Items, app.evObj.trigType))
                app.TriggertypeDropDown.Value = app.TriggertypeDropDown.Items{ ...
                    find(strcmpi(app.TriggertypeDropDown.Items, app.evObj.trigType), 1, 'first')};
            end

            if any(strcmpi(app.PolarityDropDown.Items, app.evObj.trigPolarity))
                app.PolarityDropDown.Value = app.PolarityDropDown.Items{ ...
                    find(strcmpi(app.PolarityDropDown.Items, app.evObj.trigPolarity), 1, 'first')};
            end

            if ~isempty(app.evObj.minInterStim)
                app.MininterstimSpinner.Value = max(double(app.evObj.minInterStim), eps);
            end

            if ~isempty(app.evObj.baselinePeriod)
                app.BaselinesEditField.Value = max(double(app.evObj.baselinePeriod), eps);
            end

            if ~isempty(app.evObj.EventFileParseMethod) && ...
                    any(strcmpi(app.ConditionsourceDropDown.Items, app.evObj.EventFileParseMethod))
                app.ConditionsourceDropDown.Value = lower(app.evObj.EventFileParseMethod);
            else
                app.ConditionsourceDropDown.Value = 'none';
            end

            % Condition file controls.
            if ~isempty(app.evObj.EventFileName)
                app.ConditionFilePath = app.evObj.EventFileName;
                [~, fileName, fileExt] = fileparts(app.ConditionFilePath);
                app.ChoosefileButton.Text = [fileName fileExt];
                app.ChoosefileButton.Tooltip = app.ConditionFilePath;
            else
                app.ConditionFilePath = '';
                app.ChoosefileButton.Text = 'Choose file...';
                app.ChoosefileButton.Tooltip = '';
            end

            app.ConditionsourceDropDownValueChanged([]);

            if ~isempty(app.ConditionFilePath) && isfile(app.ConditionFilePath) && ...
                    ~strcmpi(app.ConditionsourceDropDown.Value, 'none') && ...
                    ~isempty(app.evObj.timestamps)
                app.UpdateConditionsButton.Enable = 'on';
            end

            % Save state.
            app.HasUnsavedChanges = false;
            app.SaveeventsmatButton.Enable = 'off';
        end

        function storeOriginalEventNameList(app)
            %STOREORIGINALEVENTNAMELIST Store current condition names for restoration.

            if isempty(app.evObj) || isempty(app.evObj.eventNameList)
                app.OriginalEventNameList = {};
                app.RestoreConditionNamesButton.Text = 'Restore';
                app.RestoreConditionNamesButton.Enable = 'off';
                return
            end

            app.OriginalEventNameList = app.evObj.eventNameList(:).';
            app.RestoreConditionNamesButton.Text = 'Restore';
            app.RestoreConditionNamesButton.Enable = 'off';
        end

        function enterConditionNameEditMode(app)
            %ENTERCONDITIONNAMEEDITMODE Enable editing of condition names.

            if isempty(app.ConditionSummaryTable.Data) || ...
                    ~istable(app.ConditionSummaryTable.Data) || ...
                    ~ismember('ConditionName', app.ConditionSummaryTable.Data.Properties.VariableNames)
                uialert(app.UIFigure, ...
                    'Condition summary table does not contain editable condition names.', ...
                    'Cannot edit condition names', ...
                    'Icon', 'warning');
                return
            end

            app.IsEditingConditionNames = true;

            app.ConditionSummaryTable.ColumnEditable = [true false false];
            app.EditConditionNamesButton.Text = 'Apply';
            app.RestoreConditionNamesButton.Text = 'Cancel';
            app.RestoreConditionNamesButton.Enable = 'on';

            % Block operations that could conflict with pending name edits.
            app.SaveeventsmatButton.Enable = 'off';
            app.DetectButton.Enable = 'off';
            app.UpdateConditionsButton.Enable = 'off';

            try
                removeStyle(app.ConditionSummaryTable);
            catch
            end

            style = uistyle('BackgroundColor', [1.0 0.96 0.65]);
            addStyle(app.ConditionSummaryTable, style, 'column', 1);
            app.ConditionNameColumnStyle = style;

            app.setLocalStatus('Edit condition names, then click Apply. Click Cancel to discard edits.');
        end
        function applyEditedConditionNames(app)
            %APPLYEDITEDCONDITIONNAMES Validate and apply edited condition names.

            try
                T = app.ConditionSummaryTable.Data;

                if ~istable(T) || ~ismember('ConditionName', T.Properties.VariableNames)
                    error('Condition summary table does not contain ConditionName.');
                end

                newNames = strtrim(string(T.ConditionName));
                app.evObj.setEventNameList(newNames);

                app.IsEditingConditionNames = false;
                app.ConditionSummaryTable.ColumnEditable = [false false false];
                app.EditConditionNamesButton.Text = 'Edit';
                app.RestoreConditionNamesButton.Text = 'Restore';

                try
                    removeStyle(app.ConditionSummaryTable);
                catch
                end
                app.ConditionNameColumnStyle = matlab.ui.style.Style.empty;

                app.HasUnsavedChanges = true;
                app.SaveeventsmatButton.Enable = 'on';
                app.DetectButton.Enable = 'on';

                app.ConditionsourceDropDownValueChanged([]);

                app.refreshDetectionSummary();
                app.refreshConditionSummaryTable();
                app.refreshConditionNameEditButtonState();
                app.refreshPreviewAxes();

                app.setLocalStatus('Condition names applied. Save events.mat to keep this change.');

            catch ME
                uialert(app.UIFigure, ...
                    getReport(ME, 'basic', 'hyperlinks', 'off'), ...
                    'Invalid condition names', ...
                    'Icon', 'warning');
            end
        end

        function exitConditionNameEditMode(app)
            %EXITCONDITIONNAMEEDITMODE Return condition-name table to read-only state.

            app.IsEditingConditionNames = false;

            app.ConditionSummaryTable.ColumnEditable = [false false false];
            app.EditConditionNamesButton.Text = 'Edit';
            app.RestoreConditionNamesButton.Text = 'Restore';

            try
                removeStyle(app.ConditionSummaryTable);
            catch
            end

            app.ConditionNameColumnStyle = matlab.ui.style.Style.empty;

            app.DetectButton.Enable = 'on';

            % Restore condition-file button state from the current source/file state.
            app.ConditionsourceDropDownValueChanged([]);

            app.refreshConditionNameEditButtonState();
        end

        function refreshConditionRestoreButtonState(app)
            %REFRESHCONDITIONRESTOREBUTTONSTATE Enable Restore when names differ.

            if isempty(app.RestoreConditionNamesButton) || ~isvalid(app.RestoreConditionNamesButton)
                return
            end

            if app.IsEditingConditionNames
                app.RestoreConditionNamesButton.Text = 'Cancel';
                app.RestoreConditionNamesButton.Enable = 'on';
                return
            end

            app.RestoreConditionNamesButton.Text = 'Restore';

            if isempty(app.evObj) || isempty(app.evObj.eventNameList) || ...
                    isempty(app.OriginalEventNameList)
                app.RestoreConditionNamesButton.Enable = 'off';
                return
            end

            currentNames = string(app.evObj.eventNameList(:));
            originalNames = string(app.OriginalEventNameList(:));

            if numel(currentNames) ~= numel(originalNames)
                app.RestoreConditionNamesButton.Enable = 'off';
                return
            end

            if ~isequal(currentNames, originalNames)
                app.RestoreConditionNamesButton.Enable = 'on';
            else
                app.RestoreConditionNamesButton.Enable = 'off';
            end
        end


        function refreshConditionNameEditButtonState(app)
            %REFRESHCONDITIONNAMEEDITBUTTONSTATE Enable condition-name edit controls.

            if isempty(app.evObj) || isempty(app.evObj.eventNameList)
                app.EditConditionNamesButton.Enable = 'off';
                app.RestoreConditionNamesButton.Text = 'Restore';
                app.RestoreConditionNamesButton.Enable = 'off';
                return
            end

            app.EditConditionNamesButton.Enable = 'on';
            app.refreshConditionRestoreButtonState();
        end
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, parentApp, saveFolder, rawFolder)

            %STARTUPFCN Initialize modal Events Manager utility.
            %
            %   DataViewer_EventsManager_exported(parentApp, saveFolder, rawFolder)
            %   opens a modal utility for detecting triggers from ai_*.bin files and
            %   saving events.mat.

            if nargin < 2
                parentApp = [];
            end

            if nargin < 3 || isempty(saveFolder)
                saveFolder = '';
            end

            if nargin < 4 || isempty(rawFolder)
                rawFolder = '';
            end

            app.ParentApp = parentApp;
            app.SaveFolder = char(string(saveFolder));
            app.RawFolder = char(string(rawFolder));

            app.UIFigure.Visible = 'off';
            app.UIFigure.WindowStyle = 'modal';
            app.UIFigure.CloseRequestFcn = @(src, evt) app.closeEventsManagerApp();

            app.configureInitialGUIState();
            app.createOrRefreshEventsObject();

            if exist('placeAppInsideCaller', 'file') == 2 && ~isempty(parentApp)
                placeAppInsideCaller(parentApp, app.UIFigure, 'center');
            else
                movegui(app.UIFigure, 'center');
            end

            app.UIFigure.Visible = 'on';

        end

        % Button pushed function: LoadAnalogInputsButton
        function LoadAnalogInputsButtonPushed(app, event)

            %LOADANALOGINPUTSBUTTONPUSHED Load ai_*.bin analog input files.

            try
                if isempty(app.RawFolder) || ~isfolder(app.RawFolder)
                    startFolder = app.SaveFolder;
                    if isempty(startFolder) || ~isfolder(startFolder)
                        startFolder = pwd;
                    end

                    selectedFolder = uigetdir(startFolder, ...
                        'Select RawFolder containing ai_*.bin files');

                    if isequal(selectedFolder, 0)
                        app.setLocalStatus('Analog input loading cancelled.');
                        return
                    end

                    app.RawFolder = selectedFolder;
                end

                if isempty(app.evObj) || ~isa(app.evObj, 'EventsManager') || ~isvalid(app.evObj)
                    app.createOrRefreshEventsObject();
                end

                app.evObj.RawFolder = app.RawFolder;
                
                if ~isempty(app.ParentApp) && isvalid(app.ParentApp)
                    app.ParentApp.setRawFolderFromChildTool(app.RawFolder, app.SaveFolder);
                end

                if isempty(app.evObj.AnalogIN)
                    w = uiprogressdlg(app.UIFigure, ...
                        'Title', 'Loading analog inputs', ...
                        'Message', 'Reading ai_*.bin files...', ...
                        'Indeterminate', 'on');

                    cleanupObj = onCleanup(@() delete(w));
                    app.evObj.setAnalogIN();
                end

                if isempty(app.evObj.AnalogIN)
                    uialert(app.UIFigure, ...
                        sprintf('No analog input data were loaded from:\n\n%s', app.RawFolder), ...
                        'Analog input not loaded', ...
                        'Icon', 'warning');
                    app.setLocalStatus('Analog input loading failed.');
                    return
                end

                app.configureLimitsAfterAnalogLoad();
                app.populateTriggerChannelCheckboxes();
                app.refreshSourceStatusLabel();
                app.refreshDetectionSummary();
                app.refreshPreviewAxes();

                app.DetectButton.Enable = 'on';
                app.setLocalStatus('Analog inputs loaded.');

            catch ME
                app.setLocalStatus(sprintf('Analog input loading failed: %s', ME.message));
                uialert(app.UIFigure, getReport(ME, 'basic', 'hyperlinks', 'off'), ...
                    'Analog input loading failed', 'Icon', 'error');
            end

        end

        % Button pushed function: AutoBaselineButton
        function AutoBaselineButtonPushed(app, event)
            %AUTOBASELINEBUTTONPUSHED Estimate baseline period from detected events.

            try
                if isempty(app.evObj) || isempty(app.evObj.timestamps)
                    uialert(app.UIFigure, ...
                        'Detect events before estimating the baseline period.', ...
                        'No detected events', ...
                        'Icon', 'warning');
                    return
                end

                app.evObj.setBaselinePeriod();
                app.BaselinesEditField.Value = double(app.evObj.baselinePeriod);

                app.HasUnsavedChanges = true;
                app.SaveeventsmatButton.Enable = 'on';
                app.SetBaselinePeriodButton.Enable = 'off';

                app.refreshDetectionSummary();
                app.setLocalStatus('Baseline period estimated automatically. Save events.mat to keep this change.');

            catch ME
                uialert(app.UIFigure, getReport(ME, 'basic', 'hyperlinks', 'off'), ...
                    'Auto baseline failed', 'Icon', 'error');
            end

        end

        % Button pushed function: DetectButton
        function DetectButtonPushed(app, event)
            %DETECTBUTTONPUSHED Detect triggers from selected analog input channels.

            try
                app.applyDetectionParametersToEventsObject();

                w = uiprogressdlg(app.UIFigure, ...
                    'Title', 'Detecting events', ...
                    'Message', 'Detecting triggers...', ...
                    'Indeterminate', 'on');

                cleanupObj = onCleanup(@() delete(w)); %#ok<NASGU>

                app.evObj.getTriggers(app.getSelectedTriggerChannels(), false, ...
                    'FilterFreq', app.LowpassHzEditField.Value);

                if isempty(app.evObj.timestamps)
                    app.SaveeventsmatButton.Enable = 'off';
                    app.refreshPreviewAxes();
                    app.refreshDetectionSummary();
                    app.refreshConditionNameEditButtonState();

                    uialert(app.UIFigure, ...
                        'No triggers were detected in the selected channel(s).', ...
                        'Trigger detection failed', ...
                        'Icon', 'warning');

                    app.setLocalStatus('Trigger detection failed.');
                    return
                end

                if isempty(app.evObj.baselinePeriod)
                    app.evObj.setBaselinePeriod();
                end

                app.BaselinesEditField.Value = double(app.evObj.baselinePeriod);

                if ~strcmpi(app.ConditionsourceDropDown.Value, 'none') && ...
                        ~isempty(app.ConditionFilePath)
                    app.applyConditionFileToDetectedEvents();
                end

                app.HasUnsavedChanges = true;
                app.SaveeventsmatButton.Enable = 'on';

                % Detection or condition-file parsing defines the new original names
                % used by Restore. Manual renames do not update this baseline.
                app.storeOriginalEventNameList();

                app.refreshSourceStatusLabel();
                app.refreshPreviewAxes();
                app.refreshDetectionSummary();
                app.refreshConditionSummaryTable();
                app.refreshConditionNameEditButtonState();

                app.AutoBaselineButton.Enable = 'on';
                app.SetBaselinePeriodButton.Enable = 'off';
                app.setLocalStatus('Trigger detection completed.');

            catch ME
                app.setLocalStatus(sprintf('Trigger detection failed: %s', ME.message));
                uialert(app.UIFigure, getReport(ME, 'basic', 'hyperlinks', 'off'), ...
                    'Trigger detection failed', 'Icon', 'error');
            end

        end

        % Value changed function: ConditionsourceDropDown
        function ConditionsourceDropDownValueChanged(app, event)

            %CONDITIONSOURCEDROPDOWNVALUECHANGED Enable condition-file controls.
            if isempty(app.evObj)
                return
            end
            sourceType = lower(app.ConditionsourceDropDown.Value);

            app.evObj.EventFileParseMethod = sourceType;

            switch sourceType
                case 'none'
                    app.ConditionFilePath = '';
                    app.ChoosefileButton.Enable = 'off';
                    app.UpdateConditionsButton.Enable = 'off';
                    app.CSVColumnsTable.Enable = 'off';
                    app.CSVColumnsLabel.Enable = 'off';
                    app.CSVColumnsTable.Data = table(false(0, 1), strings(0, 1), ...
                        'VariableNames', {'Use', 'Column'});
                    app.ChoosefileButton.Text = 'Choose file...';

                case 'csv'
                    app.ChoosefileButton.Enable = 'on';
                    app.UpdateConditionsButton.Enable = 'off';
                    app.CSVColumnsTable.Enable = 'off';
                    app.CSVColumnsLabel.Enable = 'off';

                case 'vpixx'
                    app.ChoosefileButton.Enable = 'on';
                    app.UpdateConditionsButton.Enable = 'off';
                    app.CSVColumnsTable.Enable = 'off';
                    app.CSVColumnsLabel.Enable = 'off';
            end

            app.refreshDetectionSummary();

        end

        % Button pushed function: ChoosefileButton
        function ChoosefileButtonPushed(app, event)

            %CHOOSEFILEBUTTONPUSHED Select CSV or VPixx condition file.

            sourceType = lower(app.ConditionsourceDropDown.Value);

            if strcmpi(sourceType, 'none')
                return
            end

            startFolder = app.RawFolder;
            if isempty(startFolder) || ~isfolder(startFolder)
                startFolder = pwd;
            end

            switch sourceType
                case 'csv'
                    filterSpec = {'*.csv', 'CSV files (*.csv)'};
                case 'vpixx'
                    filterSpec = {'*.vpixx;*.txt', 'VPixx files (*.vpixx, *.txt)'};
                otherwise
                    return
            end

            [fileName, folderName] = uigetfile(filterSpec, ...
                'Select condition file', startFolder);

            if isequal(fileName, 0)
                app.setLocalStatus('Condition file selection cancelled.');
                return
            end

            app.ConditionFilePath = fullfile(folderName, fileName);
            app.ChoosefileButton.Text = fileName;
            app.ChoosefileButton.Tooltip = app.ConditionFilePath;

            if strcmpi(sourceType, 'csv')
                opts = detectImportOptions(app.ConditionFilePath);

                varNames = string(opts.VariableNames(:));
                app.CSVColumnsTable.Data = table(true(numel(varNames), 1), varNames, ...
                    'VariableNames', {'Use', 'Column'});
                app.CSVColumnsTable.ColumnName = {'Use', 'Column'};
                app.CSVColumnsTable.ColumnEditable = [true false];
                app.CSVColumnsTable.Enable = 'on';
                app.CSVColumnsLabel.Enable = 'on';
            else
                app.CSVColumnsTable.Enable = 'off';
                app.CSVColumnsLabel.Enable = 'off';
            end

            if ~isempty(app.evObj) && ~isempty(app.evObj.timestamps)
                app.UpdateConditionsButton.Enable = 'on';
            end

            app.refreshDetectionSummary();
            app.setLocalStatus('Condition file selected.');

        end

        % Button pushed function: UpdateConditionsButton
        function UpdateConditionsButtonPushed(app, event)
            %UPDATECONDITIONSBUTTONPUSHED Apply condition file to detected events.

            try
                if isempty(app.evObj) || isempty(app.evObj.timestamps)
                    uialert(app.UIFigure, ...
                        'Detect events before updating conditions.', ...
                        'No detected events', ...
                        'Icon', 'warning');
                    return
                end

                app.applyConditionFileToDetectedEvents();
                app.HasUnsavedChanges = true;

                % The condition file defines the new original names used by Restore.
                app.storeOriginalEventNameList();

                app.refreshPreviewAxes();
                app.refreshDetectionSummary();
                app.refreshConditionSummaryTable();
                app.refreshConditionNameEditButtonState();

                app.setLocalStatus('Condition list updated.');

            catch ME
                app.setLocalStatus(sprintf('Condition update failed: %s', ME.message));
                uialert(app.UIFigure, getReport(ME, 'basic', 'hyperlinks', 'off'), ...
                    'Condition update failed', 'Icon', 'error');
            end

        end

        % Button pushed function: PlotAllButton
        function PlotAllButtonPushed(app, event)

            %PLOTALLBUTTONPUSHED Plot all AI channels in a separate EventsManager figure.

            try
                if isempty(app.evObj) || isempty(app.evObj.AnalogIN)
                    uialert(app.UIFigure, ...
                        'Load analog inputs before plotting all channels.', ...
                        'No analog inputs', ...
                        'Icon', 'warning');
                    return
                end

                app.evObj.plot(app.evObj.AIChanList);
                app.setLocalStatus('Plotted all analog channels in a separate figure.');

            catch ME
                uialert(app.UIFigure, getReport(ME, 'basic', 'hyperlinks', 'off'), ...
                    'Plot all failed', 'Icon', 'error');
            end

        end

        % Button pushed function: SaveeventsmatButton
        function SaveeventsmatButtonPushed(app, event)
            %SAVEEVENTSMATBUTTONPUSHED Save detected event metadata to events.mat.

            try
                if isempty(app.evObj) || isempty(app.evObj.timestamps)
                    uialert(app.UIFigure, ...
                        'No detected events are available to save.', ...
                        'No events', ...
                        'Icon', 'warning');
                    return
                end

                if isempty(app.SaveFolder) || ~isfolder(app.SaveFolder)
                    startFolder = app.RawFolder;
                    if isempty(startFolder) || ~isfolder(startFolder)
                        startFolder = pwd;
                    end

                    selectedFolder = uigetdir(startFolder, ...
                        'Select folder to save events.mat');

                    if isequal(selectedFolder, 0)
                        app.setLocalStatus('Save events.mat cancelled.');
                        return
                    end

                    app.SaveFolder = selectedFolder;
                end

                eventsPath = fullfile(app.SaveFolder, 'events.mat');

                if isfile(eventsPath)
                    answer = uiconfirm(app.UIFigure, ...
                        sprintf('Overwrite existing events.mat in:\n\n%s', app.SaveFolder), ...
                        'Overwrite events.mat?', ...
                        'Options', {'Overwrite', 'Cancel'}, ...
                        'DefaultOption', 'Cancel', ...
                        'CancelOption', 'Cancel', ...
                        'Icon', 'warning');

                    if strcmpi(answer, 'Cancel')
                        app.setLocalStatus('Save events.mat cancelled.');
                        return
                    end
                end

                app.evObj.saveEvents(app.SaveFolder);
                app.evObj.loadEvents(app.SaveFolder);

                app.HasUnsavedChanges = false;
                app.SaveeventsmatButton.Enable = 'off';

                % Do not call storeOriginalEventNameList() here. Keeping the previous
                % original-name baseline allows Restore to undo manual renames even
                % after saving events.mat in the same session.
                app.refreshDetectionSummary();
                app.refreshConditionSummaryTable();
                app.refreshConditionNameEditButtonState();
                app.refreshParentDataViewerEvents();

                uialert(app.UIFigure, ...
                    sprintf('events.mat saved in:\n\n%s', app.SaveFolder), ...
                    'Events saved', ...
                    'Icon', 'success');

                app.setLocalStatus('events.mat saved.');

            catch ME
                app.setLocalStatus(sprintf('Save events.mat failed: %s', ME.message));
                uialert(app.UIFigure, getReport(ME, 'basic', 'hyperlinks', 'off'), ...
                    'Save events.mat failed', 'Icon', 'error');
            end

        end

        % Callback function: CloseButton, UIFigure
        function UIFigureCloseRequest(app, event)
            app.closeEventsManagerApp();

        end

        % Button pushed function: SetBaselinePeriodButton
        function SetBaselinePeriodButtonPushed(app, event)

            %SETBASELINEPERIODBUTTONPUSHED Apply manual baseline period.
            %
            %   The baseline period is saved in events.mat and later used by
            %   DataViewer to split data into event-aligned trials. It is not a
            %   trigger-detection parameter.

            try
                if isempty(app.evObj) || isempty(app.evObj.timestamps)
                    uialert(app.UIFigure, ...
                        'Detect events before setting the baseline period.', ...
                        'No detected events', ...
                        'Icon', 'warning');
                    return
                end

                baselinePeriod = double(app.BaselinesEditField.Value);

                if ~isscalar(baselinePeriod) || ~isfinite(baselinePeriod) || baselinePeriod <= 0
                    uialert(app.UIFigure, ...
                        'Baseline period must be a finite positive value in seconds.', ...
                        'Invalid baseline period', ...
                        'Icon', 'warning');
                    return
                end

                app.evObj.setBaselinePeriod(baselinePeriod);

                app.BaselinesEditField.Value = double(app.evObj.baselinePeriod);
                app.HasUnsavedChanges = true;
                app.SaveeventsmatButton.Enable = 'on';
                app.SetBaselinePeriodButton.Enable = 'off';

                app.refreshDetectionSummary();
                app.setLocalStatus('Baseline period set. Save events.mat to keep this change.');

            catch ME
                app.setLocalStatus(sprintf('Set baseline period failed: %s', ME.message));
                uialert(app.UIFigure, ...
                    getReport(ME, 'basic', 'hyperlinks', 'off'), ...
                    'Set baseline period failed', ...
                    'Icon', 'error');
            end
        end

        % Value changed function: BaselinesEditField
        function BaselinesEditFieldValueChanged(app, event)
            %BASELINESEDITFIELDVALUECHANGED Enable saving pending baseline changes.

            if isempty(app.evObj) || isempty(app.evObj.timestamps)
                app.SetBaselinePeriodButton.Enable = 'off';
                return
            end

            if isempty(app.evObj.baselinePeriod)
                currentBaseline = NaN;
            else
                currentBaseline = double(app.evObj.baselinePeriod);
            end

            newBaseline = double(app.BaselinesEditField.Value);

            if isfinite(newBaseline) && newBaseline > 0 && ...
                    ~isequaln(newBaseline, currentBaseline)
                app.SetBaselinePeriodButton.Enable = 'on';
                app.setLocalStatus('Baseline period changed. Click Set to apply.');
            else
                app.SetBaselinePeriodButton.Enable = 'off';
            end

        end

        % Button pushed function: EditConditionNamesButton
        function EditConditionNamesButtonPushed(app, event)
            %EDITCONDITIONNAMESBUTTONPUSHED Toggle/apply condition-name editing.
            %
            %   First click enters edit mode. The ConditionName column becomes
            %   editable and highlighted, and the button text changes to Apply.
            %   Second click validates and applies the edited names to EventsManager.

            if isempty(app.evObj) || isempty(app.evObj.eventNameList)
                uialert(app.UIFigure, ...
                    'Detect events before editing condition names.', ...
                    'No conditions', ...
                    'Icon', 'warning');
                return
            end

            if ~app.IsEditingConditionNames
                app.enterConditionNameEditMode();
                return
            end

            app.applyEditedConditionNames();


        end

        % Button pushed function: RestoreConditionNamesButton
        function RestoreConditionNamesButtonPushed(app, event)
            %RESTORECONDITIONNAMESBUTTONPUSHED Cancel edits or restore original names.

            try
                if app.IsEditingConditionNames
                    % Cancel pending table edits without changing EventsManager.
                    app.IsEditingConditionNames = false;
                    app.refreshConditionSummaryTable();
                    app.exitConditionNameEditMode();
                    app.setLocalStatus('Condition name editing cancelled.');
                    return
                end

                if isempty(app.evObj) || isempty(app.evObj.eventNameList)
                    return
                end

                if isempty(app.OriginalEventNameList)
                    uialert(app.UIFigure, ...
                        'No original condition names are available to restore.', ...
                        'Restore unavailable', ...
                        'Icon', 'info');
                    return
                end

                app.evObj.setEventNameList(app.OriginalEventNameList);

                app.HasUnsavedChanges = true;
                app.SaveeventsmatButton.Enable = 'on';

                app.refreshDetectionSummary();
                app.refreshConditionSummaryTable();
                app.refreshConditionNameEditButtonState();
                app.refreshPreviewAxes();

                app.setLocalStatus('Original condition names restored. Save events.mat to keep this change.');

            catch ME
                uialert(app.UIFigure, ...
                    getReport(ME, 'basic', 'hyperlinks', 'off'), ...
                    'Restore condition names failed', ...
                    'Icon', 'error');
            end

        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1191 862];
            app.UIFigure.Name = 'Events Manager app';
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, @UIFigureCloseRequest, true);

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {300, '1x', '1x'};
            app.GridLayout.RowHeight = {60, 210, 100, '1x', '1x', 40};
            app.GridLayout.ColumnSpacing = 5;

            % Create AnalogInputsPanel
            app.AnalogInputsPanel = uipanel(app.GridLayout);
            app.AnalogInputsPanel.BorderType = 'none';
            app.AnalogInputsPanel.Title = 'Analog Inputs ';
            app.AnalogInputsPanel.Layout.Row = 1;
            app.AnalogInputsPanel.Layout.Column = [1 3];
            app.AnalogInputsPanel.FontWeight = 'bold';

            % Create GridAnalogInputs
            app.GridAnalogInputs = uigridlayout(app.AnalogInputsPanel);
            app.GridAnalogInputs.ColumnWidth = {120, '1x'};
            app.GridAnalogInputs.RowHeight = {'1x'};
            app.GridAnalogInputs.ColumnSpacing = 5;
            app.GridAnalogInputs.RowSpacing = 0;
            app.GridAnalogInputs.Padding = [0 0 0 5];

            % Create LoadAnalogInputsButton
            app.LoadAnalogInputsButton = uibutton(app.GridAnalogInputs, 'push');
            app.LoadAnalogInputsButton.ButtonPushedFcn = createCallbackFcn(app, @LoadAnalogInputsButtonPushed, true);
            app.LoadAnalogInputsButton.Layout.Row = 1;
            app.LoadAnalogInputsButton.Layout.Column = 1;
            app.LoadAnalogInputsButton.Text = 'Load Analog Inputs';

            % Create SourceStatusLabel
            app.SourceStatusLabel = uilabel(app.GridAnalogInputs);
            app.SourceStatusLabel.Layout.Row = 1;
            app.SourceStatusLabel.Layout.Column = 2;
            app.SourceStatusLabel.Text = 'SourceStatus';

            % Create DetectionParametersPanel
            app.DetectionParametersPanel = uipanel(app.GridLayout);
            app.DetectionParametersPanel.Title = 'Detection Parameters';
            app.DetectionParametersPanel.Layout.Row = 2;
            app.DetectionParametersPanel.Layout.Column = 1;
            app.DetectionParametersPanel.FontWeight = 'bold';

            % Create GridDetectionParameters
            app.GridDetectionParameters = uigridlayout(app.DetectionParametersPanel);
            app.GridDetectionParameters.ColumnWidth = {80, '1x', 40};
            app.GridDetectionParameters.RowHeight = {25, 25, 25, 25, 25};

            % Create LowpassHzLabel
            app.LowpassHzLabel = uilabel(app.GridDetectionParameters);
            app.LowpassHzLabel.Layout.Row = 5;
            app.LowpassHzLabel.Layout.Column = 1;
            app.LowpassHzLabel.Text = 'Low-pass Hz: ';

            % Create LowpassHzEditField
            app.LowpassHzEditField = uieditfield(app.GridDetectionParameters, 'numeric');
            app.LowpassHzEditField.Layout.Row = 5;
            app.LowpassHzEditField.Layout.Column = [2 3];

            % Create MininterstimLabel
            app.MininterstimLabel = uilabel(app.GridDetectionParameters);
            app.MininterstimLabel.Layout.Row = 4;
            app.MininterstimLabel.Layout.Column = 1;
            app.MininterstimLabel.Text = 'Min interstim:';

            % Create MininterstimSpinner
            app.MininterstimSpinner = uispinner(app.GridDetectionParameters);
            app.MininterstimSpinner.Layout.Row = 4;
            app.MininterstimSpinner.Layout.Column = [2 3];

            % Create PolarityLabel
            app.PolarityLabel = uilabel(app.GridDetectionParameters);
            app.PolarityLabel.Layout.Row = 3;
            app.PolarityLabel.Layout.Column = 1;
            app.PolarityLabel.Text = 'Polarity:';

            % Create PolarityDropDown
            app.PolarityDropDown = uidropdown(app.GridDetectionParameters);
            app.PolarityDropDown.Items = {'Positive', 'Negative'};
            app.PolarityDropDown.Layout.Row = 3;
            app.PolarityDropDown.Layout.Column = [2 3];
            app.PolarityDropDown.Value = 'Positive';

            % Create TriggertypeDropDownLabel
            app.TriggertypeDropDownLabel = uilabel(app.GridDetectionParameters);
            app.TriggertypeDropDownLabel.Layout.Row = 2;
            app.TriggertypeDropDownLabel.Layout.Column = 1;
            app.TriggertypeDropDownLabel.Text = 'Trigger type:';

            % Create TriggertypeDropDown
            app.TriggertypeDropDown = uidropdown(app.GridDetectionParameters);
            app.TriggertypeDropDown.Items = {'EdgeSet', 'EdgeToggle'};
            app.TriggertypeDropDown.Layout.Row = 2;
            app.TriggertypeDropDown.Layout.Column = [2 3];
            app.TriggertypeDropDown.Value = 'EdgeSet';

            % Create ThresholdEditFieldLabel
            app.ThresholdEditFieldLabel = uilabel(app.GridDetectionParameters);
            app.ThresholdEditFieldLabel.Layout.Row = 1;
            app.ThresholdEditFieldLabel.Layout.Column = 1;
            app.ThresholdEditFieldLabel.Text = 'Threshold:';

            % Create ThresholdEditField
            app.ThresholdEditField = uieditfield(app.GridDetectionParameters, 'text');
            app.ThresholdEditField.HorizontalAlignment = 'right';
            app.ThresholdEditField.Layout.Row = 1;
            app.ThresholdEditField.Layout.Column = [2 3];
            app.ThresholdEditField.Value = 'auto';

            % Create TriggerChannelPanel
            app.TriggerChannelPanel = uipanel(app.GridLayout);
            app.TriggerChannelPanel.Title = 'Trigger channels';
            app.TriggerChannelPanel.Layout.Row = 4;
            app.TriggerChannelPanel.Layout.Column = 1;
            app.TriggerChannelPanel.FontWeight = 'bold';

            % Create TriggerChannelGrid
            app.TriggerChannelGrid = uigridlayout(app.TriggerChannelPanel);
            app.TriggerChannelGrid.ColumnWidth = {'1x'};
            app.TriggerChannelGrid.RowHeight = {'1x'};
            app.TriggerChannelGrid.Scrollable = 'on';

            % Create EmptyLabel
            app.EmptyLabel = uilabel(app.TriggerChannelGrid);
            app.EmptyLabel.HorizontalAlignment = 'center';
            app.EmptyLabel.FontAngle = 'italic';
            app.EmptyLabel.Layout.Row = 1;
            app.EmptyLabel.Layout.Column = 1;
            app.EmptyLabel.Text = 'Empty';

            % Create PreviewPanel
            app.PreviewPanel = uipanel(app.GridLayout);
            app.PreviewPanel.Title = 'Preview';
            app.PreviewPanel.Layout.Row = [2 4];
            app.PreviewPanel.Layout.Column = [2 3];
            app.PreviewPanel.FontWeight = 'bold';

            % Create GridPreviewPanel
            app.GridPreviewPanel = uigridlayout(app.PreviewPanel);
            app.GridPreviewPanel.ColumnWidth = {80, '1x'};
            app.GridPreviewPanel.RowHeight = {30, '1x'};

            % Create UIAxes
            app.UIAxes = uiaxes(app.GridPreviewPanel);
            title(app.UIAxes, 'Title')
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.Layout.Row = 2;
            app.UIAxes.Layout.Column = [1 2];

            % Create PlotAllButton
            app.PlotAllButton = uibutton(app.GridPreviewPanel, 'push');
            app.PlotAllButton.ButtonPushedFcn = createCallbackFcn(app, @PlotAllButtonPushed, true);
            app.PlotAllButton.Layout.Row = 1;
            app.PlotAllButton.Layout.Column = 1;
            app.PlotAllButton.Text = 'Plot All';

            % Create ConditionFilePanel
            app.ConditionFilePanel = uipanel(app.GridLayout);
            app.ConditionFilePanel.Title = 'Condition File';
            app.ConditionFilePanel.Layout.Row = 5;
            app.ConditionFilePanel.Layout.Column = 1;
            app.ConditionFilePanel.FontWeight = 'bold';

            % Create GridConditionFile
            app.GridConditionFile = uigridlayout(app.ConditionFilePanel);
            app.GridConditionFile.ColumnWidth = {100, 80, '1x'};
            app.GridConditionFile.RowHeight = {30, 20, '1x', 30};
            app.GridConditionFile.RowSpacing = 5;

            % Create ConditionsourceDropDownLabel
            app.ConditionsourceDropDownLabel = uilabel(app.GridConditionFile);
            app.ConditionsourceDropDownLabel.Layout.Row = 1;
            app.ConditionsourceDropDownLabel.Layout.Column = 1;
            app.ConditionsourceDropDownLabel.Text = 'Condition source:';

            % Create ConditionsourceDropDown
            app.ConditionsourceDropDown = uidropdown(app.GridConditionFile);
            app.ConditionsourceDropDown.Items = {'none', 'csv', 'vpixx'};
            app.ConditionsourceDropDown.ValueChangedFcn = createCallbackFcn(app, @ConditionsourceDropDownValueChanged, true);
            app.ConditionsourceDropDown.Layout.Row = 1;
            app.ConditionsourceDropDown.Layout.Column = 2;
            app.ConditionsourceDropDown.Value = 'none';

            % Create CSVColumnsLabel
            app.CSVColumnsLabel = uilabel(app.GridConditionFile);
            app.CSVColumnsLabel.FontWeight = 'bold';
            app.CSVColumnsLabel.Layout.Row = 2;
            app.CSVColumnsLabel.Layout.Column = 1;
            app.CSVColumnsLabel.Text = 'CSV Columns';

            % Create CSVColumnsTable
            app.CSVColumnsTable = uitable(app.GridConditionFile);
            app.CSVColumnsTable.ColumnName = {'Use'; 'Column'};
            app.CSVColumnsTable.RowName = {};
            app.CSVColumnsTable.Layout.Row = 3;
            app.CSVColumnsTable.Layout.Column = [1 3];

            % Create UpdateConditionsButton
            app.UpdateConditionsButton = uibutton(app.GridConditionFile, 'push');
            app.UpdateConditionsButton.ButtonPushedFcn = createCallbackFcn(app, @UpdateConditionsButtonPushed, true);
            app.UpdateConditionsButton.Layout.Row = 4;
            app.UpdateConditionsButton.Layout.Column = [1 2];
            app.UpdateConditionsButton.Text = 'Update Conditions';

            % Create ChoosefileButton
            app.ChoosefileButton = uibutton(app.GridConditionFile, 'push');
            app.ChoosefileButton.ButtonPushedFcn = createCallbackFcn(app, @ChoosefileButtonPushed, true);
            app.ChoosefileButton.Layout.Row = 1;
            app.ChoosefileButton.Layout.Column = 3;
            app.ChoosefileButton.Text = 'Choose file...';

            % Create GridFooter
            app.GridFooter = uigridlayout(app.GridLayout);
            app.GridFooter.ColumnWidth = {'1x', 100, 120, 100};
            app.GridFooter.RowHeight = {'1x'};
            app.GridFooter.RowSpacing = 0;
            app.GridFooter.Padding = [0 0 0 0];
            app.GridFooter.Layout.Row = 6;
            app.GridFooter.Layout.Column = [1 3];

            % Create DetectButton
            app.DetectButton = uibutton(app.GridFooter, 'push');
            app.DetectButton.ButtonPushedFcn = createCallbackFcn(app, @DetectButtonPushed, true);
            app.DetectButton.Layout.Row = 1;
            app.DetectButton.Layout.Column = 2;
            app.DetectButton.Text = 'Detect';

            % Create SaveeventsmatButton
            app.SaveeventsmatButton = uibutton(app.GridFooter, 'push');
            app.SaveeventsmatButton.ButtonPushedFcn = createCallbackFcn(app, @SaveeventsmatButtonPushed, true);
            app.SaveeventsmatButton.Layout.Row = 1;
            app.SaveeventsmatButton.Layout.Column = 3;
            app.SaveeventsmatButton.Text = 'Save events.mat';

            % Create CloseButton
            app.CloseButton = uibutton(app.GridFooter, 'push');
            app.CloseButton.ButtonPushedFcn = createCallbackFcn(app, @UIFigureCloseRequest, true);
            app.CloseButton.Layout.Row = 1;
            app.CloseButton.Layout.Column = 4;
            app.CloseButton.Text = 'Close';

            % Create StatusLabel
            app.StatusLabel = uilabel(app.GridFooter);
            app.StatusLabel.FontAngle = 'italic';
            app.StatusLabel.Layout.Row = 1;
            app.StatusLabel.Layout.Column = 1;
            app.StatusLabel.Text = 'Status';

            % Create TrialSplittingPanel
            app.TrialSplittingPanel = uipanel(app.GridLayout);
            app.TrialSplittingPanel.Title = 'Trial Splitting';
            app.TrialSplittingPanel.Layout.Row = 3;
            app.TrialSplittingPanel.Layout.Column = 1;
            app.TrialSplittingPanel.FontWeight = 'bold';

            % Create GridTrialSplitting
            app.GridTrialSplitting = uigridlayout(app.TrialSplittingPanel);
            app.GridTrialSplitting.ColumnWidth = {60, '1x', 60};
            app.GridTrialSplitting.RowHeight = {'1x', 30};
            app.GridTrialSplitting.Padding = [5 5 5 5];

            % Create BaselinesEditFieldLabel
            app.BaselinesEditFieldLabel = uilabel(app.GridTrialSplitting);
            app.BaselinesEditFieldLabel.Layout.Row = 1;
            app.BaselinesEditFieldLabel.Layout.Column = 1;
            app.BaselinesEditFieldLabel.Text = 'Baseline s:';

            % Create BaselinesEditField
            app.BaselinesEditField = uieditfield(app.GridTrialSplitting, 'numeric');
            app.BaselinesEditField.ValueChangedFcn = createCallbackFcn(app, @BaselinesEditFieldValueChanged, true);
            app.BaselinesEditField.Layout.Row = 1;
            app.BaselinesEditField.Layout.Column = 2;

            % Create AutoBaselineButton
            app.AutoBaselineButton = uibutton(app.GridTrialSplitting, 'push');
            app.AutoBaselineButton.ButtonPushedFcn = createCallbackFcn(app, @AutoBaselineButtonPushed, true);
            app.AutoBaselineButton.Layout.Row = 1;
            app.AutoBaselineButton.Layout.Column = 3;
            app.AutoBaselineButton.Text = 'Auto';

            % Create SetBaselinePeriodButton
            app.SetBaselinePeriodButton = uibutton(app.GridTrialSplitting, 'push');
            app.SetBaselinePeriodButton.ButtonPushedFcn = createCallbackFcn(app, @SetBaselinePeriodButtonPushed, true);
            app.SetBaselinePeriodButton.Layout.Row = 2;
            app.SetBaselinePeriodButton.Layout.Column = [1 3];
            app.SetBaselinePeriodButton.Text = 'Set';

            % Create DetectionSummaryPanel
            app.DetectionSummaryPanel = uipanel(app.GridLayout);
            app.DetectionSummaryPanel.Title = 'Detection Summary';
            app.DetectionSummaryPanel.Layout.Row = 5;
            app.DetectionSummaryPanel.Layout.Column = [2 3];
            app.DetectionSummaryPanel.FontWeight = 'bold';

            % Create GridDetectionSummary
            app.GridDetectionSummary = uigridlayout(app.DetectionSummaryPanel);
            app.GridDetectionSummary.ColumnWidth = {260, 120, 70, 70, '1x'};
            app.GridDetectionSummary.RowHeight = {30, '1x'};
            app.GridDetectionSummary.ColumnSpacing = 5;
            app.GridDetectionSummary.RowSpacing = 5;
            app.GridDetectionSummary.Padding = [10 3 10 10];

            % Create DetectionSummaryTextArea
            app.DetectionSummaryTextArea = uitextarea(app.GridDetectionSummary);
            app.DetectionSummaryTextArea.Editable = 'off';
            app.DetectionSummaryTextArea.BackgroundColor = [0.9412 0.9412 0.9412];
            app.DetectionSummaryTextArea.Layout.Row = [1 2];
            app.DetectionSummaryTextArea.Layout.Column = 1;

            % Create ConditionSummaryTable
            app.ConditionSummaryTable = uitable(app.GridDetectionSummary);
            app.ConditionSummaryTable.ColumnName = {'Column 1'; 'Column 2'; 'Column 3'; 'Column 4'};
            app.ConditionSummaryTable.RowName = {};
            app.ConditionSummaryTable.Layout.Row = 2;
            app.ConditionSummaryTable.Layout.Column = [2 5];

            % Create EditConditionNamesButton
            app.EditConditionNamesButton = uibutton(app.GridDetectionSummary, 'push');
            app.EditConditionNamesButton.ButtonPushedFcn = createCallbackFcn(app, @EditConditionNamesButtonPushed, true);
            app.EditConditionNamesButton.Layout.Row = 1;
            app.EditConditionNamesButton.Layout.Column = 3;
            app.EditConditionNamesButton.Text = 'Edit';

            % Create RestoreConditionNamesButton
            app.RestoreConditionNamesButton = uibutton(app.GridDetectionSummary, 'push');
            app.RestoreConditionNamesButton.ButtonPushedFcn = createCallbackFcn(app, @RestoreConditionNamesButtonPushed, true);
            app.RestoreConditionNamesButton.Layout.Row = 1;
            app.RestoreConditionNamesButton.Layout.Column = 4;
            app.RestoreConditionNamesButton.Text = 'Restore';

            % Create EditConditionNamesLabel
            app.EditConditionNamesLabel = uilabel(app.GridDetectionSummary);
            app.EditConditionNamesLabel.Layout.Row = 1;
            app.EditConditionNamesLabel.Layout.Column = 2;
            app.EditConditionNamesLabel.Text = 'Edit condition names:';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = DataViewer_EventsManager_exported(varargin)

            runningApp = getRunningApp(app);

            % Check for running singleton app
            if isempty(runningApp)

                % Create UIFigure and components
                createComponents(app)

                % Register the app with App Designer
                registerApp(app, app.UIFigure)

                % Execute the startup function
                runStartupFcn(app, @(app)startupFcn(app, varargin{:}))
            else

                % Focus the running singleton app
                figure(runningApp.UIFigure)

                app = runningApp;
            end

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