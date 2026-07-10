classdef PipelineManagerTool_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                      matlab.ui.Figure
        FileMenu                      matlab.ui.container.Menu
        LoadPipelineMenu              matlab.ui.container.Menu
        SavePipelineMenu              matlab.ui.container.Menu
        GenerateScriptMenu            matlab.ui.container.Menu
        DataMenu                      matlab.ui.container.Menu
        SelectDataFoldersMenu         matlab.ui.container.Menu
        ExecutionMenu                 matlab.ui.container.Menu
        RunPipelineMenu               matlab.ui.container.Menu
        ReuseExistingFilesMenu        matlab.ui.container.Menu
        ExecutionSettingsMenu         matlab.ui.container.Menu
        AdvancedRAMSettingsMenu       matlab.ui.container.Menu
        ToolsMenu                     matlab.ui.container.Menu
        ShowExecutionPlanMenu         matlab.ui.container.Menu
        RemoveInvalidStepsMenu        matlab.ui.container.Menu
        ShowLatestRunSummaryMenu      matlab.ui.container.Menu
        ViewFolderPipelineLogMenu     matlab.ui.container.Menu
        ExportErrorLogMenu            matlab.ui.container.Menu
        GridLayout                    matlab.ui.container.GridLayout
        ActivityLogPanel              matlab.ui.container.Panel
        GridActivityLog               matlab.ui.container.GridLayout
        StatusLabel                   matlab.ui.control.Label
        ActivityLogTextArea           matlab.ui.control.TextArea
        ExecutionOrderPanel           matlab.ui.container.Panel
        GridExecutionOrder            matlab.ui.container.GridLayout
        ExecGraphHTML                 matlab.ui.control.HTML
        SelectedstepdetailsPanel      matlab.ui.container.Panel
        GridSelectedStepDetails       matlab.ui.container.GridLayout
        FunctionNameField             matlab.ui.control.EditField
        ParametersTextArea            matlab.ui.control.TextArea
        ParametersLabel               matlab.ui.control.Label
        StepSummaryTextArea           matlab.ui.control.TextArea
        StepSummaryLabel              matlab.ui.control.Label
        StepTagField                  matlab.ui.control.EditField
        StepTagLabel                  matlab.ui.control.Label
        StepNameLabel                 matlab.ui.control.Label
        PipelineGraphPanel            matlab.ui.container.Panel
        GridPipelineGraph             matlab.ui.container.GridLayout
        ExecutionProgressPanel        matlab.ui.container.Panel
        GridExecutionProgress         matlab.ui.container.GridLayout
        StepProgressOuterPanel        matlab.ui.container.Panel
        FolderProgressOuterPanel      matlab.ui.container.Panel
        ClearPipelineButton           matlab.ui.control.Button
        GraphHTML                     matlab.ui.control.HTML
        AvailableFunctionsPanel       matlab.ui.container.Panel
        GridAvailableFunctions        matlab.ui.container.GridLayout
        AddasNewBranchButton          matlab.ui.control.Button
        AddButton                     matlab.ui.control.Button
        FunctionTree                  matlab.ui.container.Tree
        SearchFunctionEditField       matlab.ui.control.EditField
        GridLayoutTop                 matlab.ui.container.GridLayout
        ReuseExistingFilesCheckBox    matlab.ui.control.CheckBox
        Sep2                          matlab.ui.container.Panel
        Sep1                          matlab.ui.container.Panel
        SafetyLabel                   matlab.ui.control.Label
        RAMLabel                      matlab.ui.control.Label
        AutoSaveFinalOutputsLabel     matlab.ui.control.Label
        RAMSafetyPolicyDropDown       matlab.ui.control.DropDown
        RAMModeDropDown               matlab.ui.control.DropDown
        AutoSaveFinalOutputsDropDown  matlab.ui.control.DropDown
        RunButton                     matlab.ui.control.Button
        SelectDataFoldersButton       matlab.ui.control.Button
        SaveButton                    matlab.ui.control.Button
        LoadButton                    matlab.ui.control.Button
    end


    properties (Access = private)
        pm = []
        graphHtmlFile char = ''
        execHtmlFile char = ''
        selectedNodeID double = NaN
        selectedStepTag char = ''
        selectedFunctionName char = ''
        bGraphHtmlReady logical = false
        bExecHtmlReady logical = false

        bReconnectMode logical = false
        reconnectTargetNodeID double = NaN
        reconnectTargetStepTag char = ''
        reconnectTargetInputName char = ''
        reconnectCandidateNodeIDs double = []
        FolderPairTable table = table()
        ExecutionPlanFigure = []
        PipelineSummaryFigure = []
        LastExecutionResult = []

        bPipelineRunning logical = false
        bCancelRequested logical = false
        RunButtonDefaultBackgroundColor double = []
        RunButtonDefaultFontColor double = []
        CurrentFolderIndex double = NaN
        NumExecutionFolders double = NaN
        CurrentStepIndex double = NaN
        NumExecutionSteps double = NaN
        CurrentExecutionFolder string = ""
        CurrentExecutionStep string = ""
        FolderProgressHTML = []
        StepProgressHTML = []
        RuntimeCurrentStepID double = NaN
        RuntimeCompletedStepIDs double = []
        RuntimeSkippedStepIDs double = []
        RuntimeFailedStepIDs double = []
        RuntimeLastEventType string = ""

        bDataViewerMode logical = false
        DataViewerSaveFolder char = ''
        DataViewerRawFolder char = 'Missing'
        DataViewerCurrentFilePath char = ''
        DataViewerFileSourceNodeID double = NaN
        DataViewerExecutionFinishedFcn = []
        DataViewerToolClosedFcn = []
        bDataViewerExecutionResultNotified logical = false


    end

    methods (Access = private)

        function parseStartupArguments(app, varargin)
            %PARSESTARTUPARGUMENTS Parse startup name-value inputs.

            if isempty(varargin)
                return
            end

            p = inputParser;
            addParameter(p, 'Mode', '', @(x)ischar(x)||isstring(x));
            addParameter(p, 'SaveFolder', '', @(x)ischar(x)||isstring(x));
            addParameter(p, 'RawFolder', 'Missing', @(x)ischar(x)||isstring(x));
            addParameter(p, 'CurrentFilePath', '', @(x)ischar(x)||isstring(x));
            addParameter(p, 'PipelineManager', [], @(x)isempty(x)||isa(x,'PipelineManager'));
            addParameter(p, 'DataViewerFileSourceNodeID', NaN, @(x)isnumeric(x)&&isscalar(x));
            addParameter(p, 'LastExecutionResult', [], @(x)true);
            addParameter(p, 'ExecutionFinishedFcn', [], @(x)isempty(x)||isa(x,'function_handle'));
            addParameter(p, 'ToolClosedFcn', [], @(x)isempty(x)||isa(x,'function_handle'));
            parse(p, varargin{:});

            modeValue = lower(char(string(p.Results.Mode)));

            if strcmp(modeValue, 'dataviewer')
                app.bDataViewerMode = true;
                app.DataViewerSaveFolder = char(string(p.Results.SaveFolder));
                app.DataViewerRawFolder = char(string(p.Results.RawFolder));
                app.DataViewerCurrentFilePath = char(string(p.Results.CurrentFilePath));
                app.DataViewerFileSourceNodeID = double(p.Results.DataViewerFileSourceNodeID);
                app.LastExecutionResult = p.Results.LastExecutionResult;
                app.DataViewerExecutionFinishedFcn = p.Results.ExecutionFinishedFcn;
                app.DataViewerToolClosedFcn = p.Results.ToolClosedFcn;
                app.bDataViewerExecutionResultNotified = false;

                if ~isempty(p.Results.PipelineManager)
                    app.pm = p.Results.PipelineManager;
                end

                if isempty(app.DataViewerRawFolder) || ...
                        (~strcmpi(app.DataViewerRawFolder, 'Missing') && ~isfolder(app.DataViewerRawFolder))
                    app.DataViewerRawFolder = 'Missing';
                end
            end
        end

        function initializeUiState(app)
            %INITIALIZEUISTATE Reset passive UI content and local state.

            app.UIFigure.Name = 'Pipeline Manager Tool';

            app.pm = [];
            app.LastExecutionResult = [];


            app.selectedNodeID = NaN;
            app.selectedStepTag = '';
            app.selectedFunctionName = '';
            app.bGraphHtmlReady = false;
            app.bExecHtmlReady = false;

            app.bReconnectMode = false;
            app.reconnectTargetNodeID = NaN;
            app.reconnectTargetStepTag = '';
            app.reconnectTargetInputName = '';
            app.reconnectCandidateNodeIDs = [];

            app.resetRuntimeExecutionState(true);

            app.FolderPairTable = app.makeEmptyFolderPairTable();

            app.StepTagField.Value = '';
            app.FunctionNameField.Value = '';
            app.StepSummaryTextArea.Value = {'<no step selected>'};
            app.ParametersTextArea.Value = {'<no step selected>'};
            app.ActivityLogTextArea.Value = {''};

            app.AddButton.Enable = 'off';
            app.AddasNewBranchButton.Enable = 'off';
            app.StepNameLabel.Text = 'Function';

            app.configureExecutionControlDropDowns();
            app.closeExecutionPlanWindow();
            app.closePipelineSummaryWindow();

            try
                app.ClearPipelineButton.Tooltip = ...
                    'Remove all steps from the current pipeline while keeping the selected data folders and execution settings.';
            catch
            end

            app.setPipelineControlsAvailable(false);
            app.updateFolderSelectionButtonSummary();
            app.updateLatestRunSummaryMenuState();

            app.setStatus('Select data folders to start.');

            try
                delete(app.FunctionTree.Children);
            catch
            end
        end

        function ensureHtmlFilesExist(app)
            %ENSUREHTMLFILESEXIST Resolve the HTML files used by the uihtml controls.

            appFolder = fileparts(mfilename('fullpath'));

            graphTemplate = fullfile(appFolder, 'pipeline_manager_graph.html');
            execTemplate  = fullfile(appFolder, 'pipeline_manager_execution.html');

            if ~isfile(graphTemplate)
                error('PipelineManagerTool:ensureHtmlFilesExist:MissingGraphHtml', ...
                    ['The graph HTML file was not found. Expected file:\n%s\n\n' ...
                    'Place pipeline_manager_graph.html in the same folder as this app.'], ...
                    graphTemplate);
            end

            if ~isfile(execTemplate)
                error('PipelineManagerTool:ensureHtmlFilesExist:MissingExecutionHtml', ...
                    ['The execution-order HTML file was not found. Expected file:\n%s\n\n' ...
                    'Place pipeline_manager_execution.html in the same folder as this app.'], ...
                    execTemplate);
            end

            app.graphHtmlFile = graphTemplate;
            app.execHtmlFile = execTemplate;

            app.GraphHTML.HTMLSource = app.graphHtmlFile;
            app.ExecGraphHTML.HTMLSource = app.execHtmlFile;
        end

        function populateFunctionTree(app)
            %POPULATEFUNCTIONTREE Populate the available-function tree.

            app.AddButton.Enable = 'off';
            app.selectedFunctionName = '';

            try
                delete(app.FunctionTree.Children);
            catch
            end

            if isempty(app.pm) || ~isprop(app.pm, 'funcList') || isempty(app.pm.funcList)
                placeholderNode = uitreenode(app.FunctionTree);
                placeholderNode.Text = 'Select data folders first';
                placeholderNode.NodeData = struct('kind', 'placeholder');
                return
            end

            funcs = app.pm.funcList;

            % In DataViewer mode, the input file is already provided by a protected
            % FileSource node. Data Import functions are therefore intentionally
            % hidden from the function tree instead of disabled. This keeps the UI
            % unambiguous and prevents accidental creation of a second import entry
            % point while preserving user-added parallel FileSource branches.
            if app.bDataViewerMode && ~isempty(funcs)
                keepMask = true(numel(funcs), 1);
                for iFunc = 1:numel(funcs)
                    keepMask(iFunc) = ~app.shouldHideFunctionInDataViewerMode(funcs(iFunc));
                end
                funcs = funcs(keepMask);
            end

            filterText = strtrim(char(string(app.SearchFunctionEditField.Value)));
            if ~isempty(filterText)
                keepMask = false(numel(funcs), 1);
                filterLower = lower(filterText);

                for iFunc = 1:numel(funcs)
                    fName = lower(char(string(funcs(iFunc).name)));
                    fFolder = '';
                    if isfield(funcs(iFunc), 'folder') && ~isempty(funcs(iFunc).folder)
                        fFolder = lower(char(string(funcs(iFunc).folder)));
                    end
                    keepMask(iFunc) = contains(fName, filterLower) || contains(fFolder, filterLower);
                end

                funcs = funcs(keepMask);
            end

            if isempty(funcs)
                placeholderNode = uitreenode(app.FunctionTree);
                if app.bDataViewerMode
                    placeholderNode.Text = 'No matching functions available in DataViewer mode';
                else
                    placeholderNode.Text = 'No matching functions';
                end
                placeholderNode.NodeData = struct('kind', 'placeholder');
                return
            end

            folderList = strings(numel(funcs), 1);
            for iFunc = 1:numel(funcs)
                if isfield(funcs(iFunc), 'folder') && ~isempty(funcs(iFunc).folder)
                    folderList(iFunc) = string(funcs(iFunc).folder);
                else
                    folderList(iFunc) = "Uncategorized";
                end
            end

            % Inline common-folder-prefix logic. This avoids keeping three small
            % path-format helpers that are only used for this tree display.
            commonRoot = '';
            if ~isempty(folderList)
                partsCell = cell(numel(folderList), 1);
                minLen = inf;

                for i = 1:numel(folderList)
                    parts = regexp(char(folderList(i)), ['[\\' filesep '/]+'], 'split');
                    parts = parts(~cellfun(@isempty, parts));
                    partsCell{i} = parts;
                    minLen = min(minLen, numel(parts));
                end

                commonParts = {};
                for iPart = 1:minLen
                    token = partsCell{1}{iPart};
                    isShared = true;

                    for iList = 2:numel(partsCell)
                        if ~strcmpi(partsCell{iList}{iPart}, token)
                            isShared = false;
                            break
                        end
                    end

                    if ~isShared
                        break
                    end

                    commonParts{end+1} = token; %#ok<AGROW>
                end

                if ~isempty(commonParts)
                    commonRoot = strjoin(commonParts, filesep);
                end
            end

            nodeMap = containers.Map('KeyType', 'char', 'ValueType', 'any');

            for iFunc = 1:numel(funcs)
                folderStr = char(folderList(iFunc));
                relFolder = folderStr;

                if ~isempty(commonRoot)
                    prefixWithSep = [commonRoot filesep];

                    if startsWith(folderStr, prefixWithSep)
                        relFolder = extractAfter(string(folderStr), strlength(prefixWithSep));
                        relFolder = char(relFolder);
                    elseif strcmpi(folderStr, commonRoot)
                        relFolder = '';
                    end
                end

                if isempty(relFolder)
                    folderParts = {};
                else
                    folderParts = regexp(relFolder, ['[\\' filesep '/]+'], 'split');
                    folderParts = folderParts(~cellfun(@isempty, folderParts));
                end

                if isempty(folderParts)
                    folderParts = {'Uncategorized'};
                end

                parentHandle = app.FunctionTree;
                currentKey = '';

                for iPart = 1:numel(folderParts)
                    thisPart = folderParts{iPart};

                    if isempty(currentKey)
                        currentKey = thisPart;
                    else
                        currentKey = [currentKey '|' thisPart]; %#ok<AGROW>
                    end

                    if ~isKey(nodeMap, currentKey)
                        catNode = uitreenode(parentHandle);
                        catNode.Text = thisPart;
                        catNode.NodeData = struct('kind', 'category', 'pathKey', currentKey);
                        nodeMap(currentKey) = catNode;
                    end

                    parentHandle = nodeMap(currentKey);
                end

                funcNode = uitreenode(parentHandle);
                funcNode.Text = char(string(funcs(iFunc).name));
                funcNode.NodeData = struct( ...
                    'kind', 'function', ...
                    'name', char(string(funcs(iFunc).name)), ...
                    'folder', folderStr);
            end

            rootNodes = app.FunctionTree.Children;
            for iRoot = 1:numel(rootNodes)
                try
                    expand(rootNodes(iRoot));
                catch
                end
            end
        end

        function tf = shouldHideFunctionInDataViewerMode(app, funcInfo)
            %SHOULDHIDEFUNCTIONINDATAVIEWERMODE True for Data Import functions.
            %
            % In DataViewer mode, the currently displayed file is already the input
            % source. Pipeline-level import functions are hidden so the user can
            % build processing/export branches without accidentally creating a new
            % import entry point.

            tf = false;

            if ~app.bDataViewerMode
                return
            end

            funcName = "";
            folderText = "";

            try
                if isfield(funcInfo, 'name') && ~isempty(funcInfo.name)
                    funcName = string(funcInfo.name);
                end
            catch
                funcName = "";
            end

            try
                if isfield(funcInfo, 'folder') && ~isempty(funcInfo.folder)
                    folderText = string(funcInfo.folder);
                end
            catch
                folderText = "";
            end

            normalizedFolder = lower(regexprep(char(folderText), '[^a-zA-Z0-9]', ''));
            if contains(normalizedFolder, 'dataimport')
                tf = true;
                return
            end

            % Defensive fallback for function libraries where the folder/category
            % metadata is missing or flattened.
            normalizedName = lower(regexprep(char(funcName), '[^a-zA-Z0-9]', ''));
            importFunctionNames = [ ...
                "imagesclassification", ...
                "runimagesclassification", ...
                "importfromtif", ...
                "runimportfromtif", ...
                "importfrombin", ...
                "runimportfrombin"];

            tf = any(strcmp(normalizedName, importFunctionNames));
        end

        function refreshAllViews(app)
            %REFRESHALLVIEWS Refresh graph, execution view, and selected-step details.

            app.refreshGraphView();
            app.refreshExecutionView();
            app.refreshSelectedStepPanel();
        end

        function refreshGraphView(app)
            if ~app.bGraphHtmlReady
                return
            end
            payload = app.buildGraphPayload();
            payload.command = 'refresh';
            app.GraphHTML.Data = payload;
        end

        function refreshExecutionView(app)
            if ~app.bExecHtmlReady
                return
            end
            payload = app.buildExecutionPayload();
            app.ExecGraphHTML.Data = payload;
        end

        function payload = buildGraphPayload(app)
            %BUILDGRAPHPAYLOAD Build the payload rendered by the graph HTML view.

            payload = struct();
            payload.nodes = struct('id', {}, 'stepTag', {}, 'funcName', {}, 'kind', {}, ...
                'x', {}, 'y', {}, 'width', {}, 'height', {}, 'topoIndex', {}, ...
                'scheduleIndex', {}, 'runIndex', {}, 'isExecutable', {}, 'executionRole', {}, ...
                'runtimeState', {}, 'runtimeIcon', {}, ...
                'hasSavedOutputs', {}, 'hasWarning', {}, 'hasError', {}, ...
                'hasInvalidNonblocking', {}, 'hoverText', {}, 'isReconnectTarget', {}, ...
                'isReconnectCandidate', {}, 'isReconnectMuted', {}, ...
                'reconnectChoiceLabels', {}, 'reconnectChoiceRefs', {}, ...
                'isProtectedNode', {}, 'protectedReason', {});

            payload.edges = struct('sourceNodeID', {}, 'targetNodeID', {}, ...
                'sourceLabel', {}, 'targetLabel', {});
            payload.selectedNodeID = app.selectedNodeID;
            payload.reconnectMode = app.bReconnectMode;
            payload.runtime = app.buildHtmlRuntimePayload();

            payload.reconnectBanner = '';
            payload.reconnectGhost = struct([]);

            targetInputDef = [];
            targetNodeLocal = [];

            if app.bReconnectMode && ~isempty(app.reconnectTargetStepTag) && ~isempty(app.reconnectTargetInputName)
                payload.reconnectBanner = sprintf('Reconnect input "%s" of "%s": click a highlighted source or Use File....', ...
                    app.reconnectTargetInputName, app.reconnectTargetStepTag);
            end

            if isempty(app.pm) || isempty(app.pm.nodes)
                return
            end

            plan = app.getCurrentExecutionPlan();
            if isempty(plan)
                return
            end

            nodesLocal = app.pm.nodes;
            connectionsLocal = app.pm.connections;
            topoOrder = app.getTopoOrderFromPlan(plan);
            if isempty(topoOrder)
                topoOrder = [nodesLocal.id];
            end

            levelMap = app.computeNodeLevels(nodesLocal, connectionsLocal, topoOrder);
            nodeTable = table();
            outputPlan = table();
            nodeStatusTable = table();
            edgeTable = table();

            if isfield(plan, 'nodeTable') && istable(plan.nodeTable)
                nodeTable = plan.nodeTable;
            end
            if isfield(plan, 'outputPlan') && istable(plan.outputPlan)
                outputPlan = plan.outputPlan;
            end
            if isfield(plan, 'edgeTable') && istable(plan.edgeTable)
                edgeTable = plan.edgeTable;
            end
            if isfield(plan, 'validationReport') && isstruct(plan.validationReport) && ...
                    isfield(plan.validationReport, 'nodeStatus') && istable(plan.validationReport.nodeStatus)
                nodeStatusTable = plan.validationReport.nodeStatus;
            end

            if app.bReconnectMode
                tgtIdx = find([nodesLocal.id] == app.reconnectTargetNodeID, 1, 'first');
                if ~isempty(tgtIdx)
                    targetNodeLocal = nodesLocal(tgtIdx);
                    if isfield(targetNodeLocal, 'info') && isfield(targetNodeLocal.info, 'inputs') && ~isempty(targetNodeLocal.info.inputs)
                        dataInputsLocal = targetNodeLocal.info.inputs(arrayfun(@(x) isfield(x, 'isData') && x.isData, targetNodeLocal.info.inputs));
                        for iIn = 1:numel(dataInputsLocal)
                            if strcmpi(char(string(dataInputsLocal(iIn).name)), app.reconnectTargetInputName)
                                targetInputDef = dataInputsLocal(iIn);
                                break
                            end
                        end
                    end
                end
            end

            rowCounterMap = containers.Map('KeyType', 'char', 'ValueType', 'double');

            for iTopo = 1:numel(topoOrder)
                nodeID = topoOrder(iTopo);
                nodeIdx = find([nodesLocal.id] == nodeID, 1, 'first');
                if isempty(nodeIdx)
                    continue
                end
                nodeLocal = nodesLocal(nodeIdx);

                if isKey(levelMap, nodeID)
                    levelVal = levelMap(nodeID);
                else
                    levelVal = 1;
                end
                levelKey = sprintf('L%d', levelVal);

                if ~isKey(rowCounterMap, levelKey)
                    rowCounterMap(levelKey) = 0;
                end

                rowVal = rowCounterMap(levelKey) + 1;
                rowCounterMap(levelKey) = rowVal;

                x = 40 + (levelVal - 1) * 240;
                y = 40 + (rowVal - 1) * 120;
                widthVal = 190;
                heightVal = 76;

                nodeRow = app.getNodeTableRowByID(nodeTable, nodeLocal.id);
                orderInfo = app.getNodeOrderMetadataFromPlan(plan, nodeLocal.id);

                if isempty(orderInfo.topoIndex)
                    orderInfo.topoIndex = iTopo;
                end

                funcName = app.getNodeFunctionName(nodeLocal);
                if ~isempty(nodeRow) && ismember('functionName', nodeRow.Properties.VariableNames)
                    funcName = char(string(nodeRow.functionName(1)));
                end

                hasSavedOutputs = app.getNodePersistenceFlagFromPlan(nodeLocal.id, outputPlan);
                statusInfo = app.getNodeStatusInfo(nodeLocal, nodeStatusTable);
                statusInfo = app.addRawFolderWarningToStatusInfo(nodeLocal, statusInfo);
                hoverText = app.buildNodeHoverText(nodeLocal, orderInfo, statusInfo, plan);

                isReconnectTarget = app.bReconnectMode && nodeLocal.id == app.reconnectTargetNodeID;
                isReconnectCandidate = app.bReconnectMode && any(app.reconnectCandidateNodeIDs == nodeLocal.id);
                isReconnectMuted = app.bReconnectMode && ~(isReconnectTarget || isReconnectCandidate);

                reconnectChoiceLabels = {};
                reconnectChoiceRefs = {};

                if isReconnectCandidate && ~isempty(targetInputDef) && ~isempty(targetNodeLocal)
                    dstTypes = {};
                    if isfield(targetInputDef, 'type') && ~isempty(targetInputDef.type)
                        if ischar(targetInputDef.type) || (isstring(targetInputDef.type) && isscalar(targetInputDef.type))
                            dstTypes = {char(string(targetInputDef.type))};
                        else
                            dstTypes = cellstr(string(targetInputDef.type(:)));
                        end
                    end

                    [choiceLabels, choiceRefs] = app.buildExplicitSourceChoices(dstTypes, ...
                        'leafOnly', false, ...
                        'excludeNodeID', targetNodeLocal.id, ...
                        'candidateNodeIDs', nodeLocal.id);

                    reconnectChoiceLabels = choiceLabels;
                    reconnectChoiceRefs = choiceRefs;
                end

                runtimeState = app.getRuntimeStepState(nodeLocal.id);
                runtimeIcon = app.getRuntimeStateIcon(runtimeState);

                isProtectedNode = app.isDataViewerProtectedNodeID(nodeLocal.id);
                protectedReason = '';
                if isProtectedNode
                    protectedReason = 'DataViewer input';
                end

                payload.nodes(end+1) = struct( ...
                    'id', nodeLocal.id, ...
                    'stepTag', char(string(nodeLocal.name)), ...
                    'funcName', funcName, ...
                    'kind', char(string(nodeLocal.kind)), ...
                    'x', x, ...
                    'y', y, ...
                    'width', widthVal, ...
                    'height', heightVal, ...
                    'topoIndex', app.emptyToNaN(orderInfo.topoIndex), ...
                    'scheduleIndex', app.emptyToNaN(orderInfo.scheduleIndex), ...
                    'runIndex', app.emptyToNaN(orderInfo.runIndex), ...
                    'isExecutable', logical(orderInfo.isExecutable), ...
                    'executionRole', char(string(orderInfo.executionRole)), ...
                    'runtimeState', runtimeState, ...
                    'runtimeIcon', runtimeIcon, ...
                    'hasSavedOutputs', hasSavedOutputs, ...
                    'hasWarning', statusInfo.hasWarning, ...
                    'hasError', statusInfo.hasError, ...
                    'hasInvalidNonblocking', statusInfo.hasInvalidNonblocking, ...
                    'hoverText', hoverText, ...
                    'isReconnectTarget', isReconnectTarget, ...
                    'isReconnectCandidate', isReconnectCandidate, ...
                    'isReconnectMuted', isReconnectMuted, ...
                    'reconnectChoiceLabels', {reconnectChoiceLabels}, ...
                    'reconnectChoiceRefs', {reconnectChoiceRefs}, ...
                    'isProtectedNode', isProtectedNode, ...
                    'protectedReason', protectedReason);

            end

            if app.bReconnectMode
                ghostY = 40;
                if ~isempty(targetNodeLocal)
                    try
                        targetLevelVal = levelMap(targetNodeLocal.id);
                        rowCounterTmp = 0;
                        for iTopo = 1:numel(topoOrder)
                            curID = topoOrder(iTopo);
                            if curID == targetNodeLocal.id
                                break
                            end
                            if isKey(levelMap, curID) && levelMap(curID) == targetLevelVal
                                rowCounterTmp = rowCounterTmp + 1;
                            end
                        end
                        ghostY = 40 + rowCounterTmp * 120;
                    catch
                        ghostY = 40;
                    end
                end

                payload.reconnectGhost = struct( ...
                    'x', 20, ...
                    'y', ghostY, ...
                    'width', 170, ...
                    'height', 64, ...
                    'stepTag', 'Use File...', ...
                    'funcName', 'File Source', ...
                    'hoverText', sprintf('Select a file from SaveFolder for input "%s".', app.reconnectTargetInputName));
            end

            if istable(edgeTable) && ~isempty(edgeTable) && ...
                    all(ismember({'sourceNodeID','targetNodeID','sourceOutput','selectedFile','targetInput'}, edgeTable.Properties.VariableNames))

                for iEdge = 1:height(edgeTable)
                    srcLabel = char(string(edgeTable.sourceOutput(iEdge)));
                    if strlength(strtrim(string(edgeTable.selectedFile(iEdge)))) > 0
                        srcLabel = char(string(edgeTable.selectedFile(iEdge)));
                    end

                    payload.edges(end+1) = struct( ...
                        'sourceNodeID', edgeTable.sourceNodeID(iEdge), ...
                        'targetNodeID', edgeTable.targetNodeID(iEdge), ...
                        'sourceLabel', srcLabel, ...
                        'targetLabel', char(string(edgeTable.targetInput(iEdge))));
                end
            else
                for iConn = 1:numel(connectionsLocal)
                    connLocal = connectionsLocal(iConn);

                    srcLabel = char(string(connLocal.sourceOutputName));
                    if isfield(connLocal, 'selectedFile') && ~isempty(connLocal.selectedFile)
                        srcLabel = char(string(connLocal.selectedFile));
                    end

                    payload.edges(end+1) = struct( ...
                        'sourceNodeID', connLocal.sourceNodeID, ...
                        'targetNodeID', connLocal.targetNodeID, ...
                        'sourceLabel', srcLabel, ...
                        'targetLabel', char(string(connLocal.targetInputName)));
                end
            end
        end

        function info = getNodeStatusInfo(~, nodeLocal, nodeStatusTable)
            %GETNODESTATUSINFO Return structured node status from diagnosePipeline.

            info = struct( ...
                'status', "valid", ...
                'severity', "none", ...
                'messages', strings(0,1), ...
                'hasWarning', false, ...
                'hasInvalidNonblocking', false, ...
                'hasError', false);

            if ~istable(nodeStatusTable) || isempty(nodeStatusTable)
                return
            end

            if ~ismember('nodeID', nodeStatusTable.Properties.VariableNames) || ...
                    ~ismember('status', nodeStatusTable.Properties.VariableNames)
                return
            end

            idx = find(nodeStatusTable.nodeID == nodeLocal.id, 1, 'first');
            if isempty(idx)
                return
            end

            info.status = string(nodeStatusTable.status(idx));

            if ismember('severity', nodeStatusTable.Properties.VariableNames)
                info.severity = string(nodeStatusTable.severity(idx));
            end

            if ismember('messages', nodeStatusTable.Properties.VariableNames) && iscell(nodeStatusTable.messages)
                msgLocal = string(nodeStatusTable.messages{idx});
                info.messages = msgLocal(:);
            end

            info.hasWarning = info.status == "warning";
            info.hasInvalidNonblocking = info.status == "invalid_nonblocking";
            info.hasError = info.status == "invalid_blocking";
        end

        function hoverText = buildNodeHoverText(app, nodeLocal, orderInfo, statusInfo, plan)
            %BUILDNODEHOVERTEXT Build one compact hover summary string for a node.

            if nargin < 5
                plan = struct();
            end

            lines = strings(0,1);

            stepTag = char(string(nodeLocal.name));
            funcName = app.getNodeFunctionName(nodeLocal);

            lines(end+1,1) = "Step: " + string(stepTag);
            lines(end+1,1) = "Function: " + string(funcName);

            isExecutable = false;
            if isstruct(orderInfo) && isfield(orderInfo, 'isExecutable') && ~isempty(orderInfo.isExecutable)
                isExecutable = logical(orderInfo.isExecutable);
            elseif isfield(nodeLocal, 'kind')
                isExecutable = strcmpi(char(string(nodeLocal.kind)), 'stream');
            end

            if ~isExecutable
                lines(end+1,1) = "Role: Virtual input source";
                lines(end+1,1) = "Run order: <not executable>";
            else
                lines(end+1,1) = "Role: Analysis step";
                lines(end+1,1) = "Run order: " + string(app.valueOrPlaceholder(orderInfo.runIndex));
            end

            lines(end+1,1) = "Schedule index: " + string(app.valueOrPlaceholder(orderInfo.scheduleIndex));
            lines(end+1,1) = "Topological order: " + string(app.valueOrPlaceholder(orderInfo.topoIndex));

            if isstruct(statusInfo)
                lines(end+1,1) = "Status: " + string(statusInfo.status);
                lines(end+1,1) = "Severity: " + string(statusInfo.severity);
            end

            inputLines = string(app.getDataInputLines(nodeLocal));
            outputLines = string(app.getOutputLines(nodeLocal));
            outputPlanLines = string(app.getOutputPlanLinesForNode(plan, nodeLocal.id));
            saveLines = string(app.getSavedOutputLines(nodeLocal));

            if ~isempty(inputLines)
                lines(end+1,1) = "";
                lines(end+1,1) = "Inputs:";
                lines = [lines; inputLines(:)];
            end

            if ~isempty(outputLines)
                lines(end+1,1) = "";
                lines(end+1,1) = "Outputs:";
                lines = [lines; outputLines(:)];
            end

            if ~isempty(outputPlanLines)
                lines(end+1,1) = "";
                lines(end+1,1) = "Output handling plan:";
                lines = [lines; outputPlanLines(:)];
            elseif ~isempty(saveLines)
                lines(end+1,1) = "";
                lines(end+1,1) = "Saved outputs:";
                lines = [lines; saveLines(:)];
            end

            if isstruct(statusInfo) && isfield(statusInfo, 'messages') && ~isempty(statusInfo.messages)
                lines(end+1,1) = "";
                lines(end+1,1) = "Issues:";
                lines = [lines; "- " + string(statusInfo.messages(:))];
            end

            hoverText = char(join(lines, newline));
        end

        function payload = buildExecutionPayload(app)
            %BUILDEXECUTIONPAYLOAD Build the payload rendered by the execution HTML view.

            payload = struct();
            payload.version = 2;
            payload.sources = struct('id', {}, 'sourceIndex', {}, 'fileName', {}, ...
                'semanticType', {}, 'feeds', {}, 'runtimeState', {}, 'runtimeIcon', {});
            payload.steps = struct('id', {}, 'stepTag', {}, 'funcName', {}, ...
                'runIndex', {}, 'scheduleIndex', {}, 'topoIndex', {}, ...
                'inputSummary', {}, 'outputSummary', {}, ...
                'runtimeState', {}, 'runtimeIcon', {});
            payload.selectedNodeID = app.selectedNodeID;
            payload.runtime = app.buildHtmlRuntimePayload();
            payload.note = ['Run order lists executable analysis steps. ' ...
                'File sources are virtual inputs and are not loaded until consumed.'];

            if isempty(app.pm) || isempty(app.pm.nodes)
                return
            end

            plan = app.getCurrentExecutionPlan();
            if isempty(plan)
                return
            end

            if isfield(plan, 'leafOutputPolicy') && strlength(string(plan.leafOutputPolicy)) > 0
                payload.note = sprintf(['Run order lists executable analysis steps. ' ...
                    'File sources are virtual inputs and are not loaded until consumed. Output policy: %s'], ...
                    char(string(plan.leafOutputPolicy)));
            end

            nodeTable = table();
            edgeTable = table();
            if isfield(plan, 'nodeTable') && istable(plan.nodeTable)
                nodeTable = plan.nodeTable;
            end
            if isfield(plan, 'edgeTable') && istable(plan.edgeTable)
                edgeTable = plan.edgeTable;
            end

            % -------------------------------------------------------------
            % Virtual FileSource nodes
            % -------------------------------------------------------------
            sourceRows = app.getVirtualSourceRowsFromNodeTable(nodeTable);
            if ~isempty(sourceRows)
                [~, ord] = sort(sourceRows.scheduleIndex, 'ascend', 'MissingPlacement', 'last');
                sourceRows = sourceRows(ord, :);
            end

            for iSource = 1:height(sourceRows)
                nodeID = sourceRows.nodeID(iSource);
                fileName = char(string(sourceRows.stepName(iSource)));
                feeds = "";
                semanticType = "";

                if istable(edgeTable) && ~isempty(edgeTable) && ismember('sourceNodeID', edgeTable.Properties.VariableNames)
                    eRows = edgeTable(edgeTable.sourceNodeID == nodeID, :);

                    if ~isempty(eRows) && ismember('targetStep', eRows.Properties.VariableNames)
                        feeds = strjoin(unique(string(eRows.targetStep), 'stable'), ', ');
                    end

                    if ~isempty(eRows) && ismember('sourceTypes', eRows.Properties.VariableNames)
                        typeVals = strings(0,1);
                        for iEdge = 1:height(eRows)
                            try
                                typeVals = [typeVals; string(eRows.sourceTypes{iEdge}(:))]; %#ok<AGROW>
                            catch
                            end
                        end
                        typeVals = unique(typeVals(strlength(strtrim(typeVals)) > 0), 'stable');
                        if ~isempty(typeVals)
                            semanticType = strjoin(typeVals, ' / ');
                        end
                    end
                end

                runtimeState = 'source';
                runtimeIcon = '';

                payload.sources(end+1) = struct( ... %#ok<AGROW>
                    'id', nodeID, ...
                    'sourceIndex', iSource, ...
                    'fileName', fileName, ...
                    'semanticType', char(semanticType), ...
                    'feeds', char(feeds), ...
                    'runtimeState', runtimeState, ...
                    'runtimeIcon', runtimeIcon);
            end

            % -------------------------------------------------------------
            % Executable stream nodes only
            % -------------------------------------------------------------
            runOrder = app.getRunOrderFromPlan(plan);
            if isempty(runOrder) && ~isempty(nodeTable)
                runRows = app.getExecutableRowsFromNodeTable(nodeTable);
                if ~isempty(runRows)
                    [~, ord] = sort(runRows.runIndex, 'ascend', 'MissingPlacement', 'last');
                    runOrder = runRows.nodeID(ord).';
                end
            end

            for iStep = 1:numel(runOrder)
                nodeID = runOrder(iStep);
                nodeIdx = find([app.pm.nodes.id] == nodeID, 1, 'first');
                if isempty(nodeIdx)
                    continue
                end

                nodeLocal = app.pm.nodes(nodeIdx);
                nodeRow = app.getNodeTableRowByID(nodeTable, nodeID);
                orderInfo = app.getNodeOrderMetadataFromPlan(plan, nodeID);

                if isempty(orderInfo.runIndex)
                    orderInfo.runIndex = iStep;
                end

                funcName = app.getNodeFunctionName(nodeLocal);
                if ~isempty(nodeRow) && ismember('functionName', nodeRow.Properties.VariableNames)
                    funcName = char(string(nodeRow.functionName(1)));
                end

                inputLines = string(app.getDataInputLines(nodeLocal));
                inputLines = inputLines(~strcmp(inputLines, "<no data inputs>") & ~strcmp(inputLines, "<no inputs>"));
                inputSummary = strjoin(inputLines, ' | ');

                outputPlanLines = string(app.getOutputPlanLinesForNode(plan, nodeID));
                if isempty(outputPlanLines)
                    outputLines = string(app.getOutputLines(nodeLocal));
                    outputLines = outputLines(~strcmp(outputLines, "<no outputs>"));
                    outputSummary = strjoin(outputLines, ' | ');
                else
                    outputSummary = strjoin(outputPlanLines, ' | ');
                end

                runtimeState = app.getRuntimeStepState(nodeLocal.id);
                runtimeIcon = app.getRuntimeStateIcon(runtimeState);

                payload.steps(end+1) = struct( ... %#ok<AGROW>
                    'id', nodeLocal.id, ...
                    'stepTag', char(string(nodeLocal.name)), ...
                    'funcName', funcName, ...
                    'runIndex', app.emptyToNaN(orderInfo.runIndex), ...
                    'scheduleIndex', app.emptyToNaN(orderInfo.scheduleIndex), ...
                    'topoIndex', app.emptyToNaN(orderInfo.topoIndex), ...
                    'inputSummary', char(inputSummary), ...
                    'outputSummary', char(outputSummary), ...
                    'runtimeState', runtimeState, ...
                    'runtimeIcon', runtimeIcon);
            end
        end

        function topoOrder = computeTopoOrderFromGraph(~, nodesLocal, connectionsLocal)
            %COMPUTETOPOORDERFROMGRAPH Compute a topological order from nodes/connections.

            if isempty(nodesLocal)
                topoOrder = [];
                return
            end

            nodeIDs = [nodesLocal.id];

            if isempty(connectionsLocal)
                topoOrder = nodeIDs;
                return
            end

            inDegree = zeros(1, numel(nodeIDs));

            for iConn = 1:numel(connectionsLocal)
                tgtIdx = find(nodeIDs == connectionsLocal(iConn).targetNodeID, 1, 'first');
                if ~isempty(tgtIdx)
                    inDegree(tgtIdx) = inDegree(tgtIdx) + 1;
                end
            end

            queue = nodeIDs(inDegree == 0);
            topoOrder = [];

            while ~isempty(queue)
                currentID = queue(1);
                queue(1) = [];
                topoOrder(end+1) = currentID; %#ok<AGROW>

                outMask = [connectionsLocal.sourceNodeID] == currentID;
                outgoing = connectionsLocal(outMask);

                for iOut = 1:numel(outgoing)
                    tgtIdx = find(nodeIDs == outgoing(iOut).targetNodeID, 1, 'first');
                    if isempty(tgtIdx)
                        continue
                    end

                    inDegree(tgtIdx) = inDegree(tgtIdx) - 1;

                    if inDegree(tgtIdx) == 0
                        queue(end+1) = nodeIDs(tgtIdx); %#ok<AGROW>
                    end
                end
            end

            if numel(topoOrder) ~= numel(nodeIDs)
                topoOrder = nodeIDs;
            end
        end

        function levelMap = computeNodeLevels(~, nodesLocal, connectionsLocal, topoOrder)
            %COMPUTENODELEVELS Compute a simple layered layout index for each node.

            levelMap = containers.Map('KeyType', 'double', 'ValueType', 'double');

            if isempty(nodesLocal)
                return
            end

            % No connections -> every node is a root-level node
            if isempty(connectionsLocal)
                for iTopo = 1:numel(topoOrder)
                    levelMap(topoOrder(iTopo)) = 1;
                end
                return
            end

            for iTopo = 1:numel(topoOrder)
                nodeID = topoOrder(iTopo);

                inMask = [connectionsLocal.targetNodeID] == nodeID;
                incoming = connectionsLocal(inMask);

                if isempty(incoming)
                    levelMap(nodeID) = 1;
                    continue
                end

                parentLevels = zeros(1, numel(incoming));
                for iIn = 1:numel(incoming)
                    srcID = incoming(iIn).sourceNodeID;
                    if isKey(levelMap, srcID)
                        parentLevels(iIn) = levelMap(srcID);
                    else
                        parentLevels(iIn) = 1;
                    end
                end

                levelMap(nodeID) = max(parentLevels) + 1;
            end
        end

        function leafIDs = computeLeafNodeIDs(~, nodesLocal, connectionsLocal)
            %COMPUTELEAFNODEIDS Return node IDs with no outgoing edges.

            if isempty(nodesLocal)
                leafIDs = [];
                return
            end

            nodeIDs = [nodesLocal.id];
            if isempty(connectionsLocal)
                leafIDs = nodeIDs;
                return
            end

            srcIDs = [connectionsLocal.sourceNodeID];
            leafIDs = nodeIDs(~ismember(nodeIDs, srcIDs));
        end

        function refreshSelectedStepPanel(app)
            %REFRESHSELECTEDSTEPPANEL Refresh the selected-step details panel.

            if isempty(app.pm) || isempty(app.pm.nodes) || isnan(app.selectedNodeID)
                app.clearSelectedStepPanel();
                return
            end

            nodeIdx = find([app.pm.nodes.id] == app.selectedNodeID, 1, 'first');
            if isempty(nodeIdx)
                app.clearSelectedStepPanel();
                return
            end

            nodeLocal = app.pm.nodes(nodeIdx);

            app.StepTagField.Value = char(string(nodeLocal.name));
            app.FunctionNameField.Value = app.getNodeFunctionName(nodeLocal);
            app.StepSummaryTextArea.Value = app.buildStepSummaryText(nodeLocal);
            app.ParametersTextArea.Value = app.buildParameterText(nodeLocal);
        end

        function clearSelectedStepPanel(app)
            %CLEARSELECTEDSTEPPANEL Clear the right-hand details panel.

            app.StepTagField.Value = '';
            app.FunctionNameField.Value = '';
            app.StepSummaryTextArea.Value = {'<no step selected>'};
            app.ParametersTextArea.Value = {'<no step selected>'};
        end

        function outLines = buildStepSummaryText(app, nodeLocal)
            %BUILDSTEPSUMMARYTEXT Build the Step Summary text for one selected node.

            outLines = strings(0, 1);

            if isempty(app.pm)
                outLines = "<no step selected>";
                outLines = cellstr(outLines);
                return
            end

            plan = app.getCurrentExecutionPlan();
            orderInfo = app.getNodeOrderMetadataFromPlan(plan, nodeLocal.id);

            if isfield(nodeLocal, 'kind') && strcmpi(char(string(nodeLocal.kind)), 'folder')
                outLines(end+1, 1) = "Role: Virtual input source"; %#ok<AGROW>
                outLines(end+1, 1) = "Run order: <not executable>"; %#ok<AGROW>
            else
                outLines(end+1, 1) = "Role: Analysis step"; %#ok<AGROW>
                outLines(end+1, 1) = "Run order: " + string(app.valueOrPlaceholder(orderInfo.runIndex)); %#ok<AGROW>
            end

            outLines(end+1, 1) = "Schedule index: " + string(app.valueOrPlaceholder(orderInfo.scheduleIndex)); %#ok<AGROW>
            outLines(end+1, 1) = "Topological order: " + string(app.valueOrPlaceholder(orderInfo.topoIndex)); %#ok<AGROW>
            outLines(end+1, 1) = ""; %#ok<AGROW>

            outLines(end+1, 1) = "Inputs:"; %#ok<AGROW>
            dataInputLines = app.getDataInputLines(nodeLocal);
            outLines = [outLines; string(dataInputLines(:))]; %#ok<AGROW>
            outLines(end+1, 1) = ""; %#ok<AGROW>

            outLines(end+1, 1) = "Outputs:"; %#ok<AGROW>
            outputLines = app.getOutputLines(nodeLocal);
            outLines = [outLines; string(outputLines(:))]; %#ok<AGROW>

            outputPlanLines = app.getOutputPlanLinesForNode(plan, nodeLocal.id);
            if ~isempty(outputPlanLines)
                outLines(end+1, 1) = ""; %#ok<AGROW>
                outLines(end+1, 1) = "Output handling plan:"; %#ok<AGROW>
                outLines = [outLines; string(outputPlanLines(:))]; %#ok<AGROW>
            else
                saveLines = app.getSavedOutputLines(nodeLocal);
                if ~isempty(saveLines)
                    outLines(end+1, 1) = ""; %#ok<AGROW>
                    outLines(end+1, 1) = "Saved files:"; %#ok<AGROW>
                    outLines = [outLines; string(saveLines(:))]; %#ok<AGROW>
                end
            end

            outLines = cellstr(outLines);
        end

        function outLines = getDataInputLines(app, nodeLocal)
            %GETDATAINPUTLINES Build the input lines shown in the selected-step summary.

            outLines = strings(0, 1);

            if strcmpi(nodeLocal.kind, 'folder')
                outLines(end+1, 1) = "<file source node>"; %#ok<AGROW>
                outLines = cellstr(outLines);
                return
            end

            if ~isfield(nodeLocal, 'info') || ~isfield(nodeLocal.info, 'inputs') || isempty(nodeLocal.info.inputs)
                outLines(end+1, 1) = "<no inputs>"; %#ok<AGROW>
                outLines = cellstr(outLines);
                return
            end

            for iIn = 1:numel(nodeLocal.info.inputs)
                inDef = nodeLocal.info.inputs(iIn);

                if ~isfield(inDef, 'isData') || ~inDef.isData
                    continue
                end

                connIdx = find([app.pm.connections.targetNodeID] == nodeLocal.id & ...
                    strcmpi({app.pm.connections.targetInputName}, char(string(inDef.name))), 1, 'first');

                if isempty(connIdx)
                    outLines(end+1, 1) = "- " + string(inDef.name) + " <= <missing>"; %#ok<AGROW>
                    continue
                end

                connLocal = app.pm.connections(connIdx);
                srcIdx = find([app.pm.nodes.id] == connLocal.sourceNodeID, 1, 'first');
                srcNode = app.pm.nodes(srcIdx);

                srcRef = char(string(connLocal.sourceOutputName));
                if isfield(connLocal, 'selectedFile') && ~isempty(connLocal.selectedFile)
                    srcRef = char(string(connLocal.selectedFile));
                end

                outLines(end+1, 1) = "- " + string(inDef.name) + " <= " + string(srcNode.name) + ":" + string(srcRef); %#ok<AGROW>
            end

            if isempty(outLines)
                outLines(end+1, 1) = "<no data inputs>"; %#ok<AGROW>
            end

            outLines = cellstr(outLines);
        end

        function outLines = getOutputLines(app, nodeLocal)
            %GETOUTPUTLINES Build the output lines shown in the selected-step summary.

            outLines = strings(0, 1);

            if ~isfield(nodeLocal, 'info') || ~isfield(nodeLocal.info, 'outputs') || isempty(nodeLocal.info.outputs)
                outLines(end+1, 1) = "<no outputs>"; %#ok<AGROW>
                outLines = cellstr(outLines);
                return
            end

            for iOut = 1:numel(nodeLocal.info.outputs)
                outDef = nodeLocal.info.outputs(iOut);
                outName = char(string(outDef.name));
                outType = app.typeListToText(outDef.type);

                outMode = 'data';
                if isfield(outDef, 'outputMode') && ~isempty(outDef.outputMode)
                    outMode = char(string(outDef.outputMode));
                end

                outLines(end+1, 1) = "- " + string(outName) + " (" + string(outType) + ", " + string(outMode) + ")"; %#ok<AGROW>
            end

            outLines = cellstr(outLines);
        end

        function outLines = getSavedOutputLines(app, nodeLocal)
            %GETSAVEDOUTPUTLINES Return manager-save-target lines for one step.
            %
            % Only saveable outputs are listed:
            %   isData == true AND outputMode ~= 'file'

            outLines = strings(0, 1);

            if ~isfield(nodeLocal, 'info') || ~isfield(nodeLocal.info, 'outputs') || isempty(nodeLocal.info.outputs)
                outLines = cellstr(outLines);
                return
            end

            for iOut = 1:numel(nodeLocal.info.outputs)
                outDef = nodeLocal.info.outputs(iOut);

                outMode = 'data';
                if isfield(outDef, 'outputMode') && ~isempty(outDef.outputMode)
                    outMode = char(string(outDef.outputMode));
                end

                isSaveable = isfield(outDef, 'isData') && outDef.isData && ~strcmpi(outMode, 'file');

                if isSaveable && isfield(outDef, 'saveFileName') && ~isempty(outDef.saveFileName)
                    try
                        if ischar(outDef.saveFileName) || (isstring(outDef.saveFileName) && isscalar(outDef.saveFileName))
                            saveNames = string(outDef.saveFileName);
                        elseif isstring(outDef.saveFileName)
                            saveNames = outDef.saveFileName(:);
                        elseif iscell(outDef.saveFileName)
                            saveNames = string(outDef.saveFileName(:));
                        else
                            saveNames = string(outDef.saveFileName);
                            saveNames = saveNames(:);
                        end
                    catch
                        saveNames = strings(0,1);
                    end

                    saveNames = strip(saveNames);
                    saveNames = saveNames(strlength(saveNames) > 0);
                    saveNames = saveNames(~ismissing(saveNames));

                    for iName = 1:numel(saveNames)
                        outLines(end+1, 1) = "- " + string(outDef.name) + " -> " + saveNames(iName); %#ok<AGROW>
                    end
                end
            end

            outLines = cellstr(outLines);
        end

        function outLines = buildParameterText(app, nodeLocal)
            %BUILDPARAMETERTEXT Build the Parameters text for one selected node.

            outLines = strings(0, 1);

            if ~isfield(nodeLocal, 'info') || ~isfield(nodeLocal.info, 'parameters') || isempty(nodeLocal.info.parameters)
                outLines(end+1, 1) = "<no parameters>"; %#ok<AGROW>
                outLines = cellstr(outLines);
                return
            end

            for iParam = 1:numel(nodeLocal.info.parameters)
                p = nodeLocal.info.parameters(iParam);

                if ~isfield(p, 'name') || isempty(p.name)
                    continue
                end

                pValue = [];
                if isfield(p, 'value')
                    pValue = p.value;
                elseif isfield(p, 'default')
                    pValue = p.default;
                end

                outLines(end+1, 1) = string(p.name) + ": " + string(app.formatValueForDisplay(pValue)); %#ok<AGROW>
            end

            if isempty(outLines)
                outLines(end+1, 1) = "<no parameters>"; %#ok<AGROW>
            end

            outLines = cellstr(outLines);
        end

        function txt = formatValueForDisplay(~, valueIn)
            %FORMATVALUEFORDISPLAY Convert a MATLAB value to short readable text.

            if isempty(valueIn)
                txt = '[]';
                return
            end

            if ischar(valueIn) || (isstring(valueIn) && isscalar(valueIn))
                txt = char(string(valueIn));
                return
            end

            if isnumeric(valueIn) || islogical(valueIn)
                if isscalar(valueIn) || (isvector(valueIn) && numel(valueIn) <= 6)
                    txt = mat2str(valueIn);
                else
                    txt = sprintf('[%s %s]', strjoin(string(size(valueIn)), 'x'), class(valueIn));
                end
                return
            end

            if iscell(valueIn)
                txt = sprintf('[%s cell]', strjoin(string(size(valueIn)), 'x'));
                return
            end

            if isstruct(valueIn)
                txt = sprintf('[%s struct]', strjoin(string(size(valueIn)), 'x'));
                return
            end

            if isa(valueIn, 'function_handle')
                txt = ['@' func2str(valueIn)];
                return
            end

            txt = class(valueIn);
        end

        function txt = getNodeFunctionName(~, nodeLocal)
            %GETNODEFUNCTIONNAME Return the callable function name for one node.

            if isfield(nodeLocal, 'kind') && strcmpi(nodeLocal.kind, 'folder')
                txt = 'File Source';
                return
            end

            txt = char(string(nodeLocal.name));

            if isfield(nodeLocal, 'info') && isfield(nodeLocal.info, 'name') && ~isempty(nodeLocal.info.name)
                txt = char(string(nodeLocal.info.name));
            end
        end

        function txt = typeListToText(~, typeIn)
            %TYPELISTTOTEXT Convert type metadata to one display string.

            if isempty(typeIn)
                txt = 'Unknown';
                return
            end

            if ischar(typeIn) || (isstring(typeIn) && isscalar(typeIn))
                txt = char(string(typeIn));
                return
            end

            txt = strjoin(string(typeIn), ' / ');
            txt = char(txt);
        end

        function outVal = valueOrPlaceholder(~, val)
            %VALUEORPLACEHOLDER Return a string-friendly value or placeholder.

            if isempty(val)
                outVal = '<n/a>';
            else
                outVal = val;
            end
        end

        function appendDiagnostic(app, msg)
            %APPENDDIAGNOSTIC Append one timestamped line to the activity log.

            stamp = char(datetime('now', 'Format', 'HH:mm:ss'));
            newLine = sprintf('[%s] %s', stamp, msg);

            currentVal = string(app.ActivityLogTextArea.Value);
            if isempty(currentVal) || (numel(currentVal) == 1 && strlength(currentVal(1)) == 0)
                app.ActivityLogTextArea.Value = {newLine};
            else
                app.ActivityLogTextArea.Value = cellstr([currentVal; string(newLine)]);
            end
        end

        function setStatus(app, msg)
            %SETSTATUS Update the bottom status label.

            app.StatusLabel.Text = ['Status: ' char(string(msg))];
        end

        function funcName = getSelectedFunctionName(app)
            %GETSELECTEDFUNCTIONNAME Return the currently selected function-tree leaf name.

            funcName = '';

            selectedNodes = app.FunctionTree.SelectedNodes;
            if isempty(selectedNodes)
                return
            end

            selectedNode = selectedNodes(1);
            nodeData = selectedNode.NodeData;

            if isempty(nodeData) || ~isstruct(nodeData) || ~isfield(nodeData, 'kind')
                return
            end

            if strcmpi(nodeData.kind, 'function') && isfield(nodeData, 'name')
                funcName = char(string(nodeData.name));
            end
        end

        function tryAddSelectedFunction(app, isNewBranch)
            %TRYADDSELECTEDFUNCTION Add the selected function using GUI-resolved inputs.
            %
            %   tryAddSelectedFunction(app, false) -> standard add
            %   tryAddSelectedFunction(app, true)  -> add as new branch
            %
            % Standard add:
            %   - selected-node-first continuation
            %   - fallback to compatible leaf sources
            %
            % New-branch add:
            %   - always starts from a freshly selected file input
            %   - does not auto-use existing FileSource nodes
            %
            % In both modes, all DATA-input steps are added using explicit
            % 'input' references to avoid MATLAB CLI prompting.

            if nargin < 2 || isempty(isNewBranch)
                isNewBranch = false;
            end

            if isempty(app.pm)
                app.appendDiagnostic('No PipelineManager is attached.');
                app.setStatus('No PipelineManager attached.');
                return
            end

            funcName = app.getSelectedFunctionName();
            if isempty(funcName)
                app.appendDiagnostic('No function is selected.');
                app.setStatus('No function selected.');
                return
            end

            try
                info = app.pm.getPipelineInfo(funcName);

                dataInputs = struct([]);
                if isfield(info, 'inputs') && ~isempty(info.inputs)
                    dataMask = arrayfun(@(x) isfield(x, 'isData') && x.isData, info.inputs);
                    dataInputs = info.inputs(dataMask);
                end

                nDataInputs = numel(dataInputs);

                % -------------------------------------------------------------
                % No DATA inputs -> safe direct add
                % -------------------------------------------------------------
                if nDataInputs == 0
                    stepTag = app.pm.addStep(funcName);

                else
                    inputRefs = app.promptExplicitInputsForAdd(funcName, dataInputs, isNewBranch);

                    if isempty(inputRefs)
                        if isNewBranch
                            app.appendDiagnostic(sprintf('Add-as-new-branch operation for "%s" was cancelled.', funcName));
                            app.setStatus('Add-as-new-branch cancelled.');
                        else
                            app.appendDiagnostic(sprintf('Add-step operation for "%s" was cancelled.', funcName));
                            app.setStatus('Add-step cancelled.');
                        end
                        return
                    end

                    if nDataInputs == 1
                        stepTag = app.pm.addStep(funcName, 'input', inputRefs{1});
                    else
                        stepTag = app.pm.addStep(funcName, 'input', inputRefs);
                    end
                end

                if isempty(stepTag)
                    if isNewBranch
                        app.appendDiagnostic(sprintf('Add-as-new-branch operation for "%s" was cancelled.', funcName));
                        app.setStatus('Add-as-new-branch cancelled.');
                    else
                        app.appendDiagnostic(sprintf('Add-step operation for "%s" was cancelled.', funcName));
                        app.setStatus('Add-step cancelled.');
                    end
                    return
                end

                newIdx = find(strcmpi(stepTag, {app.pm.nodes.name}), 1, 'first');
                if ~isempty(newIdx)
                    app.selectedNodeID = app.pm.nodes(newIdx).id;
                    app.selectedStepTag = char(string(app.pm.nodes(newIdx).name));
                else
                    app.selectedNodeID = NaN;
                    app.selectedStepTag = '';
                end

                app.refreshAllViews();

                if isNewBranch
                    app.appendDiagnostic(sprintf('Added function "%s" as new-branch step "%s".', funcName, stepTag));
                    app.setStatus('New branch step added.');
                else
                    app.appendDiagnostic(sprintf('Added function "%s" as step "%s".', funcName, stepTag));
                    app.setStatus('Step added.');
                end

            catch ME
                if isNewBranch
                    app.appendDiagnostic(sprintf('Failed to add "%s" as new branch: %s', funcName, ME.message));
                    app.setStatus('Failed to add new branch step.');
                else
                    app.appendDiagnostic(sprintf('Failed to add function "%s": %s', funcName, ME.message));
                    app.setStatus('Failed to add step.');
                end
            end
        end

        function inputRefs = promptExplicitInputsForAdd(app, funcName, dataInputs, isNewBranch)
            %PROMPTEXPLICITINPUTSFORADD Resolve explicit input references for addStep.
            %
            % Standard mode:
            %   1) selected node only
            %   2) compatible leaf nodes
            %   3) browse file
            %
            % New-branch mode:
            %   - always prompt for a new file selection for each DATA input
            %   - do not auto-use existing FileSource nodes
            %
            % Returned refs are suitable for:
            %   PipelineManager.addStep(..., 'input', inputRefs)
            %
            % Important:
            %   FileSource/folder nodes return full file paths, not "StepTag:data".
            %   This forces PipelineManager.resolveInputReference to use filename
            %   semantics and populate connection.selectedFile.

            if nargin < 4 || isempty(isNewBranch)
                isNewBranch = false;
            end

            inputRefs = [];

            if nargin < 3 || isempty(dataInputs)
                inputRefs = {};
                return
            end

            inputRefs = cell(1, numel(dataInputs));

            for iIn = 1:numel(dataInputs)

                inputDef = dataInputs(iIn);
                inputName = char(string(inputDef.name));
                chosenRef = '';

                % -------------------------------------------------------------
                % MODE A: Add as new branch
                % -------------------------------------------------------------
                if isNewBranch
                    chosenRef = app.promptFileInputReference(inputName);

                    if isempty(chosenRef)
                        inputRefs = [];
                        return
                    end

                    app.appendDiagnostic(sprintf( ...
                        'Selected new file source "%s" for input "%s" of "%s".', ...
                        chosenRef, inputName, funcName));

                    inputRefs{iIn} = chosenRef;
                    continue
                end

                % -------------------------------------------------------------
                % MODE B: Standard add, selected-node-first
                % -------------------------------------------------------------
                [selectedDisplay, selectedRefs, selectedDefaultIdx] = app.getAddChoicesForInput( ...
                    inputDef, 'selectedNodeOnly', true);

                if ~isempty(selectedRefs)

                    if numel(selectedRefs) == 1
                        chosenRef = selectedRefs{1};

                        app.appendDiagnostic(sprintf( ...
                            ['Auto-selected source "%s" for input "%s" of "%s" ' ...
                            'from the currently selected step.'], ...
                            chosenRef, inputName, funcName));

                        inputRefs{iIn} = chosenRef;
                        continue
                    end

                    chosenRef = app.promptSourceChoiceDialog( ...
                        inputName, funcName, selectedDisplay, selectedRefs, selectedDefaultIdx);

                    if isempty(chosenRef)
                        inputRefs = [];
                        return
                    end

                    app.appendDiagnostic(sprintf( ...
                        ['Selected source "%s" for input "%s" of "%s" ' ...
                        'from the currently selected step.'], ...
                        chosenRef, inputName, funcName));

                    inputRefs{iIn} = chosenRef;
                    continue
                end

                % -------------------------------------------------------------
                % Second pass: compatible leaf sources
                % -------------------------------------------------------------
                [displayList, refList, defaultIdx] = app.getAddChoicesForInput( ...
                    inputDef, 'selectedNodeOnly', false);

                if isempty(refList)
                    chosenRef = app.promptFileInputReference(inputName);

                    if isempty(chosenRef)
                        inputRefs = [];
                        return
                    end

                    app.appendDiagnostic(sprintf( ...
                        'Selected file "%s" for input "%s" of "%s".', ...
                        chosenRef, inputName, funcName));

                    inputRefs{iIn} = chosenRef;
                    continue
                end

                if numel(refList) == 1
                    chosenRef = refList{1};

                    app.appendDiagnostic(sprintf( ...
                        'Auto-selected source "%s" for input "%s" of "%s".', ...
                        chosenRef, inputName, funcName));

                    inputRefs{iIn} = chosenRef;
                    continue
                end

                chosenRef = app.promptSourceChoiceDialog( ...
                    inputName, funcName, displayList, refList, defaultIdx);

                if isempty(chosenRef)
                    inputRefs = [];
                    return
                end

                app.appendDiagnostic(sprintf( ...
                    'Selected source "%s" for input "%s" of "%s".', ...
                    chosenRef, inputName, funcName));

                inputRefs{iIn} = chosenRef;
            end
        end

        function [displayList, refList, defaultIdx] = getAddChoicesForInput(app, inputDef, varargin)
            %GETADDCHOICESFORINPUT Build explicit source choices for one add-step input.
            %
            % Standard add policy:
            %   - selectedNodeOnly=true: only the currently selected graph node is used
            %   - selectedNodeOnly=false: compatible leaf nodes are used
            %
            % Source reference policy:
            %   - FileSource/folder nodes:
            %         fullfile(SaveFolder, filename)
            %     so PipelineManager resolves it as a filename and sets selectedFile.
            %
            %   - Stream data outputs:
            %         StepTag:OutputPortName
            %
            %   - Stream file outputs:
            %         StepTag:filename.ext

            p = inputParser;
            addParameter(p, 'selectedNodeOnly', false, @(x) islogical(x) && isscalar(x));
            parse(p, varargin{:});

            selectedNodeOnly = p.Results.selectedNodeOnly;

            displayList = {};
            refList = {};
            defaultIdx = [];

            dstTypes = {};
            if isfield(inputDef, 'type') && ~isempty(inputDef.type)
                if ischar(inputDef.type) || (isstring(inputDef.type) && isscalar(inputDef.type))
                    dstTypes = {char(string(inputDef.type))};
                else
                    dstTypes = cellstr(string(inputDef.type(:)));
                end
            end

            if selectedNodeOnly
                if isnan(app.selectedNodeID)
                    return
                end

                [displayList, refList, defaultIdx] = app.buildExplicitSourceChoices(dstTypes, ...
                    'leafOnly', false, ...
                    'candidateNodeIDs', app.selectedNodeID);
            else
                [displayList, refList, defaultIdx] = app.buildExplicitSourceChoices(dstTypes, ...
                    'leafOnly', true, ...
                    'preferredNodeID', app.selectedNodeID);
            end
        end

        function [displayList, refList, defaultIdx] = buildExplicitSourceChoices(app, dstTypes, varargin)
            %BUILDEXPLICITSOURCECHOICES Build explicit selectable source references.
            %
            % Output:
            %   displayList - Cellstr shown in list dialogs.
            %   refList     - Explicit references suitable for addStep/reconnectStep.
            %
            % Reference policy:
            %   FileSource/folder node:
            %       fullfile(SaveFolder, filename)
            %
            %       This avoids ambiguity with an existing graph node named "red.dat".
            %       PipelineManager.resolveInputReference normalizes this path to
            %       selectedFile="red.dat".
            %
            %   Stream data output:
            %       StepTag:OutputPortName
            %
            %   Stream file output:
            %       StepTag:filename.ext

            p = inputParser;
            addParameter(p, 'leafOnly', false, @(x) islogical(x) && isscalar(x));
            addParameter(p, 'excludeNodeID', NaN, @(x) isnumeric(x) && isscalar(x));
            addParameter(p, 'preferredNodeID', NaN, @(x) isnumeric(x) && isscalar(x));
            addParameter(p, 'candidateNodeIDs', [], @(x) isnumeric(x) || isempty(x));
            parse(p, varargin{:});

            leafOnly = p.Results.leafOnly;
            excludeNodeID = p.Results.excludeNodeID;
            preferredNodeID = p.Results.preferredNodeID;
            candidateNodeIDs = p.Results.candidateNodeIDs;

            displayList = {};
            refList = {};
            defaultIdx = [];

            if isempty(app.pm) || isempty(app.pm.nodes)
                return
            end

            candidateNodes = app.pm.nodes;

            % -------------------------------------------------------------
            % Restrict to explicit node IDs or leaves
            % -------------------------------------------------------------
            if ~isempty(candidateNodeIDs)
                candidateNodeIDs = unique(candidateNodeIDs(:).', 'stable');
                candidateNodes = candidateNodes(ismember([candidateNodes.id], candidateNodeIDs));
            elseif leafOnly
                leafIDs = app.computeLeafNodeIDs(app.pm.nodes, app.pm.connections);
                candidateNodes = candidateNodes(ismember([candidateNodes.id], leafIDs));
            end

            if ~isnan(excludeNodeID)
                candidateNodes = candidateNodes([candidateNodes.id] ~= excludeNodeID);
            end

            if isempty(candidateNodes)
                return
            end

            % -------------------------------------------------------------
            % Prefer selected node when applicable
            % -------------------------------------------------------------
            if isempty(candidateNodeIDs) && ~isnan(preferredNodeID) && any([candidateNodes.id] == preferredNodeID)
                selMask = [candidateNodes.id] == preferredNodeID;
                candidateNodes = [candidateNodes(selMask), candidateNodes(~selMask)];
            end

            saveFolder = '';
            try
                if ~isempty(app.pm.SaveFolderList)
                    saveFolder = app.pm.SaveFolderList{1};
                end
            catch
                saveFolder = '';
            end

            for iNode = 1:numel(candidateNodes)

                srcNode = candidateNodes(iNode);

                if ~isfield(srcNode, 'info') || ~isfield(srcNode.info, 'outputs') || isempty(srcNode.info.outputs)
                    continue
                end

                % -----------------------------------------------------------------
                % FileSource/folder node special case
                % -----------------------------------------------------------------
                if isfield(srcNode, 'kind') && strcmpi(char(string(srcNode.kind)), 'folder')

                    sourceFileName = char(string(srcNode.name));

                    % Prefer declared filename when present.
                    try
                        outDef0 = srcNode.info.outputs(1);
                        if isfield(outDef0, 'defOutfilename') && ~isempty(outDef0.defOutfilename)
                            if iscell(outDef0.defOutfilename)
                                defList = outDef0.defOutfilename;
                                defList = defList(~cellfun(@isempty, defList));
                                if ~isempty(defList)
                                    sourceFileName = char(string(defList{1}));
                                end
                            else
                                sourceFileName = char(string(outDef0.defOutfilename));
                            end
                        end
                    catch
                    end

                    [~, sourceBase, sourceExt] = fileparts(sourceFileName);
                    sourceFileName = [sourceBase sourceExt];

                    if isempty(sourceFileName)
                        continue
                    end

                    outDef = srcNode.info.outputs(1);

                    srcTypes = {};
                    if isfield(outDef, 'type') && ~isempty(outDef.type)
                        if ischar(outDef.type) || (isstring(outDef.type) && isscalar(outDef.type))
                            srcTypes = {char(string(outDef.type))};
                        else
                            srcTypes = cellstr(string(outDef.type(:)));
                        end
                    end

                    if ~app.areTypesCompatible(srcTypes, dstTypes)
                        continue
                    end

                    typeText = app.typeListToText(srcTypes);

                    if isempty(saveFolder)
                        sourceRef = sourceFileName;
                    else
                        sourceRef = fullfile(saveFolder, sourceFileName);
                    end

                    displayList{end+1,1} = sprintf('%s [File Source, %s]', ...
                        sourceFileName, typeText); %#ok<AGROW>

                    refList{end+1,1} = sourceRef; %#ok<AGROW>

                    continue
                end

                % -----------------------------------------------------------------
                % Normal stream outputs
                % -----------------------------------------------------------------
                for iOut = 1:numel(srcNode.info.outputs)

                    outDef = srcNode.info.outputs(iOut);

                    if ~isfield(outDef, 'isData') || ~outDef.isData
                        continue
                    end

                    srcTypes = {};
                    if isfield(outDef, 'type') && ~isempty(outDef.type)
                        if ischar(outDef.type) || (isstring(outDef.type) && isscalar(outDef.type))
                            srcTypes = {char(string(outDef.type))};
                        else
                            srcTypes = cellstr(string(outDef.type(:)));
                        end
                    end

                    if ~app.areTypesCompatible(srcTypes, dstTypes)
                        continue
                    end

                    outName = char(string(outDef.name));
                    outMode = 'data';
                    if isfield(outDef, 'outputMode') && ~isempty(outDef.outputMode)
                        outMode = char(string(outDef.outputMode));
                    end

                    typeText = app.typeListToText(srcTypes);

                    if strcmpi(outMode, 'file')

                        declaredFiles = {};
                        if isfield(outDef, 'defOutfilename') && ~isempty(outDef.defOutfilename)
                            if ischar(outDef.defOutfilename) || ...
                                    (isstring(outDef.defOutfilename) && isscalar(outDef.defOutfilename))
                                declaredFiles = {char(string(outDef.defOutfilename))};
                            elseif iscell(outDef.defOutfilename)
                                declaredFiles = cellfun(@(x) char(string(x)), ...
                                    outDef.defOutfilename, 'UniformOutput', false);
                            elseif isstring(outDef.defOutfilename)
                                declaredFiles = cellstr(outDef.defOutfilename(:));
                            end
                        end

                        declaredFiles = declaredFiles(~cellfun(@isempty, declaredFiles));
                        declaredFiles = unique(declaredFiles, 'stable');

                        if isempty(declaredFiles)
                            displayList{end+1,1} = sprintf('%s:%s [%s, file]', ...
                                srcNode.name, outName, typeText); %#ok<AGROW>
                            refList{end+1,1} = sprintf('%s:%s', ...
                                srcNode.name, outName); %#ok<AGROW>
                        else
                            for iFile = 1:numel(declaredFiles)
                                displayList{end+1,1} = sprintf('%s:%s [%s, file]', ...
                                    srcNode.name, declaredFiles{iFile}, typeText); %#ok<AGROW>
                                refList{end+1,1} = sprintf('%s:%s', ...
                                    srcNode.name, declaredFiles{iFile}); %#ok<AGROW>
                            end
                        end

                    else
                        displayList{end+1,1} = sprintf('%s:%s [%s]', ...
                            srcNode.name, outName, typeText); %#ok<AGROW>
                        refList{end+1,1} = sprintf('%s:%s', ...
                            srcNode.name, outName); %#ok<AGROW>
                    end
                end
            end

            if isempty(displayList)
                return
            end

            defaultIdx = 1;
        end

        function chosenRef = promptSourceChoiceDialog(app, inputName, funcName, displayList, refList, defaultIdx)
            %PROMPTSOURCECHOICEDIALOG Show the GUI dialog for ambiguous source selection.
            %
            % Output:
            %   chosenRef - Explicit source reference string, or '' if cancelled.

            chosenRef = '';

            if isempty(displayList) || isempty(refList)
                return
            end

            browseToken = '__BROWSE_FILE__';
            promptDisplay = [displayList(:); {sprintf('[Browse SaveFolder file...]  -> %s', inputName)}];
            promptRefs    = [refList(:); {browseToken}];

            if isempty(defaultIdx) || defaultIdx < 1 || defaultIdx > numel(displayList)
                defaultIdx = 1;
            end

            [selIdx, tf] = listdlg( ...
                'PromptString', sprintf('Select source for input "%s" of "%s":', inputName, funcName), ...
                'SelectionMode', 'single', ...
                'ListString', promptDisplay, ...
                'InitialValue', defaultIdx, ...
                'ListSize', [700 320]);

            if ~tf || isempty(selIdx)
                return
            end

            chosenRef = promptRefs{selIdx};

            if strcmp(chosenRef, browseToken)
                chosenRef = app.promptFileInputReference(inputName);
            end
        end

        function [displayList, refList] = getReconnectChoices(app, targetNode, inputDef)
            %GETRECONNECTCHOICES Build explicit reconnect choices for one target input.
            %
            % Reconnect policy:
            %   - Show all compatible source nodes except the target node itself
            %   - Do not restrict to leaf nodes

            displayList = {};
            refList = {};

            dstTypes = {};
            if isfield(inputDef, 'type') && ~isempty(inputDef.type)
                if ischar(inputDef.type) || (isstring(inputDef.type) && isscalar(inputDef.type))
                    dstTypes = {char(string(inputDef.type))};
                else
                    dstTypes = cellstr(string(inputDef.type(:)));
                end
            end

            [displayList, refList] = app.buildExplicitSourceChoices(dstTypes, ...
                'leafOnly', false, ...
                'excludeNodeID', targetNode.id);
        end

        function inputRef = promptFileInputReference(app, inputName)
            %PROMPTFILEINPUTREFERENCE Select an existing file from SaveFolder for one input.
            %
            % Output:
            %   inputRef - Full file path suitable for PipelineManager.resolveInputReference,
            %              or '' if cancelled.
            %
            % Returning the full path avoids ambiguity when a FileSource node already exists
            % with the same name as the file, for example node.name == "red.dat".

            inputRef = '';

            startFolder = pwd;
            if ~isempty(app.pm) && isprop(app.pm, 'SaveFolderList') && ~isempty(app.pm.SaveFolderList)
                startFolder = app.pm.SaveFolderList{1};
            end

            [fileName, folderName] = uigetfile('*.*', ...
                sprintf('Select file for input "%s"', inputName), ...
                startFolder);

            if isequal(fileName, 0)
                return
            end

            if isempty(folderName)
                folderName = startFolder;
            end

            inputRef = fullfile(folderName, fileName);

            app.appendDiagnostic(sprintf('Selected file "%s" for input "%s".', inputRef, inputName));
        end

        function reconnectStepInteractive(app, stepTag)
            %RECONNECTSTEPINTERACTIVE Start reconnect mode for one selected step.

            if isempty(app.pm)
                app.appendDiagnostic('Reconnect requested without an attached PipelineManager.');
                app.setStatus('No PipelineManager attached.');
                return
            end

            nodeIdx = find(strcmpi(stepTag, {app.pm.nodes.name}), 1, 'first');
            if isempty(nodeIdx)
                app.appendDiagnostic(sprintf('Reconnect failed. Step "%s" was not found.', stepTag));
                app.setStatus('Reconnect failed.');
                return
            end

            nodeLocal = app.pm.nodes(nodeIdx);

            if strcmpi(nodeLocal.kind, 'folder')
                app.appendDiagnostic('File source steps cannot be reconnected.');
                app.setStatus('Reconnect not allowed.');
                return
            end

            if ~isfield(nodeLocal, 'info') || ~isfield(nodeLocal.info, 'inputs') || isempty(nodeLocal.info.inputs)
                app.appendDiagnostic(sprintf('Step "%s" has no inputs to reconnect.', stepTag));
                app.setStatus('No inputs to reconnect.');
                return
            end

            dataInputs = nodeLocal.info.inputs(arrayfun(@(x) isfield(x, 'isData') && x.isData, nodeLocal.info.inputs));

            if isempty(dataInputs)
                app.appendDiagnostic(sprintf('Step "%s" has no data inputs to reconnect.', stepTag));
                app.setStatus('No data inputs to reconnect.');
                return
            end

            if numel(dataInputs) == 1
                targetInputName = char(string(dataInputs(1).name));
            else
                inputLabels = cell(numel(dataInputs), 1);
                for iIn = 1:numel(dataInputs)
                    inputLabels{iIn} = sprintf('%s [%s]', ...
                        char(string(dataInputs(iIn).name)), ...
                        app.typeListToText(dataInputs(iIn).type));
                end

                [selIdx, tf] = listdlg( ...
                    'PromptString', sprintf('Select input to reconnect for "%s":', stepTag), ...
                    'SelectionMode', 'single', ...
                    'ListString', inputLabels, ...
                    'ListSize', [450 260]);

                if ~tf || isempty(selIdx)
                    app.appendDiagnostic('Reconnect cancelled by user.');
                    app.setStatus('Reconnect cancelled.');
                    return
                end

                targetInputName = char(string(dataInputs(selIdx).name));
            end

            targetInputDef = [];
            for iIn = 1:numel(dataInputs)
                if strcmpi(char(string(dataInputs(iIn).name)), targetInputName)
                    targetInputDef = dataInputs(iIn);
                    break
                end
            end

            if isempty(targetInputDef)
                app.appendDiagnostic('Reconnect failed. Target input definition was not found.');
                app.setStatus('Reconnect failed.');
                return
            end

            [~, refList] = app.getReconnectChoices(nodeLocal, targetInputDef);
            candidateNodeIDs = app.refListToNodeIDs(refList);

            if isempty(candidateNodeIDs)
                app.appendDiagnostic(sprintf('No compatible source steps were found for input "%s".', targetInputName));
                app.setStatus('Reconnect aborted.');
                return
            end

            app.bReconnectMode = true;
            app.reconnectTargetNodeID = nodeLocal.id;
            app.reconnectTargetStepTag = stepTag;
            app.reconnectTargetInputName = targetInputName;
            app.reconnectCandidateNodeIDs = candidateNodeIDs;

            app.selectedNodeID = nodeLocal.id;
            app.selectedStepTag = stepTag;

            app.refreshGraphView();
            app.appendDiagnostic(sprintf('Reconnect mode started for input "%s" of "%s".', targetInputName, stepTag));
            app.setStatus(sprintf('Reconnect mode: click a new source for %s', targetInputName));
        end

        function inputRef = getCurrentExplicitInputReference(app, targetNode, inputDef)
            %GETCURRENTEXPLICITINPUTREFERENCE Return the current explicit input reference string.

            inputRef = '';

            if isempty(app.pm) || isempty(app.pm.connections)
                return
            end

            inputName = char(string(inputDef.name));

            connIdx = find([app.pm.connections.targetNodeID] == targetNode.id & ...
                strcmpi({app.pm.connections.targetInputName}, inputName), 1, 'first');

            if isempty(connIdx)
                return
            end

            connLocal = app.pm.connections(connIdx);
            srcIdx = find([app.pm.nodes.id] == connLocal.sourceNodeID, 1, 'first');
            if isempty(srcIdx)
                return
            end

            srcNode = app.pm.nodes(srcIdx);

            srcRef = char(string(connLocal.sourceOutputName));
            if isfield(connLocal, 'selectedFile') && ~isempty(connLocal.selectedFile)
                srcRef = char(string(connLocal.selectedFile));
            end

            inputRef = sprintf('%s:%s', char(string(srcNode.name)), srcRef);
        end

        function nodeIDs = refListToNodeIDs(app, refList)
            %REFLISTTONODEIDS Convert explicit source refs to unique node IDs.

            nodeIDs = [];

            if isempty(refList) || isempty(app.pm) || isempty(app.pm.nodes)
                return
            end

            nodeIDs = zeros(0,1);

            for iRef = 1:numel(refList)
                refStr = char(string(refList{iRef}));
                tok = regexp(refStr, '^([^:]+):', 'tokens', 'once');
                if isempty(tok)
                    continue
                end

                stepTag = tok{1};
                idx = find(strcmpi(stepTag, {app.pm.nodes.name}), 1, 'first');
                if ~isempty(idx)
                    nodeIDs(end+1,1) = app.pm.nodes(idx).id; %#ok<AGROW>
                end
            end

            nodeIDs = unique(nodeIDs(:).', 'stable');
        end

        function exitReconnectMode(app)
            %EXITRECONNECTMODE Leave reconnect mode and clear transient state.

            app.bReconnectMode = false;
            app.reconnectTargetNodeID = NaN;
            app.reconnectTargetStepTag = '';
            app.reconnectTargetInputName = '';
            app.reconnectCandidateNodeIDs = [];

            app.refreshGraphView();
        end

        function tf = areTypesCompatible(~, srcTypes, dstTypes)
            %ARETYPESCOMPATIBLE Return true when two type lists are compatible.

            if isempty(srcTypes) || isempty(dstTypes)
                tf = true;
                return
            end

            srcTypes = cellstr(string(srcTypes(:)));
            dstTypes = cellstr(string(dstTypes(:)));

            srcHasWildcard = any(strcmp(srcTypes, 'UnknownDataType'));
            dstHasWildcard = any(strcmp(dstTypes, 'UnknownDataType'));

            tf = srcHasWildcard || dstHasWildcard || ~isempty(intersect(srcTypes, dstTypes));
        end

        function configureSaveTargetsInteractive(app, stepTag)
            %CONFIGURESAVETARGETSINTERACTIVE Configure PipelineManager-managed save targets.

            if isempty(app.pm)
                app.appendDiagnostic('Configure-save requested without an attached PipelineManager.');
                app.setStatus('No PipelineManager attached.');
                return
            end

            nodeIdx = find(strcmpi(stepTag, {app.pm.nodes.name}), 1, 'first');
            if isempty(nodeIdx)
                app.appendDiagnostic(sprintf('Configure-save failed. Step "%s" was not found.', stepTag));
                app.setStatus('Configure-save failed.');
                return
            end

            nodeLocal = app.pm.nodes(nodeIdx);

            if strcmpi(nodeLocal.kind, 'folder')
                app.appendDiagnostic('File source steps cannot have PipelineManager-managed save targets.');
                app.setStatus('Configure-save not allowed.');
                return
            end

            saveableOutputs = app.getSaveableOutputDefsForNode(nodeLocal);

            if isempty(saveableOutputs)
                app.appendDiagnostic(sprintf('Step "%s" has no manager-saveable outputs.', stepTag));
                app.setStatus('No saveable outputs.');
                return
            end

            promptList = cell(numel(saveableOutputs), 1);
            defaultVals = cell(numel(saveableOutputs), 1);

            for iOut = 1:numel(saveableOutputs)
                outDef = saveableOutputs(iOut);

                outName = char(string(outDef.name));
                outType = app.typeListToText(outDef.type);

                defaultTarget = '';
                if isfield(outDef, 'saveFileName') && ~isempty(outDef.saveFileName)
                    defaultTarget = char(string(outDef.saveFileName));
                elseif isfield(outDef, 'defOutfilename') && ~isempty(outDef.defOutfilename)
                    if iscell(outDef.defOutfilename)
                        defaultTarget = char(string(outDef.defOutfilename{1}));
                    else
                        defaultTarget = char(string(outDef.defOutfilename));
                    end
                end

                promptList{iOut} = sprintf('Save filename for output "%s" [%s]:', outName, outType);
                defaultVals{iOut} = defaultTarget;
            end

            answer = inputdlg( ...
                promptList, ...
                sprintf('Save Outputs - %s', stepTag), ...
                repmat([1 70], numel(promptList), 1), ...
                defaultVals);

            if isempty(answer)
                app.appendDiagnostic('Configure-save action cancelled by user.');
                app.setStatus('Configure-save cancelled.');
                return
            end

            answer = cellfun(@(x) strtrim(char(string(x))), answer, 'UniformOutput', false);
            isBlank = cellfun(@isempty, answer);

            try
                if numel(saveableOutputs) == 1

                    if isBlank(1)
                        app.pm.setSaveTargets(stepTag, '');
                        app.appendDiagnostic(sprintf('Cleared save target for "%s".', stepTag));
                        app.setStatus('Save target cleared.');
                    else
                        app.pm.setSaveTargets(stepTag, answer{1});
                        app.appendDiagnostic(sprintf('Updated save target for "%s" -> %s', stepTag, answer{1}));
                        app.setStatus('Save target updated.');
                    end

                else
                    if all(isBlank)
                        app.pm.setSaveTargets(stepTag, '');
                        app.appendDiagnostic(sprintf('Cleared all save targets for "%s".', stepTag));
                        app.setStatus('Save targets cleared.');

                    elseif any(isBlank)
                        app.appendDiagnostic(['Partial save configuration is not supported yet for multi-output steps. ' ...
                            'Fill every saveable output filename or leave all blank to clear all.']);
                        app.setStatus('Invalid save-target selection.');
                        return

                    else
                        app.pm.setSaveTargets(stepTag, answer);
                        app.appendDiagnostic(sprintf('Updated %d save targets for "%s".', numel(answer), stepTag));
                        app.setStatus('Save targets updated.');
                    end
                end

                app.refreshAllViews();

            catch ME
                app.appendDiagnostic(sprintf('Failed to configure save outputs for "%s": %s', stepTag, ME.message));
                app.setStatus('Failed to configure save outputs.');
            end
        end

        function saveableOutputs = getSaveableOutputDefsForNode(~, nodeLocal)
            %GETSAVEABLEOUTPUTDEFSFORNODE Return saveable outputs for one node.
            %
            % Saveable outputs are:
            %   isData == true AND outputMode ~= 'file'

            saveableOutputs = struct([]);

            if ~isfield(nodeLocal, 'info') || ~isfield(nodeLocal.info, 'outputs') || isempty(nodeLocal.info.outputs)
                return
            end

            mask = false(1, numel(nodeLocal.info.outputs));
            for iOut = 1:numel(nodeLocal.info.outputs)

                outDef = nodeLocal.info.outputs(iOut);

                outMode = 'data';
                if isfield(outDef, 'outputMode') && ~isempty(outDef.outputMode)
                    outMode = char(string(outDef.outputMode));
                end

                mask(iOut) = isfield(outDef, 'isData') && outDef.isData && ~strcmpi(outMode, 'file');
            end

            saveableOutputs = nodeLocal.info.outputs(mask);
        end

        function completeReconnectWithRef(app, chosenRef)
            %COMPLETERECONNECTWITHREF Complete reconnect using one explicit reference.

            if ~app.bReconnectMode || isempty(app.pm)
                return
            end

            targetIdx = find([app.pm.nodes.id] == app.reconnectTargetNodeID, 1, 'first');
            if isempty(targetIdx)
                app.appendDiagnostic('Reconnect failed. Target step no longer exists.');
                app.exitReconnectMode();
                app.setStatus('Reconnect failed.');
                return
            end

            targetNode = app.pm.nodes(targetIdx);

            if ~isfield(targetNode, 'info') || ~isfield(targetNode.info, 'inputs') || isempty(targetNode.info.inputs)
                app.appendDiagnostic('Reconnect failed. Target step has no inputs.');
                app.exitReconnectMode();
                app.setStatus('Reconnect failed.');
                return
            end

            dataInputs = targetNode.info.inputs(arrayfun(@(x) isfield(x, 'isData') && x.isData, targetNode.info.inputs));

            inputRefs = cell(1, numel(dataInputs));

            for iIn = 1:numel(dataInputs)
                thisInputName = char(string(dataInputs(iIn).name));

                if strcmpi(thisInputName, app.reconnectTargetInputName)
                    inputRefs{iIn} = chosenRef;
                else
                    inputRefs{iIn} = app.getCurrentExplicitInputReference(targetNode, dataInputs(iIn));
                    if isempty(inputRefs{iIn})
                        error('Current input "%s" is not connected, so partial reconnect cannot proceed.', thisInputName);
                    end
                end
            end

            try
                targetStepTag = app.reconnectTargetStepTag;
                targetInputName = app.reconnectTargetInputName;

                app.pm.reconnectStep(targetStepTag, 'input', inputRefs);

                app.exitReconnectMode();
                app.refreshAllViews();
                app.appendDiagnostic(sprintf('Reconnected input "%s" of "%s".', targetInputName, targetStepTag));
                app.setStatus('Step reconnected.');

            catch ME
                app.appendDiagnostic(sprintf('Reconnect failed for "%s": %s', app.reconnectTargetStepTag, ME.message));
                app.setStatus('Reconnect failed.');
            end
        end

        function plan = getCurrentExecutionPlan(app)
            %GETCURRENTEXECUTIONPLAN Return the backend execution plan when available.
            %
            % The PipelineManager backend owns graph order, execution order,
            % validation status, data-flow edges, and output-persistence planning.
            % The app only falls back to a local display-only plan when the backend
            % cannot build a full plan from an already-invalid loaded graph.

            plan = [];

            if isempty(app.pm)
                return
            end

            try
                plan = app.pm.getExecutionPlan();
            catch
                plan = app.buildFallbackExecutionPlan();
            end
        end

        function plan = buildFallbackExecutionPlan(app)
            %BUILDFALLBACKEXECUTIONPLAN Build a display-only plan for invalid graphs.
            %
            % This fallback is intentionally narrow. It keeps the graph visible when
            % a loaded pipeline is structurally broken enough that getExecutionPlan
            % cannot compute backend order metadata. Normal valid/repairable graphs
            % should use PipelineManager.getExecutionPlan().

            plan = struct();
            plan.createdOn = datetime('now');
            plan.leafOutputPolicy = "";
            plan.tempOutputPrefix = "";
            plan.tempSessionTag = "";
            plan.saveFolder = "";

            if isempty(app.pm) || isempty(app.pm.nodes)
                plan.nodeOrders = struct('topoOrder', [], 'scheduleOrder', [], 'runOrder', [], 'nodeTable', table());
                plan.nodeTable = table();
                plan.edgeTable = table();
                plan.outputPlan = table();
                plan.validationReport = struct('isValid', true, 'errors', {{}}, 'warnings', {{}}, 'nodeStatus', table());
                return
            end

            try
                validationReport = app.pm.diagnosePipeline('verbose', false);
            catch
                validationReport = struct('isValid', false, 'errors', {{}}, 'warnings', {{}}, 'nodeStatus', table());
            end

            nodesLocal = app.pm.nodes;
            connectionsLocal = app.pm.connections;
            topoOrder = app.computeTopoOrderFromGraph(nodesLocal, connectionsLocal);
            scheduleOrder = topoOrder;
            runOrder = zeros(1,0);

            for iNodeID = 1:numel(scheduleOrder)
                nodeIdx = find([nodesLocal.id] == scheduleOrder(iNodeID), 1, 'first');
                if ~isempty(nodeIdx) && strcmpi(nodesLocal(nodeIdx).kind, 'stream')
                    runOrder(end+1) = scheduleOrder(iNodeID); %#ok<AGROW>
                end
            end

            nodeTable = app.buildFallbackNodeTable(nodesLocal, connectionsLocal, topoOrder, scheduleOrder, runOrder);
            edgeTable = app.buildFallbackEdgeTable(nodesLocal, connectionsLocal);

            plan.nodeOrders = struct();
            plan.nodeOrders.topoOrder = topoOrder;
            plan.nodeOrders.scheduleOrder = scheduleOrder;
            plan.nodeOrders.runOrder = runOrder;
            plan.nodeOrders.nodeTable = nodeTable;
            plan.nodeTable = nodeTable;
            plan.edgeTable = edgeTable;
            plan.outputPlan = table();
            plan.validationReport = validationReport;
        end

        function nodeTable = buildFallbackNodeTable(app, nodesLocal, connectionsLocal, topoOrder, scheduleOrder, runOrder)
            %BUILDFALLBACKNODETABLE Build a minimal backend-shaped node table.

            if nargin < 5 || isempty(scheduleOrder)
                scheduleOrder = topoOrder;
            end
            if nargin < 6
                runOrder = [];
            end

            if isempty(nodesLocal)
                nodeTable = table();
                return
            end

            nodeID = zeros(numel(nodesLocal), 1);
            stepName = strings(numel(nodesLocal), 1);
            functionName = strings(numel(nodesLocal), 1);
            kind = strings(numel(nodesLocal), 1);
            topoIndex = nan(numel(nodesLocal), 1);
            scheduleIndex = nan(numel(nodesLocal), 1);
            runIndex = nan(numel(nodesLocal), 1);
            isExecutable = false(numel(nodesLocal), 1);
            executionRole = strings(numel(nodesLocal), 1);
            isRoot = false(numel(nodesLocal), 1);
            isLeaf = false(numel(nodesLocal), 1);

            for iNode = 1:numel(nodesLocal)
                nodeLocal = nodesLocal(iNode);
                nodeID(iNode) = nodeLocal.id;
                stepName(iNode) = string(nodeLocal.name);
                functionName(iNode) = string(app.getNodeFunctionName(nodeLocal));
                kind(iNode) = string(nodeLocal.kind);
                isExecutable(iNode) = strcmpi(char(kind(iNode)), 'stream');

                if isExecutable(iNode)
                    executionRole(iNode) = "analysis_step";
                else
                    executionRole(iNode) = "virtual_source";
                end

                pos = find(topoOrder == nodeLocal.id, 1, 'first');
                if ~isempty(pos)
                    topoIndex(iNode) = pos;
                end

                pos = find(scheduleOrder == nodeLocal.id, 1, 'first');
                if ~isempty(pos)
                    scheduleIndex(iNode) = pos;
                end

                pos = find(runOrder == nodeLocal.id, 1, 'first');
                if ~isempty(pos)
                    runIndex(iNode) = pos;
                end

                if isempty(connectionsLocal)
                    isRoot(iNode) = true;
                    isLeaf(iNode) = true;
                else
                    isRoot(iNode) = ~any([connectionsLocal.targetNodeID] == nodeLocal.id);
                    isLeaf(iNode) = ~any([connectionsLocal.sourceNodeID] == nodeLocal.id);
                end
            end

            nodeTable = table(nodeID, stepName, functionName, kind, topoIndex, ...
                scheduleIndex, runIndex, isExecutable, executionRole, isRoot, isLeaf);
            [~, ord] = sort(nodeTable.topoIndex, 'ascend');
            nodeTable = nodeTable(ord, :);
        end

        function edgeTable = buildFallbackEdgeTable(app, nodesLocal, connectionsLocal)
            %BUILDFALLBACKEDGETABLE Build a minimal backend-shaped edge table.

            if isempty(connectionsLocal)
                edgeTable = table();
                return
            end

            nConn = numel(connectionsLocal);
            sourceNodeID = zeros(nConn, 1);
            sourceStep = strings(nConn, 1);
            sourceOutput = strings(nConn, 1);
            selectedFile = strings(nConn, 1);
            targetNodeID = zeros(nConn, 1);
            targetStep = strings(nConn, 1);
            targetInput = strings(nConn, 1);
            sourceTypes = cell(nConn, 1);

            for iConn = 1:nConn
                connLocal = connectionsLocal(iConn);
                sourceNodeID(iConn) = connLocal.sourceNodeID;
                targetNodeID(iConn) = connLocal.targetNodeID;
                sourceOutput(iConn) = string(connLocal.sourceOutputName);
                targetInput(iConn) = string(connLocal.targetInputName);

                if isfield(connLocal, 'selectedFile') && ~isempty(connLocal.selectedFile)
                    selectedFile(iConn) = string(connLocal.selectedFile);
                end

                srcIdx = find([nodesLocal.id] == connLocal.sourceNodeID, 1, 'first');
                if isempty(srcIdx)
                    sourceStep(iConn) = "<missing>";
                else
                    sourceStep(iConn) = string(nodesLocal(srcIdx).name);
                end

                dstIdx = find([nodesLocal.id] == connLocal.targetNodeID, 1, 'first');
                if isempty(dstIdx)
                    targetStep(iConn) = "<missing>";
                else
                    targetStep(iConn) = string(nodesLocal(dstIdx).name);
                end

                if isfield(connLocal, 'sourceOutputType') && ~isempty(connLocal.sourceOutputType)
                    sourceTypes{iConn} = cellstr(string(connLocal.sourceOutputType(:)));
                else
                    sourceTypes{iConn} = {};
                end
            end

            edgeTable = table(sourceNodeID, sourceStep, sourceOutput, selectedFile, ...
                targetNodeID, targetStep, targetInput, sourceTypes);
        end

        function topoOrder = getTopoOrderFromPlan(~, plan)
            %GETTOPOORDERFROMPLAN Extract backend topological order from a plan.

            topoOrder = [];

            if isstruct(plan) && isfield(plan, 'nodeOrders') && isstruct(plan.nodeOrders) && ...
                    isfield(plan.nodeOrders, 'topoOrder') && ~isempty(plan.nodeOrders.topoOrder)
                topoOrder = plan.nodeOrders.topoOrder(:).';
                return
            end

            if isstruct(plan) && isfield(plan, 'nodeTable') && istable(plan.nodeTable) && ...
                    ~isempty(plan.nodeTable) && all(ismember({'nodeID','topoIndex'}, plan.nodeTable.Properties.VariableNames))
                nodeTable = plan.nodeTable;
                [~, ord] = sort(nodeTable.topoIndex, 'ascend');
                topoOrder = nodeTable.nodeID(ord).';
            end
        end

        function row = getNodeTableRowByID(~, nodeTable, nodeID)
            %GETNODETABLEROWBYID Return one nodeTable row by node ID.

            row = table();

            if ~istable(nodeTable) || isempty(nodeTable) || ~ismember('nodeID', nodeTable.Properties.VariableNames)
                return
            end

            idx = find(nodeTable.nodeID == nodeID, 1, 'first');
            if ~isempty(idx)
                row = nodeTable(idx, :);
            end
        end

        function tf = getNodePersistenceFlagFromPlan(~, nodeID, outputPlan)
            %GETNODEPERSISTENCEFLAGFROMPLAN True when a step has planned persisted outputs.

            tf = false;

            if ~istable(outputPlan) || isempty(outputPlan) || ...
                    ~all(ismember({'nodeID','plannedPersistence'}, outputPlan.Properties.VariableNames))
                return
            end

            rows = outputPlan(outputPlan.nodeID == nodeID, :);
            if isempty(rows)
                return
            end

            persistedStates = ["already_file", "saved_explicit", "saved_final", "saved_temp", "saved_leaf"];
            tf = any(ismember(string(rows.plannedPersistence), persistedStates));
        end

        function outLines = getOutputPlanLinesForNode(app, plan, nodeID)
            %GETOUTPUTPLANLINESFORNODE Return output-persistence plan lines for one step.

            outLines = strings(0,1);

            if ~isstruct(plan) || ~isfield(plan, 'outputPlan') || ~istable(plan.outputPlan) || isempty(plan.outputPlan)
                outLines = cellstr(outLines);
                return
            end

            outputPlan = plan.outputPlan;
            requiredVars = {'nodeID','outputName','plannedPersistence','plannedFileName','reason'};
            if ~all(ismember(requiredVars, outputPlan.Properties.VariableNames))
                outLines = cellstr(outLines);
                return
            end

            rows = outputPlan(outputPlan.nodeID == nodeID, :);
            if isempty(rows)
                outLines = cellstr(outLines);
                return
            end

            for iRow = 1:height(rows)
                fileName = string(rows.plannedFileName(iRow));
                persistenceText = app.formatPersistenceStateForDisplay(rows.plannedPersistence(iRow));

                if strlength(strtrim(fileName)) > 0
                    outLines(end+1,1) = "- " + string(rows.outputName(iRow)) + ": " + ...
                        persistenceText + " -> " + fileName; %#ok<AGROW>
                else
                    outLines(end+1,1) = "- " + string(rows.outputName(iRow)) + ": " + ...
                        persistenceText; %#ok<AGROW>
                end

                reasonText = app.formatOutputPlanReasonForDisplay(rows.reason(iRow));
                if strlength(strtrim(reasonText)) > 0
                    outLines(end+1,1) = "  " + reasonText; %#ok<AGROW>
                end
            end

            outLines = cellstr(outLines);
        end

        function lineList = formatDiagnosticReport(~, report)
            %FORMATDIAGNOSTICREPORT Convert a diagnosePipeline report to text lines.

            lineList = strings(0,1);

            if isempty(report) || ~isstruct(report)
                lineList = "<no diagnostic report available>";
                lineList = cellstr(lineList);
                return
            end

            if isfield(report, 'isValid') && report.isValid
                lineList(end+1,1) = "Pipeline status: VALID"; %#ok<AGROW>
            else
                lineList(end+1,1) = "Pipeline status: INVALID"; %#ok<AGROW>
            end

            if isfield(report, 'graphIsValid')
                lineList(end+1,1) = "Graph valid: " + string(logical(report.graphIsValid)); %#ok<AGROW>
            end
            if isfield(report, 'paramsAreValid')
                lineList(end+1,1) = "Parameters valid: " + string(logical(report.paramsAreValid)); %#ok<AGROW>
            end

            errList = {};
            warnList = {};
            if isfield(report, 'errors') && ~isempty(report.errors)
                errList = report.errors;
            end
            if isfield(report, 'warnings') && ~isempty(report.warnings)
                warnList = report.warnings;
            end

            lineList(end+1,1) = ""; %#ok<AGROW>
            lineList(end+1,1) = "Errors (" + string(numel(errList)) + "):"; %#ok<AGROW>
            if isempty(errList)
                lineList(end+1,1) = "  <none>"; %#ok<AGROW>
            else
                for iErr = 1:numel(errList)
                    lineList(end+1,1) = "  " + string(iErr) + ") " + string(errList{iErr}); %#ok<AGROW>
                end
            end

            lineList(end+1,1) = ""; %#ok<AGROW>
            lineList(end+1,1) = "Warnings (" + string(numel(warnList)) + "):"; %#ok<AGROW>
            if isempty(warnList)
                lineList(end+1,1) = "  <none>"; %#ok<AGROW>
            else
                for iWarn = 1:numel(warnList)
                    lineList(end+1,1) = "  " + string(iWarn) + ") " + string(warnList{iWarn}); %#ok<AGROW>
                end
            end

            if isfield(report, 'nodeStatus') && istable(report.nodeStatus) && ~isempty(report.nodeStatus)
                badRows = report.nodeStatus(report.nodeStatus.status ~= "valid", :);
                lineList(end+1,1) = ""; %#ok<AGROW>
                lineList(end+1,1) = "Step status (" + string(height(badRows)) + " non-valid):"; %#ok<AGROW>
                if isempty(badRows)
                    lineList(end+1,1) = "  <all steps valid>"; %#ok<AGROW>
                else
                    for iRow = 1:height(badRows)
                        lineList(end+1,1) = "  " + string(badRows.stepName(iRow)) + ...
                            " | " + string(badRows.status(iRow)) + ...
                            " | " + string(badRows.severity(iRow)); %#ok<AGROW>
                    end
                end
            end

            lineList = cellstr(lineList);
        end

        function removeInvalidStepsInteractive(app)
            %REMOVEINVALIDSTEPSINTERACTIVE Remove invalid graph steps through the backend.

            if isempty(app.pm)
                app.appendDiagnostic('Remove-invalid requested without an attached PipelineManager.');
                app.setStatus('No PipelineManager attached.');
                return
            end

            try
                report = app.pm.diagnosePipeline('verbose', false);
                app.appendValidationReportLogSummary(report, 'Validation before invalid-step cleanup');
            catch ME
                app.appendDiagnostic(sprintf('Could not diagnose pipeline before cleanup: %s', ME.message));
                app.setStatus('Remove invalid failed.');
                return
            end

            invalidRows = table();
            if isfield(report, 'nodeStatus') && istable(report.nodeStatus) && ~isempty(report.nodeStatus)
                invalidRows = report.nodeStatus(ismember(report.nodeStatus.status, ...
                    ["invalid_blocking", "invalid_nonblocking"]), :);
            end

            if isempty(invalidRows)
                app.appendDiagnostic('No invalid steps were found.');
                app.setStatus('No invalid steps found.');
                app.refreshAllViews();
                return
            end

            nBlocking = sum(invalidRows.status == "invalid_blocking");
            nNonBlocking = sum(invalidRows.status == "invalid_nonblocking");

            msg = sprintf(['Remove invalid steps from the current pipeline?\n\n' ...
                'Blocking invalid steps: %d\n' ...
                'Non-blocking invalid steps: %d\n\n' ...
                'This cannot be undone from this dialog.'], ...
                nBlocking, nNonBlocking);

            choice = uiconfirm(app.UIFigure, msg, 'Remove Invalid Steps', ...
                'Options', {'Remove Invalid', 'Cancel'}, ...
                'DefaultOption', 2, ...
                'CancelOption', 2);

            if ~strcmp(choice, 'Remove Invalid')
                app.appendDiagnostic('Remove-invalid action cancelled.');
                app.setStatus('Remove invalid cancelled.');
                return
            end

            try
                removed = app.pm.removeInvalidNodes();

                app.selectedNodeID = NaN;
                app.selectedStepTag = '';
                app.exitReconnectMode();

                reportAfter = app.pm.diagnosePipeline('verbose', false);
                app.appendValidationReportLogSummary(reportAfter, 'Validation after invalid-step cleanup');
                app.refreshAllViews();

                if isempty(removed)
                    app.appendDiagnostic('Invalid-step cleanup completed. No steps were removed.');
                else
                    app.appendDiagnostic(sprintf('Removed invalid step ID(s): %s', strjoin(string(removed), ', ')));
                end

                if isfield(reportAfter, 'isValid') && reportAfter.isValid
                    app.setStatus('Invalid steps removed. Pipeline is valid.');
                else
                    app.setStatus('Invalid-step cleanup completed. Pipeline still needs attention.');
                end

            catch ME
                app.appendDiagnostic(sprintf('Failed to remove invalid steps: %s', ME.message));
                app.setStatus('Failed to remove invalid steps.');
            end
        end

        function txt = formatPersistenceStateForDisplay(~, state)
            %FORMATPERSISTENCESTATEFORDISPLAY Convert backend state to user-facing text.

            state = lower(strtrim(string(state)));

            switch state
                case "already_file"
                    txt = "already file-backed";
                case "saved_explicit"
                    txt = "saved by user setting";
                case {"saved_final", "saved_leaf"}
                    txt = "save final output";
                case "saved_temp"
                    txt = "temporary output for viewer";
                case "transient_explicit"
                    txt = "temporary in memory";
                case "unsaved_reported"
                    txt = "not saved";
                case "unsupported_save"
                    txt = "cannot be saved";
                case "not_final_intermediate"
                    txt = "intermediate output";
                case "not_leaf_intermediate"
                    txt = "intermediate output";
                case "not_produced"
                    txt = "not produced";
                case "missing_expected_file"
                    txt = "missing expected file";
                otherwise
                    txt = string(state);
            end
        end

        function txt = formatOutputPlanReasonForDisplay(~, reasonText)
            %FORMATOUTPUTPLANREASONFORDISPLAY Remove backend jargon from plan reasons.

            txt = string(reasonText);

            if strlength(strtrim(txt)) == 0
                return
            end

            replacements = {
                "LeafOutputPolicy=saveLeaves.", ...
                "Output policy=saveLeaves.";
                "LeafOutputPolicy=viewerTemp and output is a viewer-compatible candidate.", ...
                "Output policy=viewerTemp and output can be opened by DataViewer.";
                "LeafOutputPolicy=viewerTemp but output is not a DataViewer-compatible candidate.", ...
                "Output policy=viewerTemp, but output is not a DataViewer-compatible candidate.";
                "No explicit save target and LeafOutputPolicy=explicit.", ...
                "No explicit save target and output policy=explicit.";
                "Output feeds downstream node(s) or is not a DATA leaf.", ...
                "Output feeds another step or is not DATA."
                };

            for iRep = 1:size(replacements, 1)
                txt = replace(txt, replacements{iRep,1}, replacements{iRep,2});
            end
        end

        function configureExecutionControlDropDowns(app)
            %CONFIGUREEXECUTIONCONTROLDROPDOWNS Configure UI labels and backend values.

            app.AutoSaveFinalOutputsDropDown.Items = { ...
                'Disabled', ...
                'Enabled', ...
                'Enabled as temporary files'};
            app.AutoSaveFinalOutputsDropDown.ItemsData = { ...
                'explicit', ...
                'saveLeaves', ...
                'viewerTemp'};

            app.RAMModeDropDown.Items = { ...
                'Auto', ...
                'RAM-safe'};
            app.RAMModeDropDown.ItemsData = { ...
                'auto', ...
                'ramsafe'};

            app.RAMSafetyPolicyDropDown.Items = { ...
                'Strict', ...
                'Best effort'};
            app.RAMSafetyPolicyDropDown.ItemsData = { ...
                'strict', ...
                'bestEffort'};

            if isempty(app.AutoSaveFinalOutputsDropDown.Value) || ...
                    ~ismember(app.AutoSaveFinalOutputsDropDown.Value, app.AutoSaveFinalOutputsDropDown.ItemsData)
                app.AutoSaveFinalOutputsDropDown.Value = 'saveLeaves';
            end

            if isempty(app.RAMModeDropDown.Value) || ...
                    ~ismember(app.RAMModeDropDown.Value, app.RAMModeDropDown.ItemsData)
                app.RAMModeDropDown.Value = 'auto';
            end

            if isempty(app.RAMSafetyPolicyDropDown.Value) || ...
                    ~ismember(app.RAMSafetyPolicyDropDown.Value, app.RAMSafetyPolicyDropDown.ItemsData)
                app.RAMSafetyPolicyDropDown.Value = 'bestEffort';
            end

            app.setReuseExistingFilesState(app.getReuseExistingFilesUiValue(), false);
            app.updateRamSafetyControlState();
        end

        function syncExecutionControlsFromManager(app)
            %SYNCEXECUTIONCONTROLSFROMMANAGER Reflect PipelineManager settings in the UI.

            app.configureExecutionControlDropDowns();

            if isempty(app.pm)
                return
            end

            try
                if isprop(app.pm, 'leafOutputPolicy') && ...
                        ismember(char(string(app.pm.leafOutputPolicy)), app.AutoSaveFinalOutputsDropDown.ItemsData)
                    app.AutoSaveFinalOutputsDropDown.Value = char(string(app.pm.leafOutputPolicy));
                end
            catch
            end

            try
                if isprop(app.pm, 'ramMode') && ...
                        ismember(char(string(app.pm.ramMode)), app.RAMModeDropDown.ItemsData)
                    app.RAMModeDropDown.Value = char(string(app.pm.ramMode));
                end
            catch
            end

            try
                if isprop(app.pm, 'ramSafePolicy') && ...
                        ismember(char(string(app.pm.ramSafePolicy)), app.RAMSafetyPolicyDropDown.ItemsData)
                    app.RAMSafetyPolicyDropDown.Value = char(string(app.pm.ramSafePolicy));
                end
            catch
            end

            try
                if isprop(app.pm, 'b_skipSteps')
                    app.setReuseExistingFilesState(logical(app.pm.b_skipSteps), false);
                end
            catch
            end

            app.updateRamSafetyControlState();
        end

        function applyExecutionControlsToManager(app)
            %APPLYEXECUTIONCONTROLSTOMANAGER Push visible execution settings to backend.

            if isempty(app.pm)
                return
            end

            if app.bDataViewerMode
                try
                    app.pm.setLeafOutputPolicy('viewerTemp');
                    app.AutoSaveFinalOutputsDropDown.Value = 'viewerTemp';
                catch ME
                    app.appendDiagnostic(sprintf('Failed to enforce DataViewer output policy: %s', ME.message));
                end
            else
                try
                    app.pm.setLeafOutputPolicy(app.AutoSaveFinalOutputsDropDown.Value);
                catch ME
                    app.appendDiagnostic(sprintf('Failed to update output policy: %s', ME.message));
                end
            end

            try
                app.pm.ramMode = char(string(app.RAMModeDropDown.Value));
            catch ME
                app.appendDiagnostic(sprintf('Failed to update RAM mode: %s', ME.message));
            end

            try
                app.pm.ramSafePolicy = char(string(app.RAMSafetyPolicyDropDown.Value));
            catch ME
                app.appendDiagnostic(sprintf('Failed to update RAM safety policy: %s', ME.message));
            end

            try
                if isprop(app.pm, 'b_skipSteps')
                    app.pm.b_skipSteps = app.getReuseExistingFilesUiValue();
                end
            catch ME
                app.appendDiagnostic(sprintf('Failed to update reuse-existing-files setting: %s', ME.message));
            end

            app.updateRamSafetyControlState();
            app.lockDataViewerModeControls();
        end

        function updateRamSafetyControlState(app)
            %UPDATERAMSAFETYCONTROLSTATE Enable safety dropdown only in RAM-safe mode.

            try
                if strcmpi(char(string(app.RAMModeDropDown.Value)), 'ramsafe')
                    app.RAMSafetyPolicyDropDown.Enable = 'on';
                    app.SafetyLabel.Enable = 'on';
                else
                    app.RAMSafetyPolicyDropDown.Enable = 'off';
                    app.SafetyLabel.Enable = 'off';
                end
            catch
            end
        end

        function runPipelineFromGui(app)
            %RUNPIPELINEFROMGUI Validate and execute the current pipeline.

            if app.bPipelineRunning
                app.requestPipelineStop();
                return
            end

            if isempty(app.pm)
                app.appendDiagnostic('Run requested without an attached PipelineManager.');
                app.setStatus('No PipelineManager attached.');
                return
            end

            app.applyExecutionControlsToManager();

            try
                report = app.pm.diagnosePipeline('verbose', false);
                app.appendValidationReportLogSummary(report, 'Pre-run validation');
                app.refreshGraphView();
            catch ME
                app.appendDiagnostic(sprintf('Could not validate pipeline before execution: %s', ME.message));
                app.setStatus('Run aborted.');
                return
            end

            hasBlockingInvalid = false;
            if isfield(report, 'nodeStatus') && istable(report.nodeStatus) && ...
                    ismember('status', report.nodeStatus.Properties.VariableNames)
                hasBlockingInvalid = any(string(report.nodeStatus.status) == "invalid_blocking");
            elseif isfield(report, 'isValid')
                hasBlockingInvalid = ~logical(report.isValid);
            end

            if hasBlockingInvalid
                uialert(app.UIFigure, ...
                    'The pipeline has blocking invalid steps. Fix them or remove invalid steps before running.', ...
                    'Cannot Run Pipeline');
                app.appendDiagnostic('Run aborted because blocking invalid steps were found.');
                app.setStatus('Run blocked by invalid steps.');
                return
            end

            hasNonBlockingInvalid = false;
            if isfield(report, 'hasInvalidNonblocking')
                hasNonBlockingInvalid = logical(report.hasInvalidNonblocking);
            elseif isfield(report, 'nodeStatus') && istable(report.nodeStatus) && ...
                    ismember('status', report.nodeStatus.Properties.VariableNames)
                hasNonBlockingInvalid = any(string(report.nodeStatus.status) == "invalid_nonblocking");
            end

            if hasNonBlockingInvalid
                choice = uiconfirm(app.UIFigure, ...
                    ['The pipeline contains non-blocking invalid steps, usually unused stale steps. ' ...
                    'Execution can continue, but the graph should be cleaned when convenient.'], ...
                    'Run Pipeline', ...
                    'Options', {'Run Anyway', 'Cancel'}, ...
                    'DefaultOption', 1, ...
                    'CancelOption', 2);

                if ~strcmp(choice, 'Run Anyway')
                    app.appendDiagnostic('Run cancelled by user.');
                    app.setStatus('Run cancelled.');
                    return
                end
            end

            cleanupRunState = [];

            try
                app.prepareGuiForPipelineRun();
                cleanupRunState = onCleanup(@() app.restoreGuiAfterPipelineRun()); %#ok<NASGU>

                app.appendDiagnostic('Pipeline execution started.');
                app.setStatus('Running pipeline...');
                drawnow;

                result = app.pm.executePipeline( ...
                    'LeafOutputPolicy', app.getLockedOrSelectedLeafOutputPolicy(), ...
                    'PrintSummary', false, ...
                    'ProgressFcn', @(evt) app.onPipelineProgress(evt), ...
                    'CancelFcn', @() app.isPipelineCancelRequested());

                app.renderExecutionResult(result);
                app.refreshAllViews();

            catch ME
                app.appendDiagnostic(sprintf('Pipeline execution failed: %s', ME.message));
                app.setStatus('Pipeline execution failed.');
            end
        end

        function policy = getLockedOrSelectedLeafOutputPolicy(app)
            %GETLOCKEDORSELECTEDLEAFOUTPUTPOLICY Return active output policy.

            if app.bDataViewerMode
                policy = 'viewerTemp';
            else
                policy = char(string(app.AutoSaveFinalOutputsDropDown.Value));
            end
        end

        function renderExecutionResult(app, result)
            %RENDEREXECUTIONRESULT Store and summarize the latest execution result.
            %
            % Detailed execution tables are shown in the Pipeline Summary window.
            % The Activity Log receives only compact chronological messages.

            app.LastExecutionResult = result;
            app.updateLatestRunSummaryMenuState();

            app.appendExecutionResultLogSummary(result);

            try
                app.showPipelineSummaryInteractive(result);
            catch ME
                app.appendDiagnostic(sprintf('Failed to open pipeline summary window: %s', ME.message));
            end

            if isstruct(result) && isfield(result, 'status')
                app.setStatus(sprintf('Run finished: %s', char(string(result.status))));
            else
                app.setStatus('Run finished.');
            end

            app.notifyDataViewerExecutionFinished(result);
        end

        function T = getExecutionResultTable(~, result, fieldNames)
            %GETEXECUTIONRESULTTABLE Return the first available non-empty result table.

            T = table();

            if ~isstruct(result)
                return
            end

            for iField = 1:numel(fieldNames)
                fieldName = fieldNames{iField};
                if isfield(result, fieldName) && istable(result.(fieldName)) && ~isempty(result.(fieldName))
                    T = result.(fieldName);
                    return
                end
            end
        end

        function value = getTableStringValue(~, T, varName, rowIdx, defaultValue)
            %GETTABLESTRINGVALUE Return a string scalar from a table variable.

            value = string(defaultValue);

            if ~istable(T) || ~ismember(varName, T.Properties.VariableNames) || rowIdx > height(T)
                return
            end

            try
                rawValue = T.(varName)(rowIdx);
                if iscell(rawValue)
                    value = string(rawValue{1});
                else
                    value = string(rawValue);
                end
            catch
                value = string(defaultValue);
            end
        end

        function generateScriptInteractive(app)
            %GENERATESCRIPTINTERACTIVE Export current pipeline as a MATLAB script.

            if isempty(app.pm)
                app.appendDiagnostic('Generate-script requested without an attached PipelineManager.');
                app.setStatus('No PipelineManager attached.');
                return
            end

            [fileName, folderName] = uiputfile({'*.m', 'MATLAB script (*.m)'}, ...
                'Generate pipeline script', fullfile(pwd, 'generated_pipeline.m'));

            if isequal(fileName, 0)
                app.appendDiagnostic('Generate-script action cancelled by user.');
                app.setStatus('Generate script cancelled.');
                return
            end

            [~,~,ext] = fileparts(fileName);
            if isempty(ext)
                fileName = [fileName '.m'];
            end

            scriptFile = fullfile(folderName, fileName);

            try
                app.applyExecutionControlsToManager();
                app.pm.generateScript(scriptFile);
                app.appendDiagnostic(sprintf('Generated pipeline script: %s', scriptFile));
                app.setStatus('Pipeline script generated.');
            catch ME
                app.appendDiagnostic(sprintf('Failed to generate pipeline script: %s', ME.message));
                app.setStatus('Failed to generate script.');
            end
        end

        function exportErrorLogInteractive(app)
            %EXPORTERRORLOGINTERACTIVE Export latest pipeline error log.

            if isempty(app.pm)
                app.appendDiagnostic('Export-error-log requested without an attached PipelineManager.');
                app.setStatus('No PipelineManager attached.');
                return
            end

            startFolder = pwd;
            try
                if isprop(app.pm, 'SaveFolderList') && ~isempty(app.pm.SaveFolderList)
                    startFolder = app.pm.SaveFolderList{1};
                end
            catch
            end

            outputFolder = uigetdir(startFolder, 'Select folder for pipeline error log');
            if isequal(outputFolder, 0)
                app.appendDiagnostic('Export-error-log action cancelled by user.');
                app.setStatus('Export error log cancelled.');
                return
            end

            try
                outFile = app.pm.exportPipeErrorLog(outputFolder);
                app.appendDiagnostic(sprintf('Exported pipeline error log: %s', outFile));
                app.setStatus('Pipeline error log exported.');
            catch ME
                app.appendDiagnostic(sprintf('Failed to export pipeline error log: %s', ME.message));
                app.setStatus('Failed to export error log.');
            end
        end

        function clearPipelineInteractive(app)
            %CLEARPIPELINEINTERACTIVE Clear current pipeline with confirmation.

            if isempty(app.pm)
                app.appendDiagnostic('Clear-pipeline requested without an attached PipelineManager.');
                app.setStatus('No PipelineManager attached.');
                return
            end
            if isempty(app.pm.nodes)
                app.setStatus('Pipeline is empty.');
                return
            end
            choice = uiconfirm(app.UIFigure, ...
                sprintf(['Clear all steps from the current pipeline?\n\n' ...
                'The selected SaveFolder/RawFolder context and backend settings are kept.']), ...
                'Clear Pipeline', ...
                'Options', {'Clear Pipeline', 'Cancel'}, ...
                'DefaultOption', 2, ...
                'CancelOption', 2);

            if ~strcmp(choice, 'Clear Pipeline')
                app.appendDiagnostic('Clear-pipeline action cancelled.');
                app.setStatus('Clear pipeline cancelled.');
                return
            end

            try
                app.pm.clearPipeline();
                app.selectedNodeID = NaN;
                app.selectedStepTag = '';
                app.selectedFunctionName = '';
                app.exitReconnectMode();
                app.refreshAllViews();
                app.appendDiagnostic('Pipeline cleared.');
                app.setStatus('Pipeline cleared.');
            catch ME
                app.appendDiagnostic(sprintf('Failed to clear pipeline: %s', ME.message));
                app.setStatus('Failed to clear pipeline.');
            end
        end

        function tf = selectDataFoldersInteractive(app)
            %SELECTDATAFOLDERSINTERACTIVE Select SaveFolder/RawFolder pairs.
            %
            % This implementation keeps folder-list management local to the
            % PipelineManagerTool. It supports manual one-by-one folder selection
            % and Java-free multi-folder selection from one parent folder.

            tf = false;

            initialTable = app.getFolderPairTableFromPipelineManager();
            if isempty(initialTable)
                initialTable = app.makeEmptyFolderPairTable();
            end
            initialTable = app.validateFolderPairTable(initialTable);

            dlgW = 1120;
            dlgH = 540;

            mainPos = app.UIFigure.Position;
            dlgX = mainPos(1) + max(20, (mainPos(3) - dlgW) / 2);
            dlgY = mainPos(2) + max(20, (mainPos(4) - dlgH) / 2);

            dlg = uifigure( ...
                'Name', 'Select Data Folders', ...
                'Position', [dlgX dlgY dlgW dlgH], ...
                'Resize', 'off', ...
                'WindowStyle', 'modal');

            dlg.UserData = struct('selectedRows', []);

            mainGrid = uigridlayout(dlg);
            mainGrid.RowHeight = {36, '1x', 42};
            mainGrid.ColumnWidth = {'1x'};
            mainGrid.Padding = [10 10 10 10];
            mainGrid.RowSpacing = 8;

            headerLabel = uilabel(mainGrid);
            headerLabel.Layout.Row = 1;
            headerLabel.Layout.Column = 1;
            headerLabel.Text = ['Build the SaveFolder / RawFolder list used to create the pipeline context. ' ...
                'RawFolder may be Missing for pipelines that do not need raw data.'];
            headerLabel.FontWeight = 'bold';

            tbl = uitable(mainGrid);
            tbl.Layout.Row = 2;
            tbl.Layout.Column = 1;
            tbl.Data = initialTable;
            tbl.ColumnName = {'Use', 'SaveFolder', 'RawFolder', 'Status', 'Notes'};
            tbl.ColumnEditable = [true false false false true];
            tbl.ColumnFormat = {'logical', 'char', 'char', 'char', 'char'};
            tbl.CellSelectionCallback = @onTableSelection;
            tbl.CellEditCallback = @onTableEdited;

            bottomGrid = uigridlayout(mainGrid);
            bottomGrid.Layout.Row = 3;
            bottomGrid.Layout.Column = 1;
            bottomGrid.RowHeight = {'1x'};
            bottomGrid.ColumnWidth = {72, 108, 148, 102, 142, 82, 70, '1x', 76, 76};
            bottomGrid.ColumnSpacing = 6;
            bottomGrid.Padding = [0 0 0 0];

            addPairButton = uibutton(bottomGrid, 'push');
            addPairButton.Text = 'Add Pair';
            addPairButton.Tooltip = 'Add one SaveFolder/RawFolder pair.';
            addPairButton.Layout.Row = 1;
            addPairButton.Layout.Column = 1;
            addPairButton.ButtonPushedFcn = @(~,~) addPair();

            addSaveButton = uibutton(bottomGrid, 'push');
            addSaveButton.Text = 'Add SaveFolder';
            addSaveButton.Tooltip = 'Add one SaveFolder and mark RawFolder as Missing.';
            addSaveButton.Layout.Row = 1;
            addSaveButton.Layout.Column = 2;
            addSaveButton.ButtonPushedFcn = @(~,~) addSaveFolderOnly();

            addSaveParentButton = uibutton(bottomGrid, 'push');
            addSaveParentButton.Text = 'Save From Parent...';
            addSaveParentButton.Tooltip = 'Scan one parent folder and add selected folders as SaveFolders.';
            addSaveParentButton.Layout.Row = 1;
            addSaveParentButton.Layout.Column = 3;
            addSaveParentButton.ButtonPushedFcn = @(~,~) addSaveFoldersFromParent();

            setRawButton = uibutton(bottomGrid, 'push');
            setRawButton.Text = 'Set RawFolder';
            setRawButton.Tooltip = 'Set one RawFolder for the selected row(s).';
            setRawButton.Layout.Row = 1;
            setRawButton.Layout.Column = 4;
            setRawButton.ButtonPushedFcn = @(~,~) setRawFolderForSelection();

            addRawParentButton = uibutton(bottomGrid, 'push');
            addRawParentButton.Text = 'Raw From Parent...';
            addRawParentButton.Tooltip = 'Scan one parent folder and assign selected folders as RawFolders to selected rows.';
            addRawParentButton.Layout.Row = 1;
            addRawParentButton.Layout.Column = 5;
            addRawParentButton.ButtonPushedFcn = @(~,~) addRawFoldersFromParent();

            rawEqualsSaveButton = uibutton(bottomGrid, 'push');
            rawEqualsSaveButton.Text = 'Raw = Save';
            rawEqualsSaveButton.Tooltip = 'Set RawFolder equal to SaveFolder for the selected row(s).';
            rawEqualsSaveButton.Layout.Row = 1;
            rawEqualsSaveButton.Layout.Column = 6;
            rawEqualsSaveButton.ButtonPushedFcn = @(~,~) setRawEqualsSave();

            removeButton = uibutton(bottomGrid, 'push');
            removeButton.Text = 'Remove';
            removeButton.Tooltip = 'Remove selected row(s).';
            removeButton.Layout.Row = 1;
            removeButton.Layout.Column = 7;
            removeButton.ButtonPushedFcn = @(~,~) removeSelectedRows();

            applyButton = uibutton(bottomGrid, 'push');
            applyButton.Text = 'Apply';
            applyButton.FontWeight = 'bold';
            applyButton.Layout.Row = 1;
            applyButton.Layout.Column = 9;
            applyButton.ButtonPushedFcn = @(~,~) applySelection();

            cancelButton = uibutton(bottomGrid, 'push');
            cancelButton.Text = 'Cancel';
            cancelButton.Layout.Row = 1;
            cancelButton.Layout.Column = 10;
            cancelButton.ButtonPushedFcn = @(~,~) cancelSelection();

            dlg.CloseRequestFcn = @(~,~) cancelSelection();

            uiwait(dlg);

            if isvalid(dlg)
                delete(dlg);
            end

            % -------------------------------------------------------------
            % Nested callbacks
            % -------------------------------------------------------------

            function onTableSelection(~, event)
                if isempty(event.Indices)
                    dlg.UserData = struct('selectedRows', []);
                    return
                end

                selectedRows = unique(event.Indices(:,1), 'stable');
                dlg.UserData = struct('selectedRows', selectedRows(:).');
            end

            function onTableEdited(~, ~)
                T = tbl.Data;
                tbl.Data = app.validateFolderPairTable(T);
            end

            function rows = getSelectedRows()
                rows = [];
                if isempty(dlg.UserData) || ~isstruct(dlg.UserData) || ~isfield(dlg.UserData, 'selectedRows')
                    return
                end
                rows = dlg.UserData.selectedRows;
                rows = rows(rows >= 1 & rows <= height(tbl.Data));
            end

            function addPair()
                saveFolder = uigetdir(app.getDefaultFolderSelectionStartFolder('save'), 'Select SaveFolder');
                if isequal(saveFolder, 0)
                    return
                end

                rawFolder = uigetdir(saveFolder, 'Select RawFolder - Cancel to mark as Missing');
                if isequal(rawFolder, 0)
                    rawFolder = 'Missing';
                end

                T = tbl.Data;
                newRow = table(true, string(saveFolder), string(rawFolder), "", "", ...
                    'VariableNames', {'Use','SaveFolder','RawFolder','Status','Notes'});
                T = [T; newRow]; %#ok<AGROW>
                tbl.Data = app.validateFolderPairTable(T);
            end

            function addSaveFolderOnly()
                saveFolder = uigetdir(app.getDefaultFolderSelectionStartFolder('save'), 'Select SaveFolder');
                if isequal(saveFolder, 0)
                    return
                end

                T = tbl.Data;
                newRow = table(true, string(saveFolder), "Missing", "", "", ...
                    'VariableNames', {'Use','SaveFolder','RawFolder','Status','Notes'});
                T = [T; newRow]; %#ok<AGROW>
                tbl.Data = app.validateFolderPairTable(T);
            end

            function addSaveFoldersFromParent()
                selectedFolders = app.selectFoldersFromParentInteractive( ...
                    dlg, ...
                    app.getDefaultFolderSelectionStartFolder('save'), ...
                    'Add SaveFolders From Parent');

                if isempty(selectedFolders)
                    return
                end

                selectedFolders = string(selectedFolders(:));
                T = tbl.Data;

                for iFolder = 1:numel(selectedFolders)
                    newRow = table(true, selectedFolders(iFolder), "Missing", "", "", ...
                        'VariableNames', {'Use','SaveFolder','RawFolder','Status','Notes'});
                    T = [T; newRow]; %#ok<AGROW>
                end

                tbl.Data = app.validateFolderPairTable(T);
            end

            function setRawFolderForSelection()
                rows = getSelectedRows();
                if isempty(rows)
                    uialert(dlg, 'Select one or more table rows first.', 'No Row Selected');
                    return
                end

                rawFolder = uigetdir(app.getDefaultFolderSelectionStartFolder('raw'), 'Select RawFolder');
                if isequal(rawFolder, 0)
                    return
                end

                T = tbl.Data;
                for iRow = rows
                    T.RawFolder(iRow) = string(rawFolder);
                end
                tbl.Data = app.validateFolderPairTable(T);
            end

            function addRawFoldersFromParent()
                rows = getSelectedRows();

                if isempty(rows)
                    uialert(dlg, ...
                        ['Select the SaveFolder row(s) that should receive RawFolder assignments ' ...
                        'before using Raw From Parent.'], ...
                        'No Row Selected');
                    return
                end

                selectedFolders = app.selectFoldersFromParentInteractive( ...
                    dlg, ...
                    app.getDefaultFolderSelectionStartFolder('raw'), ...
                    'Add RawFolders From Parent');

                if isempty(selectedFolders)
                    return
                end

                selectedFolders = string(selectedFolders(:));

                if numel(selectedFolders) ~= numel(rows)
                    uialert(dlg, sprintf(['The number of selected RawFolders must match the number of selected ' ...
                        'SaveFolder rows.\n\nSelected SaveFolder rows: %d\nSelected RawFolders: %d'], ...
                        numel(rows), numel(selectedFolders)), ...
                        'RawFolder Count Mismatch');
                    return
                end

                T = tbl.Data;
                rows = rows(:);
                rows = rows(rows >= 1 & rows <= height(T));

                for iRow = 1:numel(rows)
                    T.RawFolder(rows(iRow)) = selectedFolders(iRow);
                end

                tbl.Data = app.validateFolderPairTable(T);
            end

            function setRawEqualsSave()
                rows = getSelectedRows();
                if isempty(rows)
                    uialert(dlg, 'Select one or more table rows first.', 'No Row Selected');
                    return
                end

                T = tbl.Data;
                for iRow = rows
                    T.RawFolder(iRow) = T.SaveFolder(iRow);
                end
                tbl.Data = app.validateFolderPairTable(T);
            end

            function removeSelectedRows()
                rows = getSelectedRows();
                if isempty(rows)
                    uialert(dlg, 'Select one or more table rows first.', 'No Row Selected');
                    return
                end

                T = tbl.Data;
                T(rows, :) = [];
                dlg.UserData = struct('selectedRows', []);
                tbl.Data = app.validateFolderPairTable(T);
            end

            function applySelection()
                T = app.validateFolderPairTable(tbl.Data);

                if isempty(T) || height(T) == 0
                    uialert(dlg, 'Add at least one SaveFolder/RawFolder pair before applying.', ...
                        'No Folder Pairs');
                    return
                end

                enabledRows = T(T.Use, :);

                if isempty(enabledRows)
                    uialert(dlg, 'Enable at least one folder pair before applying.', ...
                        'No Enabled Folder Pairs');
                    return
                end

                blockingRows = enabledRows(contains(enabledRows.Status, "Missing SaveFolder") | ...
                    contains(enabledRows.Status, "Duplicate SaveFolder"), :);

                if ~isempty(blockingRows)
                    uialert(dlg, ['One or more enabled rows cannot be applied.' newline newline ...
                        'Fix missing or duplicate SaveFolder entries before continuing.'], ...
                        'Invalid Folder Pairs');
                    tbl.Data = T;
                    return
                end

                if ~isempty(app.pm) && isprop(app.pm, 'nodes') && ~isempty(app.pm.nodes)
                    choice = uiconfirm(dlg, ...
                        ['Applying a new folder selection will create a new empty PipelineManager context.' newline newline ...
                        'The current pipeline steps will be replaced in this first implementation.'], ...
                        'Replace Current Pipeline Context', ...
                        'Options', {'Replace Context', 'Cancel'}, ...
                        'DefaultOption', 2, ...
                        'CancelOption', 2);

                    if ~strcmp(choice, 'Replace Context')
                        return
                    end
                end

                try
                    app.applyFolderPairTable(T);
                    tf = true;
                    uiresume(dlg);

                catch ME
                    uialert(dlg, ME.message, 'Failed to Apply Folder Pairs');
                end
            end

            function cancelSelection()
                tf = false;
                uiresume(dlg);
            end
        end

        function T = makeEmptyFolderPairTable(~)
            %MAKEEMPTYFOLDERPAIRTABLE Return an empty folder-pair table.

            T = table( ...
                false(0,1), ...
                strings(0,1), ...
                strings(0,1), ...
                strings(0,1), ...
                strings(0,1), ...
                'VariableNames', {'Use','SaveFolder','RawFolder','Status','Notes'});
        end

        function T = getFolderPairTableFromPipelineManager(app)
            %GETFOLDERPAIRTABLEFROMPIPELINEMANAGER Build a folder-pair table from app.pm.

            T = app.makeEmptyFolderPairTable();

            if isempty(app.pm)
                return
            end

            if ~isprop(app.pm, 'SaveFolderList') || isempty(app.pm.SaveFolderList)
                return
            end

            saveFolders = string(app.pm.SaveFolderList(:));

            if isprop(app.pm, 'RawFolderList') && ~isempty(app.pm.RawFolderList)
                rawFolders = string(app.pm.RawFolderList(:));
            else
                rawFolders = repmat("Missing", numel(saveFolders), 1);
            end

            if numel(rawFolders) ~= numel(saveFolders)
                rawFolders = repmat("Missing", numel(saveFolders), 1);
            end

            T = table( ...
                true(numel(saveFolders),1), ...
                saveFolders(:), ...
                rawFolders(:), ...
                strings(numel(saveFolders),1), ...
                strings(numel(saveFolders),1), ...
                'VariableNames', {'Use','SaveFolder','RawFolder','Status','Notes'});

            T = app.validateFolderPairTable(T);
        end

        function T = validateFolderPairTable(~, T)
            %VALIDATEFOLDERPAIRTABLE Validate folder-pair rows for UI display.

            if isempty(T)
                T = table( ...
                    false(0,1), ...
                    strings(0,1), ...
                    strings(0,1), ...
                    strings(0,1), ...
                    strings(0,1), ...
                    'VariableNames', {'Use','SaveFolder','RawFolder','Status','Notes'});
                return
            end

            requiredVars = {'Use','SaveFolder','RawFolder','Status','Notes'};
            for iVar = 1:numel(requiredVars)
                if ~ismember(requiredVars{iVar}, T.Properties.VariableNames)
                    switch requiredVars{iVar}
                        case 'Use'
                            T.Use = true(height(T),1);
                        otherwise
                            T.(requiredVars{iVar}) = strings(height(T),1);
                    end
                end
            end

            T.Use = logical(T.Use);
            T.SaveFolder = string(T.SaveFolder);
            T.RawFolder = string(T.RawFolder);
            T.Status = strings(height(T),1);
            T.Notes = string(T.Notes);

            emptyRaw = strlength(strtrim(T.RawFolder)) == 0;
            T.RawFolder(emptyRaw) = "Missing";

            enabledRows = find(T.Use);
            enabledSaveFolders = strtrim(T.SaveFolder(enabledRows));
            duplicateEnabledMask = false(height(T),1);

            for i = 1:numel(enabledRows)
                rowIdx = enabledRows(i);
                sf = enabledSaveFolders(i);
                if strlength(sf) == 0
                    continue
                end
                isDup = sum(strcmpi(enabledSaveFolders, sf)) > 1;
                duplicateEnabledMask(rowIdx) = isDup;
            end

            for iRow = 1:height(T)
                rowIssues = strings(0,1);

                saveFolder = strtrim(T.SaveFolder(iRow));
                rawFolder = strtrim(T.RawFolder(iRow));

                if ~T.Use(iRow)
                    T.Status(iRow) = "Disabled";
                    continue
                end

                if strlength(saveFolder) == 0 || ~isfolder(char(saveFolder))
                    rowIssues(end+1,1) = "Missing SaveFolder"; %#ok<AGROW>
                end

                if duplicateEnabledMask(iRow)
                    rowIssues(end+1,1) = "Duplicate SaveFolder"; %#ok<AGROW>
                end

                if strlength(rawFolder) == 0
                    T.RawFolder(iRow) = "Missing";
                    rawFolder = "Missing";
                end

                if ~strcmpi(rawFolder, "Missing") && ~isfolder(char(rawFolder))
                    rowIssues(end+1,1) = "Missing RawFolder"; %#ok<AGROW>
                end

                if isempty(rowIssues)
                    if strcmpi(rawFolder, "Missing")
                        T.Status(iRow) = "Ready - RawFolder Missing";
                    else
                        T.Status(iRow) = "Ready";
                    end
                else
                    T.Status(iRow) = strjoin(rowIssues, "; ");
                end
            end
        end

        function applyFolderPairTable(app, T)
            %APPLYFOLDERPAIRTABLE Create a new PipelineManager from enabled folder pairs.

            T = app.validateFolderPairTable(T);
            enabledRows = T(T.Use, :);

            if isempty(enabledRows)
                error('PipelineManagerTool:applyFolderPairTable:NoEnabledRows', ...
                    'No enabled folder pairs were provided.');
            end

            blockingRows = enabledRows(contains(enabledRows.Status, "Missing SaveFolder") | ...
                contains(enabledRows.Status, "Duplicate SaveFolder"), :);

            if ~isempty(blockingRows)
                error('PipelineManagerTool:applyFolderPairTable:InvalidRows', ...
                    'One or more enabled rows have invalid SaveFolder entries.');
            end

            saveFolders = cellstr(enabledRows.SaveFolder);
            rawFolders = cellstr(enabledRows.RawFolder);

            rawFolders(cellfun(@isempty, rawFolders)) = {'Missing'};

            % Create the backend object first.
            app.pm = PipelineManager(saveFolders, rawFolders);
            app.FolderPairTable = T;

            % Critical: the visible top-panel settings are the source of truth
            % for a newly created backend context. Do this immediately after
            % object creation, before graph/output planning is refreshed.
            app.configureExecutionControlDropDowns();
            app.applyExecutionControlsToManager();

            app.selectedNodeID = NaN;
            app.selectedStepTag = '';
            app.selectedFunctionName = '';

            app.bReconnectMode = false;
            app.reconnectTargetNodeID = NaN;
            app.reconnectTargetStepTag = '';
            app.reconnectTargetInputName = '';
            app.reconnectCandidateNodeIDs = [];

            app.setPipelineControlsAvailable(true);
            app.updateFolderSelectionButtonSummary();

            app.populateFunctionTree();
            app.refreshAllViews();

            nEnabled = height(enabledRows);
            app.appendDiagnostic(sprintf(['Selected %d folder pair(s) and created a new PipelineManager context. ' ...
                'Execution settings applied from the top panel.'], nEnabled));
            app.setStatus(sprintf('Data folders selected: %d pair(s).', nEnabled));
        end

        function updateFolderSelectionButtonSummary(app)
            %UPDATEFOLDERSELECTIONBUTTONSUMMARY Update Select Data Folders text/tooltip.

            nPairs = 0;
            firstSave = "";
            firstRaw = "";

            if ~isempty(app.FolderPairTable) && istable(app.FolderPairTable) && height(app.FolderPairTable) > 0
                enabledRows = app.FolderPairTable(app.FolderPairTable.Use, :);
                nPairs = height(enabledRows);
                if nPairs > 0
                    firstSave = string(enabledRows.SaveFolder(1));
                    firstRaw = string(enabledRows.RawFolder(1));
                end
            end

            if nPairs == 0
                app.SelectDataFoldersButton.Text = 'Select Data Folders';
                app.SelectDataFoldersButton.Tooltip = { ...
                    'Choose or edit the SaveFolder/RawFolder pairs used to create the pipeline context.'};
            elseif nPairs == 1
                app.SelectDataFoldersButton.Text = 'Data Folders (1)';
                app.SelectDataFoldersButton.Tooltip = {sprintf('1 folder pair selected.\nSaveFolder: %s\nRawFolder: %s', ...
                    char(firstSave), char(firstRaw))};
            else
                app.SelectDataFoldersButton.Text = sprintf('Data Folders (%d)', nPairs);
                app.SelectDataFoldersButton.Tooltip = {sprintf('%d folder pairs selected.\nFirst SaveFolder: %s\nFirst RawFolder: %s', ...
                    nPairs, char(firstSave), char(firstRaw))};
            end
        end

        function selectedFolders = selectFoldersFromParentInteractive(app, parentFig, startFolder, titleText)
            %SELECTFOLDERSFROMPARENTINTERACTIVE Select multiple folders without Java.
            %
            % The user selects one parent folder, optionally scans recursively, then
            % checks one or more candidate folders from a table. A filter field uses
            % case-insensitive contains matching to pre-check and highlight rows.

            selectedFolders = {};

            if nargin < 2 || isempty(parentFig) || ~isvalid(parentFig)
                parentFig = app.UIFigure;
            end

            if nargin < 3 || isempty(startFolder) || ~isfolder(char(string(startFolder)))
                startFolder = pwd;
            end

            if nargin < 4 || isempty(titleText)
                titleText = 'Select Folders From Parent';
            end

            parentFolder = uigetdir(char(string(startFolder)), titleText);
            if isequal(parentFolder, 0)
                return
            end

            dlgW = 980;
            dlgH = 560;

            try
                parentPos = parentFig.Position;
            catch
                parentPos = app.UIFigure.Position;
            end

            dlgX = parentPos(1) + max(20, (parentPos(3) - dlgW) / 2);
            dlgY = parentPos(2) + max(20, (parentPos(4) - dlgH) / 2);

            dlg = uifigure( ...
                'Name', titleText, ...
                'Position', [dlgX dlgY dlgW dlgH], ...
                'Resize', 'off', ...
                'WindowStyle', 'modal');

            dlg.UserData = struct( ...
                'parentFolder', string(parentFolder), ...
                'matchedRows', []);

            mainGrid = uigridlayout(dlg);
            mainGrid.RowHeight = {36, 36, '1x', 42};
            mainGrid.ColumnWidth = {'1x'};
            mainGrid.Padding = [10 10 10 10];
            mainGrid.RowSpacing = 8;

            parentGrid = uigridlayout(mainGrid);
            parentGrid.Layout.Row = 1;
            parentGrid.Layout.Column = 1;
            parentGrid.RowHeight = {'1x'};
            parentGrid.ColumnWidth = {90, '1x', 110, 95};
            parentGrid.Padding = [0 0 0 0];
            parentGrid.ColumnSpacing = 6;

            parentLabel = uilabel(parentGrid);
            parentLabel.Text = 'Parent folder:';
            parentLabel.Layout.Row = 1;
            parentLabel.Layout.Column = 1;

            parentField = uieditfield(parentGrid, 'text');
            parentField.Value = char(string(parentFolder));
            parentField.Editable = 'off';
            parentField.Layout.Row = 1;
            parentField.Layout.Column = 2;

            browseButton = uibutton(parentGrid, 'push');
            browseButton.Text = 'Browse...';
            browseButton.Layout.Row = 1;
            browseButton.Layout.Column = 3;
            browseButton.ButtonPushedFcn = @(~,~) browseParentFolder();

            recursiveCheckBox = uicheckbox(parentGrid);
            recursiveCheckBox.Text = 'Recursive';
            recursiveCheckBox.Value = false;
            recursiveCheckBox.Layout.Row = 1;
            recursiveCheckBox.Layout.Column = 4;
            recursiveCheckBox.ValueChangedFcn = @(~,~) refreshCandidates();

            filterGrid = uigridlayout(mainGrid);
            filterGrid.Layout.Row = 2;
            filterGrid.Layout.Column = 1;
            filterGrid.RowHeight = {'1x'};
            filterGrid.ColumnWidth = {60, '1x', 56, 90, 80};
            filterGrid.Padding = [0 0 0 0];
            filterGrid.ColumnSpacing = 6;

            filterLabel = uilabel(filterGrid);
            filterLabel.Text = 'Filter:';
            filterLabel.Layout.Row = 1;
            filterLabel.Layout.Column = 1;

            filterField = uieditfield(filterGrid, 'text');
            filterField.Placeholder = 'Text contained in folder name or path';
            filterField.Layout.Row = 1;
            filterField.Layout.Column = 2;
            filterField.ValueChangedFcn = @(~,~) previewFilterMatches();

            filterOkButton = uibutton(filterGrid, 'push');
            filterOkButton.Text = 'OK';
            filterOkButton.Tooltip = 'Apply filter selection and remove highlight.';
            filterOkButton.Layout.Row = 1;
            filterOkButton.Layout.Column = 3;
            filterOkButton.ButtonPushedFcn = @(~,~) applyFilterAndClearHighlight();

            selectAllButton = uibutton(filterGrid, 'push');
            selectAllButton.Text = 'Select All';
            selectAllButton.Layout.Row = 1;
            selectAllButton.Layout.Column = 4;
            selectAllButton.ButtonPushedFcn = @(~,~) selectAllCandidates();

            clearButton = uibutton(filterGrid, 'push');
            clearButton.Text = 'Clear';
            clearButton.Layout.Row = 1;
            clearButton.Layout.Column = 5;
            clearButton.ButtonPushedFcn = @(~,~) clearCandidateSelection();

            candTable = uitable(mainGrid);
            candTable.Layout.Row = 3;
            candTable.Layout.Column = 1;
            candTable.ColumnName = {'Use', 'FolderName', 'FullPath', 'Status'};
            candTable.ColumnEditable = [true false false false];
            candTable.ColumnFormat = {'logical', 'char', 'char', 'char'};

            bottomGrid = uigridlayout(mainGrid);
            bottomGrid.Layout.Row = 4;
            bottomGrid.Layout.Column = 1;
            bottomGrid.RowHeight = {'1x'};
            bottomGrid.ColumnWidth = {'1x', 80, 80};
            bottomGrid.Padding = [0 0 0 0];
            bottomGrid.ColumnSpacing = 8;

            infoLabel = uilabel(bottomGrid);
            infoLabel.Layout.Row = 1;
            infoLabel.Layout.Column = 1;
            infoLabel.Text = '';
            infoLabel.FontColor = [0.35 0.35 0.35];

            applyButton = uibutton(bottomGrid, 'push');
            applyButton.Text = 'Apply';
            applyButton.FontWeight = 'bold';
            applyButton.Layout.Row = 1;
            applyButton.Layout.Column = 2;
            applyButton.ButtonPushedFcn = @(~,~) applyCandidateSelection();

            cancelButton = uibutton(bottomGrid, 'push');
            cancelButton.Text = 'Cancel';
            cancelButton.Layout.Row = 1;
            cancelButton.Layout.Column = 3;
            cancelButton.ButtonPushedFcn = @(~,~) cancelCandidateSelection();

            dlg.CloseRequestFcn = @(~,~) cancelCandidateSelection();

            refreshCandidates();
            uiwait(dlg);

            if isvalid(dlg)
                delete(dlg);
            end

            % -------------------------------------------------------------
            % Nested callbacks
            % -------------------------------------------------------------

            function browseParentFolder()
                newParent = uigetdir(char(string(dlg.UserData.parentFolder)), titleText);
                if isequal(newParent, 0)
                    return
                end

                data = dlg.UserData;
                data.parentFolder = string(newParent);
                data.matchedRows = [];
                dlg.UserData = data;

                parentField.Value = char(string(newParent));
                filterField.Value = '';
                refreshCandidates();
            end

            function refreshCandidates()
                data = dlg.UserData;
                parentFolderLocal = char(string(data.parentFolder));

                T = app.scanFolderCandidates(parentFolderLocal, recursiveCheckBox.Value);
                candTable.Data = T;
                filterField.Value = '';
                data.matchedRows = [];
                dlg.UserData = data;

                app.clearTableStylesSafe(candTable);
                updateInfoLabel();
            end

            function previewFilterMatches()
                T = candTable.Data;
                if isempty(T)
                    return
                end

                app.clearTableStylesSafe(candTable);

                filterText = strtrim(string(filterField.Value));
                if strlength(filterText) == 0
                    data = dlg.UserData;
                    data.matchedRows = [];
                    dlg.UserData = data;
                    updateInfoLabel();
                    return
                end

                folderName = lower(string(T.FolderName));
                fullPath = lower(string(T.FullPath));
                token = lower(filterText);

                matchMask = contains(folderName, token) | contains(fullPath, token);

                T.Use(:) = false;
                T.Use(matchMask) = true;
                candTable.Data = T;

                matchedRows = find(matchMask);
                data = dlg.UserData;
                data.matchedRows = matchedRows(:).';
                dlg.UserData = data;

                app.highlightTableRowsSafe(candTable, matchedRows);
                updateInfoLabel();
            end

            function applyFilterAndClearHighlight()
                previewFilterMatches();
                app.clearTableStylesSafe(candTable);

                data = dlg.UserData;
                data.matchedRows = [];
                dlg.UserData = data;

                updateInfoLabel();
            end

            function selectAllCandidates()
                T = candTable.Data;
                if isempty(T)
                    return
                end

                T.Use(:) = true;
                candTable.Data = T;
                app.clearTableStylesSafe(candTable);

                data = dlg.UserData;
                data.matchedRows = [];
                dlg.UserData = data;

                updateInfoLabel();
            end

            function clearCandidateSelection()
                T = candTable.Data;
                if isempty(T)
                    return
                end

                T.Use(:) = false;
                candTable.Data = T;
                filterField.Value = '';
                app.clearTableStylesSafe(candTable);

                data = dlg.UserData;
                data.matchedRows = [];
                dlg.UserData = data;

                updateInfoLabel();
            end

            function applyCandidateSelection()
                T = candTable.Data;

                if isempty(T) || height(T) == 0
                    selectedFolders = {};
                    uiresume(dlg);
                    return
                end

                selectedRows = T(T.Use, :);
                if isempty(selectedRows)
                    uialert(dlg, 'Select at least one folder before applying.', 'No Folders Selected');
                    return
                end

                selectedFolders = cellstr(string(selectedRows.FullPath));
                app.clearTableStylesSafe(candTable);
                uiresume(dlg);
            end

            function cancelCandidateSelection()
                selectedFolders = {};
                uiresume(dlg);
            end

            function updateInfoLabel()
                T = candTable.Data;
                if isempty(T)
                    infoLabel.Text = 'No candidate folders found.';
                    return
                end

                nSelected = sum(T.Use);
                data = dlg.UserData;
                nMatched = 0;
                if isstruct(data) && isfield(data, 'matchedRows')
                    nMatched = numel(data.matchedRows);
                end

                if nMatched > 0
                    infoLabel.Text = sprintf('%d candidate folder(s), %d matched by filter, %d selected.', ...
                        height(T), nMatched, nSelected);
                else
                    infoLabel.Text = sprintf('%d candidate folder(s), %d selected.', ...
                        height(T), nSelected);
                end
            end
        end

        function T = scanFolderCandidates(~, parentFolder, recursiveFlag)
            %SCANFOLDERCANDIDATES Return candidate subfolders under one parent.
            %
            % Output table columns:
            %   Use, FolderName, FullPath, Status

            if nargin < 3 || isempty(recursiveFlag)
                recursiveFlag = false;
            end

            parentFolder = char(string(parentFolder));

            if isempty(parentFolder) || ~isfolder(parentFolder)
                T = table( ...
                    false(0,1), ...
                    strings(0,1), ...
                    strings(0,1), ...
                    strings(0,1), ...
                    'VariableNames', {'Use','FolderName','FullPath','Status'});
                return
            end

            folderPaths = strings(0,1);

            if recursiveFlag
                pending = string(parentFolder);

                while ~isempty(pending)
                    currentFolder = pending(1);
                    pending(1) = [];
                    pending = pending(:);

                    d = dir(char(currentFolder));
                    d = d([d.isdir]);

                    if isempty(d)
                        continue
                    end

                    names = string({d.name})';
                    keep = names ~= "." & names ~= "..";
                    names = names(keep);

                    if isempty(names)
                        continue
                    end

                    fullPaths = strings(numel(names), 1);
                    for iName = 1:numel(names)
                        fullPaths(iName) = string(fullfile(char(currentFolder), char(names(iName))));
                    end

                    folderPaths = [folderPaths; fullPaths]; %#ok<AGROW>
                    pending = [pending; fullPaths]; %#ok<AGROW>
                end

            else
                d = dir(parentFolder);
                d = d([d.isdir]);

                if ~isempty(d)
                    names = string({d.name})';
                    keep = names ~= "." & names ~= "..";
                    names = names(keep);

                    folderPaths = strings(numel(names), 1);
                    for iName = 1:numel(names)
                        folderPaths(iName) = string(fullfile(parentFolder, char(names(iName))));
                    end
                end
            end

            if isempty(folderPaths)
                T = table( ...
                    false(0,1), ...
                    strings(0,1), ...
                    strings(0,1), ...
                    strings(0,1), ...
                    'VariableNames', {'Use','FolderName','FullPath','Status'});
                return
            end

            folderPaths = unique(folderPaths(:), 'stable');
            [~, ord] = sort(lower(folderPaths));
            folderPaths = folderPaths(ord);

            folderNames = strings(numel(folderPaths), 1);
            status = strings(numel(folderPaths), 1);

            for iFolder = 1:numel(folderPaths)
                [~, folderNames(iFolder)] = fileparts(folderPaths(iFolder));

                if isfolder(char(folderPaths(iFolder)))
                    status(iFolder) = "Ready";
                else
                    status(iFolder) = "Missing";
                end
            end

            T = table( ...
                false(numel(folderPaths), 1), ...
                folderNames(:), ...
                folderPaths(:), ...
                status(:), ...
                'VariableNames', {'Use','FolderName','FullPath','Status'});
        end

        function startFolder = getDefaultFolderSelectionStartFolder(app, mode)
            %GETDEFAULTFOLDERSELECTIONSTARTFOLDER Return a useful folder-picker start path.

            if nargin < 2 || isempty(mode)
                mode = 'save';
            end

            startFolder = pwd;

            try
                if ~isempty(app.FolderPairTable) && istable(app.FolderPairTable) && height(app.FolderPairTable) > 0
                    enabledRows = app.FolderPairTable(app.FolderPairTable.Use, :);

                    if isempty(enabledRows)
                        enabledRows = app.FolderPairTable;
                    end

                    switch lower(char(string(mode)))
                        case 'raw'
                            rawFolders = string(enabledRows.RawFolder);
                            rawFolders = rawFolders(strlength(strtrim(rawFolders)) > 0 & ~strcmpi(rawFolders, "Missing"));
                            if ~isempty(rawFolders) && isfolder(char(rawFolders(1)))
                                startFolder = char(rawFolders(1));
                                return
                            end

                        otherwise
                            saveFolders = string(enabledRows.SaveFolder);
                            saveFolders = saveFolders(strlength(strtrim(saveFolders)) > 0);
                            if ~isempty(saveFolders) && isfolder(char(saveFolders(1)))
                                startFolder = char(saveFolders(1));
                                return
                            end
                    end
                end
            catch
            end

            try
                if ~isempty(app.pm) && isprop(app.pm, 'SaveFolderList') && ~isempty(app.pm.SaveFolderList) && ...
                        isfolder(app.pm.SaveFolderList{1})
                    startFolder = app.pm.SaveFolderList{1};
                end
            catch
            end
        end

        function highlightTableRowsSafe(~, tbl, rows)
            %HIGHLIGHTTABLEROWSSAFE Highlight table rows when uistyle is available.

            if isempty(rows)
                return
            end

            try
                removeStyle(tbl);
            catch
            end

            try
                styleObj = uistyle('BackgroundColor', [1.00 0.95 0.65]);
                addStyle(tbl, styleObj, 'row', rows);
            catch
                % Styling is cosmetic. Ignore unsupported cases.
            end
        end

        function clearTableStylesSafe(~, tbl)
            %CLEARTABLESTYLESSAFE Remove table styles when supported.

            try
                removeStyle(tbl);
            catch
            end
        end

        function statusInfo = addRawFolderWarningToStatusInfo(app, stepLocal, statusInfo)
            %ADDRAWFOLDERWARNINGTOSTATUSINFO Add GUI warning for missing RawFolder.
            %
            % This is a display-only warning. Missing RawFolder can be valid for
            % pipelines that do not use RawFolder. The warning is attached only to
            % steps that explicitly declare a RawFolder input.

            if ~app.hasAnyMissingRawFolder()
                return
            end

            if ~app.stepUsesRawFolder(stepLocal)
                return
            end

            msg = "This step uses RawFolder, but at least one selected folder pair has RawFolder set to Missing.";

            if ~isfield(statusInfo, 'messages') || isempty(statusInfo.messages)
                statusInfo.messages = strings(0,1);
            end

            if ~any(strcmpi(string(statusInfo.messages), msg))
                statusInfo.messages(end+1,1) = msg;
            end

            % Preserve stronger backend states. The graph renderer already gives
            % blocking/non-blocking invalid states priority over warnings.
            if ~isfield(statusInfo, 'hasError') || ~statusInfo.hasError
                if ~isfield(statusInfo, 'hasInvalidNonblocking') || ~statusInfo.hasInvalidNonblocking
                    statusInfo.hasWarning = true;

                    if ~isfield(statusInfo, 'status') || strcmpi(string(statusInfo.status), "valid")
                        statusInfo.status = "warning";
                    end

                    if ~isfield(statusInfo, 'severity') || ...
                            strcmpi(string(statusInfo.severity), "none") || ...
                            strlength(strtrim(string(statusInfo.severity))) == 0
                        statusInfo.severity = "warning";
                    end
                end
            end
        end

        function tf = hasAnyMissingRawFolder(app)
            %HASANYMISSINGRAWFOLDER True if any active RawFolder is Missing.

            tf = false;

            rawFolders = strings(0,1);

            % Prefer the current folder-pair table because it tracks enabled rows.
            try
                if isprop(app, 'FolderPairTable') && istable(app.FolderPairTable) && ...
                        ~isempty(app.FolderPairTable) && height(app.FolderPairTable) > 0

                    T = app.FolderPairTable;
                    if ismember('Use', T.Properties.VariableNames)
                        T = T(logical(T.Use), :);
                    end

                    if ismember('RawFolder', T.Properties.VariableNames)
                        rawFolders = string(T.RawFolder(:));
                    end
                end
            catch
                rawFolders = strings(0,1);
            end

            % Fallback for externally attached PipelineManager objects.
            if isempty(rawFolders)
                try
                    if ~isempty(app.pm) && isprop(app.pm, 'RawFolderList') && ~isempty(app.pm.RawFolderList)
                        rawFolders = string(app.pm.RawFolderList(:));
                    end
                catch
                    rawFolders = strings(0,1);
                end
            end

            if isempty(rawFolders)
                return
            end

            rawFolders = strtrim(rawFolders);
            tf = any(strcmpi(rawFolders, "Missing") | strlength(rawFolders) == 0);
        end

        function tf = stepUsesRawFolder(~, stepLocal)
            %STEPUSESRAWFOLDER True when a step declares a RawFolder input.

            tf = false;

            if isempty(stepLocal) || ~isstruct(stepLocal)
                return
            end

            if ~isfield(stepLocal, 'kind') || ~strcmpi(char(string(stepLocal.kind)), 'stream')
                return
            end

            if ~isfield(stepLocal, 'info') || ~isfield(stepLocal.info, 'inputs') || ...
                    isempty(stepLocal.info.inputs)
                return
            end

            for iIn = 1:numel(stepLocal.info.inputs)
                inDef = stepLocal.info.inputs(iIn);

                inputName = "";
                if isfield(inDef, 'name') && ~isempty(inDef.name)
                    inputName = lower(strtrim(string(inDef.name)));
                end

                inputTypes = strings(0,1);
                if isfield(inDef, 'type') && ~isempty(inDef.type)
                    if ischar(inDef.type) || (isstring(inDef.type) && isscalar(inDef.type))
                        inputTypes = string(inDef.type);
                    else
                        inputTypes = string(inDef.type(:));
                    end
                end

                if strcmp(inputName, "rawfolder") || any(strcmpi(inputTypes, "RawFolder"))
                    tf = true;
                    return
                end
            end
        end

        function cleanupAuxiliaryFigures(app)
            %CLEANUPAUXILIARYFIGURES Delete non-modal child windows owned by the app.
            %
            % This method intentionally avoids restoring control states because it
            % can be called while the main app is being destroyed.

            try
                if ~isempty(app.ExecutionPlanFigure) && isvalid(app.ExecutionPlanFigure)
                    fig = app.ExecutionPlanFigure;
                    app.ExecutionPlanFigure = [];
                    delete(fig);
                else
                    app.ExecutionPlanFigure = [];
                end
            catch
                app.ExecutionPlanFigure = [];
            end

            try
                if ~isempty(app.PipelineSummaryFigure) && isvalid(app.PipelineSummaryFigure)
                    fig = app.PipelineSummaryFigure;
                    app.PipelineSummaryFigure = [];
                    delete(fig);
                else
                    app.PipelineSummaryFigure = [];
                end
            catch
                app.PipelineSummaryFigure = [];
            end
        end

        function closeExecutionPlanWindow(app)
            %CLOSEEXECUTIONPLANWINDOW Close the execution-plan viewer.

            try
                if ~isempty(app.ExecutionPlanFigure) && isvalid(app.ExecutionPlanFigure)
                    fig = app.ExecutionPlanFigure;
                    app.ExecutionPlanFigure = [];
                    delete(fig);
                else
                    app.ExecutionPlanFigure = [];
                end
            catch
                app.ExecutionPlanFigure = [];
            end
        end

        function showExecutionPlanInteractive(app)
            %SHOWEXECUTIONPLANINTERACTIVE Show backend execution plan in UI tables.

            if isempty(app.pm)
                app.appendDiagnostic('Show-execution-plan requested without an attached PipelineManager.');
                app.setStatus('No PipelineManager attached.');
                return
            end

            try
                app.applyExecutionControlsToManager();
                plan = app.getCurrentExecutionPlan();

                if isempty(plan)
                    error('Execution plan is empty.');
                end

                app.closeExecutionPlanWindow();

                mainPos = app.UIFigure.Position;
                figW = 1050;
                figH = 620;
                x = mainPos(1) + max(20, (mainPos(3) - figW) / 2);
                y = mainPos(2) + max(20, (mainPos(4) - figH) / 2);

                app.ExecutionPlanFigure = uifigure( ...
                    'Name', 'Execution Plan', ...
                    'Position', [x y figW figH], ...
                    'Resize', 'on', ...
                    'WindowStyle', 'normal', ...
                    'CloseRequestFcn', @(src, event) app.closeExecutionPlanWindow());

                mainGrid = uigridlayout(app.ExecutionPlanFigure);
                mainGrid.RowHeight = {32, '1x'};
                mainGrid.ColumnWidth = {'1x'};
                mainGrid.Padding = [10 10 10 10];
                mainGrid.RowSpacing = 8;

                nRunSteps = 0;
                nSources = 0;
                nLinks = 0;
                nOutputs = 0;
                policyText = "";
                runTable = table();
                sourceTable = table();

                if isstruct(plan)
                    if isfield(plan, 'nodeTable') && istable(plan.nodeTable)
                        runTable = app.buildRunOrderDisplayTableFromPlan(plan);
                        sourceTable = app.buildVirtualInputDisplayTableFromPlan(plan);
                        nRunSteps = height(runTable);
                        nSources = height(sourceTable);
                    end
                    if isfield(plan, 'edgeTable') && istable(plan.edgeTable)
                        nLinks = height(plan.edgeTable);
                    end
                    if isfield(plan, 'outputPlan') && istable(plan.outputPlan)
                        nOutputs = height(plan.outputPlan);
                    end
                    if isfield(plan, 'leafOutputPolicy')
                        policyText = string(plan.leafOutputPolicy);
                    end
                end

                summaryLabel = uilabel(mainGrid);
                summaryLabel.Layout.Row = 1;
                summaryLabel.Layout.Column = 1;
                summaryLabel.FontWeight = 'bold';
                if strlength(policyText) > 0
                    summaryLabel.Text = sprintf(['Execution plan: %d executable step(s), %d virtual input source(s), ' ...
                        '%d data link(s), %d planned output(s). Output policy: %s.'], ...
                        nRunSteps, nSources, nLinks, nOutputs, char(policyText));
                else
                    summaryLabel.Text = sprintf(['Execution plan: %d executable step(s), %d virtual input source(s), ' ...
                        '%d data link(s), %d planned output(s).'], ...
                        nRunSteps, nSources, nLinks, nOutputs);
                end

                tabGroup = uitabgroup(mainGrid);
                tabGroup.Layout.Row = 2;
                tabGroup.Layout.Column = 1;

                app.createExecutionPlanTableTab(tabGroup, 'Run Order', runTable, 'runOrder');
                app.createExecutionPlanTableTab(tabGroup, 'Virtual Inputs', sourceTable, 'virtualInputs');

                if isfield(plan, 'edgeTable') && istable(plan.edgeTable)
                    app.createExecutionPlanTableTab(tabGroup, 'Data Flow', plan.edgeTable, 'dataFlow');
                else
                    app.createExecutionPlanTableTab(tabGroup, 'Data Flow', table(), 'dataFlow');
                end

                if isfield(plan, 'outputPlan') && istable(plan.outputPlan)
                    app.createExecutionPlanTableTab(tabGroup, 'Output Handling', plan.outputPlan, 'outputHandling');
                else
                    app.createExecutionPlanTableTab(tabGroup, 'Output Handling', table(), 'outputHandling');
                end

                if isfield(plan, 'validationReport') && isstruct(plan.validationReport) && ...
                        isfield(plan.validationReport, 'nodeStatus') && istable(plan.validationReport.nodeStatus)
                    app.createExecutionPlanTableTab(tabGroup, 'Diagnostics', plan.validationReport.nodeStatus, 'diagnostics');
                else
                    app.createExecutionPlanTableTab(tabGroup, 'Diagnostics', table(), 'diagnostics');
                end

                app.appendDiagnostic('Opened execution plan window.');
                app.setStatus('Execution plan shown.');

            catch ME
                app.appendDiagnostic(sprintf('Failed to show execution plan: %s', ME.message));
                app.setStatus('Failed to show execution plan.');
            end
        end

        function createExecutionPlanTableTab(app, tabGroup, tabTitle, T, tableKind)
            %CREATEEXECUTIONPLANTABLETAB Create one table tab in the plan viewer.

            tab = uitab(tabGroup);
            tab.Title = tabTitle;

            grid = uigridlayout(tab);
            grid.RowHeight = {'1x'};
            grid.ColumnWidth = {'1x'};
            grid.Padding = [8 8 8 8];

            if ~istable(T) || isempty(T)
                label = uilabel(grid);
                label.Layout.Row = 1;
                label.Layout.Column = 1;
                label.Text = sprintf('No %s information available.', lower(tabTitle));
                label.HorizontalAlignment = 'center';
                label.FontColor = [0.45 0.45 0.45];
                return
            end

            [displayTable, columnNames] = app.prepareExecutionPlanDisplayTable(T, tableKind);

            tbl = uitable(grid);
            tbl.Layout.Row = 1;
            tbl.Layout.Column = 1;
            tbl.Data = displayTable;
            tbl.ColumnName = columnNames;
            tbl.ColumnEditable = false(1, width(displayTable));
            tbl.RowName = 'numbered';
        end

        function [displayTable, columnNames] = prepareExecutionPlanDisplayTable(app, T, tableKind)
            %PREPAREEXECUTIONPLANDISPLAYTABLE Convert backend tables to UI-friendly tables.

            if nargin < 3
                tableKind = '';
            end

            displayTable = T;

            for iVar = 1:width(displayTable)
                varName = displayTable.Properties.VariableNames{iVar};
                values = displayTable.(varName);

                if iscell(values)
                    displayTable.(varName) = string(cellfun( ...
                        @(x) app.formatTableValueForDisplay(x), ...
                        values, ...
                        'UniformOutput', false));

                elseif isstruct(values)
                    displayTable.(varName) = arrayfun( ...
                        @(x) string(app.formatTableValueForDisplay(x)), ...
                        values);

                elseif isnumeric(values) || islogical(values) || isstring(values) || ...
                        isdatetime(values) || isduration(values) || iscategorical(values)
                    % These are supported by uitable/table display.
                else
                    displayTable.(varName) = arrayfun( ...
                        @(x) string(app.formatTableValueForDisplay(x)), ...
                        values);
                end
            end

            if strcmpi(tableKind, 'outputHandling')
                if ismember('plannedPersistence', displayTable.Properties.VariableNames)
                    displayTable.plannedPersistence = arrayfun( ...
                        @(x) string(app.formatPersistenceStateForDisplay(x)), ...
                        string(displayTable.plannedPersistence));
                end

                if ismember('reason', displayTable.Properties.VariableNames)
                    displayTable.reason = arrayfun( ...
                        @(x) string(app.formatOutputPlanReasonForDisplay(x)), ...
                        string(displayTable.reason));
                end
            end

            if ismember('executionRole', displayTable.Properties.VariableNames)
                displayTable.executionRole = arrayfun( ...
                    @(x) string(app.formatExecutionRoleForDisplay(x)), ...
                    string(displayTable.executionRole));
            end

            columnNames = app.getExecutionPlanColumnNames(displayTable.Properties.VariableNames, tableKind);
        end

        function txt = formatTableValueForDisplay(app, valueIn)
            %FORMATTABLEVALUEFORDISPLAY Convert one table element to readable text.

            if isempty(valueIn)
                txt = "";
                return
            end

            if ischar(valueIn) || (isstring(valueIn) && isscalar(valueIn))
                txt = char(string(valueIn));
                return
            end

            if iscell(valueIn)
                if isempty(valueIn)
                    txt = "";
                    return
                end

                try
                    asText = string(valueIn(:));
                    asText = asText(strlength(strtrim(asText)) > 0);
                    txt = char(strjoin(asText, ', '));
                catch
                    txt = sprintf('[%s cell]', strjoin(string(size(valueIn)), 'x'));
                end
                return
            end

            if isnumeric(valueIn) || islogical(valueIn)
                if isscalar(valueIn) || (isvector(valueIn) && numel(valueIn) <= 8)
                    txt = mat2str(valueIn);
                else
                    txt = sprintf('[%s %s]', strjoin(string(size(valueIn)), 'x'), class(valueIn));
                end
                return
            end

            if isdatetime(valueIn) || isduration(valueIn) || iscategorical(valueIn)
                txt = char(string(valueIn));
                return
            end

            if isstruct(valueIn)
                if isscalar(valueIn)
                    txt = sprintf('[struct: %s]', strjoin(string(fieldnames(valueIn)), ', '));
                else
                    txt = sprintf('[%s struct]', strjoin(string(size(valueIn)), 'x'));
                end
                return
            end

            txt = app.formatValueForDisplay(valueIn);
        end

        function columnNames = getExecutionPlanColumnNames(~, varNames, tableKind)
            %GETEXECUTIONPLANCOLUMNNAMES Return user-facing column labels.

            if nargin < 3
                tableKind = '';
            end

            columnNames = cellstr(string(varNames));

            for iName = 1:numel(columnNames)
                nameLower = lower(columnNames{iName});

                switch nameLower
                    case 'nodeid'
                        columnNames{iName} = 'Step ID';
                    case {'stepname', 'steptag'}
                        columnNames{iName} = 'Step';
                    case 'functionname'
                        columnNames{iName} = 'Function';
                    case 'kind'
                        columnNames{iName} = 'Kind';
                    case 'topoindex'
                        columnNames{iName} = 'Topological order';
                    case 'scheduleindex'
                        columnNames{iName} = 'Schedule index';
                    case 'runindex'
                        columnNames{iName} = 'Run order';
                    case 'isexecutable'
                        columnNames{iName} = 'Executable';
                    case 'executionrole'
                        columnNames{iName} = 'Role';
                    case 'sourceindex'
                        columnNames{iName} = 'Source';
                    case 'sourcefile'
                        columnNames{iName} = 'Source file';
                    case 'feedsstep'
                        columnNames{iName} = 'Feeds step';
                    case 'sourcetype'
                        columnNames{iName} = 'Source type';
                    case 'isroot'
                        columnNames{iName} = 'Starts branch';
                    case 'isleaf'
                        columnNames{iName} = 'Final step';
                    case 'isfinaloutput'
                        columnNames{iName} = 'Final output';
                    case 'sourcenodeid'
                        columnNames{iName} = 'Source step ID';
                    case 'sourcestep'
                        columnNames{iName} = 'Source step';
                    case 'sourceoutput'
                        columnNames{iName} = 'Source output';
                    case 'selectedfile'
                        columnNames{iName} = 'Selected file';
                    case 'targetnodeid'
                        columnNames{iName} = 'Target step ID';
                    case 'targetstep'
                        columnNames{iName} = 'Target step';
                    case 'targetinput'
                        columnNames{iName} = 'Target input';
                    case 'sourcetypes'
                        columnNames{iName} = 'Source type';
                    case 'outputname'
                        columnNames{iName} = 'Output';
                    case 'outputmode'
                        columnNames{iName} = 'Output mode';
                    case 'plannedpersistence'
                        columnNames{iName} = 'Handling';
                    case 'plannedfilename'
                        columnNames{iName} = 'Planned file';
                    case 'reason'
                        columnNames{iName} = 'Reason';
                    case 'status'
                        columnNames{iName} = 'Status';
                    case 'severity'
                        columnNames{iName} = 'Severity';
                    case 'messages'
                        columnNames{iName} = 'Messages';
                    otherwise
                        % Keep backend name if no user-facing alias is known.
                end
            end

            %#ok<NASGU> tableKind is kept for future table-specific aliases.
        end

        function value = emptyToNaN(~, value)
            %EMPTYTONAN Convert [] to NaN for HTML payload scalars.

            if isempty(value)
                value = NaN;
            end
        end

        function scheduleOrder = getScheduleOrderFromPlan(~, plan)
            %GETSCHEDULEORDERFROMPLAN Extract backend scheduler/visit order from a plan.

            scheduleOrder = [];

            if isstruct(plan) && isfield(plan, 'nodeOrders') && isstruct(plan.nodeOrders) && ...
                    isfield(plan.nodeOrders, 'scheduleOrder') && ~isempty(plan.nodeOrders.scheduleOrder)
                scheduleOrder = plan.nodeOrders.scheduleOrder(:).';
                return
            end

            if isstruct(plan) && isfield(plan, 'nodeTable') && istable(plan.nodeTable) && ...
                    ~isempty(plan.nodeTable) && all(ismember({'nodeID','scheduleIndex'}, plan.nodeTable.Properties.VariableNames))
                nodeTable = plan.nodeTable;
                keep = ~isnan(nodeTable.scheduleIndex);
                nodeTable = nodeTable(keep, :);
                [~, ord] = sort(nodeTable.scheduleIndex, 'ascend');
                scheduleOrder = nodeTable.nodeID(ord).';
            end
        end

        function runOrder = getRunOrderFromPlan(~, plan)
            %GETRUNORDERFROMPLAN Extract executable stream-step order from a plan.

            runOrder = [];

            if isstruct(plan) && isfield(plan, 'nodeOrders') && isstruct(plan.nodeOrders) && ...
                    isfield(plan.nodeOrders, 'runOrder') && ~isempty(plan.nodeOrders.runOrder)
                runOrder = plan.nodeOrders.runOrder(:).';
                return
            end

            if isstruct(plan) && isfield(plan, 'nodeTable') && istable(plan.nodeTable) && ...
                    ~isempty(plan.nodeTable) && all(ismember({'nodeID','runIndex'}, plan.nodeTable.Properties.VariableNames))
                nodeTable = plan.nodeTable;
                if ismember('isExecutable', nodeTable.Properties.VariableNames)
                    nodeTable = nodeTable(logical(nodeTable.isExecutable), :);
                end
                keep = ~isnan(nodeTable.runIndex);
                nodeTable = nodeTable(keep, :);
                [~, ord] = sort(nodeTable.runIndex, 'ascend');
                runOrder = nodeTable.nodeID(ord).';
            end
        end

        function orderInfo = getNodeOrderMetadataFromPlan(app, plan, nodeID)
            %GETNODEORDERMETADATAFROMPLAN Return order metadata for one node.

            orderInfo = struct( ...
                'topoIndex', [], ...
                'scheduleIndex', [], ...
                'runIndex', [], ...
                'isExecutable', false, ...
                'executionRole', "");

            if ~isstruct(plan) || ~isfield(plan, 'nodeTable')
                return
            end

            row = app.getNodeTableRowByID(plan.nodeTable, nodeID);
            if isempty(row)
                return
            end

            if ismember('topoIndex', row.Properties.VariableNames) && ~isnan(row.topoIndex(1))
                orderInfo.topoIndex = row.topoIndex(1);
            end
            if ismember('scheduleIndex', row.Properties.VariableNames) && ~isnan(row.scheduleIndex(1))
                orderInfo.scheduleIndex = row.scheduleIndex(1);
            end
            if ismember('runIndex', row.Properties.VariableNames) && ~isnan(row.runIndex(1))
                orderInfo.runIndex = row.runIndex(1);
            end
            if ismember('isExecutable', row.Properties.VariableNames)
                orderInfo.isExecutable = logical(row.isExecutable(1));
            end
            if ismember('executionRole', row.Properties.VariableNames)
                orderInfo.executionRole = string(row.executionRole(1));
            elseif orderInfo.isExecutable
                orderInfo.executionRole = "analysis_step";
            else
                orderInfo.executionRole = "virtual_source";
            end
        end

        function T = buildRunOrderDisplayTableFromPlan(app, plan)
            %BUILDRUNORDERDISPLAYTABLEFROMPLAN Build the user-facing run-order table.

            T = table();

            if ~isstruct(plan) || ~isfield(plan, 'nodeTable') || ~istable(plan.nodeTable) || isempty(plan.nodeTable)
                return
            end

            rows = app.getExecutableRowsFromNodeTable(plan.nodeTable);
            if isempty(rows)
                return
            end

            [~, ord] = sort(rows.runIndex, 'ascend', 'MissingPlacement', 'last');
            rows = rows(ord, :);

            keepVars = {'runIndex','nodeID','stepName','functionName','scheduleIndex','topoIndex','isLeaf'};
            keepVars = keepVars(ismember(keepVars, rows.Properties.VariableNames));
            T = rows(:, keepVars);
        end

        function T = buildVirtualInputDisplayTableFromPlan(app, plan)
            %BUILDVIRTUALINPUTDISPLAYTABLEFROMPLAN Build the virtual-input table.

            T = table();

            if ~isstruct(plan) || ~isfield(plan, 'nodeTable') || ~istable(plan.nodeTable) || isempty(plan.nodeTable)
                return
            end

            sourceRows = app.getVirtualSourceRowsFromNodeTable(plan.nodeTable);
            if isempty(sourceRows)
                return
            end

            [~, ord] = sort(sourceRows.scheduleIndex, 'ascend', 'MissingPlacement', 'last');
            sourceRows = sourceRows(ord, :);

            edgeTable = table();
            if isfield(plan, 'edgeTable') && istable(plan.edgeTable)
                edgeTable = plan.edgeTable;
            end

            sourceIndex = (1:height(sourceRows)).';
            nodeID = sourceRows.nodeID;
            sourceFile = sourceRows.stepName;
            scheduleIndex = sourceRows.scheduleIndex;
            feedsStep = strings(height(sourceRows), 1);
            sourceType = strings(height(sourceRows), 1);

            for iRow = 1:height(sourceRows)
                if istable(edgeTable) && ~isempty(edgeTable) && ...
                        all(ismember({'sourceNodeID','targetStep'}, edgeTable.Properties.VariableNames))

                    eRows = edgeTable(edgeTable.sourceNodeID == sourceRows.nodeID(iRow), :);

                    if ~isempty(eRows)
                        feedsStep(iRow) = strjoin(unique(string(eRows.targetStep), 'stable'), ', ');

                        if ismember('sourceTypes', eRows.Properties.VariableNames)
                            typeVals = strings(0,1);
                            for iEdge = 1:height(eRows)
                                try
                                    typeVals = [typeVals; string(eRows.sourceTypes{iEdge}(:))]; %#ok<AGROW>
                                catch
                                end
                            end
                            typeVals = unique(typeVals(strlength(strtrim(typeVals)) > 0), 'stable');
                            if ~isempty(typeVals)
                                sourceType(iRow) = strjoin(typeVals, ' / ');
                            end
                        end
                    end
                end
            end

            T = table(sourceIndex, nodeID, sourceFile, scheduleIndex, feedsStep, sourceType);
        end

        function rows = getExecutableRowsFromNodeTable(~, nodeTable)
            %GETEXECUTABLEROWSFROMNODETABLE Return stream/executable rows from nodeTable.

            rows = table();
            if ~istable(nodeTable) || isempty(nodeTable)
                return
            end

            if ismember('isExecutable', nodeTable.Properties.VariableNames)
                rows = nodeTable(logical(nodeTable.isExecutable), :);
            elseif ismember('kind', nodeTable.Properties.VariableNames)
                rows = nodeTable(strcmpi(string(nodeTable.kind), "stream"), :);
            end
        end

        function rows = getVirtualSourceRowsFromNodeTable(~, nodeTable)
            %GETVIRTUALSOURCEROWSFROMNODETABLE Return virtual FileSource rows.

            rows = table();
            if ~istable(nodeTable) || isempty(nodeTable)
                return
            end

            if ismember('executionRole', nodeTable.Properties.VariableNames)
                rows = nodeTable(strcmpi(string(nodeTable.executionRole), "virtual_source"), :);
            elseif ismember('isExecutable', nodeTable.Properties.VariableNames)
                rows = nodeTable(~logical(nodeTable.isExecutable), :);
            elseif ismember('kind', nodeTable.Properties.VariableNames)
                rows = nodeTable(strcmpi(string(nodeTable.kind), "folder"), :);
            end
        end

        function txt = formatExecutionRoleForDisplay(~, roleText)
            %FORMATEXECUTIONROLEFORDISPLAY Convert backend role to user-facing text.

            roleText = lower(strtrim(string(roleText)));

            switch roleText
                case "virtual_source"
                    txt = "virtual input source";
                case "analysis_step"
                    txt = "analysis step";
                otherwise
                    txt = string(roleText);
            end
        end

        function setMenuControlsAvailable(app, tf)
            %SETMENUCONTROLSAVAILABLE Enable/disable menu actions requiring backend context.

            if tf
                state = 'on';
            else
                state = 'off';
            end

            menuItems = { ...
                app.LoadPipelineMenu, ...
                app.SavePipelineMenu, ...
                app.GenerateScriptMenu, ...
                app.RunPipelineMenu, ...
                app.ReuseExistingFilesMenu, ...
                app.ExecutionSettingsMenu, ...
                app.AdvancedRAMSettingsMenu, ...
                app.ShowExecutionPlanMenu, ...
                app.ExportErrorLogMenu, ...
                app.RemoveInvalidStepsMenu, ...
                app.ViewFolderPipelineLogMenu};


            for iMenu = 1:numel(menuItems)
                try
                    menuItems{iMenu}.Enable = state;
                catch
                end
            end

            try
                if app.bDataViewerMode
                    app.SelectDataFoldersMenu.Enable = 'off';
                else
                    app.SelectDataFoldersMenu.Enable = 'on';
                end
            catch
            end


            app.updateLatestRunSummaryMenuState();
        end

        function tf = getReuseExistingFilesUiValue(app)
            %GETREUSEEXISTINGFILESUIVALUE Return the reuse-existing-files UI state.

            tf = true;

            try
                tf = logical(app.ReuseExistingFilesCheckBox.Value);
                return
            catch
            end

            try
                tf = strcmpi(app.ReuseExistingFilesMenu.Checked, 'on');
            catch
            end
        end

        function setReuseExistingFilesState(app, tf, pushToManager)
            %SETREUSEEXISTINGFILESSTATE Synchronize checkbox, menu, and backend flag.
            %
            % User-facing label:
            %   Reuse Existing Files
            %
            % Backend mapping:
            %   PipelineManager.b_skipSteps

            if nargin < 3 || isempty(pushToManager)
                pushToManager = true;
            end

            tf = logical(tf);

            try
                app.ReuseExistingFilesCheckBox.Value = tf;
            catch
            end

            try
                if tf
                    app.ReuseExistingFilesMenu.Checked = 'on';
                else
                    app.ReuseExistingFilesMenu.Checked = 'off';
                end
            catch
            end

            if pushToManager && ~isempty(app.pm)
                try
                    if isprop(app.pm, 'b_skipSteps')
                        app.pm.b_skipSteps = tf;
                    end
                catch ME
                    app.appendDiagnostic(sprintf('Failed to update reuse-existing-files setting: %s', ME.message));
                end
            end
        end

        function showExecutionSettingsInteractive(app)
            %SHOWEXECUTIONSETTINGSINTERACTIVE Edit non-RAM execution behavior.

            if isempty(app.pm)
                app.appendDiagnostic('Execution-settings requested without an attached PipelineManager.');
                app.setStatus('No PipelineManager attached.');
                return
            end

            reuseExisting = true;
            avoidOverwrite = true;
            saveBeforeFail = false;

            try
                if isprop(app.pm, 'b_skipSteps')
                    reuseExisting = logical(app.pm.b_skipSteps);
                else
                    reuseExisting = app.getReuseExistingFilesUiValue();
                end
            catch
                reuseExisting = app.getReuseExistingFilesUiValue();
            end

            try
                if isprop(app.pm, 'b_avoidOverwrite')
                    avoidOverwrite = logical(app.pm.b_avoidOverwrite);
                end
            catch
            end

            try
                if isprop(app.pm, 'b_saveDataBeforeFail')
                    saveBeforeFail = logical(app.pm.b_saveDataBeforeFail);
                end
            catch
            end

            mainPos = app.UIFigure.Position;
            figW = 470;
            figH = 230;
            x = mainPos(1) + max(20, (mainPos(3) - figW) / 2);
            y = mainPos(2) + max(20, (mainPos(4) - figH) / 2);

            dlg = uifigure( ...
                'Name', 'Execution Settings', ...
                'Position', [x y figW figH], ...
                'Resize', 'off', ...
                'WindowStyle', 'modal');

            grid = uigridlayout(dlg);
            grid.RowHeight = {32, 34, 34, 34, '1x', 36};
            grid.ColumnWidth = {'1x', 80, 80};
            grid.Padding = [12 12 12 12];
            grid.RowSpacing = 6;
            grid.ColumnSpacing = 8;

            titleLabel = uilabel(grid);
            titleLabel.Text = 'Execution behavior';
            titleLabel.FontWeight = 'bold';
            titleLabel.Layout.Row = 1;
            titleLabel.Layout.Column = [1 3];

            reuseBox = uicheckbox(grid);
            reuseBox.Text = 'Reuse existing files';
            reuseBox.Value = reuseExisting;
            reuseBox.Tooltip = 'When enabled, data files that already exist and match pipeline history may be reused. Disable to force re-execution.';
            reuseBox.Layout.Row = 2;
            reuseBox.Layout.Column = [1 3];

            overwriteBox = uicheckbox(grid);
            overwriteBox.Text = 'Avoid overwriting existing files';
            overwriteBox.Value = avoidOverwrite;
            overwriteBox.Tooltip = 'When enabled, PipelineManager avoids overwriting existing output files when possible.';
            overwriteBox.Layout.Row = 3;
            overwriteBox.Layout.Column = [1 3];

            saveFailBox = uicheckbox(grid);
            saveFailBox.Text = 'Save recovery data on failure';
            saveFailBox.Value = saveBeforeFail;
            saveFailBox.Tooltip = 'When enabled, PipelineManager saves the latest available data before a step failure when possible.';
            saveFailBox.Layout.Row = 4;
            saveFailBox.Layout.Column = [1 3];

            applyButton = uibutton(grid, 'push');
            applyButton.Text = 'Apply';
            applyButton.FontWeight = 'bold';
            applyButton.Layout.Row = 6;
            applyButton.Layout.Column = 2;
            applyButton.ButtonPushedFcn = @(~,~) applySettings();

            cancelButton = uibutton(grid, 'push');
            cancelButton.Text = 'Cancel';
            cancelButton.Layout.Row = 6;
            cancelButton.Layout.Column = 3;
            cancelButton.ButtonPushedFcn = @(~,~) cancelSettings();

            uiwait(dlg);

            if isvalid(dlg)
                delete(dlg);
            end

            function applySettings()
                try
                    app.setReuseExistingFilesState(reuseBox.Value, true);

                    if isprop(app.pm, 'b_avoidOverwrite')
                        app.pm.b_avoidOverwrite = logical(overwriteBox.Value);
                    end

                    if isprop(app.pm, 'b_saveDataBeforeFail')
                        app.pm.b_saveDataBeforeFail = logical(saveFailBox.Value);
                    end

                    app.refreshAllViews();
                    app.appendDiagnostic('Execution settings updated.');
                    app.setStatus('Execution settings updated.');
                    uiresume(dlg);

                catch ME
                    uialert(dlg, ME.message, 'Failed to Apply Execution Settings');
                end
            end

            function cancelSettings()
                app.appendDiagnostic('Execution-settings dialog cancelled.');
                app.setStatus('Execution settings cancelled.');
                uiresume(dlg);
            end
        end

        function showAdvancedRamSettingsInteractive(app)
            %SHOWADVANCEDRAMSETTINGSINTERACTIVE Edit advanced RAM planning settings.

            if isempty(app.pm)
                app.appendDiagnostic('Advanced-RAM-settings requested without an attached PipelineManager.');
                app.setStatus('No PipelineManager attached.');
                return
            end

            ramLookaheadLength = app.getNumericManagerProperty('ramLookaheadLength', 1);
            ramCooldownSteps = app.getNumericManagerProperty('ramCooldownSteps', 0);
            ramSafetyOverheadPct = app.getNumericManagerProperty('ramSafetyOverheadPct', 0.20);
            ramSafetyDataMultiplier = app.getNumericManagerProperty('ramSafetyDataMultiplier', 1.00);

            mainPos = app.UIFigure.Position;
            figW = 520;
            figH = 290;
            x = mainPos(1) + max(20, (mainPos(3) - figW) / 2);
            y = mainPos(2) + max(20, (mainPos(4) - figH) / 2);

            dlg = uifigure( ...
                'Name', 'Advanced RAM Settings', ...
                'Position', [x y figW figH], ...
                'Resize', 'off', ...
                'WindowStyle', 'modal');

            grid = uigridlayout(dlg);
            grid.RowHeight = {32, 34, 34, 34, 34, '1x', 36};
            grid.ColumnWidth = {190, '1x', 80, 80};
            grid.Padding = [12 12 12 12];
            grid.RowSpacing = 6;
            grid.ColumnSpacing = 8;

            titleLabel = uilabel(grid);
            titleLabel.Text = 'Advanced RAM planning settings';
            titleLabel.FontWeight = 'bold';
            titleLabel.Layout.Row = 1;
            titleLabel.Layout.Column = [1 4];

            lookaheadLabel = uilabel(grid);
            lookaheadLabel.Text = 'Lookahead length';
            lookaheadLabel.Layout.Row = 2;
            lookaheadLabel.Layout.Column = 1;

            lookaheadField = uieditfield(grid, 'numeric');
            lookaheadField.Value = ramLookaheadLength;
            lookaheadField.Limits = [0 Inf];
            lookaheadField.RoundFractionalValues = 'on';
            lookaheadField.Tooltip = 'Number of downstream steps considered during RAM planning.';
            lookaheadField.Layout.Row = 2;
            lookaheadField.Layout.Column = [2 4];

            cooldownLabel = uilabel(grid);
            cooldownLabel.Text = 'Cooldown steps';
            cooldownLabel.Layout.Row = 3;
            cooldownLabel.Layout.Column = 1;

            cooldownField = uieditfield(grid, 'numeric');
            cooldownField.Value = ramCooldownSteps;
            cooldownField.Limits = [0 Inf];
            cooldownField.RoundFractionalValues = 'on';
            cooldownField.Tooltip = 'Number of steps retained before intermediate RAM can be released.';
            cooldownField.Layout.Row = 3;
            cooldownField.Layout.Column = [2 4];

            overheadLabel = uilabel(grid);
            overheadLabel.Text = 'Safety overhead fraction';
            overheadLabel.Layout.Row = 4;
            overheadLabel.Layout.Column = 1;

            overheadField = uieditfield(grid, 'numeric');
            overheadField.Value = ramSafetyOverheadPct;
            overheadField.Limits = [0 Inf];
            overheadField.Tooltip = 'Extra RAM fraction reserved as a safety margin, for example 0.20 for 20%.';
            overheadField.Layout.Row = 4;
            overheadField.Layout.Column = [2 4];

            multiplierLabel = uilabel(grid);
            multiplierLabel.Text = 'Data size multiplier';
            multiplierLabel.Layout.Row = 5;
            multiplierLabel.Layout.Column = 1;

            multiplierField = uieditfield(grid, 'numeric');
            multiplierField.Value = ramSafetyDataMultiplier;
            multiplierField.Limits = [0 Inf];
            multiplierField.Tooltip = 'Multiplier applied to estimated data size during RAM planning.';
            multiplierField.Layout.Row = 5;
            multiplierField.Layout.Column = [2 4];

            applyButton = uibutton(grid, 'push');
            applyButton.Text = 'Apply';
            applyButton.FontWeight = 'bold';
            applyButton.Layout.Row = 7;
            applyButton.Layout.Column = 3;
            applyButton.ButtonPushedFcn = @(~,~) applyRamSettings();

            cancelButton = uibutton(grid, 'push');
            cancelButton.Text = 'Cancel';
            cancelButton.Layout.Row = 7;
            cancelButton.Layout.Column = 4;
            cancelButton.ButtonPushedFcn = @(~,~) cancelRamSettings();

            uiwait(dlg);

            if isvalid(dlg)
                delete(dlg);
            end

            function applyRamSettings()
                try
                    app.setNumericManagerProperty('ramLookaheadLength', max(0, round(lookaheadField.Value)));
                    app.setNumericManagerProperty('ramCooldownSteps', max(0, round(cooldownField.Value)));
                    app.setNumericManagerProperty('ramSafetyOverheadPct', max(0, overheadField.Value));
                    app.setNumericManagerProperty('ramSafetyDataMultiplier', max(0, multiplierField.Value));

                    app.refreshAllViews();
                    app.appendDiagnostic('Advanced RAM settings updated.');
                    app.setStatus('Advanced RAM settings updated.');
                    uiresume(dlg);

                catch ME
                    uialert(dlg, ME.message, 'Failed to Apply Advanced RAM Settings');
                end
            end

            function cancelRamSettings()
                app.appendDiagnostic('Advanced-RAM-settings dialog cancelled.');
                app.setStatus('Advanced RAM settings cancelled.');
                uiresume(dlg);
            end
        end

        function value = getNumericManagerProperty(app, propName, defaultValue)
            %GETNUMERICMANAGERPROPERTY Read one numeric PipelineManager property.

            value = defaultValue;

            try
                if ~isempty(app.pm) && isprop(app.pm, propName)
                    candidate = app.pm.(propName);
                    if isnumeric(candidate) && isscalar(candidate) && isfinite(candidate)
                        value = double(candidate);
                    end
                end
            catch
                value = defaultValue;
            end
        end

        function setNumericManagerProperty(app, propName, value)
            %SETNUMERICMANAGERPROPERTY Write one numeric PipelineManager property.

            if isempty(app.pm) || ~isprop(app.pm, propName)
                return
            end

            if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value)
                error('PipelineManagerTool:setNumericManagerProperty:InvalidValue', ...
                    'Invalid value for %s.', propName);
            end

            app.pm.(propName) = value;
        end

        function restoreRunControlsAfterExecution(app)
            %RESTORERUNCONTROLSAFTEREXECUTION Backward-compatible wrapper.

            app.restoreGuiAfterPipelineRun();
        end

        function updateLatestRunSummaryMenuState(app)
            %UPDATELATESTRUNSUMMARYMENUSTATE Enable summary menu only when a run exists.

            try
                if ~isempty(app.LastExecutionResult)
                    app.ShowLatestRunSummaryMenu.Enable = 'on';
                else
                    app.ShowLatestRunSummaryMenu.Enable = 'off';
                end
            catch
            end
        end

        function closePipelineSummaryWindow(app, varargin)
            %CLOSEPIPELINESUMMARYWINDOW Close the latest-run summary viewer.
            %
            % Name-Value:
            %   CloseToolIfDataViewer - when true, closing the summary also closes this
            %                           modal PM Tool and notifies DataViewer.

            p = inputParser;
            addParameter(p, 'CloseToolIfDataViewer', false, @(x)islogical(x)&&isscalar(x));
            parse(p, varargin{:});

            closeToolIfDataViewer = logical(p.Results.CloseToolIfDataViewer);

            try
                if ~isempty(app.PipelineSummaryFigure) && isvalid(app.PipelineSummaryFigure)
                    fig = app.PipelineSummaryFigure;
                    app.PipelineSummaryFigure = [];
                    delete(fig);
                else
                    app.PipelineSummaryFigure = [];
                end
            catch
                app.PipelineSummaryFigure = [];
            end

            if closeToolIfDataViewer && app.bDataViewerMode
                app.closeDataViewerManagedTool();
            end
        end

        function closeDataViewerManagedTool(app)
            %CLOSEDATAVIEWERMANAGEDTOOL Close PM Tool and notify DataViewer.

            if ~app.bDataViewerMode
                return
            end

            toolState = struct();
            toolState.status = "closed";
            toolState.closedOn = datetime('now');
            toolState.lastExecutionResult = app.LastExecutionResult;
            toolState.dataViewerFileSourceNodeID = app.DataViewerFileSourceNodeID;

            try
                if ~isempty(app.DataViewerToolClosedFcn)
                    app.DataViewerToolClosedFcn(toolState);
                end
            catch ME
                try
                    app.appendDiagnostic(sprintf('DataViewer close callback failed: %s', ME.message));
                catch
                end
            end

            try
                app.UIFigure.WindowStyle = 'normal';
            catch
            end

            try
                delete(app.UIFigure);
            catch
            end
        end

        function showPipelineSummaryInteractive(app, result)
            %SHOWPIPELINESUMMARYINTERACTIVE Show latest execution result in UI tables.

            if nargin < 2 || isempty(result)
                result = app.LastExecutionResult;
            end

            if isempty(result)
                app.appendDiagnostic('No pipeline run summary is available yet.');
                app.setStatus('No latest run summary.');
                return
            end

            app.LastExecutionResult = result;
            app.updateLatestRunSummaryMenuState();

            app.closePipelineSummaryWindow('CloseToolIfDataViewer', false);

            mainPos = app.UIFigure.Position;
            figW = 1120;
            figH = 650;
            x = mainPos(1) + max(20, (mainPos(3) - figW) / 2);
            y = mainPos(2) + max(20, (mainPos(4) - figH) / 2);

            app.PipelineSummaryFigure = uifigure( ...
                'Name', 'Pipeline Summary', ...
                'Position', [x y figW figH], ...
                'Resize', 'on', ...
                'WindowStyle', 'modal', ...
                'CloseRequestFcn', @(~,~) app.closePipelineSummaryWindow( ...
                'CloseToolIfDataViewer', true));

            mainGrid = uigridlayout(app.PipelineSummaryFigure);
            mainGrid.RowHeight = {32, '1x'};
            mainGrid.ColumnWidth = {'1x'};
            mainGrid.Padding = [10 10 10 10];
            mainGrid.RowSpacing = 8;

            headerGrid = uigridlayout(mainGrid);
            headerGrid.RowHeight = {'1x'};
            headerGrid.ColumnWidth = {'1x', 110};
            headerGrid.Padding = [0 0 0 0];
            headerGrid.ColumnSpacing = 8;
            headerGrid.Layout.Row = 1;
            headerGrid.Layout.Column = 1;

            titleLabel = uilabel(headerGrid);
            titleLabel.FontWeight = 'bold';
            titleLabel.Text = app.buildPipelineSummaryHeaderText(result);
            titleLabel.Layout.Row = 1;
            titleLabel.Layout.Column = 1;

            closeButton = uibutton(headerGrid, 'push');
            closeButton.Text = 'Close';
            closeButton.Layout.Row = 1;
            closeButton.Layout.Column = 2;
            closeButton.ButtonPushedFcn = @(~,~) app.closePipelineSummaryWindow( ...
                'CloseToolIfDataViewer', true);

            tabGroup = uitabgroup(mainGrid);
            tabGroup.Layout.Row = 2;
            tabGroup.Layout.Column = 1;

            % ---------------------------------------------------------------------
            % Global log tab
            % ---------------------------------------------------------------------
            logTab = uitab(tabGroup);
            logTab.Title = 'Run Log';

            logGrid = uigridlayout(logTab);
            logGrid.RowHeight = {'1x'};
            logGrid.ColumnWidth = {'1x'};

            logTable = uitable(logGrid);
            logTable.Layout.Row = 1;
            logTable.Layout.Column = 1;
            logTable.RowName = {};

            if isfield(result, 'globalPipeLog') && istable(result.globalPipeLog)
                logTable.Data = result.globalPipeLog;
            else
                logTable.Data = table();
            end

            % ---------------------------------------------------------------------
            % Created files tab
            % ---------------------------------------------------------------------
            createdTab = uitab(tabGroup);
            createdTab.Title = 'Created Files';

            createdGrid = uigridlayout(createdTab);
            createdGrid.RowHeight = {'1x'};
            createdGrid.ColumnWidth = {'1x'};

            createdTable = uitable(createdGrid);
            createdTable.Layout.Row = 1;
            createdTable.Layout.Column = 1;
            createdTable.RowName = {};

            if isfield(result, 'createdFiles') && istable(result.createdFiles)
                createdTable.Data = result.createdFiles;
            else
                createdTable.Data = table();
            end

            % ---------------------------------------------------------------------
            % Final outputs tab
            % ---------------------------------------------------------------------
            finalTab = uitab(tabGroup);
            finalTab.Title = 'Final Outputs';

            finalGrid = uigridlayout(finalTab);
            finalGrid.RowHeight = {'1x'};
            finalGrid.ColumnWidth = {'1x'};

            finalTable = uitable(finalGrid);
            finalTable.Layout.Row = 1;
            finalTable.Layout.Column = 1;
            finalTable.RowName = {};

            if isfield(result, 'outputManifest') && istable(result.outputManifest)
                finalTable.Data = result.outputManifest;
            else
                finalTable.Data = table();
            end

            % ---------------------------------------------------------------------
            % Progress tab
            % ---------------------------------------------------------------------
            progressTab = uitab(tabGroup);
            progressTab.Title = 'Progress';

            progressGrid = uigridlayout(progressTab);
            progressGrid.RowHeight = {'1x'};
            progressGrid.ColumnWidth = {'1x'};

            progressTable = uitable(progressGrid);
            progressTable.Layout.Row = 1;
            progressTable.Layout.Column = 1;
            progressTable.RowName = {};

            if isfield(result, 'progressLog') && ~isempty(result.progressLog)
                try
                    progressTable.Data = struct2table(result.progressLog);
                catch
                    progressTable.Data = table();
                end
            else
                progressTable.Data = table();
            end
        end

        function txt = buildPipelineSummaryHeaderText(~, result)
            %BUILDPIPELINESUMMARYHEADERTEXT Compact summary heading.

            statusText = "unknown";
            durationText = "";

            try
                if isfield(result, 'status')
                    statusText = string(result.status);
                end
            catch
            end

            try
                if isfield(result, 'duration_s') && ~isempty(result.duration_s) && isfinite(result.duration_s)
                    durationText = sprintf(' | %.2f s', double(result.duration_s));
                end
            catch
                durationText = "";
            end

            txt = sprintf('Pipeline Summary | Status: %s%s', char(statusText), char(durationText));
        end

        function T = buildPipelineSummaryOverviewTable(app, result)
            %BUILDPIPELINESUMMARYOVERVIEWTABLE Build key/value overview table.

            fieldNames = strings(0,1);
            values = strings(0,1);

            addRow("Run status", app.getResultScalarString(result, 'status', '<unknown>'));

            if isstruct(result) && isfield(result, 'duration_s') && isnumeric(result.duration_s) && isscalar(result.duration_s)
                addRow("Duration (s)", string(round(result.duration_s, 3)));
            end

            addRow("Output policy", app.getResultScalarString(result, 'leafOutputPolicy', '<unknown>'));
            addRow("RAM mode", app.getResultScalarString(result, 'ramMode', '<unknown>'));
            addRow("RAM safety", app.getResultScalarString(result, 'ramSafePolicy', '<unknown>'));

            nFolders = 0;
            nCompleted = 0;
            nSkipped = 0;
            nFailed = 0;

            if isstruct(result) && isfield(result, 'globalPipeLog') && istable(result.globalPipeLog)
                logTable = result.globalPipeLog;
                nFolders = height(logTable);

                if ismember('Status', logTable.Properties.VariableNames)
                    statusVals = lower(string(logTable.Status));
                    nCompleted = sum(contains(statusVals, "completed"));
                    nSkipped = sum(contains(statusVals, "skipped"));
                    nFailed = sum(contains(statusVals, "failed") | contains(statusVals, "error"));
                end
            end

            addRow("Folder runs", string(nFolders));
            addRow("Completed folders", string(nCompleted));
            addRow("Skipped folders", string(nSkipped));
            addRow("Failed folders", string(nFailed));

            nCreated = 0;
            nTemporary = 0;
            if isstruct(result) && isfield(result, 'createdFiles') && istable(result.createdFiles)
                createdTable = result.createdFiles;
                nCreated = height(createdTable);
                if ismember('IsTemporary', createdTable.Properties.VariableNames)
                    nTemporary = sum(logical(createdTable.IsTemporary));
                end
            end

            addRow("Created files", string(nCreated));
            addRow("Temporary files", string(nTemporary));

            viewerTable = app.getExecutionResultTable(result, ...
                {'viewerCompatibleFinalOutputs', 'viewerCompatibleLeafOutputs'});
            unsavedTable = app.getExecutionResultTable(result, ...
                {'unsavedFinalOutputs', 'unsavedLeafOutputs'});

            addRow("Viewer-compatible final outputs", string(height(viewerTable)));
            addRow("Unsaved final outputs", string(height(unsavedTable)));

            T = table(fieldNames, values, 'VariableNames', {'Field','Value'});

            function addRow(name, value)
                fieldNames(end+1,1) = string(name); %#ok<AGROW>
                values(end+1,1) = string(value); %#ok<AGROW>
            end
        end

        function T = buildFinalOutputSummaryTable(app, result)
            %BUILDFINALOUTPUTSUMMARYTABLE Build a compact final-output manifest table.

            category = strings(0,1);
            saveFolder = strings(0,1);
            stepTag = strings(0,1);
            outputName = strings(0,1);
            handling = strings(0,1);
            fileName = strings(0,1);
            filePath = strings(0,1);

            viewerTable = app.getExecutionResultTable(result, ...
                {'viewerCompatibleFinalOutputs', 'viewerCompatibleLeafOutputs'});
            unsavedTable = app.getExecutionResultTable(result, ...
                {'unsavedFinalOutputs', 'unsavedLeafOutputs'});

            appendRows(viewerTable, "Viewer-compatible");
            appendRows(unsavedTable, "Not saved");

            T = table(category, saveFolder, stepTag, outputName, handling, fileName, filePath, ...
                'VariableNames', {'Category','SaveFolder','Step','Output','Handling','FileName','FilePath'});

            function appendRows(srcTable, categoryText)
                if ~istable(srcTable) || isempty(srcTable)
                    return
                end

                for iRow = 1:height(srcTable)
                    category(end+1,1) = categoryText; %#ok<AGROW>
                    saveFolder(end+1,1) = app.getFirstTableStringValue(srcTable, ...
                        {'SaveFolder','saveFolder'}, iRow, ""); %#ok<AGROW>
                    stepTag(end+1,1) = app.getFirstTableStringValue(srcTable, ...
                        {'StepTag','stepTag','Step','stepName'}, iRow, ""); %#ok<AGROW>
                    outputName(end+1,1) = app.getFirstTableStringValue(srcTable, ...
                        {'OutputName','outputName','Output'}, iRow, ""); %#ok<AGROW>

                    stateText = app.getFirstTableStringValue(srcTable, ...
                        {'ActualPersistence','actualPersistence','plannedPersistence','Persistence'}, iRow, "");
                    if strlength(strtrim(stateText)) > 0
                        stateText = app.formatPersistenceStateForDisplay(stateText);
                    end
                    handling(end+1,1) = stateText; %#ok<AGROW>

                    pathText = app.getFirstTableStringValue(srcTable, ...
                        {'FilePath','filePath','plannedFilePath','FullPath','fullPath'}, iRow, "");

                    nameText = app.getFirstTableStringValue(srcTable, ...
                        {'FileName','fileName','plannedFileName','OutputFile','outputFile'}, iRow, "");

                    % Some backend result tables only provide FilePath. For display,
                    % derive FileName from FilePath so the compact column remains useful.
                    if strlength(strtrim(nameText)) == 0 && strlength(strtrim(pathText)) > 0
                        try
                            [~, f, e] = fileparts(char(pathText));
                            nameText = string([f e]);
                        catch
                            nameText = "";
                        end
                    end

                    fileName(end+1,1) = nameText; %#ok<AGROW>
                    filePath(end+1,1) = pathText; %#ok<AGROW>
                end
            end
        end

        function T = buildPipelineIssuesTable(app, result)
            %BUILDPIPELINEISSUESTABLE Build a compact issues table from run result.

            category = strings(0,1);
            item = strings(0,1);
            detail = strings(0,1);

            if isstruct(result) && isfield(result, 'globalPipeLog') && istable(result.globalPipeLog)
                logTable = result.globalPipeLog;

                for iRow = 1:height(logTable)
                    statusText = app.getFirstTableStringValue(logTable, {'Status'}, iRow, "");
                    msgText = app.getFirstTableStringValue(logTable, {'Messages_short','Messages'}, iRow, "");
                    saveFolderText = app.getFirstTableStringValue(logTable, {'SaveFolder'}, iRow, "");

                    hasIssueStatus = strlength(statusText) > 0 && ...
                        ~(contains(lower(statusText), "completed") || contains(lower(statusText), "skipped"));

                    hasIssueMessage = strlength(strtrim(msgText)) > 0 && ...
                        ~strcmpi(strtrim(msgText), "No Errors") && ...
                        ~strcmpi(strtrim(msgText), "No Error");

                    if hasIssueStatus || hasIssueMessage
                        category(end+1,1) = "Folder run"; %#ok<AGROW>
                        item(end+1,1) = saveFolderText; %#ok<AGROW>
                        detail(end+1,1) = strtrim(statusText + " | " + msgText); %#ok<AGROW>
                    end
                end
            end

            unsavedTable = app.getExecutionResultTable(result, ...
                {'unsavedFinalOutputs', 'unsavedLeafOutputs'});
            if istable(unsavedTable) && ~isempty(unsavedTable)
                for iRow = 1:height(unsavedTable)
                    stepText = app.getFirstTableStringValue(unsavedTable, ...
                        {'StepTag','stepTag','Step','stepName'}, iRow, "<step>");
                    outputText = app.getFirstTableStringValue(unsavedTable, ...
                        {'OutputName','outputName','Output'}, iRow, "<output>");

                    category(end+1,1) = "Unsaved final output"; %#ok<AGROW>
                    item(end+1,1) = stepText + ":" + outputText; %#ok<AGROW>
                    detail(end+1,1) = "Final output was reported but not saved."; %#ok<AGROW>
                end
            end

            T = table(category, item, detail, 'VariableNames', {'Category','Item','Detail'});
        end

        function createPipelineSummaryTableTab(app, tabGroup, tabTitle, T, tableKind)
            %CREATEPIPELINESUMMARYTABLETAB Create one table tab in the summary viewer.

            tab = uitab(tabGroup);
            tab.Title = tabTitle;

            grid = uigridlayout(tab);
            grid.RowHeight = {'1x'};
            grid.ColumnWidth = {'1x'};
            grid.Padding = [8 8 8 8];

            if ~istable(T) || isempty(T)
                label = uilabel(grid);
                label.Layout.Row = 1;
                label.Layout.Column = 1;
                label.Text = sprintf('No %s information available.', lower(tabTitle));
                label.HorizontalAlignment = 'center';
                label.FontColor = [0.45 0.45 0.45];
                return
            end

            [displayTable, columnNames] = app.preparePipelineSummaryDisplayTable(T, tableKind);

            tbl = uitable(grid);
            tbl.Layout.Row = 1;
            tbl.Layout.Column = 1;
            tbl.Data = displayTable;
            tbl.ColumnName = columnNames;
            tbl.ColumnEditable = false(1, width(displayTable));
            tbl.RowName = 'numbered';
        end

        function [displayTable, columnNames] = preparePipelineSummaryDisplayTable(app, T, tableKind)
            %PREPAREPIPELINESUMMARYDISPLAYTABLE Convert run-result tables to UI-friendly tables.

            if nargin < 3
                tableKind = '';
            end

            displayTable = T;

            for iVar = 1:width(displayTable)
                varName = displayTable.Properties.VariableNames{iVar};
                values = displayTable.(varName);

                if iscell(values)
                    displayTable.(varName) = string(cellfun( ...
                        @(x) app.formatTableValueForDisplay(x), ...
                        values, ...
                        'UniformOutput', false));

                elseif isstruct(values)
                    displayTable.(varName) = arrayfun( ...
                        @(x) string(app.formatTableValueForDisplay(x)), ...
                        values);

                elseif isnumeric(values) || islogical(values) || isstring(values) || ...
                        isdatetime(values) || isduration(values) || iscategorical(values)
                    % These are supported by uitable/table display.
                else
                    displayTable.(varName) = arrayfun( ...
                        @(x) string(app.formatTableValueForDisplay(x)), ...
                        values);
                end
            end

            persistenceVars = {'Handling','ActualPersistence','actualPersistence','plannedPersistence','Persistence'};
            for iVar = 1:numel(persistenceVars)
                varName = persistenceVars{iVar};
                if ismember(varName, displayTable.Properties.VariableNames)
                    displayTable.(varName) = arrayfun( ...
                        @(x) string(app.formatPersistenceStateForDisplay(x)), ...
                        string(displayTable.(varName)));
                end
            end

            columnNames = app.getPipelineSummaryColumnNames(displayTable.Properties.VariableNames, tableKind);
        end

        function columnNames = getPipelineSummaryColumnNames(~, varNames, tableKind)
            %GETPIPELINESUMMARYCOLUMNNAMES Return user-facing summary table labels.

            if nargin < 3
                tableKind = '';
            end

            columnNames = cellstr(string(varNames));

            for iName = 1:numel(columnNames)
                nameLower = lower(columnNames{iName});

                switch nameLower
                    case 'savefolder'
                        columnNames{iName} = 'SaveFolder';
                    case 'rawfolder'
                        columnNames{iName} = 'RawFolder';
                    case 'status'
                        columnNames{iName} = 'Status';
                    case 'messages_short'
                        columnNames{iName} = 'Message';
                    case 'startedon'
                        columnNames{iName} = 'Started';
                    case 'finishedon'
                        columnNames{iName} = 'Finished';
                    case 'duration_s'
                        columnNames{iName} = 'Duration (s)';
                    case 'outputfiles'
                        columnNames{iName} = 'Output files';
                    case 'filepath'
                        columnNames{iName} = 'File path';
                    case 'filename'
                        columnNames{iName} = 'File name';
                    case 'istemporary'
                        columnNames{iName} = 'Temporary';
                    case 'viewcompatible'
                        columnNames{iName} = 'Viewer-compatible';
                    case 'step'
                        columnNames{iName} = 'Step';
                    case {'steptag','stepname'}
                        columnNames{iName} = 'Step';
                    case {'output','outputname'}
                        columnNames{iName} = 'Output';
                    case {'handling','actualpersistence','plannedpersistence','persistence'}
                        columnNames{iName} = 'Handling';
                    case 'category'
                        columnNames{iName} = 'Category';
                    case 'item'
                        columnNames{iName} = 'Item';
                    case 'detail'
                        columnNames{iName} = 'Detail';
                    case 'field'
                        columnNames{iName} = 'Field';
                    case 'value'
                        columnNames{iName} = 'Value';
                    otherwise
                        % Keep backend name if no user-facing alias is known.
                end
            end

            %#ok<NASGU> tableKind is kept for future table-specific aliases.
        end

        function appendExecutionResultLogSummary(app, result)
            %APPENDEXECUTIONRESULTLOGSUMMARY Append one compact run summary line.

            statusText = "<unknown>";
            durationText = "";
            nFolders = 0;
            nCompleted = 0;
            nSkipped = 0;
            nFailed = 0;
            nCreated = 0;

            if isstruct(result)
                if isfield(result, 'status')
                    statusText = string(result.status);
                end

                if isfield(result, 'duration_s') && isnumeric(result.duration_s) && isscalar(result.duration_s)
                    durationText = sprintf(', duration %.3f s', result.duration_s);
                end

                if isfield(result, 'globalPipeLog') && istable(result.globalPipeLog)
                    logTable = result.globalPipeLog;
                    nFolders = height(logTable);

                    if ismember('Status', logTable.Properties.VariableNames)
                        statusVals = lower(string(logTable.Status));
                        nCompleted = sum(contains(statusVals, "completed"));
                        nSkipped = sum(contains(statusVals, "skipped"));
                        nFailed = sum(contains(statusVals, "failed") | contains(statusVals, "error"));
                    end
                end

                if isfield(result, 'createdFiles') && istable(result.createdFiles)
                    nCreated = height(result.createdFiles);
                end
            end

            app.appendDiagnostic(sprintf(['Pipeline execution finished: %s%s. ' ...
                'Folders: %d total, %d completed, %d skipped, %d failed. Created files: %d.'], ...
                char(statusText), durationText, nFolders, nCompleted, nSkipped, nFailed, nCreated));
        end

        function appendValidationReportLogSummary(app, report, labelText)
            %APPENDVALIDATIONREPORTLOGSUMMARY Append compact validation summary to Activity Log.

            if nargin < 3 || isempty(labelText)
                labelText = 'Pipeline validation';
            end

            if isempty(report) || ~isstruct(report)
                app.appendDiagnostic(sprintf('%s: no validation report available.', labelText));
                return
            end

            isValidText = "<unknown>";
            nErrors = 0;
            nWarnings = 0;
            nBlocking = 0;
            nNonBlocking = 0;

            if isfield(report, 'isValid')
                isValidText = string(logical(report.isValid));
            end

            if isfield(report, 'errors') && ~isempty(report.errors)
                nErrors = numel(report.errors);
            end

            if isfield(report, 'warnings') && ~isempty(report.warnings)
                nWarnings = numel(report.warnings);
            end

            if isfield(report, 'nodeStatus') && istable(report.nodeStatus) && ...
                    ismember('status', report.nodeStatus.Properties.VariableNames)
                statusVals = string(report.nodeStatus.status);
                nBlocking = sum(statusVals == "invalid_blocking");
                nNonBlocking = sum(statusVals == "invalid_nonblocking");
            end

            app.appendDiagnostic(sprintf('%s: valid=%s, errors=%d, warnings=%d, blocking invalid steps=%d, non-blocking invalid steps=%d.', ...
                labelText, char(isValidText), nErrors, nWarnings, nBlocking, nNonBlocking));
        end

        function txt = getResultScalarString(app, result, fieldName, defaultValue)
            %GETRESULTSCALARSTRING Return one scalar result field as string.
            %
            % The execution result does not always duplicate backend execution
            % settings. If the result field is missing or empty, fall back to the
            % active PipelineManager property with the same name.

            txt = string(defaultValue);

            try
                if isstruct(result) && isfield(result, fieldName)
                    value = result.(fieldName);
                    if ~isempty(value)
                        if ischar(value) || (isstring(value) && isscalar(value)) || ...
                                (isnumeric(value) && isscalar(value)) || ...
                                (islogical(value) && isscalar(value))
                            txt = string(value);
                        else
                            txt = string(app.formatTableValueForDisplay(value));
                        end
                    end
                end
            catch
                txt = string(defaultValue);
            end

            if ~strcmp(txt, string(defaultValue)) && strlength(strtrim(txt)) > 0
                return
            end

            try
                if ~isempty(app.pm) && isprop(app.pm, fieldName)
                    value = app.pm.(fieldName);
                    if ~isempty(value)
                        if ischar(value) || (isstring(value) && isscalar(value)) || ...
                                (isnumeric(value) && isscalar(value)) || ...
                                (islogical(value) && isscalar(value))
                            txt = string(value);
                        else
                            txt = string(app.formatTableValueForDisplay(value));
                        end
                    end
                end
            catch
                txt = string(defaultValue);
            end
        end

        function value = getFirstTableStringValue(~, T, varNames, rowIdx, defaultValue)
            %GETFIRSTTABLESTRINGVALUE Return first available table value from a list.

            value = string(defaultValue);

            if ~istable(T) || rowIdx < 1 || rowIdx > height(T)
                return
            end

            varNames = cellstr(string(varNames));

            for iVar = 1:numel(varNames)
                varName = varNames{iVar};

                if ~ismember(varName, T.Properties.VariableNames)
                    continue
                end

                try
                    rawValue = T.(varName)(rowIdx);
                    if iscell(rawValue)
                        value = string(rawValue{1});
                    else
                        value = string(rawValue);
                    end
                    return
                catch
                    value = string(defaultValue);
                end
            end
        end

        function showFolderPipelineLogInteractive(app)
            %SHOWFOLDERPIPELINELOGINTERACTIVE View each SaveFolder pipeLog.mat table.
            %
            % This is a read-only modal viewer. It lists enabled SaveFolders and
            % renders the table stored as variable "pipeLog" in each folder's
            % pipeLog.mat file.

            if isempty(app.pm) && (isempty(app.FolderPairTable) || ~istable(app.FolderPairTable) || height(app.FolderPairTable) == 0)
                app.appendDiagnostic('Folder pipeline-log viewer requested without selected data folders.');
                app.setStatus('No data folders selected.');
                return
            end

            folderIndex = buildFolderIndexTable();

            if isempty(folderIndex) || height(folderIndex) == 0
                app.appendDiagnostic('Folder pipeline-log viewer found no enabled SaveFolders.');
                app.setStatus('No enabled SaveFolders.');
                return
            end

            mainPos = app.UIFigure.Position;
            figW = 1100;
            figH = 680;
            x = mainPos(1) + max(20, (mainPos(3) - figW) / 2);
            y = mainPos(2) + max(20, (mainPos(4) - figH) / 2);

            dlg = uifigure( ...
                'Name', 'Folder Pipeline Log', ...
                'Position', [x y figW figH], ...
                'Resize', 'on', ...
                'WindowStyle', 'modal');

            mainGrid = uigridlayout(dlg);
            mainGrid.RowHeight = {32, 160, 28, '1x', 36};
            mainGrid.ColumnWidth = {'1x'};
            mainGrid.Padding = [10 10 10 10];
            mainGrid.RowSpacing = 8;

            titleLabel = uilabel(mainGrid);
            titleLabel.Layout.Row = 1;
            titleLabel.Layout.Column = 1;
            titleLabel.FontWeight = 'bold';
            titleLabel.Text = sprintf('Folder pipeline logs: %d folder(s). Select a row to load pipeLog.mat.', height(folderIndex));

            folderTable = uitable(mainGrid);
            folderTable.Layout.Row = 2;
            folderTable.Layout.Column = 1;
            folderTable.Data = folderIndex;
            folderTable.ColumnName = {'Folder','SaveFolder','Has pipeLog','Status'};
            folderTable.ColumnEditable = false(1, width(folderIndex));
            folderTable.RowName = 'numbered';
            folderTable.CellSelectionCallback = @(~, event) selectFolderRow(event);

            try
                folderTable.ColumnWidth = {160, 'auto', 72, 170};
            catch
            end


            selectedLabel = uilabel(mainGrid);
            selectedLabel.Layout.Row = 3;
            selectedLabel.Layout.Column = 1;
            selectedLabel.Text = 'Selected folder: <none>';
            selectedLabel.FontColor = [0.35 0.35 0.35];

            logTable = uitable(mainGrid);
            logTable.Layout.Row = 4;
            logTable.Layout.Column = 1;
            logTable.Data = table();
            logTable.ColumnEditable = false;

            bottomGrid = uigridlayout(mainGrid);
            bottomGrid.Layout.Row = 5;
            bottomGrid.Layout.Column = 1;
            bottomGrid.RowHeight = {'1x'};
            bottomGrid.ColumnWidth = {'1x', 90, 90};
            bottomGrid.Padding = [0 0 0 0];
            bottomGrid.ColumnSpacing = 8;

            infoLabel = uilabel(bottomGrid);
            infoLabel.Layout.Row = 1;
            infoLabel.Layout.Column = 1;
            infoLabel.Text = '';
            infoLabel.FontColor = [0.35 0.35 0.35];

            refreshButton = uibutton(bottomGrid, 'push');
            refreshButton.Text = 'Refresh';
            refreshButton.Tooltip = 'Refresh the folder list and reload the selected folder log.';
            refreshButton.Layout.Row = 1;
            refreshButton.Layout.Column = 2;
            refreshButton.ButtonPushedFcn = @(~,~) refreshViewer();

            closeButton = uibutton(bottomGrid, 'push');
            closeButton.Text = 'Close';
            closeButton.Layout.Row = 1;
            closeButton.Layout.Column = 3;
            closeButton.ButtonPushedFcn = @(~,~) closeViewer();

            selectedRow = find(folderIndex.HasPipeLog, 1, 'first');
            if isempty(selectedRow)
                selectedRow = 1;
            end

            loadSelectedFolderLog(selectedRow);

            app.appendDiagnostic('Opened folder pipeline-log viewer.');
            app.setStatus('Folder pipeline-log viewer opened.');

            uiwait(dlg);

            if isvalid(dlg)
                delete(dlg);
            end

            % -------------------------------------------------------------
            % Nested helpers
            % -------------------------------------------------------------

            function T = buildFolderIndexTable()
                saveFolders = strings(0,1);

                try
                    if ~isempty(app.FolderPairTable) && istable(app.FolderPairTable) && height(app.FolderPairTable) > 0
                        Tpairs = app.FolderPairTable;
                        if ismember('Use', Tpairs.Properties.VariableNames)
                            Tpairs = Tpairs(logical(Tpairs.Use), :);
                        end

                        if ismember('SaveFolder', Tpairs.Properties.VariableNames)
                            saveFolders = string(Tpairs.SaveFolder(:));
                        end
                    end
                catch
                    saveFolders = strings(0,1);
                end

                if isempty(saveFolders)
                    try
                        if ~isempty(app.pm) && isprop(app.pm, 'SaveFolderList') && ~isempty(app.pm.SaveFolderList)
                            saveFolders = string(app.pm.SaveFolderList(:));
                        end
                    catch
                        saveFolders = strings(0,1);
                    end
                end

                saveFolders = saveFolders(strlength(strtrim(saveFolders)) > 0);
                saveFolders = unique(saveFolders(:), 'stable');

                folderName = strings(numel(saveFolders), 1);
                hasPipeLog = false(numel(saveFolders), 1);
                status = strings(numel(saveFolders), 1);


                for iFolder = 1:numel(saveFolders)
                    sf = saveFolders(iFolder);

                    try
                        [~, folderName(iFolder)] = fileparts(char(sf));
                    catch
                        folderName(iFolder) = "";
                    end

                    if ~isfolder(char(sf))
                        status(iFolder) = "SaveFolder missing";
                        hasPipeLog(iFolder) = false;
                        continue
                    end

                    logFile = fullfile(char(sf), 'pipeLog.mat');
                    hasPipeLog(iFolder) = isfile(logFile);

                    if hasPipeLog(iFolder)
                        status(iFolder) = "Ready";
                    else
                        status(iFolder) = "pipeLog.mat not found";
                    end
                end

                T = table(folderName, saveFolders, hasPipeLog, status, ...
                    'VariableNames', {'Folder','SaveFolder','HasPipeLog','Status'});
            end

            function selectFolderRow(event)
                if isempty(event.Indices)
                    return
                end

                rowIdx = event.Indices(1,1);
                if rowIdx < 1 || rowIdx > height(folderTable.Data)
                    return
                end

                selectedRow = rowIdx;
                loadSelectedFolderLog(selectedRow);
            end

            function refreshViewer()
                oldRow = selectedRow;

                folderIndex = buildFolderIndexTable();
                folderTable.Data = folderIndex;

                if isempty(folderIndex) || height(folderIndex) == 0
                    selectedRow = [];
                    selectedLabel.Text = 'Selected folder: <none>';
                    logTable.Data = table();
                    logTable.ColumnName = {};
                    infoLabel.Text = 'No enabled SaveFolders.';
                    return
                end

                selectedRow = min(max(1, oldRow), height(folderIndex));
                loadSelectedFolderLog(selectedRow);
            end

            function loadSelectedFolderLog(rowIdx)
                if isempty(folderIndex) || rowIdx < 1 || rowIdx > height(folderIndex)
                    selectedLabel.Text = 'Selected folder: <none>';
                    logTable.Data = table();
                    logTable.ColumnName = {};
                    infoLabel.Text = 'No folder selected.';
                    return
                end

                saveFolder = string(folderIndex.SaveFolder(rowIdx));
                folderNameLocal = string(folderIndex.Folder(rowIdx));
                selectedLabel.Text = sprintf('Selected folder: %s', char(saveFolder));

                if ~isfolder(char(saveFolder))
                    showMessageTable("SaveFolder is missing or inaccessible.");
                    infoLabel.Text = 'SaveFolder missing.';
                    return
                end

                logFile = fullfile(char(saveFolder), 'pipeLog.mat');

                if ~isfile(logFile)
                    showMessageTable("pipeLog.mat not found in selected SaveFolder.");
                    infoLabel.Text = 'pipeLog.mat not found.';
                    return
                end

                try
                    S = load(logFile, '-mat', 'pipeLog');
                catch ME
                    showMessageTable("Failed to load pipeLog.mat: " + string(ME.message));
                    infoLabel.Text = 'Failed to load pipeLog.mat.';
                    return
                end

                if ~isfield(S, 'pipeLog')
                    showMessageTable("Variable pipeLog was not found in pipeLog.mat.");
                    infoLabel.Text = 'pipeLog variable missing.';
                    return
                end

                if ~istable(S.pipeLog)
                    showMessageTable("Variable pipeLog exists but is not a table.");
                    infoLabel.Text = 'pipeLog is not a table.';
                    return
                end

                displayTable = prepareDisplayTable(S.pipeLog);

                logTable.Data = displayTable;
                logTable.ColumnName = displayTable.Properties.VariableNames;
                logTable.ColumnEditable = false(1, width(displayTable));
                logTable.RowName = 'numbered';

                try
                    sortableVars = {'Status', 'StartedOn', 'FinishedOn'};
                    logTable.ColumnSortable = ismember(displayTable.Properties.VariableNames, sortableVars);
                catch
                    % ColumnSortable is cosmetic and version-dependent.
                end


                if isempty(displayTable)
                    infoLabel.Text = sprintf('Loaded empty pipeLog table for folder: %s', char(folderNameLocal));
                else
                    infoLabel.Text = sprintf('Loaded pipeLog table: %d row(s), %d column(s).', ...
                        height(displayTable), width(displayTable));
                end

                app.appendDiagnostic(sprintf('Loaded pipeLog.mat for folder: %s', char(folderNameLocal)));
            end

            function displayTable = prepareDisplayTable(T)
                displayTable = T;

                for iVar = 1:width(displayTable)
                    varName = displayTable.Properties.VariableNames{iVar};
                    values = displayTable.(varName);

                    if iscell(values)
                        displayTable.(varName) = string(cellfun( ...
                            @(x) formatOneValue(x), ...
                            values, ...
                            'UniformOutput', false));

                    elseif isstruct(values)
                        displayTable.(varName) = arrayfun( ...
                            @(x) string(formatOneValue(x)), ...
                            values);

                    elseif isnumeric(values) || islogical(values) || isstring(values) || ...
                            isdatetime(values) || isduration(values) || iscategorical(values)
                        % uitable supports these directly.

                    else
                        try
                            displayTable.(varName) = string(values);
                        catch
                            displayTable.(varName) = repmat("<unsupported>", height(displayTable), 1);
                        end
                    end
                end
            end

            function txt = formatOneValue(valueIn)
                if isempty(valueIn)
                    txt = "";
                    return
                end

                if ischar(valueIn) || (isstring(valueIn) && isscalar(valueIn))
                    txt = char(string(valueIn));
                    return
                end

                if iscell(valueIn)
                    try
                        asText = string(valueIn(:));
                        asText = asText(strlength(strtrim(asText)) > 0);
                        txt = char(strjoin(asText, ', '));
                    catch
                        txt = sprintf('[%s cell]', strjoin(string(size(valueIn)), 'x'));
                    end
                    return
                end

                if isnumeric(valueIn) || islogical(valueIn)
                    if isscalar(valueIn) || (isvector(valueIn) && numel(valueIn) <= 8)
                        txt = mat2str(valueIn);
                    else
                        txt = sprintf('[%s %s]', strjoin(string(size(valueIn)), 'x'), class(valueIn));
                    end
                    return
                end

                if isdatetime(valueIn) || isduration(valueIn) || iscategorical(valueIn)
                    txt = char(string(valueIn));
                    return
                end

                if isstruct(valueIn)
                    if isscalar(valueIn)
                        txt = sprintf('[struct: %s]', strjoin(string(fieldnames(valueIn)), ', '));
                    else
                        txt = sprintf('[%s struct]', strjoin(string(size(valueIn)), 'x'));
                    end
                    return
                end

                txt = class(valueIn);
            end

            function showMessageTable(msg)
                message = string(msg);
                logTable.Data = table(message, 'VariableNames', {'Message'});
                logTable.ColumnName = {'Message'};
                logTable.ColumnEditable = false;
                logTable.RowName = 'numbered';

                try
                    logTable.ColumnSortable = false;
                catch
                end

            end

            function closeViewer()
                if isvalid(dlg)
                    uiresume(dlg);
                end
            end
        end

        function cacheRunButtonDefaultStyle(app)
            %CACHERUNBUTTONDEFAULTSTYLE Store the original Run button appearance.

            try
                if isempty(app.RunButtonDefaultBackgroundColor)
                    app.RunButtonDefaultBackgroundColor = app.RunButton.BackgroundColor;
                end
            catch
                app.RunButtonDefaultBackgroundColor = [];
            end

            try
                if isempty(app.RunButtonDefaultFontColor)
                    app.RunButtonDefaultFontColor = app.RunButton.FontColor;
                end
            catch
                app.RunButtonDefaultFontColor = [];
            end
        end

        function prepareGuiForPipelineRun(app)
            %PREPAREGUIFORPIPELINERUN Put the GUI in running state.

            app.bPipelineRunning = true;
            app.bCancelRequested = false;

            app.CurrentFolderIndex = NaN;
            app.NumExecutionFolders = NaN;
            app.CurrentStepIndex = NaN;
            app.NumExecutionSteps = NaN;
            app.CurrentExecutionFolder = "";
            app.CurrentExecutionStep = "";
            app.resetRuntimeExecutionState(true);

            app.resetExecutionProgressPanel();
            app.setRunButtonState('running');
            app.refreshRuntimeHtmlViews();

            controlsToDisable = { ...
                app.SelectDataFoldersButton, ...
                app.LoadButton, ...
                app.SaveButton, ...
                app.AutoSaveFinalOutputsDropDown, ...
                app.RAMModeDropDown, ...
                app.RAMSafetyPolicyDropDown, ...
                app.ReuseExistingFilesCheckBox, ...
                app.ClearPipelineButton, ...
                app.AddButton, ...
                app.AddasNewBranchButton, ...
                app.SearchFunctionEditField, ...
                app.FunctionTree};

            for iCtrl = 1:numel(controlsToDisable)
                try
                    controlsToDisable{iCtrl}.Enable = 'off';
                catch
                end
            end

            menuItemsToDisable = { ...
                app.LoadPipelineMenu, ...
                app.SavePipelineMenu, ...
                app.GenerateScriptMenu, ...
                app.SelectDataFoldersMenu, ...
                app.RunPipelineMenu, ...
                app.ReuseExistingFilesMenu, ...
                app.ExecutionSettingsMenu, ...
                app.AdvancedRAMSettingsMenu, ...
                app.ShowExecutionPlanMenu, ...
                app.RemoveInvalidStepsMenu, ...
                app.ShowLatestRunSummaryMenu, ...
                app.ViewFolderPipelineLogMenu, ...
                app.ExportErrorLogMenu};

            for iMenu = 1:numel(menuItemsToDisable)
                try
                    menuItemsToDisable{iMenu}.Enable = 'off';
                catch
                end
            end

            app.updateRamSafetyControlState();
        end

        function restoreGuiAfterPipelineRun(app)
            %RESTOREGUIAFTERPIPELINERUN Restore GUI controls after executePipeline returns.

            app.bPipelineRunning = false;
            app.bCancelRequested = false;

            app.setRunButtonState('idle');

            try
                app.setPipelineControlsAvailable(~isempty(app.pm));
            catch
                try
                    app.RunButton.Enable = 'on';
                catch
                end
            end

            app.refreshRuntimeHtmlViews();

            app.updateLatestRunSummaryMenuState();
        end

        function setRunButtonState(app, stateName)
            %SETRUNBUTTONSTATE Update Run button text, color, and enabled state.

            stateName = lower(char(string(stateName)));
            app.cacheRunButtonDefaultStyle();

            switch stateName
                case 'idle'
                    app.RunButton.Text = 'Run';
                    app.RunButton.Enable = 'on';

                    if ~isempty(app.RunButtonDefaultBackgroundColor)
                        try
                            app.RunButton.BackgroundColor = app.RunButtonDefaultBackgroundColor;
                        catch
                        end
                    end

                    if ~isempty(app.RunButtonDefaultFontColor)
                        try
                            app.RunButton.FontColor = app.RunButtonDefaultFontColor;
                        catch
                        end
                    end

                    try
                        app.RunButton.Tooltip = ...
                            {'Execute the current pipeline using the selected folder pairs and execution settings.'};
                    catch
                    end

                case 'running'
                    app.RunButton.Text = 'Stop';
                    app.RunButton.Enable = 'on';

                    try
                        app.RunButton.BackgroundColor = [0.85 0.20 0.20];
                        app.RunButton.FontColor = [1 1 1];
                    catch
                    end

                    try
                        app.RunButton.Tooltip = ...
                            {'Request pipeline stop. The current function will finish before execution stops.'};
                    catch
                    end

                case 'stopping'
                    app.RunButton.Text = 'Stopping...';
                    app.RunButton.Enable = 'off';

                    try
                        app.RunButton.BackgroundColor = [0.60 0.20 0.20];
                        app.RunButton.FontColor = [1 1 1];
                    catch
                    end

                    try
                        app.RunButton.Tooltip = ...
                            {'Stop requested. Waiting for the current function to finish.'};
                    catch
                    end

                otherwise
                    app.RunButton.Text = 'Run';
                    app.RunButton.Enable = 'on';
            end
        end

        function requestPipelineStop(app)
            %REQUESTPIPELINESTOP Request cooperative cancellation.

            if ~app.bPipelineRunning
                return
            end

            if app.bCancelRequested
                return
            end

            app.bCancelRequested = true;
            app.setRunButtonState('stopping');

            app.setStatus('Stop requested. Waiting for the current function to finish...');
            app.appendDiagnostic(['Stop requested. The current function will finish before ' ...
                'PipelineManager stops before the next folder, step, or file load.']);

            app.setProgressLabelsForStopping();
            drawnow;
        end

        function tf = isPipelineCancelRequested(app)
            %ISPIPELINECANCELREQUESTED Return true when the user requested stop.

            tf = logical(app.bCancelRequested);
        end

        function onPipelineProgress(app, evt)
            %ONPIPELINEPROGRESS Receive one PipelineManager progress event.

            if isempty(evt) || ~isstruct(evt)
                return
            end

            app.updatePipelineProgressPanel(evt);

            evtType = lower(char(app.getProgressEventField(evt, 'type', "")));
            app.updateRuntimeExecutionStateFromProgressEvent(evt);

            switch evtType
                case 'runstarted'
                    app.appendDiagnostic('Pipeline execution started.');

                case 'folderstarted'
                    folderName = app.getProgressEventDisplayFolder(evt);
                    app.appendDiagnostic(sprintf('Started folder: %s', folderName));

                case 'folderfinished'
                    folderName = app.getProgressEventDisplayFolder(evt);
                    app.appendDiagnostic(sprintf('Finished folder: %s', folderName));

                case 'folderskipped'
                    folderName = app.getProgressEventDisplayFolder(evt);
                    msg = char(app.getProgressEventField(evt, 'message', "folder skipped"));
                    app.appendDiagnostic(sprintf('Skipped folder: %s (%s)', folderName, msg));

                case 'stepfailed'
                    stepTag = char(app.getProgressEventField(evt, 'stepTag', "<step>"));
                    msg = char(app.getProgressEventField(evt, 'message', "step failed"));
                    app.appendDiagnostic(sprintf('Step failed: %s (%s)', stepTag, msg));

                case 'cancelrequested'
                    app.setProgressLabelsForStopping();

                case 'runcancelled'
                    app.appendDiagnostic('Pipeline execution cancelled.');
                    app.setStatus('Pipeline execution cancelled.');

                case 'runfinished'
                    app.appendDiagnostic('Pipeline execution finished.');
            end

            app.refreshRuntimeHtmlViews();

            drawnow limitrate
        end

        function updatePipelineProgressPanel(app, evt)
            %UPDATEPIPELINEPROGRESSPANEL Update compact folder/step progress bars.

            evtType = lower(char(app.getProgressEventField(evt, 'type', "")));

            folderIndex = app.getProgressEventNumericField(evt, 'folderIndex', app.CurrentFolderIndex);
            numFolders = app.getProgressEventNumericField(evt, 'numFolders', app.NumExecutionFolders);
            stepIndex = app.getProgressEventNumericField(evt, 'stepIndex', app.CurrentStepIndex);
            numSteps = app.getProgressEventNumericField(evt, 'numSteps', app.NumExecutionSteps);

            if ~isnan(folderIndex); app.CurrentFolderIndex = folderIndex; end
            if ~isnan(numFolders); app.NumExecutionFolders = numFolders; end
            if ~isnan(stepIndex); app.CurrentStepIndex = stepIndex; end
            if ~isnan(numSteps); app.NumExecutionSteps = numSteps; end

            switch evtType
                case 'runstarted'
                    nF = app.getProgressEventNumericField(evt, 'numFolders', NaN);
                    if isnan(nF)
                        nF = app.NumExecutionFolders;
                    end
                    app.NumExecutionFolders = nF;
                    app.CurrentFolderIndex = 0;
                    app.CurrentStepIndex = 0;
                    app.NumExecutionSteps = NaN;
                    app.setFolderProgressTextAndFraction(sprintf('Folders: 0/%s - running', app.formatProgressCount(nF)), 0);
                    app.setStepProgressTextAndFraction('Steps: waiting', 0);

                case 'folderstarted'
                    folderName = app.getProgressEventDisplayFolder(evt);
                    app.CurrentExecutionFolder = string(folderName);
                    frac = app.getCompletedFraction(app.CurrentFolderIndex - 1, app.NumExecutionFolders);
                    app.setFolderProgressTextAndFraction( ...
                        sprintf('Folders: %s/%s - %s', ...
                        app.formatProgressCount(app.CurrentFolderIndex), ...
                        app.formatProgressCount(app.NumExecutionFolders), ...
                        folderName), ...
                        frac);
                    app.setStepProgressTextAndFraction('Steps: waiting', 0);

                case 'folderfinished'
                    folderName = app.getProgressEventDisplayFolder(evt);
                    frac = app.getCompletedFraction(app.CurrentFolderIndex, app.NumExecutionFolders);
                    app.setFolderProgressTextAndFraction( ...
                        sprintf('Folders: %s/%s - finished %s', ...
                        app.formatProgressCount(app.CurrentFolderIndex), ...
                        app.formatProgressCount(app.NumExecutionFolders), ...
                        folderName), ...
                        frac);

                case 'folderskipped'
                    folderName = app.getProgressEventDisplayFolder(evt);
                    frac = app.getCompletedFraction(app.CurrentFolderIndex, app.NumExecutionFolders);
                    app.setFolderProgressTextAndFraction( ...
                        sprintf('Folders: %s/%s - skipped %s', ...
                        app.formatProgressCount(app.CurrentFolderIndex), ...
                        app.formatProgressCount(app.NumExecutionFolders), ...
                        folderName), ...
                        frac);

                case 'stepstarted'
                    stepTag = char(app.getProgressEventField(evt, 'stepTag', "<step>"));
                    app.CurrentExecutionStep = string(stepTag);
                    frac = app.getCompletedFraction(app.CurrentStepIndex - 1, app.NumExecutionSteps);
                    app.setStepProgressTextAndFraction( ...
                        sprintf('Steps: %s/%s - %s', ...
                        app.formatProgressCount(app.CurrentStepIndex), ...
                        app.formatProgressCount(app.NumExecutionSteps), ...
                        stepTag), ...
                        frac);
                    app.setStatus(sprintf('Running: %s', stepTag));

                case {'stepfinished','stepskipped'}
                    stepTag = char(app.getProgressEventField(evt, 'stepTag', "<step>"));
                    frac = app.getCompletedFraction(app.CurrentStepIndex, app.NumExecutionSteps);
                    app.setStepProgressTextAndFraction( ...
                        sprintf('Steps: %s/%s - finished %s', ...
                        app.formatProgressCount(app.CurrentStepIndex), ...
                        app.formatProgressCount(app.NumExecutionSteps), ...
                        stepTag), ...
                        frac);

                case 'stepfailed'
                    stepTag = char(app.getProgressEventField(evt, 'stepTag', "<step>"));
                    frac = app.getCompletedFraction(app.CurrentStepIndex, app.NumExecutionSteps);
                    app.setStepProgressTextAndFraction( ...
                        sprintf('Steps: %s/%s - failed %s', ...
                        app.formatProgressCount(app.CurrentStepIndex), ...
                        app.formatProgressCount(app.NumExecutionSteps), ...
                        stepTag), ...
                        frac);

                case 'cancelrequested'
                    app.setProgressLabelsForStopping();

                case 'runcancelled'
                    app.setProgressLabelsForStopping();
                    app.setStatus('Pipeline execution cancelled.');

                case 'runfinished'
                    app.setFolderProgressTextAndFraction('Folders: completed', 1);
                    app.setStepProgressTextAndFraction('Steps: completed', 1);
            end
        end

        function resetExecutionProgressPanel(app)
            %RESETEXECUTIONPROGRESSPANEL Reset compact progress display.

            app.CurrentFolderIndex = NaN;
            app.NumExecutionFolders = NaN;
            app.CurrentStepIndex = NaN;
            app.NumExecutionSteps = NaN;
            app.CurrentExecutionFolder = "";
            app.CurrentExecutionStep = "";

            app.ensureExecutionProgressHtmlBars();

            app.setFolderProgressTextAndFraction('Folders: ready', 0);
            app.setStepProgressTextAndFraction('Steps: ready', 0);
        end

        function setProgressLabelsForStopping(app)
            %SETPROGRESSLABELSFORSTOPPING Show cooperative stop status.

            if strlength(app.CurrentExecutionFolder) > 0
                folderText = sprintf('Folders: %s/%s - stopping after current function', ...
                    app.formatProgressCount(app.CurrentFolderIndex), ...
                    app.formatProgressCount(app.NumExecutionFolders));
            else
                folderText = 'Folders: stopping after current function';
            end

            if strlength(app.CurrentExecutionStep) > 0
                stepText = sprintf('Steps: %s/%s - current function must finish before stopping', ...
                    app.formatProgressCount(app.CurrentStepIndex), ...
                    app.formatProgressCount(app.NumExecutionSteps));
            else
                stepText = 'Steps: current function must finish before stopping';
            end

            app.setFolderProgressTextAndFraction(folderText, ...
                app.getCompletedFraction(app.CurrentFolderIndex - 1, app.NumExecutionFolders));
            app.setStepProgressTextAndFraction(stepText, ...
                app.getCompletedFraction(app.CurrentStepIndex - 1, app.NumExecutionSteps));
        end

        function setFolderProgressTextAndFraction(app, txt, frac)
            %SETFOLDERPROGRESSTEXTANDFRACTION Update folder progress bar.

            app.ensureExecutionProgressHtmlBars();

            app.setProgressHtmlBar( ...
                app.FolderProgressHTML, ...
                app.FolderProgressOuterPanel, ...
                txt, ...
                frac, ...
                '#5A99E6');

            % Keep the legacy label synchronized in case the HTML component is
            % unavailable in a deployed environment.
            try
                app.FolderProgressLabel.Text = char(string(txt));
            catch
            end
        end

        function setStepProgressTextAndFraction(app, txt, frac)
            %SETSTEPPROGRESSTEXTANDFRACTION Update step progress bar.

            app.ensureExecutionProgressHtmlBars();

            app.setProgressHtmlBar( ...
                app.StepProgressHTML, ...
                app.StepProgressOuterPanel, ...
                txt, ...
                frac, ...
                '#59BF73');

            % Keep the legacy label synchronized in case the HTML component is
            % unavailable in a deployed environment.
            try
                app.StepProgressLabel.Text = char(string(txt));
            catch
            end
        end

        function frac = getCompletedFraction(~, indexValue, totalValue)
            %GETCOMPLETEDFRACTION Return bounded completed-unit fraction.

            if isempty(indexValue) || isempty(totalValue) || isnan(indexValue) || isnan(totalValue) || totalValue <= 0
                frac = 0;
                return
            end

            frac = max(0, min(1, double(indexValue) ./ double(totalValue)));
        end

        function txt = formatProgressCount(~, value)
            %FORMATPROGRESSCOUNT Return compact progress count text.

            if isempty(value) || isnan(value)
                txt = '?';
            else
                txt = sprintf('%d', round(double(value)));
            end
        end

        function value = getProgressEventField(~, evt, fieldName, defaultValue)
            %GETPROGRESSEVENTFIELD Read one progress event field safely.

            value = string(defaultValue);

            try
                if isstruct(evt) && isfield(evt, fieldName) && ~isempty(evt.(fieldName))
                    raw = evt.(fieldName);
                    if ischar(raw) || (isstring(raw) && isscalar(raw)) || ...
                            (isnumeric(raw) && isscalar(raw)) || ...
                            (islogical(raw) && isscalar(raw))
                        value = string(raw);
                    else
                        value = string(defaultValue);
                    end
                end
            catch
                value = string(defaultValue);
            end
        end

        function value = getProgressEventNumericField(~, evt, fieldName, defaultValue)
            %GETPROGRESSEVENTNUMERICFIELD Read one numeric progress event field.

            value = defaultValue;

            try
                if isstruct(evt) && isfield(evt, fieldName) && ~isempty(evt.(fieldName))
                    raw = evt.(fieldName);
                    if isnumeric(raw) && isscalar(raw) && isfinite(raw)
                        value = double(raw);
                    elseif isstring(raw) || ischar(raw)
                        tmp = str2double(char(string(raw)));
                        if isfinite(tmp)
                            value = tmp;
                        end
                    end
                end
            catch
                value = defaultValue;
            end
        end

        function folderName = getProgressEventDisplayFolder(app, evt)
            %GETPROGRESSEVENTDISPLAYFOLDER Return compact folder name from event.

            saveFolder = char(app.getProgressEventField(evt, 'saveFolder', ""));
            if isempty(strtrim(saveFolder))
                folderName = '<folder>';
                return
            end

            try
                [~, folderName] = fileparts(saveFolder);
                if isempty(folderName)
                    folderName = saveFolder;
                end
            catch
                folderName = saveFolder;
            end
        end

        function ensureExecutionProgressHtmlBars(app)
            %ENSUREEXECUTIONPROGRESSHTMLBARS Create HTML progress bars if needed.
            %
            % The original panel+label overlay is unreliable because uilabel is
            % not transparent and can hide the fill panel. The HTML bars keep the
            % text centered while showing a true filled background.

            try
                % Hide legacy overlay elements. Keep them in the app for Designer
                % compatibility, but do not use them for rendering.
                app.FolderProgressFillPanel.Visible = 'off';
                app.FolderProgressLabel.Visible = 'off';
            catch
            end

            try
                app.StepProgressFillPanel.Visible = 'off';
                app.StepProgressLabel.Visible = 'off';
            catch
            end

            try
                if isempty(app.FolderProgressHTML) || ~isvalid(app.FolderProgressHTML)
                    app.FolderProgressHTML = uihtml(app.FolderProgressOuterPanel);
                    app.FolderProgressHTML.Position = app.getInnerPanelPosition(app.FolderProgressOuterPanel);
                end
            catch
                app.FolderProgressHTML = [];
            end

            try
                if isempty(app.StepProgressHTML) || ~isvalid(app.StepProgressHTML)
                    app.StepProgressHTML = uihtml(app.StepProgressOuterPanel);
                    app.StepProgressHTML.Position = app.getInnerPanelPosition(app.StepProgressOuterPanel);
                end
            catch
                app.StepProgressHTML = [];
            end
        end

        function setProgressHtmlBar(app, htmlObj, outerPanel, txt, frac, fillColorHex)
            %SETPROGRESSHTMLBAR Render one compact progress bar with centered text.

            try
                if isempty(htmlObj) || ~isvalid(htmlObj)
                    return
                end

                frac = max(0, min(1, double(frac)));
                if isnan(frac)
                    frac = 0;
                end

                pct = round(100 * frac);
                txt = app.escapeHtmlText(char(string(txt)));

                htmlObj.Position = app.getInnerPanelPosition(outerPanel);

                htmlObj.HTMLSource = sprintf([ ...
                    '<!DOCTYPE html>' ...
                    '<html><head><meta charset="utf-8">' ...
                    '<style>' ...
                    'html,body{margin:0;padding:0;width:100%%;height:100%%;overflow:hidden;' ...
                    'font-family:Arial,Helvetica,sans-serif;font-size:10px;}' ...
                    '.bar{position:relative;width:100%%;height:100%%;' ...
                    'box-sizing:border-box;border:1px solid #B8B8B8;background:#EEEEEE;}' ...
                    '.fill{position:absolute;left:0;top:0;height:100%%;width:%d%%;' ...
                    'background:%s;}' ...
                    '.text{position:absolute;left:0;top:0;width:100%%;height:100%%;' ...
                    'display:flex;align-items:center;justify-content:center;' ...
                    'white-space:nowrap;overflow:hidden;text-overflow:ellipsis;' ...
                    'color:#111111;font-weight:normal;}' ...
                    '</style></head><body>' ...
                    '<div class="bar"><div class="fill"></div><div class="text">%s</div></div>' ...
                    '</body></html>'], ...
                    pct, fillColorHex, txt);

            catch
            end
        end

        function pos = getInnerPanelPosition(~, outerPanel)
            %GETINNERPANELPOSITION Return a child position that fills one panel.

            try
                p = outerPanel.Position;
                w = max(1, p(3));
                h = max(1, p(4));
                pos = [1 1 w h];
            catch
                pos = [1 1 100 12];
            end
        end

        function txt = escapeHtmlText(~, txt)
            %ESCAPEHTMLTEXT Escape text for inline HTML display.

            txt = char(string(txt));
            txt = strrep(txt, '&', '&amp;');
            txt = strrep(txt, '<', '&lt;');
            txt = strrep(txt, '>', '&gt;');
            txt = strrep(txt, '"', '&quot;');
            txt = strrep(txt, '''', '&#39;');
        end

        function resetRuntimeExecutionState(app, clearFinalState)
            %RESETRUNTIMEEXECUTIONSTATE Reset live execution state sent to HTML views.

            if nargin < 2 || isempty(clearFinalState)
                clearFinalState = false;
            end

            app.RuntimeCurrentStepID = NaN;
            app.RuntimeLastEventType = "";

            if clearFinalState
                app.RuntimeCompletedStepIDs = [];
                app.RuntimeSkippedStepIDs = [];
                app.RuntimeFailedStepIDs = [];
            end
        end

        function updateRuntimeExecutionStateFromProgressEvent(app, evt)
            %UPDATERUNTIMEEXECUTIONSTATEFROMPROGRESSEVENT Track step state for HTML views.

            eventType = lower(char(app.getProgressEventField(evt, 'type', "")));
            app.RuntimeLastEventType = string(eventType);

            switch eventType
                case 'runstarted'
                    app.resetRuntimeExecutionState(true);

                case 'stepstarted'
                    stepID = app.getNumericProgressEventField(evt, 'stepID', NaN);
                    app.RuntimeCurrentStepID = stepID;

                    % A re-run can visit a step that was present in stale runtime state.
                    app.RuntimeCompletedStepIDs = app.removeIDFromList(app.RuntimeCompletedStepIDs, stepID);
                    app.RuntimeSkippedStepIDs = app.removeIDFromList(app.RuntimeSkippedStepIDs, stepID);
                    app.RuntimeFailedStepIDs = app.removeIDFromList(app.RuntimeFailedStepIDs, stepID);

                case 'stepfinished'
                    stepID = app.getNumericProgressEventField(evt, 'stepID', app.RuntimeCurrentStepID);
                    app.RuntimeCompletedStepIDs = app.appendUniqueID(app.RuntimeCompletedStepIDs, stepID);
                    if isequaln(app.RuntimeCurrentStepID, stepID)
                        app.RuntimeCurrentStepID = NaN;
                    end

                case 'stepskipped'
                    stepID = app.getNumericProgressEventField(evt, 'stepID', app.RuntimeCurrentStepID);
                    app.RuntimeSkippedStepIDs = app.appendUniqueID(app.RuntimeSkippedStepIDs, stepID);
                    if isequaln(app.RuntimeCurrentStepID, stepID)
                        app.RuntimeCurrentStepID = NaN;
                    end

                case 'stepfailed'
                    stepID = app.getNumericProgressEventField(evt, 'stepID', app.RuntimeCurrentStepID);
                    app.RuntimeFailedStepIDs = app.appendUniqueID(app.RuntimeFailedStepIDs, stepID);
                    if isequaln(app.RuntimeCurrentStepID, stepID)
                        app.RuntimeCurrentStepID = NaN;
                    end

                case {'runcancelled','runfinished','runfailed'}
                    app.RuntimeCurrentStepID = NaN;
            end
        end

        function payload = buildHtmlRuntimePayload(app)
            %BUILDHTMLRUNTIMEPAYLOAD Build runtime state shared by both HTML views.

            payload = struct();
            payload.isRunning = logical(app.bPipelineRunning);
            payload.isStopping = logical(app.bCancelRequested);
            payload.interactionDisabled = logical(app.bPipelineRunning);
            payload.currentStepID = app.nanToEmptyNumeric(app.RuntimeCurrentStepID);
            payload.completedStepIDs = app.numericVectorForPayload(app.RuntimeCompletedStepIDs);
            payload.skippedStepIDs = app.numericVectorForPayload(app.RuntimeSkippedStepIDs);
            payload.failedStepIDs = app.numericVectorForPayload(app.RuntimeFailedStepIDs);
            payload.lastEventType = char(string(app.RuntimeLastEventType));
        end

        function state = getRuntimeStepState(app, stepID)
            %GETRUNTIMESTEPSTATE Return one display state for GraphHTML/ExecutionHTML.

            state = 'idle';

            if isempty(stepID) || isnan(stepID)
                return
            end

            if any(app.RuntimeFailedStepIDs == stepID)
                state = 'failed';
                return
            end

            if any(app.RuntimeSkippedStepIDs == stepID)
                state = 'skipped';
                return
            end

            if any(app.RuntimeCompletedStepIDs == stepID)
                state = 'passed';
                return
            end

            if ~isnan(app.RuntimeCurrentStepID) && app.RuntimeCurrentStepID == stepID
                if app.bCancelRequested
                    state = 'stopping';
                else
                    state = 'running';
                end
                return
            end

            if app.bPipelineRunning
                state = 'pending';
            end
        end

        function icon = getRuntimeStateIcon(~, state)
            %GETRUNTIMESTATEICON Return compact status icon for HTML views.

            switch lower(char(string(state)))
                case 'passed'
                    icon = '✓';
                case 'failed'
                    icon = '✕';
                case 'skipped'
                    icon = '↷';
                case 'running'
                    icon = '●';
                case 'stopping'
                    icon = '■';
                otherwise
                    icon = '';
            end
        end

        function refreshRuntimeHtmlViews(app)
            %REFRESHRUNTIMEHTMLVIEWS Refresh only HTML views affected by live state.
            %
            % This intentionally does not refresh the selected-step details panel.

            try
                app.refreshGraphView();
            catch
            end

            try
                app.refreshExecutionView();
            catch
            end
        end

        function value = getNumericProgressEventField(~, evt, fieldName, defaultValue)
            %GETNUMERICPROGRESSEVENTFIELD Read one numeric progress event field.

            value = defaultValue;

            try
                if isstruct(evt) && isfield(evt, fieldName) && ~isempty(evt.(fieldName))
                    raw = evt.(fieldName);
                    if isnumeric(raw) && isscalar(raw)
                        value = double(raw);
                    elseif islogical(raw) && isscalar(raw)
                        value = double(raw);
                    elseif ischar(raw) || (isstring(raw) && isscalar(raw))
                        tmp = str2double(char(string(raw)));
                        if ~isnan(tmp)
                            value = tmp;
                        end
                    end
                end
            catch
                value = defaultValue;
            end
        end

        function out = appendUniqueID(~, in, idValue)
            %APPENDUNIQUEID Append one numeric ID to a row vector if valid and absent.

            out = in;

            if isempty(idValue) || isnan(idValue)
                return
            end

            idValue = double(idValue);

            if isempty(out)
                out = idValue;
            elseif ~any(out == idValue)
                out(end+1) = idValue;
            end
        end

        function out = removeIDFromList(~, in, idValue)
            %REMOVEIDFROMLIST Remove one numeric ID from a vector.

            out = in;

            if isempty(out) || isempty(idValue) || isnan(idValue)
                return
            end

            out = out(out ~= double(idValue));
        end

        function value = nanToEmptyNumeric(~, value)
            %NANTOEMPTYNUMERIC Return [] for NaN so JSON payload is cleaner.

            if isempty(value) || ~isnumeric(value) || ~isscalar(value) || isnan(value)
                value = [];
            end
        end

        function vec = numericVectorForPayload(~, vec)
            %NUMERICVECTORFORPAYLOAD Return a finite row vector for uihtml payloads.

            if isempty(vec)
                vec = [];
                return
            end

            vec = double(vec(:).');
            vec = vec(isfinite(vec));
            vec = unique(vec, 'stable');
        end

        function configureForDataViewerMode(app)
            %CONFIGUREFORDATAVIEWERMODE Attach managed DataViewer execution context.
            %
            % DataViewer mode uses a normal window so the app menu bar remains visible.
            % DataViewer itself is hidden while this tool is open, so the user cannot
            % modify the viewer behind the PM Tool.

            if isempty(app.DataViewerSaveFolder) || ~isfolder(app.DataViewerSaveFolder)
                app.appendDiagnostic('DataViewer mode failed: invalid SaveFolder.');
                app.setStatus('Invalid DataViewer SaveFolder.');
                return
            end

            if isempty(app.DataViewerRawFolder) || ...
                    (~strcmpi(app.DataViewerRawFolder, 'Missing') && ~isfolder(app.DataViewerRawFolder))
                app.DataViewerRawFolder = 'Missing';
            end

            try
                app.UIFigure.WindowStyle = 'normal';
                figure(app.UIFigure);
            catch
            end

            try
                if isempty(app.pm)
                    app.pm = PipelineManager(app.DataViewerSaveFolder, app.DataViewerRawFolder);
                end

                app.pm.setLeafOutputPolicy('viewerTemp');

                try
                    app.pm.b_inputFromDataViewer = true;
                catch
                end

                app.FolderPairTable = app.getFolderPairTableFromPipelineManager();

                app.ensureDataViewerPresetFileSource();
                app.syncExecutionControlsFromManager();
                app.lockDataViewerModeControls();

                app.setPipelineControlsAvailable(true);
                app.updateFolderSelectionButtonSummary();
                app.updateLatestRunSummaryMenuState();
                app.populateFunctionTree();
                app.refreshAllViews();

                app.appendDiagnostic(sprintf('Opened from DataViewer. SaveFolder: %s', app.DataViewerSaveFolder));
                app.appendDiagnostic(sprintf('RawFolder: %s', app.DataViewerRawFolder));
                app.appendDiagnostic('Final output policy locked to temporary viewer files.');
                app.setStatus('DataViewer mode active. Pipeline context is locked.');

            catch ME
                app.appendDiagnostic(sprintf('Failed to configure DataViewer mode: %s', ME.message));
                app.setStatus('Failed to configure DataViewer mode.');
            end
        end

        function ensureDataViewerPresetFileSource(app)
            %ENSUREDATAVIEWERPRESETFILESOURCE Update only the DataViewer-managed FileSource.
            %
            % DataViewer mode treats the persistent pipeline as a reusable recipe:
            %
            %   DataViewer input -> Step 1 -> Step 2 -> ...
            %
            % On each PM Tool opening, only the FileSource previously registered as
            % the DataViewer input is replaced/reconnected to the current DataViewer
            % file. User-added parallel FileSource nodes are left untouched.
            %
            % The GUI must not mutate app.pm.nodes/app.pm.connections directly. It
            % uses only PipelineManager public APIs:
            %   1) ensure/reuse the current DataViewer FileSource,
            %   2) reconnect consumers of the stored DataViewer FileSource only,
            %   3) remove only that old stored FileSource if it becomes disconnected.

            if isempty(app.pm) || isempty(app.DataViewerCurrentFilePath) || ...
                    ~isfile(app.DataViewerCurrentFilePath)
                return
            end

            try
                semanticTypes = app.getDataViewerCurrentFileSemanticTypes(app.DataViewerCurrentFilePath);

                % ---------------------------------------------------------------------
                % 1) Ensure the current DataViewer file exists as a backend FileSource.
                % ---------------------------------------------------------------------
                [newStepTag, newNodeID] = app.pm.ensureFileSourceStep( ...
                    app.DataViewerCurrentFilePath, ...
                    'SemanticTypes', semanticTypes);

                newNodeID = double(newNodeID);
                newStepTag = char(string(newStepTag));

                if isempty(newStepTag) || ~isfinite(newNodeID)
                    error('PipelineManagerTool:DataViewerFileSourceNotCreated', ...
                        'Could not create or find the current DataViewer FileSource node.');
                end

                % ---------------------------------------------------------------------
                % 2) Reconnect only the previously registered DataViewer input source.
                % ---------------------------------------------------------------------
                oldNodeID = iGetStoredDataViewerFileSourceNodeID(newNodeID);

                if ~isempty(oldNodeID) && isfinite(oldNodeID) && oldNodeID ~= newNodeID
                    oldStepTag = iGetNodeStepTag(oldNodeID);

                    if ~isempty(oldStepTag)
                        % Move only downstream DATA inputs currently sourced by the old
                        % DataViewer FileSource to the new current FileSource. Other
                        % FileSource nodes, including user-added parallel branches, are
                        % never inferred or touched here.
                        iReconnectDownstreamConsumers(oldNodeID, newStepTag);

                        % Remove only the old DataViewer FileSource, and only if it is
                        % now fully disconnected. Do not clean unrelated FileSources.
                        iRemoveFolderNodeIfDisconnected(oldNodeID, oldStepTag);
                    end
                end

                app.DataViewerFileSourceNodeID = newNodeID;
                app.selectedNodeID = newNodeID;
                app.selectedStepTag = newStepTag;
                app.selectedFunctionName = '';

                try
                    app.pm.autoValidate();
                catch
                end

                app.appendDiagnostic(sprintf('DataViewer input FileSource set to: %s', app.DataViewerCurrentFilePath));

            catch ME
                app.appendDiagnostic(sprintf('Could not preset DataViewer file source: %s', ME.message));
            end

            % =====================================================================
            % Local helpers
            % =====================================================================
            function nodeIDOut = iGetStoredDataViewerFileSourceNodeID(newNodeIDIn)
                %IGETSTOREDDATAVIEWERFILESOURCENODEID Return only the stored DV input node.
                %
                % Important: do not fallback to "first connected FileSource". That
                % fallback can steal user-created parallel FileSource branches when
                % reopening the PM Tool.

                nodeIDOut = NaN;

                if isempty(app.pm) || isempty(app.pm.nodes)
                    return
                end

                storedID = double(app.DataViewerFileSourceNodeID);

                if ~isfinite(storedID) || storedID == newNodeIDIn
                    return
                end

                if ~iNodeExists(storedID)
                    return
                end

                if ~iIsFolderNode(storedID)
                    return
                end

                nodeIDOut = storedID;
            end

            function iReconnectDownstreamConsumers(oldNodeIDIn, newStepTagIn)
                if isempty(app.pm.connections)
                    return
                end

                oldNodeIDIn = double(oldNodeIDIn);

                connList = app.pm.connections([app.pm.connections.sourceNodeID] == oldNodeIDIn);
                if isempty(connList)
                    return
                end

                targetIDs = unique([connList.targetNodeID], 'stable');

                for iTarget = 1:numel(targetIDs)
                    targetID = targetIDs(iTarget);
                    targetIdx = find([app.pm.nodes.id] == targetID, 1, 'first');

                    if isempty(targetIdx)
                        continue
                    end

                    targetNode = app.pm.nodes(targetIdx);

                    if strcmpi(char(string(targetNode.kind)), 'folder')
                        continue
                    end

                    if ~isfield(targetNode, 'info') || ~isfield(targetNode.info, 'inputs') || ...
                            isempty(targetNode.info.inputs)
                        continue
                    end

                    allInputs = targetNode.info.inputs;
                    dataMask = arrayfun(@(x) isfield(x, 'isData') && x.isData, allInputs);
                    dataInputs = allInputs(dataMask);

                    if isempty(dataInputs)
                        continue
                    end

                    inputRefs = cell(1, numel(dataInputs));

                    for iIn = 1:numel(dataInputs)
                        inputName = char(string(dataInputs(iIn).name));
                        connIdx = iFindInputConnection(targetID, inputName);

                        if isempty(connIdx)
                            error('PipelineManagerTool:MissingExistingInput', ...
                                'Cannot reconnect "%s": input "%s" is not currently connected.', ...
                                char(string(targetNode.name)), inputName);
                        end

                        connLocal = app.pm.connections(connIdx);

                        if double(connLocal.sourceNodeID) == oldNodeIDIn
                            inputRefs{iIn} = sprintf('%s:data', newStepTagIn);
                        else
                            inputRefs{iIn} = iConnectionToExplicitRef(connLocal);
                        end
                    end

                    app.pm.reconnectStep(char(string(targetNode.name)), 'input', inputRefs);
                end
            end

            function idx = iFindInputConnection(targetNodeID, inputName)
                idx = [];

                if isempty(app.pm.connections)
                    return
                end

                idx = find([app.pm.connections.targetNodeID] == targetNodeID & ...
                    strcmpi({app.pm.connections.targetInputName}, inputName), 1, 'first');
            end

            function inputRef = iConnectionToExplicitRef(connLocal)
                inputRef = '';

                srcIdx = find([app.pm.nodes.id] == connLocal.sourceNodeID, 1, 'first');
                if isempty(srcIdx)
                    error('PipelineManagerTool:InvalidExistingConnection', ...
                        'Existing source node %d was not found.', connLocal.sourceNodeID);
                end

                srcNode = app.pm.nodes(srcIdx);
                srcStepTag = char(string(srcNode.name));

                srcToken = char(string(connLocal.sourceOutputName));
                if isfield(connLocal, 'selectedFile') && ~isempty(connLocal.selectedFile)
                    selectedFile = char(string(connLocal.selectedFile));
                    if ~isempty(strtrim(selectedFile))
                        srcToken = selectedFile;
                    end
                end

                inputRef = sprintf('%s:%s', srcStepTag, srcToken);
            end

            function iRemoveFolderNodeIfDisconnected(nodeIDIn, stepTagIn)
                if isempty(stepTagIn) || ~isfinite(nodeIDIn)
                    return
                end

                if ~iIsFolderNode(nodeIDIn)
                    return
                end

                if ~isempty(app.pm.connections) && any([app.pm.connections.sourceNodeID] == nodeIDIn)
                    return
                end

                if ~isempty(app.pm.connections) && any([app.pm.connections.targetNodeID] == nodeIDIn)
                    return
                end

                try
                    app.pm.rmStep(stepTagIn);
                    app.appendDiagnostic(sprintf('Removed previous DataViewer FileSource: %s', stepTagIn));
                catch ME
                    app.appendDiagnostic(sprintf('Could not remove previous DataViewer FileSource "%s": %s', stepTagIn, ME.message));
                end
            end

            function tf = iNodeExists(nodeIDIn)
                tf = false;

                if isempty(app.pm) || isempty(app.pm.nodes) || ~isfinite(nodeIDIn)
                    return
                end

                tf = any([app.pm.nodes.id] == double(nodeIDIn));
            end

            function tf = iIsFolderNode(nodeIDIn)
                tf = false;

                if isempty(app.pm) || isempty(app.pm.nodes) || ~isfinite(nodeIDIn)
                    return
                end

                idx = find([app.pm.nodes.id] == double(nodeIDIn), 1, 'first');
                if isempty(idx)
                    return
                end

                try
                    tf = isfield(app.pm.nodes(idx), 'kind') && ...
                        strcmpi(char(string(app.pm.nodes(idx).kind)), 'folder');
                catch
                    tf = false;
                end
            end

            function stepTag = iGetNodeStepTag(nodeIDIn)
                stepTag = '';

                if isempty(app.pm) || isempty(app.pm.nodes) || ~isfinite(nodeIDIn)
                    return
                end

                idx = find([app.pm.nodes.id] == double(nodeIDIn), 1, 'first');
                if isempty(idx)
                    return
                end

                stepTag = char(string(app.pm.nodes(idx).name));
            end
        end

        function tf = isDataViewerProtectedNodeID(app, nodeID)
            %ISDATAVIEWERPROTECTEDNODEID PM Tool-only protection for DataViewer input.

            tf = app.bDataViewerMode && ...
                ~isempty(app.DataViewerFileSourceNodeID) && ...
                isfinite(app.DataViewerFileSourceNodeID) && ...
                isfinite(nodeID) && ...
                double(nodeID) == double(app.DataViewerFileSourceNodeID);
        end

        function semanticTypes = getDataViewerCurrentFileSemanticTypes(~, filePath)
            %GETDATAVIEWERCURRENTFILESEMANTICTYPES Return permissive source types.
            %
            % UnknownDataType is intentionally included so the current DataViewer file
            % can feed functions that declare ImageTimeSeries, Image, or ProcessedData.

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

        function lockDataViewerModeControls(app)
            %LOCKDATAVIEWERMODECONTROLS Lock DataViewer-owned controls.

            if ~app.bDataViewerMode
                return
            end

            try
                app.AutoSaveFinalOutputsDropDown.Value = 'viewerTemp';
                app.AutoSaveFinalOutputsDropDown.Enable = 'off';
            catch
            end

            try
                app.AutoSaveFinalOutputsLabel.Enable = 'off';
            catch
            end

            try
                app.SelectDataFoldersButton.Enable = 'off';
            catch
            end

            try
                app.SelectDataFoldersMenu.Enable = 'off';
            catch
            end
        end

        function notifyDataViewerExecutionFinished(app, result)
            %NOTIFYDATAVIEWEREXECUTIONFINISHED Send execution result back to DataViewer.

            if ~app.bDataViewerMode || isempty(app.DataViewerExecutionFinishedFcn)
                return
            end

            try
                app.bDataViewerExecutionResultNotified = true;
                app.DataViewerExecutionFinishedFcn(result);
                app.appendDiagnostic('Execution result returned to DataViewer.');
            catch ME
                app.appendDiagnostic(sprintf('Failed to notify DataViewer: %s', ME.message));
            end
        end

        function notifyDataViewerToolClosed(app)
            %NOTIFYDATAVIEWERTOOLCLOSED Inform DataViewer that the PM tool closed.

            if ~app.bDataViewerMode || isempty(app.DataViewerToolClosedFcn)
                return
            end

            toolState = struct();
            toolState.status = "closed";
            toolState.closedOn = datetime('now');
            toolState.lastExecutionResult = app.LastExecutionResult;
            toolState.wasRunning = app.bPipelineRunning;
            toolState.dataViewerFileSourceNodeID = app.DataViewerFileSourceNodeID;

            try
                app.DataViewerToolClosedFcn(toolState);
            catch
            end
        end

        function setPipelineControlsAvailable(app, tf)
            %SETPIPELINECONTROLSAVAILABLE Enable/disable controls requiring a PipelineManager.
            %
            % In standalone mode, folder selection remains available because it creates
            % or replaces the PipelineManager context.
            %
            % In DataViewer mode, DataViewer owns SaveFolder/RawFolder and the leaf
            % output policy. Therefore:
            %   - data-folder selection is locked
            %   - final-output policy is locked to viewerTemp
            %   - normal pipeline editing/execution controls still follow tf

            if tf
                state = 'on';
            else
                state = 'off';
            end

            bDataViewerMode = false;
            try
                bDataViewerMode = isprop(app, 'bDataViewerMode') && logical(app.bDataViewerMode);
            catch
                bDataViewerMode = false;
            end

            % ---------------------------------------------------------------------
            % Top controls that require an attached PipelineManager context
            % ---------------------------------------------------------------------
            topControls = { ...
                app.LoadButton, ...
                app.SaveButton, ...
                app.AutoSaveFinalOutputsDropDown, ...
                app.RAMModeDropDown, ...
                app.RAMSafetyPolicyDropDown, ...
                app.RunButton, ...
                app.ReuseExistingFilesCheckBox};

            for iCtrl = 1:numel(topControls)
                try
                    topControls{iCtrl}.Enable = state;
                catch
                end
            end

            try
                app.ClearPipelineButton.Enable = state;
            catch
            end

            topLabels = { ...
                app.AutoSaveFinalOutputsLabel, ...
                app.RAMLabel, ...
                app.SafetyLabel};

            for iLabel = 1:numel(topLabels)
                try
                    topLabels{iLabel}.Enable = state;
                catch
                end
            end

            % ---------------------------------------------------------------------
            % Available-function controls
            % ---------------------------------------------------------------------
            sideControls = { ...
                app.SearchFunctionEditField, ...
                app.FunctionTree};

            for iCtrl = 1:numel(sideControls)
                try
                    sideControls{iCtrl}.Enable = state;
                catch
                end
            end

            if ~tf
                try
                    app.AddButton.Enable = 'off';
                catch
                end

                try
                    app.AddasNewBranchButton.Enable = 'off';
                catch
                end
            end

            % ---------------------------------------------------------------------
            % Menus
            % ---------------------------------------------------------------------
            try
                app.setMenuControlsAvailable(tf);
            catch
            end

            % ---------------------------------------------------------------------
            % Data-folder selection
            % ---------------------------------------------------------------------
            if bDataViewerMode
                % DataViewer owns SaveFolder/RawFolder in this mode.
                try
                    app.SelectDataFoldersButton.Enable = 'off';
                catch
                end

                try
                    app.SelectDataFoldersMenu.Enable = 'off';
                catch
                end
            else
                % Standalone mode: always allow folder selection because it creates or
                % replaces the execution context.
                try
                    app.SelectDataFoldersButton.Enable = 'on';
                catch
                end

                try
                    app.SelectDataFoldersMenu.Enable = 'on';
                catch
                end
            end

            % ---------------------------------------------------------------------
            % Leaf-output policy lock for DataViewer mode
            % ---------------------------------------------------------------------
            if bDataViewerMode
                try
                    app.AutoSaveFinalOutputsDropDown.Value = 'viewerTemp';
                catch
                end

                try
                    app.AutoSaveFinalOutputsDropDown.Enable = 'off';
                catch
                end

                try
                    app.AutoSaveFinalOutputsLabel.Enable = 'off';
                catch
                end

                try
                    if ~isempty(app.pm)
                        app.pm.setLeafOutputPolicy('viewerTemp');
                    end
                catch ME
                    app.appendDiagnostic(sprintf( ...
                        'Failed to enforce DataViewer output policy: %s', ME.message));
                end
            end

            % ---------------------------------------------------------------------
            % RAM safety dependent enable/disable
            % ---------------------------------------------------------------------
            try
                app.updateRamSafetyControlState();
            catch
            end

            % Re-apply DataViewer locks after updateRamSafetyControlState because that
            % method may modify dropdown state.
            if bDataViewerMode
                try
                    app.AutoSaveFinalOutputsDropDown.Value = 'viewerTemp';
                    app.AutoSaveFinalOutputsDropDown.Enable = 'off';
                    app.AutoSaveFinalOutputsLabel.Enable = 'off';
                catch
                end

                try
                    app.SelectDataFoldersButton.Enable = 'off';
                    app.SelectDataFoldersMenu.Enable = 'off';
                catch
                end
            end

            % ---------------------------------------------------------------------
            % Running state
            % ---------------------------------------------------------------------
            if app.bPipelineRunning
                % Do not call prepareGuiForPipelineRun here because that resets progress.
                % Only enforce the running UI state.
                try
                    app.setRunButtonState('running');
                catch
                end

                try
                    app.LoadButton.Enable = 'off';
                    app.SaveButton.Enable = 'off';
                    app.ClearPipelineButton.Enable = 'off';
                    app.SearchFunctionEditField.Enable = 'off';
                    app.FunctionTree.Enable = 'off';
                    app.AddButton.Enable = 'off';
                    app.AddasNewBranchButton.Enable = 'off';
                    app.AutoSaveFinalOutputsDropDown.Enable = 'off';
                    app.RAMModeDropDown.Enable = 'off';
                    app.RAMSafetyPolicyDropDown.Enable = 'off';
                    app.ReuseExistingFilesCheckBox.Enable = 'off';
                    app.SelectDataFoldersButton.Enable = 'off';
                    app.SelectDataFoldersMenu.Enable = 'off';
                catch
                end
            end
        end

    end

    methods (Access = public)

    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, varargin)
            app.initializeUiState();

            app.parseStartupArguments(varargin{:});

            try
                app.ensureHtmlFilesExist();
                app.appendDiagnostic(sprintf('Loaded graph HTML: %s', app.graphHtmlFile));
                app.appendDiagnostic(sprintf('Loaded execution HTML: %s', app.execHtmlFile));
            catch ME
                app.appendDiagnostic(sprintf('Failed to initialize HTML components: %s', ME.message));
                app.setStatus('Failed to initialize HTML components.');
            end

            if app.bDataViewerMode
                app.configureForDataViewerMode();
                return
            end

            app.populateFunctionTree();
            app.refreshAllViews();

            app.setPipelineControlsAvailable(false);
            app.updateFolderSelectionButtonSummary();

            app.appendDiagnostic('App initialized. Select data folders to start.');
            app.setStatus('Select data folders to start.');



        end

        % Callback function: LoadButton, LoadPipelineMenu
        function LoadButtonPushed(app, event)
            if isempty(app.pm)
                app.appendDiagnostic('Load requested before selecting data folders.');
                app.setStatus('Select data folders before loading a pipeline.');
                return
            end

            [fileName, folderName] = uigetfile({'*.pipe', 'Pipeline files (*.pipe)'}, 'Load pipeline');
            if isequal(fileName, 0)
                app.appendDiagnostic('Load pipeline cancelled by user.');
                app.setStatus('Load cancelled.');
                return
            end

            pipeFile = fullfile(folderName, fileName);

            try
                app.pm.loadPipe(pipeFile);
                if app.bDataViewerMode
                    app.pm.setLeafOutputPolicy('viewerTemp');
                    app.ensureDataViewerPresetFileSource();
                    app.lockDataViewerModeControls();
                end

                if ~app.bDataViewerMode
                    app.LastExecutionResult = [];
                end
                app.updateLatestRunSummaryMenuState();

                app.selectedNodeID = NaN;
                app.selectedStepTag = '';
                app.exitReconnectMode();

                app.FolderPairTable = app.getFolderPairTableFromPipelineManager();
                app.setPipelineControlsAvailable(true);
                app.updateFolderSelectionButtonSummary();

                try
                    app.syncExecutionControlsFromManager();
                catch ME
                    app.appendDiagnostic(sprintf('Loaded pipeline, but failed to sync execution controls: %s', ME.message));
                end

                app.populateFunctionTree();
                app.refreshAllViews();

                try
                    report = app.pm.diagnosePipeline('verbose', false);
                    app.appendValidationReportLogSummary(report, 'Loaded-pipeline validation');
                catch ME
                    app.appendDiagnostic(sprintf('Loaded pipeline, but validation summary failed: %s', ME.message));
                end

                app.appendDiagnostic(sprintf('Loaded pipeline: %s', pipeFile));
                app.setStatus('Pipeline loaded.');

            catch ME
                app.appendDiagnostic(sprintf('Failed to load pipeline: %s', ME.message));
                app.setStatus('Failed to load pipeline.');
            end

        end

        % Callback function: SaveButton, SavePipelineMenu
        function SaveButtonPushed(app, event)
            if isempty(app.pm)
                app.appendDiagnostic('Save requested without an attached PipelineManager.');
                app.setStatus('No PipelineManager attached.');
                return
            end

            [fileName, folderName] = uiputfile({'*.pipe', 'Pipeline files (*.pipe)'}, ...
                'Save pipeline as', fullfile(pwd, 'MyPipeline.pipe'));

            if isequal(fileName, 0)
                app.appendDiagnostic('Save pipeline cancelled by user.');
                app.setStatus('Save cancelled.');
                return
            end

            [~,~,ext] = fileparts(fileName);
            if isempty(ext)
                fileName = [fileName '.pipe'];
            end

            pipeFile = fullfile(folderName, fileName);

            try
                app.applyExecutionControlsToManager();
                app.pm.savePipe(pipeFile);
                app.appendDiagnostic(sprintf('Saved pipeline to %s', pipeFile));
                app.setStatus('Pipeline saved.');

            catch ME
                app.appendDiagnostic(sprintf('Failed to save pipeline: %s', ME.message));
                app.setStatus('Failed to save pipeline.');
            end

        end

        % Value changed function: SearchFunctionEditField
        function SearchFunctionEditFieldValueChanged(app, event)
            app.populateFunctionTree();

            if isempty(strtrim(char(string(app.SearchFunctionEditField.Value))))
                app.appendDiagnostic('Function search cleared.');
            else
                app.appendDiagnostic(sprintf('Function search filter updated: %s', char(string(app.SearchFunctionEditField.Value))));
            end

            app.setStatus('Function tree updated.');


        end

        % Button pushed function: AddButton
        function AddButtonPushed(app, event)
            app.tryAddSelectedFunction(false);
        end

        % Button pushed function: AddasNewBranchButton
        function AddasNewBranchButtonPushed(app, event)
            app.tryAddSelectedFunction(true);
        end

        % Data changed function: GraphHTML
        function GraphHTMLDataChanged(app, event)
            % This callback is currently not used.
            % Graph interactions are handled through HTMLEventReceived.


        end

        % Event received function: GraphHTML
        function GraphHTMLEventReceived(app, event)

            if app.bPipelineRunning
                app.appendDiagnostic('Graph interactions are disabled while the pipeline is running.');
                app.setStatus('Graph interactions disabled during execution.');
                return
            end

            eventName = lower(char(string(event.HTMLEventName)));
            eventData = event.HTMLEventData;

            switch eventName

                case 'htmlready'
                    if isstruct(eventData) && isfield(eventData, 'component') && strcmpi(char(string(eventData.component)), 'graph')
                        app.bGraphHtmlReady = true;
                        app.refreshGraphView();
                        app.appendDiagnostic('HTML graph is ready.');
                        app.setStatus('Graph view ready.');
                    end

                case 'graphrendered'
                    % No action needed.

                case 'nodeselected'
                    if isstruct(eventData) && isfield(eventData, 'nodeID')
                        app.selectedNodeID = double(eventData.nodeID);

                        if isfield(eventData, 'stepTag')
                            app.selectedStepTag = char(string(eventData.stepTag));
                        else
                            app.selectedStepTag = '';
                        end

                        app.refreshSelectedStepPanel();
                        app.refreshGraphView();
                        app.refreshExecutionView();

                        if ~isempty(app.selectedStepTag)
                            app.setStatus(['Selected step: ' app.selectedStepTag]);
                        else
                            app.setStatus('Step selected.');
                        end
                    end

                case 'reconnectsourcechosen'
                    if isstruct(eventData) && isfield(eventData, 'ref')
                        chosenRef = char(string(eventData.ref));
                        if ~isempty(chosenRef)
                            app.completeReconnectWithRef(chosenRef);
                        end
                    end

                case 'reconnectusefile'
                    if app.bReconnectMode
                        chosenRef = app.promptFileInputReference(app.reconnectTargetInputName);
                        if isempty(chosenRef)
                            app.appendDiagnostic('Reconnect file selection cancelled.');
                            app.setStatus('Reconnect still active.');
                        else
                            app.completeReconnectWithRef(chosenRef);
                        end
                    end

                case 'cancelreconnect'
                    if app.bReconnectMode
                        app.exitReconnectMode();
                        app.appendDiagnostic('Reconnect mode cancelled.');
                        app.setStatus('Reconnect cancelled.');
                    end

                case 'contextaction'
                    if ~isstruct(eventData) || ~isfield(eventData, 'action')
                        return
                    end

                    action = lower(char(string(eventData.action)));
                    stepTag = '';

                    if isfield(eventData, 'stepTag')
                        stepTag = char(string(eventData.stepTag));
                    end

                    switch action

                        case 'removestep'
                            if isempty(stepTag) || isempty(app.pm)
                                return
                            end
                            nodeIDToRemove = NaN;
                            try
                                nodeIdxToRemove = find(strcmpi(stepTag, {app.pm.nodes.name}), 1, 'first');
                                if ~isempty(nodeIdxToRemove)
                                    nodeIDToRemove = app.pm.nodes(nodeIdxToRemove).id;
                                end
                            catch
                                nodeIDToRemove = NaN;
                            end

                            if app.isDataViewerProtectedNodeID(nodeIDToRemove)
                                app.appendDiagnostic('DataViewer input FileSource cannot be removed.');
                                app.setStatus('DataViewer input FileSource is protected.');
                                app.refreshGraphView();
                                return
                            end



                            choice = uiconfirm(app.UIFigure, ...
                                sprintf('Remove step "%s"?', stepTag), ...
                                'Remove Step', ...
                                'Options', {'Remove', 'Cancel'}, ...
                                'DefaultOption', 2, ...
                                'CancelOption', 2);

                            if ~strcmp(choice, 'Remove')
                                app.appendDiagnostic('Remove-step action cancelled.');
                                app.setStatus('Remove-step cancelled.');
                                return
                            end

                            try
                                app.pm.rmStep(stepTag);

                                if strcmpi(app.selectedStepTag, stepTag)
                                    app.selectedNodeID = NaN;
                                    app.selectedStepTag = '';
                                end

                                app.exitReconnectMode();
                                app.refreshAllViews();
                                app.appendDiagnostic(sprintf('Removed step "%s".', stepTag));
                                app.setStatus('Step removed.');

                            catch ME
                                app.appendDiagnostic(sprintf('Failed to remove step "%s": %s', stepTag, ME.message));
                                app.setStatus('Failed to remove step.');
                            end

                        case 'editparameters'
                            if isempty(stepTag) || isempty(app.pm)
                                return
                            end

                            try
                                app.pm.setParameters(stepTag);
                                app.refreshAllViews();
                                app.appendDiagnostic(sprintf('Parameter dialog completed for "%s".', stepTag));
                                app.setStatus('Parameters updated.');

                            catch ME
                                app.appendDiagnostic(sprintf('Failed to edit parameters for "%s": %s', stepTag, ME.message));
                                app.setStatus('Failed to edit parameters.');
                            end

                        case 'reconnectstep'
                            if isempty(stepTag)
                                return
                            end
                            app.reconnectStepInteractive(stepTag);

                        case 'configuresaveoutputs'
                            if isempty(stepTag)
                                return
                            end
                            app.configureSaveTargetsInteractive(stepTag);

                        otherwise
                            app.appendDiagnostic(sprintf('Unknown graph context action: %s', action));
                            app.setStatus('Unknown graph action.');
                    end
            end

        end

        % Data changed function: ExecGraphHTML
        function ExecGraphHTMLDataChanged(app, event)

            if app.bPipelineRunning
                app.appendDiagnostic('Execution-view interactions are disabled while the pipeline is running.');
                app.setStatus('Execution-view interactions disabled during execution.');
                return
            end


        end

        % Event received function: ExecGraphHTML
        function ExecGraphHTMLEventReceived(app, event)

            eventName = lower(char(string(event.HTMLEventName)));
            eventData = event.HTMLEventData;

            switch eventName

                case 'htmlready'
                    if isstruct(eventData) && isfield(eventData, 'component') && strcmpi(char(string(eventData.component)), 'execution')
                        app.bExecHtmlReady = true;
                        app.refreshExecutionView();
                        app.appendDiagnostic('Execution-order HTML view is ready.');
                        app.setStatus('Execution view ready.');
                    end

                case 'stepselected'
                    if isstruct(eventData) && isfield(eventData, 'nodeID')
                        app.selectedNodeID = double(eventData.nodeID);

                        if isfield(eventData, 'stepTag')
                            app.selectedStepTag = char(string(eventData.stepTag));
                        else
                            app.selectedStepTag = '';
                        end

                        app.refreshSelectedStepPanel();
                        app.refreshGraphView();
                        app.refreshExecutionView();

                        if ~isempty(app.selectedStepTag)
                            app.setStatus(['Selected step: ' app.selectedStepTag]);
                        else
                            app.setStatus('Execution step selected.');
                        end
                    end
            end


        end

        % Close request function: UIFigure
        function UIFigureCloseRequest(app, event)
            app.notifyDataViewerToolClosed();
            app.cleanupAuxiliaryFigures();
            delete(app)

        end

        % Selection changed function: FunctionTree
        function FunctionTreeSelectionChanged(app, event)

            selectedNodes = app.FunctionTree.SelectedNodes;

            if isempty(selectedNodes)
                app.selectedFunctionName = '';
                app.AddButton.Enable = 'off';
                app.AddasNewBranchButton.Enable = 'off';
                return
            end

            selectedNode = selectedNodes(1);
            nodeData = selectedNode.NodeData;

            if isempty(nodeData) || ~isstruct(nodeData) || ~isfield(nodeData, 'kind')
                app.selectedFunctionName = '';
                app.AddButton.Enable = 'off';
                app.AddasNewBranchButton.Enable = 'off';
                return
            end

            if strcmpi(nodeData.kind, 'function') && isfield(nodeData, 'name')
                app.selectedFunctionName = char(string(nodeData.name));
                app.AddButton.Enable = 'on';
                app.AddasNewBranchButton.Enable = 'on';
                app.appendDiagnostic(sprintf('Selected function: %s', app.selectedFunctionName));
                app.setStatus('Function selected.');
            else
                app.selectedFunctionName = '';
                app.AddButton.Enable = 'off';
                app.AddasNewBranchButton.Enable = 'off';
            end

        end

        % Callback function: SelectDataFoldersButton, SelectDataFoldersMenu
        function SelectDataFoldersButtonPushed(app, event)
            if app.bDataViewerMode
                app.appendDiagnostic('Data folders are controlled by DataViewer in this mode.');
                app.setStatus('Data folder selection is locked.');
                return
            end

            app.selectDataFoldersInteractive();


        end

        % Value changed function: AutoSaveFinalOutputsDropDown
        function AutoSaveFinalOutputsDropDownValueChanged(app, event)
            value = app.AutoSaveFinalOutputsDropDown.Value;

            if app.bDataViewerMode
                app.AutoSaveFinalOutputsDropDown.Value = 'viewerTemp';
                app.appendDiagnostic('Final output policy is controlled by DataViewer and locked to temporary viewer files.');
                app.setStatus('Output policy locked by DataViewer.');
                return
            end

            if isempty(app.pm)
                app.appendDiagnostic(sprintf('Auto-save final outputs setting changed: %s', char(string(value))));
                app.setStatus('Output policy updated.');
                return
            end

            try
                app.pm.setLeafOutputPolicy(value);
                app.refreshAllViews();
                app.appendDiagnostic(sprintf('Auto-save final outputs setting changed: %s', char(string(value))));
                app.setStatus('Output policy updated.');
            catch ME
                app.appendDiagnostic(sprintf('Failed to update output policy: %s', ME.message));
                app.syncExecutionControlsFromManager();
                app.setStatus('Failed to update output policy.');
            end


        end

        % Value changed function: RAMModeDropDown
        function RAMModeDropDownValueChanged(app, event)
            value = app.RAMModeDropDown.Value;
            app.updateRamSafetyControlState();

            if isempty(app.pm)
                app.appendDiagnostic(sprintf('RAM mode changed: %s', char(string(value))));
                app.setStatus('RAM mode updated.');
                return
            end

            try
                app.pm.ramMode = char(string(value));
                app.refreshAllViews();
                app.appendDiagnostic(sprintf('RAM mode changed: %s', char(string(value))));
                app.setStatus('RAM mode updated.');
            catch ME
                app.appendDiagnostic(sprintf('Failed to update RAM mode: %s', ME.message));
                app.syncExecutionControlsFromManager();
                app.setStatus('Failed to update RAM mode.');
            end

        end

        % Value changed function: RAMSafetyPolicyDropDown
        function RAMSafetyPolicyDropDownValueChanged(app, event)
            value = app.RAMSafetyPolicyDropDown.Value;

            if isempty(app.pm)
                app.appendDiagnostic(sprintf('RAM safety policy changed: %s', char(string(value))));
                app.setStatus('RAM safety updated.');
                return
            end

            try
                app.pm.ramSafePolicy = char(string(value));
                app.refreshAllViews();
                app.appendDiagnostic(sprintf('RAM safety policy changed: %s', char(string(value))));
                app.setStatus('RAM safety updated.');
            catch ME
                app.appendDiagnostic(sprintf('Failed to update RAM safety policy: %s', ME.message));
                app.syncExecutionControlsFromManager();
                app.setStatus('Failed to update RAM safety.');
            end

        end

        % Callback function: RunButton, RunPipelineMenu
        function RunButtonPushed(app, event)

            if app.bPipelineRunning
                app.requestPipelineStop();
            else
                app.runPipelineFromGui();
            end

        end

        % Button pushed function: ClearPipelineButton
        function ClearPipelineButtonPushed(app, event)

            %CLEARPIPELINEBUTTONPUSHED Clear graph while preserving execution context.

            if isempty(app.pm)
                app.appendDiagnostic('No PipelineManager is attached.');
                app.setStatus('No PipelineManager attached.');
                return
            end

            try
                app.pm.clearPipeline();

                app.selectedNodeID = NaN;
                app.selectedStepTag = '';
                app.selectedFunctionName = '';

                app.bReconnectMode = false;
                app.reconnectTargetNodeID = NaN;
                app.reconnectTargetStepTag = '';
                app.reconnectTargetInputName = '';
                app.reconnectCandidateNodeIDs = [];

                app.resetRuntimeExecutionState(true);

                if app.bDataViewerMode
                    app.ensureDataViewerPresetFileSource();
                    app.lockDataViewerModeControls();
                    app.appendDiagnostic('Pipeline cleared. DataViewer file source was restored.');
                    app.setStatus('Pipeline cleared. DataViewer source restored.');
                else
                    app.appendDiagnostic('Pipeline cleared.');
                    app.setStatus('Pipeline cleared.');
                end

                app.populateFunctionTree();
                app.refreshAllViews();

            catch ME
                app.appendDiagnostic(sprintf('Failed to clear pipeline: %s', ME.message));
                app.setStatus('Failed to clear pipeline.');
            end

        end

        % Menu selected function: ReuseExistingFilesMenu
        function ReuseExistingFilesMenuSelected(app, event)
            currentState = false;
            try
                currentState = strcmpi(app.ReuseExistingFilesMenu.Checked, 'on');
            catch
            end

            newState = ~currentState;
            app.setReuseExistingFilesState(newState, true);
            app.refreshAllViews();

            if newState
                app.appendDiagnostic('Reuse Existing Files enabled.');
                app.setStatus('Reuse existing files enabled.');
            else
                app.appendDiagnostic('Reuse Existing Files disabled. Pipeline run will force re-execution where possible.');
                app.setStatus('Reuse existing files disabled.');
            end


        end

        % Menu selected function: ExecutionSettingsMenu
        function ExecutionSettingsMenuSelected(app, event)
            app.showExecutionSettingsInteractive();

        end

        % Menu selected function: AdvancedRAMSettingsMenu
        function AdvancedRAMSettingsMenuSelected(app, event)
            app.showAdvancedRamSettingsInteractive();

        end

        % Menu selected function: ShowExecutionPlanMenu
        function ShowExecutionPlanMenuSelected(app, event)
            app.showExecutionPlanInteractive();

        end

        % Menu selected function: ExportErrorLogMenu
        function ExportErrorLogMenuSelected(app, event)
            app.exportErrorLogInteractive();

        end

        % Menu selected function: RemoveInvalidStepsMenu
        function RemoveInvalidStepsMenuSelected(app, event)
            app.removeInvalidStepsInteractive();

        end

        % Menu selected function: GenerateScriptMenu
        function GenerateScriptMenuSelected(app, event)
            app.generateScriptInteractive();
        end

        % Value changed function: ReuseExistingFilesCheckBox
        function ReuseExistingFilesCheckBoxValueChanged(app, event)
            value = logical(app.ReuseExistingFilesCheckBox.Value);
            app.setReuseExistingFilesState(value, true);
            app.refreshAllViews();

            if value
                app.appendDiagnostic('Reuse Existing Files enabled.');
                app.setStatus('Reuse existing files enabled.');
            else
                app.appendDiagnostic('Reuse Existing Files disabled. Pipeline run will force re-execution where possible.');
                app.setStatus('Reuse existing files disabled.');
            end

        end

        % Menu selected function: ShowLatestRunSummaryMenu
        function ShowLatestRunSummaryMenuSelected(app, event)
            app.showPipelineSummaryInteractive();

        end

        % Menu selected function: ViewFolderPipelineLogMenu
        function ViewFolderPipelineLogMenuSelected(app, event)
            app.showFolderPipelineLogInteractive();
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1253 871];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, @UIFigureCloseRequest, true);

            % Create FileMenu
            app.FileMenu = uimenu(app.UIFigure);
            app.FileMenu.Text = 'File';

            % Create LoadPipelineMenu
            app.LoadPipelineMenu = uimenu(app.FileMenu);
            app.LoadPipelineMenu.MenuSelectedFcn = createCallbackFcn(app, @LoadButtonPushed, true);
            app.LoadPipelineMenu.Text = 'Load Pipeline';

            % Create SavePipelineMenu
            app.SavePipelineMenu = uimenu(app.FileMenu);
            app.SavePipelineMenu.MenuSelectedFcn = createCallbackFcn(app, @SaveButtonPushed, true);
            app.SavePipelineMenu.Text = 'Save Pipeline';

            % Create GenerateScriptMenu
            app.GenerateScriptMenu = uimenu(app.FileMenu);
            app.GenerateScriptMenu.MenuSelectedFcn = createCallbackFcn(app, @GenerateScriptMenuSelected, true);
            app.GenerateScriptMenu.Text = 'Generate Script';

            % Create DataMenu
            app.DataMenu = uimenu(app.UIFigure);
            app.DataMenu.Text = 'Data';

            % Create SelectDataFoldersMenu
            app.SelectDataFoldersMenu = uimenu(app.DataMenu);
            app.SelectDataFoldersMenu.MenuSelectedFcn = createCallbackFcn(app, @SelectDataFoldersButtonPushed, true);
            app.SelectDataFoldersMenu.Text = 'Select Data Folders';

            % Create ExecutionMenu
            app.ExecutionMenu = uimenu(app.UIFigure);
            app.ExecutionMenu.Text = 'Execution';

            % Create RunPipelineMenu
            app.RunPipelineMenu = uimenu(app.ExecutionMenu);
            app.RunPipelineMenu.MenuSelectedFcn = createCallbackFcn(app, @RunButtonPushed, true);
            app.RunPipelineMenu.Text = 'Run Pipeline';

            % Create ReuseExistingFilesMenu
            app.ReuseExistingFilesMenu = uimenu(app.ExecutionMenu);
            app.ReuseExistingFilesMenu.MenuSelectedFcn = createCallbackFcn(app, @ReuseExistingFilesMenuSelected, true);
            app.ReuseExistingFilesMenu.Checked = 'on';
            app.ReuseExistingFilesMenu.Text = 'Reuse Existing Files';

            % Create ExecutionSettingsMenu
            app.ExecutionSettingsMenu = uimenu(app.ExecutionMenu);
            app.ExecutionSettingsMenu.MenuSelectedFcn = createCallbackFcn(app, @ExecutionSettingsMenuSelected, true);
            app.ExecutionSettingsMenu.Separator = 'on';
            app.ExecutionSettingsMenu.Text = 'Execution Settings...';

            % Create AdvancedRAMSettingsMenu
            app.AdvancedRAMSettingsMenu = uimenu(app.ExecutionMenu);
            app.AdvancedRAMSettingsMenu.MenuSelectedFcn = createCallbackFcn(app, @AdvancedRAMSettingsMenuSelected, true);
            app.AdvancedRAMSettingsMenu.Text = 'Advanced RAM Settings...';

            % Create ToolsMenu
            app.ToolsMenu = uimenu(app.UIFigure);
            app.ToolsMenu.Text = 'Tools';

            % Create ShowExecutionPlanMenu
            app.ShowExecutionPlanMenu = uimenu(app.ToolsMenu);
            app.ShowExecutionPlanMenu.MenuSelectedFcn = createCallbackFcn(app, @ShowExecutionPlanMenuSelected, true);
            app.ShowExecutionPlanMenu.Text = 'Show Execution Plan';

            % Create RemoveInvalidStepsMenu
            app.RemoveInvalidStepsMenu = uimenu(app.ToolsMenu);
            app.RemoveInvalidStepsMenu.MenuSelectedFcn = createCallbackFcn(app, @RemoveInvalidStepsMenuSelected, true);
            app.RemoveInvalidStepsMenu.Text = 'Remove Invalid Steps';

            % Create ShowLatestRunSummaryMenu
            app.ShowLatestRunSummaryMenu = uimenu(app.ToolsMenu);
            app.ShowLatestRunSummaryMenu.MenuSelectedFcn = createCallbackFcn(app, @ShowLatestRunSummaryMenuSelected, true);
            app.ShowLatestRunSummaryMenu.Text = 'Show Latest Run Summary';

            % Create ViewFolderPipelineLogMenu
            app.ViewFolderPipelineLogMenu = uimenu(app.ToolsMenu);
            app.ViewFolderPipelineLogMenu.MenuSelectedFcn = createCallbackFcn(app, @ViewFolderPipelineLogMenuSelected, true);
            app.ViewFolderPipelineLogMenu.Separator = 'on';
            app.ViewFolderPipelineLogMenu.Text = 'View Folder Pipeline Log';

            % Create ExportErrorLogMenu
            app.ExportErrorLogMenu = uimenu(app.ToolsMenu);
            app.ExportErrorLogMenu.MenuSelectedFcn = createCallbackFcn(app, @ExportErrorLogMenuSelected, true);
            app.ExportErrorLogMenu.Text = 'Export Error Log';

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {280, '1x', 280};
            app.GridLayout.RowHeight = {50, '1x', 180};

            % Create GridLayoutTop
            app.GridLayoutTop = uigridlayout(app.GridLayout);
            app.GridLayoutTop.ColumnWidth = {160, 90, 90, 3, 130, '1x', 30, 90, 40, 90, 3, 130, 90};
            app.GridLayoutTop.RowHeight = {'1x'};
            app.GridLayoutTop.Layout.Row = 1;
            app.GridLayoutTop.Layout.Column = [1 3];

            % Create LoadButton
            app.LoadButton = uibutton(app.GridLayoutTop, 'push');
            app.LoadButton.ButtonPushedFcn = createCallbackFcn(app, @LoadButtonPushed, true);
            app.LoadButton.Tooltip = {'Load a saved pipeline definition from a .pipe file.'};
            app.LoadButton.Layout.Row = 1;
            app.LoadButton.Layout.Column = 2;
            app.LoadButton.Text = 'Load';

            % Create SaveButton
            app.SaveButton = uibutton(app.GridLayoutTop, 'push');
            app.SaveButton.ButtonPushedFcn = createCallbackFcn(app, @SaveButtonPushed, true);
            app.SaveButton.Tooltip = {'Save the current pipeline definition to a .pipe file.'};
            app.SaveButton.Layout.Row = 1;
            app.SaveButton.Layout.Column = 3;
            app.SaveButton.Text = 'Save';

            % Create SelectDataFoldersButton
            app.SelectDataFoldersButton = uibutton(app.GridLayoutTop, 'push');
            app.SelectDataFoldersButton.ButtonPushedFcn = createCallbackFcn(app, @SelectDataFoldersButtonPushed, true);
            app.SelectDataFoldersButton.Tooltip = {'Choose or edit the SaveFolder/RawFolder pairs used to create the pipeline context.'};
            app.SelectDataFoldersButton.Layout.Row = 1;
            app.SelectDataFoldersButton.Layout.Column = 1;
            app.SelectDataFoldersButton.Text = 'Select Data Folders';

            % Create RunButton
            app.RunButton = uibutton(app.GridLayoutTop, 'push');
            app.RunButton.ButtonPushedFcn = createCallbackFcn(app, @RunButtonPushed, true);
            app.RunButton.Tooltip = {'Execute the current pipeline using the selected folder pairs and execution settings.'};
            app.RunButton.Layout.Row = 1;
            app.RunButton.Layout.Column = 13;
            app.RunButton.Text = 'Run';

            % Create AutoSaveFinalOutputsDropDown
            app.AutoSaveFinalOutputsDropDown = uidropdown(app.GridLayoutTop);
            app.AutoSaveFinalOutputsDropDown.Items = {'Disabled', 'Enabled', 'Enabled as temporary files'};
            app.AutoSaveFinalOutputsDropDown.ValueChangedFcn = createCallbackFcn(app, @AutoSaveFinalOutputsDropDownValueChanged, true);
            app.AutoSaveFinalOutputsDropDown.Tooltip = {'Controls whether final pipeline outputs are automatically saved. Manually selected outputs are always saved.'};
            app.AutoSaveFinalOutputsDropDown.Placeholder = 'Output Policy';
            app.AutoSaveFinalOutputsDropDown.Layout.Row = 1;
            app.AutoSaveFinalOutputsDropDown.Layout.Column = 6;
            app.AutoSaveFinalOutputsDropDown.Value = 'Disabled';

            % Create RAMModeDropDown
            app.RAMModeDropDown = uidropdown(app.GridLayoutTop);
            app.RAMModeDropDown.Items = {'Auto', 'RAM-safe'};
            app.RAMModeDropDown.ValueChangedFcn = createCallbackFcn(app, @RAMModeDropDownValueChanged, true);
            app.RAMModeDropDown.Tooltip = {'Choose how pipeline data is passed between steps: automatically or using file-backed RAM-safe execution.'};
            app.RAMModeDropDown.Placeholder = 'RAM Mode';
            app.RAMModeDropDown.Layout.Row = 1;
            app.RAMModeDropDown.Layout.Column = 8;
            app.RAMModeDropDown.Value = 'Auto';

            % Create RAMSafetyPolicyDropDown
            app.RAMSafetyPolicyDropDown = uidropdown(app.GridLayoutTop);
            app.RAMSafetyPolicyDropDown.Items = {'Strict', 'Best effort'};
            app.RAMSafetyPolicyDropDown.ValueChangedFcn = createCallbackFcn(app, @RAMSafetyPolicyDropDownValueChanged, true);
            app.RAMSafetyPolicyDropDown.Tooltip = {'Choose how strictly RAM-safe mode handles steps that cannot receive file-backed inputs.'};
            app.RAMSafetyPolicyDropDown.Placeholder = 'RAM safety';
            app.RAMSafetyPolicyDropDown.Layout.Row = 1;
            app.RAMSafetyPolicyDropDown.Layout.Column = 10;
            app.RAMSafetyPolicyDropDown.Value = 'Best effort';

            % Create AutoSaveFinalOutputsLabel
            app.AutoSaveFinalOutputsLabel = uilabel(app.GridLayoutTop);
            app.AutoSaveFinalOutputsLabel.Layout.Row = 1;
            app.AutoSaveFinalOutputsLabel.Layout.Column = 5;
            app.AutoSaveFinalOutputsLabel.Text = 'Auto-save final outputs:';

            % Create RAMLabel
            app.RAMLabel = uilabel(app.GridLayoutTop);
            app.RAMLabel.Layout.Row = 1;
            app.RAMLabel.Layout.Column = 7;
            app.RAMLabel.Text = 'RAM:';

            % Create SafetyLabel
            app.SafetyLabel = uilabel(app.GridLayoutTop);
            app.SafetyLabel.Layout.Row = 1;
            app.SafetyLabel.Layout.Column = 9;
            app.SafetyLabel.Text = 'Safety:';

            % Create Sep1
            app.Sep1 = uipanel(app.GridLayoutTop);
            app.Sep1.BorderType = 'none';
            app.Sep1.BackgroundColor = [0.8 0.8 0.8];
            app.Sep1.Layout.Row = 1;
            app.Sep1.Layout.Column = 4;

            % Create Sep2
            app.Sep2 = uipanel(app.GridLayoutTop);
            app.Sep2.BorderType = 'none';
            app.Sep2.BackgroundColor = [0.8 0.8 0.8];
            app.Sep2.Layout.Row = 1;
            app.Sep2.Layout.Column = 11;

            % Create ReuseExistingFilesCheckBox
            app.ReuseExistingFilesCheckBox = uicheckbox(app.GridLayoutTop);
            app.ReuseExistingFilesCheckBox.ValueChangedFcn = createCallbackFcn(app, @ReuseExistingFilesCheckBoxValueChanged, true);
            app.ReuseExistingFilesCheckBox.Tooltip = {'When enabled, completed outputs that already exist and match pipeline history may be reused. Disable to force re-execution.'};
            app.ReuseExistingFilesCheckBox.Text = 'Reuse Existing Files';
            app.ReuseExistingFilesCheckBox.Layout.Row = 1;
            app.ReuseExistingFilesCheckBox.Layout.Column = 12;

            % Create AvailableFunctionsPanel
            app.AvailableFunctionsPanel = uipanel(app.GridLayout);
            app.AvailableFunctionsPanel.Title = 'Available Functions';
            app.AvailableFunctionsPanel.Layout.Row = [2 3];
            app.AvailableFunctionsPanel.Layout.Column = 1;
            app.AvailableFunctionsPanel.FontWeight = 'bold';

            % Create GridAvailableFunctions
            app.GridAvailableFunctions = uigridlayout(app.AvailableFunctionsPanel);
            app.GridAvailableFunctions.RowHeight = {30, '1x', 40};

            % Create SearchFunctionEditField
            app.SearchFunctionEditField = uieditfield(app.GridAvailableFunctions, 'text');
            app.SearchFunctionEditField.ValueChangedFcn = createCallbackFcn(app, @SearchFunctionEditFieldValueChanged, true);
            app.SearchFunctionEditField.Placeholder = 'Search function ...';
            app.SearchFunctionEditField.Layout.Row = 1;
            app.SearchFunctionEditField.Layout.Column = [1 2];

            % Create FunctionTree
            app.FunctionTree = uitree(app.GridAvailableFunctions);
            app.FunctionTree.SelectionChangedFcn = createCallbackFcn(app, @FunctionTreeSelectionChanged, true);
            app.FunctionTree.Layout.Row = 2;
            app.FunctionTree.Layout.Column = [1 2];

            % Create AddButton
            app.AddButton = uibutton(app.GridAvailableFunctions, 'push');
            app.AddButton.ButtonPushedFcn = createCallbackFcn(app, @AddButtonPushed, true);
            app.AddButton.Enable = 'off';
            app.AddButton.Layout.Row = 3;
            app.AddButton.Layout.Column = 1;
            app.AddButton.Text = 'Add ';

            % Create AddasNewBranchButton
            app.AddasNewBranchButton = uibutton(app.GridAvailableFunctions, 'push');
            app.AddasNewBranchButton.ButtonPushedFcn = createCallbackFcn(app, @AddasNewBranchButtonPushed, true);
            app.AddasNewBranchButton.Enable = 'off';
            app.AddasNewBranchButton.Layout.Row = 3;
            app.AddasNewBranchButton.Layout.Column = 2;
            app.AddasNewBranchButton.Text = 'Add as New Branch';

            % Create PipelineGraphPanel
            app.PipelineGraphPanel = uipanel(app.GridLayout);
            app.PipelineGraphPanel.Title = 'Pipeline Graph ';
            app.PipelineGraphPanel.Layout.Row = 2;
            app.PipelineGraphPanel.Layout.Column = 2;
            app.PipelineGraphPanel.FontWeight = 'bold';

            % Create GridPipelineGraph
            app.GridPipelineGraph = uigridlayout(app.PipelineGraphPanel);
            app.GridPipelineGraph.ColumnWidth = {'1x', 100};
            app.GridPipelineGraph.RowHeight = {30, '1x'};
            app.GridPipelineGraph.ColumnSpacing = 2;
            app.GridPipelineGraph.RowSpacing = 2;

            % Create GraphHTML
            app.GraphHTML = uihtml(app.GridPipelineGraph);
            app.GraphHTML.DataChangedFcn = createCallbackFcn(app, @GraphHTMLDataChanged, true);
            app.GraphHTML.HTMLEventReceivedFcn = createCallbackFcn(app, @GraphHTMLEventReceived, true);
            app.GraphHTML.Layout.Row = 2;
            app.GraphHTML.Layout.Column = [1 2];

            % Create ClearPipelineButton
            app.ClearPipelineButton = uibutton(app.GridPipelineGraph, 'push');
            app.ClearPipelineButton.ButtonPushedFcn = createCallbackFcn(app, @ClearPipelineButtonPushed, true);
            app.ClearPipelineButton.Layout.Row = 1;
            app.ClearPipelineButton.Layout.Column = 2;
            app.ClearPipelineButton.Text = 'Clear Pipeline';

            % Create ExecutionProgressPanel
            app.ExecutionProgressPanel = uipanel(app.GridPipelineGraph);
            app.ExecutionProgressPanel.BorderType = 'none';
            app.ExecutionProgressPanel.Layout.Row = 1;
            app.ExecutionProgressPanel.Layout.Column = 1;

            % Create GridExecutionProgress
            app.GridExecutionProgress = uigridlayout(app.ExecutionProgressPanel);
            app.GridExecutionProgress.ColumnWidth = {'1x'};
            app.GridExecutionProgress.RowSpacing = 2;
            app.GridExecutionProgress.Padding = [0 0 0 0];

            % Create FolderProgressOuterPanel
            app.FolderProgressOuterPanel = uipanel(app.GridExecutionProgress);
            app.FolderProgressOuterPanel.BorderType = 'none';
            app.FolderProgressOuterPanel.Layout.Row = 1;
            app.FolderProgressOuterPanel.Layout.Column = 1;

            % Create StepProgressOuterPanel
            app.StepProgressOuterPanel = uipanel(app.GridExecutionProgress);
            app.StepProgressOuterPanel.BorderType = 'none';
            app.StepProgressOuterPanel.Layout.Row = 2;
            app.StepProgressOuterPanel.Layout.Column = 1;

            % Create SelectedstepdetailsPanel
            app.SelectedstepdetailsPanel = uipanel(app.GridLayout);
            app.SelectedstepdetailsPanel.Title = 'Selected step details';
            app.SelectedstepdetailsPanel.Layout.Row = 2;
            app.SelectedstepdetailsPanel.Layout.Column = 3;
            app.SelectedstepdetailsPanel.FontWeight = 'bold';

            % Create GridSelectedStepDetails
            app.GridSelectedStepDetails = uigridlayout(app.SelectedstepdetailsPanel);
            app.GridSelectedStepDetails.ColumnWidth = {80, '1x'};
            app.GridSelectedStepDetails.RowHeight = {30, 30, 30, '1x', 30, '1x'};
            app.GridSelectedStepDetails.RowSpacing = 5;

            % Create StepNameLabel
            app.StepNameLabel = uilabel(app.GridSelectedStepDetails);
            app.StepNameLabel.Layout.Row = 1;
            app.StepNameLabel.Layout.Column = 1;
            app.StepNameLabel.Text = 'Step';

            % Create StepTagLabel
            app.StepTagLabel = uilabel(app.GridSelectedStepDetails);
            app.StepTagLabel.Layout.Row = 2;
            app.StepTagLabel.Layout.Column = 1;
            app.StepTagLabel.Text = 'Step Tag';

            % Create StepTagField
            app.StepTagField = uieditfield(app.GridSelectedStepDetails, 'text');
            app.StepTagField.Editable = 'off';
            app.StepTagField.BackgroundColor = [0.9412 0.9412 0.9412];
            app.StepTagField.Layout.Row = 2;
            app.StepTagField.Layout.Column = 2;

            % Create StepSummaryLabel
            app.StepSummaryLabel = uilabel(app.GridSelectedStepDetails);
            app.StepSummaryLabel.Layout.Row = 3;
            app.StepSummaryLabel.Layout.Column = 1;
            app.StepSummaryLabel.Text = 'Step Summary';

            % Create StepSummaryTextArea
            app.StepSummaryTextArea = uitextarea(app.GridSelectedStepDetails);
            app.StepSummaryTextArea.Editable = 'off';
            app.StepSummaryTextArea.BackgroundColor = [0.9412 0.9412 0.9412];
            app.StepSummaryTextArea.Layout.Row = 4;
            app.StepSummaryTextArea.Layout.Column = [1 2];

            % Create ParametersLabel
            app.ParametersLabel = uilabel(app.GridSelectedStepDetails);
            app.ParametersLabel.Layout.Row = 5;
            app.ParametersLabel.Layout.Column = 1;
            app.ParametersLabel.Text = 'Parameters';

            % Create ParametersTextArea
            app.ParametersTextArea = uitextarea(app.GridSelectedStepDetails);
            app.ParametersTextArea.Editable = 'off';
            app.ParametersTextArea.BackgroundColor = [0.9412 0.9412 0.9412];
            app.ParametersTextArea.Layout.Row = 6;
            app.ParametersTextArea.Layout.Column = [1 2];

            % Create FunctionNameField
            app.FunctionNameField = uieditfield(app.GridSelectedStepDetails, 'text');
            app.FunctionNameField.Editable = 'off';
            app.FunctionNameField.BackgroundColor = [0.9412 0.9412 0.9412];
            app.FunctionNameField.Layout.Row = 1;
            app.FunctionNameField.Layout.Column = 2;

            % Create ExecutionOrderPanel
            app.ExecutionOrderPanel = uipanel(app.GridLayout);
            app.ExecutionOrderPanel.Title = 'Execution Order';
            app.ExecutionOrderPanel.Layout.Row = 3;
            app.ExecutionOrderPanel.Layout.Column = 2;

            % Create GridExecutionOrder
            app.GridExecutionOrder = uigridlayout(app.ExecutionOrderPanel);
            app.GridExecutionOrder.ColumnWidth = {'1x'};
            app.GridExecutionOrder.RowHeight = {'1x'};

            % Create ExecGraphHTML
            app.ExecGraphHTML = uihtml(app.GridExecutionOrder);
            app.ExecGraphHTML.DataChangedFcn = createCallbackFcn(app, @ExecGraphHTMLDataChanged, true);
            app.ExecGraphHTML.HTMLEventReceivedFcn = createCallbackFcn(app, @ExecGraphHTMLEventReceived, true);
            app.ExecGraphHTML.Layout.Row = 1;
            app.ExecGraphHTML.Layout.Column = 1;

            % Create ActivityLogPanel
            app.ActivityLogPanel = uipanel(app.GridLayout);
            app.ActivityLogPanel.Title = 'Activity Log';
            app.ActivityLogPanel.Layout.Row = 3;
            app.ActivityLogPanel.Layout.Column = 3;

            % Create GridActivityLog
            app.GridActivityLog = uigridlayout(app.ActivityLogPanel);
            app.GridActivityLog.ColumnWidth = {'1x'};
            app.GridActivityLog.RowHeight = {'1x', 30};
            app.GridActivityLog.RowSpacing = 5;

            % Create ActivityLogTextArea
            app.ActivityLogTextArea = uitextarea(app.GridActivityLog);
            app.ActivityLogTextArea.Editable = 'off';
            app.ActivityLogTextArea.BackgroundColor = [0.9412 0.9412 0.9412];
            app.ActivityLogTextArea.Layout.Row = 1;
            app.ActivityLogTextArea.Layout.Column = 1;

            % Create StatusLabel
            app.StatusLabel = uilabel(app.GridActivityLog);
            app.StatusLabel.Layout.Row = 2;
            app.StatusLabel.Layout.Column = 1;
            app.StatusLabel.Text = 'Status';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = PipelineManagerTool_exported(varargin)

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