classdef DataViewer_Coreg2Cams_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        GridLayout                      matlab.ui.container.GridLayout
        StatusLabel                     matlab.ui.control.Label
        GrifLayoutLeft                  matlab.ui.container.GridLayout
        CalibrationSourcePanel          matlab.ui.container.Panel
        GridCalibrationSource           matlab.ui.container.GridLayout
        LoadCalibrationDataButton       matlab.ui.control.Button
        CalibrationfilePanel            matlab.ui.container.Panel
        GridGeometricPanel              matlab.ui.container.GridLayout
        ArchiveActiveCalibrationButton  matlab.ui.control.Button
        LoadBackupButton                matlab.ui.control.Button
        CalibrationFileInfoTextArea     matlab.ui.control.TextArea
        ManualfinetuningPanel           matlab.ui.container.Panel
        GridManualFineTuning            matlab.ui.container.GridLayout
        ResetButton                     matlab.ui.control.Button
        ApplyStepSpinner                matlab.ui.control.Spinner
        ApplyStepSpinnerLabel           matlab.ui.control.Label
        StepsizeEditField               matlab.ui.control.NumericEditField
        StepsizeEditFieldLabel          matlab.ui.control.Label
        TransformationDropDown          matlab.ui.control.DropDown
        TransformationDropDownLabel     matlab.ui.control.Label
        ViewmodeButtonGroup             matlab.ui.container.ButtonGroup
        FalseColorOverlayButton         matlab.ui.control.RadioButton
        AlternateViewButton             matlab.ui.control.RadioButton
        SaveCalibrationButton           matlab.ui.control.Button
        CoregistrationparametersPanel   matlab.ui.container.Panel
        GridCoregistrationPanel         matlab.ui.container.GridLayout
        FlipHorizontallyCheckBox        matlab.ui.control.CheckBox
        UseenhancedcontrastCheckBox     matlab.ui.control.CheckBox
        RunAutomaticcoregistrationButton  matlab.ui.control.Button
        UIAxes                          matlab.ui.control.UIAxes
    end


    properties (Access = private)
        % Parent app/context.
        MainApp = []
        appMode char = 'CalibrationUtility'

        % Calibration source images. These must represent unregistered camera
        % images. CurrentMovingImage is transformed into CurrentRegisteredMovingImage.
        CurrentFixedImage = []
        CurrentMovingImage = []
        CurrentRegisteredMovingImage = []
        EnhancedFixedImage = []
        EnhancedMovingImage = []
        HasCalibrationSource logical = false
        CalibrationSourceFolder char = ''
        CalibrationSourceDescription char = ''

        % Active calibration loaded from coreg2cam_tform.mat.
        ActiveTform = []
        ActiveTformInfo = struct()

        % Candidate calibration currently previewed in the app.
        CurrentTform = []
        CurrentTformInfo = struct()
        CandidateSource char = ''
        HasUnsavedCalibration logical = false

        % Manual adjustment state. These are applied on top of BaseTform.
        BaseTform = []
        ManualDx double = 0
        ManualDy double = 0
        ManualRotationDeg double = 0
        ManualScale double = 1

        % Backup management.
        BackupFiles = struct([])
        MaxBackupFiles double = 10

        % Alternate-view timer state.
        AlternateViewTimer = timer.empty
        AlternateViewImageHandle = []
        AlternateViewIndex double = 1
    end


    properties (Access = public)
        filename = 'coreg2cam_tform' % Active calibration filename without extension.
        tformFile                 % Scalar dir struct for active calibration file, or empty.
        DataFolder                % Optional context folder from the parent app.
        LabeoFolder               % .../Documents/LabeoTech/Config/umIT/tformFiles
    end

    methods (Access = {?DataViewer})

        function lookForTformFile(app)
            %LOOKFORTFORMFILE Refresh active calibration-file location and metadata.
            %
            %   Ensures the LabeoTech calibration folder exists and updates
            %   app.tformFile with the active coreg2cam_tform.mat file when found.

            app.LabeoFolder = getUmitFolder('tformFiles');
           
            app.tformFile = dir(fullfile(app.LabeoFolder, [app.filename '.mat']));

            if numel(app.tformFile) > 1
                [~, idxNewest] = max([app.tformFile.datenum]);
                app.tformFile = app.tformFile(idxNewest);
            end

            app.BackupFiles = app.listBackupFiles();
        end

        function saveTform(app, tform, tformInfo)
            %SAVETFORM Save active dual-camera calibration transform.
            %
            %   saveTform(app, tform, tformInfo) saves tform and tformInfo as the
            %   active coreg2cam_tform.mat file. If an active file already exists,
            %   it is backed up automatically before being overwritten. Older backup
            %   files are pruned after save.

            if isempty(app.LabeoFolder) || ~isfolder(app.LabeoFolder)
                app.lookForTformFile();
            end

            activeFile = fullfile(app.LabeoFolder, [app.filename '.mat']);

            if isfile(activeFile)
                app.backupActiveCalibration();
            end

            save(activeFile, 'tform', 'tformInfo');

            app.pruneBackupFiles();
            app.lookForTformFile();
            app.loadActiveCalibrationFromDisk();
        end


        function switchAppMode(app)
            % Switches the application between Coregistration Mode and Results View Mode.
            %
            % This function toggles the application mode between Coregistration Mode and Results View Mode.
            % In Coregistration Mode, it adjusts the UI to display coregistration-related information.
            % In Results View Mode, it adjusts the UI to display registered camera images and related controls.
            %

            app.loadCameraImages;  % Load composite images from each camera.
            app.ViewmodeButtonGroup.Visible = 'on';
            if strcmpi(app.appMode, 'CoregistrationMode')
                % In Coregistration Mode:

                % Hide ResultsViewMode panel:
                app.GrifLayoutLeft.RowHeight{4} = 0;
                %
                app.UseenhancedcontrastCheckBoxValueChanged('');
            else
                % In Results View Mode:

                % Hide CoregistrationMode panel and Execution button:
                app.GrifLayoutLeft.RowHeight([2,3]) = {0,0};
                % Check if there are backup files:
                bkpFiles = dir(fullfile(app.LabeoFolder, [app.filename '_*_bkp.mat']));
                if ~isempty(bkpFiles)
                    app.LoadBackupButton.Enable = 'on';
                else
                    app.LoadBackupButton.Enable = 'off';
                end
                % Check if there is a tformFile:
                if isempty(app.tformFile)
                    set(app.DiscardButton,'Enable','off','Text','No File!')
                    app.ShowinfolderButton.Enable = 'off';
                else
                    set(app.DiscardButton,'Enable','on','Text','Discard')
                    app.ShowinfolderButton.Enable = 'on';
                end
                app.ViewmodeButtonGroupSelectionChanged('');
            end
            set(app.UIAxes, 'XTick',[],'YTick',[],'DataAspectRatio',[1 1 1], 'XLim',[1 size(app.camImages{1},2)], 'YLim',[1 size(app.camImages{1},1)], 'Colormap',gray)
            set(app.UIAxes2, 'XTick',[],'YTick',[],'DataAspectRatio',[1 1 1], 'XLim',[1 size(app.camImages{2},2)], 'YLim',[1 size(app.camImages{2},1)], 'Colormap',gray)
        end
    end

    methods (Access = private)
        function loadCameraImages(app)
            %LOADCAMERAIMAGES Load unregistered calibration-source camera images.
            %
            %   The calibration source folder must contain AcqInfos.mat and the
            %   per-channel .dat files generated without applying dual-camera
            %   coregistration. This function does not import raw data by itself.

            if isempty(app.CalibrationSourceFolder) || ~isfolder(app.CalibrationSourceFolder)
                error('DataViewerCoreg2Cams:MissingCalibrationSource', ...
                    'Select a valid calibration-source folder first.');
            end
                

            %%%% Import Calibration recording %%%%

            % Get list of current files in the folder to preserve after
            % cleanup           
            fileInfoCalibFolder = dir(app.CalibrationSourceFolder);
            fileInfoCalibFolder([fileInfoCalibFolder.isdir]) = [];
            % Import data without applying the registration by calling
            % ImagesClassification directly.
            bkpFile = ['backup_' datestr(now(),'yyyymmdd_hhMMss')];
            calibFiles = ImagesClassification(app.CalibrationSourceFolder,app.CalibrationSourceFolder,1,1,'backupOpts',bkpFile);
            bkpFile = [bkpFile '.zip'];
            % Get AcqInfo:
            S = load(fullfile(app.CalibrationSourceFolder,'AcqInfos.mat'));
            AcqInfoStream = S.AcqInfoStream;
            illumFields = fieldnames(AcqInfoStream);
            illumFields = illumFields(startsWith(illumFields, 'Illumination'));

            if isempty(illumFields)
                error('DataViewerCoreg2Cams:NoIlluminations', ...
                    'No Illumination fields were found in AcqInfoStream.');
            end

            cam1List = {};
            cam2List = {};

            for iField = 1:numel(illumFields)
                thisInfo = AcqInfoStream.(illumFields{iField});

                if ~isfield(thisInfo, 'CamIdx') || ~isfield(thisInfo, 'Color')
                    continue
                end

                idx = thisInfo.CamIdx;
                chan = lower(char(string(thisInfo.Color)));
                % Manage Fluo Channel Name
                if contains(chan, 'fluo') 
                    tok= regexp(chan, '(\d+)\s*nm\b*', 'tokens');
                    if ~isempty(tok)
                        wavTag = tok{:}{:};
                        chan = ['fluo_' wavTag];
                    else
                        chan = 'fluo';
                    end
                end                                  
                % Manage Yellow Channel name
                if contains(chan, 'amber')
                    chan = 'yellow';
                end

                datName = [chan '.dat'];
                if idx == 1
                    cam1List{end+1} = datName; %#ok<AGROW>
                elseif idx == 2
                    cam2List{end+1} = datName; %#ok<AGROW>
                end
            end

            if isempty(cam1List) || isempty(cam2List)
                error('DataViewerCoreg2Cams:MissingCameraChannels', ...
                    'Could not identify at least one .dat channel for each camera.');
            end
            assert(all(ismember([cam1List, cam2List],calibFiles)), ...
                'There is a mismatch between the existing .dat files and the channel names extracted from AcqInfos.mat file.');

            app.CurrentFixedImage = app.buildCompositeCameraImage(cam1List);
            app.CurrentMovingImage = app.buildCompositeCameraImage(cam2List);

            if app.FlipHorizontallyCheckBox.Value
                app.CurrentMovingImage = fliplr(app.CurrentMovingImage);
            end

            app.CurrentFixedImage = app.normalizeImage(app.CurrentFixedImage);
            app.CurrentMovingImage = app.normalizeImage(app.CurrentMovingImage);

            app.EnhancedFixedImage = adapthisteq(app.CurrentFixedImage);
            app.EnhancedMovingImage = adapthisteq(app.CurrentMovingImage);

            app.HasCalibrationSource = true;
            app.CalibrationSourceDescription = app.CalibrationSourceFolder;
            % Restore previous files
            unzip(fullfile(app.CalibrationSourceFolder,bkpFile),app.CalibrationSourceFolder);
            delete(fullfile(app.CalibrationSourceFolder,bkpFile));
        end

        function configureInitialState(app)
            %CONFIGUREINITIALSTATE Initialize UI state before file/source loading.

            app.UIFigure.Visible = 'off';
            app.UIFigure.WindowStyle = 'modal';

            app.CalibrationFileInfoTextArea.Editable = 'off';
            app.CalibrationFileInfoTextArea.Value = {'No calibration information loaded.'};

            app.StepsizeEditField.Limits = [eps Inf];
            app.StepsizeEditField.Value = 1;
            app.updateStepSizeLabel();

            app.ApplyStepSpinner.Limits = [-1 1];
            app.ApplyStepSpinner.Step = 1;
            app.ApplyStepSpinner.Value = 0;

            app.RunAutomaticcoregistrationButton.Enable = 'off';
            app.SaveCalibrationButton.Enable = 'off';
            app.LoadBackupButton.Enable = 'off';
            app.ArchiveActiveCalibrationButton.Enable = 'off';
            app.setManualControlsEnabled(false);

            cla(app.UIAxes);
            title(app.UIAxes, 'Load calibration data to start');
            app.UIAxes.XTick = [];
            app.UIAxes.YTick = [];
            axis(app.UIAxes, 'image');
            colormap(app.UIAxes, gray);

            app.setLocalStatus('Ready. Load calibration data to preview or create a calibration.');
        end

        function setLocalStatus(app, msg)
            %SETLOCALSTATUS Update status label text.

            app.StatusLabel.Text = char(string(msg));
            drawnow limitrate
        end

        function setManualControlsEnabled(app, tf)
            %SETMANUALCONTROLSENABLED Enable or disable manual fine-tuning controls.

            state = matlab.lang.OnOffSwitchState(tf);
            app.TransformationDropDown.Enable = state;
            app.StepsizeEditField.Enable = state;
            app.ApplyStepSpinner.Enable = state;
            app.ResetButton.Enable = state;
        end

        function refreshButtonStates(app)
            %REFRESHBUTTONSTATES Enable controls from current calibration state.

            hasSource = app.HasCalibrationSource && ...
                ~isempty(app.CurrentFixedImage) && ~isempty(app.CurrentMovingImage);
            hasCandidate = ~isempty(app.CurrentTform);
            hasActive = ~isempty(app.tformFile) && isfile(fullfile(app.tformFile.folder, app.tformFile.name));
            hasBackups = ~isempty(app.BackupFiles);

            app.LoadCalibrationDataButton.Enable = 'on';
            app.RunAutomaticcoregistrationButton.Enable = matlab.lang.OnOffSwitchState(hasSource);
            app.LoadBackupButton.Enable = matlab.lang.OnOffSwitchState(hasBackups);
            app.ArchiveActiveCalibrationButton.Enable = matlab.lang.OnOffSwitchState(hasActive);
            app.SaveCalibrationButton.Enable = matlab.lang.OnOffSwitchState(~isempty(app.CurrentTform) && app.HasUnsavedCalibration);
            app.setManualControlsEnabled(hasSource && hasCandidate);
        end

        function refreshCalibrationFileInfo(app)
            %REFRESHCALIBRATIONFILEINFO Update calibration/source status text area.

            lines = strings(0, 1);

            activeFile = fullfile(app.LabeoFolder, [app.filename '.mat']);
            if isfile(activeFile)
                d = dir(activeFile);
                lines(end+1) = "Active calibration: found";
                lines(end+1) = "File: " + string(d.name);
                lines(end+1) = "Modified: " + string(datetime(d.datenum, 'ConvertFrom', 'datenum'));
            else
                lines(end+1) = "Active calibration: not found";
            end

            lines(end+1) = "Folder: " + string(app.LabeoFolder);
            lines(end+1) = "Backups: " + string(numel(app.BackupFiles)) + " / " + string(app.MaxBackupFiles);
            lines(end+1) = "";

            if app.HasCalibrationSource
                lines(end+1) = "Calibration source: loaded";
                lines(end+1) = "Source folder: " + string(app.CalibrationSourceFolder);
            else
                lines(end+1) = "Calibration source: not loaded";
            end

            lines(end+1) = "";
            if isempty(app.CandidateSource)
                lines(end+1) = "Preview candidate: none";
            else
                lines(end+1) = "Preview candidate: " + string(app.CandidateSource);
            end

            if app.HasUnsavedCalibration
                lines(end+1) = "Unsaved candidate: yes";
            else
                lines(end+1) = "Unsaved candidate: no";
            end

            if ~isempty(app.CurrentTform)
                lines(end+1) = sprintf('Manual adjustment: dx=%.3g px, dy=%.3g px, rot=%.3g deg, scale=%.6g', ...
                    app.ManualDx, app.ManualDy, app.ManualRotationDeg, app.ManualScale);
            end

            app.CalibrationFileInfoTextArea.Value = cellstr(lines);
        end

        function refreshAllState(app)
            %REFRESHALLSTATE Refresh all status, buttons, and preview.

            app.lookForTformFile();
            app.refreshCalibrationFileInfo();
            app.refreshButtonStates();
            app.refreshPreview();
        end

        function files = listBackupFiles(app)
            %LISTBACKUPFILES Return backup files sorted newest first.

            if isempty(app.LabeoFolder) || ~isfolder(app.LabeoFolder)
                files = struct([]);
                return
            end

            files = dir(fullfile(app.LabeoFolder, [app.filename '_*_bkp.mat']));
            if isempty(files)
                return
            end

            [~, idx] = sort([files.datenum], 'descend');
            files = files(idx);
        end

        function backupPath = backupActiveCalibration(app)
            %BACKUPACTIVECALIBRATION Copy active tform to timestamped backup.

            backupPath = '';
            activeFile = fullfile(app.LabeoFolder, [app.filename '.mat']);

            if ~isfile(activeFile)
                return
            end

            stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
            backupPath = fullfile(app.LabeoFolder, sprintf('%s_%s_bkp.mat', app.filename, stamp));
            copyfile(activeFile, backupPath);
        end

        function pruneBackupFiles(app)
            %PRUNEBACKUPFILES Keep only the newest MaxBackupFiles backup files.

            files = app.listBackupFiles();
            if numel(files) <= app.MaxBackupFiles
                return
            end

            for iFile = app.MaxBackupFiles+1:numel(files)
                thisPath = fullfile(files(iFile).folder, files(iFile).name);
                if isfile(thisPath)
                    delete(thisPath);
                end
            end
        end

        function loadActiveCalibrationFromDisk(app)
            %LOADACTIVECALIBRATIONFROMDISK Load active tform into active/candidate state.

            app.ActiveTform = [];
            app.ActiveTformInfo = struct();

            if isempty(app.tformFile)
                return
            end

            activePath = fullfile(app.tformFile.folder, app.tformFile.name);
            if ~isfile(activePath)
                return
            end

            S = load(activePath, 'tform', 'tformInfo');
            if ~isfield(S, 'tform')
                error('DataViewerCoreg2Cams:InvalidTformFile', ...
                    'The active calibration file does not contain variable "tform".');
            end

            app.ActiveTform = S.tform;
            if isfield(S, 'tformInfo')
                app.ActiveTformInfo = S.tformInfo;
            else
                app.ActiveTformInfo = struct();
            end

            app.setCandidateTransform(app.ActiveTform, app.ActiveTformInfo, ...
                'Active calibration', false, true);
        end

        function setCandidateTransform(app, tform, tformInfo, sourceLabel, isUnsaved, resetManual)
            %SETCANDIDATETRANSFORM Set current preview transform and metadata.

            if nargin < 6
                resetManual = true;
            end

            app.BaseTform = tform;
            app.CurrentTform = tform;
            app.CurrentTformInfo = tformInfo;
            app.CandidateSource = char(string(sourceLabel));
            app.HasUnsavedCalibration = logical(isUnsaved);

            if resetManual
                app.ManualDx = 0;
                app.ManualDy = 0;
                app.ManualRotationDeg = 0;
                app.ManualScale = 1;
            end

            app.updateRegisteredMovingImage();
            app.refreshCalibrationFileInfo();
            app.refreshButtonStates();
            app.refreshPreview();
        end

        function updateRegisteredMovingImage(app)
            %UPDATEREGISTEREDMOVINGIMAGE Apply current transform to moving image.

            app.CurrentRegisteredMovingImage = [];

            if isempty(app.CurrentTform) || isempty(app.CurrentMovingImage) || isempty(app.CurrentFixedImage)
                return
            end

            outRef = imref2d(size(app.CurrentFixedImage));
            app.CurrentRegisteredMovingImage = imwarp(app.CurrentMovingImage, app.CurrentTform, ...
                'OutputView', outRef, 'FillValues', 0);
            app.CurrentRegisteredMovingImage = app.normalizeImage(app.CurrentRegisteredMovingImage);
        end

        function img = getDisplayFixedImage(app)
            %GETDISPLAYFIXEDIMAGE Return fixed image with selected contrast mode.

            if app.UseenhancedcontrastCheckBox.Value && ~isempty(app.EnhancedFixedImage)
                img = app.EnhancedFixedImage;
            else
                img = app.CurrentFixedImage;
            end
        end

        function img = getDisplayMovingImage(app)
            %GETDISPLAYMOVINGIMAGE Return registered/unregistered moving image for display.

            if ~isempty(app.CurrentRegisteredMovingImage)
                img = app.CurrentRegisteredMovingImage;
            else
                img = app.CurrentMovingImage;
            end

            if app.UseenhancedcontrastCheckBox.Value
                img = adapthisteq(app.normalizeImage(img));
            else
                img = app.normalizeImage(img);
            end
        end

        function refreshPreview(app)
            %REFRESHPREVIEW Draw the current preview mode.

            app.stopAlternateViewTimer();
            cla(app.UIAxes);

            if ~app.HasCalibrationSource || isempty(app.CurrentFixedImage) || isempty(app.CurrentMovingImage)
                title(app.UIAxes, 'Load calibration data to preview calibration');
                app.UIAxes.XTick = [];
                app.UIAxes.YTick = [];
                axis(app.UIAxes, 'image');
                return
            end

            fixedImg = app.getDisplayFixedImage();
            movingImg = app.getDisplayMovingImage();

            if startsWith(app.ViewmodeButtonGroup.SelectedObject.Text, 'False', 'IgnoreCase', true)
                imshowpair(fixedImg, movingImg, 'falsecolor', 'Scaling', 'joint', 'Parent', app.UIAxes);
                title(app.UIAxes, app.getPreviewTitle('False-color overlay'), 'Interpreter', 'none');
                app.UIAxes.XTick = [];
                app.UIAxes.YTick = [];
                axis(app.UIAxes, 'image');
            else
                app.startAlternateView(fixedImg, movingImg);
            end
        end

        function titleText = getPreviewTitle(app, modeText)
            %GETPREVIEWTITLE Compose preview title.

            if isempty(app.CandidateSource)
                candidateText = 'unregistered source';
            else
                candidateText = app.CandidateSource;
            end

            titleText = sprintf('%s | %s', modeText, candidateText);
        end

        function startAlternateView(app, fixedImg, movingImg)
            %STARTALTERNATEVIEW Start 1 Hz flicker comparison on one axes.

            app.stopAlternateViewTimer();
            app.AlternateViewIndex = 1;

            app.AlternateViewImageHandle = imagesc(app.UIAxes, fixedImg, [0 1]);
            colormap(app.UIAxes, gray);
            axis(app.UIAxes, 'image');
            app.UIAxes.XTick = [];
            app.UIAxes.YTick = [];
            title(app.UIAxes, app.getPreviewTitle('Alternate view: Camera 1'), 'Interpreter', 'none');

            app.AlternateViewTimer = timer( ...
                'ExecutionMode', 'fixedRate', ...
                'Period', 1, ...
                'BusyMode', 'drop', ...
                'TimerFcn', @(~, ~) app.onAlternateViewTimerTick(fixedImg, movingImg));

            start(app.AlternateViewTimer);
        end

        function onAlternateViewTimerTick(app, fixedImg, movingImg)
            %ONALTERNATEVIEWTIMERTICK Alternate fixed and moving images in UIAxes.

            if isempty(app) || ~isvalid(app) || isempty(app.UIFigure) || ~isvalid(app.UIFigure)
                return
            end

            if isempty(app.AlternateViewImageHandle) || ~isvalid(app.AlternateViewImageHandle)
                return
            end

            if app.AlternateViewIndex == 1
                app.AlternateViewImageHandle.CData = movingImg;
                app.AlternateViewIndex = 2;
                title(app.UIAxes, app.getPreviewTitle('Alternate view: Camera 2'), 'Interpreter', 'none');
            else
                app.AlternateViewImageHandle.CData = fixedImg;
                app.AlternateViewIndex = 1;
                title(app.UIAxes, app.getPreviewTitle('Alternate view: Camera 1'), 'Interpreter', 'none');
            end

            drawnow limitrate
        end

        function stopAlternateViewTimer(app)
            %STOPALTERNATEVIEWTIMER Stop and delete alternate-view timer.

            if ~isempty(app.AlternateViewTimer)
                try
                    if isvalid(app.AlternateViewTimer)
                        stop(app.AlternateViewTimer);
                        delete(app.AlternateViewTimer);
                    end
                catch
                end
            end

            app.AlternateViewTimer = timer.empty;
            app.AlternateViewImageHandle = [];
        end

        function updateStepSizeLabel(app)
            %UPDATESTEPSIZELABEL Update step-size label unit from selected transform.

            switch app.TransformationDropDown.Value
                case {'Horizontal (transl.)', 'Vertical (transl.)'}
                    app.StepsizeEditFieldLabel.Text = 'Step size (px)';
                    if app.StepsizeEditField.Value <= 0
                        app.StepsizeEditField.Value = 1;
                    end
                case 'Rotation'
                    app.StepsizeEditFieldLabel.Text = 'Step size (deg)';
                    if app.StepsizeEditField.Value <= 0
                        app.StepsizeEditField.Value = 0.1;
                    end
                case 'Scaling'
                    app.StepsizeEditFieldLabel.Text = 'Step size (%)';
                    if app.StepsizeEditField.Value <= 0
                        app.StepsizeEditField.Value = 1;
                    end
            end
        end

        function applyManualStep(app, signedStep)
            %APPLYMANUALSTEP Apply one manual nudge to the candidate transform.

            if isempty(app.BaseTform) || isempty(app.CurrentTform)
                return
            end

            stepSize = double(app.StepsizeEditField.Value);
            if ~isfinite(stepSize) || stepSize <= 0
                error('DataViewerCoreg2Cams:InvalidStepSize', ...
                    'Step size must be a finite positive value.');
            end

            delta = signedStep * stepSize;

            switch app.TransformationDropDown.Value
                case 'Horizontal (transl.)'
                    app.ManualDx = app.ManualDx + delta;
                case 'Vertical (transl.)'
                    app.ManualDy = app.ManualDy + delta;
                case 'Rotation'
                    app.ManualRotationDeg = app.ManualRotationDeg + delta;
                case 'Scaling'
                    app.ManualScale = app.ManualScale * (1 + delta / 100);
                    app.ManualScale = max(app.ManualScale, eps);
            end

            app.recomposeCurrentTransform();
            app.HasUnsavedCalibration = true;
            app.CandidateSource = app.appendManualCandidateLabel(app.CandidateSource);

            app.updateRegisteredMovingImage();
            app.refreshCalibrationFileInfo();
            app.refreshButtonStates();
            app.refreshPreview();
            app.setLocalStatus('Manual adjustment applied. Save Calibration to make it active.');
        end

        function label = appendManualCandidateLabel(app, labelIn) %#ok<INUSL>
            %APPENDMANUALCANDIDATELABEL Mark source label as manually adjusted.

            label = char(string(labelIn));
            if isempty(label)
                label = 'Manual adjustment';
            elseif ~contains(label, '+ manual adjustment')
                label = [label '+ manual adjustment'];
            end
        end

        function recomposeCurrentTransform(app)
            %RECOMPOSECURRENTTRANSFORM Compose base transform plus manual adjustment.

            if isempty(app.BaseTform)
                return
            end

            fixedSize = size(app.CurrentFixedImage);
            cx = (fixedSize(2) + 1) / 2;
            cy = (fixedSize(1) + 1) / 2;

            TtoOrigin = [1 0 0; 0 1 0; -cx -cy 1];
            TfromOrigin = [1 0 0; 0 1 0; cx cy 1];

            theta = deg2rad(app.ManualRotationDeg);
            R = [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];
            S = [app.ManualScale 0 0; 0 app.ManualScale 0; 0 0 1];
            Tr = [1 0 0; 0 1 0; app.ManualDx app.ManualDy 1];

            manualT = TtoOrigin * S * R * TfromOrigin * Tr;
            app.CurrentTform = affine2d(app.BaseTform.T * manualT);

            app.CurrentTformInfo = app.updateTformInfoForManualAdjustment(app.CurrentTformInfo);
        end

        function info = updateTformInfoForManualAdjustment(app, info)
            %UPDATETFORMINFOFORMANUALADJUSTMENT Add manual-adjustment metadata.

            if isempty(info) || ~isstruct(info)
                info = struct();
            end

            info.ManualAdjustment = struct( ...
                'dx_px', app.ManualDx, ...
                'dy_px', app.ManualDy, ...
                'rotation_deg', app.ManualRotationDeg, ...
                'scale', app.ManualScale, ...
                'updatedOn', datetime('now'));

            info.CandidateSource = app.CandidateSource;
        end

        function img = buildCompositeCameraImage(app, datFileList)
            %BUILDCOMPOSITECAMERAIMAGE Create one normalized mean image from channels.

            if isempty(datFileList)
                error('DataViewerCoreg2Cams:EmptyDatList', ...
                    'No .dat files were provided for camera image construction.');
            end

            firstPath = fullfile(app.CalibrationSourceFolder, datFileList{1});
            if ~isfile(firstPath)
                error('DataViewerCoreg2Cams:MissingDatFile', ...
                    'Missing required channel file: %s', firstPath);
            end

            datMap = mapDat(firstPath);
            img = zeros([numel(datFileList), size(datMap.Data.data, 1), size(datMap.Data.data, 2)], 'single');
            clear datMap

            for iFile = 1:numel(datFileList)
                datPath = fullfile(app.CalibrationSourceFolder, datFileList{iFile});
                if ~isfile(datPath)
                    error('DataViewerCoreg2Cams:MissingDatFile', ...
                        'Missing required channel file: %s', datPath);
                end

                datMap = mapDat(datPath);
                nFrames = size(datMap.Data.data, 3);

                if nFrames > 100
                    frameIdx = round(linspace(1, nFrames, 100));
                    dat = mean(datMap.Data.data(:, :, frameIdx), 3);
                else
                    dat = mean(datMap.Data.data, 3);
                end

                if mean(dat(:), 'omitnan') > 1e3
                    dat = dat ./ mean(dat(:), 'omitnan');
                else
                    dat = zeros(size(dat), 'single');
                end

                img(iFile, :, :) = single(dat);
                clear datMap
            end

            img = squeeze(sum(img, 1));
        end

        function img = normalizeImage(app, img) %#ok<INUSL>
            %NORMALIZEIMAGE Normalize image to [0 1] safely.

            img = single(img);
            imgMin = min(img(:), [], 'omitnan');
            imgMax = max(img(:), [], 'omitnan');

            if ~isfinite(imgMin) || ~isfinite(imgMax) || imgMax <= imgMin
                img = zeros(size(img), 'single');
            else
                img = (img - imgMin) ./ (imgMax - imgMin);
            end
        end

        function ensureCurrentTformInfoHasRegisteredImages(app)
            %ENSURECURRENTTFORMINFOHASREGISTEREDIMAGES Store preview images in tformInfo.

            if isempty(app.CurrentTformInfo) || ~isstruct(app.CurrentTformInfo)
                app.CurrentTformInfo = struct();
            end

            if ~isempty(app.CurrentFixedImage) && ~isempty(app.CurrentRegisteredMovingImage)
                app.CurrentTformInfo.RegisteredImages = { ...
                    app.CurrentFixedImage, ...
                    app.CurrentRegisteredMovingImage};
            end

            app.CurrentTformInfo.SourceDataFolder = app.CalibrationSourceFolder;
            app.CurrentTformInfo.CandidateSource = app.CandidateSource;
            app.CurrentTformInfo.SavedOn = datetime('now');
        end

        function tf = hasCandidateTransform(app)
            %HASCANDIDATETRANSFORM Return true when a candidate transform exists.

            tf = ~isempty(app.CurrentTform);
        end

        function tf = hasActiveCalibration(app)
            %HASACTIVECALIBRATION Return true when active calibration exists on disk.

            activeFile = fullfile(app.LabeoFolder, [app.filename '.mat']);
            tf = isfile(activeFile);
        end

        function cleanupCoregistrationApp(app)
            %CLEANUPCOREGISTRATIONAPP Stop timers and release transient UI resources.

            try
                app.stopAlternateViewTimer();
            catch
            end
        end
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, parentApp, dataFolder, appMode)
            %STARTUPFCN Initialize the dual-camera calibration utility.
            %
            %   DataViewer_Coreg2Cams(parentApp, dataFolder)
            %   DataViewer_Coreg2Cams(parentApp, dataFolder, appMode)
            %
            %   The dataFolder argument is treated as context only. Calibration
            %   source data must be explicitly loaded by the user because processed
            %   .dat files may already have been coregistered during import.

            if nargin < 2
                parentApp = [];
            end

            if nargin < 3 || isempty(dataFolder)
                dataFolder = '';
            end

            if nargin < 4 || isempty(appMode)
                appMode = 'CalibrationUtility';
            end

            app.MainApp = parentApp;
            app.DataFolder = char(string(dataFolder));
            app.appMode = char(string(appMode));

            app.configureInitialState();
            app.lookForTformFile();
            app.loadActiveCalibrationFromDisk();

            app.refreshCalibrationFileInfo();
            app.refreshButtonStates();
            app.refreshPreview();

            if exist('placeAppInsideCaller', 'file') == 2 && ~isempty(parentApp)
                try
                    placeAppInsideCaller(parentApp, app.UIFigure, 'center');
                catch
                    movegui(app.UIFigure, 'center');
                end
            else
                movegui(app.UIFigure, 'center');
            end

            app.UIFigure.Visible = 'on';

            %             if strcmpi(AppMode, 'CoregistrationMode')
            %                 app.UseenhancedcontrastCheckBoxValueChanged('');
            %             end
        end

        % Close request function: UIFigure
        function UIFigureCloseRequest(app, event)

            %UIFIGURECLOSEREQUEST Close app after stopping preview timers.

            if app.HasUnsavedCalibration
                answer = uiconfirm(app.UIFigure, ...
                    'There is an unsaved calibration candidate. Close without saving?', ...
                    'Unsaved calibration', ...
                    'Options', {'Close without saving', 'Cancel'}, ...
                    'DefaultOption', 'Cancel', ...
                    'CancelOption', 'Cancel', ...
                    'Icon', 'warning');

                if strcmpi(answer, 'Cancel')
                    return
                end
            end

            app.cleanupCoregistrationApp();
            delete(app)

        end

        % Button pushed function: RunAutomaticcoregistrationButton
        function RunAutomaticcoregistrationButtonPushed(app, event)
            %RUNAUTOMATICCOREGISTRATIONBUTTONPUSHED Generate automatic candidate tform.

            if ~app.HasCalibrationSource
                uialert(app.UIFigure, ...
                    'Load unregistered calibration data before running automatic coregistration.', ...
                    'No calibration source', ...
                    'Icon', 'warning');
                return
            end

            try
                w = uiprogressdlg(app.UIFigure, ...
                    'Message', 'Calculating geometric transformation...', ...
                    'Title', 'Please wait', ...
                    'Indeterminate', 'on');
                cleanupObj = onCleanup(@() delete(w));

                [tform, tformInfo, warnmsg] = genTform2Cams(app.CalibrationSourceFolder, ...
                    app.UseenhancedcontrastCheckBox.Value, ...
                    app.FlipHorizontallyCheckBox.Value);

                if ~isempty(warnmsg)
                    uialert(app.UIFigure, warnmsg, ...
                        'Failed to generate geometric transformation file', ...
                        'Icon', 'warning');
                    app.setLocalStatus('Automatic coregistration failed.');
                    return
                end

                tformInfo.CandidateSource = 'Automatic coregistration';
                tformInfo.SourceDataFolder = app.CalibrationSourceFolder;
                tformInfo.GeneratedOn = datetime('now');

                app.setCandidateTransform(tform, tformInfo, ...
                    'Automatic coregistration', true, true);

                app.setLocalStatus('Automatic coregistration completed. Review and save the calibration if acceptable.');

            catch ME
                app.setLocalStatus(sprintf('Automatic coregistration failed: %s', ME.message));
                uialert(app.UIFigure, getReport(ME, 'basic', 'hyperlinks', 'off'), ...
                    'Automatic coregistration failed', ...
                    'Icon', 'error');
            end

        end

        % Selection changed function: ViewmodeButtonGroup
        function ViewmodeButtonGroupSelectionChanged(app, event)
            %VIEWMODEBUTTONGROUPSELECTIONCHANGED Refresh preview mode.

            app.refreshPreview();

        end

        % Value changed function: UseenhancedcontrastCheckBox
        function UseenhancedcontrastCheckBoxValueChanged(app, event)
            %USEENHANCEDCONTRASTCHECKBOXVALUECHANGED Refresh preview contrast mode.

            app.refreshPreview();

        end

        % Value changed function: FlipHorizontallyCheckBox
        function FlipHorizontallyCheckBoxValueChanged(app, event)
            %FLIPHORIZONTALLYCHECKBOXVALUECHANGED Flip moving source image before calibration.

            if app.HasCalibrationSource && ~isempty(app.CurrentMovingImage)
                app.CurrentMovingImage = fliplr(app.CurrentMovingImage);
                app.EnhancedMovingImage = adapthisteq(app.normalizeImage(app.CurrentMovingImage));
                app.updateRegisteredMovingImage();
                app.refreshPreview();
                app.setLocalStatus('Moving camera image flipped horizontally. Re-run automatic coregistration if needed.');
            end

        end

        % Value changed function: TransformationDropDown
        function TransformationDropDownValueChanged(app, event)
            %TRANSFORMATIONDROPDOWNVALUECHANGED Update manual-step units.

            app.updateStepSizeLabel();

        end

        % Value changed function: StepsizeEditField
        function StepsizeEditFieldValueChanged(app, event)
            %STEPSIZEEDITFIELDVALUECHANGED Validate manual adjustment step size.

            if ~isfinite(app.StepsizeEditField.Value) || app.StepsizeEditField.Value <= 0
                app.StepsizeEditField.Value = 1;
                uialert(app.UIFigure, ...
                    'Step size must be a finite positive value.', ...
                    'Invalid step size', ...
                    'Icon', 'warning');
            end

        end

        % Value changed function: ApplyStepSpinner
        function ApplyStepSpinnerValueChanged(app, event)
            %APPLYSTEPSPINNERVALUECHANGED Apply relative manual transform step.

            try
                if isempty(app.CurrentTform)
                    app.ApplyStepSpinner.Value = 0;
                    return
                end

                signedStep = event.Value - event.PreviousValue;
                if signedStep == 0
                    return
                end

                app.applyManualStep(signedStep);

            catch ME
                uialert(app.UIFigure, getReport(ME, 'basic', 'hyperlinks', 'off'), ...
                    'Manual adjustment failed', ...
                    'Icon', 'error');
            end

            app.ApplyStepSpinner.Value = 0;

        end

        % Button pushed function: SaveCalibrationButton
        function SaveCalibrationButtonPushed(app, event)
            %SAVECALIBRATIONBUTTONPUSHED Save current candidate as active calibration.

            if isempty(app.CurrentTform)
                uialert(app.UIFigure, ...
                    'No candidate calibration is available to save.', ...
                    'No calibration candidate', ...
                    'Icon', 'warning');
                return
            end

            answer = uiconfirm(app.UIFigure, ...
                ['Save the currently previewed calibration as the active calibration file?' newline newline ...
                'If an active calibration exists, it will be backed up automatically.'], ...
                'Save active calibration?', ...
                'Options', {'Save Calibration', 'Cancel'}, ...
                'DefaultOption', 'Save Calibration', ...
                'CancelOption', 'Cancel', ...
                'Icon', 'question');

            if strcmpi(answer, 'Cancel')
                return
            end

            try
                app.ensureCurrentTformInfoHasRegisteredImages();
                tform = app.CurrentTform;
                tformInfo = app.CurrentTformInfo;

                app.saveTform(tform, tformInfo);

                app.HasUnsavedCalibration = false;
                app.setCandidateTransform(app.ActiveTform, app.ActiveTformInfo, ...
                    'Active calibration', false, true);

                app.refreshAllState();
                app.setLocalStatus('Calibration saved as active coreg2cam_tform.mat.');

            catch ME
                app.setLocalStatus(sprintf('Save calibration failed: %s', ME.message));
                uialert(app.UIFigure, getReport(ME, 'basic', 'hyperlinks', 'off'), ...
                    'Save calibration failed', ...
                    'Icon', 'error');
            end

        end

        % Button pushed function: ArchiveActiveCalibrationButton
        function ArchiveActiveCalibrationButtonPushed(app, event)
            %ARCHIVEACTIVECALIBRATIONBUTTONPUSHED Backup and remove active tform.

            activeFile = fullfile(app.LabeoFolder, [app.filename '.mat']);
            if ~isfile(activeFile)
                app.refreshAllState();
                return
            end

            answer = uiconfirm(app.UIFigure, ...
                ['This will remove the active calibration used by future dual-camera imports.' newline newline ...
                'A backup copy will be kept. Future imports will not be coregistered until a new calibration is saved.'], ...
                'Archive active calibration?', ...
                'Options', {'Archive', 'Cancel'}, ...
                'DefaultOption', 'Cancel', ...
                'CancelOption', 'Cancel', ...
                'Icon', 'warning');

            if strcmpi(answer, 'Cancel')
                return
            end

            try
                app.backupActiveCalibration();
                delete(activeFile);

                app.ActiveTform = [];
                app.ActiveTformInfo = struct();

                if strcmpi(app.CandidateSource, 'Active calibration')
                    app.CurrentTform = [];
                    app.CurrentTformInfo = struct();
                    app.CandidateSource = '';
                    app.CurrentRegisteredMovingImage = [];
                    app.HasUnsavedCalibration = false;
                end

                app.pruneBackupFiles();
                app.refreshAllState();
                app.setLocalStatus('Active calibration archived. Future imports will not be coregistered until a new calibration is saved.');

            catch ME
                app.setLocalStatus(sprintf('Archive failed: %s', ME.message));
                uialert(app.UIFigure, getReport(ME, 'basic', 'hyperlinks', 'off'), ...
                    'Archive active calibration failed', ...
                    'Icon', 'error');
            end

        end

        % Button pushed function: LoadBackupButton
        function LoadBackupButtonPushed(app, event)
            %LOADBACKUPBUTTONPUSHED Load backup transform as current preview candidate.

            backupFiles = app.listBackupFiles();
            if isempty(backupFiles)
                uialert(app.UIFigure, ...
                    'No calibration backup files were found.', ...
                    'No backups available', ...
                    'Icon', 'info');
                return
            end

            [fileName, folderName] = uigetfile( ...
                fullfile(app.LabeoFolder, [app.filename '_*_bkp.mat']), ...
                'Select backup calibration file', ...
                app.LabeoFolder);

            if isequal(fileName, 0)
                return
            end

            backupPath = fullfile(folderName, fileName);

            try
                S = load(backupPath, 'tform', 'tformInfo');
                if ~isfield(S, 'tform')
                    error('Selected backup does not contain variable "tform".');
                end

                if isfield(S, 'tformInfo')
                    tformInfo = S.tformInfo;
                else
                    tformInfo = struct();
                end

                tformInfo.CandidateSource = 'Backup';
                tformInfo.SourceBackupFile = backupPath;
                tformInfo.LoadedFromBackupOn = datetime('now');

                app.setCandidateTransform(S.tform, tformInfo, ...
                    ['Backup: ' fileName], true, true);

                if app.HasCalibrationSource
                    app.setLocalStatus('Backup loaded as preview candidate. Review it, optionally fine tune small drifts, then Save Calibration to make it active.');
                else
                    app.setLocalStatus('Backup loaded as candidate. Load calibration data to preview it before saving.');
                end

            catch ME
                app.setLocalStatus(sprintf('Load backup failed: %s', ME.message));
                uialert(app.UIFigure, getReport(ME, 'basic', 'hyperlinks', 'off'), ...
                    'Load backup failed', ...
                    'Icon', 'error');
            end

        end

        % Button pushed function: ResetButton
        function ResetButtonPushed(app, event)
            %RESETBUTTONPUSHED Reset manual fine-tuning adjustment.

            if isempty(app.BaseTform)
                return
            end

            app.ManualDx = 0;
            app.ManualDy = 0;
            app.ManualRotationDeg = 0;
            app.ManualScale = 1;
            app.CurrentTform = app.BaseTform;
            app.CurrentTformInfo = app.updateTformInfoForManualAdjustment(app.CurrentTformInfo);
            app.updateRegisteredMovingImage();
            app.HasUnsavedCalibration = true;

            app.refreshCalibrationFileInfo();
            app.refreshButtonStates();
            app.refreshPreview();
            app.setLocalStatus('Manual adjustment reset. Save Calibration to keep the reset candidate.');

        end

        % Button pushed function: LoadCalibrationDataButton
        function LoadCalibrationDataButtonPushed(app, event)
            %LOADCALIBRATIONDATABUTTONPUSHED Load unregistered source images.
            %
            %   This implementation expects a folder containing AcqInfos.mat and
            %   unregistered .dat channel files. Do not select processed data that
            %   were imported while an active coreg2cam_tform.mat existed.

            startFolder = app.DataFolder;
            if isempty(startFolder) || ~isfolder(startFolder)
                startFolder = pwd;
            end

            answer = uiconfirm(app.UIFigure, ...
                ['Select a folder containing UNREGISTERED dual-camera calibration .dat files.' newline newline ...
                'Do not select processed data that may already have been coregistered during import.'], ...
                'Load calibration source data', ...
                'Options', {'Select Folder', 'Cancel'}, ...
                'DefaultOption', 'Select Folder', ...
                'CancelOption', 'Cancel', ...
                'Icon', 'warning');

            if strcmpi(answer, 'Cancel')
                return
            end

            selectedFolder = uigetdir(startFolder, 'Select unregistered calibration-source folder');
            if isequal(selectedFolder, 0)
                return
            end

            try
                app.stopAlternateViewTimer();
                app.CalibrationSourceFolder = selectedFolder;
                app.CurrentFixedImage = [];
                app.CurrentMovingImage = [];
                app.CurrentRegisteredMovingImage = [];
                app.EnhancedFixedImage = [];
                app.EnhancedMovingImage = [];
                app.HasCalibrationSource = false;

                w = uiprogressdlg(app.UIFigure, ...
                    'Title', 'Loading calibration data', ...
                    'Message', 'Building camera preview images...', ...
                    'Indeterminate', 'on');
                cleanupObj = onCleanup(@() delete(w));

                app.loadCameraImages();

                if ~isempty(app.CurrentTform)
                    app.updateRegisteredMovingImage();
                elseif ~isempty(app.ActiveTform)
                    app.setCandidateTransform(app.ActiveTform, app.ActiveTformInfo, ...
                        'Active calibration', false, true);
                end

                app.refreshCalibrationFileInfo();
                app.refreshButtonStates();
                app.refreshPreview();
                app.setLocalStatus('Calibration source loaded. Run automatic coregistration or preview the current candidate.');

            catch ME
                app.HasCalibrationSource = false;
                app.CalibrationSourceFolder = '';
                app.refreshCalibrationFileInfo();
                app.refreshButtonStates();
                app.refreshPreview();
                app.setLocalStatus(sprintf('Load calibration data failed: %s', ME.message));
                uialert(app.UIFigure, getReport(ME, 'basic', 'hyperlinks', 'off'), ...
                    'Load calibration data failed', ...
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
            app.UIFigure.Position = [100 100 1202 855];
            app.UIFigure.Name = 'Dual Camera coregistration';
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, @UIFigureCloseRequest, true);

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {280, '1x'};
            app.GridLayout.RowHeight = {'1x', 30};

            % Create UIAxes
            app.UIAxes = uiaxes(app.GridLayout);
            title(app.UIAxes, 'Title')
            app.UIAxes.XTick = [];
            app.UIAxes.YTick = [];
            app.UIAxes.Layout.Row = 1;
            app.UIAxes.Layout.Column = 2;

            % Create GrifLayoutLeft
            app.GrifLayoutLeft = uigridlayout(app.GridLayout);
            app.GrifLayoutLeft.ColumnWidth = {'1x'};
            app.GrifLayoutLeft.RowHeight = {70, 80, 80, 40, '1x', 40, '1x'};
            app.GrifLayoutLeft.Padding = [0 0 0 0];
            app.GrifLayoutLeft.Layout.Row = 1;
            app.GrifLayoutLeft.Layout.Column = 1;

            % Create RunAutomaticcoregistrationButton
            app.RunAutomaticcoregistrationButton = uibutton(app.GrifLayoutLeft, 'push');
            app.RunAutomaticcoregistrationButton.ButtonPushedFcn = createCallbackFcn(app, @RunAutomaticcoregistrationButtonPushed, true);
            app.RunAutomaticcoregistrationButton.BackgroundColor = [0.902 0.902 0.902];
            app.RunAutomaticcoregistrationButton.Layout.Row = 4;
            app.RunAutomaticcoregistrationButton.Layout.Column = 1;
            app.RunAutomaticcoregistrationButton.Text = 'Run Automatic coregistration';

            % Create CoregistrationparametersPanel
            app.CoregistrationparametersPanel = uipanel(app.GrifLayoutLeft);
            app.CoregistrationparametersPanel.Title = 'Coregistration parameters';
            app.CoregistrationparametersPanel.Layout.Row = 3;
            app.CoregistrationparametersPanel.Layout.Column = 1;
            app.CoregistrationparametersPanel.FontWeight = 'bold';

            % Create GridCoregistrationPanel
            app.GridCoregistrationPanel = uigridlayout(app.CoregistrationparametersPanel);
            app.GridCoregistrationPanel.ColumnWidth = {'1x'};

            % Create UseenhancedcontrastCheckBox
            app.UseenhancedcontrastCheckBox = uicheckbox(app.GridCoregistrationPanel);
            app.UseenhancedcontrastCheckBox.ValueChangedFcn = createCallbackFcn(app, @UseenhancedcontrastCheckBoxValueChanged, true);
            app.UseenhancedcontrastCheckBox.Tooltip = {'Check this box to enhance the image contrast to improve coregistration'};
            app.UseenhancedcontrastCheckBox.Text = 'Use enhanced contrast';
            app.UseenhancedcontrastCheckBox.Layout.Row = 1;
            app.UseenhancedcontrastCheckBox.Layout.Column = 1;
            app.UseenhancedcontrastCheckBox.Value = true;

            % Create FlipHorizontallyCheckBox
            app.FlipHorizontallyCheckBox = uicheckbox(app.GridCoregistrationPanel);
            app.FlipHorizontallyCheckBox.ValueChangedFcn = createCallbackFcn(app, @FlipHorizontallyCheckBoxValueChanged, true);
            app.FlipHorizontallyCheckBox.Tooltip = {'Some Imaging systems will have a mirrored image in camera #2. Check this box if this is the case.'};
            app.FlipHorizontallyCheckBox.Text = 'Flip Horizontally';
            app.FlipHorizontallyCheckBox.Layout.Row = 2;
            app.FlipHorizontallyCheckBox.Layout.Column = 1;

            % Create SaveCalibrationButton
            app.SaveCalibrationButton = uibutton(app.GrifLayoutLeft, 'push');
            app.SaveCalibrationButton.ButtonPushedFcn = createCallbackFcn(app, @SaveCalibrationButtonPushed, true);
            app.SaveCalibrationButton.Layout.Row = 6;
            app.SaveCalibrationButton.Layout.Column = 1;
            app.SaveCalibrationButton.Text = 'Save Calibration';

            % Create ViewmodeButtonGroup
            app.ViewmodeButtonGroup = uibuttongroup(app.GrifLayoutLeft);
            app.ViewmodeButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @ViewmodeButtonGroupSelectionChanged, true);
            app.ViewmodeButtonGroup.Title = 'View mode';
            app.ViewmodeButtonGroup.Layout.Row = 2;
            app.ViewmodeButtonGroup.Layout.Column = 1;
            app.ViewmodeButtonGroup.FontWeight = 'bold';

            % Create AlternateViewButton
            app.AlternateViewButton = uiradiobutton(app.ViewmodeButtonGroup);
            app.AlternateViewButton.Text = 'Alternate View';
            app.AlternateViewButton.Position = [11 30 99 22];
            app.AlternateViewButton.Value = true;

            % Create FalseColorOverlayButton
            app.FalseColorOverlayButton = uiradiobutton(app.ViewmodeButtonGroup);
            app.FalseColorOverlayButton.Text = 'False Color Overlay';
            app.FalseColorOverlayButton.Position = [11 8 128 22];

            % Create ManualfinetuningPanel
            app.ManualfinetuningPanel = uipanel(app.GrifLayoutLeft);
            app.ManualfinetuningPanel.Title = 'Manual fine tuning';
            app.ManualfinetuningPanel.Layout.Row = 5;
            app.ManualfinetuningPanel.Layout.Column = 1;
            app.ManualfinetuningPanel.FontWeight = 'bold';

            % Create GridManualFineTuning
            app.GridManualFineTuning = uigridlayout(app.ManualfinetuningPanel);
            app.GridManualFineTuning.ColumnWidth = {90, '1x'};
            app.GridManualFineTuning.RowHeight = {35, 35, 35, 40};

            % Create TransformationDropDownLabel
            app.TransformationDropDownLabel = uilabel(app.GridManualFineTuning);
            app.TransformationDropDownLabel.Layout.Row = 1;
            app.TransformationDropDownLabel.Layout.Column = 1;
            app.TransformationDropDownLabel.Text = 'Transformation';

            % Create TransformationDropDown
            app.TransformationDropDown = uidropdown(app.GridManualFineTuning);
            app.TransformationDropDown.Items = {'Horizontal (transl.)', 'Vertical (transl.)', 'Rotation', 'Scaling'};
            app.TransformationDropDown.ValueChangedFcn = createCallbackFcn(app, @TransformationDropDownValueChanged, true);
            app.TransformationDropDown.Layout.Row = 1;
            app.TransformationDropDown.Layout.Column = 2;
            app.TransformationDropDown.Value = 'Horizontal (transl.)';

            % Create StepsizeEditFieldLabel
            app.StepsizeEditFieldLabel = uilabel(app.GridManualFineTuning);
            app.StepsizeEditFieldLabel.Layout.Row = 2;
            app.StepsizeEditFieldLabel.Layout.Column = 1;
            app.StepsizeEditFieldLabel.Text = 'Step size';

            % Create StepsizeEditField
            app.StepsizeEditField = uieditfield(app.GridManualFineTuning, 'numeric');
            app.StepsizeEditField.ValueChangedFcn = createCallbackFcn(app, @StepsizeEditFieldValueChanged, true);
            app.StepsizeEditField.Layout.Row = 2;
            app.StepsizeEditField.Layout.Column = 2;

            % Create ApplyStepSpinnerLabel
            app.ApplyStepSpinnerLabel = uilabel(app.GridManualFineTuning);
            app.ApplyStepSpinnerLabel.Layout.Row = 3;
            app.ApplyStepSpinnerLabel.Layout.Column = 1;
            app.ApplyStepSpinnerLabel.Text = 'Apply Step';

            % Create ApplyStepSpinner
            app.ApplyStepSpinner = uispinner(app.GridManualFineTuning);
            app.ApplyStepSpinner.ValueChangedFcn = createCallbackFcn(app, @ApplyStepSpinnerValueChanged, true);
            app.ApplyStepSpinner.Layout.Row = 3;
            app.ApplyStepSpinner.Layout.Column = 2;

            % Create ResetButton
            app.ResetButton = uibutton(app.GridManualFineTuning, 'push');
            app.ResetButton.ButtonPushedFcn = createCallbackFcn(app, @ResetButtonPushed, true);
            app.ResetButton.Layout.Row = 4;
            app.ResetButton.Layout.Column = [1 2];
            app.ResetButton.Text = 'Reset';

            % Create CalibrationfilePanel
            app.CalibrationfilePanel = uipanel(app.GrifLayoutLeft);
            app.CalibrationfilePanel.Title = 'Calibration file';
            app.CalibrationfilePanel.Layout.Row = 7;
            app.CalibrationfilePanel.Layout.Column = 1;

            % Create GridGeometricPanel
            app.GridGeometricPanel = uigridlayout(app.CalibrationfilePanel);
            app.GridGeometricPanel.ColumnWidth = {'1x'};
            app.GridGeometricPanel.RowHeight = {'1x', 40, 40};

            % Create CalibrationFileInfoTextArea
            app.CalibrationFileInfoTextArea = uitextarea(app.GridGeometricPanel);
            app.CalibrationFileInfoTextArea.BackgroundColor = [0.9412 0.9412 0.9412];
            app.CalibrationFileInfoTextArea.Layout.Row = 1;
            app.CalibrationFileInfoTextArea.Layout.Column = 1;

            % Create LoadBackupButton
            app.LoadBackupButton = uibutton(app.GridGeometricPanel, 'push');
            app.LoadBackupButton.ButtonPushedFcn = createCallbackFcn(app, @LoadBackupButtonPushed, true);
            app.LoadBackupButton.Tooltip = {'Replace the curent Tform file with a backup Tform file.'};
            app.LoadBackupButton.Layout.Row = 2;
            app.LoadBackupButton.Layout.Column = 1;
            app.LoadBackupButton.Text = 'Load Backup';

            % Create ArchiveActiveCalibrationButton
            app.ArchiveActiveCalibrationButton = uibutton(app.GridGeometricPanel, 'push');
            app.ArchiveActiveCalibrationButton.ButtonPushedFcn = createCallbackFcn(app, @ArchiveActiveCalibrationButtonPushed, true);
            app.ArchiveActiveCalibrationButton.Layout.Row = 3;
            app.ArchiveActiveCalibrationButton.Layout.Column = 1;
            app.ArchiveActiveCalibrationButton.Text = 'Archive Active Calibration';

            % Create CalibrationSourcePanel
            app.CalibrationSourcePanel = uipanel(app.GrifLayoutLeft);
            app.CalibrationSourcePanel.Title = 'Calibration Source';
            app.CalibrationSourcePanel.Layout.Row = 1;
            app.CalibrationSourcePanel.Layout.Column = 1;
            app.CalibrationSourcePanel.FontWeight = 'bold';

            % Create GridCalibrationSource
            app.GridCalibrationSource = uigridlayout(app.CalibrationSourcePanel);
            app.GridCalibrationSource.ColumnWidth = {'1x'};
            app.GridCalibrationSource.RowHeight = {'1x'};
            app.GridCalibrationSource.Padding = [5 5 5 5];

            % Create LoadCalibrationDataButton
            app.LoadCalibrationDataButton = uibutton(app.GridCalibrationSource, 'push');
            app.LoadCalibrationDataButton.ButtonPushedFcn = createCallbackFcn(app, @LoadCalibrationDataButtonPushed, true);
            app.LoadCalibrationDataButton.Layout.Row = 1;
            app.LoadCalibrationDataButton.Layout.Column = 1;
            app.LoadCalibrationDataButton.Text = 'Load Calibration Data....';

            % Create StatusLabel
            app.StatusLabel = uilabel(app.GridLayout);
            app.StatusLabel.FontAngle = 'italic';
            app.StatusLabel.Layout.Row = 2;
            app.StatusLabel.Layout.Column = [1 2];
            app.StatusLabel.Text = 'Status';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = DataViewer_Coreg2Cams_exported(varargin)

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