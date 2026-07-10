classdef ImageAlignmentTool_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                       matlab.ui.Figure
        GridLayout                     matlab.ui.container.GridLayout
        GridLayoutBottom               matlab.ui.container.GridLayout
        StatusLabel                    matlab.ui.control.Label
        ApplyRegistrationtoFolderButton  matlab.ui.control.Button
        CloseButton                    matlab.ui.control.Button
        GridLayoutLeft                 matlab.ui.container.GridLayout
        RegistrationoptionsPanel       matlab.ui.container.Panel
        GridRegistrationOptions        matlab.ui.container.GridLayout
        ResetRegistrationButton        matlab.ui.control.Button
        FreehandAlignmentButton        matlab.ui.control.Button
        AutomaticcoregistrationButton  matlab.ui.control.Button
        DataSelectionPanel             matlab.ui.container.Panel
        GridDataSelectionPanel         matlab.ui.container.GridLayout
        TextArea                       matlab.ui.control.TextArea
        SelectMovingsourceButton       matlab.ui.control.Button
        SelectImageReferenceButton     matlab.ui.control.Button
        ManualfinetuningPanel          matlab.ui.container.Panel
        GridManualFineTuning           matlab.ui.container.GridLayout
        ResetFineTuningButton          matlab.ui.control.Button
        ApplyStepSpinner               matlab.ui.control.Spinner
        ApplyStepSpinnerLabel          matlab.ui.control.Label
        StepsizeEditField              matlab.ui.control.NumericEditField
        StepsizeEditFieldLabel         matlab.ui.control.Label
        TransformationDropDown         matlab.ui.control.DropDown
        TransformationDropDownLabel    matlab.ui.control.Label
        CoregistrationparametersPanel  matlab.ui.container.Panel
        GridCoregistrationPanel        matlab.ui.container.GridLayout
        FlipHorizontallyCheckBox       matlab.ui.control.CheckBox
        UseenhancedcontrastCheckBox    matlab.ui.control.CheckBox
        ViewmodeButtonGroup            matlab.ui.container.ButtonGroup
        FalseColorOverlayButton        matlab.ui.control.RadioButton
        AlternateViewButton            matlab.ui.control.RadioButton
        UIAxes                         matlab.ui.control.UIAxes
    end


    properties (Access = private)
        % Parent app/context.
        MainApp = []
        AppMode char = 'standalone'
        DataFolder char = ''
        TargetFolder char = ''

        % Fixed reference image from ImageReferenceManager.
        CurrentImageReference = []
        ReferenceFileAbs char = ''
        ReferenceFileRel char = ''
        CurrentFixedImage = []
        EnhancedFixedImage = []
        HasReference logical = false

        % Moving image source. CurrentMovingImage is the unflipped averaged
        % source image. Horizontal flipping is applied dynamically when
        % building previews/registration inputs.
        CurrentMovingImage = []
        EnhancedMovingImage = []
        MovingSourceFile char = ''
        MovingSourceFolder char = ''
        MovingSourceDescription char = ''
        MovingSourceSize double = []
        HasMovingSource logical = false

        % Current preview image after applying candidate transform to the
        % prepared moving image.
        CurrentRegisteredMovingImage = []

        % Candidate transform state.
        BaseTform = []
        CurrentTform = []
        CurrentTformInfo = struct()
        CandidateSource char = ''
        HasUnsavedRegistration logical = false
        RegistrationApplied logical = false

        % Manual adjustment state. These are applied on top of BaseTform.
        ManualDx double = 0
        ManualDy double = 0
        ManualRotationDeg double = 0
        ManualScale double = 1

        % Alternate-view timer state.
        AlternateViewTimer = timer.empty
        AlternateViewImageHandle = matlab.graphics.primitive.Image.empty
        AlternateViewIndex double = 1

        % Last backend report returned by applyImageAlignmentToFolder.
        LastApplyReport = []
        MovingSourceLocked logical = false

        % Free-hand alignment runtime state.
        IsFreehandEditing logical = false
        FreehandBoxHandle = []
        FreehandBoxListeners cell = {}
        FreehandOverlayHandle = matlab.graphics.primitive.Image.empty


    end

    properties (Access = public)
        OutputReport = []
        OutputTargetFolder char = ''
        OutputTform = []
    end

    methods (Access = private)

        function configureInitialState(app)
            %CONFIGUREINITIALSTATE Initialize UI and runtime state.

            app.UIFigure.Visible = 'off';
            app.UIFigure.WindowStyle = 'normal';
            app.UIFigure.CloseRequestFcn = @(src, evt) app.UIFigureCloseRequest(evt);

            app.TextArea.Editable = 'off';
            app.TextArea.Value = {'No ImageReference or moving source selected.'};

            app.StepsizeEditField.Limits = [eps Inf];
            app.StepsizeEditField.Value = 1;
            app.updateStepSizeLabel();

            app.ApplyStepSpinner.Limits = [-1 1];
            app.ApplyStepSpinner.Step = 1;
            app.ApplyStepSpinner.Value = 0;

            app.ViewmodeButtonGroup.SelectedObject = app.FalseColorOverlayButton;

            app.clearCandidateTransform();

            if app.MovingSourceLocked
                app.SelectMovingsourceButton.Enable = 'off';
            end


            cla(app.UIAxes);
            title(app.UIAxes, 'Select ImageReference and moving source', 'Interpreter', 'none');
            app.UIAxes.XTick = [];
            app.UIAxes.YTick = [];
            axis(app.UIAxes, 'image');
            colormap(app.UIAxes, gray);

            app.setLocalStatus('Ready. Select an ImageReference and a moving source.');
        end

        function parseStartupOptions(app, nvArgs)
            %PARSESTARTUPOPTIONS Parse optional name-value startup inputs.
            %
            %   nvArgs is the varargin cell from startupFcn. Do not expand it when
            %   calling this method.

            if nargin < 2 || isempty(nvArgs)
                nvArgs = {};
            end

            p = inputParser;
            p.FunctionName = 'ImageAlignmentTool.startupFcn';

            addParameter(p, 'movingImage', [], ...
                @(x) isempty(x) || isnumeric(x) || islogical(x));

            addParameter(p, 'movingSourceFile', '', ...
                @(x) ischar(x) || (isstring(x) && isscalar(x)));

            addParameter(p, 'movingSourceDescription', '', ...
                @(x) ischar(x) || (isstring(x) && isscalar(x)));

            addParameter(p, 'targetFolder', '', ...
                @(x) ischar(x) || (isstring(x) && isscalar(x)));

            addParameter(p, 'lockMovingSource', false, ...
                @(x) islogical(x) && isscalar(x));

            parse(p, nvArgs{:});
            R = p.Results;

            if ~isempty(R.targetFolder)
                app.TargetFolder = char(string(R.targetFolder));
            end

            app.MovingSourceLocked = R.lockMovingSource;

            if ~isempty(R.movingImage)
                sourceFile = char(string(R.movingSourceFile));
                sourceFolder = app.TargetFolder;

                if ~isempty(sourceFile)
                    sourceFolder = fileparts(sourceFile);
                end

                if isempty(app.TargetFolder) && ~isempty(sourceFolder)
                    app.TargetFolder = sourceFolder;
                end

                description = char(string(R.movingSourceDescription));

                if isempty(description)
                    description = 'Startup moving image';
                end

                app.setMovingSourceFromImage( ...
                    R.movingImage, ...
                    sourceFile, ...
                    sourceFolder, ...
                    description);
            end
        end

        function setMovingSourceFromImage(app, movingImage, sourceFile, sourceFolder, sourceDescription)
            %SETMOVINGSOURCEFROMIMAGE Store averaged moving image and context.

            movingImage = app.validate2DImage(movingImage, 'movingImage');

            app.CurrentMovingImage = single(movingImage);
            app.CurrentMovingImage(~isfinite(app.CurrentMovingImage)) = 0;
            app.MovingSourceSize = size(app.CurrentMovingImage);

            app.EnhancedMovingImage = app.enhanceImageForRegistration( ...
                app.normalizeImage(app.CurrentMovingImage));

            app.MovingSourceFile = char(string(sourceFile));
            app.MovingSourceFolder = char(string(sourceFolder));
            app.MovingSourceDescription = char(string(sourceDescription));

            app.TargetFolder = app.MovingSourceFolder;

            app.HasMovingSource = true;
            app.clearCandidateTransform();
        end

        function img = readDatAverageImage(app, datFile, maxFramesToAverage)
            %READDATAVERAGEIMAGE Read a .dat file as an averaged 2D moving image.
            %
            %   The function averages up to maxFramesToAverage evenly spaced frames.
            %   It first tries DatImageSource and falls back to mapDat for compatibility.

            if nargin < 3 || isempty(maxFramesToAverage)
                maxFramesToAverage = 100;
            end

            if ~isfile(datFile)
                error('ImageAlignmentTool:MissingDatFile', ...
                    'DAT file not found: %s', datFile);
            end

            if exist('DatImageSource', 'class') == 8
                src = DatImageSource(char(datFile));
                sz = src.getSize();

                if numel(sz) < 3
                    nFrames = 1;
                else
                    nFrames = double(sz(3));
                end

                frameIdx = app.makeAverageFrameIndex(nFrames, maxFramesToAverage);

                img = [];
                n = 0;

                for iFrame = 1:numel(frameIdx)
                    thisFrame = single(squeeze(src.getFrame(frameIdx(iFrame))));

                    if isempty(img)
                        img = zeros(size(thisFrame), 'single');
                    elseif ~isequal(size(img), size(thisFrame))
                        error('ImageAlignmentTool:InconsistentDatFrameSize', ...
                            'DAT frames do not have consistent spatial size.');
                    end

                    n = n + 1;
                    img = img + (thisFrame - img) ./ n;
                end

                img = app.validate2DImage(img, 'dat average image');
                return
            end

            if exist('mapDat', 'file') ~= 2
                error('ImageAlignmentTool:MissingDatReader', ...
                    'Neither DatImageSource nor mapDat is available on the MATLAB path.');
            end

            datMap = mapDat(char(datFile));
            nFrames = size(datMap.Data.data, 3);
            frameIdx = app.makeAverageFrameIndex(nFrames, maxFramesToAverage);

            img = mean(single(datMap.Data.data(:, :, frameIdx)), 3, 'omitnan');
            img = app.validate2DImage(img, 'dat average image');
        end

        function frameIdx = makeAverageFrameIndex(app, nFrames, maxFramesToAverage) %#ok<INUSL>
            %MAKEAVERAGEFRAMEINDEX Return evenly spaced frame indices.

            nFrames = max(1, round(double(nFrames)));
            nRead = min(nFrames, maxFramesToAverage);

            if nFrames > nRead
                frameIdx = unique(round(linspace(1, nFrames, nRead)));
            else
                frameIdx = 1:nFrames;
            end
        end

        function setCandidateTransform(app, tform, tformInfo, sourceLabel, isUnsaved, resetManual)
            %SETCANDIDATETRANSFORM Store candidate transform and refresh state.

            if nargin < 6
                resetManual = true;
            end

            app.BaseTform = tform;
            app.CurrentTform = tform;

            if isempty(tformInfo) || ~isstruct(tformInfo)
                tformInfo = struct();
            end

            app.CurrentTformInfo = tformInfo;
            app.CandidateSource = char(string(sourceLabel));
            app.HasUnsavedRegistration = logical(isUnsaved);
            app.RegistrationApplied = false;

            if resetManual
                app.ManualDx = 0;
                app.ManualDy = 0;
                app.ManualRotationDeg = 0;
                app.ManualScale = 1;
            end

            app.updateRegisteredMovingImage();
            app.refreshDataSelectionInfo();
            app.refreshButtonStates();
            app.refreshPreview();
        end

        function clearCandidateTransform(app)
            %CLEARCANDIDATETRANSFORM Clear transform and manual adjustment state.

            app.stopFreehandAlignmentRuntime();

            app.BaseTform = [];
            app.CurrentTform = [];
            app.CurrentTformInfo = struct();
            app.CandidateSource = '';
            app.CurrentRegisteredMovingImage = [];
            app.HasUnsavedRegistration = false;
            app.RegistrationApplied = false;

            app.ManualDx = 0;
            app.ManualDy = 0;
            app.ManualRotationDeg = 0;
            app.ManualScale = 1;
        end

        function startFreehandAlignment(app)
            %STARTFREEHANDALIGNMENT Enter interactive rectangle-based alignment mode.

            if ~app.HasReference || isempty(app.CurrentFixedImage)
                uialert(app.UIFigure, ...
                    'Select an ImageReference before using free-hand alignment.', ...
                    'Missing ImageReference', ...
                    'Icon', 'warning');
                return
            end

            if ~app.HasMovingSource || isempty(app.CurrentMovingImage)
                uialert(app.UIFigure, ...
                    'Select a moving source before using free-hand alignment.', ...
                    'Missing moving source', ...
                    'Icon', 'warning');
                return
            end

            app.stopAlternateViewTimer();
            app.stopFreehandAlignmentRuntime();

            app.IsFreehandEditing = true;
            app.ViewmodeButtonGroup.SelectedObject = app.FalseColorOverlayButton;

            [initialPosition, initialRotationDeg] = app.getInitialFreehandBoxGeometry();

            initialTform = app.calculateFreehandTformFromGeometry( ...
                initialPosition, initialRotationDeg);

            app.CurrentTform = initialTform;
            app.BaseTform = initialTform;
            app.CandidateSource = 'Free-hand alignment preview';
            app.HasUnsavedRegistration = true;
            app.RegistrationApplied = false;

            app.ManualDx = 0;
            app.ManualDy = 0;
            app.ManualRotationDeg = 0;
            app.ManualScale = 1;

            app.updateRegisteredMovingImage();

            % Draw overlay first, then set image-bounded manual limits. Do not use
            % padded limits here. Users can zoom out manually when they need access
            % to rectangle handles outside the visible image area.
            app.renderFreehandPreview();
            app.setFreehandAxesToImageLimits();

            app.FreehandBoxHandle = drawrectangle(app.UIAxes, ...
                'Position', initialPosition, ...
                'Rotatable', true, ...
                'DrawingArea','unlimited',...
                'InteractionsAllowed', 'all');

            if isprop(app.FreehandBoxHandle, 'RotationAngle')
                app.FreehandBoxHandle.RotationAngle = initialRotationDeg;
            end

            if isprop(app.FreehandBoxHandle, 'FixedAspectRatio')
                app.FreehandBoxHandle.FixedAspectRatio = true;
            end

            if isprop(app.FreehandBoxHandle, 'Color')
                app.FreehandBoxHandle.Color = [1 1 0];
            end

            if isprop(app.FreehandBoxHandle, 'FaceAlpha')
                app.FreehandBoxHandle.FaceAlpha = 0;
            end

            app.addFreehandBoxListeners();

            app.refreshDataSelectionInfo();
            app.refreshButtonStates();
            app.setLocalStatus('Free-hand alignment active. Zoom out if needed, then move/rotate/resize the box and click Confirm.');
        end

        function confirmFreehandAlignment(app)
            %CONFIRMFREEHANDALIGNMENT Commit current free-hand rectangle transform.

            if ~app.IsFreehandEditing || ~app.isUsableFreehandBox()
                error('ImageAlignmentTool:InvalidFreehandState', ...
                    'No active free-hand alignment box is available.');
            end

            freehandTform = app.calculateFreehandTformFromBox();

            tformInfo = struct();
            tformInfo.Method = 'freehand';
            tformInfo.TransformType = 'similarity';
            tformInfo.GeneratedOn = datetime('now');
            tformInfo.UseEnhancedContrast = logical(app.UseenhancedcontrastCheckBox.Value);
            tformInfo.PreFlipHorizontally = logical(app.FlipHorizontallyCheckBox.Value);
            tformInfo.ReferenceFile = app.ReferenceFileRel;
            tformInfo.MovingSourceFile = app.MovingSourceFile;
            tformInfo.TargetFolder = app.TargetFolder;
            tformInfo.FreehandBoxPosition = app.FreehandBoxHandle.Position;

            if isprop(app.FreehandBoxHandle, 'RotationAngle')
                tformInfo.FreehandRotationDeg = app.FreehandBoxHandle.RotationAngle;
            else
                tformInfo.FreehandRotationDeg = 0;
            end

            app.stopFreehandAlignmentRuntime();

            app.setCandidateTransform( ...
                freehandTform, ...
                tformInfo, ...
                'Free-hand similarity alignment', ...
                true, ...
                true);

            app.setLocalStatus('Free-hand alignment confirmed. Review or fine tune if needed.');
        end

        function stopFreehandAlignmentRuntime(app)
            %STOPFREEHANDALIGNMENTRUNTIME Delete free-hand ROI/listeners/overlay.

            for iListener = 1:numel(app.FreehandBoxListeners)
                try
                    delete(app.FreehandBoxListeners{iListener});
                catch
                end
            end

            app.FreehandBoxListeners = {};

            try
                if ~isempty(app.FreehandBoxHandle) && isvalid(app.FreehandBoxHandle)
                    delete(app.FreehandBoxHandle);
                end
            catch
            end

            app.FreehandBoxHandle = [];

            try
                if ~isempty(app.FreehandOverlayHandle) && isvalid(app.FreehandOverlayHandle)
                    delete(app.FreehandOverlayHandle);
                end
            catch
            end

            app.FreehandOverlayHandle = matlab.graphics.primitive.Image.empty;

            app.IsFreehandEditing = false;

            try
                app.FreehandAlignmentButton.Text = 'Free-hand Alignment';
            catch
            end
        end

        function addFreehandBoxListeners(app)
            %ADDFREEHANDBOXLISTENERS Attach ROI listeners for live free-hand preview.

            app.FreehandBoxListeners = {};

            if ~app.isUsableFreehandBox()
                return
            end

            try
                app.FreehandBoxListeners{end+1} = addlistener( ...
                    app.FreehandBoxHandle, ...
                    'MovingROI', ...
                    @(src, evt) app.onFreehandBoxChanged(src, evt));
            catch
            end

            try
                app.FreehandBoxListeners{end+1} = addlistener( ...
                    app.FreehandBoxHandle, ...
                    'ROIMoved', ...
                    @(src, evt) app.onFreehandBoxChanged(src, evt));
            catch
            end
        end

        function onFreehandBoxChanged(app, src, evt) %#ok<INUSD>
            %ONFREEHANDBOXCHANGED Update preview while user edits the rectangle.

            if ~app.IsFreehandEditing || ~app.isUsableFreehandBox()
                return
            end

            try
                app.CurrentTform = app.calculateFreehandTformFromBox();
                app.BaseTform = app.CurrentTform;
                app.CandidateSource = 'Free-hand alignment preview';
                app.HasUnsavedRegistration = true;
                app.RegistrationApplied = false;

                app.updateRegisteredMovingImage();
                app.renderFreehandPreview();

                drawnow limitrate

            catch ME
                app.setLocalStatus(sprintf('Free-hand preview update failed: %s', ME.message));
            end
        end

        function renderFreehandPreview(app)
            %RENDERFREEHANDPREVIEW Draw/update false-color overlay without clearing ROI.
            %
            %   Important:
            %       Do not reset XLim/YLim during live updates. The user may zoom out
            %       during free-hand editing to reach rectangle edges/handles.

            if ~app.HasReference || isempty(app.CurrentFixedImage) || ...
                    ~app.HasMovingSource || isempty(app.CurrentMovingImage)
                return
            end

            fixedImg = app.getDisplayFixedImage();

            if isempty(app.CurrentRegisteredMovingImage)
                movingImg = zeros(size(fixedImg), 'single');
            else
                movingImg = app.CurrentRegisteredMovingImage;
            end

            overlayImg = imfuse(fixedImg, movingImg, ...
                'falsecolor', ...
                'Scaling', 'joint');

            if isempty(app.FreehandOverlayHandle) || ~isvalid(app.FreehandOverlayHandle)
                cla(app.UIAxes);
                app.FreehandOverlayHandle = image(app.UIAxes, overlayImg);
                app.UIAxes.YDir = 'reverse';
                app.UIAxes.DataAspectRatio = [1 1 1];
                app.UIAxes.XTick = [];
                app.UIAxes.YTick = [];
                title(app.UIAxes, 'Free-hand alignment: adjust moving-image box', 'Interpreter', 'none');
            else
                app.FreehandOverlayHandle.CData = overlayImg;
            end
        end

        function setFreehandAxesToImageLimits(app)
            %SETFREEHANDAXESTOIMAGELIMITS Initialize free-hand axes to image bounds.
            %
            %   The limits are manual so user zoom/pan interactions can change them.
            %   The app should not repeatedly reset them while the user is editing.

            fixedSize = size(app.CurrentFixedImage);

            if numel(fixedSize) < 2
                return
            end

            fixedH = fixedSize(1);
            fixedW = fixedSize(2);

            app.UIAxes.XLim = [0.5, fixedW + 0.5];
            app.UIAxes.YLim = [0.5, fixedH + 0.5];
            app.UIAxes.XLimMode = 'manual';
            app.UIAxes.YLimMode = 'manual';
            app.UIAxes.DataAspectRatio = [1 1 1];
            app.UIAxes.YDir = 'reverse';
        end

        function tform = calculateFreehandTformFromBox(app)
            %CALCULATEFREEHANDTFORMFROMBOX Return transform represented by ROI box.

            if ~app.isUsableFreehandBox()
                error('ImageAlignmentTool:InvalidFreehandBox', ...
                    'Free-hand alignment box is missing or invalid.');
            end

            position = app.FreehandBoxHandle.Position;

            if isprop(app.FreehandBoxHandle, 'RotationAngle')
                rotationDeg = app.FreehandBoxHandle.RotationAngle;
            else
                rotationDeg = 0;
            end

            tform = app.calculateFreehandTformFromGeometry(position, rotationDeg);
        end

        function tform = calculateFreehandTformFromGeometry(app, position, rotationDeg)
            %CALCULATEFREEHANDTFORMFROMGEOMETRY Convert box geometry to similarity transform.
            %
            %   The rectangle represents where the prepared moving image should land
            %   in reference-image coordinates.
            %
            %   The ROI RotationAngle visual convention is opposite to the row-vector
            %   affine matrix convention used here with YDir = reverse. Therefore,
            %   the sign is inverted when building the affine transform.
            %
            %   Row-vector convention:
            %       [x y 1] * TtoOrigin * S * R * TtoEditedCenter

            movingSize = size(app.CurrentMovingImage);
            movingH = movingSize(1);
            movingW = movingSize(2);

            if numel(position) ~= 4 || any(~isfinite(position)) || ...
                    position(3) <= 0 || position(4) <= 0
                error('ImageAlignmentTool:InvalidFreehandPosition', ...
                    'Free-hand rectangle position must be [x y width height].');
            end

            oldCenter = [(movingW + 1) / 2, (movingH + 1) / 2];
            newCenter = [position(1) + position(3) / 2, position(2) + position(4) / 2];

            scaleX = position(3) / movingW;
            scaleY = position(4) / movingH;
            scale = mean([scaleX scaleY]);

            if abs(scaleX - scaleY) > 1e-3 * max(1, abs(scale))
                error('ImageAlignmentTool:NonUniformFreehandScale', ...
                    'Free-hand rectangle scaling is not uniform. Keep the rectangle aspect ratio fixed.');
            end

            % Sign flip fixes the visual mismatch:
            % turning the ROI box clockwise should rotate the moving image clockwise.
            theta = deg2rad(-rotationDeg);

            TtoOrigin = [1 0 0; 0 1 0; -oldCenter(1) -oldCenter(2) 1];
            S = [scale 0 0; 0 scale 0; 0 0 1];
            R = [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];
            TtoEditedCenter = [1 0 0; 0 1 0; newCenter(1) newCenter(2) 1];

            tform = affine2d(TtoOrigin * S * R * TtoEditedCenter);
        end

        function [position, rotationDeg] = getInitialFreehandBoxGeometry(app)
            %GETINITIALFREEHANDBOXGEOMETRY Initialize free-hand box from current candidate.
            %
            %   If a candidate transform exists, initialize the box from the current
            %   transformed moving-image footprint. Otherwise, use the unregistered
            %   prepared moving image footprint.
            %
            %   calculateFreehandTformFromGeometry applies a sign flip to convert ROI
            %   RotationAngle into row-vector affine rotation. Therefore this method
            %   applies the inverse sign when converting an existing transform back
            %   into the ROI rectangle angle.

            movingSize = size(app.CurrentMovingImage);
            movingH = movingSize(1);
            movingW = movingSize(2);

            if isempty(app.CurrentTform)
                T = eye(3);
            else
                T = app.getTransformMatrix(app.CurrentTform);
            end

            corners = [ ...
                0.5,          0.5; ...
                movingW+0.5,  0.5; ...
                movingW+0.5,  movingH+0.5; ...
                0.5,          movingH+0.5];

            [x, y] = transformPointsForward(affine2d(T), corners(:,1), corners(:,2));
            pts = [x(:), y(:)];

            edgeW = pts(2,:) - pts(1,:);
            edgeH = pts(4,:) - pts(1,:);

            boxW = norm(edgeW);
            boxH = norm(edgeH);

            if ~isfinite(boxW) || ~isfinite(boxH) || boxW <= 0 || boxH <= 0
                position = [0.5 0.5 movingW movingH];
                rotationDeg = 0;
                return
            end

            centerXY = mean(pts, 1);

            % Invert sign so existing transform -> ROI angle -> transform is stable.
            rotationDeg = -atan2d(edgeW(2), edgeW(1));

            position = [centerXY(1) - boxW/2, centerXY(2) - boxH/2, boxW, boxH];
        end

        function tf = isUsableFreehandBox(app)
            %ISUSABLEFREEHANDBOX True if free-hand ROI object is usable.

            tf = false;

            try
                tf = ~isempty(app.FreehandBoxHandle) && ...
                    isvalid(app.FreehandBoxHandle) && ...
                    isprop(app.FreehandBoxHandle, 'Position') && ...
                    ~isempty(app.FreehandBoxHandle.Position);
            catch
                tf = false;
            end
        end

        function updateRegisteredMovingImage(app)
            %UPDATEREGISTEREDMOVINGIMAGE Apply current transform to moving preview image.

            app.CurrentRegisteredMovingImage = [];

            if isempty(app.CurrentTform) || isempty(app.CurrentMovingImage) || ...
                    isempty(app.CurrentFixedImage)
                return
            end

            movingImg = app.getPreparedMovingImage();
            outRef = imref2d(size(app.CurrentFixedImage));

            registered = imwarp(movingImg, app.CurrentTform, ...
                'OutputView', outRef, ...
                'FillValues', 0, ...
                'InterpolationMethod', 'linear');

            app.CurrentRegisteredMovingImage = app.normalizeImage(registered);
        end

        function fixedImg = getRegistrationFixedImage(app)
            %GETREGISTRATIONFIXEDIMAGE Return fixed image used by imregtform.

            if app.UseenhancedcontrastCheckBox.Value && ~isempty(app.EnhancedFixedImage)
                fixedImg = app.EnhancedFixedImage;
            else
                fixedImg = app.normalizeImage(app.CurrentFixedImage);
            end
        end

        function movingImg = getRegistrationMovingImage(app)
            %GETREGISTRATIONMOVINGIMAGE Return moving image used by imregtform.

            movingImg = app.getPreparedMovingImage();

            if app.UseenhancedcontrastCheckBox.Value
                movingImg = app.enhanceImageForRegistration(movingImg);
            else
                movingImg = app.normalizeImage(movingImg);
            end
        end

        function movingImg = getPreparedMovingImage(app)
            %GETPREPAREDMOVINGIMAGE Return normalized moving image with optional flip.

            movingImg = app.normalizeImage(app.CurrentMovingImage);

            if app.FlipHorizontallyCheckBox.Value
                movingImg = fliplr(movingImg);
            end
        end

        function img = getDisplayFixedImage(app)
            %GETDISPLAYFIXEDIMAGE Return fixed image for preview.

            img = app.getRegistrationFixedImage();
        end

        function img = getDisplayMovingImage(app)
            %GETDISPLAYMOVINGIMAGE Return moving or registered moving image for preview.

            if ~isempty(app.CurrentRegisteredMovingImage)
                img = app.CurrentRegisteredMovingImage;
            else
                img = app.getPreparedMovingImage();
            end

            if app.UseenhancedcontrastCheckBox.Value
                img = app.enhanceImageForRegistration(img);
            else
                img = app.normalizeImage(img);
            end
        end

        function refreshPreview(app)
            %REFRESHPREVIEW Draw current preview mode.

            if app.IsFreehandEditing
                app.renderFreehandPreview();
                return
            end

            app.stopAlternateViewTimer();
            cla(app.UIAxes);

            hasFixed = app.HasReference && ~isempty(app.CurrentFixedImage);
            hasMoving = app.HasMovingSource && ~isempty(app.CurrentMovingImage);

            if ~hasFixed && ~hasMoving
                title(app.UIAxes, 'Select ImageReference and moving source', 'Interpreter', 'none');
                app.formatPreviewAxes();
                return
            end

            if hasFixed && ~hasMoving
                imagesc(app.UIAxes, app.getDisplayFixedImage(), [0 1]);
                colormap(app.UIAxes, gray);
                title(app.UIAxes, 'Selected ImageReference', 'Interpreter', 'none');
                app.formatPreviewAxes();
                return
            end

            if ~hasFixed && hasMoving
                imagesc(app.UIAxes, app.getDisplayMovingImage(), [0 1]);
                colormap(app.UIAxes, gray);
                title(app.UIAxes, 'Selected moving source', 'Interpreter', 'none');
                app.formatPreviewAxes();
                return
            end

            fixedImg = app.getDisplayFixedImage();
            movingImg = app.getDisplayMovingImage();

            if isempty(app.CurrentTform) && ~isequal(size(fixedImg), size(movingImg))
                imagesc(app.UIAxes, movingImg, [0 1]);
                colormap(app.UIAxes, gray);
                title(app.UIAxes, ...
                    'Moving source loaded. Run automatic alignment to overlay with reference.', ...
                    'Interpreter', 'none');
                app.formatPreviewAxes();
                return
            end

            if startsWith(app.ViewmodeButtonGroup.SelectedObject.Text, 'False', 'IgnoreCase', true)
                imshowpair(fixedImg, movingImg, ...
                    'falsecolor', ...
                    'Scaling', 'joint', ...
                    'Parent', app.UIAxes);
                title(app.UIAxes, app.getPreviewTitle('False-color overlay'), 'Interpreter', 'none');
                app.formatPreviewAxes();
            else
                app.startAlternateView(fixedImg, movingImg);
            end
        end

        function formatPreviewAxes(app)
            %FORMATPREVIEWAXES Apply common preview axis formatting.

            axis(app.UIAxes, 'image');
            app.UIAxes.XTick = [];
            app.UIAxes.YTick = [];
        end

        function titleText = getPreviewTitle(app, modeText)
            %GETPREVIEWTITLE Compose preview title.

            if isempty(app.CandidateSource)
                candidateText = 'unregistered moving source';
            else
                candidateText = app.CandidateSource;
            end

            titleText = sprintf('%s | %s', modeText, candidateText);
        end

        function startAlternateView(app, fixedImg, movingImg)
            %STARTALTERNATEVIEW Start 1 Hz flicker comparison on UIAxes.

            app.stopAlternateViewTimer();
            app.AlternateViewIndex = 1;

            app.AlternateViewImageHandle = imagesc(app.UIAxes, fixedImg, [0 1]);
            colormap(app.UIAxes, gray);
            app.formatPreviewAxes();
            title(app.UIAxes, app.getPreviewTitle('Alternate view: reference'), 'Interpreter', 'none');

            app.AlternateViewTimer = timer( ...
                'ExecutionMode', 'fixedRate', ...
                'Period', 1, ...
                'BusyMode', 'drop', ...
                'Name', 'ImageAlignmentToolAlternateViewTimer', ...
                'TimerFcn', @(~, ~) app.onAlternateViewTimerTick(fixedImg, movingImg));

            start(app.AlternateViewTimer);
        end

        function onAlternateViewTimerTick(app, fixedImg, movingImg)
            %ONALTERNATEVIEWTIMERTICK Alternate fixed and moving images.

            if isempty(app) || ~isvalid(app) || isempty(app.UIFigure) || ~isvalid(app.UIFigure)
                return
            end

            if isempty(app.AlternateViewImageHandle) || ~isvalid(app.AlternateViewImageHandle)
                return
            end

            if app.AlternateViewIndex == 1
                app.AlternateViewImageHandle.CData = movingImg;
                app.AlternateViewIndex = 2;
                title(app.UIAxes, app.getPreviewTitle('Alternate view: moving'), 'Interpreter', 'none');
            else
                app.AlternateViewImageHandle.CData = fixedImg;
                app.AlternateViewIndex = 1;
                title(app.UIAxes, app.getPreviewTitle('Alternate view: reference'), 'Interpreter', 'none');
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
            app.AlternateViewImageHandle = matlab.graphics.primitive.Image.empty;
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
                error('ImageAlignmentTool:InvalidStepSize', ...
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

            app.HasUnsavedRegistration = true;
            app.CandidateSource = app.appendManualCandidateLabel(app.CandidateSource);

            app.updateRegisteredMovingImage();
            app.refreshDataSelectionInfo();
            app.refreshButtonStates();
            app.refreshPreview();

            app.setLocalStatus('Manual adjustment applied.');
        end

        function label = appendManualCandidateLabel(app, labelIn) %#ok<INUSL>
            %APPENDMANUALCANDIDATELABEL Mark candidate as manually adjusted.

            label = char(string(labelIn));

            if isempty(label)
                label = 'Manual adjustment';
            elseif ~contains(label, '+ manual adjustment')
                label = [label '+ manual adjustment'];
            end
        end

        function recomposeCurrentTransform(app)
            %RECOMPOSECURRENTTRANSFORM Compose base transform plus manual adjustment.
            %
            %   BaseTform maps the prepared moving image into reference space.
            %
            %   Manual rotation/scaling are applied around the moving image center in
            %   moving-image coordinates. Manual translation is applied afterward in
            %   reference-image coordinates so horizontal/vertical nudges behave like
            %   screen-space shifts in the overlay.

            if isempty(app.BaseTform)
                return
            end

            movingSize = size(app.CurrentMovingImage);
            cx = (movingSize(2) + 1) / 2;
            cy = (movingSize(1) + 1) / 2;

            % Row-vector transform convention: [x y 1] * T.
            TtoOrigin = [1 0 0; 0 1 0; -cx -cy 1];
            TfromOrigin = [1 0 0; 0 1 0; cx cy 1];

            theta = deg2rad(app.ManualRotationDeg);
            R = [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];
            S = [app.ManualScale 0 0; 0 app.ManualScale 0; 0 0 1];

            % Rotation/scaling around the moving-image center.
            manualMovingT = TtoOrigin * S * R * TfromOrigin;

            % Translation remains in reference/output coordinates for intuitive nudging.
            manualReferenceTranslationT = [1 0 0; 0 1 0; app.ManualDx app.ManualDy 1];

            baseT = app.getTransformMatrix(app.BaseTform);

            % Apply moving-centered rotation/scale before the automatic transform,
            % then apply the translation in reference coordinates.
            app.CurrentTform = affine2d(manualMovingT * baseT * manualReferenceTranslationT);

            app.CurrentTformInfo = app.updateTformInfoForManualAdjustment(app.CurrentTformInfo);
        end

        function info = updateTformInfoForManualAdjustment(app, info)
            %UPDATETFORMINFOFORMANUALADJUSTMENT Add manual-adjustment metadata.

            if isempty(info) || ~isstruct(info)
                info = struct();
            end

            info.ManualAdjustment = app.getManualAdjustmentStruct();
            info.CandidateSource = app.CandidateSource;
            info.UpdatedOn = datetime('now');
        end

        function s = getManualAdjustmentStruct(app)
            %GETMANUALADJUSTMENTSTRUCT Return current manual adjustment metadata.

            s = struct( ...
                'dx_px', app.ManualDx, ...
                'dy_px', app.ManualDy, ...
                'rotation_deg', app.ManualRotationDeg, ...
                'scale', app.ManualScale, ...
                'updatedOn', datetime('now'));
        end

        function finalTform = getFinalTformForRawData(app)
            %GETFINALTFORMFORRAWDATA Return transform from raw moving data to reference.
            %
            %   CurrentTform maps the prepared moving image to reference space. If
            %   horizontal flip is enabled, the prepared moving image is fliplr(raw).
            %   Therefore, the final raw-data transform must include the flip before
            %   the candidate transform.
            %
            %   This method uses getTransformMatrix so it works with simtform2d,
            %   affinetform2d, affine2d, and numeric 3 x 3 matrices.

            if isempty(app.CurrentTform)
                finalTform = [];
                return
            end

            currentT = app.getTransformMatrix(app.CurrentTform);

            if ~app.FlipHorizontallyCheckBox.Value
                finalTform = affine2d(currentT);
                return
            end

            movingSize = size(app.CurrentMovingImage);
            Nx = movingSize(2);

            % Row-vector transform convention: [x y 1] * T.
            Tflip = [-1 0 0; 0 1 0; Nx + 1 0 1];

            finalTform = affine2d(Tflip * currentT);
        end

        function refreshButtonStates(app)
            %REFRESHBUTTONSTATES Enable controls from current app state.

            hasReference = app.HasReference && ~isempty(app.CurrentFixedImage);
            hasMoving = app.HasMovingSource && ~isempty(app.CurrentMovingImage);
            hasCandidate = ~isempty(app.CurrentTform);
            hasTargetFolder = ~isempty(app.TargetFolder) && isfolder(app.TargetFolder);

            if app.IsFreehandEditing
                app.SelectImageReferenceButton.Enable = 'off';
                app.SelectMovingsourceButton.Enable = 'off';
                app.AutomaticcoregistrationButton.Enable = 'off';
                app.FreehandAlignmentButton.Enable = 'on';
                app.FreehandAlignmentButton.Text = 'Confirm';
                app.ResetRegistrationButton.Enable = 'on';
                app.ApplyRegistrationtoFolderButton.Enable = 'off';
                app.setManualControlsEnabled(false);
                return
            end

            app.SelectImageReferenceButton.Enable = 'on';
            app.SelectMovingsourceButton.Enable = app.onOff(~app.MovingSourceLocked);

            app.AutomaticcoregistrationButton.Enable = app.onOff(hasReference && hasMoving);
            app.FreehandAlignmentButton.Enable = app.onOff(hasReference && hasMoving);
            app.FreehandAlignmentButton.Text = 'Free-hand Alignment';

            app.ResetRegistrationButton.Enable = app.onOff(hasCandidate);
            app.ApplyRegistrationtoFolderButton.Enable = app.onOff(hasCandidate && hasTargetFolder);

            app.setManualControlsEnabled(hasCandidate);
        end

        function setManualControlsEnabled(app, tf)
            %SETMANUALCONTROLSENABLED Enable or disable manual fine-tuning controls.

            state = app.onOff(tf);

            app.TransformationDropDown.Enable = state;
            app.StepsizeEditField.Enable = state;
            app.ApplyStepSpinner.Enable = state;
            app.ResetFineTuningButton.Enable = state;
        end

        function refreshDataSelectionInfo(app)
            %REFRESHDATASELECTIONINFO Update data-selection summary text.

            lines = strings(0, 1);

            lines(end+1) = "Image Alignment Tool";
            lines(end+1) = "";

            if app.HasReference && ~isempty(app.CurrentImageReference)
                lines(end+1) = "ImageReference: selected";

                if isfield(app.CurrentImageReference, 'name')
                    lines(end+1) = "Name: " + string(app.CurrentImageReference.name);
                end

                if ~isempty(app.ReferenceFileRel)
                    lines(end+1) = "Reference file: " + string(app.ReferenceFileRel);
                elseif ~isempty(app.ReferenceFileAbs)
                    lines(end+1) = "Reference file: " + string(app.ReferenceFileAbs);
                end

                lines(end+1) = "Reference size: " + app.sizeToText(size(app.CurrentFixedImage));
            else
                lines(end+1) = "ImageReference: not selected";
            end

            lines(end+1) = "";

            if app.HasMovingSource
                lines(end+1) = "Moving source: selected";
                lines(end+1) = "Source: " + string(app.MovingSourceDescription);

                if app.MovingSourceLocked
                    lines(end+1) = "Moving source is locked by launch context.";
                end

                if ~isempty(app.MovingSourceFile)
                    lines(end+1) = "File: " + string(app.MovingSourceFile);
                end

                lines(end+1) = "Moving size: " + app.sizeToText(size(app.CurrentMovingImage));
            else
                lines(end+1) = "Moving source: not selected";
            end

            lines(end+1) = "";

            if isempty(app.TargetFolder)
                lines(end+1) = "Target folder: not set";
            else
                lines(end+1) = "Target folder:";
                lines(end+1) = string(app.TargetFolder);
            end

            lines(end+1) = "";

            lines = [lines, app.getTransformSummaryLines()];

            if app.RegistrationApplied
                lines(end+1) = "";
                lines(end+1) = "Registration applied: yes";
            end

            app.TextArea.Value = cellstr(lines);
        end

        function img = enhanceImageForRegistration(app, img)
            %ENHANCEIMAGEFORREGISTRATION Normalize and contrast-enhance one image.

            img = app.normalizeImage(img);

            try
                img = adapthisteq(img);
            catch
                img = imadjust(img);
            end

            img = app.normalizeImage(img);
        end

        function img = normalizeImage(app, img) %#ok<INUSL>
            %NORMALIZEIMAGE Normalize image to [0 1] safely.

            img = single(img);
            img(~isfinite(img)) = 0;

            imgMin = min(img(:), [], 'omitnan');
            imgMax = max(img(:), [], 'omitnan');

            if ~isfinite(imgMin) || ~isfinite(imgMax) || imgMax <= imgMin
                img = zeros(size(img), 'single');
            else
                img = (img - imgMin) ./ (imgMax - imgMin);
                img = min(max(img, 0), 1);
            end
        end

        function img = validate2DImage(app, img, varName) %#ok<INUSL>
            %VALIDATE2DIMAGE Validate one 2D numeric/logical image.

            if ~(isnumeric(img) || islogical(img)) || ndims(img) ~= 2 || isempty(img)
                error('ImageAlignmentTool:InvalidImage', ...
                    '%s must be a non-empty 2D numeric or logical image.', varName);
            end

            img = single(img);
        end

        function txt = sizeToText(app, sz) %#ok<INUSL>
            %SIZETOTEXT Return compact [Y X] size text.

            if isempty(sz) || numel(sz) < 2
                txt = "[]";
            else
                txt = sprintf('%d x %d', sz(1), sz(2));
            end
        end

        function state = onOff(app, tf) %#ok<INUSL>
            %ONOFF Convert logical to MATLAB on/off state char.

            if tf
                state = 'on';
            else
                state = 'off';
            end
        end

        function value = getObjectPropertyIfPresent(app, obj, propName, defaultValue) %#ok<INUSL>
            %GETOBJECTPROPERTYIFPRESENT Return property value or default.

            value = defaultValue;

            try
                if isprop(obj, propName)
                    value = obj.(propName);
                end
            catch
                value = defaultValue;
            end
        end

        function deleteAppIfValid(app, appObj) %#ok<INUSL>
            %DELETEAPPIFVALID Delete app object if still valid.

            try
                if ~isempty(appObj) && isvalid(appObj)
                    delete(appObj);
                end
            catch
            end
        end

        function deleteIfValid(app, obj) %#ok<INUSL>
            %DELETEIFVALID Delete graphics/UI object if valid.

            try
                if ~isempty(obj) && isvalid(obj)
                    delete(obj);
                end
            catch
            end
        end

        function setLocalStatus(app, msg)
            %SETLOCALSTATUS Update status label text.

            app.StatusLabel.Text = char(string(msg));
            drawnow limitrate
        end

        function cleanupImageAlignmentTool(app)
            %CLEANUPIMAGEALIGNMENTTOOL Stop transient runtime resources.

            try
                app.stopAlternateViewTimer();
            catch
            end

            try
                app.stopFreehandAlignmentRuntime();
            catch
            end
        end

        function requestClose(app)
            %REQUESTCLOSE Close the tool after optional confirmation.

            if app.HasUnsavedRegistration && ~app.RegistrationApplied
                answer = uiconfirm(app.UIFigure, ...
                    'There is an unapplied registration candidate. Close without applying it?', ...
                    'Unapplied registration', ...
                    'Options', {'Close without applying', 'Cancel'}, ...
                    'DefaultOption', 'Cancel', ...
                    'CancelOption', 'Cancel', ...
                    'Icon', 'warning');

                if strcmpi(answer, 'Cancel')
                    return
                end
            end

            app.cleanupImageAlignmentTool();
            delete(app);
        end

        function info = get2DTransformInfo(app, tform)
            %GET2DTRANSFORMINFO Extract compact transform info.
            %
            %   Works with older and newer MATLAB transform objects through
            %   getTransformMatrix.

            if isempty(tform)
                info = struct();
                info.dx = NaN;
                info.dy = NaN;
                info.rotationDeg = NaN;
                info.scale = NaN;
                info.scaleX = NaN;
                info.scaleY = NaN;
                info.determinant = NaN;
                info.hasReflection = false;
                return
            end

            T = app.getTransformMatrix(tform);
            A = T(1:2, 1:2);

            scaleX = hypot(A(1,1), A(1,2));
            scaleY = hypot(A(2,1), A(2,2));

            info = struct();
            info.dx = T(3,1);
            info.dy = T(3,2);
            info.rotationDeg = atan2d(A(1,2), A(1,1));
            info.scale = mean([scaleX scaleY]);
            info.scaleX = scaleX;
            info.scaleY = scaleY;
            info.determinant = det(A);
            info.hasReflection = info.determinant < 0;
        end

        function lines = getTransformSummaryLines(app)
            %GETTRANSFORMSUMMARYLINES Return automatic/manual/final transform summary.
            %
            %   The summary distinguishes:
            %       - automatic transform estimated from the prepared moving image
            %       - manual delta applied on top of the automatic transform
            %       - final raw-data transform sent to the backend

            lines = strings(0, 1);

            if isempty(app.CurrentTform)
                lines(end+1) = "Candidate transform: none";
                return
            end

            lines(end+1) = "Candidate transform: ready";
            lines(end+1) = "Transform type: similarity";

            if ~isempty(app.CandidateSource)
                lines(end+1) = "Source: " + string(app.CandidateSource);
            end

            lines(end+1) = "";

            if ~isempty(app.BaseTform)
                baseInfo = app.get2DTransformInfo(app.BaseTform);
                lines(end+1) = sprintf( ...
                    'Automatic prepared-image transform: dx=%.3g px, dy=%.3g px, rot=%.3g deg, scale=%.6g', ...
                    baseInfo.dx, baseInfo.dy, baseInfo.rotationDeg, baseInfo.scale);
            else
                lines(end+1) = "Automatic prepared-image transform: none";
            end

            lines(end+1) = sprintf( ...
                'Manual delta: dx=%.3g px, dy=%.3g px, rot=%.3g deg, scale=%.6g', ...
                app.ManualDx, app.ManualDy, app.ManualRotationDeg, app.ManualScale);

            finalTform = app.getFinalTformForRawData();
            finalInfo = app.get2DTransformInfo(finalTform);

            if finalInfo.hasReflection
                reflectionText = "yes";
            else
                reflectionText = "no";
            end

            lines(end+1) = sprintf( ...
                'Final raw-data transform: dx=%.3g px, dy=%.3g px, rot=%.3g deg, scale=%.6g, reflection=%s', ...
                finalInfo.dx, finalInfo.dy, finalInfo.rotationDeg, finalInfo.scale, reflectionText);

            if app.FlipHorizontallyCheckBox.Value
                lines(end+1) = "Final raw-data transform includes horizontal pre-flip.";
            end
        end

        function T = getTransformMatrix(app, tform) %#ok<INUSL>
            %GETTRANSFORMMATRIX Return 3 x 3 transform matrix.

            if isnumeric(tform)
                T = tform;
                return
            end

            if isprop(tform, 'T')
                T = tform.T;
                return
            end

            if isprop(tform, 'A')
                T = tform.A;
                return
            end

            error('ImageAlignmentTool:MissingTransformMatrix', ...
                'Could not extract transform matrix from class "%s".', class(tform));
        end

        function impact = getAlignmentImpactSummary(app)
            %GETALIGNMENTIMPACTSUMMARY Summarize files found before registration.
            %
            %   This is a GUI-side folder listing. The backend may still skip
            %   non-image UMT files or unsupported legacy ROI files after schema
            %   validation, so the wording intentionally says "files found" rather
            %   than "files transformed".

            targetFolder = app.TargetFolder;

            datFiles = dir(fullfile(targetFolder, '*.dat'));
            umtFiles = dir(fullfile(targetFolder, '*.umt'));
            roiFiles = dir(fullfile(targetFolder, '*.roi'));

            impact = struct();
            impact.datFiles = string({datFiles.name}).';
            impact.umtFiles = string({umtFiles.name}).';
            impact.roiFiles = string({roiFiles.name}).';

            impact.lines = strings(0, 1);
            impact.lines(end+1) = "Files found in target folder:";

            impact.lines(end+1) = sprintf("  DAT files: %d", numel(datFiles));
            impact.lines = [impact.lines, app.previewFileList(impact.datFiles, 6)];

            impact.lines(end+1) = sprintf("  UMT files: %d", numel(umtFiles));
            impact.lines = [impact.lines, app.previewFileList(impact.umtFiles, 6)];

            impact.lines(end+1) = sprintf("  ROI files: %d", numel(roiFiles));
            impact.lines = [impact.lines, app.previewFileList(impact.roiFiles, 6)];
        end

        function lines = previewFileList(app, fileList, maxItems) %#ok<INUSL>
            %PREVIEWFILELIST Return compact indented file list.

            lines = strings(0, 1);

            if isempty(fileList)
                lines(end+1) = "    <none>";
                return
            end

            nShow = min(numel(fileList), maxItems);

            for iFile = 1:nShow
                lines(end+1) = "    - " + fileList(iFile);
            end

            nMore = numel(fileList) - nShow;
            if nMore > 0
                lines(end+1) = sprintf("    ... and %d more", nMore);
            end
        end

        function tf = confirmRoiImpactBeforeApply(app, impact)
            %CONFIRMROIIMPACTBEFOREAPPLY Warn about ROI effects before transformation.
            %
            %   This method avoids accessing DataViewer private properties directly.
            %   Until DataViewer exposes getLoadedROICountForExternalTool(), loaded
            %   in-memory ROIs are reported as zero.

            tf = true;

            hasRoiFiles = isfield(impact, 'roiFiles') && ~isempty(impact.roiFiles);
            loadedRoiCount = app.getLoadedRoiCountForWarning();
            hasLoadedRois = loadedRoiCount > 0;

            if ~hasRoiFiles && ~hasLoadedRois
                return
            end

            msg = sprintf([ ...
                'This registration applies a geometric transform to the image data.\n\n' ...
                'ROI impact:\n' ...
                '- ROI files in the target folder will be transformed when supported.\n' ...
                '- Unsupported legacy ROI files will be skipped by the backend with a warning.\n' ...
                '- ROIs currently loaded in DataViewer may become stale unless DataViewer reloads after registration.\n\n' ...
                'ROI files found in target folder: %d\n' ...
                'Current DataViewer ROIs loaded: %d\n\n' ...
                'Continue?'], ...
                numel(impact.roiFiles), loadedRoiCount);

            choice = uiconfirm(app.UIFigure, msg, ...
                'ROI geometry will be affected', ...
                'Options', {'Continue registration', 'Cancel'}, ...
                'DefaultOption', 'Cancel', ...
                'CancelOption', 'Cancel', ...
                'Icon', 'warning');

            tf = strcmp(choice, 'Continue registration');
        end

        function n = getLoadedRoiCountForWarning(app)
            %GETLOADEDROICOUNTFORWARNING Return number of ROIs loaded in parent app.
            %
            %   The ImageAlignmentTool must not access DataViewer private ROIList
            %   directly. DataViewer will later expose getLoadedROICountForExternalTool().
            %   Until then, this method returns zero.

            n = 0;

            try
                if isempty(app.MainApp) || ~isvalid(app.MainApp)
                    return
                end

                if ismethod(app.MainApp, 'getLoadedROICountForExternalTool')
                    n = app.MainApp.getLoadedROICountForExternalTool();
                end
            catch
                n = 0;
            end
        end

        function refreshParentAfterSuccessfulRegistration(app)
            %REFRESHPARENTAFTERSUCCESSFULREGISTRATION Refresh launching DataViewer.
            %
            %   This method is intentionally safe before DataViewer is patched. Once
            %   DataViewer exposes refreshAfterImageAlignment(targetFolder), this
            %   method will call it automatically.

            try
                if isempty(app.MainApp) || ~isvalid(app.MainApp)
                    return
                end

                if ismethod(app.MainApp, 'refreshAfterImageAlignment')
                    app.MainApp.refreshAfterImageAlignment(app.TargetFolder);
                end

            catch ME
                warning('ImageAlignmentTool:ParentRefreshFailed', ...
                    'Registration succeeded, but parent app refresh failed: %s', ME.message);
            end
        end

    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, parentApp, dataFolder, appMode, varargin)
            %STARTUPFCN Initialize ImageAlignmentTool.
            %
            %   ImageAlignmentTool()
            %   ImageAlignmentTool(parentApp, dataFolder, appMode)
            %   ImageAlignmentTool(parentApp, dataFolder, appMode, 'Name', Value, ...)
            %
            %   This startup function receives all constructor inputs as varargin.
            %   parseStartupOptions receives the varargin cell directly.


            app.UIFigure.Visible = 'off';


            if nargin < 2
                parentApp = [];
            end

            if nargin < 3 || isempty(dataFolder)
                dataFolder = '';
            end

            if nargin < 4 || isempty(appMode)
                appMode = 'standalone';
            end

            % Store the fixed launch context here.
            app.MainApp = parentApp;
            app.DataFolder = char(string(dataFolder));
            app.AppMode = char(string(appMode));

            if ~isempty(app.DataFolder) && isfolder(app.DataFolder)
                app.TargetFolder = app.DataFolder;
            else
                app.TargetFolder = '';
            end

            app.parseStartupOptions(varargin);

            app.configureInitialState();

            app.refreshDataSelectionInfo();
            app.refreshButtonStates();
            app.refreshPreview();


            if ~isempty(app.MainApp)
                placeAppInsideCaller(app.MainApp, app.UIFigure, 'center');
            end



            app.UIFigure.Visible = 'on';

        end

        % Button pushed function: SelectImageReferenceButton
        function SelectImageReferenceButtonPushed(app, event)
            %SELECTIMAGEREFERENCEBUTTONPUSHED Select fixed ImageReference.

            try
                if exist('ImageReferenceManager', 'file') ~= 2 && ...
                        exist('ImageReferenceManager', 'class') ~= 8
                    error('ImageAlignmentTool:MissingImageReferenceManager', ...
                        'ImageReferenceManager was not found on the MATLAB path.');
                end

                refApp = ImageReferenceManager(app, 'select');
                cleanupObj = onCleanup(@() app.deleteAppIfValid(refApp));

                uiwait(refApp.UIFigure);

                if ~isvalid(refApp) || ~refApp.WasSelectionConfirmed || ...
                        isempty(refApp.OutputImageReference)
                    app.setLocalStatus('ImageReference selection cancelled.');
                    return
                end

                ImageReference = refApp.OutputImageReference;

                if ~isfield(ImageReference, 'image') || isempty(ImageReference.image)
                    error('ImageAlignmentTool:InvalidImageReference', ...
                        'The selected ImageReference does not contain an image.');
                end

                app.CurrentImageReference = ImageReference;
                app.ReferenceFileAbs = refApp.OutputReferenceFileAbs;
                app.ReferenceFileRel = refApp.OutputReferenceFileRel;

                app.CurrentFixedImage = app.normalizeImage(ImageReference.image);
                app.EnhancedFixedImage = app.enhanceImageForRegistration(app.CurrentFixedImage);
                app.HasReference = true;

                app.clearCandidateTransform();
                app.refreshDataSelectionInfo();
                app.refreshButtonStates();
                app.refreshPreview();

                app.setLocalStatus('ImageReference selected.');

            catch ME
                app.setLocalStatus(sprintf('ImageReference selection failed: %s', ME.message));
                uialert(app.UIFigure, ...
                    getReport(ME, 'basic', 'hyperlinks', 'off'), ...
                    'ImageReference selection failed', ...
                    'Icon', 'error');
            end

        end

        % Button pushed function: SelectMovingsourceButton
        function SelectMovingsourceButtonPushed(app, event)
            %SELECTMOVINGSOURCEBUTTONPUSHED Select a .dat moving source.

            startFolder = app.TargetFolder;

            if isempty(startFolder) || ~isfolder(startFolder)
                if ~isempty(app.DataFolder) && isfolder(app.DataFolder)
                    startFolder = app.DataFolder;
                else
                    startFolder = pwd;
                end
            end

            [fileName, folderName] = uigetfile( ...
                {'*.dat', 'DAT image time series (*.dat)'}, ...
                'Select moving source .dat file', ...
                startFolder);

            if isequal(fileName, 0)
                app.setLocalStatus('Moving source selection cancelled.');
                return
            end

            datFile = fullfile(folderName, fileName);

            try
                w = uiprogressdlg(app.UIFigure, ...
                    'Title', 'Loading moving source', ...
                    'Message', 'Building averaged moving image...', ...
                    'Indeterminate', 'on');
                cleanupObj = onCleanup(@() app.deleteIfValid(w));

                movingImage = app.readDatAverageImage(datFile, 100);

                app.setMovingSourceFromImage( ...
                    movingImage, ...
                    datFile, ...
                    folderName, ...
                    sprintf('DAT average: %s', fileName));

                app.refreshDataSelectionInfo();
                app.refreshButtonStates();
                app.refreshPreview();

                app.setLocalStatus('Moving source loaded.');
                
            catch ME
                app.setLocalStatus(sprintf('Moving source loading failed: %s', ME.message));
                uialert(app.UIFigure, ...
                    getReport(ME, 'basic', 'hyperlinks', 'off'), ...
                    'Moving source loading failed', ...
                    'Icon', 'error');
            end
            % Set focus back to figure
            figure(app.UIFigure);
        end

        % Selection changed function: ViewmodeButtonGroup
        function ViewmodeButtonGroupSelectionChanged(app, event)
            %VIEWMODEBUTTONGROUPSELECTIONCHANGED Refresh preview mode.

            app.refreshPreview();

        end

        % Value changed function: UseenhancedcontrastCheckBox
        function UseenhancedcontrastCheckBoxValueChanged(app, event)
            %USEENHANCEDCONTRASTCHECKBOXVALUECHANGED Refresh preview contrast.

            app.refreshPreview();


        end

        % Value changed function: FlipHorizontallyCheckBox
        function FlipHorizontallyCheckBoxValueChanged(app, event)
            %FLIPHORIZONTALLYCHECKBOXVALUECHANGED Update moving-source preparation.
            %
            %   The candidate transform is cleared if it exists because the current
            %   transform was estimated from the previous moving-image orientation.

            if ~isempty(app.CurrentTform)
                app.clearCandidateTransform();
                app.setLocalStatus('Horizontal flip changed. Candidate transform cleared; run automatic alignment again.');
            else
                app.setLocalStatus('Horizontal flip setting changed.');
            end

            app.refreshDataSelectionInfo();
            app.refreshButtonStates();
            app.refreshPreview();


        end

        % Button pushed function: AutomaticcoregistrationButton
        function AutomaticcoregistrationButtonPushed(app, event)
            %RUNAUTOMATICCOREGISTRATIONBUTTONPUSHED Estimate similarity transform.

            if ~app.HasReference || isempty(app.CurrentFixedImage)
                uialert(app.UIFigure, ...
                    'Select an ImageReference before running automatic alignment.', ...
                    'Missing ImageReference', ...
                    'Icon', 'warning');
                return
            end

            if ~app.HasMovingSource || isempty(app.CurrentMovingImage)
                uialert(app.UIFigure, ...
                    'Select a moving source before running automatic alignment.', ...
                    'Missing moving source', ...
                    'Icon', 'warning');
                return
            end

            try
                w = uiprogressdlg(app.UIFigure, ...
                    'Title', 'Automatic alignment', ...
                    'Message', 'Estimating similarity transform...', ...
                    'Indeterminate', 'on');
                cleanupObj = onCleanup(@() app.deleteIfValid(w));

                fixedImg = app.getRegistrationFixedImage();
                movingImg = app.getRegistrationMovingImage();

                [optimizer, metric] = imregconfig('multimodal');

                try
                    optimizer.MaximumIterations = 300;
                catch
                end

                try
                    optimizer.InitialRadius = optimizer.InitialRadius / 3;
                catch
                end

                tform = imregtform(movingImg, fixedImg, ...
                    'similarity', ...
                    optimizer, ...
                    metric, ...
                    'PyramidLevels', 3);

                tformInfo = struct();
                tformInfo.Method = 'imregtform';
                tformInfo.TransformType = 'similarity';
                tformInfo.GeneratedOn = datetime('now');
                tformInfo.UseEnhancedContrast = logical(app.UseenhancedcontrastCheckBox.Value);
                tformInfo.PreFlipHorizontally = logical(app.FlipHorizontallyCheckBox.Value);
                tformInfo.ReferenceFile = app.ReferenceFileRel;
                tformInfo.MovingSourceFile = app.MovingSourceFile;
                tformInfo.TargetFolder = app.TargetFolder;

                app.setCandidateTransform(tform, tformInfo, ...
                    'Automatic similarity alignment', true, true);

                app.setLocalStatus('Automatic similarity alignment completed. Review and fine tune if needed.');

            catch ME
                app.setLocalStatus(sprintf('Automatic alignment failed: %s', ME.message));
                uialert(app.UIFigure, ...
                    getReport(ME, 'basic', 'hyperlinks', 'off'), ...
                    'Automatic alignment failed', ...
                    'Icon', 'error');
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

            app.updateStepSizeLabel();


        end

        % Value changed function: ApplyStepSpinner
        function ApplyStepSpinnerValueChanged(app, event)
            %APPLYSTEPSPINNERVALUECHANGED Apply one relative manual transform step.

            try
                if isempty(app.CurrentTform)
                    app.ApplyStepSpinner.Value = 0;
                    return
                end

                if isprop(event, 'PreviousValue')
                    signedStep = event.Value - event.PreviousValue;
                else
                    signedStep = app.ApplyStepSpinner.Value;
                end

                if signedStep == 0
                    return
                end

                app.applyManualStep(signedStep);

            catch ME
                uialert(app.UIFigure, ...
                    getReport(ME, 'basic', 'hyperlinks', 'off'), ...
                    'Manual adjustment failed', ...
                    'Icon', 'error');
            end

            app.ApplyStepSpinner.Value = 0;


        end

        % Button pushed function: ResetFineTuningButton
        function ResetFineTuningButtonPushed(app, event)
            %RESETBUTTONPUSHED Reset manual adjustment to automatic/base transform.

            if isempty(app.BaseTform)
                return
            end

            app.ManualDx = 0;
            app.ManualDy = 0;
            app.ManualRotationDeg = 0;
            app.ManualScale = 1;

            app.CurrentTform = app.BaseTform;
            app.CurrentTformInfo = app.updateTformInfoForManualAdjustment(app.CurrentTformInfo);
            app.HasUnsavedRegistration = true;

            app.updateRegisteredMovingImage();
            app.refreshDataSelectionInfo();
            app.refreshButtonStates();
            app.refreshPreview();

            app.setLocalStatus('Manual adjustment reset.');

        end

        % Button pushed function: ApplyRegistrationtoFolderButton
        function ApplyRegistrationtoFolderButtonPushed(app, event)
            %APPLYREGISTRATIONTOFOLDERBUTTONPUSHED Apply candidate registration to folder.
            %
            %   This callback validates GUI state and delegates the destructive
            %   folder-level work to applyImageAlignmentToFolder.

            if isempty(app.CurrentTform)
                uialert(app.UIFigure, ...
                    'No candidate transform is available. Run automatic alignment first.', ...
                    'No registration candidate', ...
                    'Icon', 'warning');
                return
            end

            if isempty(app.TargetFolder) || ~isfolder(app.TargetFolder)
                uialert(app.UIFigure, ...
                    'Target folder is missing or invalid.', ...
                    'Invalid target folder', ...
                    'Icon', 'warning');
                return
            end

            if exist('applyImageAlignmentToFolder', 'file') ~= 2
                uialert(app.UIFigure, ...
                    ['The backend function applyImageAlignmentToFolder.m was not found ' ...
                    'on the MATLAB path. No files were modified.'], ...
                    'Backend not found', ...
                    'Icon', 'error');
                return
            end

            if ~app.HasReference || isempty(app.CurrentImageReference) || ...
                    ~isfield(app.CurrentImageReference, 'image') || ...
                    isempty(app.CurrentImageReference.image)
                uialert(app.UIFigure, ...
                    'The selected ImageReference is missing or invalid.', ...
                    'Invalid ImageReference', ...
                    'Icon', 'warning');
                return
            end

            referenceSize = size(app.CurrentFixedImage);
            finalTform = app.getFinalTformForRawData();

            if isempty(finalTform)
                uialert(app.UIFigure, ...
                    'Could not build final transform for raw data.', ...
                    'Invalid transform', ...
                    'Icon', 'error');
                return
            end

            impact = app.getAlignmentImpactSummary();

            if ~app.confirmRoiImpactBeforeApply(impact)
                app.setLocalStatus('Registration cancelled.');
                return
            end

            messageText = sprintf([ ...
                'Apply the current registration to all supported files in this folder?\n\n' ...
                '%s\n\n' ...
                'Output image size will be: %d x %d\n\n' ...
                '%s\n\n' ...
                'A managed-folder backup will be created before files are modified.'], ...
                app.TargetFolder, referenceSize(1), referenceSize(2), ...
                strjoin(cellstr(impact.lines), newline));

            answer = uiconfirm(app.UIFigure, ...
                messageText, ...
                'Apply registration to folder?', ...
                'Options', {'Apply Registration', 'Cancel'}, ...
                'DefaultOption', 'Cancel', ...
                'CancelOption', 'Cancel', ...
                'Icon', 'warning');

            if strcmpi(answer, 'Cancel')
                app.setLocalStatus('Registration cancelled.');
                return
            end

            try
                w = uiprogressdlg(app.UIFigure, ...
                    'Title', 'Applying registration', ...
                    'Message', 'Creating backup and transforming folder data...', ...
                    'Indeterminate', 'on');
                cleanupObj = onCleanup(@() app.deleteIfValid(w)); %#ok<NASGU>

                report = applyImageAlignmentToFolder( ...
                    app.TargetFolder, ...
                    finalTform, ...
                    'imageReference', app.CurrentImageReference, ...
                    'referenceFile', app.ReferenceFileRel, ...
                    'referenceImage', app.CurrentImageReference.image, ...
                    'baseTform', app.BaseTform, ...
                    'previewTform', app.CurrentTform, ...
                    'manualAdjustment', app.getManualAdjustmentStruct(), ...
                    'movingSourceFile', app.MovingSourceFile, ...
                    'preFlipHorizontally', app.FlipHorizontallyCheckBox.Value, ...
                    'createdBy', 'ImageAlignmentTool');

                app.LastApplyReport = report;
                app.OutputReport = report;
                app.OutputTargetFolder = app.TargetFolder;
                app.OutputTform = finalTform;
                app.RegistrationApplied = true;
                app.HasUnsavedRegistration = false;

                app.refreshDataSelectionInfo();
                app.refreshButtonStates();

                app.refreshParentAfterSuccessfulRegistration();

                app.setLocalStatus('Registration applied to folder.');

                uialert(app.UIFigure, ...
                    'Registration was applied successfully.', ...
                    'Registration complete', ...
                    'Icon', 'success');

            catch ME
                app.setLocalStatus(sprintf('Apply registration failed: %s', ME.message));
                uialert(app.UIFigure, ...
                    getReport(ME, 'basic', 'hyperlinks', 'off'), ...
                    'Apply registration failed', ...
                    'Icon', 'error');
            end

        end

        % Button pushed function: CloseButton
        function CloseButtonPushed(app, event)
            %CLOSEBUTTONPUSHED Close the tool.

            app.requestClose();

        end

        % Close request function: UIFigure
        function UIFigureCloseRequest(app, event)
            %UIFIGURECLOSEREQUEST Stop timers and close the tool.

            app.requestClose();


        end

        % Button pushed function: FreehandAlignmentButton
        function FreehandAlignmentButtonPushed(app, event)
            %FREEHANDALIGNMENTBUTTONPUSHED Toggle interactive free-hand alignment.

            try
                if app.IsFreehandEditing
                    app.confirmFreehandAlignment();
                else
                    app.startFreehandAlignment();
                end

            catch ME
                app.setLocalStatus(sprintf('Free-hand alignment failed: %s', ME.message));
                uialert(app.UIFigure, ...
                    getReport(ME, 'basic', 'hyperlinks', 'off'), ...
                    'Free-hand alignment failed', ...
                    'Icon', 'error');
            end

        end

        % Button pushed function: ResetRegistrationButton
        function ResetRegistrationButtonPushed(app, event)
            %RESETREGISTRATIONBUTTONPUSHED Clear automatic/free-hand registration.
            %
            %   This resets the full registration candidate. It is different from
            %   Reset Fine Tuning, which keeps BaseTform and only resets manual nudges.

            app.clearCandidateTransform();

            app.refreshDataSelectionInfo();
            app.refreshButtonStates();
            app.refreshPreview();

            app.setLocalStatus('Registration candidate reset.');

        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1300 944];
            app.UIFigure.Name = 'Image Alignment Tool';
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, @UIFigureCloseRequest, true);
            app.UIFigure.WindowStyle = 'modal';

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {280, '1x'};
            app.GridLayout.RowHeight = {'1x', 40};

            % Create UIAxes
            app.UIAxes = uiaxes(app.GridLayout);
            title(app.UIAxes, 'Title')
            app.UIAxes.XTick = [];
            app.UIAxes.YTick = [];
            app.UIAxes.Layout.Row = 1;
            app.UIAxes.Layout.Column = 2;

            % Create GridLayoutLeft
            app.GridLayoutLeft = uigridlayout(app.GridLayout);
            app.GridLayoutLeft.ColumnWidth = {'1x'};
            app.GridLayoutLeft.RowHeight = {'1x', 80, 80, 200, 185};
            app.GridLayoutLeft.Layout.Row = 1;
            app.GridLayoutLeft.Layout.Column = 1;

            % Create ViewmodeButtonGroup
            app.ViewmodeButtonGroup = uibuttongroup(app.GridLayoutLeft);
            app.ViewmodeButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @ViewmodeButtonGroupSelectionChanged, true);
            app.ViewmodeButtonGroup.Title = 'View mode';
            app.ViewmodeButtonGroup.Layout.Row = 2;
            app.ViewmodeButtonGroup.Layout.Column = 1;

            % Create AlternateViewButton
            app.AlternateViewButton = uiradiobutton(app.ViewmodeButtonGroup);
            app.AlternateViewButton.Text = 'Alternate View';
            app.AlternateViewButton.Position = [11 30 99 22];
            app.AlternateViewButton.Value = true;

            % Create FalseColorOverlayButton
            app.FalseColorOverlayButton = uiradiobutton(app.ViewmodeButtonGroup);
            app.FalseColorOverlayButton.Text = 'False Color Overlay';
            app.FalseColorOverlayButton.Position = [11 8 128 22];

            % Create CoregistrationparametersPanel
            app.CoregistrationparametersPanel = uipanel(app.GridLayoutLeft);
            app.CoregistrationparametersPanel.Title = 'Coregistration parameters';
            app.CoregistrationparametersPanel.Layout.Row = 3;
            app.CoregistrationparametersPanel.Layout.Column = 1;

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

            % Create ManualfinetuningPanel
            app.ManualfinetuningPanel = uipanel(app.GridLayoutLeft);
            app.ManualfinetuningPanel.Title = 'Manual fine tuning';
            app.ManualfinetuningPanel.Layout.Row = 5;
            app.ManualfinetuningPanel.Layout.Column = 1;

            % Create GridManualFineTuning
            app.GridManualFineTuning = uigridlayout(app.ManualfinetuningPanel);
            app.GridManualFineTuning.ColumnWidth = {90, '1x'};
            app.GridManualFineTuning.RowHeight = {30, 30, 30, 40};
            app.GridManualFineTuning.RowSpacing = 5;

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

            % Create ResetFineTuningButton
            app.ResetFineTuningButton = uibutton(app.GridManualFineTuning, 'push');
            app.ResetFineTuningButton.ButtonPushedFcn = createCallbackFcn(app, @ResetFineTuningButtonPushed, true);
            app.ResetFineTuningButton.Layout.Row = 4;
            app.ResetFineTuningButton.Layout.Column = [1 2];
            app.ResetFineTuningButton.Text = 'Reset Fine Tuning';

            % Create DataSelectionPanel
            app.DataSelectionPanel = uipanel(app.GridLayoutLeft);
            app.DataSelectionPanel.Title = 'Data Selection';
            app.DataSelectionPanel.Layout.Row = 1;
            app.DataSelectionPanel.Layout.Column = 1;
            app.DataSelectionPanel.FontWeight = 'bold';

            % Create GridDataSelectionPanel
            app.GridDataSelectionPanel = uigridlayout(app.DataSelectionPanel);
            app.GridDataSelectionPanel.ColumnWidth = {'1x'};
            app.GridDataSelectionPanel.RowHeight = {40, 40, '1x'};

            % Create SelectImageReferenceButton
            app.SelectImageReferenceButton = uibutton(app.GridDataSelectionPanel, 'push');
            app.SelectImageReferenceButton.ButtonPushedFcn = createCallbackFcn(app, @SelectImageReferenceButtonPushed, true);
            app.SelectImageReferenceButton.Layout.Row = 1;
            app.SelectImageReferenceButton.Layout.Column = 1;
            app.SelectImageReferenceButton.Text = 'Select Image Reference ...';

            % Create SelectMovingsourceButton
            app.SelectMovingsourceButton = uibutton(app.GridDataSelectionPanel, 'push');
            app.SelectMovingsourceButton.ButtonPushedFcn = createCallbackFcn(app, @SelectMovingsourceButtonPushed, true);
            app.SelectMovingsourceButton.Layout.Row = 2;
            app.SelectMovingsourceButton.Layout.Column = 1;
            app.SelectMovingsourceButton.Text = 'Select Moving source...';

            % Create TextArea
            app.TextArea = uitextarea(app.GridDataSelectionPanel);
            app.TextArea.BackgroundColor = [0.9412 0.9412 0.9412];
            app.TextArea.Placeholder = 'Data Selection Info';
            app.TextArea.Layout.Row = 3;
            app.TextArea.Layout.Column = 1;

            % Create RegistrationoptionsPanel
            app.RegistrationoptionsPanel = uipanel(app.GridLayoutLeft);
            app.RegistrationoptionsPanel.Title = 'Registration options';
            app.RegistrationoptionsPanel.Layout.Row = 4;
            app.RegistrationoptionsPanel.Layout.Column = 1;
            app.RegistrationoptionsPanel.FontWeight = 'bold';

            % Create GridRegistrationOptions
            app.GridRegistrationOptions = uigridlayout(app.RegistrationoptionsPanel);
            app.GridRegistrationOptions.ColumnWidth = {'1x'};
            app.GridRegistrationOptions.RowHeight = {50, 50, 40};

            % Create AutomaticcoregistrationButton
            app.AutomaticcoregistrationButton = uibutton(app.GridRegistrationOptions, 'push');
            app.AutomaticcoregistrationButton.ButtonPushedFcn = createCallbackFcn(app, @AutomaticcoregistrationButtonPushed, true);
            app.AutomaticcoregistrationButton.Layout.Row = 1;
            app.AutomaticcoregistrationButton.Layout.Column = 1;
            app.AutomaticcoregistrationButton.Text = 'Automatic coregistration';

            % Create FreehandAlignmentButton
            app.FreehandAlignmentButton = uibutton(app.GridRegistrationOptions, 'push');
            app.FreehandAlignmentButton.ButtonPushedFcn = createCallbackFcn(app, @FreehandAlignmentButtonPushed, true);
            app.FreehandAlignmentButton.Layout.Row = 2;
            app.FreehandAlignmentButton.Layout.Column = 1;
            app.FreehandAlignmentButton.Text = 'Free-hand Alignment';

            % Create ResetRegistrationButton
            app.ResetRegistrationButton = uibutton(app.GridRegistrationOptions, 'push');
            app.ResetRegistrationButton.ButtonPushedFcn = createCallbackFcn(app, @ResetRegistrationButtonPushed, true);
            app.ResetRegistrationButton.Layout.Row = 3;
            app.ResetRegistrationButton.Layout.Column = 1;
            app.ResetRegistrationButton.Text = 'Reset Registration';

            % Create GridLayoutBottom
            app.GridLayoutBottom = uigridlayout(app.GridLayout);
            app.GridLayoutBottom.ColumnWidth = {'1x', '1x', '1x', '1x', 180, 80};
            app.GridLayoutBottom.RowHeight = {'1x'};
            app.GridLayoutBottom.Padding = [0 0 0 0];
            app.GridLayoutBottom.Layout.Row = 2;
            app.GridLayoutBottom.Layout.Column = [1 2];

            % Create CloseButton
            app.CloseButton = uibutton(app.GridLayoutBottom, 'push');
            app.CloseButton.ButtonPushedFcn = createCallbackFcn(app, @CloseButtonPushed, true);
            app.CloseButton.Layout.Row = 1;
            app.CloseButton.Layout.Column = 6;
            app.CloseButton.Text = 'Close';

            % Create ApplyRegistrationtoFolderButton
            app.ApplyRegistrationtoFolderButton = uibutton(app.GridLayoutBottom, 'push');
            app.ApplyRegistrationtoFolderButton.ButtonPushedFcn = createCallbackFcn(app, @ApplyRegistrationtoFolderButtonPushed, true);
            app.ApplyRegistrationtoFolderButton.Layout.Row = 1;
            app.ApplyRegistrationtoFolderButton.Layout.Column = 5;
            app.ApplyRegistrationtoFolderButton.Text = 'Apply Registration to Folder';

            % Create StatusLabel
            app.StatusLabel = uilabel(app.GridLayoutBottom);
            app.StatusLabel.FontAngle = 'italic';
            app.StatusLabel.Layout.Row = 1;
            app.StatusLabel.Layout.Column = [1 4];
            app.StatusLabel.Text = 'Status';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = ImageAlignmentTool_exported(varargin)

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