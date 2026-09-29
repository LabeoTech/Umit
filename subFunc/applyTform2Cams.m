function [status, warnmsg] = applyTform2Cams(DataFolder, tform, tformInfo, bRAMSafeMode)
%APPLYTFORM2CAMS Apply geometric transformation to Camera 2 .dat files.
%
%   [status, warnmsg] = applyTform2Cams(DataFolder, tform, tformInfo)
%   [status, warnmsg] = applyTform2Cams(DataFolder, tform, tformInfo, bRAMSafeMode)
%
%   This helper applies the geometric transformation (tform) to the
%   processed .dat files belonging to Camera 2 from a multi-camera OiS200 LightTrack
%   Imaging System (LabeoTech). The transformed data are written back to the
%   original .dat files.
%
%   Inputs:
%       DataFolder    - Folder containing AcqInfos.mat and channel .dat files.
%       tform         - affine2d geometric transformation.
%       tformInfo     - Struct containing extra parameters from the saved
%                       tform file.
%       bRAMSafeMode  - Logical scalar. If true, process one frame at a
%                       time using a temporary file. Default: false
%
%   Outputs:
%       status        - True when the operation succeeds.
%       warnmsg       - Warning/error message when status is false.
%
%   Notes:
%       - This function is intended for already-classified channel .dat
%         files from dual-camera acquisitions.
%       - Frame size is resolved from AcqInfos.mat directly.
%       - In RAM-safe mode, each file is rewritten through a temporary file
%         and replaced only after successful completion.

if nargin < 4 || isempty(bRAMSafeMode)
    bRAMSafeMode = false;
end

warnmsg = '';
status = false;

if ~(ischar(DataFolder) || (isstring(DataFolder) && isscalar(DataFolder)))
    warnmsg = 'DataFolder must be a character vector or string scalar.';
    return
end
DataFolder = char(string(DataFolder));

if ~isfolder(DataFolder)
    warnmsg = sprintf('Data folder not found: "%s".', DataFolder);
    return
end

if ~exist(fullfile(DataFolder, 'AcqInfos.mat'), 'file')
    warnmsg = ['AcqInfos.mat file not found in "' DataFolder '". ' ...
        'Run ImagesClassification and try again.'];
    return
end

S = load(fullfile(DataFolder, 'AcqInfos.mat'), 'AcqInfoStream');
if ~isfield(S, 'AcqInfoStream') || ~isstruct(S.AcqInfoStream) || ~isscalar(S.AcqInfoStream)
    warnmsg = 'AcqInfos.mat does not contain a valid scalar AcqInfoStream struct.';
    return
end
AcqInfo = S.AcqInfoStream;
clear S

assertBinningMetadata(AcqInfo, 'AcqInfos.mat', ...
    'Re-import the original raw data with the current data importer.');

% Get list of channels for each camera:
if AcqInfo.MultiCam
    NbIllum = sum(cellfun(@(X) contains(X, 'Illumination'), fieldnames(AcqInfo)));
    Cam1List = {};
    Cam2List = {};
    for ind = 1:NbIllum
        idx = AcqInfo.("Illumination" + int2str(ind)).CamIdx;
        chan = lower(AcqInfo.("Illumination" + int2str(ind)).Color);
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
       
        if contains(chan, 'amber')
            chan = 'yellow';
        end
        if idx == 1
            Cam1List{end + 1} = [chan '.dat']; %#ok<AGROW>
        else
            Cam2List{end + 1} = [chan '.dat']; %#ok<AGROW>
        end
    end
    clear NbIllum ind idx chan
else
    disp('Only one camera was used. No need to coregister images');
    status = true;
    return;
end

if isempty(Cam2List)
    warnmsg = 'No Camera 2 files were found from AcqInfoStream illumination metadata.';
    return
end

assertBinningMetadata(tformInfo, 'camera-coregistration transform', ...
    'Regenerate the camera-coregistration transform with the current calibration workflow.');

% Check for existence of rotation and X/Y offset fields
acqFieldNames = fieldnames(AcqInfo);
tformFieldNames = fieldnames(tformInfo);

fNames = {'Rotation', 'X_Offset', 'Y_Offset'};
defaults = [0 0 0];
for ii = 1:length(fNames)
    if ~ismember(fNames{ii}, acqFieldNames)
        AcqInfo.(fNames{ii}) = defaults(ii);
    end
    if ~ismember(fNames{ii}, tformFieldNames)
        tformInfo.(fNames{ii}) = defaults(ii);
    end
end

if ~isfield(AcqInfo, 'Height') || ~isfield(AcqInfo, 'Width')
    warnmsg = 'AcqInfoStream is missing Height and/or Width.';
    return
end

frameSizeYX = [double(AcqInfo.Height), double(AcqInfo.Width)];

% Account for rotation in acquisition software:
rot_diff = 90 * double(AcqInfo.Rotation) - 90 * double(tformInfo.Rotation);

% Update TFORM to account for differences in rotation, binning and ROI offset:
tform = updateTForm(tform, tformInfo, AcqInfo, frameSizeYX, rot_diff);
RA = imref2d(frameSizeYX);

% Apply tform to data from Camera 2:
for ii = 1:length(Cam2List)
    fprintf('----------------------------------\n')
    fprintf('Coregistration of file: "%s"\n', Cam2List{ii})

    datPath = fullfile(DataFolder, Cam2List{ii});
    if ~isfile(datPath)
        warnmsg = sprintf('Camera 2 file not found: "%s".', datPath);
        return
    end

    fprintf('\t- Loading metadata...\n')
    md = loadMetaData(datPath);

    if ~isfield(md, 'dimNames') || ~isfield(md, 'dimSizes') || datAxisSize(md, 'T') == 0
        warnmsg = sprintf('Could not resolve the Y, X, and T sizes of "%s".', datPath);
        return
    end
    if ~strcmp(md.dataClass, 'single')
        error('Umitoolbox:applyTform2Cams:unsupportedDataClass', ...
            'applyTform2Cams supports single-precision .dat files; "%s" stores %s.', ...
            datPath, md.dataClass);
    end

    ny = datAxisSize(md, 'Y');
    nx = datAxisSize(md, 'X');
    nt = datAxisSize(md, 'T');

    if ny ~= frameSizeYX(1) || nx ~= frameSizeYX(2)
        warnmsg = sprintf(['File "%s" has frame size [%d %d], which does not match ' ...
            'AcqInfos.mat frame size [%d %d].'], Cam2List{ii}, ny, nx, frameSizeYX(1), frameSizeYX(2));
        return
    end

    % Write to a temporary headered file and replace the original only
    % after a fully successful write, so the only copy of the data is never
    % truncated. The header keeps the input's class, sizes, rate, exposure,
    % and name.
    tmpPath = [datPath '.tmp'];
    if isfile(tmpPath)
        delete(tmpPath);
    end
    [~, channelName] = fileparts(datPath);
    outHeader = datHeaderFromInfo(md, channelName);

    if ~bRAMSafeMode
        % Standard mode
        fprintf('\t- Loading data...\n')
        dat = loadData(datPath);
        fprintf('\t- Applying geometric transformation...\n')
        dat = imwarp(dat, RA, tform, 'nearest', 'OutputView', RA, 'FillValues', 0);

        fprintf('\t- Writing transformed data to a temporary file...\n')
        try
            slabOut = spatialSlabIO('create', tmpPath, outHeader);
            cleanupOut = onCleanup(@() spatialSlabIO('close', slabOut));
            spatialSlabIO('write', slabOut, 1:nx, dat);
            spatialSlabIO('finalize', slabOut);
            clear cleanupOut
        catch ME
            clear cleanupOut
            iDeleteIfExists(tmpPath);
            rethrow(ME);
        end
    else
        % RAM-safe mode
        fprintf('\t- Applying geometric transformation in RAM-safe mode...\n')
        try
            slabIn = spatialSlabIO('open', datPath, 'Info', md);
            cleanupIn = onCleanup(@() spatialSlabIO('close', slabIn));
            slabOut = spatialSlabIO('create', tmpPath, outHeader);
            cleanupOut = onCleanup(@() spatialSlabIO('close', slabOut));

            for t = 1:nt
                frame = reshape(spatialSlabIO('read', slabIn, 1:nx, t), ny, nx);
                frame = imwarp(frame, RA, tform, 'nearest', 'OutputView', RA, 'FillValues', 0);
                spatialSlabIO('write', slabOut, 1:nx, frame, t);
            end
            spatialSlabIO('finalize', slabOut);
            clear cleanupIn cleanupOut
        catch ME
            clear cleanupIn cleanupOut
            iDeleteIfExists(tmpPath);
            rethrow(ME);
        end
    end

    fprintf('\t- Replacing original .DAT file...\n')
    delete(datPath);
    movefile(tmpPath, datPath, 'f');

    fprintf('Done.\n')
    fprintf('----------------------------------\n')
end

status = true;

% Save copy of updated tform in Data folder:
save(fullfile(DataFolder, 'tformDualCam.mat'), 'tform');
end

% Local function

function iDeleteIfExists(filePath)
if isfile(filePath)
    delete(filePath);
end
end

function newtform = updateTForm(tform, tf_info, acqInfo, frameSizeYX, ang)
%UPDATETFORM Update a geometric transformation to account for processed geometry.
%
% Inputs:
%   tform       - Original geometric transformation.
%   tf_info     - Transformation info struct.
%   acqInfo     - Processed acquisition info struct from AcqInfos.mat.
%   frameSizeYX - Processed frame size [Height Width].
%   ang         - Rotation angle in degrees.
%
% Output:
%   newtform    - Updated affine2d transformation.

% 1. Process Spatial Binning
AcqBinFactor = double(acqInfo.Binning) * double(acqInfo.BinningSpatial);
TformBinFactor = double(tf_info.Binning) * double(tf_info.BinningSpatial);
binFactor = AcqBinFactor / TformBinFactor;
binningMat = [binFactor 0 0; 0 binFactor 0; 0 0 1];

% 2. Process XY Offset
Xoffset = double(acqInfo.X_Offset) - double(tf_info.X_Offset);
Yoffset = double(acqInfo.Y_Offset) - double(tf_info.Y_Offset);
offsetMat = [1 0 0; 0 1 0; Xoffset Yoffset 1];

% 3. Create Rotation and Centering Matrices
frSize = frameSizeYX;
centerImg = [1 0 0; 0 1 0; -frSize(2)/2 -frSize(1)/2 1];
centerImgFlip = [1 0 0; 0 1 0; -frSize(1)/2 -frSize(2)/2 1];
rot = [cosd(ang) sind(ang) 0; -sind(ang) cosd(ang) 0; 0 0 1];

% 4. Create New tform
switch abs(ang)
    case 0
        newMat = binningMat * offsetMat * tform.T * inv(offsetMat) * inv(binningMat);
    case {90, 270}
        newMat = centerImg * rot * inv(centerImgFlip) * binningMat * offsetMat * ...
            tform.T * inv(offsetMat) * inv(binningMat) * centerImgFlip * inv(rot) * inv(centerImg);
    case 180
        newMat = centerImg * rot * inv(centerImg) * binningMat * offsetMat * ...
            tform.T * inv(offsetMat) * inv(binningMat) * centerImg * inv(rot) * inv(centerImg);
    otherwise
        error('Umitoolbox:applyTform2Cams:InvalidRotation', ...
            'Unsupported rotation difference: %g degrees.', ang);
end

newMat(:,3) = [0; 0; 1];
newtform = affine2d(newMat);
end

function assertBinningMetadata(info, sourceName, recoveryInstruction)
%ASSERTBINNINGMETADATA Require unambiguous hardware and software binning.

requiredFields = {'Binning', 'BinningSpatial'};
hasRequiredFields = isstruct(info) && isscalar(info) && ...
    all(isfield(info, requiredFields));

if hasRequiredFields
    for iField = 1:numel(requiredFields)
        value = info.(requiredFields{iField});
        if ~(isnumeric(value) && isscalar(value) && isreal(value) && ...
                isfinite(value) && value > 0)
            hasRequiredFields = false;
            break
        end
    end
end

if ~hasRequiredFields
    error('Umitoolbox:applyTform2Cams:MissingBinningMetadata', ...
        ['%s must contain positive scalar Binning (hardware) and ' ...
         'BinningSpatial (software classification) fields. The metadata ' ...
         'are legacy or incomplete. %s'], sourceName, recoveryInstruction);
end
end
