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

% Check for existence of BinningSpatial/Binning, rotation and X/Y offset fields
acqFieldNames = fieldnames(AcqInfo);
tformFieldNames = fieldnames(tformInfo);

if ~isfield(AcqInfo, 'BinningSpatial')
    if isfield(AcqInfo, 'Binning')
        AcqInfo.BinningSpatial = AcqInfo.Binning;
    else
        AcqInfo.BinningSpatial = 1;
    end
end
if ~isfield(tformInfo, 'BinningSpatial')
    if isfield(tformInfo, 'Binning')
        tformInfo.BinningSpatial = tformInfo.Binning;
    else
        tformInfo.BinningSpatial = 1;
    end
end

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

    if ~isfield(md, 'Height') || ~isfield(md, 'Width') || ~isfield(md, 'Length')
        warnmsg = sprintf('Could not resolve Height/Width/Length for "%s".', datPath);
        return
    end

    ny = double(md.Height);
    nx = double(md.Width);
    nt = double(md.Length);

    if ny ~= frameSizeYX(1) || nx ~= frameSizeYX(2)
        warnmsg = sprintf(['File "%s" has frame size [%d %d], which does not match ' ...
            'AcqInfos.mat frame size [%d %d].'], Cam2List{ii}, ny, nx, frameSizeYX(1), frameSizeYX(2));
        return
    end

    if ~bRAMSafeMode
        % Standard mode
        fprintf('\t- Loading data...\n')
        dat = loadData(datPath);
        fprintf('\t- Applying geometric transformation...\n')
        dat = imwarp(dat, RA, tform, 'nearest', 'OutputView', RA, 'FillValues', 0);

        % Write to a temporary file and replace the original only after a
        % fully successful write. Opening datPath directly with 'w' would
        % truncate it immediately, destroying the only copy of the data if
        % the write subsequently failed for any reason.
        fprintf('\t- Writing transformed data to a temporary file...\n')
        tmpPath = [datPath '.tmp'];
        if isfile(tmpPath)
            delete(tmpPath);
        end

        fid = fopen(tmpPath, 'w');
        if fid == -1
            warnmsg = sprintf('Could not open temporary file "%s" for writing.', tmpPath);
            return
        end
        cleanupFid = onCleanup(@() safeFclose(fid));

        nWritten = fwrite(fid, dat, 'single');
        clear cleanupFid % triggers fclose(fid) now, before movefile

        if nWritten ~= numel(dat)
            if isfile(tmpPath)
                delete(tmpPath);
            end
            warnmsg = sprintf(['Failed to write all data to temporary file for "%s" ' ...
                '(%d/%d elements written). Original file was not modified.'], ...
                datPath, nWritten, numel(dat));
            return
        end

        fprintf('\t- Replacing original .DAT file...\n')
        delete(datPath);
        movefile(tmpPath, datPath, 'f');
    else
        % RAM-safe mode
        fprintf('\t- Applying geometric transformation in RAM-safe mode...\n')
        tmpPath = [datPath '.tmp'];

        fidIn = fopen(datPath, 'r');

        if fidIn == -1
            warnmsg = sprintf('Could not open "%s" for reading.', datPath);
            return
        end

        fidOut = fopen(tmpPath, 'w');
        if fidOut == -1
            fclose(fidIn);
            warnmsg = sprintf('Could not open temporary file "%s" for writing.', tmpPath);
            return
        end

        cleanupIn = onCleanup(@() safeFclose(fidIn));
        cleanupOut = onCleanup(@() safeFclose(fidOut));

        bytesPerFrame = ny * nx * 4; % single precision

        for t = 1:nt
            fseek(fidIn, (t-1) * bytesPerFrame, 'bof');
            frame = fread(fidIn, ny * nx, '*single');
            if numel(frame) ~= ny * nx
                warnmsg = sprintf('Could not read frame %d from "%s".', t, datPath);
                if isfile(tmpPath)
                    delete(tmpPath);
                end
                return
            end
            frame = reshape(frame, ny, nx);
            frame = imwarp(frame, RA, tform, 'nearest', 'OutputView', RA, 'FillValues', 0);
            nWritten = fwrite(fidOut, frame, 'single');
            if nWritten ~= numel(frame)
                warnmsg = sprintf(['Failed to write frame %d to temporary file for "%s" ' ...
                    '(%d/%d elements written). Original file was not modified.'], ...
                    t, datPath, nWritten, numel(frame));
                clear cleanupIn cleanupOut
                if isfile(tmpPath)
                    delete(tmpPath);
                end
                return
            end
        end

        clear cleanupIn cleanupOut
        delete(datPath);
        movefile(tmpPath, datPath, 'f');
    end

    fprintf('Done.\n')
    fprintf('----------------------------------\n')
end

status = true;

% Save copy of updated tform in Data folder:
save(fullfile(DataFolder, 'tformDualCam.mat'), 'tform');
end

% Local function

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
AcqBinFactor = double(acqInfo.BinningSpatial);
AcqBinFactor = AcqBinFactor * double(acqInfo.Width) / frameSizeYX(2);
binFactor = AcqBinFactor / double(tf_info.BinningSpatial);
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
