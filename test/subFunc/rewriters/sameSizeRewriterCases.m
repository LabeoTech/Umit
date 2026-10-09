function [cases, utils] = sameSizeRewriterCases()
%SAMESIZEREWRITERCASES Characterization cases for the same-size rewriters.
%
%   [cases, utils] = sameSizeRewriterCases()
%
%   .dat header Phase 4e-1. Shared by makeSameSizeRewriterReference (which
%   recorded the committed functions' outputs) and
%   TestSameSizeRewritersCharacterization. Each case has fields:
%       name    - Case name.
%       build   - @(folder, bHeadered) creates the input folder and returns
%                 a cleanup function handle. The .dat inputs are headerless
%                 with a legacy sidecar .mat (writeTestDat 'legacySidecar';
%                 until .dat header Phase 5a they were headerless files
%                 described only by AcqInfos.mat, with the same values).
%                 With bHeadered, every .dat is given a header right after
%                 it is written (same values, rate and exposure from
%                 loadMetaData of the sidecar file) and its sidecar removed.
%       run     - @(folder) runs the function; returns extra results (the
%                 correctMotionArtifact shifts) or [].
%       outputs - Output .dat file names whose values are recorded.
%       rewritten - The outputs the function rewrites (and so writes with a
%                 header); the others (e.g. Camera 1 for applyTform2Cams)
%                 are recorded to show they are left unchanged.
%   utils.headerize(folder) and utils.hash(values) are the helpers.

cases = struct('name', {}, 'build', {}, 'run', {}, 'outputs', {}, 'rewritten', {});

cases(end+1) = struct('name', 'applyTform2Cams_standard', ...
    'build', @(f, b) iBuildDualCamera(f, b), ...
    'run', @(f) iRunTform2Cams(f, false), ...
    'outputs', {{'green.dat', 'red.dat'}}, 'rewritten', {{'green.dat'}});
cases(end+1) = struct('name', 'applyTform2Cams_RAMsafe', ...
    'build', @(f, b) iBuildDualCamera(f, b), ...
    'run', @(f) iRunTform2Cams(f, true), ...
    'outputs', {{'green.dat', 'red.dat'}}, 'rewritten', {{'green.dat'}});
cases(end+1) = struct('name', 'applyRegistrationTformOnFolder', ...
    'build', @(f, b) iBuildRegistrationFolder(f, b), ...
    'run', @(f) iRunRegistration(f), ...
    'outputs', {{'green.dat', 'red.dat'}}, 'rewritten', {{'green.dat', 'red.dat'}});

utils = struct('headerize', @iHeaderize, 'hash', @iHash);
end

% =========================================================================
% Dual-camera folder (applyTform2Cams)
% =========================================================================
function cleanupFcn = iBuildDualCamera(folder, bHeadered)
mkdir(folder);
frameSize = [8, 8];
nt = 5;
AcqInfoStream = struct();
AcqInfoStream.Width = frameSize(2);
AcqInfoStream.Height = frameSize(1);
AcqInfoStream.Length = nt;
AcqInfoStream.FrameRateHz = 10;
AcqInfoStream.ExposureMsec = 3;
AcqInfoStream.Datatype = 'single';
AcqInfoStream.MultiCam = 1;
AcqInfoStream.Binning = 4;
AcqInfoStream.BinningSpatial = 1;
AcqInfoStream.Illumination1 = struct('Color', 'red', 'CamIdx', 1, 'FrameIdx', 1);
AcqInfoStream.Illumination2 = struct('Color', 'green', 'CamIdx', 2, 'FrameIdx', 1);
save(fullfile(folder, 'AcqInfos.mat'), 'AcqInfoStream');

rs = RandStream('mt19937ar', 'Seed', 7);
iWriteRaw(fullfile(folder, 'red.dat'), rand(rs, [frameSize nt], 'single'), 10, 3);
[yy, xx] = ndgrid(1:frameSize(1), 1:frameSize(2));
pattern = single(exp(-((yy - 4).^2 + (xx - 5).^2) / 6));
green = pattern .* reshape(single(1:nt), 1, 1, []) + 0.01 * rand(rs, [frameSize nt], 'single');
iWriteRaw(fullfile(folder, 'green.dat'), green, 10, 3);
if bHeadered
    iHeaderize(folder);
end
cleanupFcn = @() [];
end

function out = iRunTform2Cams(folder, bRAMSafe)
tform = affine2d([1 0 0; 0 1 0; 4 -2 1]);
tformInfo = struct('Binning', 2, 'BinningSpatial', 1, 'Rotation', 0, ...
    'X_Offset', 0, 'Y_Offset', 0);
[status, warnmsg] = applyTform2Cams(folder, tform, tformInfo, bRAMSafe);
assert(status, 'applyTform2Cams failed: %s', warnmsg);
out = [];
end

% =========================================================================
% Shifted-stack folder (correctMotionArtifact), as in TestCorrectMotionArtifact
% =========================================================================
function cleanupFcn = iBuildShiftedFolder(folder, bHeadered)
mkdir(folder);
ny = 40; nx = 40; nt = 4;
[xg, yg] = meshgrid(1:nx, 1:ny);
refPattern = single(exp(-((xg-15).^2 + (yg-18).^2)/40) + ...
    0.5*exp(-((xg-28).^2 + (yg-25).^2)/60));
shifts = [0 0; 2 -3; -1 4; 3 1];
iWriteRaw(fullfile(folder, 'green.dat'), iShiftedStack(refPattern, shifts), 10, 20);
iWriteRaw(fullfile(folder, 'red.dat'), iShiftedStack(refPattern * 0.6, shifts), 10, 20);
AcqInfoStream = struct('Width', nx, 'Height', ny, 'Length', nt, ...
    'FrameRateHz', 10, 'ExposureMsec', 20);
save(fullfile(folder, 'AcqInfos.mat'), 'AcqInfoStream');
if bHeadered
    iHeaderize(folder);
end
cleanupFcn = @() [];
end

function stack = iShiftedStack(pattern, expectedShifts)
nt = size(expectedShifts, 1);
stack = zeros([size(pattern), nt], 'single');
for t = 1:nt
    stack(:,:,t) = imtranslate(pattern, [-expectedShifts(t, 2), -expectedShifts(t, 1)], ...
        'cubic', 'FillValues', 0);
end
end

% =========================================================================
% Registration folder (applyRegistrationTformOnFolder)
% =========================================================================
function cleanupFcn = iBuildRegistrationFolder(folder, bHeadered)
fixture = setupSyntheticRegistrationFolder(folder, 'DatFormat', 'legacySidecar');
if bHeadered
    iHeaderize(folder);
end
evalc('createRegistrationTform(folder, ''ShowFigure'', false)');
cleanupFcn = @() iRemoveFolders({fixture.ProjectRoot, fixture.RigRoot});
end

function out = iRunRegistration(folder)
applyRegistrationTformOnFolder(folder, 'RequireUserConfirmation', false, 'OpenQCFigure', false);
out = [];
end

function iRemoveFolders(folders)
for k = 1:numel(folders)
    if isfolder(folders{k})
        rmdir(folders{k}, 's');
    end
end
end

% =========================================================================
% Helpers
% =========================================================================
function iWriteRaw(filePath, data, frameRateHz, exposureMsec)
% Headerless single .dat with a legacy sidecar (the values the committed
% functions were characterized with).
writeTestDat(filePath, single(data), frameRateHz, exposureMsec, 'Format', 'legacySidecar');
end

function iHeaderize(folder)
% Put a v1 header in front of every .dat, with the values, rate, and
% exposure loadMetaData gives the headerless file.
listing = dir(fullfile(folder, '*.dat'));
for k = 1:numel(listing)
    f = fullfile(folder, listing(k).name);
    if isDatWithHeader(f)
        continue
    end
    [~, info] = evalc('loadMetaData(f)');
    [~, values] = evalc('loadData(f)');
    [~, base] = fileparts(f);
    hdr = datHeaderFromInfo(info, base);
    hdr.writeComplete = true;
    fid = fopen(f, 'w', 'ieee-le');
    fwrite(fid, encodeDatHeader(hdr), 'uint8');
    fwrite(fid, values, info.dataClass);
    fclose(fid);
    sidecar = fullfile(folder, [base '.mat']);
    if isfile(sidecar)
        delete(sidecar);
    end
end
end

function h = iHash(values)
% SHA-256 of the values' bytes in memory (column-major) order.
md = java.security.MessageDigest.getInstance('SHA-256');
md.update(typecast(values(:), 'uint8'));
h = lower(reshape(dec2hex(typecast(md.digest(), 'uint8'))', 1, []));
end
