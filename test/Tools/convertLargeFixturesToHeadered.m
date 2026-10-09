function summary = convertLargeFixturesToHeadered(backupRoot)
%CONVERTLARGEFIXTURESTOHEADERED Shrink and header the large test fixtures (run once).
%
%   summary = convertLargeFixturesToHeadered(backupRoot)
%
%   .dat header Phase 5a. Must run while loadMetaData still reads
%   AcqInfos-bound .dat files (before Phase 5b). For each fixture below it:
%     1) reads the original (AcqInfos-bound) file in frame chunks;
%     2) bins it by block mean (spatial s x s, temporal t; trailing pixels
%        or frames that do not fill a block are dropped);
%     3) moves the original .dat and AcqInfos.mat to BACKUPROOT (same
%        relative path under the repository), deleting nothing;
%     4) writes the binned data with saveData (header: frame rate / t, the
%        original file's own exposure, channelName = file name);
%     5) spot-checks 300 random output voxels against block means computed
%        directly from the moved original;
%     6) writes the folder's AcqInfos.mat with Height/Width/Length/
%        FrameRateHz (and ImportedChannels Length/FrameRateHz,
%        BinningSpatial, BinningTemp when present) matching the new files.
%
%   Factors (spatial x temporal) and results:
%     TestingData_retinotopy/fluo.dat           4 x 2   256x256x6180 @3.33 Hz -> 64x64x3090 @1.67 Hz
%     TestingData_speckle/speckle.dat           4 x 1   512x512x603  @5 Hz    -> 128x128x603 @5 Hz
%     TestingData_with_events/{green,red,yellow} 2 x 1  112x112x1056 @20 Hz   -> 56x56x1056 @20 Hz
%     GSR/green.dat                              1 x 1  64x64x203 @2.5 Hz (header only)
%
%   Raw acquisition files (ai_*.bin, img_*.bin, info.txt), events.mat, and
%   snapshots are not touched.

repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
fixtures = struct( ...
    'folder', {'test/Analysis/TestingData_retinotopy', 'test/Analysis/TestingData_speckle', ...
               'test/Analysis/TestingData_with_events', 'test/Analysis/GSR'}, ...
    'files', {{'fluo.dat'}, {'speckle.dat'}, {'green.dat', 'red.dat', 'yellow.dat'}, {'green.dat'}}, ...
    'spatial', {4, 4, 2, 1}, ...
    'temporal', {2, 1, 1, 1});

summary = struct('file', {}, 'before', {}, 'after', {}, 'rateHz', {}, 'exposureMsec', {}, 'maxRelError', {});
for iFix = 1:numel(fixtures)
    fx = fixtures(iFix);
    folder = fullfile(repoRoot, fx.folder);
    backupFolder = fullfile(backupRoot, fx.folder);
    if ~isfolder(backupFolder)
        mkdir(backupFolder);
    end

    S = load(fullfile(folder, 'AcqInfos.mat'));
    acq = S.AcqInfoStream;
    converted = struct('size', {}, 'rate', {});

    for iFile = 1:numel(fx.files)
        src = fullfile(folder, fx.files{iFile});
        assert(~isDatWithHeader(src), '%s is already headered.', src);
        info = loadMetaData(src);
        [binned, rateHz] = iBinFile(src, info, fx.spatial, fx.temporal);

        backupFile = fullfile(backupFolder, fx.files{iFile});
        assert(~isfile(backupFile), 'Backup %s already exists.', backupFile);
        [ok, msg] = movefile(src, backupFile);
        assert(ok, 'Could not move %s: %s', src, msg);

        saveData(src, binned, 'FrameRateHz', rateHz, 'Info', struct('exposureMsec', info.exposureMsec), ...
            'DimNames', {'Y', 'X', 'T'});
        written = loadData(src);
        assert(isequal(written, binned), 'Written values differ from the binned array for %s.', src);
        maxRel = iSpotCheck(backupFile, info, written, fx.spatial, fx.temporal);
        assert(maxRel < 1e-5, 'Spot check failed for %s (max relative error %g).', src, maxRel);

        summary(end+1) = struct('file', src, 'before', info.dimSizes, 'after', size(written), ...
            'rateHz', rateHz, 'exposureMsec', info.exposureMsec, 'maxRelError', maxRel); %#ok<AGROW>
        converted(end+1) = struct('size', size(written), 'rate', rateHz); %#ok<AGROW>
        fprintf('%s: %s -> %s, %g Hz, exposure %g ms (spot check max rel. error %.2g)\n', ...
            src, mat2str(info.dimSizes), mat2str(size(written)), rateHz, info.exposureMsec, maxRel);
    end

    % AcqInfos.mat matching the converted files (the original is moved to the backup).
    [ok, msg] = movefile(fullfile(folder, 'AcqInfos.mat'), fullfile(backupFolder, 'AcqInfos.mat'));
    assert(ok, 'Could not move AcqInfos.mat of %s: %s', folder, msg);
    newSize = converted(1).size;
    acq.Height = newSize(1);
    acq.Width = newSize(2);
    acq.Length = newSize(3);
    acq.FrameRateHz = converted(1).rate;
    if isfield(acq, 'ImportedChannels')
        for k = 1:numel(acq.ImportedChannels)
            acq.ImportedChannels(k).Length = floor(double(acq.ImportedChannels(k).Length) / fx.temporal);
            acq.ImportedChannels(k).FrameRateHz = double(acq.ImportedChannels(k).FrameRateHz) / fx.temporal;
        end
    end
    if isfield(acq, 'BinningSpatial')
        acq.BinningSpatial = double(acq.BinningSpatial) * fx.spatial;
    end
    if isfield(acq, 'BinningTemp')
        acq.BinningTemp = double(acq.BinningTemp) * fx.temporal;
    end
    S.AcqInfoStream = acq;
    save(fullfile(folder, 'AcqInfos.mat'), '-struct', 'S');
end
end

function [binned, rateHz] = iBinFile(src, info, s, t)
% Block-mean binning, read in frame chunks through spatialSlabIO.
ny = datAxisSize(info, 'Y'); nx = datAxisSize(info, 'X'); nt = datAxisSize(info, 'T');
nyB = floor(ny / s); nxB = floor(nx / s); ntB = floor(nt / t);
binned = zeros(nyB, nxB, ntB, 'single');
h = spatialSlabIO('open', src, 'Info', info);
cleanupObj = onCleanup(@() spatialSlabIO('close', h));
chunk = max(1, floor(2e8 / (4 * ny * nx * t))) * t;       % ~200 MB per read
for f1 = 1:chunk:ntB * t
    f2 = min(f1 + chunk - 1, ntB * t);
    block = double(spatialSlabIO('read', h, 1:nx, f1:f2));
    block = block(1:nyB * s, 1:nxB * s, :);
    nF = size(block, 3) / t;
    block = reshape(block, s, nyB, s, nxB, t, nF);
    block = squeeze(mean(mean(mean(block, 1), 3), 5));
    binned(:, :, (f1 - 1) / t + (1:nF)) = single(reshape(block, nyB, nxB, nF));
end
rateHz = double(info.frameRateHz) / t;
end

function maxRel = iSpotCheck(original, info, written, s, t)
% Recompute 300 random output voxels directly from the original file.
infoOrig = info;
infoOrig.filePath = original;
h = spatialSlabIO('open', original, 'Info', infoOrig);
cleanupObj = onCleanup(@() spatialSlabIO('close', h));
rs = RandStream('mt19937ar', 'Seed', 5);
maxRel = 0;
for k = 1:300
    iy = randi(rs, size(written, 1)); ix = randi(rs, size(written, 2)); it = randi(rs, size(written, 3));
    cols = (ix - 1) * s + (1:s);
    frames = (it - 1) * t + (1:t);
    v = double(spatialSlabIO('read', h, cols, frames));
    v = v((iy - 1) * s + (1:s), :, :);
    expected = mean(v(:));
    rel = abs(double(written(iy, ix, it)) - expected) / max(abs(expected), eps);
    maxRel = max(maxRel, rel);
end
end
