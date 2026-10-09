function summary = trimLargeFixtures(backupRoot, action)
%TRIMLARGEFIXTURES Back up and shrink the raw payloads of the large test fixtures.
%
%   summary = trimLargeFixtures(backupRoot)
%   summary = trimLargeFixtures(backupRoot, 'restore')
%
%   Companion to convertLargeFixturesToHeadered. The processed .dat files
%   were already binned; this trims what is left: raw acquisition files
%   that no test reads in full, and one oversized events.mat. Nothing is
%   deleted: every touched file is first moved to BACKUPROOT (same relative
%   path under the repository). ACTION 'restore' moves every backed-up file
%   back, overwriting its trimmed replacement.
%
%   Trims (tests copy only the .dat, AcqInfos.mat and events.mat of the
%   speckle and retinotopy folders):
%     TestingData_speckle/img_0000[0-2].bin   keep img_00000.bin, first 3 frames
%     TestingData_speckle/ai_00000.bin        keep the first 1 s block
%     TestingData_retinotopy/ai_0000[0-7].bin keep ai_00000.bin, first 1 s block
%     TestingData_retinotopy/events.mat       AnalogIN decimated 10x, sr / 10
%   TestingData_with_events is left whole. TestImagesClassification's
%   testPixelValuesMatchHeaderlessReference compares a SHA-256 of every
%   channel built from the full img_00000.bin, and TestGetEvents builds
%   events from the full 54 s ai_00000.bin trace.
%
%   File layouts (as read by IOIAnalysis/ImagesClassification and
%   EventsManager.setAnalogIN):
%     img_*.bin  5 int32 header, then frames of [3 uint64; nx*ny uint16]
%     ai_*.bin   5 int32 header, then 1e4-sample blocks of nChan doubles
%
%   Each trim is verified: the kept bytes equal the original's first bytes
%   and the size satisfies the reader's block/frame arithmetic.

if nargin < 2
    action = 'trim';
end
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));

if strcmpi(action, 'restore')
    summary = iRestore(backupRoot, repoRoot);
    return
end
assert(strcmpi(action, 'trim'), 'Unknown action "%s".', action);

headerBytes = 20;
aiBlockSamples = 1e4;
decimation = 10;

plan = struct('file', {}, 'mode', {}, 'keep', {});
plan(end+1) = struct('file', 'test/Analysis/TestingData_speckle/img_00000.bin', 'mode', 'img', 'keep', 3);
plan(end+1) = struct('file', 'test/Analysis/TestingData_speckle/img_00001.bin', 'mode', 'drop', 'keep', 0);
plan(end+1) = struct('file', 'test/Analysis/TestingData_speckle/img_00002.bin', 'mode', 'drop', 'keep', 0);
plan(end+1) = struct('file', 'test/Analysis/TestingData_speckle/ai_00000.bin', 'mode', 'ai', 'keep', 1);
plan(end+1) = struct('file', 'test/Analysis/TestingData_retinotopy/ai_00000.bin', 'mode', 'ai', 'keep', 1);
for k = 1:7
    plan(end+1) = struct('file', sprintf('test/Analysis/TestingData_retinotopy/ai_%05d.bin', k), ...
        'mode', 'drop', 'keep', 0); %#ok<AGROW>
end
plan(end+1) = struct('file', 'test/Analysis/TestingData_retinotopy/events.mat', 'mode', 'events', 'keep', decimation);

% Fail before touching anything.
for k = 1:numel(plan)
    assert(isfile(fullfile(repoRoot, plan(k).file)), 'Missing %s.', plan(k).file);
    assert(~isfile(fullfile(backupRoot, plan(k).file)), 'Backup %s already exists.', plan(k).file);
end

summary = struct('file', {}, 'mode', {}, 'before', {}, 'after', {});
for k = 1:numel(plan)
    item = plan(k);
    src = fullfile(repoRoot, item.file);
    backupFile = fullfile(backupRoot, item.file);
    if ~isfolder(fileparts(backupFile))
        mkdir(fileparts(backupFile));
    end
    before = dir(src).bytes;
    fprintf('%-8s %s (%.1f MB)\n', item.mode, item.file, before / 1e6);

    [ok, msg] = movefile(src, backupFile);
    assert(ok, 'Could not move %s: %s', src, msg);

    switch item.mode
        case 'drop'
            % Stays in the backup only.
        case 'img'
            frameBytes = iImageFrameBytes(backupFile);
            keepBytes = headerBytes + item.keep * frameBytes;
            assert(mod(before - headerBytes, frameBytes) == 0 && before >= keepBytes, ...
                'Unexpected frame arithmetic in %s.', item.file);
            iCopyPrefix(backupFile, src, keepBytes);
            iVerifyPrefix(backupFile, src, keepBytes);
        case 'ai'
            nChan = iInfoValue(fileparts(src), 'AINChannels');
            keepBytes = headerBytes + item.keep * aiBlockSamples * nChan * 8;
            assert(before >= keepBytes, 'File %s is shorter than the kept blocks.', item.file);
            iCopyPrefix(backupFile, src, keepBytes);
            iVerifyPrefix(backupFile, src, keepBytes);
        case 'events'
            iDecimateEvents(backupFile, src, item.keep);
    end

    if isfile(src)
        after = dir(src).bytes;
    else
        after = 0;
    end
    summary(end+1) = struct('file', item.file, 'mode', item.mode, 'before', before, 'after', after); %#ok<AGROW>
end

fprintf('\nTrimmed %d files: %.1f MB -> %.1f MB\n', numel(summary), ...
    sum([summary.before]) / 1e6, sum([summary.after]) / 1e6);
end

function frameBytes = iImageFrameBytes(file)
fid = fopen(file, 'r');
cleaner = onCleanup(@() fclose(fid));
header = fread(fid, 5, 'int32');
frameBytes = header(2) * header(3) * 2 + 3 * 8;
end

function value = iInfoValue(folder, key)
txt = fileread(fullfile(folder, 'info.txt'));
tok = regexp(txt, ['(?m)^' key ':\s*(\d+)'], 'tokens', 'once');
assert(~isempty(tok), '%s not found in info.txt of %s.', key, folder);
value = str2double(tok{1});
end

function iCopyPrefix(src, dst, nBytes)
fin = fopen(src, 'r');
cleanIn = onCleanup(@() fclose(fin));
bytes = fread(fin, nBytes, '*uint8');
assert(numel(bytes) == nBytes, 'Short read from %s.', src);
fout = fopen(dst, 'w');
cleanOut = onCleanup(@() fclose(fout));
assert(fwrite(fout, bytes, 'uint8') == nBytes, 'Short write to %s.', dst);
end

function iVerifyPrefix(original, trimmed, nBytes)
assert(dir(trimmed).bytes == nBytes, 'Trimmed %s has the wrong size.', trimmed);
fa = fopen(original, 'r');
cleanA = onCleanup(@() fclose(fa));
fb = fopen(trimmed, 'r');
cleanB = onCleanup(@() fclose(fb));
assert(isequal(fread(fa, nBytes, '*uint8'), fread(fb, nBytes, '*uint8')), ...
    'Trimmed %s differs from its original.', trimmed);
end

function iDecimateEvents(original, dst, factor)
fid = fopen(original, 'r');
magic = fread(fid, 20, '*char')';
fclose(fid);
versionFlag = '-v7';
if contains(magic, '7.3')
    versionFlag = '-v7.3';
end

S = load(original);
nOrig = size(S.AnalogIN, 1);
srOrig = double(S.sr);
S.AnalogIN = S.AnalogIN(1:factor:end, :);
S.sr = single(srOrig / factor);
save(dst, '-struct', 'S', versionFlag);

T = load(dst);
assert(size(T.AnalogIN, 1) == ceil(nOrig / factor), 'Decimated AnalogIN has the wrong length.');
assert(abs(size(T.AnalogIN, 1) / double(T.sr) - nOrig / srOrig) < 1 / double(T.sr), ...
    'Decimated AnalogIN changed the recording duration.');
others = setdiff(fieldnames(S), {'AnalogIN', 'sr'});
Orig = load(original, others{:});
for k = 1:numel(others)
    assert(isequaln(T.(others{k}), Orig.(others{k})), 'Field %s changed.', others{k});
end
end

function summary = iRestore(backupRoot, repoRoot)
listing = dir(fullfile(backupRoot, '**', '*'));
listing = listing(~[listing.isdir]);
summary = struct('file', {}, 'mode', {}, 'before', {}, 'after', {});
for k = 1:numel(listing)
    backupFile = fullfile(listing(k).folder, listing(k).name);
    rel = backupFile(numel(backupRoot) + 2:end);
    dst = fullfile(repoRoot, rel);
    if isfile(dst)
        delete(dst);
    end
    [ok, msg] = movefile(backupFile, dst);
    assert(ok, 'Could not restore %s: %s', rel, msg);
    summary(end+1) = struct('file', rel, 'mode', 'restore', 'before', 0, 'after', listing(k).bytes); %#ok<AGROW>
end
fprintf('Restored %d files.\n', numel(summary));
end
