function [cases, utils] = mergeRecordingsCases()
%MERGERECORDINGSCASES Characterization cases for mergeRecordings.
%
%   [cases, utils] = mergeRecordingsCases()
%
%   .dat header Phase 4e-3. Shared by makeMergeRecordingsReference (which
%   recorded the committed function's outputs) and
%   TestMergeRecordingsCharacterization. utils.build(root, bHeadered)
%   creates three recording folders rec1..rec3 under ROOT, each with a
%   single green.dat (8 x 6 x T, T = 12, 9, 15, 10 Hz), its LabeoTech
%   sidecar green.mat, events.mat, and AcqInfos.mat. With bHeadered, each
%   green.dat is given a header (values, rate, and exposure from
%   loadMetaData) and its sidecar is deleted. Each case has fields:
%       name  - Case name.
%       prep  - @(root) adjusts the folders before the run (or does nothing).
%       run   - @(root, outFile) calls mergeRecordings.
%   utils.outFile(root) is the merged file path; utils.hash(values) and
%   utils.readEvents(folder) are the recording helpers.

cases = struct('name', {}, 'prep', {}, 'run', {});
cases(end+1) = struct('name', 'mergedEvents', 'prep', @(r) [], ...
    'run', @(r, out) mergeRecordings(out, iFolders(r), 'green'));
cases(end+1) = struct('name', 'ignoreEvents', 'prep', @(r) [], ...
    'run', @(r, out) mergeRecordings(out, iFolders(r), 'green', [], {}, true));
cases(end+1) = struct('name', 'trialNamesPermuted', 'prep', @(r) [], ...
    'run', @(r, out) mergeRecordings(out, iFolders(r), 'green', [3 1 2], {'A', 'B', 'C'}));
cases(end+1) = struct('name', 'missingFileSkipped', ...
    'prep', @(r) delete(fullfile(r, 'rec2', 'green.dat')), ...
    'run', @(r, out) mergeRecordings(out, iFolders(r), 'green'));

utils = struct('build', @iBuild, 'outFile', @(r) fullfile(r, 'out', 'merged.dat'), ...
    'folders', @iFolders, 'hash', @iHash, 'readEvents', @iReadEvents);
end

function folders = iFolders(root)
folders = {fullfile(root, 'rec1'), fullfile(root, 'rec2'), fullfile(root, 'rec3')};
end

function iBuild(root, bHeadered)
lengths = [12 9 15];
ny = 8; nx = 6; rate = 10;
names = {{'Air', 'Odor'}, {'Odor'}, {'Tone', 'Air'}};
mkdir(fullfile(root, 'out'));
for k = 1:3
    folder = fullfile(root, sprintf('rec%d', k));
    mkdir(folder);
    nt = lengths(k);
    data = single(reshape(1:ny * nx * nt, ny, nx, nt)) + single(1000 * k);
    fid = fopen(fullfile(folder, 'green.dat'), 'w');
    fwrite(fid, data, 'single');
    fclose(fid);

    meta = struct('Freq', rate, 'datName', 'data', 'datLength', nt, 'FirstDim', 'Y', ...
        'dim_names', {{'Y', 'X', 'T'}}, 'Datatype', 'single', 'datSize', [ny nx]);
    save(fullfile(folder, 'green.mat'), '-struct', 'meta');

    nEv = numel(names{k});
    eventID = uint16(repelem(1:nEv, 2).');
    state = logical(repmat([1; 0], nEv, 1));
    timestamps = single((0:2 * nEv - 1).' * 0.2 + 0.1);
    saveEventsFile(folder, eventID, timestamps, state, names{k});

    AcqInfoStream = struct('Height', ny, 'Width', nx, 'Length', nt, 'FrameRateHz', rate, ...
        'ExposureMsec', 5, 'RecordingIndex', k);
    save(fullfile(folder, 'AcqInfos.mat'), 'AcqInfoStream');

    if bHeadered
        f = fullfile(folder, 'green.dat');
        [~, info] = evalc('loadMetaData(f)');
        [~, values] = evalc('loadData(f)');
        hdr = datHeaderFromInfo(info, 'green');
        hdr.writeComplete = true;
        fid = fopen(f, 'w', 'ieee-le');
        fwrite(fid, encodeDatHeader(hdr), 'uint8');
        fwrite(fid, values, 'single');
        fclose(fid);
        delete(fullfile(folder, 'green.mat'));
    end
end
end

function ev = iReadEvents(folder)
S = load(fullfile(folder, 'events.mat'));
ev = struct('eventID', S.eventID, 'state', S.state, 'timestamps', S.timestamps, ...
    'eventNameList', {S.eventNameList});
end

function h = iHash(values)
md = java.security.MessageDigest.getInstance('SHA-256');
md.update(typecast(values(:), 'uint8'));
h = lower(reshape(dec2hex(typecast(md.digest(), 'uint8'))', 1, []));
end
