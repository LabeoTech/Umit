function ref = makeMergeRecordingsReference()
%MAKEMERGERECORDINGSREFERENCE Record mergeRecordings' outputs (run once).
%
%   ref = makeMergeRecordingsReference()
%
%   .dat header Phase 4e-3 characterization. Run against the committed
%   (sidecar-based, headerless) mergeRecordings (commit 325137f) BEFORE it
%   is changed. For each case of mergeRecordingsCases, on legacy sidecar
%   folders, it records the merged size and value hash, events.mat, and the
%   copied AcqInfos.mat in fixtures/mergeRecordingsReference.mat. The file
%   is then frozen: never regenerate it after the function changes.
%
%   The committed function ends by calling genDataHistory for its output
%   sidecar, but genDataHistory was removed from dev (commit 5fd1ae4), so
%   every call failed there after writing the merged .dat, events.mat, and
%   AcqInfos.mat. A stub that returns an empty struct is put on the path
%   while recording; it only affects the sidecar, which Phase 4e-3 drops.

stubFolder = tempname;
mkdir(stubFolder);
fid = fopen(fullfile(stubFolder, 'genDataHistory.m'), 'w');
fprintf(fid, ['function dH = genDataHistory(varargin)\n' ...
    '%%Recording stub (see makeMergeRecordingsReference).\n' ...
    'dH = struct();\n' ...
    'end\n']);
fclose(fid);
addpath(stubFolder);
cleanupStub = onCleanup(@() iRemoveStub(stubFolder));

[cases, utils] = mergeRecordingsCases();
ref = struct('name', {}, 'dimSizes', {}, 'sha256', {}, 'events', {}, 'acqInfos', {});
for k = 1:numel(cases)
    c = cases(k);
    root = fullfile(tempdir, ['mergeRef_' char(java.util.UUID.randomUUID)]);
    mkdir(root);
    try
        utils.build(root, false);
        c.prep(root);
        outFile = utils.outFile(root);
        evalc('c.run(root, outFile)');
        [~, values] = evalc('loadData(outFile)');
        acq = load(fullfile(fileparts(outFile), 'AcqInfos.mat'));
        ref(end+1) = struct('name', c.name, 'dimSizes', size(values), ...
            'sha256', utils.hash(values), 'events', utils.readEvents(fileparts(outFile)), ...
            'acqInfos', acq); %#ok<AGROW>
    catch ME
        rmdir(root, 's');
        rethrow(ME);
    end
    rmdir(root, 's');
end

fixtureFolder = fullfile(fileparts(mfilename('fullpath')), 'fixtures');
if ~isfolder(fixtureFolder)
    mkdir(fixtureFolder);
end
save(fullfile(fixtureFolder, 'mergeRecordingsReference.mat'), 'ref');
clear cleanupStub
end

function iRemoveStub(stubFolder)
rmpath(stubFolder);
rmdir(stubFolder, 's');
end
