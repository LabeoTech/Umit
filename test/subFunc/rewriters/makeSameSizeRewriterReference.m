function ref = makeSameSizeRewriterReference()
%MAKESAMESIZEREWRITERREFERENCE Record the same-size rewriters' outputs (run once).
%
%   ref = makeSameSizeRewriterReference()
%
%   .dat header Phase 4e-1 characterization. Run against the committed
%   (headerless) applyTform2Cams, correctMotionArtifact, and
%   applyRegistrationTformOnFolder (commit 7df18d5) BEFORE they are
%   changed. For each case of sameSizeRewriterCases, the function runs on
%   a headerless folder and records every output's size and value hash,
%   plus the extra results (correctMotionArtifact shifts), in
%   fixtures/sameSizeRewriterReference.mat. The file is then frozen: never
%   regenerate it after the functions change.

[cases, utils] = sameSizeRewriterCases();
ref = struct('name', {}, 'outputs', {}, 'dimSizes', {}, 'sha256', {}, 'extra', {});
for k = 1:numel(cases)
    c = cases(k);
    folder = fullfile(tempdir, ['sameSizeRef_' char(java.util.UUID.randomUUID)]);
    cleanupFcn = c.build(folder, false);
    try
        [~, extra] = evalc('c.run(folder)');
        sizes = cell(1, numel(c.outputs));
        hashes = cell(1, numel(c.outputs));
        for j = 1:numel(c.outputs)
            [~, values] = evalc('loadData(fullfile(folder, c.outputs{j}))');
            sizes{j} = size(values);
            hashes{j} = utils.hash(values);
        end
    catch ME
        cleanupFcn();
        rmdir(folder, 's');
        rethrow(ME);
    end
    cleanupFcn();
    rmdir(folder, 's');
    ref(end+1) = struct('name', c.name, 'outputs', {c.outputs}, 'dimSizes', {sizes}, ...
        'sha256', {hashes}, 'extra', extra); %#ok<AGROW>
end

fixtureFolder = fullfile(fileparts(mfilename('fullpath')), 'fixtures');
if ~isfolder(fixtureFolder)
    mkdir(fixtureFolder);
end
save(fullfile(fixtureFolder, 'sameSizeRewriterReference.mat'), 'ref');
end
