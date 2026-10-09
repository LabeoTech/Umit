function ref = makeImageAlignmentReference()
%MAKEIMAGEALIGNMENTREFERENCE Record applyImageAlignmentToFolder's outputs (run once).
%
%   ref = makeImageAlignmentReference()
%
%   .dat header Phase 4e-2 characterization. Run against the committed
%   (headerless) applyImageAlignmentToFolder (commit 296abb7) BEFORE it is
%   changed. It runs imageAlignmentCase on a headerless folder and records
%   each .dat output's size and value hash, the image .umt entry's size and
%   hash, DataParams.view.imageSizeYX, and the AcqInfos.mat Height/Width
%   before the run, in fixtures/imageAlignmentReference.mat. The file is
%   then frozen: never regenerate it after the function changes.

c = imageAlignmentCase();
folder = fullfile(tempdir, ['imageAlignRef_' char(java.util.UUID.randomUUID)]);
cleanupFcn = c.build(folder, 'legacySidecar');   % recorded with 'headerless' (AcqInfos-bound) inputs before Phase 5a
try
    acqBefore = load(fullfile(folder, 'AcqInfos.mat'), 'AcqInfoStream');
    evalc('c.run(folder)');

    ref = struct();
    ref.datOutputs = c.datOutputs;
    ref.datSizes = cell(1, numel(c.datOutputs));
    ref.datSha256 = cell(1, numel(c.datOutputs));
    for j = 1:numel(c.datOutputs)
        [~, values] = evalc('loadData(fullfile(folder, c.datOutputs{j}))');
        ref.datSizes{j} = size(values);
        ref.datSha256{j} = c.hash(values);
    end
    U = load(fullfile(folder, c.umtOutput), '-mat');
    ref.umtSize = size(U.umt.data.map.value);
    ref.umtSha256 = c.hash(single(U.umt.data.map.value));
    D = load(fullfile(folder, 'DataParams.mat'), 'DataParams');
    ref.viewImageSizeYX = D.DataParams.view.imageSizeYX;
    ref.acqHeightWidthBefore = [acqBefore.AcqInfoStream.Height, acqBefore.AcqInfoStream.Width];
catch ME
    cleanupFcn();
    rmdir(folder, 's');
    rethrow(ME);
end
cleanupFcn();
rmdir(folder, 's');

fixtureFolder = fullfile(fileparts(mfilename('fullpath')), 'fixtures');
if ~isfolder(fixtureFolder)
    mkdir(fixtureFolder);
end
save(fullfile(fixtureFolder, 'imageAlignmentReference.mat'), 'ref');
end
