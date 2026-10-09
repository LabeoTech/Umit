function c = imageAlignmentCase()
%IMAGEALIGNMENTCASE Characterization case for applyImageAlignmentToFolder.
%
%   c = imageAlignmentCase()
%
%   .dat header Phase 4e-2. Shared by makeImageAlignmentReference (which
%   recorded the committed function's outputs) and
%   TestImageAlignmentCharacterization. Fields:
%       build(folder, dataClass) - Synthetic folder (setupSyntheticRegistrationFolder
%                        plus one image-kind map.umt). dataClass is
%                        'legacySidecar' (headerless with legacy sidecars;
%                        the reference was recorded from the same values in
%                        headerless files described by AcqInfos.mat before
%                        .dat header Phase 5a), 'single' (headered single
%                        copy), or an integer class such as 'uint16'
%                        (headered, values scaled and rounded). Returns a
%                        cleanup function handle.
%       run(folder)    - Runs applyImageAlignmentToFolder with the case's
%                        similarity transform and reference image.
%       tform, referenceImage - The transform and reference used by run.
%       datOutputs     - The .dat files transformed.
%       umtOutput      - The image .umt transformed.
%       hash(values)   - SHA-256 of the values' bytes.
%       integerScale   - Scale applied before rounding for integer classes.

s = 0.85;
theta = deg2rad(8);
c.tform = affine2d([s*cos(theta), s*sin(theta), 0; ...
    -s*sin(theta), s*cos(theta), 0; ...
    3, -2, 1]);
[yy, xx] = ndgrid(1:50, 1:56);
c.referenceImage = single(yy + 2 * xx);
c.datOutputs = {'green.dat', 'red.dat'};
c.umtOutput = 'map.umt';
c.integerScale = 1000;
c.build = @iBuild;
c.run = @(folder) applyImageAlignmentToFolder(folder, c.tform, ...
    'referenceImage', c.referenceImage);
c.hash = @iHash;
end

function cleanupFcn = iBuild(folder, dataClass)
fixture = setupSyntheticRegistrationFolder(folder, 'DatFormat', 'legacySidecar');
cleanupFcn = @() iRemoveFolders({fixture.ProjectRoot, fixture.RigRoot});

map = reshape(single(1:fixture.Ny * fixture.Nx), fixture.Ny, fixture.Nx);
% ImageAlignmentTool reads a UMT stored as one variable in the MAT file
% (saveData writes the struct fields at top level, which it skips; see
% DFR-20260929-005), so the map is saved in the format it transforms.
umt = genUMTStruct(map, 'kind', 'image', 'entryName', 'map', 'dimNames', {'Y', 'X'}); %#ok<NASGU>
save(fullfile(folder, 'map.umt'), 'umt', '-mat');

if strcmp(dataClass, 'legacySidecar')
    return
end
listing = dir(fullfile(folder, '*.dat'));
for k = 1:numel(listing)
    f = fullfile(folder, listing(k).name);
    [~, info] = evalc('loadMetaData(f)');
    [~, values] = evalc('loadData(f)');
    if ~strcmp(dataClass, 'single')
        values = cast(round(double(values) * 1000), dataClass);
    end
    [~, base] = fileparts(f);
    hdr = datHeaderFromInfo(info, base, 'dataClass', dataClass);
    hdr.writeComplete = true;
    fid = fopen(f, 'w', 'ieee-le');
    fwrite(fid, encodeDatHeader(hdr), 'uint8');
    fwrite(fid, values, dataClass);
    fclose(fid);
    delete(fullfile(folder, [base '.mat']));
end
end

function iRemoveFolders(folders)
for k = 1:numel(folders)
    if isfolder(folders{k})
        rmdir(folders{k}, 's');
    end
end
end

function h = iHash(values)
md = java.security.MessageDigest.getInstance('SHA-256');
md.update(typecast(values(:), 'uint8'));
h = lower(reshape(dec2hex(typecast(md.digest(), 'uint8'))', 1, []));
end
