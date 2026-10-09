function makeLegacySidecarFixture(outFolder, varargin)
%MAKELEGACYSIDECARFIXTURE Generate the legacy sidecar .dat reference file.
%
%   makeLegacySidecarFixture()
%   makeLegacySidecarFixture(outFolder)
%   makeLegacySidecarFixture(outFolder, 'Overwrite', true)
%
%   Documents how the legacy reference file used by TestDatFormatLoading
%   was produced. A legacy .dat file is headerless and described by its own
%   sidecar .mat (the Astrocyte-era format). The folder deliberately has
%   no AcqInfos.mat: legacy files must open from their sidecar alone.
%
%   The files are generated once and must never be regenerated after the
%   .dat header Phase 2 task is accepted; the function refuses to overwrite
%   existing files unless 'Overwrite' is true.
%
%   Default outFolder: the "legacy_sidecar" folder next to this file.
%
%   Files written:
%       green.dat  headerless, Y=6 X=5 T=4, single, values 1:120,
%                  little-endian, column-major
%       green.mat  datSize = [6 5], datLength = 4, Freq = 10,
%                  Datatype = 'single', dim_names = {'Y','X','T'}

p = inputParser;
addOptional(p, 'outFolder', fullfile(fileparts(mfilename('fullpath')), 'legacy_sidecar'), ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'Overwrite', false, @(x) islogical(x) && isscalar(x));
if nargin < 1
    parse(p);
else
    parse(p, outFolder, varargin{:});
end
outFolder = char(p.Results.outFolder);

targets = {fullfile(outFolder, 'green.dat'), fullfile(outFolder, 'green.mat')};
existing = targets(cellfun(@isfile, targets));
if ~isempty(existing) && ~p.Results.Overwrite
    error('Umitoolbox:makeLegacySidecarFixture:fixturesExist', ...
        ['Reference files already exist (e.g. %s). They must not be regenerated; ' ...
         'pass ''Overwrite'', true only if the task explicitly allows it.'], existing{1});
end

if ~isfolder(outFolder)
    mkdir(outFolder);
end

fid = fopen(targets{1}, 'w', 'ieee-le');
assert(fid >= 0, 'Cannot create legacy sidecar fixture.');
fwrite(fid, single(1:120), 'single');
fclose(fid);

datSize = [6 5];
datLength = 4;
Freq = 10;
Datatype = 'single';
dim_names = {'Y', 'X', 'T'};
save(targets{2}, 'datSize', 'datLength', 'Freq', 'Datatype', 'dim_names');

end
