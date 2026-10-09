function makeDatFixturesV1(outFolder, varargin)
%MAKEDATFIXTURESV1 Generate the version-1 .dat header reference files.
%
%   makeDatFixturesV1()
%   makeDatFixturesV1(outFolder)
%   makeDatFixturesV1(outFolder, 'Overwrite', true)
%
%   Documents how the reference files used by TestDatHeader and TestMapDat
%   were produced. The files are generated once and must never be
%   regenerated after the .dat header Phase 1 task is accepted; the
%   function refuses to overwrite existing files unless 'Overwrite' is
%   true.
%
%   Default outFolder: the "v1" folder next to this file.
%
%   Files written (all data little-endian, column-major, deterministic):
%       yxt_single.dat      Y=6 X=5 T=4, single, 1:120, 10 Hz, 5 ms, 'green'
%       yx_uint16.dat       Y=6 X=5, uint16, 1:30, NaN Hz, NaN ms, 'amplitude'
%       yxte_single.dat     Y=4 X=3 T=5 E=2, single, 1:120, 10 Hz, 5 ms, 'green'
%       yxe_single.dat      Y=4 X=3 E=2 (no T), single, 1:24, NaN Hz, NaN ms, 'green'
%       yxt_incomplete.dat  as yxt_single.dat, write-complete bit 0
%       legacy/green.dat    headerless, Y=6 X=5 T=4, single, 1:120
%       legacy/AcqInfos.mat AcqInfoStream as written by importFromTif
%
%   Header bytes come from encodeDatHeader; data follows as a plain fwrite.

p = inputParser;
addOptional(p, 'outFolder', fullfile(fileparts(mfilename('fullpath')), 'v1'), ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'Overwrite', false, @(x) islogical(x) && isscalar(x));
if nargin < 1
    parse(p);
else
    parse(p, outFolder, varargin{:});
end
outFolder = char(p.Results.outFolder);
legacyFolder = fullfile(outFolder, 'legacy');

names = {'yxt_single.dat', 'yx_uint16.dat', 'yxte_single.dat', ...
    'yxe_single.dat', 'yxt_incomplete.dat'};
targets = [fullfile(outFolder, names), ...
    {fullfile(legacyFolder, 'green.dat'), fullfile(legacyFolder, 'AcqInfos.mat')}];
existing = targets(cellfun(@isfile, targets));
if ~isempty(existing) && ~p.Results.Overwrite
    error('Umitoolbox:makeDatFixturesV1:fixturesExist', ...
        ['Reference files already exist (e.g. %s). They must not be regenerated; ' ...
         'pass ''Overwrite'', true only if the task explicitly allows it.'], existing{1});
end

if ~isfolder(legacyFolder)
    mkdir(legacyFolder);
end

% --- new-format files ------------------------------------------------------
yxt = struct('dataClass', 'single', 'frameRateHz', 10, 'exposureMsec', 5, ...
    'channelName', 'green', 'dimNames', {{'Y', 'X', 'T'}}, 'dimSizes', [6 5 4], ...
    'writeComplete', true);
iWriteDat(fullfile(outFolder, 'yxt_single.dat'), yxt, single(1:120));

yx = struct('dataClass', 'uint16', 'frameRateHz', NaN, 'exposureMsec', NaN, ...
    'channelName', 'amplitude', 'dimNames', {{'Y', 'X'}}, 'dimSizes', [6 5], ...
    'writeComplete', true);
iWriteDat(fullfile(outFolder, 'yx_uint16.dat'), yx, uint16(1:30));

yxte = struct('dataClass', 'single', 'frameRateHz', 10, 'exposureMsec', 5, ...
    'channelName', 'green', 'dimNames', {{'Y', 'X', 'T', 'E'}}, 'dimSizes', [4 3 5 2], ...
    'writeComplete', true);
iWriteDat(fullfile(outFolder, 'yxte_single.dat'), yxte, single(1:120));

yxe = struct('dataClass', 'single', 'frameRateHz', NaN, 'exposureMsec', NaN, ...
    'channelName', 'green', 'dimNames', {{'Y', 'X', 'E'}}, 'dimSizes', [4 3 2], ...
    'writeComplete', true);
iWriteDat(fullfile(outFolder, 'yxe_single.dat'), yxe, single(1:24));

incomplete = yxt;
incomplete.writeComplete = false;
iWriteDat(fullfile(outFolder, 'yxt_incomplete.dat'), incomplete, single(1:120));

% --- legacy headerless file --------------------------------------------------
fid = fopen(fullfile(legacyFolder, 'green.dat'), 'w', 'ieee-le');
assert(fid >= 0, 'Cannot create legacy fixture.');
fwrite(fid, single(1:120), 'single');
fclose(fid);

AcqInfoStream = struct();
AcqInfoStream.Width = 5;
AcqInfoStream.Height = 6;
AcqInfoStream.Length = 4;
AcqInfoStream.FrameRateHz = 10;
AcqInfoStream.ExposureMsec = 5;
AcqInfoStream.ImportedChannels = struct( ...
    'DatFile', 'green.dat', 'Tag', 'green', 'Color', 'green', ...
    'Length', 4, 'FrameRateHz', 10, 'ExposureMsec', 5, 'CamIdx', 1);
save(fullfile(legacyFolder, 'AcqInfos.mat'), 'AcqInfoStream');

end

% =========================================================================
function iWriteDat(filePath, hdr, values)
headerBytes = encodeDatHeader(hdr);
fid = fopen(filePath, 'w', 'ieee-le');
assert(fid >= 0, 'Cannot create fixture %s.', filePath);
cleanupObj = onCleanup(@() fclose(fid));
fwrite(fid, headerBytes, 'uint8');
fwrite(fid, values, class(values));
clear cleanupObj
end
