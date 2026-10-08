function outFile = run_ConvertToTiff(data, SaveFolder)
%RUN_CONVERTTOTIFF Export image .dat data to TIFF file(s).
%
%   outFile = run_ConvertToTiff(data, SaveFolder)
%   info    = run_ConvertToTiff('pipelineInfo')
%
%   Supported input:
%       .dat file (name or path; a bare name is also looked up in
%       SaveFolder) with layout Y-X-T, Y-X-T-E, or a single Y-X frame.
%       Arrays, UMT structs, and .umt files are not supported.
%
%   Behavior:
%       - Y-X-T (and Y-X) data produce one TIFF file with one page per
%         frame.
%       - Event-split Y-X-T-E data produce one TIFF file per E slice. The E
%         axis is labeled from the events.mat in the file's folder
%         (resolveDatEventMapping). Without an events.mat the slices are
%         labeled as repetitions 1..E of one condition; an E axis that
%         cannot be matched to events.mat is rejected, since its file names
%         would carry wrong labels. A companion text file named
%         <baseName>_info.txt is created using CSV formatting. It lists the
%         generated TIFF file name, condition name, and repetition index.
%       - Pages are written as 32-bit floating-point values. The .dat file
%         is streamed in blocks of frames, so its size is limited by disk,
%         not RAM (Low-RAM mode is always on; the block size follows the
%         available RAM). A file larger than 3.9 GB is written as BigTIFF.
%
%   Output:
%       outFile - File manifest cell array containing the generated file
%                 name(s) saved in SaveFolder. Names are derived from the
%                 input file:
%                     Y-X-T, Y-X   -> img_<inputStem>.tif
%                     Y-X-T-E      -> one img_<inputStem>_C<id>_R<rep>.tif
%                                     per E slice, plus
%                                     img_<inputStem>_info.txt

default_Output = {'img_out.tif'};

if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) ...
        && strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outFile = localPipelineInfo();
    return
end

p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'data');
addRequired(p, 'SaveFolder', @isfolder);
parse(p, data, SaveFolder);

SaveFolder = p.Results.SaveFolder;

[dataFile, baseName] = iResolveDatFile(data, SaveFolder);
datInfo = loadMetaData(dataFile);
assertDatLayout(datInfo, {{'Y','X','T'}, {'Y','X'}, {'Y','X','T','E'}}, 'run_ConvertToTiff');

dimNames = cellstr(string(datInfo.dimNames));
hasE = any(strcmpi(dimNames, 'E'));

if hasE
    mapping = resolveDatEventMapping(datInfo, fileparts(dataFile));
    assert(ismember(mapping.status, {'matched','aggregated','noEvents'}), ...
        'Umitoolbox:run_ConvertToTiff:EventMappingMismatch', ...
        'The E axis of "%s" cannot be matched to events.mat: %s', dataFile, mapping.message);
end

slabIn = spatialSlabIO('open', dataFile, 'Info', datInfo);
closeIn = onCleanup(@() spatialSlabIO('close', slabIn));

nT = datAxisSize(datInfo, 'T');
if nT < 1
    nT = 1;   % single Y-X frame: one TIFF page
end
nE = datAxisSize(datInfo, 'E');
if nE < 1
    nE = 1;
end

if ~hasE
    tifName = [baseName '.tif'];
    iStreamTiff(fullfile(SaveFolder, tifName), slabIn, 1:nT);
    outFile = {tifName};
    return
end

eventInfo = mapping.eventInfo;
outFile = cell(1, nE);
infoRows = cell(nE, 3);
for iE = 1:nE
    tifName = sprintf('%s_C%d_R%d.tif', baseName, ...
        eventInfo.eventID(iE), eventInfo.repetitionIndex(iE));
    % E is the last axis, so the frames of slice iE are one contiguous run
    % of the flattened trailing axes.
    iStreamTiff(fullfile(SaveFolder, tifName), slabIn, (iE-1)*nT + (1:nT));
    outFile{iE} = tifName;
    infoRows{iE,1} = tifName;
    infoRows{iE,2} = char(string(eventInfo.eventName(iE)));
    infoRows{iE,3} = eventInfo.repetitionIndex(iE);
end

infoName = [baseName '_info.txt'];
iWriteEventInfoText(fullfile(SaveFolder, infoName), infoRows);
outFile{end+1} = infoName;

    function info = localPipelineInfo()
        info = PipelineManager.createPipelineInfo(mfilename, ...
            'Export image .dat data to TIFF file(s).');
        info.version = '2.0.0';

        info = PipelineManager.addInput(info, 'data', ...
            {'ImageTimeSeries','ProcessedData','UnknownDataType'}, ...
            ['Image .dat file (Y-X-T or event-split Y-X-T-E) to export as TIFF. ' ...
             'PipelineManager passes the file so the source stem, and therefore ' ...
             'the exported file identity, is independent of RAM mode.'], ...
            'kind', 'input', 'position', 1, 'callType', 'positional', ...
            'isData', true, 'supportsFile', true, 'dataMode', 'file');

        info = PipelineManager.addInput(info, 'SaveFolder', 'SaveFolder', ...
            'Folder where TIFF file(s) will be saved.', ...
            'kind', 'input', 'position', 2, 'callType', 'positional', ...
            'isData', false);

        info = PipelineManager.addOutput(info, 'outFile', 'ImageTimeSeries', 'file', ...
            ['Generated TIFF file manifest saved in SaveFolder. The base name is ' ...
             '''img_<inputStem>''; event-split input writes one ' ...
             '''<baseName>_C<eventID>_R<repetition>.tif'' per E slice instead of a ' ...
             'single ''<baseName>.tif''. The declared name is the non-event case. ' ...
             'Read the returned manifest for the names actually written.'], ...
            default_Output, 1, 'isData', true, 'saveFileName', '');
    end
end

function [dataFile, baseName] = iResolveDatFile(data, SaveFolder)
%IRESOLVEDATFILE Resolve the input to an existing .dat path and base name.

assert(ischar(data) || (isstring(data) && isscalar(data)), ...
    'Umitoolbox:run_ConvertToTiff:UnsupportedInputType', ...
    'run_ConvertToTiff accepts only a .dat file name or path.');

dataFile = char(string(data));
if ~isfile(dataFile)
    altPath = fullfile(SaveFolder, dataFile);
    if isfile(altPath)
        dataFile = altPath;
    else
        error('Umitoolbox:run_ConvertToTiff:InputFileNotFound', ...
            'Input file "%s" was not found.', char(string(data)));
    end
end

[~, stem, ext] = fileparts(dataFile);
assert(strcmpi(ext, '.dat'), ...
    'Umitoolbox:run_ConvertToTiff:UnsupportedInputFile', ...
    'Unsupported input file extension "%s". Only .dat files are supported.', ext);
baseName = ['img_' stem];
end

function iStreamTiff(filePath, slabIn, frameIdx)
%ISTREAMTIFF Write the given frames of an open .dat as one TIFF stack.
%   Frames are read in blocks sized from the available RAM and written as
%   they arrive, so only one block is in memory.

Ny = slabIn.Ny;
Nx = slabIn.Nx;
nFrames = numel(frameIdx);

% Fast_Tiff_Write uses 32-bit offsets; switch to BigTIFF past ~4 GB.
if double(Ny) * Nx * nFrames * 4 < 3.9e9
    writer = Fast_Tiff_Write(filePath, 1, 0);
else
    writer = Fast_BigTiff_Write(filePath, 1, 0);
end
[~, tifName, tifExt] = fileparts(filePath);
cleanupWriter = onCleanup(@() iCloseWriter(writer));

nBlocks = calculateMaxChunkSize(double(Ny) * Nx * nFrames * 4, 2, 0.1);
framesPerBlock = max(1, ceil(nFrames / nBlocks));
nBlocks = ceil(nFrames / framesPerBlock);

for iBlock = 1:nBlocks
    sel = ((iBlock-1)*framesPerBlock + 1):min(iBlock*framesPerBlock, nFrames);
    block = single(spatialSlabIO('read', slabIn, 1:Nx, frameIdx(sel)));
    for iFrame = 1:size(block, 3)
        writer.WriteIMG(block(:,:,iFrame)');
    end
    if nBlocks > 1
        fprintf('%s%s: %d/%d frames written\n', tifName, tifExt, sel(end), nFrames);
    end
end
end

function iCloseWriter(writer)
%ICLOSEWRITER Close a TIFF writer once, whether or not the write succeeded.
if ~writer.Closed
    writer.close();
end
end

function iWriteEventInfoText(filePath, rows)
%IWRITEEVENTINFOTEXT Write event export manifest in CSV-formatted text.

fid = fopen(filePath, 'w');
assert(fid ~= -1, 'Umitoolbox:run_ConvertToTiff:FileOpenFailed', ...
    'Failed to create "%s".', filePath);
cleanupObj = onCleanup(@() fclose(fid));

fprintf(fid, 'tiffFile,conditionName,repetitionIndex\n');
for iRow = 1:size(rows,1)
    fprintf(fid, '%s,%s,%d\n', rows{iRow,1}, rows{iRow,2}, rows{iRow,3});
end
end
