function outFile = pmFreshInvalidMetadataInitializer(RawFolder, SaveFolder)
%PMFRESHINVALIDMETADATAINITIALIZER Test fixture for fresh-folder validation.

if nargin == 1 && strcmpi(char(string(RawFolder)), 'pipelineInfo')
    outFile = PipelineManager.createPipelineInfo(mfilename, ...
        'Write malformed acquisition metadata and return normally.');
    outFile.freshSaveFolderRole = 'acquisition-initializer';
    outFile = PipelineManager.addInput(outFile, 'RawFolder', 'RawFolder', ...
        'Unused fixture raw folder.', 'kind', 'input', 'position', 1, ...
        'callType', 'positional', 'isData', false);
    outFile = PipelineManager.addInput(outFile, 'SaveFolder', 'SaveFolder', ...
        'Fixture output folder.', 'kind', 'input', 'position', 2, ...
        'callType', 'positional', 'isData', false);
    outFile = PipelineManager.addOutput(outFile, 'outFile', ...
        'ImageTimeSeries', 'file', 'Fixture data file.', ...
        'partial.dat', 1, 'isData', true, 'saveFileName', '');
    return
end

assert(isfolder(RawFolder));
assert(isfolder(SaveFolder));

malformed = struct('junk', 42);
save(fullfile(SaveFolder, 'AcqInfos.mat'), '-struct', 'malformed');

fid = fopen(fullfile(SaveFolder, 'partial.dat'), 'w');
assert(fid ~= -1);
fwrite(fid, single(1), 'single');
fclose(fid);
outFile = 'partial.dat';
end
