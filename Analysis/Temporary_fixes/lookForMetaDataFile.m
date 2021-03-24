function lookForMetaDataFile(object, SaveFolder)
% LOOKFORMETADATAFILE searchs for a .MAT with experiment meta data with name containing the
% subject ID inside OBJECT.RAWFOLDER and adds to OBJECT.METADATAFILE.

%%% Arguments parsing and validation %%%
p = inputParser;
% The input of the function must be a File , RawFolder or SaveFolder
% Imaging Object:
addRequired(p, 'object', @(x) isa(x,'Modality'));
addRequired(p, 'SaveFolder', @isfolder);
% Parse inputs:
parse(p,object, SaveFolder);
% Initialize Variables:
object = p.Results.object;
%%%%
% Find Subject ID from RawFolder string:
str = object.RawFolder;
idx = strfind(str,filesep);
str = str(idx(end-1)+1:idx(end)-1);
idx = strfind(str,'_');
subjName = str(1:idx(1)-1);
%
cd(object.RawFolder)
mFiles = dir('*.mat');
idx = contains({mFiles.name}, subjName);

if sum(idx)>0
    metaFile = mFiles(idx).name;
    object.MetaDataFile = fullfile(object.RawFolder, metaFile);
else 
    disp(['Could not file Meta data File in ' object.RawFolder]);
end

end