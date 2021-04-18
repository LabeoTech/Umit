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
subjName = object.MyParent.MyParent.ID;
%
cd(object.RawFolder)
mFiles = dir('*.mat');
idx = contains({mFiles.name}, subjName);

if sum(idx)>0
    metaFile = mFiles(idx).name;
    object.MetaDataFileName = metaFile;
else 
    disp(['Could not file MetaData File in ' object.RawFolder]);
end

end