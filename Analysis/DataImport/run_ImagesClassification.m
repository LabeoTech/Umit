function outFile = run_ImagesClassification(RawFolder, SaveFolder, varargin)
% RUN_IMAGESCLASSIFICATION calls the function
% IMAGESCLASSIFICATION from the IOI library (LabeoTech).

default_Output = {'fluo_475.dat','fluo.dat', 'red.dat', 'green.dat', 'yellow.dat', 'speckle.dat'}; % This is here only as a reference for PIPELINEMANAGER.m . The real outputs will be stored in OUTFILE.
%%% Arguments parsing and validation %%%
p = inputParser;
% Raw folder:
addRequired(p, 'RawFolder', @isfolder);
% Save folder:
addRequired(p, 'SaveFolder', @isfolder);
% Optional Parameters:
% opts structure:
default_opts = struct('BinningSpatial', 1, 'BinningTemp', 1, 'b_SubROI', false, 'b_IgnoreStim', true);
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));

parse(p, RawFolder, SaveFolder, varargin{:});
%Initialize Variables:
RawFolder = p.Results.RawFolder;
SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;
outFile = {};
%%%%
% Get existing ImagesClassification files in directory:
cd(SaveFolder);
existing_ChanList  = dir('*.dat');
idxName = ismember({existing_ChanList.name}, default_Output);
existing_ChanList = existing_ChanList(idxName);
% Calls function from IOI library. Temporary for now.
ImagesClassification(RawFolder, SaveFolder, opts.BinningSpatial, opts.BinningTemp, opts.b_IgnoreStim, opts.b_SubROI)
% Get only new files created during ImagesClassification:
chanList = dir('*.dat');
idx = ismember({chanList.name}, default_Output);
chanList = chanList(idx);
idxName = ismember({chanList.name}, {existing_ChanList.name});
idxDate = ismember([chanList.datenum], [existing_ChanList.datenum]);
idxNew = ~all([idxName; idxDate],1);
chanList = {chanList(idxNew).name};
for i = 1:length(chanList)
    chanName = chanList{i};    
    MetaDataFileName = strrep(chanName, '.dat', '.mat');
    a =  matfile(fullfile(SaveFolder,MetaDataFileName), 'Writable', true);
    a.fileUUID = char(java.util.UUID.randomUUID);
    a.Datatype = 'single';
    a.datName = 'data';
    a.dim_names = {'Y', 'X', 'T'};
    % TEMPORARY FIX.
    filePath = fullfile(SaveFolder,chanName);
    newmDfile = strrep(MetaDataFileName, '.mat', '_info.mat');
    statusMat = movefile(MetaDataFileName, newmDfile);
    if ~statusMat
        disp(['Failed to rename ' MetaDataFileName]);
    else
        % Flip X,Y image axes:
        [mData, metaData] = mapDatFile(filePath);
        mData.Writable = true;
        metaData.Properties.Writable = true;
        data = mData.Data.data;
        mData.Data.data = permute(data, [2 1 3]);
        metaData.datSize = fliplr(metaData.datSize);
    end
    outFile = [outFile, chanName];
end
end