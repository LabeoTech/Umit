function outFile = importFromTif(RawFolder,SaveFolder,varargin)
% IMPORTFROMTIF is a data import function that reads single-channel imaging
% data stored in TIF format. The function looks for any .tif (e.g. fluo.tif)
% file inside the "RawFolder" and the associate info file (e.g. fluo.txt)
% containing the recording meta data such as sample rate and exposure time.

% Defaults:
default_Output = {'fluo.dat', 'red.dat', 'green.dat', 'yellow.dat', 'speckle.dat'};%#ok This is here only as a reference for PIPELINEMANAGER.m. The real outputs will be stored in OUTFILE.
default_opts = struct('BinningSpatial', 1, 'BinningTemp', 1,'FrameRateHz',1,'ExposureTimeMs',0.2,'BackupOption','Erase');
opts_values = struct('BinningSpatial', 2.^[0:4], 'BinningTemp',2.^[0:4],'FrameRateHz',[eps,Inf],'ExposureTimeMs',[eps,Inf],'BackupOption',{{'Erase','genBackup'}});%#ok  % This is here only as a reference for PIPELINEMANAGER.m.
% Arguments validation:
p = inputParser;
addRequired(p, 'RawFolder', @isfolder);
addRequired(p, 'SaveFolder', @isfolder);
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
parse(p,RawFolder,SaveFolder,varargin{:});
% Initialize Variables:
RawFolder = p.Results.RawFolder;
SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;
outFile = {};
clear p
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Control for existing files and create backup or erase them:
genBackupFolder(SaveFolder,opts.BackupOption);
% Save AcqInfo.mat file in SaveFolder with
% Read all TIF inside "RawFolder":
tif_list = dir(fullfile(RawFolder,'*.tif'));
% Find Image sequence files:\
% Here, all TIF files ending with a number will be concatenated as a single
% one in numerical ascending order (e.g. red_000, red_001...; green-1,
% green-2...). Common delimiters {-,_, } will be used to parse the file
% names.
exp = '^(.*?)[_\-\s]*(\d*).tif$';
seqFiles = regexp({tif_list.name}', exp,'tokens');
seqNames = vertcat(seqFiles{:}); seqNames = vertcat(seqNames{:});
% Treat single .TIF files as a sequence of files with N = 1:
seqNames(cellfun(@isempty,seqNames(:,2)),2) = {'00001'};
uniqNames = unique(seqNames(:,1));
%
w = waitbar(0, 'Importing frames from TIFF files...');
w.Children(1).Title.Interpreter = 'none';
% Import image sequence files:
for ii = 1:length(uniqNames)
    waitbar(0,w,['Importing frames from file sequence: ' uniqNames{ii} '...']);
    % Find numerical sequence for the current file:
    indx_name = find(strcmp(seqNames(:,1),uniqNames{ii}));
    [~,indx_num]=sort(str2double(seqNames(indx_name,2)));
    fileSeq = tif_list(indx_name); fileSeq = fileSeq(indx_num);
    % Create .dat file with first .tif file from the sequence:
    data = readTiff(fullfile(fileSeq(1).folder,fileSeq(1).name));
    % Apply binning
    data = binData(data,opts.BinningSpatial,opts.BinningTemp);
    % Create AcqInfo.mat file:
    AcqInfo = struct.empty(0,1);
    if ~isfile(fullfile(SaveFolder,'AcqInfos.mat'))
        AcqInfo = genAcqInfo(size(data,1),size(data,2),size(data,3),opts.FrameRateHz/opts.BinningTemp,opts.ExposureTimeMs);
    end
    % Save first file:
    saveFilename = [uniqNames{ii},'.dat'];
    saveData(fullfile(SaveFolder,saveFilename),data,AcqInfo);
    % Append the rest of the sequence to the existing .dat file:
    for ind = 2:length(fileSeq)
        data = readTiff(fullfile(fileSeq(ind).folder,fileSeq(ind).name));
        % Apply binning
        data = binData(data,opts.BinningSpatial,opts.BinningTemp);
        % Append the data to the .dat file
        saveData(fullfile(SaveFolder,saveFilename),data,'Append',true);
    end
    outFile = vertcat(outFile,saveFilename);
end
close(w);
disp('Finished importFromTif.')
end


% Local function(s)
function data= binData(data,SpatBin,TempBin)
% BINDATA applies spatial and temporal binning to the data.

% Spatial Binning
if ( SpatBin > 1 )
    data = imresize(data,1/SpatBin);
end

% Temporal Binning
if ( TempBin > 1 )
    data = imresize3(data, [size(data,1), size(data,2),...
        size(data,3)/TempBin], 'linear');
end
end

function out = readTiff(tif_file)
% READTIFF reads the TIF file and outputs the 3D array (single precisison).
info_tif = imfinfo(tif_file);

% Check if data has a valid format:
tmp = imread(tif_file, 'Info', info_tif(1));
errID = 'umIToolbox:importFromTif:InvalidFile';
errMsg = ['Invalid data format: This function accepts only unsigned integer data '...
    'with 8, 16 or 32 Bits.'];
assert(ismember(class(tmp),{'uint8', 'uint16', 'uint32'}),errID,errMsg);
% Import frames:
out = zeros([size(tmp),length(info_tif)],class(tmp));
for ind = 1:length(info_tif)
    out (:,:,ind) = imread(tif_file, 'Info', info_tif(ind));
end
% Transform data format to "single":
out = single(out);
end





