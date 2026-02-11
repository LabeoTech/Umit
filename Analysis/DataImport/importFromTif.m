function outFile = importFromTif(RawFolder, SaveFolder, varargin)
% IMPORTFROMTIF  Import single-channel imaging data from TIFF files.
%
% This function reads imaging data stored in TIFF format and converts it
% into raw binary (.dat) files compatible with the umIToolbox pipeline.
%
% The function supports LOW-RAM streaming by reading TIFF frames in
% time chunks and appending them sequentially to a preallocated .dat file.
%
% Supported input formats:
%   - uint8, uint16, uint32 TIFF stacks
%
% -------------------------------------------------------------------------
% Inputs:
%
%   RawFolder :
%       Folder containing TIFF (.tif/.tiff) files and associated .txt
%       metadata files.
%
%   SaveFolder :
%       Destination folder for generated .dat files.
%
%   opts (optional struct) with fields:
%       - BinningSpatial : spatial binning factor (default = 1)
%       - BinningTemp    : temporal binning factor (default = 1)
%
% -------------------------------------------------------------------------
% Outputs:
%
%   outFile :
%       Cell array with names of generated .dat files.
%
% -------------------------------------------------------------------------
% Notes:
%   - Frames are streamed in T to minimize RAM usage
%   - fwrite is used for sequential disk writes
%   - Temporal binning is applied on-the-fly
%
% -------------------------------------------------------------------------

% Defaults (for PipelineManager reference)
default_Output = {'fluo.dat', 'red.dat', 'green.dat', 'yellow.dat', 'speckle.dat'}; %#ok
default_opts   = struct('BinningSpatial', 1, 'BinningTemp', 1);

% -------------------------------------------------------------------------
% Argument parsing
% -------------------------------------------------------------------------
p = inputParser;
addRequired(p, 'RawFolder', @isfolder);
addRequired(p, 'SaveFolder', @isfolder);
addOptional(p, 'opts', default_opts, @(x) isstruct(x) && ~isempty(x));
parse(p, RawFolder, SaveFolder, varargin{:});

opts       = p.Results.opts;
outFile    = {};
clear p

% -------------------------------------------------------------------------
% Locate valid TIFF/TXT pairs
% -------------------------------------------------------------------------
tif_list = [dir(fullfile(RawFolder,'*.tif')); dir(fullfile(RawFolder,'*.tiff'))];
txt_list = dir(fullfile(RawFolder,'*.txt'));

[~, tifNames] = cellfun(@fileparts, {tif_list.name}, 'UniformOutput', false);
[~, txtNames] = cellfun(@fileparts, {txt_list.name}, 'UniformOutput', false);

[idx, loc] = ismember(tifNames, txtNames);
if ~any(idx)
    error('umIToolbox:importFromTif:FileNotFound', ...
          'No valid TIFF/TXT file pairs found.');
end

tif_list = tif_list(idx);
txt_list = txt_list(loc(loc ~= 0));

disp('Importing data from TIFF files...')
w = waitbar(0, 'Importing TIFF files...');
w.Children(1).Title.Interpreter = 'none';

% ========================================================================
%                          Main loop
% ========================================================================

for i = 1:numel(tif_list)
    w.Name = ['Importing file ' num2str(i) ' / ' num2str(numel(tif_list))];drawnow()
    tif_file = fullfile(tif_list(i).folder, tif_list(i).name);
    txt_file = fullfile(txt_list(i).folder, txt_list(i).name);
   
    % ---------------------------------------------------------------------
    % Read TIFF info and validate format
    % ---------------------------------------------------------------------
    waitbar(0, w, ['Processing file: ' tif_list(i).name]);
    info_tif = imfinfo(tif_file);
    tmp = imread(tif_file, 1);

    if ~ismember(class(tmp), {'uint8','uint16','uint32','single'})
        error('umIToolbox:importFromTif:InvalidFile', ...
              'Only uint8, uint16, or uint32 TIFF files are supported.');
    end

    % ---------------------------------------------------------------------
    % Read acquisition metadata
    % ---------------------------------------------------------------------
    acq_info = ReadInfoFile(txt_list(i).folder, txt_list(i).name);

    % Frame geometry
    Ny = size(tmp,1);
    Nx = size(tmp,2);
    Nt = numel(info_tif);

    % Apply spatial binning size
    Ny_out = round(Ny / opts.BinningSpatial);
    Nx_out = round(Nx / opts.BinningSpatial);

    % Temporal binning
    Nt_out = floor(Nt / opts.BinningTemp);

    % ---------------------------------------------------------------------
    % Estimate RAM and chunking (slice in T)
    % ---------------------------------------------------------------------    
    dataBytes     = Ny * Nx * Nt * 4; % single precision
    nChunks       = calculateMaxChunkSize(dataBytes, 1,.15);

    tChunk = ceil(Nt / nChunks);

    % ---------------------------------------------------------------------
    % Prepare output DAT file and metadata
    % ---------------------------------------------------------------------
    datFileName = [acq_info.Illumination1.Color '.dat'];
    datFilePath = fullfile(SaveFolder, datFileName);

    extraParams.tExposure = acq_info.ExposureMsec;
    extraParams.Freq      = acq_info.FrameRateHz / opts.BinningTemp;

    metaData = genMetaData(zeros(Ny_out, Nx_out, Nt_out, 'single'), ...
                           {'Y','X','T'}, extraParams);
    metaData.datFile = datFilePath;
    save(strrep(datFilePath,'.dat','.mat'),'-struct','metaData')
    
    fid = fopen(datFilePath, 'w');
    c_out = onCleanup(@() safeFclose(fid));

    % ---------------------------------------------------------------------
    % Stream frames in time
    % ---------------------------------------------------------------------
    
        
    for c = 1:nChunks                
        
        tStart = (c-1)*tChunk + 1;
        tEnd   = min(c*tChunk, Nt);

        slab = zeros(Ny, Nx, tEnd-tStart+1, 'single');
        
        for t = 1:size(slab,3)
            if mod(t,50)== 0 || t == 1 || t == size(slab,3)
                waitbar(t/size(slab,3),w, ['Reading frames : ' tif_list(i).name]);
            end
            slab(:,:,t) = single(imread(tif_file,'Info', info_tif(tStart+t-1)));
        end

        % Apply binning
        if any(opts.BinningTemp > 1 | opts.BinningSpatial>1)
            waitbar(c/nChunks, w, ['Applying binning: ' tif_list(i).name]);
            slab = imresize3(slab, [Ny_out,Nx_out,Nt_out], 'nearest');                        
        end
                
        % -------------------------------------------------------------
        % Append directly to DAT (time-contiguous write)
        % -------------------------------------------------------------
        waitbar(c/nChunks,w,['Writing to .dat file: ' datFileName]);
        fwrite(fid, slab, 'single');
           
    end

    fclose(fid);
    
    outFile = [outFile, {datFileName}];
end

close(w);
disp('Finished importFromTif.')

end
