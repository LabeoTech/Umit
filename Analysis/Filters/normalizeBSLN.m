function [outData, metaData] = normalizeBSLN(data, metaData, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NORMALIZEBSLN  Normalize imaging data by baseline (DeltaR/R0)
%
% Supports two modes:
%   1) STANDARD MODE (in-memory)
%      - Triggered when 'data' is a numeric array
%      - Baseline normalization done in RAM
%
%   2) LOW-RAM MODE (streaming / hybrid)
%      - Triggered when 'data' is a filename (.dat)
%      - Data are read in temporal chunks
%      - Output is written back to disk
%
% Inputs:
%   data : 3D (Y,X,T) or 4D (E,Y,X,T) numeric array or .dat filename
%   metaData : struct or matlab.io.MatFile with associated metadata
%   opts (optional struct):
%       - baseline_sec : 'auto' or numeric (baseline period in sec)
%       - b_centerAtOne : boolean, TRUE to center at 1 (default: false)
%
% Outputs:
%   outData : normalized array (if input was numeric) or output filename (if low-RAM)
%   metaData : updated metadata structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Defaults:
default_Output = 'normBSLN.dat'; %#ok This is here only as a reference for PIPELINEMANAGER.m.
default_opts = struct('baseline_sec','auto', 'b_centerAtOne', false);
opts_values = struct('baseline_sec', {{'auto',Inf}}, 'b_centerAtOne', [true,false]);%#ok. This is here only as a reference for PIPELINEMANAGER.m.
% Parse inputs
p = inputParser;
addRequired(p,'data',@(x) isnumeric(x) | ischar(x));
addRequired(p,'metaData', @(x) isa(x,'matlab.io.MatFile') | isstruct(x));
addOptional(p,'opts',default_opts,@(x) isstruct(x) && ~isempty(x));
parse(p,data,metaData,varargin{:});
opts     = p.Results.opts;
outData  = p.Results.data;
metaData = p.Results.metaData;
clear p data

% Validate dimensions
assert(all(ismember({'Y','X','T'}, metaData.dim_names)), ...
    'umIToolbox:normalizeBSLN:WrongInput', ...
    'Data must have dimensions "Y","X","T" (and optionally "E").');

% Determine baseline period
if any(strcmpi('E',metaData.dim_names))
    % 4D data (E,Y,X,T)
    if ~strcmpi(opts.baseline_sec,'auto')
        bsln_sec = single(opts.baseline_sec);
    else
        bsln_sec = metaData.preEventTime_sec;
    end
else
    % 3D data (Y,X,T)
    if strcmpi(opts.baseline_sec,'auto')
        bsln_sec = 0.2*metaData.datLength/metaData.Freq;
        warning('Setting baseline to first 20%% of frames (%.2f sec).', bsln_sec);
    else
        bsln_sec = single(opts.baseline_sec);
    end
end

% Execute RAM safe mode if input is filename
if ischar(outData)
    % RAM safe processing
    outData = normalizeBSLN_lowRAMmode(outData, metaData, bsln_sec, opts.b_centerAtOne);
    
else
    % Standard in-memory mode
    if any(strcmpi('E',metaData.dim_names))
        % 4D data
        bsln = median(outData(:,:,:,1:round(bsln_sec*metaData.Freq)), 4, 'omitnan');
    else
        % 3D data
        bsln = median(outData(:,:,1:round(bsln_sec*metaData.Freq)), 3, 'omitnan');
    end
    
    outData = (outData - bsln) ./ bsln;
    if opts.b_centerAtOne
        outData = outData + 1;
    end
end

end

%% ========================================================================
function outFile = normalizeBSLN_lowRAMmode(inFile, metaData, baseline_sec,b_centerAtOne)
% NORMALIZEBSLN_LOWRAMMODE  Low RAM normalization of image time series data
% Supports both continuous (3D) and event-split (4D, E,Y,X,T) data
%
% Inputs:
%   inFile   : Input .dat file containing the data
%   metaData : Associated meta data struct
%   opts     : Options struct
%       - baseline_sec : baseline period in seconds or 'auto'
%       - b_centerAtOne: logical, center data at 1
%
% Output:
%   outFile  : Output .dat file with normalized data

% Output filename
outFile = fullfile(fileparts(inFile), 'NORMALIZED_DATA.dat');


% Get dimensions
dimNames = metaData.dim_names;
dims     = [metaData.datSize metaData.datLength];
if any(strcmpi('E', dimNames))
    % Event-split data: [E,Y,X,T]
    hasEvents = true;
    Ne = dims(strcmpi('E', dimNames));
    Ny = dims(strcmpi('Y', dimNames));
    Nx = dims(strcmpi('X', dimNames));
    Nt = dims(strcmpi('T', dimNames));
else
    % Continuous data: [Y,X,T]
    hasEvents = false;
    Ny = dims(strcmpi('Y', dimNames));
    Nx = dims(strcmpi('X', dimNames));
    Nt = dims(strcmpi('T', dimNames));
end

disp('Starting RAM safe baseline normalization...');
% Preallocate output dat file
preallocateDatFile(outFile, metaData);

fidIn  = fopen(inFile,'r');
cIn = onCleanup(@() safeFclose(fidIn));
fidOut = fopen(outFile,'r+');
cOut= onCleanup(@() safeFclose(fidOut));

if hasEvents
    % ----------------------------
    % Event-split data (E,Y,X,T)
    % ----------------------------
    sz = [metaData.datSize, metaData.datLength];
    Ne = sz(1);
    Ny = sz(2);
    Nx = sz(3);
    Nt = sz(4);
    
    for e = 1:Ne
        % Read single trial:
        
        trialData = readTrial(fidIn,e,sz,'single');
        
        % Compute baseline
        if strcmpi(baseline_sec,'auto')
            bslnFrames = 1:round(0.2*Nt);  % first 20% frames
        else
            bslnFrames = 1:round(baseline_sec*metaData.Freq);
        end
        bsln = median(trialData(:,:,bslnFrames), 3, 'omitnan');
        
        % Normalize
        trialData = (trialData - bsln)./bsln;
        if b_centerAtOne
            trialData = trialData + 1;
        end
        
        % Write to output
        writeTrial_YXTE(fidOut,e,trialData,sz([2 3 4 1]),'single');
        clear trialData bsln
    end
    fclose(fidOut);
    fclose(fidIn);
    permuteDat_YXTE_to_EYXT_inplace(outFile,sz([2 3 4 1]),'single')
    
else
    % ----------------------------
    % Continuous data (Y,X,T)
    % ----------------------------
    % Estimate if baseline can fit in RAM
    if strcmpi(baseline_sec,'auto')
        nBaseFrames = round(0.2*Nt);
    else
        nBaseFrames = round(baseline_sec*metaData.Freq);
    end
    
    % Try to load baseline into memory
    try
        fseek(fidIn,0,'bof');
        baseData = fread(fidIn, Ny*Nx*nBaseFrames, metaData.Datatype);
        baseData = reshape(baseData, Ny, Nx, nBaseFrames);
        bsln = median(baseData,3,'omitnan');
        clear baseData
    catch
        warning('Low RAM available: using mean instead of median for baseline calculation');
        % Fallback: calculate mean in an iterative way
        fseek(fidIn,0,'bof');
        bsln = zeros(Ny,Nx,'single');
        for t = 1:nBaseFrames
            frame = fread(fidIn, Ny*Nx, metaData.Datatype);
            frame = reshape(frame, Ny, Nx);
            bsln = bsln + frame;
        end
        bsln = bsln / nBaseFrames;
    end
    
    % Normalize all frames
    fseek(fidIn,0,'bof');
    for t = 1:Nt
        frame = fread(fidIn, Ny*Nx, metaData.Datatype);
        frame = reshape(frame, Ny, Nx);
        frame = (frame - bsln)./bsln;
        if b_centerAtOne
            frame = frame + 1;
        end
        fwrite(fidOut, frame, 'single');
    end
    
    fclose(fidIn);
    fclose(fidOut);
end

disp('Finished RAM safe normalization by baseline.');
end

