function outData = normalizeLPF(data, metaData, varargin)
% NORMALIZELPF performs a data normalization by low-pass filtering.
% The filtering algorithm consists in creating two low-passed versions of the
% signal with a given cut-off frequency ("LowCutOffHz" and "HighCutOffHz") and
% subsequently subtracting the two filtered signals.
% Optionally, the subtracted signals can be normalized by the low cut-off signal
% to express the signal as DeltaR/R.
%
% This function is a wrapper of the IOI library function
% "NormalisationFiltering.m".
%
% Limitations:
%   The data must be an Image time series with dimensions:
%   {Y, X, T} or {E, Y, X, T}.
%
% Inputs:
%   data     : numeric array containing image time series
%   metaData : struct or matlab.io.MatFile with associated meta data
%   opts     : (optional) structure with fields:
%       - LowCutOffHz   : low cut-off frequency (Hz)
%       - HighCutOffHz  : high cut-off frequency (Hz)
%       - Normalize     : true ? ?R/R, false ? ?R
%       - bApplyExpFit  : apply exponential decay correction
%
% Outputs:
%   outData  : filtered data, same size as input
%   metaData : unchanged meta data


default_Output = 'normLPF.dat'; %#ok  Required by PipelineManager
default_opts = struct('LowCutOffHz', 0.0083, 'HighCutOffHz', 1,'Normalize', true, 'bApplyExpFit', false);% Required by PipelineManager
opts_values = struct('LowCutOffHz', [0, Inf], 'HighCutOffHz', [eps, Inf], 'Normalize', [false, true], 'bApplyExpFit', [true, false]); %#ok  Required by PipelineManager

p = inputParser;
addRequired(p, 'data', @(x) isnumeric(x) | ischar(x));
addRequired(p, 'metaData', @(x) isa(x,'matlab.io.MatFile') || isstruct(x));
addOptional(p, 'opts', default_opts, @(x) isstruct(x) && ~isempty(x));
parse(p, data, metaData, varargin{:});

outData  = p.Results.data;
metaData = p.Results.metaData;
opts     = p.Results.opts;
clear p

% -------------------------------------------------------------------------
% Validate dimensions
% -------------------------------------------------------------------------
errID  = 'umIToolbox:normalizeLPF:InvalidInput';
errMsg = 'Data must have dimensions {"Y","X","T"} or {"E","Y","X","T"}.';
assert(all(ismember({'Y','X','T'}, metaData.dim_names)), errID, errMsg);

idxE = strcmpi(metaData.dim_names, 'E');
b_HasEvents = any(idxE);

% -------------------------------------------------------------------------
% Replace NaNs with zeros (filter safety)
% -------------------------------------------------------------------------
idx_nan = isnan(outData);
outData(idx_nan) = 0;

% -------------------------------------------------------------------------
% Validate cut-off frequencies
% -------------------------------------------------------------------------
Fs = metaData.Freq;

if opts.LowCutOffHz < 0 || opts.LowCutOffHz > Fs/2
    error(errID, ...
        'LowCutOffHz must be between 0 and Nyquist frequency.');
end

if opts.HighCutOffHz <= 0 || opts.HighCutOffHz > Fs/2
    error(errID, ...
        'HighCutOffHz must be > 0 and <= Nyquist frequency.');
end

if opts.HighCutOffHz < opts.LowCutOffHz
    error(errID, ...
        'HighCutOffHz must be >= LowCutOffHz.');
end

% ========================================================================
% CASE 1 — DATA IS ON DISK (.dat)
% ========================================================================
if ischar(outData)
       
    inFile  = outData;
    outFile = fullfile(fileparts(inFile), 'NORMALIZEDDATA.dat');
    outData = outFile;
    sz = [metaData.datSize, metaData.datLength];        
    % ---------------------------------------------------------------------
    % EVENT-SPLIT DATA ON DISK
    % ---------------------------------------------------------------------
    if b_HasEvents
        Ne = sz(1);
%         Ny = sz(2);
%         Nx = sz(3);
%         Nt = sz(4);
        
        % Preallocate output file
        preallocateDatFile(outFile, metaData);
        fid_out = fopen(outFile,'r+');   
        c_out = onCleanup(@() safeFclose(fid_out));
        fid_in = fopen(inFile, 'r');
        c_in = onCleanup(@() safeFclose(fid_in));
        
        disp('Filtering event-split data (disk-based)...');
        
        for e = 1:Ne
            % Read one full trial (safe in RAM)
            trial = readTrial(fid_in,e,sz,'single');
                        
            % Run filtering on single trial
            trial = NormalisationFiltering( ...
                pwd, trial, ...
                opts.LowCutOffHz, ...
                opts.HighCutOffHz, ...
                opts.Normalize, ...
                opts.bApplyExpFit, ...
                Fs);
            
            % Write back   
            writeTrial_YXTE(fid_out,e,trial,sz([2 3 4 1]),'single')            
        end
        
        fclose(fid_in);
        fclose(fid_out);
        % Fix Permutted output file
        permuteDat_YXTE_to_EYXT_inplace(outFile,sz([2 3 4 1]),'single');
        
        % ---------------------------------------------------------------------
        % NO EVENTS — LET NormalisationFiltering HANDLE HYBRID MODE
        % ---------------------------------------------------------------------
    else               
        NormalisationFiltering( ...
            pwd, inFile, ...
            opts.LowCutOffHz, ...
            opts.HighCutOffHz, ...
            opts.Normalize, ...
            opts.bApplyExpFit, ...
            Fs,outFile);
    end
    
    return
end

% ========================================================================
% CASE 2 — DATA IS IN MEMORY
% ========================================================================
idx_nan = isnan(outData);
outData(idx_nan) = 0;

disp('Filtering data...');

if b_HasEvents
    % Reorder to E,Y,X,T
    dimorder = 1:ndims(outData);
    dimorder = [dimorder(idxE), dimorder(~idxE)];
    outData = permute(outData, dimorder);
    
    for e = 1:size(outData,1)
        outData(e,:,:,:) = NormalisationFiltering( ...
            pwd, squeeze(outData(e,:,:,:)), ...
            opts.LowCutOffHz, ...
            opts.HighCutOffHz, ...
            opts.Normalize, ...
            opts.bApplyExpFit, ...
            Fs);
    end
    
    outData = ipermute(outData, dimorder);
else
    outData = NormalisationFiltering( ...
        pwd, outData, ...
        opts.LowCutOffHz, ...
        opts.HighCutOffHz, ...
        opts.Normalize, ...
        opts.bApplyExpFit, ...
        Fs);
end

outData(idx_nan) = NaN;
disp('Finished with temporal filter.');

end