function permuteDat_YXTE_to_EYXT_inplace(srcFile, sz, datatype)
% PERMUTEDAT_YXTE_TO_EYXT_INPLACE
% Permutes a .dat file from [Y,X,T,E] layout to [E,Y,X,T] layout on disk
% using a temporary copy (safe, low-RAM).
%
% Steps:
%   1) Open source file (read-only)
%   2) Create temp copy
%   3) Write permuted data to temp file
%   4) Overwrite source with temp
%   5) Delete temp
%
% Inputs:
%   srcFile  : path to source .dat file
%   sz       : [Y X T E]
%   datatype : e.g. 'single', 'double'

% -------------------------------
% Sizes
% -------------------------------
nY = sz(1);
nX = sz(2);
nT = sz(3);
nE = sz(4);

bytes = getByteSize(datatype);

YX  = nY * nX;
YXT = YX * nT;
EYX = nE * YX;

% -------------------------------
% Temp file
% -------------------------------
[tmpPath, tmpName, tmpExt] = fileparts(srcFile);
tmpFile = fullfile(tmpPath, [tmpName '_perm_tmp' tmpExt]);

% Safety: remove stale temp file
if exist(tmpFile,'file')
    delete(tmpFile);
end

% -------------------------------
% Open files
% -------------------------------
fidIn  = fopen(srcFile,'r');
assert(fidIn > 0, 'Could not open source file');

fidOut = fopen(tmpFile,'w');
assert(fidOut > 0, 'Could not create temp file');
fprintf('Permutting file to EYXT dimension order...\n')
% cleanupObj = onCleanup(@() cleanupFiles(fidIn, fidOut, tmpFile));

% -------------------------------
% Main permutation
% -------------------------------
for t = 1:nT

    % Read one full [Y,X,E] slab
    slabYXE = zeros(nY, nX, nE, datatype);

    for e = 1:nE
        offset = ((e-1)*YXT + (t-1)*YX) * bytes;
        fseek(fidIn, offset, 'bof');
        slabYXE(:,:,e) = fread(fidIn, [nY, nX], ['*' datatype]);
    end

    % Permute to [E,Y,X]
    slabEYX = permute(slabYXE, [3 1 2]);

    % Write contiguous block
    offset = (t-1) * EYX * bytes;
    fseek(fidOut, offset, 'bof');
    fwrite(fidOut, slabEYX, datatype);
end

% -------------------------------
% Close files
% -------------------------------
fclose(fidIn);
fclose(fidOut);

% -------------------------------
% Overwrite source
% -------------------------------
movefile(tmpFile, srcFile, 'f');
disp('Finished permutation');
end

