function slab = spatialSlabIO(mode, fid, Ny, Nx, Nt, xIdx, precision, slab)
% SPATIALSLABIO  Read or write a spatial X-slab from/to a raw [Y,X,T] binary file.
%
%   slab = SPATIALSLABIO('read', fid, Ny, Nx, Nt, xIdx, precision)
%   SPATIALSLABIO('write', fid, Ny, Nx, Nt, xIdx, precision, slab)
%
% The binary file is assumed to store a MATLAB-formatted array of size
% [Y,X,T] in column-major order.
%
% Inputs:
%   mode      : 'read' or 'write'
%   fid       : file identifier from fopen
%   Ny, Nx, Nt: dimensions of on-disk data [Y,X,T]
%   xIdx      : vector of X indices to read/write
%   precision : numeric class (e.g. 'single', 'double')
%
% Additional input (WRITE mode only):
%   slab      : array of size [Ny, numel(xIdx), Nt]
%
% Output (READ mode only):
%   slab      : array of size [Ny, numel(xIdx), Nt]
%
% Notes:
%   - If xIdx is contiguous, fast block I/O is used.
%   - Otherwise, a safe column-wise fallback is used.
%   - The file must already be preallocated in WRITE mode.
%
% See also: fread, fwrite, fseek

    bytes = getByteSize(precision);
    nX = numel(xIdx);

    isRead  = strcmpi(mode, 'read');
    isWrite = strcmpi(mode, 'write');

    if ~(isRead || isWrite)
        error('spatialSlabIO:InvalidMode', ...
              'Mode must be ''read'' or ''write''.');
    end

    % ---------------------------------------------------------------------
    % Allocate / validate slab
    % ---------------------------------------------------------------------
    if isRead
        slab = zeros(Ny, nX, Nt, precision);
    else
        if nargin < 8
            error('spatialSlabIO:MissingInput', ...
                  'Slab input required in WRITE mode.');
        end

        if size(slab,1) ~= Ny || size(slab,2) ~= nX || size(slab,3) ~= Nt
            error('spatialSlabIO:SizeMismatch', ...
                  'Slab must be of size [Ny, numel(xIdx), Nt].');
        end
    end

    % ---------------------------------------------------------------------
    % Contiguous vs non-contiguous X indices
    % ---------------------------------------------------------------------
    isContiguous = all(diff(xIdx) == 1);

    if isContiguous
        % Fast path (block I/O)
        x0 = xIdx(1);

        for t = 1:Nt
            offset = ((t-1)*Nx*Ny + (x0-1)*Ny) * bytes;
            fseek(fid, offset, 'bof');

            if isRead
                slab(:,:,t) = fread(fid, [Ny, nX], precision);
            else
                fwrite(fid, slab(:,:,t), precision);
            end
        end

    else
        % Safe fallback (column-wise)
        for t = 1:Nt
            tOffset = (t-1) * Ny * Nx * bytes;

            for xi = 1:nX
                x = xIdx(xi);
                offset = tOffset + (x-1) * Ny * bytes;
                fseek(fid, offset, 'bof');

                if isRead
                    slab(:,xi,t) = fread(fid, Ny, precision);
                else
                    fwrite(fid, slab(:,xi,t), precision);
                end
            end
        end
    end
end
