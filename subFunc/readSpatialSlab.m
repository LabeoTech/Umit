function slab = readSpatialSlab(fid, Ny, Nx, Nt, xIdx, precision)
% READSPATIALSLAB  Stream a spatial X-slab from a raw [Y,X,T] binary volume.
%
%   slab = READSPATIALSLAB(fid, Ny, Nx, Nt, xIdx, precision) reads the X-columns
%   specified in xIdx for all Y and T from a raw binary file storing a
%   MATLAB-formatted array of size [Y,X,T].
%
%   This function avoids OS page-cache growth and is safe for very large files.
%
%   If xIdx is contiguous, the data is read using fast block I/O.
%   Otherwise, it falls back to a safe column-wise reader.
%
%   See also: fopen, fread, fseek

    bytes = getPrecisionBytes(precision);
    nX = numel(xIdx);

    slab = zeros(Ny, nX, Nt, precision);

    % ---------------------------------------------------------------------
    % Check whether indices are contiguous
    % ---------------------------------------------------------------------
    if all(diff(xIdx) == 1)
        % Fast path (contiguous block)

        x0 = xIdx(1);

        for t = 1:Nt
            % Offset to (y=1, x=x0, t)
            offset = ((t-1)*Nx*Ny + (x0-1)*Ny) * bytes;
            fseek(fid, offset, 'bof');

            % Read contiguous Y×X block
            slab(:,:,t) = fread(fid, [Ny, nX], precision);
        end

    else
        % Safe fallback (non-contiguous columns)

        for t = 1:Nt
            tOffset = (t-1) * Ny * Nx * bytes;

            for xi = 1:nX
                x = xIdx(xi);
                offset = tOffset + (x-1) * Ny * bytes;

                fseek(fid, offset, 'bof');
                slab(:, xi, t) = fread(fid, Ny, precision);
            end
        end
    end
end


% ======================================================================
% Local funciton
% ======================================================================
function bytes = getPrecisionBytes(precision)
% Returns number of bytes per element for a given MATLAB precision string

    switch lower(precision)
        case {'double','int64','uint64'}
            bytes = 8;

        case {'single','int32','uint32'}
            bytes = 4;

        case {'int16','uint16'}
            bytes = 2;

        case {'int8','uint8','logical'}
            bytes = 1;

        otherwise
            error('Unsupported precision "%s".', precision);
    end

end
