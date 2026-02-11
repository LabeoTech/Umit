function nBytes = getByteSize(precision)
% GETBYTESIZE Number of bytes per element for a given precision.

switch precision
    case {'single','uint32','int32'}
        nBytes = 4;
    case {'double','uint64','int64'}
        nBytes = 8;
    case {'uint16','int16'}
        nBytes = 2;
    case {'uint8','int8','char'}
        nBytes = 1;
    otherwise
        error('getByteSize:UnsupportedPrecision', ...
              'Unsupported precision: %s', precision);
end
end
