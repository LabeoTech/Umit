function [datMap, matMap] = resetMemmapHandles(Map,varargin)
% This function replaces the handles of memmapfile (Map) and the
% corresponding matfile (metaMap) with the original Writable property.
% This is used by some of the GUI apps to manage RAM usage.

if nargin == 2
    metaMap = varargin{:};
    writable_meta = metaMap.Properties.Writable;
    % Erase matfile first:
    metaMap = [];
end
% Get original writable maps
writable_map = Map.writable;

% Recreate maps:
[datMap, matMap] = mapDatFile(Map.Filename);
datMap.Writable = writable_map;
matMap.Properties.Writable = writable_meta;
end