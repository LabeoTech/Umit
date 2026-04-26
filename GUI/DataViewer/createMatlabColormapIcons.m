function createMatlabColormapIcons(iconFolder)
%CREATEMATLABCOLORMAPICONS Create PNG icons for MATLAB built-in colormaps.
%
%   createMatlabColormapIcons(iconFolder)
%
%   Creates horizontal gradient PNG files for MATLAB colormaps used in the
%   DataViewer colormap dropdown.

if nargin < 1 || isempty(iconFolder)
    iconFolder = fullfile(pwd, 'resources', 'colormaps', 'icons');
end

if ~isfolder(iconFolder)
    mkdir(iconFolder);
end

mapNames = {'parula', 'turbo', 'jet', 'hsv', 'hot', 'gray'};
nColors = 64;
iconHeight = 18;

for iMap = 1:numel(mapNames)
    mapName = mapNames{iMap};

    try
        cmap = feval(mapName, nColors);
    catch ME
        warning('createMatlabColormapIcons:UnavailableMap', ...
            'Could not create icon for "%s": %s', mapName, ME.message);
        continue
    end

    img = repmat(reshape(cmap, [1, nColors, 3]), iconHeight, 1, 1);

    outFile = fullfile(iconFolder, [mapName, '.png']);
    imwrite(img, outFile);
end

end