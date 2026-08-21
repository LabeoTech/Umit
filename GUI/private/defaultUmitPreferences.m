function preferences = defaultUmitPreferences()
%DEFAULTUMITPREFERENCES Return the current default preference structure.

dataViewerDefaults = DataViewerPreferences.defaults();
theme = dataViewerDefaults.theme;
dataViewerDefaults = rmfield(dataViewerDefaults, 'theme');

preferences = struct( ...
    'schemaVersion', 1, ...
    'theme', theme, ...
    'dataViewer', dataViewerDefaults);

end
