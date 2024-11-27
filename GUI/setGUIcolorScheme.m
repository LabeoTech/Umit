function setGUIcolorScheme(parentHandle,Mode)
% Sets the color scheme of a GUIelement and its children.
% Used to set light/dark modes in Umit GUIs
assert(ishandle(parentHandle),'Invalid Input! First argument must be a valid Graphical Element handle')
assert(ismember(lower(Mode),{'light','dark'}),'Invalid category! Second argument must be "Light" or "Dark".')
% Load Color scheme:
colorScheme = jsondecode(fileread(fullfile(fileparts(mfilename('fullpath')),'GUIcolorScheme.json')));
% Get list of all GUI element handles and respective class names:
allHandles = findobj(parentHandle);
allClasses = arrayfun(@(x) class(x),allHandles,'UniformOutput',false);
% Loop across all elements and update colors:
Mode = lower(Mode);
fieldMap = containers.Map({'light','dark'},{'LightMode','DarkMode'});
iconMap = containers.Map({'light','dark'},{'White.png','Black.png'});
iconRepMap = containers.Map({'White.png','Black.png'},{'Black.png','White.png'});
icoFields = {'Icon','ImageSource'}; % Fields containing Icons.
for ii = 1:length(allHandles)
    % Special case:
    % For buttons with icons OR Images used as buttons, look for light/dark versions of existing
    % icons:
    icoIdx = ismember(icoFields,properties(allHandles(ii)));
    if any(icoIdx)
        iconName = allHandles(ii).(icoFields{icoIdx});
        if ~iconName;continue;end
        newIcon = replace(iconName,iconMap(Mode),iconRepMap(iconMap(Mode)));
        if exist(newIcon,'file')
            allHandles(ii).(icoFields{icoIdx}) = newIcon;
        end
    end
    % Check if element is listed in color scheme
    idx = cellfun(@(x) any(strcmpi(x,allClasses{ii})),{colorScheme.ElementID});
    if ~any(idx)
        continue
    end
    % Get field information:
    fn = fieldnames(colorScheme(idx).LightMode);
    bHasField = ismember(regexprep(fn, '_.*$', ''),properties(allHandles(ii)));
    fn(~bHasField) = [];
    if isempty(fn);continue;end
    % Set fields:
    for jj = 1:length(fn)
        value = colorScheme(idx).(fieldMap(Mode)).(fn{jj});
        if isnumeric(value)
            value = repelem(value,1,3);
        end
        % Handle subfield (e.g. ax.Title.Color):
        if contains(fn{jj},'_')
            str = strsplit(fn{jj},'_');
            eval(['set(allHandles(ii).' str{1} ',"' str{2} '",value);'])
        else
            set(allHandles(ii),fn{jj},value);
        end
    end
end
end