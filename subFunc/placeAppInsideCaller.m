function placeAppInsideCaller(callerApp, mainApp, positionName, varargin)
%PLACEAPPINSIDECALLER Place one app/window relative to another app/window.
%
%   placeAppInsideCaller(callerApp, mainApp, positionName)
%   placeAppInsideCaller(..., 'margin', 20)
%
%   Positions mainApp relative to callerApp using figure Position
%   coordinates [left bottom width height].
%
%   Inputs:
%       callerApp    - App Designer app object with UIFigure property,
%                      figure/uifigure handle, or object with Position.
%
%       mainApp      - App Designer app object with UIFigure property,
%                      figure/uifigure handle, or object with Position.
%
%       positionName - Placement string:
%                         "center"
%                         "top-left"
%                         "top-right"
%                         "bottom-left"
%                         "bottom-right"
%
%   Name-Value options:
%       margin       - Pixel margin for corner positions. Default: 20.
%
%   Rules:
%       - "center" overlays the center of mainApp over the center of callerApp.
%       - "top-*" aligns top coordinates.
%       - "bottom-*" aligns bottom coordinates.
%       - "left" and "right" use the requested margin.
%       - mainApp size is preserved.

p = inputParser;
p.FunctionName = 'placeAppInsideCaller';

addRequired(p, 'callerApp');
addRequired(p, 'mainApp');
addRequired(p, 'positionName', ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));

addParameter(p, 'margin', 20, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);

parse(p, callerApp, mainApp, positionName, varargin{:});

positionName = lower(strtrim(char(string(p.Results.positionName))));
margin = double(p.Results.margin);

validPositions = { ...
    'center', ...
    'top-left', ...
    'top-right', ...
    'bottom-left', ...
    'bottom-right'};

if ~ismember(positionName, validPositions)
    error('placeAppInsideCaller:InvalidPosition', ...
        'positionName must be one of: %s.', strjoin(validPositions, ', '));
end

callerFig = iGetFigureHandle(callerApp, 'callerApp');
mainFig   = iGetFigureHandle(mainApp, 'mainApp');

drawnow limitrate

callerUnits = callerFig.Units;
mainUnits   = mainFig.Units;

cleanupObj = onCleanup(@() iRestoreUnits(callerFig, callerUnits, mainFig, mainUnits)); %#ok<NASGU>

callerFig.Units = 'pixels';
mainFig.Units = 'pixels';

callerPos = double(callerFig.Position);
mainPos   = double(mainFig.Position);

callerLeft   = callerPos(1);
callerBottom = callerPos(2);
callerWidth  = callerPos(3);
callerHeight = callerPos(4);

mainWidth  = mainPos(3);
mainHeight = mainPos(4);

callerRight = callerLeft + callerWidth;
callerTop   = callerBottom + callerHeight;

switch positionName
    case 'center'
        newLeft = callerLeft + (callerWidth - mainWidth) / 2;
        newBottom = callerBottom + (callerHeight - mainHeight) / 2;

    case 'top-left'
        newLeft = callerLeft + margin;
        newBottom = callerTop - mainHeight;

    case 'top-right'
        newLeft = callerRight - margin - mainWidth;
        newBottom = callerTop - mainHeight;

    case 'bottom-left'
        newLeft = callerLeft + margin;
        newBottom = callerBottom;

    case 'bottom-right'
        newLeft = callerRight - margin - mainWidth;
        newBottom = callerBottom;
end

mainFig.Position = [newLeft, newBottom, mainWidth, mainHeight];

drawnow limitrate

end

% =========================================================================
% Local helpers
% =========================================================================
function fig = iGetFigureHandle(appOrFig, inputName)
%IGETFIGUREHANDLE Resolve app object or figure/uifigure handle.

if isa(appOrFig, 'matlab.ui.Figure')
    fig = appOrFig;
    return
end

if isobject(appOrFig) && isprop(appOrFig, 'UIFigure') && ...
        isa(appOrFig.UIFigure, 'matlab.ui.Figure')
    fig = appOrFig.UIFigure;
    return
end

if isobject(appOrFig) && isprop(appOrFig, 'Position')
    fig = appOrFig;
    return
end

error('placeAppInsideCaller:InvalidInput', ...
    '%s must be an app object with UIFigure or a figure/uifigure handle.', ...
    inputName);

end

function iRestoreUnits(callerFig, callerUnits, mainFig, mainUnits)
%IRESTOREUNITS Restore original figure units.

try
    if isvalid(callerFig)
        callerFig.Units = callerUnits;
    end
catch
end

try
    if isvalid(mainFig)
        mainFig.Units = mainUnits;
    end
catch
end

end