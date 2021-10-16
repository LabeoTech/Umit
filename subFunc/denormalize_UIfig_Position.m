function absPos = denormalize_UIfig_Position(normPos)
% This function creates the absolute position of a UIFigure (normalized
% units not supported in Version R2019a) based on a normalized position and 
% the current screen size.
% Input 
% normPos(1x4 numerical array) = desired normalized figure position;
% Output
% absPos (1x4 numerical array) = absolute position in pixel units

% Get screen size
sz = get( 0, 'ScreenSize');
% calculate absolute position in pixels based on screen size:
absPos = normPos;
absPos([1 3]) = normPos([1 3]).*sz(3);
absPos([2 4]) = normPos([2 4]).*sz(4);
end