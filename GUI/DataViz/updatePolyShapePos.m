function updatePolyShapePos(src,evt,polyPlotHandle)
% Changes the position of polyshape objects in pShapes to translation,
% rotation and scaling of an images.roi.Rectangle.
pShapes= [polyPlotHandle.Shape];
% Positions of the center of the Rectangle:
CurrentRectCenter = [(evt.CurrentPosition(1) + (evt.CurrentPosition(3)/2)) ...
    (evt.CurrentPosition(2) + (evt.CurrentPosition(4)/2))];
PreviousRectCenter = [(evt.PreviousPosition(1) + (evt.PreviousPosition(3)/2))...
    (evt.PreviousPosition(2) + (evt.PreviousPosition(4)/2))];
% Translation
if any(CurrentRectCenter ~= PreviousRectCenter)
    pos_dif = CurrentRectCenter - PreviousRectCenter;
    
    pShapes = arrayfun(@(x) translate(x, pos_dif(1),pos_dif(2)), pShapes);
%     pShapes = translate(pShapes, pos_dif);
end
% Scaling
if evt.CurrentPosition(3) ~= evt.PreviousPosition(3)
    scale_factor = 1 - (evt.PreviousPosition(3) - evt.CurrentPosition(3)) / ...
        evt.PreviousPosition(3);
    pShapes = arrayfun(@(x) scale(x, scale_factor,CurrentRectCenter), pShapes);
%     pShapes = scale(pShapes, scale_factor,CurrentRectCenter);
end
% Rotation
if evt.PreviousRotationAngle ~= evt.CurrentRotationAngle
    ang = evt.PreviousRotationAngle - evt.CurrentRotationAngle;
    pShapes = arrayfun(@(x) rotate(x, ang,CurrentRectCenter), pShapes);
%     pShapes = rotate(pShapes, ang,CurrentRectCenter);
end
src.HandleVisibility = 'off';
% set(polyPlotHandle, 'Shape', pShapes);
arrayfun(@(x,y) set(x, 'Shape',y), polyPlotHandle, pShapes');
src.HandleVisibility = 'on';
end


