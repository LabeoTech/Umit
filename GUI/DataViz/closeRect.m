function closeRect(src,evnt)
% Callback for drawrectangle function.
% It deletes the rectangle object if the User double-clicks
% inside the rectangle.
click = evnt.SelectionType;
selFace = evnt.SelectedPart;
if strcmp(click, 'double') && strcmp(selFace, 'face')
    delete(src);
end

end