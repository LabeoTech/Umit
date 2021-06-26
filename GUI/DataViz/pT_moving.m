function pT_moving(src,evnt)
    XY_str = [num2str(round(evnt.CurrentPosition(1),2)) '; ' num2str(round(evnt.CurrentPosition(2),2))];
    src.Label = XY_str;
    


end