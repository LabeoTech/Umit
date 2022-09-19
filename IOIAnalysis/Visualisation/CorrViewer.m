function CorrViewer(cMat)

currentPixel = [1,1];
h_fig = figure('WindowButtonDownFcn', @ChangePtPos);
h_ax = axes('Parent', h_fig);
imagesc(h_ax, reshape(cMat(1,:), 192, 192), [0 1]);
axis(h_ax, 'off', 'image');
hold(h_ax, 'on');
plot(h_ax, currentPixel(1), currentPixel(2), 'or');
hold(h_ax, 'off');

    function ChangePtPos(Obj, Evnt)
        Pos = round(Obj.Children(1).CurrentPoint);
        
        Pos = Pos(1,1:2);
        if( any(Pos < 1) | any(Pos > 256) )
            return;
        end
        currentPixel = Pos;
                
        Id = floor((currentPixel(1)-1))*192 + floor(currentPixel(2));
        if( Id < 1 )
            Id = 1;
        end
        if( Id > 192*192)
            Id = 192*192;
        end
        imagesc(h_ax, reshape(cMat(:,Id), 192, 192), [-0.25 0.25]);
        axis(h_ax, 'off', 'image');
        hold(h_ax, 'on');
        plot(h_ax, currentPixel(1), currentPixel(2), 'or');
        hold(h_ax, 'off');
    end
end