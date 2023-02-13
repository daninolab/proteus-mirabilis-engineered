function colrad_cm = getColRadManual_Polar(polarim, tempdist)

    % Given an experimental image, get colony radius
    % Display

    fig = figure; imshow(imadjust(polarim));
    w = size(polarim, 2);
    % Ask user to draw circle
    disp('extend a rectangle over the colony (include flares) & double click');
    h = drawrectangle('Position', [0, 0, w, 300],'FaceAlpha', 0.7); %get max colony radius for the image
    waitfordoubleclick;
    % Get radius
    pos = h.Position;
    ht = pos(4);
    close(fig);
    if exist('tempdist', 'var')
        colrad_cm = tempdist(round(ht));
    else
        disp('No dist given')
        colrad_cm = [];
    end
    
end