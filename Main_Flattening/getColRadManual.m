function colrad = getColRadManual(tempim, tempcenter, og_im, fig_size)

    % Given an experimental image, get colony radius
    % Display
    figs = {};
    if exist('og_im', 'var') && ~isempty(og_im)
        figs{1} = figure;
        imshow(og_im);
        figs{2} = figure;
        imshow((tempim));
    else
        figs{1} = figure;
        imshow((tempim));
    end
    if exist('fig_size', 'var')
        set(gcf, 'Position', fig_size);
    end
%     fig = figure; imshow(tempim);
    % Ask user to draw circle
    disp('extend a circle over the colony (include flares) & double click');
    h = drawcircle('Center', tempcenter,'Radius', 100, 'FaceAlpha', 0.7); %get max colony radius for the image
    waitfordoubleclick;
    % Get radius
    colrad = h.Radius;
    for i = 1:length(figs)
        close(figs{i});
    end
%     close(fig);
    
end