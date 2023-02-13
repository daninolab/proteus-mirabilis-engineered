function newim = getCentCropIm(tempim, petrimask, colcent)
    % Given a cartesian colony image, the mask for the petri rim, and the
    % colony center, this function will mask out the rim, center the colony
    % in the image, and crop as much as possible in a square around it.
    rimless = im2double(rgb2gray(tempim));
    rimless( petrimask == 0) = 1;

    % crop and recenter around actual center
    % Find bounding box of petri dish
    stats = regionprops(petrimask,'Boundingbox');
    % topleft = [stats.BoundingBox(1), stats.BoundingBox(2)];
    % topleft = round(topleft);
%     xleft = round(stats.BoundingBox(1));
%     ytop = round(stats.BoundingBox(2));
%     xright = xleft+stats.BoundingBox(3);
%     ybot = ytop + stats.BoundingBox(4);

    leftcol = round(stats.BoundingBox(1));
    toprow = round(stats.BoundingBox(2));
    rightcol = leftcol+stats.BoundingBox(3);
    botrow = toprow + stats.BoundingBox(4);
    
        % Now get bounding box centered around col center which includes whole
    % petri dish
    centrow = round(colcent(2)); 
    centcol = round(colcent(1)); % cent col = center column
    wcol = round(max(abs(centcol-leftcol), abs(rightcol-centcol)));
    hrow = round(max(abs(centrow-toprow), abs(botrow-centrow)));
    % Now figure out the square's width/height
    newsqw = max(wcol, hrow); 
    newleftcol = centcol-newsqw;
    newrightcol = centcol+newsqw;
    newtoprow = centrow-newsqw;
    newbotrow = centrow+newsqw;
%     newleftcol = centcol-wcol;
%     newrightcol = centcol+wcol;
%     newtoprow = centrow-hrow;
%     newbotrow = centrow+hrow;
    
    % Check that this bounding box will fit in image, otherwise need to pad
    % image
    if newleftcol<0
        % Need to pad with columns to the left so that newleftcol becomes
        % column 1, then modify newrightcol accordingly
        rimless = padarray(rimless, [0, abs(newleftcol)] , 1, 'pre');
        newrightcol = newrightcol+abs(newleftcol);
        newleftcol = 1;
    end
    if newrightcol>size(rimless, 2)
        % righthand needs padding
        padval = newrightcol-size(rimless, 2);
        rimless = padarray(rimless, [0, padval] , 1, 'post');
    end
    if newtoprow<0
        % need to pad above top row of image
        rimless = padarray(rimless, [abs(newtoprow), 0] , 1, 'pre');
        newbotrow = newbotrow+abs(newtoprow);
        newtoprow = 1;
    end
    if newbotrow>size(rimless, 1)
        % bottom needs padding
        padval = newbotrow-size(rimless, 1);
        rimless = padarray(rimless, [padval, 0] , 1, 'post');
    end

    % Crop image; note that x refers to columns and y refers to rows
%     newim = rimless(newytop:newybot, newxleft:newxright);
    newim = rimless(newtoprow:newbotrow, newleftcol:newrightcol);
    
    
end