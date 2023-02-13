function [colregion, colmask, colarea_percent, edgemask] = ...
    getColonyRegion(rimless, petri_mask)
    % This function will use morphological operations, filtering, etc to
    % extract the colony region (including any outer faint swarm rings) from
    % the input rimless image. It will also give a rough percentage of how much
    % percentage of the agar the colony takes up.

    % Try to get the rim perfect
    [tempcenter, temprad] = imfindcircles(petri_mask, ...
        [round(length(petri_mask)/2)-200, round(length(petri_mask)/2)+200], ...
        'ObjectPolarity','bright', 'Sensitivity', 0.99);
       
    if isempty(tempcenter)
        % If auto_centers is empty, try regionprops instead
        stats = regionprops(petri_mask, 'Centroid',...
        'MajorAxisLength','MinorAxisLength', 'Eccentricity');
        if ~(isempty(stats) || isempty(stats.Centroid))
            temprad = mean([stats.MajorAxisLength, stats.MinorAxisLength], 2)/2;
            tempcenter = stats.Centroid;
        end
    end
    if ~isempty(tempcenter)
        figure;
%         set(gcf, 'Visible', 'off');
        clf('reset');
        imshow(petri_mask);
        % shrink rim a little
        h = drawcircle('Center', tempcenter, 'Radius', temprad-50);
        tempmask2 = ~createMask(h);
        tempmask2(~petri_mask)=false;
        close(gcf);
        % mask the og image
        tempim2 = rimless;
        tempim2(~tempmask2)=1;
        % or what if i just decide tempthresh arbitrarily?
        tempthresh = 0.7;
        rimless(tempim2<tempthresh) = 1;
    end

    
    
    % Use rimless image to get colony
    contrast_im2 = adapthisteq(rimless);
    entropim = entropyfilt(contrast_im2);
    test = contrast_im2; 
    t = 0.885; % set threshold
    test(test>t)=1; 
    test = 1-test; % to get the negative, so putative bacteria pixels are in white

    % Try getting the intersection between that and the entropy image
    entropim2 = rescale(entropim);
    test2 = rescale(test);
    diffim = imabsdiff(entropim, rimless);
    diffim(rimless==1) = 0;

    combine_im = rescale(1*entropim2+1.5*test2);

    morphim = imbinarize(rescale(combine_im+diffim));

    se = strel("disk", 3, 4);
    morphim1 = imopen(morphim, se);
    se = strel("disk", 3, 8); 
    morphim2 = imdilate(morphim1, se);
    morphim3 = bwareaopen(morphim2, 400, 8);
    morphim4 = imfill(morphim3, "holes");
    morphim5 = bwareaopen(morphim4, 1000, 8);
    se = strel('disk', 10, 8); 
    morphim6 = imclose(morphim5, se);
    morphim7 = imfill(morphim6, 'holes'); 
    % make mask using largest area thing
    colmask = bwareafilt(morphim7, 1);

    colregion = contrast_im2;
    colregion(colmask==0) = 1;
    colarea_percent = length(find(colmask))/length(find(petri_mask));

    edgemask = edge(colmask);



end