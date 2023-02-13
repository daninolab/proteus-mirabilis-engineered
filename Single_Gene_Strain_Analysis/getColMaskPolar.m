function [colmask, colregion, colarea_percent] = getColMaskPolar(polarim)

    % Given a polar image, get the colony region & its mask using
    % morphological operations.
    
    startim = im2double(polarim);
    tempim = startim;
    rim_pix = startim<0.7;
    rim_pix(1:400, :) = false;
    bg_val = 0.885;
    tempim(rim_pix) = bg_val; 
%     tempim(1:200, :) = startim(1:200, :);
    % fill in white space with 'agar'
    
%     plate_bg = startim==1;
    plate_bg = startim>0.99;
    tempim(plate_bg) = bg_val;

    % First filter a little
    tempim2 = medfilt2(tempim);
    h = fspecial('average', [1 6]);
    tempim3 = imfilter(tempim2, h);

    % Now let's try using entropy filtering 
    contrast_im2 = adapthisteq(tempim3);
    entropim = entropyfilt(contrast_im2);
    test = contrast_im2; 
    t = 0.885; % set threshold
    test(tempim2>t)=1; 
    test = 1-test; % to get the negative, so putative bacteria pixels are in white
    entropim2 = rescale(entropim);

    % Try emphasizing the correct areas a little extra
    diffim = imabsdiff(entropim, tempim);
    diffim(tempim==1) = 0;

    combine_im = rescale(1*entropim2+1.5*test);
    morphim = imbinarize(rescale(combine_im+diffim));

    % Dilate image a little
    se = strel("disk", 6, 4);
    morphim = imdilate(morphim, se);


    % Use morphological operations for cleanup/filling
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
    colmask(1:50, :) = true;
    % Finish filling in using mean traj
    coltraj = mean(colmask, 2);
    
    for i = 51:100
        if coltraj(i) > 0.5
            % greater than half of the pixels are filled in--fill them all
            % in
            colmask(i, :) = true;
        end
    end
    for i = 400:length(coltraj)
        if coltraj(i)<0.1
            colmask(i, :) = false;
        end
    end
    
    colregion = startim;
    colregion(colmask==0) = 1;

    % Make a petri dish mask: keep petri area to white, and areas where
    % polarim is white, set to black in the mask
    petri_mask = ones(size(polarim));
    petri_mask(polarim==1) = 0;
    
    % Calculate percent of dish
    colarea_percent = length(find(colmask))/length(find(petri_mask));


end
