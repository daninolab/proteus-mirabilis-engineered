function [x_preproc, colmask, colregion, col_area_percent] = ...
    getGroundTruthColRegion(polarim)

    % First preprocess the polar image--remove leftover rim pixels, add
    % 'agar' background to white space
    tempim = im2double(polarim);
    rim_pix = polarim<0.7;
    rim_pix(1:400, :) = false;
    tempim(rim_pix) = 1; 
    % fill in white space with 'agar'
    bg_val = 0.885;
    plate_bg = polarim==1;
    
    tempim(plate_bg) = bg_val;
    
    % NEW AS OF 07/25/21: Adding a broad low-pass filter over the image
    h = fspecial('gaussian', 50, 0.5);
    tempim = imfilter(tempim, h);
    
    x_preproc = tempim;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Getting the best estimate of the mask for the colony region
    % We will first try an active contours approach, use it in combination
    % with an entropy filter/contrast based approach, then compare the
    % result with the active contours result alone and choose which is the
    % preferred mask.
    
    % First filter a little
    tempim2 = medfilt2(tempim);
    h = fspecial('average', [1 6]);
    tempim3 = imfilter(tempim2, h);
    
    % Preparing for active contours approach
    tempim4 = imadjust(tempim3, [0.8 1]);
    % Start a mask after which do regiongrowing
    tempim5 = imbinarize(tempim4, 0.2);
    X = imadjust(tempim4);
    disp('Creating gabor features...');
    gaborX = createGaborFeatures(X);
    disp('Gabor feature creation complete');

    % Threshold image - manual threshold
    BW = X > 3.058800e-01;

    % Active contour with texture
    iterations = 250;
    disp('beginning active contours...');
    BW = activecontour(gaborX, BW, iterations, 'Chan-Vese');
    disp('active contour complete.');
    % Invert mask
    BW = imcomplement(BW);

    % Create masked image.
    maskedImage = X;
    maskedImage(~BW) = 0;

    % Now let's try alternate approach using entropy filtering 
    contrast_im2 = adapthisteq(tempim3);
    entropim = entropyfilt(contrast_im2);
    test = contrast_im2; 
    t = 0.885; % set threshold
    test(tempim2>t)=1; 
    test = 1-test; % to get the negative, so putative bacteria pixels are in white
    entropim2 = rescale(entropim);
    test2 = rescale(test);
    entropim2 = rescale(entropim);
    test2 = rescale(test);

    % Try emphasizing the correct areas a little extra
    test3 = imadd(maskedImage, test2);

    diffim = imabsdiff(entropim, polarim);
    diffim(polarim==1) = 0;

    combine_im = rescale(1*entropim2+1.5*test3);
    morphim = imbinarize(rescale(combine_im+diffim));

    % Compare morphim & activecontour mask to decide what to do next
    morphim_area = sum(sum(morphim));
    BW_area = sum(sum(BW));
    if morphim_area>BW_area
        % Check if BW_area is better
        imshowpair(maskedImage, morphim, 'montage');
        prompttext = 'Type 1 to choose lefthand (active contours mask) and nothing to choose right hand (entropy based) mask';
        resp = input(prompttext);
        if resp==1
            morphim = BW;
        end
    elseif abs(morphim_area-BW_area)<500000
        % the two images are roughly close/accurate
        % start by a dilation to connect edges etc
        se = strel("disk", 6, 4);
        morphim = imdilate(morphim, se);
    end

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
    
    colregion = contrast_im2;
    colregion(colmask==0) = 1;

    
    plate_mask = ~(polarim==1);
    col_area_percent = length(find(colmask))/length(find(plate_mask));
%     edgemask = edge(colmask);


end

function gaborFeatures = createGaborFeatures(im)

    if size(im,3) == 3
        im = prepLab(im);
    end

    im = im2single(im);

    imageSize = size(im);
    numRows = imageSize(1);
    numCols = imageSize(2);

    wavelengthMin = 4/sqrt(2);
    wavelengthMax = hypot(numRows,numCols);
    n = floor(log2(wavelengthMax/wavelengthMin));
    wavelength = 2.^(0:(n-2)) * wavelengthMin;

    deltaTheta = 45;
    orientation = 0:deltaTheta:(180-deltaTheta);

    g = gabor(wavelength,orientation);
    gabormag = imgaborfilt(im(:,:,1),g);

    for i = 1:length(g)
        sigma = 0.5*g(i).Wavelength;
        K = 3;
        gabormag(:,:,i) = imgaussfilt(gabormag(:,:,i),K*sigma);
    end

    % Increases likelihood that neighboring pixels/subregions are segmented together
    X = 1:numCols;
    Y = 1:numRows;
    [X,Y] = meshgrid(X,Y);
    featureSet = cat(3,gabormag,X);
    featureSet = cat(3,featureSet,Y);
    featureSet = reshape(featureSet,numRows*numCols,[]);

    % Normalize feature set
    featureSet = featureSet - mean(featureSet);
    featureSet = featureSet ./ std(featureSet);

    gaborFeatures = reshape(featureSet,[numRows,numCols,size(featureSet,2)]);

    % Add color/intensity into feature set
    gaborFeatures = cat(3,gaborFeatures,im);

end

function out = prepLab(in)

    % Convert L*a*b* image to range [0,1]
    out = in;
    out(:,:,1) = in(:,:,1) / 100;  % L range is [0 100].
    out(:,:,2) = (in(:,:,2) + 86.1827) / 184.4170;  % a* range is [-86.1827,98.2343].
    out(:,:,3) = (in(:,:,3) + 107.8602) / 202.3382;  % b* range is [-107.8602,94.4780].

end