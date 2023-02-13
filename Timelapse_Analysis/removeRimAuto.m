function [maskim, mask] = removeRimAuto(grayim, use_imfindcircles, threshold_on)
    disp('Starting auto rim removal function...');
    
    % Set defaults
    if ~exist('use_imfindcircles', 'var')
        use_imfindcircles = true;
    end
    if ~exist('threshold_on', 'var')
        threshold_on = false;
    end
    
    % First get image to correct format
    if  isa(grayim,'char')
        %if passed a file name, use imread
        file_name = grayim;
        grayim = imread(file_name); 
    end

    % if passed an RGB tiff
    if length(size(grayim))==3
        grayim = rgb2gray(im2double(grayim));
    end
        
    tempim = grayim;
    %Estimate threshold plate radii by taking half the width of the image
    radbound = round(length(grayim)/2);
    radbounds = [radbound-200, radbound+200];
    
    % Try doing the rim finding
    center = [];
    while isempty(center)
        % If desired, threshold image before proceeding (may be needed for cut
        % off images)
        if threshold_on
            tempim = grayim<0.7; %binarizes im and turns rim white
        end

        bigholes = convertIm(tempim);

        if use_imfindcircles
            disp('Trying imfindcircles');
            %Try  imfindcircles to get the border of the agar
            [auto_centers, auto_radii, ~] = ...
                    imfindcircles(bigholes,radbounds, 'ObjectPolarity', 'bright',...
                    'Sensitivity', 0.99);
            if ~isempty(auto_centers)
                center = auto_centers(1, :);
                radius = auto_radii(1);
            end
        end
        if ~use_imfindcircles | isempty(center)
            disp('imfindcircles did not work, trying regionprops instead');
            % If auto_centers is empty, try regionprops instead
            stats = regionprops(bigholes, 'Centroid',...
            'MajorAxisLength','MinorAxisLength', 'Eccentricity');
            if ~(isempty(stats) || isempty(stats.Centroid))
                radius = mean([stats.MajorAxisLength, stats.MinorAxisLength], 2)/2;
                center = stats.Centroid;
            end
        end
        if isempty(center) & ~threshold_on %if not already thresholded
            prompt = 'No center found. Try thresholding first? y/n: ';
            response = lower(input(prompt, 's'));
            if strcmp(response, 'y') || strcmp(response, 'yes') 
                threshold_on = true;
                continue; %go right back to the start and try thresholding
            else
                break; %stop trying to find center/rim--end while loop
            end
        end
        % If the center still is not found
        if isempty(center)
            break
        end
    end

    if isempty(center)
        disp('No center found.');
        maskim = [];
        mask = [];
    else
        disp('Masking image...');
        %Generate mask & mask out image

        maskim = grayim;
        fig = figure('Visible', 'off');
        imshow(maskim);
        h = images.roi.Circle(gca, 'Center', center, 'Radius', radius);
        mask = createMask(h);
        maskim(mask == 0) = 1;
        close(fig);
    end

end


function bigholes = convertIm(grayim)

    disp('Getting rim region...');
    % Go through steps to get inner edge of petri dish/full rim
    % Get image with edges
    edgeim = edge(grayim);
    % Dilate image to connect edges
    se = strel('disk',1);
    dilatedim = imdilate(edgeim,se);

    % Fill in holes
    filled = imfill(dilatedim, 'holes');
    % Determine which pixels were 'holes'
    holes = filled & ~dilatedim;

    % Remove all tiny holes, leaving only the agar as a 'hole'
    bigholes = bwareaopen(holes, 500000);
end