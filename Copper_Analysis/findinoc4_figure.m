function [center, threshold, petrimask] = findinoc4_figure(image, threshold)

% Using thresholding, find the inoculum of the colony.
% Function: [center, threshold, petrimask] = findinoc4(image, threshold)
% If no threshold is given, automatically uses 0.69 as threshold
% This one should close the figure unlike findinoc4

    %Find the inoculum of a given mirabilis colony image
    %Assume the input is a color image
    if  isa(image,'char')
        %if passed a file name, use imread
        file_name = image;
        I = imread(file_name); 
    else
        %if passed an image, skip imread
        I = image;
    end
    
    %use given threshold if input
    if ~exist('threshold', 'var')
        threshold = 0.69;
    end
    
    
    %mask out petri dish rim
    if length(size(I))==3
        I = rgb2gray(im2double(I));
    end
    
    [petri_no_rim, petrimask] = removeRimAuto(I, false, false);
    
    %Mask out rest of colony except center inoculum
    %Check if default thresholding works
    check_val = 0; %filter until okay with the edges
    fig = figure('Position', [500 200 900 600]);
    while check_val == 0
        temp = petri_no_rim; temp(temp>threshold)=1; bw_temp = imbinarize(temp); 
        bw_temp = ~bw_temp; 
        figure(fig);
        imshow(bw_temp);
        try
            response = input('Type Yes or y if inoculum okay, No if not okay: ','s');
            if strcmpi(response, 'yes') || strcmpi(response, 'y')
                check_val = 1;
            else
                disp('Threshold current: ');
                disp(threshold);
                threshold = input('Type new threshold val: ');
                check_val = 0;
%                 close gcf
            end
        catch
            disp('Something went wrong. Please try again.');
            continue
        end

    end
    
    %having chosen threshold, use it to find the inoculum
    inoculum = petri_no_rim;
    inoculum(petri_no_rim>threshold)=1; %0.7 somehow keeps the center inoculum woo
    
    %Find center of center inoculum
    bw_Itemp = imbinarize(inoculum); bw_Itemp = ~bw_Itemp;
    figure(fig);
%     set(fig, 'Visible', 'off')
    imshow(bw_Itemp);
    stats = regionprops('table',bw_Itemp,'Centroid',...
        'MajorAxisLength','MinorAxisLength');
    centers = stats.Centroid;
    diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
    radii = diameters/2;
    ind = radii>10; %get the inoculum
    inoc_center = [];
    
    if ~isempty(centers)
        inoc_center = [centers(ind, 1) centers(ind, 2)];
        center = [];
        % If we got more than one, choose the correct circle
        if ~isempty(inoc_center)
            if all(size(inoc_center) == [1, 2])
                center = inoc_center; 
                figure(fig);
                close(gcf);
            else %in case we got multiple circles
                set(fig, 'Visible', 'on');
                figure(fig);
%                 figure('Position', [500 200 900 600]); 
%                 figure('visible', 'on');
                imshow(petri_no_rim);
                radii_temp = 100*ones(1, length(inoc_center));
%                 all_circs = viscircles(inoc_center, radii_temp, 'Color', 'r');
                viscircles(inoc_center, radii_temp, 'Color', 'r');
                for k = 1:length(inoc_center)
                    center_curr = inoc_center(k, :);
                    imcircle = viscircles(center_curr, 100, 'Color', 'g');
                    prompt = 'Use this circle? Type yes or y: ';
                    response = input(prompt, 's');
                    if strcmpi(response, 'yes') || strcmpi(response, 'y')
                        center = center_curr;
                        delete(imcircle);
                        break
                    else
                        delete(imcircle);
                        continue;
                    end
%                     delete(imcircle);
                end
                if isempty(center)
                    % get the center manually
                    center = getCenter;
                end
                clear all_circs radii_temp center_curr k
                clear prompt response
                figure(fig);
                close(gcf);
%                 close(gcf);
%                 close(fig);
            end
        else
            % Get center manually
            center = getCenter;
        end
    else
        % Get center manually
        center = getCenter;
    end
    
    
%     close all;
end


function center = getCenter
    disp('Could not automatically detect inoculum. Draw a circle around the inoculum.');
    h = drawcircle(); 
    waitfordoubleclick;
    center = get(h, 'Center');
    delete(h);
end