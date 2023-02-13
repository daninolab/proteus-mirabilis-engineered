% 101822: Processing Metal Images

% This script is like the original flattening script but modified to: Load
% each image, split in half (draw a line?), create two new
% images each with a mirror reflection of one half concatenated, then save
% in new folder; for our metal sensor images

%% Define folders, file names, etc
% define the folder to iterate over
clear all;
clc;

imfolder = input('Write path to image directory to analyze: ', 's');
%'/Users/marianshaw/Dropbox/_Shared_Patterning/Patterning_Expts/Exploratory_Expts/Summer2022/08-17-22_GFP_cheW_HudsonRiver_spotting_allThroughout/spots'
exptdate = input('Write expt date as mm-dd-yy: ', 's');
exptname = input('Write name of expt (eg lrp1): ', 's');
% make a subfolder in it of the new images
oldFolder = cd(imfolder);
subdir_halves = 'Half_Images';
subdir_polar = 'Polar_Half_Images';
if ~isfolder(subdir_halves); mkdir(subdir_halves); end
if ~isfolder(subdir_polar); mkdir(subdir_polar); end

% % get the list of images in the folder
% fileList = dir();
% fileNames = {fileList.name}; % list of file names
% tifNames = {};
% for i = 1:length(fileNames)
%     if contains(fileNames{i}, '.tif')
%         tifNames{end+1} = fileNames{i};
%     end
% end

tifNames = convertFiles(imfolder);

centers = cell(length(tifNames), 1);
%% Define polar analysis struct

polar_analysis = struct('filename', {}, 'colcenter', {}, ...
        'colrad', {}, 'petrimask', {}, 'dist_vec', {}, ...
        'avgd', {}, 'cvs', {}, 'stdevs', {}, 'half', {});

%% iterate over the images
dpcm = 400/2.54; %~158 pix per cm
interp_val = 1000;
thresh = 0.7;
for imnum = 1:length(tifNames)
    tempimname = tifNames{imnum};

    % Check if already done (in case errored out early)
    if length(polar_analysis) >= (2*imnum-1) && ...
            contains(polar_analysis(2*imnum-1).filename, erase(tempimname, '.tif')) && ...
            ~isempty(polar_analysis(2*imnum-1).dist_vec)
        continue
    end
        
    % Show progress
    fprintf('Image %d of %d\n', imnum, length(tifNames));

    % load each image & display
    tempim = imread(tempimname);

%     % Check if this image has a splotch/shouldn't be used
%     chkval = checkIm(tempim, tempimname);
%     if chkval == 1 % Bad image
%         continue % Move to next image
%     end

    disp('Generating halved images...')
    
    [left_half_concat, right_half_concat, imcenter, grayim, petrimask, thresh] =...
        processSpotIm(tempim, thresh);
    % convert each image to polar
    figure(2);
    title(tempimname);
    imshowpair(left_half_concat, right_half_concat, 'montage');
    drawnow;
    centers{imnum} = imcenter;
%     disp('Press any button to continue')
%     waitforbuttonpress;
    % make each image polar
    disp('Generating left half polar im...')
    % Flatten image & save the flattened image file
    [left_polar_im, dist_vec_left] = flattenColonyInterp(left_half_concat, ...
        imcenter, dpcm, interp_val);
    disp('Generating right half polar im...')
    [right_polar_im, dist_vec_right] = flattenColonyInterp(right_half_concat, ...
        imcenter, dpcm, interp_val);
    imshowpair(left_polar_im, right_polar_im, 'montage');
    pause(0.5);

    % Create the new file names for saving
    tempname = erase(tempimname, '.tif');
    newname_right = strcat('Half_Images/', tempname, '_righthalf.tif');
    newname_right_polar = strcat('Polar_Half_Images/', tempname, '_righthalf_polarim.tif');
    newname_left = strcat('Half_Images/', tempname,  '_lefthalf.tif');
    newname_left_polar = strcat('Polar_Half_Images/', tempname,  '_lefthalf_polarim.tif');

    
    % Save
    disp('Saving images and storing analysis...');
    imwrite(left_half_concat, newname_left);
    imwrite(right_half_concat, newname_right);
    imwrite(left_polar_im, newname_left_polar);
    imwrite(right_polar_im, newname_right_polar);

    % Store all the info etc in the polar analysis struct
    polar_analysis(2*imnum-1).filename = newname_right_polar;
    polar_analysis(2*imnum-1).colcenter = imcenter;
    polar_analysis(2*imnum-1).petrimask = petrimask;
    polar_analysis(2*imnum-1).dist_vec = dist_vec_right;
    polar_analysis(2*imnum-1).half = 'right';
    polar_analysis(2*imnum).filename = newname_left_polar;
    polar_analysis(2*imnum).colcenter = imcenter;
    polar_analysis(2*imnum).petrimask = petrimask;
    polar_analysis(2*imnum).dist_vec = dist_vec_left;
    polar_analysis(2*imnum).half = 'left';

    % Do and store preliminary analysis
    right_polar_im = im2double(right_polar_im); %in case if not double
    avgd_right = mean(right_polar_im, 2);
    stdevs_right = std(right_polar_im, 0, 2);
    cvs_right = stdevs_right./avgd_right;
    %Store
    polar_analysis(2*imnum-1).avgd = avgd_right;
    polar_analysis(2*imnum-1).cvs = cvs_right;
    polar_analysis(2*imnum-1).stdevs = stdevs_right;

    left_polar_im = im2double(left_polar_im); %in case if not double
    avgd_left = mean(left_polar_im, 2);
    stdevs_left = std(left_polar_im, 0, 2);
    cvs_left = stdevs_left./avgd_left;
    %Store
    polar_analysis(2*imnum).avgd = avgd_left;
    polar_analysis(2*imnum).cvs = cvs_left;
    polar_analysis(2*imnum).stdevs = stdevs_left;

end % move to next image


% Save the struct of analysis
disp('Done, saving analysis...');
    
% Save the polar analysis file to the current folder
analysis_file_name = strcat(exptdate, '_', exptname, '_radial_analysis.mat');
cd(subdir_polar);
save(analysis_file_name,'polar_analysis');
cd(oldFolder);



%% functions %%%%%%%%%%%%%%%%%%%%%%%%


function chkval = checkIm(tempim, tempfile)
    chkval = 0; % default is to use the image
    imshow(tempim);
    title(tempfile, 'Interpreter', 'none');
    prompt = 'Type n or no if image is not okay to use: ';
    response = lower(input(prompt, 's'));
    if strcmpi(response, 'n') || strcmp(response, 'no')
        chkval = 1; %this will return 1 and skip image
    end
    close(gcf);
end


function [left_half_concat, right_half_concat, imcenter, grayim, petrimask, thresh] =...
    processSpotIm(tempim, thresh)
    
    % To go from the original image with spots on left and not on right to
    % two separate images (each with a mirror half concatenated) of
    % 'spotted' and 'not spotted' (the right half originally)
    figure(2);
    imshow(tempim);
    % get the rotation angle
    h = imdistline;
    disp('Move the line so one end is on the blue kan mark and one is on the inoculum center, then double click.');
    waitfordoubleclick;
    angle = getAngleFromHorizontal(h);
    disp(angle);
    
    % rotate the image at the appropriate angle
    if angle>90
        rotation_angle = -(angle-90);
    else
        rotation_angle = (90-angle);
    end
    
    % get the grayim
    grayim = im2double(rgb2gray(tempim));
    
    imshow(grayim);
    % get the center
    [tempcenter, thresh, petrimask] = findinoc4_figure(tempim, thresh);
    % remove the rim
    grayim(~petrimask) = 1;
    
    % create new image with the colony centered horizontally
    center_x = round(tempcenter(1));
    % get the half of number of columns in the image
    current_col_num = size(grayim, 2);
    current_half_point = round(current_col_num/2);
    % add the needed number of columns to center the image
    if center_x < current_half_point
        % will add the columns to the left
        num_cols_to_add = current_half_point-center_x;
        new_num_cols = current_col_num + num_cols_to_add;
        centered_im = ones(size(grayim, 1), new_num_cols);
        centered_im(:, (new_num_cols-current_col_num+1):end) = grayim;
        figure(2);
        imshow(centered_im);
    else
        % will add the columns to the right
        num_cols_to_add = center_x-current_half_point;
        new_num_cols = current_col_num + num_cols_to_add;
        centered_im = ones(size(grayim, 1), new_num_cols);
        centered_im(:, 1:current_col_num) = grayim;
        imshow(centered_im);
    end
    pause(1);
    % rotate the centered_im
    rotated_im = imrotate(centered_im,rotation_angle,"nearest", "crop");
    rotated_im(rotated_im==0)=1;
    figure(2);
    imshow(rotated_im);
    pause(1);
    % get the halves of the image
    left_half = rotated_im(:, 1:round(size(rotated_im, 2)/2));
    left_half_concat = [left_half, flip(left_half, 2)];
    right_half = rotated_im(:, round(size(rotated_im, 2)/2):end);
    right_half_concat = [flip(right_half, 2), right_half];

    % update the new center's x coordinate
    new_center_x = round(new_num_cols/2);
    imcenter = round(tempcenter);
    imcenter(1) = new_center_x;
end