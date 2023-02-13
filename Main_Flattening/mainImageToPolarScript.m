%%%%%%%%% Proteus mirabilis Colony Image to Polar%%%%%%%%%

% Developed for "Engineered bacterial swarm patterns as spatial records of
% environmental inputs", Nature Chemical Biology, 2023, by Anjali Doshi &
% Marian Shaw, in the lab of Tal Danino at Columbia University. 

% Given path to expt images of P. mirabilis scans, will flatten (ie will 
% convert them from Cartesian to polar, based on finding the center of the 
% colony, such that colony rings become rows in the new image), then will
% save polarims in a new subdirectory, will get avgs/cvs/col radii etc into a 
% structure and save the struct in a new analysis_polar.mat file

%% Get Image Directory & Info

% Set folder path and define experiment info
origFolder = input('Write path to image directory to analyze: ', 's');
exptdate = input('Write expt date as mm-dd-yy: ', 's');
exptname = input('Write name of expt (eg lrp1): ', 's');

% within folder, convert all image files to uint8 TIFF images w/ 3 channels
tifNames = convertFiles(origFolder);

% Set constants
interp_val = 1000;
resolution = input('Write resolution in dpi (eg 400): ');
dpcm = resolution/2.54;
% dpcm = 400/2.54; %~158 pix per cm

%% Set Up Struct for Analysis

tempname = strcat(exptdate, '_', exptname, '_radial_analysis.mat');

% Make a directory to store the polar images
subdir = strcat(exptdate, '_', exptname, '_polarims');
if ~isfolder(subdir)
    mkdir(subdir);
end
oldFolder = cd(subdir);

% If analysis already done, move to next expt
if isfile(tempname)
    disp('Already exists');
else
    % Make a structure to hold the interpolants/masks/etc
    polar_analysis = struct('filename', {}, 'colcenter', {}, ...
        'colrad', {}, 'petrimask', {}, 'dist_vec', {}, ...
        'avgd', {}, 'cvs', {}, 'stdevs', {});
end
    

%% Main Analysis Loop
thresh = 0.7; % default
% Iterate over each image to produce/save polarim + analysis
for imnum = 1:length(tifNames)
    % Get the current image & center
    tempfile = tifNames{imnum};
    tempim = imread(strcat('../', tempfile)); % currently in subdir
    tempname = erase(tempfile, '.tif');
    newname = strcat(tempname, '_polarim.tif');

    % Check if already done (in case errored out early)
    if length(polar_analysis) >= imnum && ...
            strcmp(polar_analysis(imnum).filename, tempfile) && ...
            ~isempty(polar_analysis(imnum).dist_vec)
        continue
    end
        
    % Show progress
    fprintf('Image %d of %d\n', imnum, length(tifNames));

    % Check if this image has a splotch/shouldn't be used
    chkval = checkIm(tempim, tempfile);
    if chkval == 1 % Bad image
        continue % Move to next image
    end
        
    %Find the inoculum center; use the latest threshold from each previous
    %loop
    [tempcenter, thresh, petrimask] = findinoc4(tempim, thresh);
%     tempcenter = findinoc3(tempim);
    polar_analysis(imnum).colcenter = tempcenter;
    % Get the radius
    polar_analysis(imnum).colrad = getColRadManual(tempim, tempcenter);
              
%     % Preprocess Image & get a masked image
%     [maskim, petrimask] = getMask(tempim);
                    
    % Save rimless petri mask
    polar_analysis(imnum).petrimask = petrimask;
    polar_analysis(imnum).filename = tempfile;

    disp('Flattening image...');
    maskim = maskImage(tempim, petrimask);
    % Flatten image & save the flattened image file
    [polar_im, dist_vec] = flattenColonyInterp(maskim, ...
        tempcenter, dpcm, interp_val);
    tempname = erase(tempfile, '.tif');
    newname = strcat(tempname, '_polarim.tif');
        
    % Save
    disp('Saving image and storing analysis...');
    imwrite(polar_im, newname);
        
    % Save analysis in polar_analysis struct

    polar_im = im2double(polar_im); %in case if not double
    avgd = mean(polar_im, 2);
    stdevs = std(polar_im, 0, 2);
    cvs = stdevs./avgd;
    %Store
    polar_analysis(imnum).avgd = avgd;
    polar_analysis(imnum).cvs = cvs;
    polar_analysis(imnum).stdevs = stdevs;
    polar_analysis(imnum).dist_vec = dist_vec;  
    clear tempfile tempim tempcenter chkval maskim petrimask;
    clear polar_im dist_vec tempname newname avgd stdevs cvs;
        
end %End iterating over that expt's images
disp('Done, saving analysis...');
    
% Save the polar analysis file to the current folder
tempname = strcat(exptdate, '_', exptname, '_radial_analysis.mat');
save(tempname,'polar_analysis');
    
%     clear imnum date_curr gene_curr subdir polar_analysis tempname
cd(oldFolder);



%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function chkval = checkIm(tempim, tempfile)
    chkval = 0; % default is to use the image
    figure; imshow(tempim);
    title(tempfile, 'Interpreter', 'none');
    prompt = 'Type n or no if image is not okay to use: ';
    response = lower(input(prompt, 's'));
    if strcmpi(response, 'n') || strcmp(response, 'no')
        chkval = 1; %this will return 1 and skip image
    end
    close(gcf);
end


function [maskim, petrimask] = getMask(tempim)
    
    use_imfindcircles = true;
    threshold_on = false;
    %Get the mask to remove rim
    [maskim, petrimask] = removeRimAuto(tempim, use_imfindcircles, threshold_on);
%     [maskim, petrimask] = removeRimAuto(tempim);
    figure; imshow(maskim); 
    prompt = 'Ok? Type n or no if not: ';
    response = lower(input(prompt, 's'));
    close(gcf);
    if strcmpi(response, 'n') || strcmp(response, 'no') 
        %Find the inoculum center
        [maskim, petrimask] = removerim_manual(tempim); 
    end    
      
%     maskim = getMaskedPlate(tempim, petrimask);
end

function maskim = maskImage(tempim, petrimask, colcenter)
        maskim = rgb2gray(im2double(tempim));
        maskim(petrimask == 0) = 1;
end

