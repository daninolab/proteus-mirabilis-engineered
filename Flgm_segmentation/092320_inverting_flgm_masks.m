% Inverting all the unet masks for flgm for fig 4

flgm_imdir = '/Users/anjali/Dropbox/Patterning_Expts_Analysis/Experiments/Single_Gene_Expts/datasets_for_col_segmentation/flgm_subset/unet_fig_testims_masks_preds/';
dir_list = dir(fullfile(flgm_imdir, '*.tif'));
im_list = {dir_list.name};

for i = 1:length(im_list)
    if contains(im_list{i}, 'skel') | contains(im_list{i}, 'uncrop') & ~contains(im_list{i}, 'invert')
        % load the image
        tempim = im2double(imread(fullfile(dir_list(i).folder, dir_list(i).name)));
        
        % convert to black
        tempim2 = ~tempim;
        % save new name
        newname = im_list{i};
        newname = strrep(newname, '.tif', '_invert.tif');
        % save
        imwrite(tempim2, fullfile(flgm_imdir, newname));
    end
end

%% Instead of inverting, let's use the mask to make two versions of each image
% Let's first make the image from polarim A
% Steps: 
% Load the mask
tempmask = im2double(imread(fullfile(dir_list(2).folder, dir_list(2).name)));
tempimname = dir_list(2).name;
tempimname = strrep(tempimname, '_skel.tif', '.tif');
im_ind = find(strcmpi({dir_list.name}, tempimname));
polarim = im2double(imread(fullfile(dir_list(im_ind).folder, dir_list(im_ind).name)));

% imshow(tempmask);
tempmask_bin = imbinarize(tempmask);

% Let's figure out how to crop this to where the regions will be split
tempmodes = zeros(1, size(tempmask, 2));
for i = 1:size(tempmask, 2)
    tempmodes(i) = sum(tempmask(:, i));
end
numbounds = mode(tempmodes);

for i = 1:size(tempmask, 2)
    if sum(tempmask(:, i)) == numbounds
        leftcolumn_val = i;
        break;
    end
end

for i = size(tempmask, 2):-1:1
    if sum(tempmask(:, i)) == numbounds
        rightcolumn_val = i;
        break;
    end
end

tempmask_bin_crop = tempmask_bin(:, leftcolumn_val:rightcolumn_val);
% Dilate with a horizontal line in order to connect boundaries...
se = strel('line', 10, 0);
tempmask2 = imdilate(tempmask_bin_crop, se);
tempmask3 = ~tempmask2;

% Now get the regions:
CC = bwconncomp(tempmask3);
L = labelmatrix(CC);

% Expand/uncrop L
L2 = zeros(size(polarim));
L2(:, leftcolumn_val:rightcolumn_val) = L;
L2(:, 1:(leftcolumn_val-1)) = repmat(L(:, 1), 1, (leftcolumn_val-1));
L2(:, (rightcolumn_val+1):end) = repmat(L(:, end), 1, size(polarim, 2)-rightcolumn_val);

% Now convert into an alternating image
output_mask1 = zeros(size(polarim));
output_mask2 = output_mask1;

output_mask1(L2==2) = 1;
output_mask2(L2==3) = 1;

% Now finally we will save the output image
polarim_cond1 = polarim;
polarim_cond1(output_mask1==0) = 1;
polarim_cond2 = polarim;
polarim_cond2(output_mask2==0) = 1;
imshowpair(polarim_cond1, polarim_cond2, 'montage');


imwrite(polarim_cond1, fullfile(flgm_imdir, strrep(tempimname, '.tif', '_region1.tif')));
imwrite(polarim_cond2, fullfile(flgm_imdir, strrep(tempimname, '.tif', '_region2.tif')));


%% Now let's repeat for test image 2: Condition b

tempimname = '07-20-21_flgM_5i_B_ben24_ben38_inc62_002_polarim_skel.tif';  
im_ind = find(strcmpi({dir_list.name}, tempimname));
tempmask = im2double(imread(fullfile(dir_list(im_ind).folder, dir_list(im_ind).name)));

tempimname = strrep(tempimname, '_skel.tif', '.tif');
im_ind = find(strcmpi({dir_list.name}, tempimname));
polarim = im2double(imread(fullfile(dir_list(im_ind).folder, dir_list(im_ind).name)));

% imshow(tempmask);
tempmask_bin = imbinarize(tempmask);

% Let's figure out how to crop this to where the regions will be split
tempmodes = zeros(1, size(tempmask, 2));
for i = 1:size(tempmask, 2)
    tempmodes(i) = sum(tempmask(:, i));
end
numbounds = mode(tempmodes);

for i = 1:size(tempmask, 2)
    if sum(tempmask(:, i)) == numbounds
        leftcolumn_val = i;
        break;
    end
end

for i = size(tempmask, 2):-1:1
    if sum(tempmask(:, i)) == numbounds
        rightcolumn_val = i;
        break;
    end
end

tempmask_bin_crop = tempmask_bin(:, leftcolumn_val:rightcolumn_val);
% Dilate with a horizontal line in order to connect boundaries...
se = strel('line', 10, 0);
tempmask2 = imdilate(tempmask_bin_crop, se);
tempmask3 = ~tempmask2;

% Now get the regions:
CC = bwconncomp(tempmask3);
L = labelmatrix(CC);

% Expand/uncrop L
L2 = zeros(size(polarim));
L2(:, leftcolumn_val:rightcolumn_val) = L;
L2(:, 1:(leftcolumn_val-1)) = repmat(L(:, 1), 1, (leftcolumn_val-1));
L2(:, (rightcolumn_val+1):end) = repmat(L(:, end), 1, size(polarim, 2)-rightcolumn_val);

% Now convert into an alternating image
output_mask1 = zeros(size(polarim));
output_mask2 = output_mask1;

output_mask1(L2==2) = 1;
output_mask2(L2==3) = 1;

% Now finally we will save the output image
polarim_cond1 = polarim;
polarim_cond1(output_mask1==0) = 1;
polarim_cond2 = polarim;
polarim_cond2(output_mask2==0) = 1;
imshowpair(polarim_cond1, polarim_cond2, 'montage');


imwrite(polarim_cond1, fullfile(flgm_imdir, strrep(tempimname, '.tif', '_region1.tif')));
imwrite(polarim_cond2, fullfile(flgm_imdir, strrep(tempimname, '.tif', '_region2.tif')));

%% Finally repeat for image 3



tempimname = '07-20-21_flgM_5i_C_ben24_inc38_inc62_002_polarim_skel.tif';   
im_ind = find(strcmpi({dir_list.name}, tempimname));
tempmask = im2double(imread(fullfile(dir_list(im_ind).folder, dir_list(im_ind).name)));

tempimname = strrep(tempimname, '_skel.tif', '.tif');
im_ind = find(strcmpi({dir_list.name}, tempimname));
polarim = im2double(imread(fullfile(dir_list(im_ind).folder, dir_list(im_ind).name)));

% imshow(tempmask);
tempmask_bin = imbinarize(tempmask);

% Let's figure out how to crop this to where the regions will be split
tempmodes = zeros(1, size(tempmask, 2));
for i = 1:size(tempmask, 2)
    tempmodes(i) = sum(tempmask(:, i));
end
numbounds = mode(tempmodes);

for i = 1:size(tempmask, 2)
    if sum(tempmask(:, i)) == numbounds
        leftcolumn_val = i;
        break;
    end
end

for i = size(tempmask, 2):-1:1
    if sum(tempmask(:, i)) == numbounds
        rightcolumn_val = i;
        break;
    end
end

tempmask_bin_crop = tempmask_bin(:, leftcolumn_val:rightcolumn_val);
% Dilate with a horizontal line in order to connect boundaries...
se = strel('line', 10, 0);
tempmask2 = imdilate(tempmask_bin_crop, se);
tempmask3 = ~tempmask2;

% Now get the regions:
CC = bwconncomp(tempmask3);
L = labelmatrix(CC);

% Expand/uncrop L
L2 = zeros(size(polarim));
L2(:, leftcolumn_val:rightcolumn_val) = L;
L2(:, 1:(leftcolumn_val-1)) = repmat(L(:, 1), 1, (leftcolumn_val-1));
L2(:, (rightcolumn_val+1):end) = repmat(L(:, end), 1, size(polarim, 2)-rightcolumn_val);

% Now convert into an alternating image
output_mask1 = zeros(size(polarim));
output_mask2 = output_mask1;

output_mask1(L2==2) = 1;
output_mask2(L2==3) = 1;

% Now finally we will save the output image
polarim_cond1 = polarim;
polarim_cond1(output_mask1==0) = 1;
polarim_cond2 = polarim;
polarim_cond2(output_mask2==0) = 1;
imshowpair(polarim_cond1, polarim_cond2, 'montage');


imwrite(polarim_cond1, fullfile(flgm_imdir, strrep(tempimname, '.tif', '_region1.tif')));
imwrite(polarim_cond2, fullfile(flgm_imdir, strrep(tempimname, '.tif', '_region2.tif')));
