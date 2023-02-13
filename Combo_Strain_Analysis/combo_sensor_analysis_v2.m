%% New Combo Sensor Plotting 2021 %%

exptdir = input('Input full path to combo expt directory: ', 's');
cd(exptdir);
%% Compile or update Combo Sensor Data
cd(comboExptDir);
%% Load latest Data
allExptData = loadLatestComboSensorData();
%% Fill in iptgs/aras, as needed
allExptData = fillInIptgAra(allExptData);

%% 2/3/21: Fixing iptgs/aras using latest function

for i = 1:length(allExptData)
   [iptgs, aras] = getConcsFromFileNames(allExptData(i).imfiles); 
   dispcell = allExptData(i).imfiles';
   dispcell(:, 2) = num2cell(iptgs');
   dispcell(:, 3) = num2cell(aras');
   disp(dispcell);
   if input('Type 1 if okay: ')==1
       allExptData(i).iptgs = iptgs;
       allExptData(i).aras = aras;
   end
end
%%
for i = 1:length(allExptData)
    allExptData(i).iptgs = cell2mat(allExptData(i).iptgs);
    allExptData(i).aras = cell2mat(allExptData(i).aras);
end


%% 2/8/21: Fill in if colony has reached edge
cd(comboExptDir);
% Save zeros for 'reached_edge' field
figure;
for i = 11:11 %1:length(allExptData)
    fprintf('Expt %g of %g\n', i, length(allExptData));
    allExptData(i).reached_edge = zeros(1, length(allExptData(i).imfiles));
    % For each expt, display images in chunks of five with indices above each
    oldFolder = cd(allExptData(i).exptfolder);
    for j = 1:length(allExptData(i).imfiles)
        fprintf('Im %g of %g\n', j, length(allExptData(i).imfiles));
        tempim = imread(allExptData(i).imfiles{j});
        imshow(tempim);
        title(allExptData(i).imfiles{j}, 'Interpreter',"none");
        resp = input('Type 1 if colony covers whole plate, 2 if image is unusable: ');
        if resp == 1
            allExptData(i).reached_edge(j) = 1;
        elseif resp == 2
            allExptData(i).reached_edge(j) = NaN;
        end
    end

    % cd back to main directory
    cd(oldFolder);

end
%% Add polarim folders to allexptdata
for i = 1:length(allExptData)
    fileinfo = dir(fullfile('**', allExptData(i).exptfile));
    allExptData(i).polfolder = fileinfo.folder;
    
end
clear fileinfo
%% Try adding putative 'percent of agar area' msmts
% Use latest version of 'getColonyRegion' function on images which are not
% colonies that went to the edge (for those use 100% as the percent of agar
% area covered); for each expt, cd to the directory, load the analysis
% struct to get the petri masks, use each petri mask + getColRegion on each
% image to get the area, store in main allExptData struct

% fill in nans to start
for i = 1:length(allExptData)
    allExptData(i).percent_area = NaN(1, length(allExptData(i).imfiles));
end

%%
% iterate over experiment
use_imfindcircles = true;
threshold_on = false;
for i = 1:length(allExptData)
    fprintf('Expt %g of %g\n', i, length(allExptData));
  
    % iterate over the images
    for j = 1:length(allExptData(i).imfiles)
        fprintf('Im %g of %g\n', j, length(allExptData(i).imfiles));
        % skip previously done ones
        if isnan(allExptData(i).mean_cv(j))
            % load the polar image
            polname = allExptData(i).polfiles{j};
            if isempty(polname)
                continue;
            end
            testinfo = dir(fullfile('**', strcat('*', polname)));
            polarim = imread(fullfile(testinfo.folder, testinfo.name));
%                 polarim = imread(fullfile(allExptData(i).polfolder, ...
%                     allExptData(i).polfiles{j}));
            % Replace the Nans with white
            polarim = im2double(polarim);
            polarim(isnan(polarim)) = 1;
            if allExptData(i).reached_edge(j) == 1
                % it reached the edge
                allExptData(i).percent_area(j) = 1;
                colregion = polarim;
            else
                
                % Calculate the area
                [colmask, colregion, colarea_percent] = getColMaskPolar(polarim);
                
            end
            
            % Use the colony region to do calculations
            temp_cvs = allExptData(i).cvs{j};
            temp_std = allExptData(i).stdevs{j};
            
            msmts = getMsmts(colregion, temp_cvs, temp_std);

            % Store these calculations
             
            allExptData(i).percent_area(j) = colarea_percent;
            allExptData(i).mean_intens(j) = msmts{1};
            allExptData(i).min_intens(j) = msmts{2};
            allExptData(i).max_intens(j) = msmts{3};
            allExptData(i).mean_cv(j) = msmts{4};
            allExptData(i).mean_std(j) = msmts{5};

        end
    end

end
clear i j;

%% Fixing 10-14-19 expt

% Load the analysis and 's' structs first
% Let's save the polarimages!
imnames = {s.img_name};
for i = 17:26
    imname = allExptData(2).imfiles{i};
    og_imind = find(strcmpi(imname, imnames));
    polarim = s(og_imind).polarim;
    % Replace the Nans with white
    polarim(isnan(polarim)) = 1;
    polname = strrep(imname, '.tif', '.jpg');
    polname = strcat('10-14-19_', polname);
    allExptData(2).polfiles{i} = polname;
%     
%     imwrite(polarim, polname);
%     
%     % Now let's do msmts on this polarim
%     
%     if allExptData(2).reached_edge(i) == 1
%         % it reached the edge
%         allExptData(2).percent_area(i) = 1;
%         colregion = polarim;
%     else
%         % Calculate the area
%         [colmask, colregion, colarea_percent] = getColMaskPolar(polarim);
%     end
%     
%     imshow(colregion);
%     waitforbuttonpress;
%     % Use the colony region to do calculations
%     mean_intens = mean(mean(colregion(colregion~=1)));
%     min_intens = min(min(colregion(colregion~=1)));
%     max_intens = max(max(colregion(colregion~=1)));
%     col_traj = mean(colregion, 2);
%     temp_cvs = allExptData(2).cvs{i};
%     temp_std = allExptData(2).stdevs{i};
%     if iscell(temp_cvs)
%         temp_cvs = temp_cvs{1};
%         temp_std = temp_std{1};
%     end
%     inds = col_traj<1;
%     inds = inds(1:length(temp_cvs));
%     mean_cv = mean(temp_cvs(inds));
%     mean_stdev = mean(temp_std(inds));
% 
%     % Store these calculations
%     allExptData(2).percent_area(i) = colarea_percent;
%     allExptData(2).mean_intens(i) = mean_intens;
%     allExptData(2).min_intens(i) = min_intens;
%     allExptData(2).max_intens(i) = max_intens;
%     allExptData(2).mean_cv(i) = mean_cv;
%     allExptData(2).mean_std(i) = mean_stdev;
%     
end


%% Fixing 9-14-19 expt
% Load the analysis and 's' structs first
% Let's save the polarimages!
imnames = {s.img_name};
for i = 1:4 %13:34
%     if isnan(allExptData(11).mean_cv(i))
        disp(i);
        imname = allExptData(11).imfiles{i};
        disp(imname);
        polname = strrep(imname, '.tif', '.jpg');
        polname = strcat('09-14-19_', polname);
        fileinfo = dir(fullfile('**', strcat('*', polname)));
        og_imind = find(strcmpi(imname, imnames));
        ogiminfo = dir(fullfile('**', strcat('*', imname)));
        ogim = imread(fullfile(ogiminfo.folder, ogiminfo.name));
        imsize = size(ogim);
        imsize = imsize(1:2);
        polarim = s(og_imind).polarim;
        % Replace the Nans with white
        polarim(isnan(polarim)) = 1;
%         if isempty(fileinfo)
%             disp(i);
%     %         imshowpair(ogim, imresize(polarim, imsize), 'montage');
%     %         imshow(ogim);
%     %         tempcenter = [allExptData(11).colcenters(2*i-1), allExptData(11).colcenters(2*i)];
%     %         h = drawpoint('Position', tempcenter);
%     %         title(polname, 'Interpreter', 'none');
%     %         waitforbuttonpress;
%     %         delete(h);
%         end
        allExptData(11).polfiles{i} = polname;

        imwrite(polarim, polname);

        % Now let's do msmts on this polarim

%         if isnan(allExptData(11).reached_edge(i))
%             imshow(ogim);
%             resp = input('Type 1 if reached edge: ');
%             if resp==1
%                 allExptData(11).reached_edge(i) = 1;
%             else
%                 allExptData(11).reached_edge(i) = 0;
%             end
%         end
% 
%         if allExptData(11).reached_edge(i) == 1
%             % it reached the edge
%             allExptData(11).percent_area(i) = 1;
%             colregion = polarim;
%         else
%             % Calculate the area
%             [colmask, colregion, colarea_percent] = getColMaskPolar(polarim);
%         end
% 
%         imshow(colregion);
%         waitforbuttonpress;
%         % Use the colony region to do calculations
%         temp_cvs = allExptData(11).cvs{i};
%         temp_std = allExptData(11).stdevs{i};
%         msmts = getMsmts(colregion, temp_cvs, temp_std);
% 
%         % Store these calculations
% 
%         allExptData(11).percent_area(i) = colarea_percent;
%         allExptData(11).mean_intens(i) = msmts{1};
%         allExptData(11).min_intens(i) = msmts{2};
%         allExptData(11).max_intens(i) = msmts{3};
%         allExptData(11).mean_cv(i) = msmts{4};
%         allExptData(11).mean_std(i) = msmts{5};


%     end
    
end

%% Fix ALL min & max intensities

% for i = 1:length(allExptData)
%     if isa(allExptData(i).min_intens, 'uint8')
%         allExptData(i).min_intens = zeros(1, length(allExptData(i).min_intens));
%         allExptData(i).max_intens = zeros(1, length(allExptData(i).max_intens));
%     end
% end

for i = 11:11 %1:length(allExptData)
    fprintf('Expt %g of %g\n', i, length(allExptData));
  
    % iterate over the images
    for j = 1:length(allExptData(i).imfiles)
        if allExptData(i).max_intens(j) >0.999
            fprintf('Im %g of %g\n', j, length(allExptData(i).imfiles));

            polname = allExptData(i).polfiles{j};
            testinfo = dir(fullfile('**', strcat('*', polname)));
            polarim = imread(fullfile(testinfo.folder, testinfo.name));
            % Replace the Nans with white
            polarim = im2double(polarim);
            polarim(isnan(polarim)) = 1;
                if allExptData(i).reached_edge(j) == 1
                    % it reached the edge
                    allExptData(i).percent_area(j) = 1;
                    colregion = polarim;
                else
                    % Calculate the area
                    [colmask, colregion, colarea_percent] = getColMaskPolar(polarim);
                end

           % get min & max intensity
           colregion(colregion==1) = NaN;
           colregion(colregion==0) = NaN;
           allExptData(i).min_intens(j) = nanmin(nanmin(colregion));
           allExptData(i).max_intens(j) = nanmax(nanmax(colregion));
%            imshow(colregion);
%            disp(nanmin(nanmin(colregion)));
%            disp(nanmax(nanmax(colregion)));
%            waitforbuttonpress;
        end
    end
end


%% 6-21-21 Fix allExptData Exptfolder & Polfolder

for i = 1:length(allExptData)
    allExptData(i).exptfolder = strrep(allExptData(i).exptfolder, '/Anjali_Updates', '');
    allExptData(i).polfolder = strrep(allExptData(i).polfolder, '/Anjali_Updates', '');
    disp(i);
    disp(isfolder(allExptData(i).exptfolder));
end

%% 6-21-21 Redo all Colony Radii Manually (get the pix measurement)
if ~isfield(allExptData, 'colrads_pix_redo')
    allExptData(1).colrads_pix_redo = [];
end
fig_size = [622    72   747   707];
for i = 1:length(allExptData)
%     if isempty(allExptData(i).colrads_cm) & isempty(allExptData(i).colrads_pix)
    if isempty(allExptData(i).colrads_pix_redo) & ~isempty(allExptData(i).exptfolder)
        fprintf('Expt %g of %g\n', i, length(allExptData));
        oldFolder = cd(allExptData(i).exptfolder);
        if isempty(colrads_pix)
            colrads_pix = zeros(1, length(allExptData(i).imfiles));
        end
%         colrads_cm = colrads_pix;
        disp(i);
        % Iterate over images
        for j = 1:length(allExptData(i).imfiles)
            fprintf('Img %g of %g\n', j, length(allExptData(i).imfiles));
            if(colrads_pix(j)==0)
               % Load the scan
               tempim = rgb2gray(im2double(imread(allExptData(i).imfiles{j})));
               % Increase contrast for ease of display
               tempim = imadjust(tempim);
               % Get the center
               tempcolcent = [allExptData(i).colcenters((2*j) - 1), ...
                allExptData(i).colcenters(2*j)];
               % Draw circle
               temprad = getColRadManual(tempim, tempcolcent, [], fig_size);

               % store in colrads_pix
               colrads_pix(j) = temprad;
               % get in cm
    %            colrads_cm(j) = allExptData(i).dists{j}(round(temprad));
            end
        end
        allExptData(i).colrads_pix_redo = colrads_pix;
        colrads_pix = [];
%         allExptData(i).colrads_cm = colrads_cm;
    end
end
cd(oldFolder);
clear colrads_pix colrads_cm tempim tempcolcent temprad;


%% 6-21-21 Redo All Colrads cm
% Incorrect to use colrad pixels to get radius in CM unless the radius was
% measured *AFTER* flattening, which for most images it wasn't. Therefore,
% in order to convert from pixels to cm, we need just the DPCM (usually
% 400/2.54). We will do a sanity check by showing the original image & a
% circle around the center of the approx radius. If it looks good, we will
% do the conversion.
% First create a way to track this as we go: 
redone = cell(1, length(allExptData));
for i = 1:length(allExptData)
    redone{i} = zeros(1, length(allExptData(i).imfiles));
end
%%
% iterate over expts
if ~isfield(allExptData, 'colrads_cm_redo')
    allExptData(1).colrads_cm_redo = [];
end
for i = 1:length(allExptData)
    if ~isempty(allExptData(i).exptfolder)
        fprintf('Expt %g of %g\n', i, length(allExptData));
        tic;
        % load the resolution
        testinfo = imfinfo(fullfile(allExptData(i).exptfolder, allExptData(i).imfiles{1})); 
        if isfield(testinfo, 'XResolution')
            res = testinfo.XResolution;
        else
            disp(allExptData(i).exptfile);
            resp = input('Type the resolution in pixels: ');
            res = resp;
        end
        for j = 1:length(allExptData(i).imfiles)
            fprintf('Im %g of %g\n', j, length(allExptData(i).polfiles));
            if ~isempty(allExptData(i).exptfolder)
                if redone{i}(j) == 0
                    % iterate over images
                    temprad = allExptData(i).colrads_pix_redo(j);
                    temprad_cm = temprad/(res/2.54);
                    % Store
                    allExptData(i).colrads_cm_redo(j) = temprad_cm;
                    redone{i}(j) = 1;
                end
            end
        end
    end
end



%% Try Sliding Window CV Method

% Instead of getting the average CV for the polarims across each row, try
% getting a CV in a certain width rectangular window that slides across the
% row, then get the average of that--so this will be local CV
if ~isfield(allExptData, 'slide_CVs')
    allExptData(1).slide_CVs = [];
end

for i = 1:length(allExptData)
    if isempty(allExptData(i).slide_CVs)
        allExptData(i).slide_CVs = cell(1, length(allExptData(i).imfiles));
        allExptData(i).slide_CVs(:) = {NaN};
        allExptData(i).slide_CV_means = zeros(1, length(allExptData(i).imfiles));
    end
end

%% Redoing sliding CV Mean Calculation using Colony Radii
cd(comboExptDir);
winwidth = 10;
for i = 1:length(allExptData)
    if ~isempty(allExptData(i).exptfolder)
        fprintf('Expt %g of %g\n', i, length(allExptData));
        tic;
        % load the image
        for j = 1:length(allExptData(i).imfiles)
        fprintf('Im %g of %g\n', j, length(allExptData(i).polfiles));
            if isnan(allExptData(i).slide_CVs{j})
                iminfo = dir(fullfile('**', allExptData(i).polfiles{j}));
                polarim = imread(fullfile(iminfo.folder, iminfo.name));
                [cvs_all, cv_traj] = getRollingCV(polarim, winwidth);
                % Store
                allExptData(i).slide_CVs{j} = cv_traj;
%                 cv_traj = allExptData(i).slide_CVs{j};
                cv_traj(cv_traj==0) = NaN;
                % Now let's check if the whole plate is not covered
                if allExptData(i).reached_edge ~= 1
                    % The whole plate is not covered. let's get the colony
                    % radius and use it to calculate the cvs
                    temprad = allExptData(i).colrads_cm_redo(j);
                    tempdist = allExptData(i).dists{j};
                    temprad_pix = find(tempdist>temprad, 1);
                    allExptData(i).slide_CV_means(j) = nanmean(cv_traj(1:temprad_pix));
                else
                    allExptData(i).slide_CV_means(j) = nanmean(cv_traj);
                end
            end
        end
        toc;
    end
end

%% Redoing Colony Min Intensity Calculation

% Get the median of the darkest 100 pixels for each colony, after masking
% out, & store
% First put in empty values
if ~isfield(allExptData, 'med_min_intens')
    allExptData(1).med_min_intens = [];
end
for i = 1:length(allExptData)
    if isempty(allExptData(i).med_min_intens)
        allExptData(i).med_min_intens = NaN(1, length(allExptData(i).imfiles));
    end
end
%%
% iterate over every image
for i = 1:length(allExptData)

    fprintf('Expt %g of %g\n', i, length(allExptData));
    if ~isempty(allExptData(i).exptfolder)
        tic;
        for j = 1:length(allExptData(i).imfiles)
        fprintf('Im %g of %g\n', j, length(allExptData(i).polfiles));
            % check if already done
            if isnan(allExptData(i).med_min_intens(j))
                iminfo = dir(fullfile('**', allExptData(i).polfiles{j}));
                polarim = imread(fullfile(iminfo.folder, iminfo.name));
                % get col region to calculate min intens over
                [colmask, colregion, colarea_percent] = getColMaskPolar(polarim);
                colregion = colregion(:)';
                colregion(colregion==0) = NaN;
                minpix = mink(colregion, 50);
                allExptData(i).med_min_intens(j) = median(minpix);
            end
        end
        toc;
    end
end

%% 7-1-21: Redoing all combo sensor image flattening

% Create a field in the struct or something to store whether image has been re-flattened
if ~isfield(allExptData, 'reflattened')
    allExptData(1).reflattened = [];
end
for i = 1:length(allExptData)
    if isempty(allExptData(i).reflattened)
        allExptData(i).reflattened = NaN(1, length(allExptData(i).imfiles));
    end
end

% Set constants
interp_val = 1000;
use_imfindcircles = true;
threshold_on = false;

% first expt on the 14th is done
allExptData(1).reflattened = ones(1, length(allExptData(1).imfiles));
%%
% Iterate over the expts

% iterate over every image
for i = 1:length(allExptData)
    fprintf('Expt %g of %g\n', i, length(allExptData));
    if ~isempty(allExptData(i).exptfolder)
        tic;
        % cd into the folder
        oldFolder = cd(allExptData(i).exptfolder);
        % make directory if it doesn't yet exist
        exptdate = allExptData(i).date;
        tempdate = string(exptdate, 'MM-dd-uuuu');
        tempdate = strrep(tempdate, '00', '');
        subdir = strcat(tempdate, '_polarims');
        if ~isfolder(subdir)
            mkdir(subdir);
        end
        
        % get or set resolution
        testinfo = imfinfo(fullfile(allExptData(i).exptfolder, allExptData(i).imfiles{1})); 
        if isfield(testinfo, 'XResolution')
            res = testinfo.XResolution;
        else
            disp(allExptData(i).exptfile);
            resp = input('Type the resolution in pixels: ');
            res = resp;
        end
        dpcm = res/2.54;
        
        % Iterate over the images
        for j = 1:length(allExptData(i).imfiles)
            fprintf('Im %g of %g\n', j, length(allExptData(i).polfiles));
            % check if already done
            if isnan(allExptData(i).reflattened(j))
               tempfile = allExptData(i).imfiles{j};
                % load the scanned image
               tempim = rgb2gray(im2double(imread(tempfile)));
                % get the colony center
                tempcenter = [allExptData(i).colcenters(2*j - 1),...
                    allExptData(i).colcenters(2*j)];
                %Get the mask to remove rim
                [maskim, petrimask] = removeRimAuto(tempim, use_imfindcircles, threshold_on);
                
                % Flatten image & save the flattened image file
                disp('Flattening...');
                [polar_im, dist_vec] = flattenColonyInterp(maskim, ...
                    tempcenter, dpcm, interp_val);
                % make the name of the polarimage--get from allExptData

                tempname = erase(tempfile, '.tif');
                polarimname = strcat(tempname, '_polarim.tif');
                % imwrite to the subdir
                imwrite(polar_im, fullfile(subdir, polarimname));
                
                % Do the basic analysis
                polar_im = im2double(polar_im); %in case if not double
                avgd = mean(polar_im, 2);
                stdevs = std(polar_im, 0, 2);
                cvs = stdevs./avgd;
                
                % save the new data
                allExptData(i).avgd{j} = avgd;
                allExptData(i).cvs{j} = cvs;
                allExptData(i).stdevs{j} = stdevs;
                allExptData(i).dists{j} = dist_vec;  
                
                % mark as done
                allExptData(i).reflattened(j) = 1;
            end
            allExptData(i).polfolder = fullfile(allExptData(i).exptfolder, subdir);
        end
        toc;
        % cd out of the folder
        cd(oldFolder);
    end
end

cd(comboExptDir);


%% Adding Inoc Msmts to Combo Sensor allExptData

load('/Users/anjali/Dropbox/Patterning_Expts_Analysis/Experiments/placcheW_pbadumod_expts/Main_Expts/cheWumoD_inoc_features_07-05-21.mat');
%%
% iterate over each folder in the thing
for i = 1:length(inocFeatures)
    tempfold = inocFeatures(i).ExptFolder;
    tempfold = erase(tempfold, '_cheWumoD');
    for j = 1:length(allExptData)
        if contains(allExptData(j).polfolder, tempfold)
            tempind = j;
            disp(tempind);
            % now get the measurements from inocfeatures and add to
            % allexptdata
            allExptData(j).inoc_meanintens = inocFeatures(i).Inoc_MeanIntensity;
            allExptData(j).inoc_minintens = inocFeatures(i).Inoc_MinIntensity;
            allExptData(j).inoc_maxintens = inocFeatures(i).Inoc_MaxIntensity;
            allExptData(j).inoc_area = inocFeatures(i).Inoc_Area;
            allExptData(j).inoc_minoraxlen = inocFeatures(i).Inoc_MinorAxisLength;

        end
    end
    
end

%% Fixing pol file names

for i = 1:length(allExptData)
    for j = 1:length(allExptData(i).polfiles)
        tempname_pol = allExptData(i).polfiles{j};
        tempfile = allExptData(i).imfiles{j};
        %get the name
        tempname = erase(tempfile, '.tif');
        polarimname = strcat(tempname, '_polarim.tif');
        % erase jpg
        % find it
        im_info = dir(fullfile('**', polarimname));
        if ~isempty(im_info)
%             disp('found');
            if ~strcmp(allExptData(i).polfiles{j}, polarimname)
                allExptData(i).polfiles{j} = polarimname;
                disp('not saved');
            end
        else
            disp(polarimname);
        end
%         dir(fullfile('**', allExptData(i).exptfile));
        % save the new name in allExpt data
        % make sure to save the latest allExptData afterwards
        
    end
end

%% Fixing all col percent areas as of August 2021 using 'getColonyMaskPolar' function

for i = 1:length(allExptData)
    disp(i);
    for j = 1:length(allExptData(i).polfiles)
        fprintf('Image %g of %g\n', j, length(allExptData(i).polfiles));
        if allExptData(i).reached_edge(j) ~= 1
            polarim = imread(fullfile(allExptData(i).polfolder, allExptData(i).polfiles{j}));
            polarim = im2double(polarim);
%             [~, colmask, colregion, ~] = getGroundTruthColRegion(polarim);
            [colmask1, colregion2, colarea_percent] = getColMaskPolar(polarim);
%             imshowpair(colregion, colregion2, 'montage');
            if abs(allExptData(i).percent_area(j)-colarea_percent)>0.1
                disp(allExptData(i).percent_area(j));
                disp(colarea_percent);
                imshowpair(polarim, colregion2, 'montage');
                title(sprintf('New area %g', colarea_percent));
%                 waitforbuttonpress;
                % If it's more accurate (the mask area), then hit enter,
                % otherwise press 1 to not update
                resp = input('Type 1 not to update the percent area: ');
                if resp ~= 1
                    allExptData(i).percent_area(j) = colarea_percent;
                end
            end
        end
    end
end




%% CHECKING what the Colony Masks Look Like at intermediate condition
figure(7);
clf('reset');
for i = 1:length(allExptData)
    disp(i);
    for j = 1:length(allExptData(i).polfiles)
        fprintf('Image %g of %g\n', j, length(allExptData(i).polfiles));
        if allExptData(i).reached_edge(j) ~= 1
            polarim = imread(fullfile(allExptData(i).polfolder, allExptData(i).polfiles{j}));
            polarim = im2double(polarim);
%             [~, colmask, colregion, ~] = getGroundTruthColRegion(polarim);
            [colmask1, colregion2, colarea_percent] = getColMaskPolar(polarim);
%             imshowpair(colregion, colregion2, 'montage');
            if allExptData(i).iptgs(j) == 2.5 && allExptData(i).aras(j) == 0.1
%                 disp(allExptData(i).percent_area(j));
%                 disp(colarea_percent);
                imshowpair(polarim, colregion2, 'montage');
                title(sprintf('New area %g', colarea_percent));
                disp(allExptData(i).polfiles{j});
                waitforbuttonpress;
                
            end
        end
    end
end

%%
% Save updated version
disp('Saving updated version...');
filename = strcat('all_combo_expt_data_', date, '.mat');
save(filename, 'allExptData', '-v7.3');


%% FUNCTIONS
function msmts = getMsmts(colregion, temp_cvs, temp_std)
% Use the colony region to do calculations
    mean_intens = nanmean(nanmean(colregion(colregion~=1)));
    min_intens = nanmin(nanmin(colregion(colregion~=1)));
    max_intens = nanmax(nanmax(colregion(colregion~=1)));
    col_traj = mean(colregion, 2);

    if iscell(temp_cvs)
        temp_cvs = temp_cvs{1};
        temp_std = temp_std{1};
    end
    inds = col_traj<1;
    inds = inds(1:length(temp_cvs));
    mean_cv = nanmean(temp_cvs(inds));
    mean_stdev = nanmean(temp_std(inds));
    
    msmts{1} = mean_intens;
    msmts{2} = min_intens;
    msmts{3} = max_intens;
    msmts{4} = mean_cv;
    msmts{5} = mean_stdev;
            
end
