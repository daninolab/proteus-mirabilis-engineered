%% New Single Gene Sensor Plotting 2021 %%

% exptdir = input('Input full path to combo expt directory: ', 's');
cd(singleExptDir);


%% Load latest Data
allExptData = loadLatestSingleGeneData();

%% Check if duplicate dates

% the 4-19 expt was duplicated; save the later one
% remove '4-9-19 analysis_polar.mat'
ind = find(strcmpi({allExptData.exptfile}, '4-9-19_WT_analysis_polar.mat'));
allExptData(ind) = [];

%% Let's also remove the 5 hour dataset--not relevant for analysis
ind = find(strcmpi({allExptData.exptfile}, '11-14-19_umoD_5hr_analysis_polar.mat'));
allExptData(ind) = [];

%% Fill in info as needed
allExptData = fillInGene(allExptData);

%% Fill in IPTGs from Previously Done Analysis-Don't redo since we have done all the analysis since then

% % Use previously done analysis to do the initial filling in of missing
% % parameters such as IPTG, petrimask
% load('all_genes_analysis_9_11_20.mat', 'all_expts_table');
% for i = 1:length(allExptData)
%     iptgs = zeros(1, length(allExptData(i).imfiles));
%     % Iterate over images
%     for j = 1:length(allExptData(i).imfiles)
%         % Find image file in all_expts_table
%         im_curr = allExptData(i).imfiles{j};
%         if any(find(strcmpi(all_expts_table.filename, im_curr)))
%         % If it is there, get the iptg
%         iptgs(j) = ...
%             all_expts_table.IPTG(find(strcmpi(all_expts_table.filename, im_curr), 1));
%         end
%     end
%     allExptData(i).iptgs = iptgs;
% end
% clear all_expts_table;

%% Fill In New Experiment IPTGs
% Find iptgs which are all 0s
allExptData = fillInIPTGs(allExptData);


%% Double check iptgs and manually fix any bad ones
test = checkiptgs(allExptData);

% To fix:
allExptData(15).iptgs(1) = 0;
allExptData(15).iptgs(9) = 5;
allExptData(15).iptgs(22) = 5;
allExptData(21).iptgs(5) = 5;
allExptData(21).iptgs(17) = 5;

%% 5/13/21: Fill in Iptgs for the New WT expts
for i = 22:23
    numims = length(allExptData(i).imfiles);
    allExptData(i).iptgs = zeros(1, numims);
end
%% Find the ones where ara isn't 0 for the 9-17-19 expt
tempfiles = allExptData(22).imfiles;
tempiptgs = zeros(1, length(tempfiles));
discardims = tempiptgs;
for i = 1:length(tempfiles)
    tempfile = tempfiles{i};
    tempfile = split(tempfile, '_');
    tempiptg = str2double(erase(tempfile{3}, 'i'));
    tempara = erase(tempfile{4}, 'a');
    tempara = strrep(tempara, 'pt', '.');
    tempara = str2double(tempara);
    
    tempiptgs(i) = tempiptg;
    if tempara ~= 0
        discardims(i) = 1;
    end
end
useims = ~logical(discardims);
allExptData(22).iptgs = tempiptgs;
%% Discard ara images from 9-17-19
% Use the discardims inds to retrieve the correct images
fnames = fieldnames(allExptData);
for i = 2:length(fnames)
    tempdata = allExptData(22).(fnames{i});
    disp(fnames{i});
    fprintf('Current data is length %g\n', length(tempdata));
    if length(tempdata)==31
        tempdata2 = tempdata(useims);
        allExptData(22).(fnames{i}) = tempdata2;
    elseif length(tempdata)==62
        tempdata2 = tempdata(repelem(useims, 2));
        allExptData(22).(fnames{i}) = tempdata2;
    else
        tempdata2 = [];
    end
    
    fprintf('New data is length %g\n', length(tempdata2));
    if strcmpi(fnames{i}, 'reached_edge');
        break;
    end
end

%% Find the 21 hour ones for 10-23-18 and discard
% Images 1-4 are the problem--remove
fnames = fieldnames(allExptData);
for i = 2:length(fnames)
    tempdata = allExptData(23).(fnames{i});
    disp(fnames{i});
    fprintf('Current data is length %g\n', length(tempdata));
    if length(tempdata)==8
        allExptData(23).(fnames{i}) = tempdata(5:8);
    elseif length(tempdata)==16
        allExptData(23).(fnames{i}) = tempdata(9:16);

    end
    if strcmpi(fnames{i}, 'reached_edge');
        break;
    end
end



%% Fill in Missing Petri Masks for 4-9-19 & fix any others desired
% exptind = 8;
for exptind = 21:21 %1:length(allExptData)
    fprintf('Expt %g of %g\n', exptind, length(allExptData));
    if isempty(allExptData(exptind).scanfolder)
        continue;
    end
    % Fix the expt masks

    for i = 1:length(allExptData(exptind).imfiles)
        tempsize = size(imread(fullfile(allExptData(exptind).scanfolder, ...
                allExptData(exptind).imfiles{i})));
        % Check if image size = petri mask size
        if ~isequal(size(allExptData(exptind).petrimasks{i}), ...
                tempsize(1:2))
            % If not, regenerate mask & store
            tempim = imread(fullfile(allExptData(exptind).scanfolder, ...
                allExptData(exptind).imfiles{i}));
            [~, allExptData(exptind).petrimasks{i}] = ...
                removeRimAuto(tempim);
        end
        clear tempim;
    end
end
clear tempsize tempmasks exptind i

%% Fill in folder with scanned images for each expt
for i = 1:length(allExptData)
    oldFolder = cd(allExptData(i).exptfolder);
    cd ../
    tempinfo = dir('*scans');
    if ~isempty(tempinfo)
        allExptData(i).scanfolder = fullfile(tempinfo.folder, tempinfo.name);
    end
    
    cd(oldFolder);
end

%% Fix scanfolders for new WT expts
for i = 22:24
    tempfolder = allExptData(i).scanfolder;
    tempfolder = strcat('/', tempfolder);
    allExptData(i).scanfolder = tempfolder;
end

%% Fill in genes for new WT expts
for i = 6 %22:24
    genes_cell = cell(1, length(allExptData(i).imfiles));
    genes_cell(:) = {'wt'};
    allExptData(i).genes = genes_cell;
end

%% Fill in outliers for new WT expts--call them all non outliers
for i = 22:24
    allExptData(i).outlier = zeros(1, length(allExptData(i).imfiles));
end

%% Fix scanfolders
for i = 24:28
    allExptData(i).scanfolder = strcat('/', allExptData(i).scanfolder);
end

%% fix chew scanfolder
allExptData(27).scanfolder = '/Users/anjali/Dropbox/Patterning_Expts_Analysis/Experiments/Single_Gene_Expts/10-1-19-pLaccheW/10-1-19_scans';

%% 2/8/21: Fill in if colony has reached edge
cd(singleExptDir);
% Save zeros for 'reached_edge' field
figure;
for i = 1:length(allExptData)
    fprintf('Expt %g of %g\n', i, length(allExptData));
    if length(allExptData(i).reached_edge) == length(allExptData(i).imfiles)
        % Expt already done
        continue;
    end
    if ~isempty(allExptData(i).scanfolder)
%         allExptData(i).reached_edge = nan(1, length(allExptData(i).imfiles));
        oldFolder = cd(allExptData(i).scanfolder);
        % Displaying
        for j = 1:length(allExptData(i).imfiles)
            % skip ones already done
            if length(allExptData(i).reached_edge) >= j+1
                continue;
            end
            fprintf('Im %g of %g\n', j, length(allExptData(i).imfiles));
            tempim = imread(allExptData(i).imfiles{j});
            tempim = rgb2gray(im2double(tempim));
            tempmask = allExptData(i).petrimasks{j};
            tempim(tempmask==0)=1;
            imshow(adapthisteq(tempim), []);
%             imshow(tempim);
            title(allExptData(i).imfiles{j}, 'Interpreter',"none");
            resp = input('Type 1 if colony covers whole plate, 2 if image is unusable: ');
            if resp == 1
                allExptData(i).reached_edge(j) = 1;
            elseif resp == 2
                allExptData(i).reached_edge(j) = NaN;
            else
                allExptData(i).reached_edge(j) = 0;
            end
        end

        % cd back to main directory
        cd(oldFolder);
    end
end

%% Add putative 'percent of agar area' field
% Use latest version of 'getColonyRegion' function on images which are not
% colonies that went to the edge (for those use 100% as the percent of agar
% area covered); for each expt, cd to the directory, load the analysis
% struct to get the petri masks, use each petri mask + getColRegion on each
% image to get the area, store in main allExptData struct

% fill in nans to start
for i = 1:length(allExptData)
    if isempty(allExptData(i).percent_area)
        allExptData(i).percent_area = NaN(1, length(allExptData(i).imfiles));

    end
end

%% Add putative 'percent of agar area' msmts from experiments
use_imfindcircles = true;
threshold_on = false;

for i = 1:length(allExptData)
    if ~isempty(allExptData(i).scanfolder)
        fprintf('Expt %g of %g\n', i, length(allExptData));
        if any(isnan(allExptData(i).percent_area))
            % cd to the directory
            oldFolder = cd(allExptData(i).exptfolder);
            % load the polar analysis struct

            % iterate over the images
            for j = 1:length(allExptData(i).imfiles)
                fprintf('Im %g of %g\n', j, length(allExptData(i).imfiles));
                % skip previously done ones
                if isnan(allExptData(i).percent_area(j))
                    if allExptData(i).reached_edge(j) == 1
                        % it reached the edge
                        allExptData(i).percent_area(j) = 1;
                    else
                        imname = allExptData(i).imfiles{j};
        %                 % load the petri mask
                        petri_mask = allExptData(i).petrimasks{j};
                        tempim = imread(fullfile(allExptData(i).scanfolder, imname));
                        % get rimless
                        rimless = rgb2gray(im2double(tempim));
                        rimless(petri_mask==0) = 1; 

                        % area in decimals (eg 50% = 0.5)
                        [~, ~, colarea_percent] = getColonyRegion(rimless, petri_mask);
                        allExptData(i).percent_area(j) = colarea_percent;
                    end
                end
            end
            % move back to main directory
            cd(oldFolder);
        end
    end
end
clear tempim imname ind petri_mask agar_area rimless colarea_percent;
clear i j;

%% Fix allexpt genes for 4-9-19 wt expt

for i = 1:length(allExptData(6).imfiles)
    allExptData(6).genes{i} = 'wt';
end

%% Fix iptgs
allExptData(5).iptgs = 0.1 * ones(1, 4);

%% Fix 4-9-19 folder
allExptData(6).exptfolder = ...
    '/Users/anjali/Dropbox/Anjali_Updates/Patterning_Expts_Analysis/Experiments/Single_Gene_Expts/4-9-19_WT_iptg/WT_4-9-19_polarims';

%% Fix lrp 9-25-19 polarim names
allExptData(21).polfiles{4} = 'pLacGFP_0mM_IPTG_012_polarim.tif';
allExptData(21).polfiles{16} = 'pLaclrp_0mM_IPTG_017_polarim.tif';
%% Get approx colony symmetry for each polar image

for i = 1:length(allExptData)
    if ~isempty(allExptData(i).scanfolder)
        fprintf('Expt %g of %g\n', i, length(allExptData));
        % cd to the directory
        oldFolder = cd(allExptData(i).exptfolder);
        tic;
        % iterate over the images
        for j = 1:length(allExptData(i).polfiles)
            fprintf('Im %g of %g\n', j, length(allExptData(i).polfiles));
            % skip previously done ones
            if length(allExptData(i).binim_cvtrajs)>(j)
                continue;
            end
            % load polarim
            polarim = imread(fullfile(allExptData(i).exptfolder, allExptData(i).polfiles{j}));
            
            % Get CV & STDtrajs
            cvtraj2 = approxColSym2(polarim);

            %Store
            allExptData(i).binim_cvtrajs{j} = cvtraj2;
            

        end
        % move back to main directory
        cd(oldFolder);
        toc;
    end
end
clear oldFolder i j polarim CVtraj STDtraj;

%% 3-17-21 Fill in Colony Radii from Previous Analysis

% Don't redo this if we have the latest single gene expt data analysis

% % Fill in IPTGs from Previously Done Analysis
% % Use previously done analysis to do the initial filling in of missing
% % parameters such as IPTG, petrimask
% load('all_genes_analysis_9_11_20.mat', 'all_expts_table');
% for i = 1:length(allExptData)
%     radii = zeros(1, length(allExptData(i).imfiles));
%     % Iterate over images
%     for j = 1:length(allExptData(i).imfiles)
%         % Find image file in all_expts_table
%         im_curr = allExptData(i).imfiles{j};
%         if any(find(strcmpi(all_expts_table.filename, im_curr)))
%             % check if it has the same colony center as the one in all
%             % expts table
%             im_inds = find(strcmpi(all_expts_table.filename, im_curr));
%             if length(im_inds)==1
%             % If it is there, get the iptg
%             radii(j) = all_expts_table.col_radii(im_inds);
% %                 all_expts_table.col_radii(find(strcmpi(all_expts_table.filename, im_curr), 1));
%                 
%             else
%                 % iterate over
%                 for k = 1:length(im_inds)
%                     tempcent = [allExptData(i).colcenters((2*j)-1), ...
%                         allExptData(i).colcenters(2*j)];
%                     tempind = im_inds(k);
%                     if all_expts_table.colcenter(tempind, :) == tempcent
%                         radii(j) = all_expts_table.col_radii(tempind);
%                     else
%                         continue;
%                     end
%                 end
%             end
%         end
%     end
%     if any(radii)
%         allExptData(i).colrads_cm = radii;
%     else
%         allExptData(i).colrads_cm = [];
%     end
% end
% clear all_expts_table;


%% For expts where we have the radius in pixels, use dist vec to get rad in cm
for i = 1:length(allExptData)
    
   if isempty(allExptData(i).colrads_cm)
       if ~isempty(allExptData(i).colrads_pix)
           temprads = allExptData(i).colrads_pix;
           tempdists = allExptData(i).dists;
           temprads_cm = zeros(size(temprads));
           % Iterate over images
           for j = 1:length(temprads)
               temprad = temprads(j);
               tempdist = tempdists{j};
               temprads_cm(j) = tempdist(round(temprad));
           end
           allExptData(i).colrads_cm = temprads_cm;
       end
   end
end
clear temprads tempdists temprads_cm temprad tempdist

%% Fill in any radii left over
% allExptData(1).colrads_pix_redo = [];
for i = 1:length(allExptData)
%     if isempty(allExptData(i).colrads_cm) & isempty(allExptData(i).colrads_pix)
    if isempty(allExptData(i).colrads_pix_redo) & ~isempty(allExptData(i).scanfolder)
        fprintf('Expt %g of %g\n', i, length(allExptData));
        oldFolder = cd(allExptData(i).scanfolder);
        colrads_pix = zeros(1, length(allExptData(i).imfiles));
%         colrads_cm = colrads_pix;
        disp(i);
        % Iterate over images
        for j = 1:length(allExptData(i).imfiles)
            fprintf('Img %g of %g\n', j, length(allExptData(i).imfiles));
           % Load the scan
           tempim = rgb2gray(im2double(imread(allExptData(i).imfiles{j})));
           % Increase contrast for ease of display
           tempim = imadjust(tempim);
           % Get the center
           tempcolcent = [allExptData(i).colcenters((2*j) - 1), ...
            allExptData(i).colcenters(2*j)];
           % Draw circle
           temprad = getColRadManual(tempim, tempcolcent);
           
           % store in colrads_pix
           colrads_pix(j) = temprad;
           % get in cm
%            colrads_cm(j) = allExptData(i).dists{j}(round(temprad));
            
        end
        allExptData(i).colrads_pix_redo = colrads_pix;
%         allExptData(i).colrads_cm = colrads_cm;
    end
end
cd(oldFolder);
clear colrads_pix colrads_cm tempim tempcolcent temprad;

%% Fix the col rads for pLaclrp_0mM_IPTG_017.tif
tempim = imread(fullfile(allExptData(21).exptfolder, allExptData(21).polfiles{16}));
tempdist = allExptData(21).dists{16};
temprad = getColRadManual_Polar(tempim, tempdist);

allExptData(21).colrads_cm(16) = temprad;

%% Fix for 10-20-20 expt
for i = 1:length(allExptData(10).imfiles)
    tempim = imread(fullfile(allExptData(10).exptfolder, allExptData(10).polfiles{i}));
    tempdist = allExptData(10).dists{i};
    temprad = getColRadManual_Polar(tempim, tempdist);

    allExptData(10).colrads_cm(i) = temprad;    
    
    
end

%% Add 'Colony Area Accuracy'
% measure how 'accurate' the colony area thing is

% iterate over every image
for i = 1:length(allExptData)
    for j = 1:length(allExptData(i).imfiles)
    % get the percent area
    percent_area = allExptData(i).percent_area(j);
    % get the petri mask
    petri_mask = allExptData(i).petrimasks{j};
    % get the radius 
    colrad_cm = allExptData(i).colrads_cm(j);
    tempdist = allExptData(i).dists{j};
    colrad = find(tempdist>colrad_cm, 1);
    % get percent of the petri dish based on that mask
    colarea = pi*(colrad^2);
    % get petri mask idsk area
%     length(find(petri_mask));
    props = regionprops(petri_mask, 'Area');
    maskarea = props.Area;
    col_percent_area = colarea/maskarea;
    
    % calculate how close the estimated one was to the measured radius by
    % subtracting the two and dividing the difference by the measured 'area'
    % store
    acc = 1-abs((percent_area-col_percent_area)/col_percent_area);
    allExptData(i).area_acc(j) = acc;
    
    end
end

%% Fix Dates

for i = 1:length(allExptData)
    newdate = replace(string(allExptData(i).date), '00', '20');
    allExptData(i).date = datetime(newdate);
end

%% Identify Outliers
figure(4);
clf('reset');
% iterate over the dataset
for i = 1:length(allExptData)
    if isempty(allExptData(i).outlier)
        allExptData(i).outlier = NaN(1, length(allExptData(i).imfiles));
    end
%     allExptData(i).outlier = NaN(1, length(allExptData(i).imfiles));
    disp(i);
    disp(allExptData(i).date);
    for j = 1:length(allExptData(i).imfiles)
        if ~isempty(allExptData(i).scanfolder)
            % check if already done
            if isnan(allExptData(i).outlier(j))
                % read image
                tempim = imread(fullfile(allExptData(i).scanfolder, allExptData(i).imfiles{j}));
                % display image
                clf('reset');
                imshow(tempim);
                title(allExptData(i).imfiles{j}, 'Interpreter', 'none');
                % type 1 if outlier
                resp = input('Is this an outlier? Type 1 if outlier: ');
                if resp==1
                    allExptData(i).outlier(j) = 1;
                else
                    allExptData(i).outlier(j) = 0;
                end
            end
        end
    end
end

%% Measure Freq Magnitude in Fourier Domain
% Transform the polar image to fourier domain, then measure the magnitude
% in a certain low-but-not-lowest frequency region.

% fill in nans to start
for i = 1:length(allExptData)
    if isempty(allExptData(i).ftregmag)
        allExptData(i).ftregmag = NaN(1, length(allExptData(i).imfiles));

    end
end

%%
% iterate over every image
for i = 1:length(allExptData)
    if isempty(allExptData(i).ftregmag)
        allExptData(i).ftregmag = NaN(1, length(allExptData(i).imfiles));
    end
    fprintf('Expt %g of %g\n', i, length(allExptData));
%     disp(i);
%     disp(allExptData(i).date);
    tic;
    for j = 1:length(allExptData(i).imfiles)
        fprintf('Im %g of %g\n', j, length(allExptData(i).polfiles));
        if ~isempty(allExptData(i).scanfolder)
            % check if already done
            if isnan(allExptData(i).ftregmag(j))
                
                % Calculate the ft region mag
                polarim = imread(fullfile(allExptData(i).exptfolder, ...
                    allExptData(i).polfiles{j}));
                [colmask, colregion] = getColMask(polarim);
                fftim = fftshift(fft2(colregion));
                mag_fftim = abs(fftim);
                ftregmag = mean(mean(mag_fftim(400:480, 400:480)));
                % Store in allexptdata:
                allExptData(i).ftregmag(j) = ftregmag;
                
            end
        end
    end
    toc;
end

%% Fixing the exptfolders
for i = 1:length(allExptData)
    if ~isempty(allExptData(i).scanfolder)
       tempfolder = allExptData(i).exptfolder;
       tempfolder2 = allExptData(i).scanfolder;
       tempfolder = erase(tempfolder, 'Anjali_Updates/');
       tempfolder2 = erase(tempfolder2, 'Anjali_Updates/');
       allExptData(i).scanfolder = tempfolder2;
       allExptData(i).exptfolder = tempfolder;
    end
end

%% Get Region Magnitudes over First Quadrant of FT Mag Im
% fill in nans to start
for i = 1:length(allExptData)
    if isempty(allExptData(i).ftregmags)
        allExptData(i).ftregmags = cell(1, length(allExptData(i).imfiles));
        allExptData(i).ftregmags(:) = {NaN};
    end
end
%%
% iterate over every image
for i = 1:length(allExptData)
%     if isempty(allExptData(i).ftregmags)
%         allExptData(i).ftregmags = cell(1, length(allExptData(i).imfiles));
%     end
    fprintf('Expt %g of %g\n', i, length(allExptData));
%     disp(i);
%     disp(allExptData(i).date);
    tic;
    for j = 1:length(allExptData(i).imfiles)
        fprintf('Im %g of %g\n', j, length(allExptData(i).polfiles));
        if ~isempty(allExptData(i).scanfolder)
            % check if already done
            if isnan(allExptData(i).ftregmags{j})
                
                % Calculate the ft region mag
                polarim = imread(fullfile(allExptData(i).exptfolder, ...
                    allExptData(i).polfiles{j}));
                [colmask, colregion] = getColMask(polarim);
                fftim = fftshift(fft2(colregion));
                mag_fftim = abs(fftim);
                regionmag_array = getFTRegionMags(mag_fftim);
                % Store in allexptdata:
                allExptData(i).ftregmags{j} = regionmag_array;
                
            end
        end
    end
    toc;
end

%% 06/03/21-Filling In missing colrads pix

% For each expt, if colrads_pix is empty, load the exptfile. Then look for
% a field for col radii?
for i = 10:10 %1:length(allExptData)
    curr_file = fullfile(allExptData(i).exptfolder, allExptData(i).exptfile);
    if exist(curr_file)==0
        % need to find it
        tempinfo = dir(fullfile('**', strcat('*', allExptData(i).exptfile)));
        curr_file = fullfile(tempinfo.folder, tempinfo.name);
    end
    listOfVariables = who('-file', curr_file);
    temp = getSingleExptData(curr_file);
    if isempty(temp{4})
        disp(i);
        disp(listOfVariables);
%     else
%         disp(temp{4}(1));
    end
end

% Looks like a lot of the experiments just don't have the col radii in
% pixels saved. Will need to redo manually

%% Manually Fill in missing colrads pix


%% Fixing ALL Col Radii (cm)

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
for i = 1:length(allExptData)
    if ~isempty(allExptData(i).scanfolder)
        fprintf('Expt %g of %g\n', i, length(allExptData));
        tic;
        % load the resolution
        testinfo = imfinfo(fullfile(allExptData(i).scanfolder, allExptData(i).imfiles{1})); 
        if isfield(testinfo, 'XResolution')
            resolution = testinfo.XResolution;
        else
            disp(allExptData(i).exptfile);
            resp = input('Type the resolution in pixels: ');
            res = resp;
        end
        for j = 1:length(allExptData(i).imfiles)
            fprintf('Im %g of %g\n', j, length(allExptData(i).polfiles));
            if ~isempty(allExptData(i).scanfolder)
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

for i = 1:length(allExptData)
    if isempty(allExptData(i).slide_CVs)
        allExptData(i).slide_CVs = cell(1, length(allExptData(i).imfiles));
        allExptData(i).slide_CVs(:) = {NaN};
        allExptData(i).slide_CV_means = zeros(1, length(allExptData(i).imfiles));
    end
end

%% iterate over expts
winwidth = 10;
for i = 1:length(allExptData)
    if ~isempty(allExptData(i).scanfolder)
        fprintf('Expt %g of %g\n', i, length(allExptData));
        tic;
        % load the image
        for j = 1:length(allExptData(i).imfiles)
        fprintf('Im %g of %g\n', j, length(allExptData(i).polfiles));
            if isnan(allExptData(i).slide_CVs{j})
                polarim = imread(fullfile(allExptData(i).exptfolder, ...
                    allExptData(i).polfiles{j}));
                [cvs_all, cv_traj] = getRollingCV(polarim, winwidth);
                % Store
                allExptData(i).slide_CVs{j} = cv_traj;
%                 cv_traj = allExptData(i).slide_CVs{j};
                cv_traj(cv_traj==0) = NaN;
                allExptData(i).slide_CV_means(j) = nanmean(cv_traj);
            end
        end
        toc;
    end
end

%% Redoing sliding CV Mean Calculation using Colony Radii

winwidth = 10;
for i = 1:length(allExptData)
    if ~isempty(allExptData(i).scanfolder)
        fprintf('Expt %g of %g\n', i, length(allExptData));
        tic;
        % load the image
        for j = 1:length(allExptData(i).imfiles)
        fprintf('Im %g of %g\n', j, length(allExptData(i).polfiles));
            if isnan(allExptData(i).slide_CVs{j})
                polarim = imread(fullfile(allExptData(i).exptfolder, ...
                    allExptData(i).polfiles{j}));
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
for i = 1:length(allExptData)
    if isempty(allExptData(i).med_min_intens)
        allExptData(i).med_min_intens = NaN(1, length(allExptData(i).imfiles));
    end
end
%%
% iterate over every image
for i = 1:length(allExptData)

    fprintf('Expt %g of %g\n', i, length(allExptData));
    if ~isempty(allExptData(i).scanfolder)
        tic;
        for j = 1:length(allExptData(i).imfiles)
        fprintf('Im %g of %g\n', j, length(allExptData(i).polfiles));
            % check if already done
            if isnan(allExptData(i).med_min_intens(j))
                
                % Calculate the ft region mag
                polarim = imread(fullfile(allExptData(i).exptfolder, ...
                    allExptData(i).polfiles{j}));
                [colmask, colregion] = getColMask(polarim);
                colregion = colregion(:)';
                colregion(colregion==0) = NaN;
                minpix = mink(colregion, 50);
                allExptData(i).med_min_intens(j) = nanmedian(minpix);
            end
        end
        toc;
    end
end

%% Redoing Colony Mean Min Intensity Calculation

% Get the mean of the darkest 100 pixels for each colony, after masking
% out, & store
% First put in empty values
if ~isfield(allExptData, 'mean_min_intens')
    allExptData(1).mean_min_intens = [];
end
for i = 1:length(allExptData)
%     if isempty(allExptData(i).mean_min_intens)
        allExptData(i).mean_min_intens = NaN(1, length(allExptData(i).imfiles));
%     end
end
%%
% iterate over every image
for i = 1:length(allExptData)

    fprintf('Expt %g of %g\n', i, length(allExptData));
    if ~isempty(allExptData(i).scanfolder)
        tic;
        for j = 1:length(allExptData(i).imfiles)
        fprintf('Im %g of %g\n', j, length(allExptData(i).polfiles));
            % check if already done
            if isnan(allExptData(i).mean_min_intens(j))
                
                % Calculate the ft region mag
                polarim = imread(fullfile(allExptData(i).exptfolder, ...
                    allExptData(i).polfiles{j}));
                [colmask, colregion] = getColMask(polarim);
                colregion = colregion(:)';
                colregion(colregion==0) = NaN;
                minpix = mink(colregion, 50);
                allExptData(i).mean_min_intens(j) = nanmean(minpix);
            end
        end
        toc;
    end
end

%% Redoing Colony MEAN & MEDIAN Intensity Calculation

% Get the median and mean intens of each colony, after masking
% out, & store
% First put in empty values
if ~isfield(allExptData, 'med_intens')
    allExptData(1).med_intens = [];
end
if ~isfield(allExptData, 'mean_intens')
    allExptData(1).mean_intens = [];
end

for i = 1:length(allExptData)
    if isempty(allExptData(i).med_intens)
        allExptData(i).med_intens = NaN(1, length(allExptData(i).imfiles));
    end
    if isempty(allExptData(i).mean_intens)
        allExptData(i).mean_intens = NaN(1, length(allExptData(i).imfiles));
    end
end
%%
% iterate over every image
for i = 1:length(allExptData)

    fprintf('Expt %g of %g\n', i, length(allExptData));
    if ~isempty(allExptData(i).scanfolder)
        tic;
        for j = 1:length(allExptData(i).imfiles)
        fprintf('Im %g of %g\n', j, length(allExptData(i).polfiles));
            % check if already done
            if isnan(allExptData(i).mean_intens(j))
                
                % Calculate the ft region mag
                polarim = imread(fullfile(allExptData(i).exptfolder, ...
                    allExptData(i).polfiles{j}));
                [colmask, colregion] = getColMaskPolar(polarim);
                colregion = colregion(:)';
                colregion(colregion==0) = NaN;
                colregion(colregion==1) = NaN;
                allExptData(i).med_intens(j) = nanmedian(colregion);
                allExptData(i).mean_intens(j) = nanmean(colregion);
            end
        end
        toc;
    end
end

%% 09/01/21 Trying Out FT Xform of 1D Radial Trajs

for i = 1:27
    if isempty(allExptData(i).scanfolder)
        continue;
    end
    disp(i);
    for j = 1:length(allExptData(i).imfiles)

        temptraj = allExptData(i).avgd{j}';
        tempdist = allExptData(i).dists{j};

        % Instead let's try to set the end of colony based on COL RADIUS
        temprad = allExptData(i).colrads_cm_redo(j);
        lastind = find(tempdist>temprad, 1);
        
        % now find inoculum border
        [~, inocedge] = min(temptraj(1:200));
        
        temptraj = temptraj(inocedge:lastind);
        tempdist = tempdist(inocedge:lastind);
        
        temptraj = temptraj-mean(temptraj);

        temptraj = -temptraj;

        L = length(temptraj);
        Y = fft(temptraj);

        P2 = abs(Y/L);
        P1 = P2(1:(L/2 + 1));
        P1(2:end-1) = 2*P1(2:end-1);
        % To get sampling frequency, we will use the 'dist' vector which tell
        % us where each of the measurements in the 1000 length vector lands in
        % terms of cm

        % Divide the # msmts over length in mm to get # msmts per 1 mm
        msmt_freq = length(temptraj)/((tempdist(end)-tempdist(1))*10); 

        Fs = msmt_freq;
        f = Fs*(0:(L/2))/L;
        
        % Disregard the first two max frequencies as they represent the
        % colony itself (ie 38 mm ring width...)
        [maxval, maxind] = max(P1(3:end));
        maxfreq = f(maxind+2);

        % Let's try zeroing out all but a certain few freqs
        [maxvals, maxinds] = maxk(Y, 10);
        test_max_freqs = zeros(1, L);
        test_max_freqs(maxinds) = Y(maxinds);

        % can we use find peaks on this and get mean peak distance?
        newvec = ifft(test_max_freqs);
        [pks, locs] = findpeaks(newvec);
        locs_cm = tempdist(locs);
        meanpeakdist = mean(diff(locs_cm))*10;
        
        % Store the max freq
        allExptData(i).max_freqs(j) = maxfreq;
        allExptData(i).meanpeakdist(j) = meanpeakdist;
    end
end

%% 09/08/21 Getting FT Values as Percent of Total for Disk Region

% Use a disk shaped region of increasing radius; get the total of the
% magnitude of FT xform in that central region; divide by total of the FT
% xform (sum) to see how much of the frequency is located in that region;
% then later we can try comparing that vs iptg vs gene to see if any have a
% change or a trend
x = 1000;
y = 1000;
for exptnum = 1:length(allExptData)
    fprintf('Expt %g of %g\n', exptnum, length(allExptData));
    
%     if size(allExptData(exptnum).FT_weightvals, 1) == length(allExptData(exptnum).imfiles)
%         continue;
%     end
    if isempty(allExptData(exptnum).scanfolder)
        continue;
    end
    
    
    allweightvals = zeros(length(allExptData(exptnum).imfiles), 20);

    for j = 1:length(allExptData(exptnum).imfiles)
        disp(j);
        % Calculate the ft region mag
        polarim = imread(fullfile(allExptData(exptnum).exptfolder, ...
            allExptData(exptnum).polfiles{j}));
        [colmask, colregion] = getColMask(polarim);
        fftim = fftshift(fft2(colregion));
        mag_fftim = abs(fftim);

        weightvals = zeros(1, 20);
        count = 1;
%         for k = 1:10:200
        for k = 1:20:400
            rval2 = k;
            mask2 = fspecial('disk', rval2) == 0;
            mask2 = imresize(padarray(mask2, [floor((x/2)-rval2) floor((y/2)-rval2)], 1, 'both'), [x y]);

            % Let's calculate how much of the values are contained in this
              vals = mag_fftim(~mask2);
              weightval = sum(vals)/(sum(sum(mag_fftim)));
              weightvals(count) = weightval;
              count = count+1;
        end
        allweightvals(j, :) = weightvals;

    end
%     allExptData(exptnum).FT_weightvals = allweightvals;
    allExptData(exptnum).FT_weightvals_wide = allweightvals;
end

%% 10/20/21 Getting Inoc Edge Mean Intens for all Images
for exptnum = 1:length(allExptData)
    fprintf('Expt %g of %g\n', exptnum, length(allExptData));
    if isempty(allExptData(exptnum).scanfolder)
        continue;
    end
    for j = 1:length(allExptData(exptnum).imfiles)
        disp(j);
        % Calculate the ft region mag
        polarim = imread(fullfile(allExptData(exptnum).exptfolder, ...
            allExptData(exptnum).polfiles{j}));
        tempim = im2double(polarim);
        inocedge_mean = getInocEdgeMean(tempim);
        allExptData(exptnum).inoc_edge_means(j) = inocedge_mean;
    end
end
        
%% 10/21/21 Try getting a range to reduce effects of outliers
for exptnum = 1:length(allExptData)
    fprintf('Expt %g of %g\n', exptnum, length(allExptData));
    if isempty(allExptData(exptnum).scanfolder)
        continue;
    end
    for j = 1:length(allExptData(exptnum).imfiles)
        disp(j);
        % Calculate the ft region mag
        polarim = imread(fullfile(allExptData(exptnum).exptfolder, ...
            allExptData(exptnum).polfiles{j}));
        tempim = im2double(polarim);
        inocedge_mean = getInocEdgeMean_Rng(tempim);
        allExptData(exptnum).inoc_edge_rng_means(j) = inocedge_mean;
    end
end

%% 10/26/21 Trying a Vertical CV
winheight = 15;
for tempexptnum = 1:length(allExptData)
    if ~isempty(allExptData(tempexptnum).scanfolder)
        disp(tempexptnum);
        roll_cvs = zeros(1, length(allExptData(tempexptnum).imfiles));
        for i = 1:length(allExptData(tempexptnum).imfiles)
            tempim = im2double(imread(fullfile(allExptData(tempexptnum).exptfolder, ...
                allExptData(tempexptnum).polfiles{i})));    
            temprad = allExptData(tempexptnum).colrads_cm_redo(i);
            tempdist = allExptData(tempexptnum).dists{i};
            [cvs_all, cv_traj] = getVertCV(tempim, winheight);
            for j = 1:1000
                [inocedge, lastind] = cropTraj(tempim(:, j), temprad, tempdist);
                cvs_all(1:inocedge, j) = NaN;
                cvs_all(lastind:end, j) = NaN;
            end
            roll_cvs(i) = nanmean(nanmean(cvs_all));
        end
        allExptData(tempexptnum).vert_cvs15 = roll_cvs;
    end
end

%% Try getting the Peaks of the Gradients
for tempexptnum = 1:length(allExptData)
    if ~isempty(allExptData(tempexptnum).scanfolder)
        disp(tempexptnum);
        roll_cvs = zeros(1, length(allExptData(tempexptnum).imfiles));
        for i = 1:length(allExptData(tempexptnum).imfiles)
            tempim = im2double(imread(fullfile(allExptData(tempexptnum).exptfolder, ...
                allExptData(tempexptnum).polfiles{i})));    
            temprad = allExptData(tempexptnum).colrads_cm_redo(i);
            tempdist = allExptData(tempexptnum).dists{i};
            % Now get the gradient
            [~, gradim] = gradient(tempim); gradim = abs(gradim);
            % Crop each column
            for j = 1:1000
                [inocedge, lastind] = cropTraj(tempim(:, j), temprad, tempdist);
                gradim(1:(inocedge+20), j) = NaN;
                gradim(lastind:end, j) = NaN;
            end
            % Get the trajectory, then smooth, then find the peaks
            tempgrad_traj = nanmean(gradim, 2);
            tempmean = nanmean(tempgrad_traj);
            tempstd = nanstd(tempgrad_traj);
            tempgrad_traj_smooth = movmean(tempgrad_traj, 10);
            [pks, locs] = findpeaks(tempgrad_traj_smooth, 'MinPeakHeight', tempmean+0.3*tempstd);
            pksmean = mean(pks);
            roll_cvs(i) = pksmean;
            if isempty(pks)
                roll_cvs(i) = NaN;
            end
        end
        allExptData(tempexptnum).grad_pks = roll_cvs;
    end
end

%% Let's try getting peak gradient near the inoculum (ie how distinct inoculum edge is)
for tempexptnum = 1:length(allExptData)
    if ~isempty(allExptData(tempexptnum).scanfolder)
        disp(tempexptnum);
        for i = 1:length(allExptData(tempexptnum).imfiles)
            temptraj = allExptData(tempexptnum).avgd{i};
            peak_inoc_grad = getPeakInocEdgeGrad(temptraj);
            allExptData(tempexptnum).inoc_pk_grads(i) = peak_inoc_grad;
        end
    end
end

%% Get the peak gradient but columnwise to reduce noise
mean_width = 25;
for tempexptnum = 1:length(allExptData)
    if ~isempty(allExptData(tempexptnum).scanfolder)
        disp(tempexptnum);
        for i = 1:length(allExptData(tempexptnum).imfiles)
            tempim = im2double(imread(fullfile(allExptData(tempexptnum).exptfolder, ...
            allExptData(tempexptnum).polfiles{i})));  
            peak_inoc_grad = getColumnPeakInocEdgeGrad(tempim, mean_width);
            allExptData(tempexptnum).inoc_pk_grads_col(i) = peak_inoc_grad;
        end
    end
end

%% Get the max-min inoc edge gradient but columnwise to reduce noise
mean_width = 25;
for tempexptnum = 1:length(allExptData)
    if ~isempty(allExptData(tempexptnum).scanfolder)
        disp(tempexptnum);
        for i = 1:length(allExptData(tempexptnum).imfiles)
            tempim = im2double(imread(fullfile(allExptData(tempexptnum).exptfolder, ...
            allExptData(tempexptnum).polfiles{i})));  
            inoc_grad_rng = getColumnInocEdgeGradRng(tempim, mean_width);
            allExptData(tempexptnum).inoc_pk_rngs_col(i) = inoc_grad_rng;
        end
    end
end

%% Get the max-min inoc edge of the TRAJ but columnwise to reduce noise
mean_width = 25;
for tempexptnum = 1:length(allExptData)
    if ~isempty(allExptData(tempexptnum).scanfolder)
        disp(tempexptnum);
        for i = 1:length(allExptData(tempexptnum).imfiles)
            tempim = im2double(imread(fullfile(allExptData(tempexptnum).exptfolder, ...
            allExptData(tempexptnum).polfiles{i})));  
            inoc_traj_rng = getColumnInocEdgeRng(tempim, mean_width);
            allExptData(tempexptnum).inoc_traj_rngs_col(i) = inoc_traj_rng;
        end
    end
end


%% SAVE Updated Single Gene Analysis
% Save updated version
disp('Saving updated version...');
filename = strcat('all_single_gene_expt_data_', date, '.mat');
save(filename, 'allExptData', '-v7.3');
disp('Done.');









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS
%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%
% getting col mask
function [colmask, colregion] = getColMask(polarim)
    startim = polarim;
    tempim = im2double(startim);
    rim_pix = startim<0.7;
    rim_pix(1:400, :) = false;
    bg_val = 0.885;
    tempim(rim_pix) = bg_val; 
%     tempim(1:200, :) = startim(1:200, :);
    % fill in white space with 'agar'
    
    plate_bg = startim==1;
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
    yikes = imbinarize(rescale(combine_im+diffim));

    % Dilate image a little
    se = strel("disk", 6, 4);
    yikes = imdilate(yikes, se);


    % Use morphological operations for cleanup/filling
    se = strel("disk", 3, 4);
    morphim1 = imopen(yikes, se);
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


end

function regionmag_array = getFTRegionMags(ftim)
    % Iterate over the ft image (magnitude image) first quadrant. Using a
    % specific patch size & stride, get the mean magnitude in each region.
    % Store the output in a new array/image
    start_ind = 80;
    end_ind = 480;
    patch_size = 50;
    num_patches = (end_ind-start_ind)/patch_size;
    regionmag_array = zeros(num_patches, num_patches);
    % Iterate over the rows
    for i = start_ind:patch_size:(end_ind-patch_size)
        % Iterate over the columns
        for j = start_ind:patch_size:(end_ind-patch_size)
            tempregion = ftim(i:i+patch_size, j:j+patch_size);
            regionmag_array(((i-start_ind)/patch_size)+1, ((j-start_ind)/patch_size)+1) = ...
                mean(mean(tempregion));
        end
    end
end

function inocedge_mean = getInocEdgeMean(tempim)
    im_width = size(tempim, 2);
    inoc_vals = zeros(1, im_width);
    for i = 1:im_width
        inoc_vals(i) = min(tempim(20:100, i));
    end
    inocedge_mean = mean(inoc_vals);

end

function inocedge_mean = getInocEdgeMean_Rng(tempim)
    im_width = size(tempim, 2);
    inoc_vals = zeros(1, im_width);
    for i = 1:im_width
        [~, inocloc] = min(tempim(20:100, i));
        inocloc = inocloc+19;
        tempvals = tempim((inocloc-10):(inocloc+10), i);
        inoc_vals(i) = mean(tempvals);
    end
    inocedge_mean = mean(inoc_vals);

end

function [inocedge, lastind] = cropTraj(temptraj, temprad, tempdist)
    lastind = find(tempdist>temprad, 1);
    % if temptraj = 1 cut it off there
    if any(temptraj==1)
        last_oneind = find(temptraj==1, 1)-10;
        lastind = min(lastind, last_oneind);
    end
    % now find inoculum border
    [~, inocedge] = min(temptraj(25:200));
    inocedge = inocedge+24;
end

function [cvs_all, cv_traj] = getVertCV(polarim, winheight)
    % This function, given a window width to use, will get windows of that
    % sliding from left to right for each row of the image and get the CV
    % of each window (ie, a local CV).
    polarim = im2double(polarim);
    
    cvs_all = zeros(size(polarim));
    % get movmean and movstd
    stdim = movstd(polarim, winheight, 0, 2);
    meanim = movmean(polarim, winheight, 2);
    cvs_all = stdim./meanim;

    cv_traj = mean(cvs_all, 2);

end

function peak_inoc_grad = getPeakInocEdgeGrad(temptraj)
    temptraj = temptraj(1:500);
    tempgradtraj = gradient(temptraj);
    peak_inoc_grad = max(tempgradtraj(25:60));    
end

function peak_inoc_grad = getColumnPeakInocEdgeGrad(tempim, mean_width)
    meanim = movmean(tempim, mean_width, 2);
    peak_inoc_grads = zeros(1, size(meanim, 2));
    for j = 1:size(meanim, 2)
        temptraj = meanim(1:500, j);
        tempgradtraj = gradient(temptraj);
        peak_inoc_grads(j) = max(tempgradtraj(15:60));  
    end
    peak_inoc_grad = mean(peak_inoc_grads);  
end

function inoc_grad_rng = getColumnInocEdgeGradRng(tempim, mean_width)
    meanim = movmean(tempim, mean_width, 2);
    peak_inoc_grads = zeros(1, size(meanim, 2));
    for j = 1:size(meanim, 2)
        temptraj = meanim(1:500, j);
        tempgradtraj = gradient(temptraj);
        inoc_grad_rngs(j) = max(tempgradtraj(15:60))-min(tempgradtraj(15:60));  
    end
    inoc_grad_rng = mean(inoc_grad_rngs);
end

function inoc_traj_rng = getColumnInocEdgeRng(tempim, mean_width)
    meanim = movmean(tempim, mean_width, 2);
    peak_inoc_grads = zeros(1, size(meanim, 2));
    for j = 1:size(meanim, 2)
        temptraj = meanim(1:500, j);
        inoc_traj_rngs(j) = max(temptraj(15:60))-min(temptraj(15:60));  
    end
    inoc_traj_rng = mean(inoc_traj_rngs);
end